// Copyright (c) 2026 Alejandro Gonzalez-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! `collapse chain`: transcript models from segmented reads by intron chain (ARK PLAN.md §5).
//!
//! Junction wobble is corrected per strand, reads are grouped by intron chain, truncated
//! subchains are absorbed into their best-supported parent, mono-exonic reads are absorbed,
//! excluded as intronic or clustered, and models pass support, tail and read-fraction
//! thresholds. Every model is a real read (its representative) tagged `#CN<n_reads>`.

use std::cmp::Reverse;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::cli::{ChainArgs, Preset};
use crate::utils::Dsu;

/// Half-open genomic interval (intron or exon).
type Iv = (u32, u32);

/// Chrom, name, strand (true = +), start, end and introns of a BED12 line.
type Bed<'a> = (&'a str, &'a str, bool, u32, u32, Vec<Iv>);

const OUTPUTS: [&str; 5] = [
    "models.bed",
    "support.tsv",
    "members.tsv",
    "excluded.bed",
    "counts.tsv",
];

/// A segmented read: its BED12 line and the iso-segment tags the caller uses.
struct Read<'a> {
    line: &'a str,
    chrom: &'a str,
    name: &'a str,
    plus: bool,
    start: u32,
    end: u32,
    introns: Vec<Iv>,
    corrected: bool,
    pa: u32,
    pr: u32,
    tc: u32,
    iy: u32,
    fg: bool,
}

impl Read<'_> {
    fn tail(&self) -> bool {
        self.pa >= 20
    }

    /// No tail and more than 5 genomic A's (iso-pas: PR - (TC + PA)).
    fn intra(&self) -> bool {
        !self.tail() && self.pr.saturating_sub(self.tc.saturating_add(self.pa)) > 5
    }

    fn five(&self) -> u32 {
        if self.plus {
            self.start
        } else {
            self.end
        }
    }

    fn three(&self) -> u32 {
        if self.plus {
            self.end
        } else {
            self.start
        }
    }
}

/// Reference transcripts of one (chrom, strand).
#[derive(Default)]
struct Ref {
    starts: Vec<u32>,
    ends: Vec<u32>,
    chains: HashSet<Vec<Iv>>,
    /// first intron in transcript orientation -> 5′ ends
    first: HashMap<Iv, Vec<u32>>,
    /// single-exon transcripts, sorted
    mono: Vec<(u32, u32, usize)>,
}

/// A root chain with the reads absorbed into it, or a mono-exonic cluster.
struct Model<'a> {
    chain: &'a [Iv],
    full: Vec<usize>,
    partial: Vec<usize>,
    mol: usize,
    five: u32,
    three: u32,
    rep: usize,
    known: bool,
    protected: bool,
    reason: Option<&'static str>,
}

impl Model<'_> {
    fn n_reads(&self) -> usize {
        self.full.len() + self.partial.len()
    }

    fn category(&self) -> &'static str {
        if self.known {
            "known"
        } else if self.chain.is_empty() {
            "mono"
        } else {
            "novel"
        }
    }
}

/// Runs `collapse chain`, writing `<prefix>.{models.bed,support.tsv,members.tsv,excluded.bed,counts.tsv}`.
pub fn run(args: ChainArgs) -> Result<(), Box<dyn Error>> {
    // ponytail: inputs are held in memory like `collapse run`; stream per chromosome if a
    // [sample, chr] ever outgrows RAM
    let read = |p: &Path| fs::read_to_string(p).map_err(|e| format!("{}: {e}", p.display()));
    let texts = args
        .bed
        .iter()
        .map(|p| read(p))
        .collect::<Result<Vec<_>, _>>()?;
    let reference = read(&args.reference)?;
    let beds: Vec<(String, &str)> = args
        .bed
        .iter()
        .map(|p| p.display().to_string())
        .zip(texts.iter().map(String::as_str))
        .collect();

    let mut out = OUTPUTS
        .iter()
        .map(|ext| File::create(format!("{}.{ext}", args.prefix)).map(BufWriter::new))
        .collect::<Result<Vec<_>, _>>()?;
    let ref_path = args.reference.display().to_string();
    chain(&beds, (&ref_path, &reference), &args, &mut out)?;
    for w in &mut out {
        w.flush()?;
    }

    Ok(())
}

/// The caller (PLAN.md §5.2 steps 1-10) over in-memory BED texts; `out` follows `OUTPUTS`.
fn chain<W: Write>(
    beds: &[(String, &str)],
    (ref_path, ref_text): (&str, &str),
    args: &ChainArgs,
    out: &mut [W],
) -> Result<(), Box<dyn Error>> {
    let [models_w, support_w, members_w, excluded_w, counts_w] = out else {
        unreachable!("one writer per output");
    };

    // presets (§5.3): novel molecules, read fraction, mono molecules, tail-switch read rate
    let (novel, fraction, mono_min, tail_min) = match args.preset {
        Preset::Sensitive => (1, 1.0, 3, f64::INFINITY),
        Preset::Balanced => (2, 0.99, 10, 0.7),
        Preset::Strict => (3, 0.95, 10, 0.1),
    };
    let novel = args.min_support_novel.unwrap_or(novel);
    let mono_min = args.min_support_mono.unwrap_or(mono_min);
    let fraction = args.min_read_fraction.unwrap_or(fraction);
    if !(fraction > 0.0 && fraction <= 1.0) {
        return Err(format!("--min-read-fraction must be in (0, 1], got {fraction}").into());
    }

    // 1. parse; tags are key-based: split once on `__`, then on `#`, two-letter keys
    let mut reads = Vec::new();
    for (path, text) in beds {
        for (n, line) in text.lines().enumerate() {
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let (chrom, name, plus, start, end, introns) =
                parse_bed(line).map_err(|e| format!("{path}:{}: {e}", n + 1))?;
            let (mut pa, mut pr, mut tc, mut iy, mut fg) = (0, 0, 0, 0, false);
            for tag in name.split_once("__").map_or("", |t| t.1).split('#') {
                let v = tag.get(2..).and_then(|v| v.parse().ok()).unwrap_or(0);
                match tag.get(..2) {
                    Some("PA") => pa = v,
                    Some("PR") => pr = v,
                    Some("TC") => tc = v,
                    Some("IY") => iy = v,
                    Some("FG") => fg = true,
                    _ => {}
                }
            }
            reads.push(Read {
                line,
                chrom,
                name,
                plus,
                start,
                end,
                introns,
                corrected: false,
                pa,
                pr,
                tc,
                iy,
                fg,
            });
        }
    }

    // reference transcripts on the reads' chromosomes only
    let chroms: HashSet<&str> = reads.iter().map(|r| r.chrom).collect();
    let mut refs: HashMap<(&str, bool), Ref> = HashMap::new();
    for (n, line) in ref_text.lines().enumerate() {
        if !chroms.contains(line.split('\t').next().unwrap_or_default()) {
            continue;
        }
        let (chrom, _, plus, start, end, introns) =
            parse_bed(line).map_err(|e| format!("{ref_path}:{}: {e}", n + 1))?;
        let r = refs.entry((chrom, plus)).or_default();
        match (introns.first(), introns.last()) {
            (Some(&a), Some(&b)) => {
                r.starts.extend(introns.iter().map(|i| i.0));
                r.ends.extend(introns.iter().map(|i| i.1));
                let (first, five) = if plus { (a, start) } else { (b, end) };
                r.first.entry(first).or_default().push(five);
                r.chains.insert(introns);
            }
            _ => r.mono.push((start, end, 0)),
        }
    }
    for r in refs.values_mut() {
        r.starts.sort_unstable();
        r.starts.dedup();
        r.ends.sort_unstable();
        r.ends.dedup();
        r.mono.sort_unstable();
    }
    let none = Ref::default();

    // tail switch (step 7): on when enough of the input reads carry a polyA tail
    let tail_rate = reads.iter().filter(|r| r.tail()).count() as f64 / reads.len().max(1) as f64;
    let tail_on = tail_rate >= tail_min;

    let mut buckets: BTreeMap<(&str, bool), Vec<usize>> = BTreeMap::new();
    for (i, r) in reads.iter().enumerate() {
        buckets.entry((r.chrom, r.plus)).or_default().push(i);
    }

    // 2. junction correction per strand; intron starts and ends are separate site classes
    let (mut sites_moved, mut reads_moved, mut reads_kept) = (0, 0, 0);
    for (key, idx) in &buckets {
        let rf = refs.get(key).unwrap_or(&none);
        let (mut starts, mut ends) = (HashMap::new(), HashMap::new());
        for &i in idx {
            for &(s, e) in &reads[i].introns {
                *starts.entry(s).or_insert(0) += 1;
                *ends.entry(e).or_insert(0) += 1;
            }
        }
        let w = args.junction_wobble;
        let (ms, me) = (moves(&starts, &rf.starts, w), moves(&ends, &rf.ends, w));
        sites_moved += ms.len() + me.len();
        if ms.is_empty() && me.is_empty() {
            continue;
        }
        for &i in idx {
            let r = &mut reads[i];
            let new: Vec<Iv> = r
                .introns
                .iter()
                .map(|&(s, e)| (*ms.get(&s).unwrap_or(&s), *me.get(&e).unwrap_or(&e)))
                .collect();
            if new == r.introns {
                continue;
            }
            // a correction that would empty an exon (or intron) keeps the read's own blocks
            if new
                .iter()
                .chain(&exons(r.start, &new, r.end))
                .all(|x| x.0 < x.1)
            {
                (r.introns, r.corrected) = (new, true);
                reads_moved += 1;
            } else {
                reads_kept += 1;
            }
        }
    }

    let (mut n_chains, mut n_absorbed) = (0, 0);
    let mut all: Vec<Model> = Vec::new();
    for (key, idx) in &buckets {
        let (rf, plus) = (refs.get(key).unwrap_or(&none), key.1);

        // 3-4. group by (strand, corrected chain); molecules are counted in `model`
        let mut by_chain: HashMap<&[Iv], Vec<usize>> = HashMap::new();
        let mut mono = Vec::new();
        for &i in idx {
            if reads[i].introns.is_empty() {
                mono.push(i);
            } else {
                by_chain.entry(&reads[i].introns[..]).or_default().push(i);
            }
        }
        let mut chains: Vec<_> = by_chain.into_iter().collect();
        chains.sort_unstable_by(|a, b| a.0.cmp(b.0));
        n_chains += chains.len();
        let mut ms: Vec<Model> = chains
            .into_iter()
            .map(|(chain, full)| Model {
                known: rf.chains.contains(chain),
                ..model(&reads, chain, full)
            })
            .collect();

        // 5. absorb subchains, most introns first, into the best-supported unabsorbed parent
        let mut by_intron: HashMap<Iv, Vec<usize>> = HashMap::new();
        for (c, m) in ms.iter().enumerate() {
            for &iv in m.chain {
                by_intron.entry(iv).or_default().push(c);
            }
        }
        let mut order: Vec<usize> = (0..ms.len()).collect();
        order.sort_by_key(|&c| Reverse(ms[c].chain.len()));
        let mut absorbed = vec![false; ms.len()];
        for a in order {
            let (ac, k) = (ms[a].chain, ms[a].chain.len());
            let (left, right) = if plus {
                (ms[a].five, ms[a].three)
            } else {
                (ms[a].three, ms[a].five)
            };
            let tail = ms[a].full.iter().any(|&i| reads[i].tail());
            let (mut best, mut blocked, mut protect) = (None::<usize>, false, None);
            // ponytail: parents come from A's rarest intron; quadratic only for a locus of
            // thousands of truncations that share every intron
            let cands = ac
                .iter()
                .map(|iv| &by_intron[iv])
                .min_by_key(|v| v.len())
                .expect("multi-exon chain");
            for &b in cands {
                let (bc, m) = (ms[b].chain, ms[b].chain.len());
                if m <= k || absorbed[b] {
                    continue;
                }
                let Ok(p) = bc.binary_search(&ac[0]) else {
                    continue;
                };
                // A's introns are a contiguous run of B's and its ends stay inside B's exons
                if p + k > m
                    || bc[p..p + k] != *ac
                    || (p > 0 && left < bc[p - 1].1)
                    || (p + k < m && right > bc[p + k].0)
                {
                    continue;
                }
                let (five_cut, three_inner) = if plus {
                    (p > 0, p + k < m)
                } else {
                    (p + k < m, p > 0)
                };
                if three_inner && tail {
                    continue;
                }
                if five_cut && *protect.get_or_insert_with(|| tss(&reads, &ms[a], rf, plus)) {
                    blocked = true;
                    continue;
                }
                if best.is_none_or(|x| ms[b].n_reads() > ms[x].n_reads()) {
                    best = Some(b);
                }
            }
            match best {
                Some(b) => {
                    absorbed[a] = true;
                    n_absorbed += 1;
                    let full = std::mem::take(&mut ms[a].full);
                    ms[b].mol += ms[a].mol;
                    ms[b].partial.extend(full);
                }
                None => ms[a].protected = blocked,
            }
        }
        let mut ms: Vec<Model> = ms
            .into_iter()
            .zip(absorbed)
            .filter_map(|(m, a)| (!a).then_some(m))
            .collect();

        // 7. support and tail thresholds; a known chain needs one molecule
        for m in ms.iter_mut().filter(|m| !m.known) {
            if m.mol < novel {
                m.reason = Some("support");
            } else if tail_on && !m.full.iter().any(|&i| reads[i].tail()) {
                m.reason = Some("tail");
            }
        }

        // 6. mono-exonic reads: inside an exon of a kept model -> absorbed by the
        // best-supported one, inside one of its introns -> intronic, else overlap clusters.
        // Terminal exons reach the model's outermost reads, so end jitter still absorbs.
        mono.sort_by_key(|&i| (reads[i].start, reads[i].end));
        let spans: Vec<Iv> = mono
            .iter()
            .map(|&i| (reads[i].start, reads[i].end))
            .collect();
        let mut rank: Vec<usize> = (0..ms.len()).filter(|&j| ms[j].reason.is_none()).collect();
        rank.sort_by_key(|&j| Reverse(ms[j].n_reads()));
        let (mut ex, mut gaps) = (Vec::new(), Vec::new());
        for (r, &j) in rank.iter().enumerate() {
            let m = &ms[j];
            let members = || m.full.iter().chain(&m.partial).map(|&i| &reads[i]);
            let lo = members().map(|x| x.start).min().unwrap_or(0);
            let hi = members().map(|x| x.end).max().unwrap_or(0);
            ex.extend(exons(lo, m.chain, hi).into_iter().map(|(s, e)| (s, e, r)));
            gaps.extend(m.chain.iter().map(|&(s, e)| (s, e, r)));
        }
        ex.sort_unstable();
        gaps.sort_unstable();
        let (mut intronic, mut free, mut mol) = (Vec::new(), Vec::new(), Vec::new());
        let hits = contained(&ex, &spans)
            .into_iter()
            .zip(contained(&gaps, &spans));
        for (&i, hit) in mono.iter().zip(hits) {
            match hit {
                (Some(r), _) => {
                    ms[rank[r]].partial.push(i);
                    mol.push((rank[r], reads[i].start, reads[i].end));
                }
                (None, Some(_)) => intronic.push(i),
                (None, None) => free.push(i),
            }
        }
        mol.sort_unstable();
        mol.dedup();
        for (j, _, _) in mol {
            ms[j].mol += 1;
        }
        for cl in clusters(&reads, &intronic) {
            ms.push(Model {
                reason: Some("intronic"),
                ..model(&reads, &[], cl)
            });
        }
        let cls = clusters(&reads, &free);
        let spans: Vec<Iv> = cls
            .iter()
            .map(|c| {
                (
                    reads[c[0]].start,
                    c.iter().map(|&i| reads[i].end).max().unwrap_or(0),
                )
            })
            .collect();
        for (cl, inside_ref) in cls.into_iter().zip(contained(&rf.mono, &spans)) {
            let m = model(&reads, &[], cl);
            let n = m.full.len();
            let intra = m.full.iter().filter(|&&i| reads[i].intra()).count();
            let tails = m.full.iter().filter(|&&i| reads[i].tail()).count();
            // 9. intra-priming (like the tail rule) applies only when tails are present
            let reason = if inside_ref.is_some() {
                None
            } else if m.mol < mono_min {
                Some("support")
            } else if tail_on && 2 * intra > n {
                Some("intraprimed")
            } else if tail_on && 2 * tails < n {
                Some("tail")
            } else {
                None
            };
            ms.push(Model {
                known: inside_ref.is_some(),
                reason,
                ..m
            });
        }

        // 7. read fraction per locus (kept models sharing >= 1 exonic bp); known are exempt
        let kept: Vec<usize> = (0..ms.len()).filter(|&j| ms[j].reason.is_none()).collect();
        let mut ex = Vec::new();
        for (j, &m) in kept.iter().enumerate() {
            let rep = &reads[ms[m].rep];
            ex.extend(
                exons(rep.start, &rep.introns, rep.end)
                    .into_iter()
                    .map(|(s, e)| (s, e, j)),
            );
        }
        ex.sort_unstable();
        let mut dsu = Dsu::new(kept.len());
        let (mut end, mut cur) = (0, 0);
        for &(s, e, j) in &ex {
            if s < end {
                dsu.union(j, cur);
            }
            if e > end {
                (end, cur) = (e, j);
            }
        }
        let mut loci: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for (j, &m) in kept.iter().enumerate() {
            loci.entry(dsu.find(j)).or_default().push(m);
        }
        for mut locus in loci.into_values() {
            let total = locus.iter().map(|&m| ms[m].n_reads()).sum::<usize>() as f64;
            locus.sort_by_key(|&m| Reverse(ms[m].n_reads()));
            let mut cum = 0;
            for m in locus {
                if !ms[m].known && cum as f64 >= fraction * total {
                    ms[m].reason = Some("fraction");
                }
                cum += ms[m].n_reads();
            }
        }
        all.extend(ms);
    }

    // 10. output sorted by position, ties by strand then representative name
    let (mut kept, mut excl): (Vec<&Model>, Vec<&Model>) =
        all.iter().partition(|m| m.reason.is_none());
    let key = |m: &&Model| {
        let r = &reads[m.rep];
        (r.chrom, r.start, r.end, !r.plus, r.name)
    };
    kept.sort_by_key(key);
    excl.sort_by_key(key);

    let prefix = &args.prefix;
    let mut fate: Vec<Result<usize, &str>> = vec![Err(""); reads.len()];
    writeln!(
        support_w,
        "model\trep\tstrand\tchain\tcategory\tn_reads\tn_molecules\tn_partial\tn_tail\tfive_prime\tthree_prime\tflags"
    )?;
    for (n, m) in kept.iter().enumerate() {
        let r = &reads[m.rep];
        let members = || m.full.iter().chain(&m.partial);
        members().for_each(|&i| fate[i] = Ok(n + 1));
        writeln!(models_w, "{}", bed_line(r, m.n_reads()))?;
        let chain: Vec<String> = m.chain.iter().map(|(s, e)| format!("{s}-{e}")).collect();
        let flags: Vec<&str> = [(r.corrected, "corrected"), (m.protected, "protected")]
            .into_iter()
            .filter_map(|(on, flag)| on.then_some(flag))
            .collect();
        let or_dot = |v: String| if v.is_empty() { ".".into() } else { v };
        writeln!(
            support_w,
            "{prefix}.{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            n + 1,
            r.name,
            if r.plus { '+' } else { '-' },
            or_dot(chain.join(",")),
            m.category(),
            m.n_reads(),
            m.mol,
            m.partial.len(),
            members().filter(|&&i| reads[i].tail()).count(),
            m.five,
            m.three,
            or_dot(flags.join(",")),
        )?;
    }
    for m in &excl {
        let why = m.reason.unwrap_or_default();
        m.full
            .iter()
            .chain(&m.partial)
            .for_each(|&i| fate[i] = Err(why));
        writeln!(
            excluded_w,
            "{}\t{why}",
            bed_line(&reads[m.rep], m.n_reads())
        )?;
    }
    writeln!(members_w, "read\tmodel")?;
    for (r, f) in reads.iter().zip(&fate) {
        match f {
            Ok(n) => writeln!(members_w, "{}\t{prefix}.{n}", r.name)?,
            Err(why) => writeln!(members_w, "{}\texcluded:{why}", r.name)?,
        }
    }

    let excluded = |why: &str| fate.iter().filter(|f| **f == Err(why)).count();
    let models = |c: &str| kept.iter().filter(|m| m.category() == c).count();
    let multi = reads.iter().filter(|r| !r.introns.is_empty()).count();
    let steps = [
        ("reads_in", reads.len()),
        ("reads_multi", multi),
        ("reads_mono", reads.len() - multi),
        ("sites_corrected", sites_moved),
        ("reads_corrected", reads_moved),
        ("reads_correction_skipped", reads_kept),
        ("chains", n_chains),
        ("chains_absorbed", n_absorbed),
        ("models_known", models("known")),
        ("models_novel", models("novel")),
        ("models_mono", models("mono")),
        ("excluded_support", excluded("support")),
        ("excluded_fraction", excluded("fraction")),
        ("excluded_tail", excluded("tail")),
        ("excluded_intronic", excluded("intronic")),
        ("excluded_intraprimed", excluded("intraprimed")),
        ("models_out", kept.len()),
    ];
    for (step, v) in steps {
        writeln!(counts_w, "{step}\t{v}")?;
    }
    writeln!(counts_w, "tail_rate\t{tail_rate:.4}")?;
    log::info!(
        "INFO: chain: {} reads -> {} models, {} excluded chains/clusters",
        reads.len(),
        kept.len(),
        excl.len()
    );

    Ok(())
}

/// Parses a BED12 line; columns past 12 are ignored.
fn parse_bed(line: &str) -> Result<Bed<'_>, &'static str> {
    let f: Vec<&str> = line.splitn(13, '\t').collect();
    if f.len() < 12 {
        return Err("expected 12 BED columns");
    }
    let plus = match f[5] {
        "+" => true,
        "-" => false,
        _ => return Err("strand must be + or -"),
    };
    let num = |s: &str| s.parse::<u32>().map_err(|_| "non-numeric BED field");
    let list = |s: &str| {
        s.split(',')
            .filter(|x| !x.is_empty())
            .map(num)
            .collect::<Result<Vec<_>, _>>()
    };
    let (start, end, n) = (num(f[1])? as u64, num(f[2])? as u64, num(f[9])? as usize);
    let (sizes, offsets) = (list(f[10])?, list(f[11])?);
    let block = |i: usize| {
        let s = start + offsets[i] as u64;
        (s, s + sizes[i] as u64)
    };
    if n == 0
        || sizes.len() != n
        || offsets.len() != n
        || offsets[0] != 0
        || sizes.contains(&0)
        || block(n - 1).1 != end
    {
        return Err("blocks do not tile chromStart-chromEnd");
    }
    let mut introns = Vec::with_capacity(n - 1);
    for i in 1..n {
        let (s, e) = (block(i - 1).1, block(i).0);
        if e < s {
            return Err("overlapping blocks");
        }
        if e > s {
            introns.push((s as u32, e as u32)); // adjacent blocks are one exon
        }
    }

    Ok((f[0], f[3], plus, start as u32, end as u32, introns))
}

/// Exons from a read's bounds and introns.
fn exons(start: u32, introns: &[Iv], end: u32) -> Vec<Iv> {
    let mut b = vec![start];
    b.extend(introns.iter().flat_map(|&(s, e)| [s, e]));
    b.push(end);
    b.chunks(2).map(|c| (c[0], c[1])).collect()
}

/// Junction correction of one site class (intron starts or ends) on one strand: site -> target.
/// A site moves to an annotated site, or one with >= 5x its reads, within the wobble; a site
/// holding >= 20% of the reads in its window never moves (NAGNAG), nor does an annotated one.
fn moves(support: &HashMap<u32, u32>, annotated: &[u32], wobble: u32) -> HashMap<u32, u32> {
    let mut sites: Vec<(u32, u32)> = support.iter().map(|(&p, &n)| (p, n)).collect();
    sites.sort_unstable();
    let ann = |p: u32| annotated.binary_search(&p).is_ok();
    let mut mv = HashMap::new();
    for &(p, n) in &sites {
        if ann(p) {
            continue;
        }
        let (lo, hi) = (p.saturating_sub(wobble), p.saturating_add(wobble));
        let near =
            &sites[sites.partition_point(|s| s.0 < lo)..sites.partition_point(|s| s.0 <= hi)];
        if 5 * n >= near.iter().map(|s| s.1).sum::<u32>() {
            continue;
        }
        let refs = &annotated
            [annotated.partition_point(|&a| a < lo)..annotated.partition_point(|&a| a <= hi)];
        // target: annotated, then highest support, then nearest, then lowest coordinate
        let target = near
            .iter()
            .copied()
            .filter(|&(q, m)| q != p && (ann(q) || m >= 5 * n))
            .chain(
                refs.iter()
                    .map(|&q| (q, support.get(&q).copied().unwrap_or(0))),
            )
            .min_by_key(|&(q, m)| (!ann(q), Reverse(m), q.abs_diff(p), q));
        if let Some((q, _)) = target {
            mv.insert(p, q);
        }
    }
    // a target that itself moved: follow it (support strictly grows, so this ends)
    let keys: Vec<u32> = mv.keys().copied().collect();
    for p in keys {
        let mut q = mv[&p];
        while let Some(&r) = mv.get(&q) {
            q = r;
        }
        mv.insert(p, q);
    }

    mv
}

/// For each query (sorted by start), the smallest tag of the intervals (sorted by start)
/// containing it; a sweep that keeps the open intervals keyed by end.
fn contained(ivs: &[(u32, u32, usize)], queries: &[Iv]) -> Vec<Option<usize>> {
    let mut open: BTreeMap<u32, usize> = BTreeMap::new();
    let mut j = 0;
    queries
        .iter()
        .map(|&(s, e)| {
            while j < ivs.len() && ivs[j].0 <= s {
                let tag = open.entry(ivs[j].1).or_insert(usize::MAX);
                *tag = (*tag).min(ivs[j].2);
                j += 1;
            }
            while open.first_key_value().is_some_and(|(&end, _)| end < s) {
                open.pop_first();
            }
            open.range(e..).map(|(_, &tag)| tag).min()
        })
        .collect()
}

/// Clusters of reads (sorted by start) linked by >= 1 bp of overlap.
fn clusters(reads: &[Read], idx: &[usize]) -> Vec<Vec<usize>> {
    let mut out: Vec<Vec<usize>> = Vec::new();
    let mut end = 0;
    for &i in idx {
        match out.last_mut() {
            Some(c) if reads[i].start < end => c.push(i),
            _ => out.push(vec![i]),
        }
        end = end.max(reads[i].end);
    }

    out
}

/// Step 8 for a read set (a chain's own reads or a mono cluster): 3′ end = mode of
/// tail-supported 3′ ends (all if none has a tail, ties to the most 3′, median if none repeats), 5′ end = most
/// upstream after trimming the outer 10%, representative = read closest to both ends
/// (ties: tail, higher IY, no FG, name), molecules = distinct (start, end).
fn model<'a>(reads: &[Read], chain: &'a [Iv], full: Vec<usize>) -> Model<'a> {
    let plus = reads[full[0]].plus;
    let upstream_first = |a: &u32, b: &u32| if plus { a.cmp(b) } else { b.cmp(a) };
    let mut f: Vec<u32> = full.iter().map(|&i| reads[i].five()).collect();
    f.sort_unstable_by(upstream_first);
    let five = f[f.len() / 10];
    let tails: Vec<u32> = full
        .iter()
        .filter(|&&i| reads[i].tail())
        .map(|&i| reads[i].three())
        .collect();
    let mut t = if tails.is_empty() {
        full.iter().map(|&i| reads[i].three()).collect()
    } else {
        tails
    };
    t.sort_unstable_by(upstream_first);
    let (mut three, mut best, mut run) = (t[0], 0, 0);
    for (j, &x) in t.iter().enumerate() {
        run = if j > 0 && t[j - 1] == x { run + 1 } else { 1 };
        if run >= best {
            (three, best) = (x, run);
        }
    }
    if best == 1 {
        three = t[t.len() / 2]; // INFO: no repeated end: the median, not the most extreme read
    }
    let rep = *full
        .iter()
        .min_by_key(|&&i| {
            let r = &reads[i];
            let d = r.five().abs_diff(five) as u64 + r.three().abs_diff(three) as u64;
            (d, !r.tail(), Reverse(r.iy), r.fg, r.name)
        })
        .expect("non-empty read set");
    let mut mol: Vec<Iv> = full
        .iter()
        .map(|&i| (reads[i].start, reads[i].end))
        .collect();
    mol.sort_unstable();
    mol.dedup();

    Model {
        chain,
        mol: mol.len(),
        five,
        three,
        rep,
        full,
        partial: Vec::new(),
        known: false,
        protected: false,
        reason: None,
    }
}

/// A 5′ truncation is kept separate when >= 10 of its 5′ ends fall within a 50-nt window, or
/// a reference shares its first intron with a 5′ end within 50 nt of its consensus 5′ end.
fn tss(reads: &[Read], m: &Model, rf: &Ref, plus: bool) -> bool {
    let mut f: Vec<u32> = m.full.iter().map(|&i| reads[i].five()).collect();
    f.sort_unstable();
    let first = if plus {
        m.chain[0]
    } else {
        m.chain[m.chain.len() - 1]
    };
    f.windows(10).any(|w| w[9] - w[0] < 50)
        || rf
            .first
            .get(&first)
            .is_some_and(|v| v.iter().any(|&x| x.abs_diff(m.five) <= 50))
}

/// The representative's first 12 columns: name without `#CN<n>`/`#SG` plus `#CN<n_reads>`
/// (and `#SG` for one read), thick = bounds, blocks rebuilt only if step 2 moved a junction.
fn bed_line(r: &Read, n: usize) -> String {
    let f: Vec<&str> = r.line.split('\t').take(12).collect();
    let mut name: Vec<&str> = f[3]
        .split('#')
        .enumerate()
        .filter(|&(j, t)| {
            let cn = t
                .strip_prefix("CN")
                .is_some_and(|d| d.bytes().all(|b| b.is_ascii_digit()));
            j == 0 || !(cn || t == "SG")
        })
        .map(|(_, t)| t)
        .collect();
    let cn = format!("CN{n}");
    name.push(&cn);
    if n == 1 {
        name.push("SG");
    }
    let (count, sizes, starts) = if r.corrected {
        let ex = exons(r.start, &r.introns, r.end);
        let join = |v: Vec<u32>| v.iter().map(u32::to_string).collect::<Vec<_>>().join(",");
        (
            ex.len().to_string(),
            join(ex.iter().map(|e| e.1 - e.0).collect()),
            join(ex.iter().map(|e| e.0 - r.start).collect()),
        )
    } else {
        (f[9].to_string(), f[10].to_string(), f[11].to_string())
    };

    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{count}\t{sizes}\t{starts}",
        f[0],
        f[1],
        f[2],
        name.join("#"),
        f[4],
        f[5],
        r.start,
        r.end,
        f[8]
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;

    /// A chr1 BED12 line with iso-segment style tags.
    fn bed(name: &str, strand: char, exons: &[Iv], tags: &str) -> String {
        let (s, e) = (exons[0].0, exons[exons.len() - 1].1);
        let join = |v: Vec<u32>| v.iter().map(u32::to_string).collect::<Vec<_>>().join(",");
        format!(
            "chr1\t{s}\t{e}\t{name}__FC0#TC0#{tags}\t0\t{strand}\t{s}\t{e}\t0\t{}\t{}\t{}\n",
            exons.len(),
            join(exons.iter().map(|x| x.1 - x.0).collect()),
            join(exons.iter().map(|x| x.0 - s).collect()),
        )
    }

    fn call(reads: &str, reference: &str, flags: &[&str]) -> Result<Vec<String>, Box<dyn Error>> {
        let cli = [
            "chain", "--bed", "a.bed", "--ref", "ref.bed", "--prefix", "t",
        ];
        let args = ChainArgs::parse_from(cli.iter().chain(flags));
        let mut out = vec![Vec::new(); 5];
        chain(
            &[("a.bed".into(), reads)],
            ("ref.bed", reference),
            &args,
            &mut out,
        )?;
        Ok(out
            .into_iter()
            .map(|o| String::from_utf8(o).unwrap())
            .collect())
    }

    #[test]
    fn chain_calls_models() {
        let (t, no) = ("PA30#PR30#IY990", "PA0#PR0#IY990"); // tail / no tail
        let (e2, e3, e4) = ((1200, 1300), (1400, 1500), (1600, 1800));
        let mut r = String::new();
        // known chain A (reference T1): a1 is a0's PCR twin, a10/a11 carry a donor 2 bp off,
        // a5 ends 5 bp past the representative
        for i in 0..12 {
            let s = if i == 1 { 995 } else { 995 + i };
            let e2 = if i >= 10 { (1200, 1302) } else { e2 };
            let e4 = if i == 5 { (1600, 1805) } else { e4 };
            r += &bed(&format!("a{i}"), '+', &[(s, 1100), e2, e3, e4], t);
        }
        r = r.replacen('\n', "\textra\n", 1); // a 13th column is dropped on output
                                              // 5′ cuts: b* are absorbed, c* form a TSS cluster and stay separate
        for i in 0..3 {
            r += &bed(&format!("b{i}"), '+', &[(1240 + 10 * i, 1300), e3, e4], t);
        }
        for i in 0..10 {
            r += &bed(&format!("c{i}"), '+', &[(1420 + i, 1500), e4], t);
        }
        let k0 = [(1000, 1100), (1402, 1500), e4]; // acceptor 2 bp off reference T2, one read
        r += &bed("k0", '+', &k0, &format!("{t}#SG"));
        // novel chains: n* hold 3 of 30 locus reads, d* are PCR twins, e* lack tails
        for i in 0..3 {
            r += &bed(&format!("n{i}"), '+', &[(1000 + i, 1100), e2, e4], t);
        }
        for i in 0..2 {
            r += &bed(&format!("d{i}"), '+', &[(1000, 1100), (1350, 1500), e4], t);
        }
        for i in 0..2 {
            r += &bed(&format!("e{i}"), '+', &[(1000 + i, 1100), (1450, 1700)], no);
        }
        r += &bed("m1", '+', &[(1650, 1803)], t); // inside A's last exon as its reads reach it
        r += &bed("m2", '+', &[(1120, 1180)], t); // inside A's first intron
        for i in 0..10 {
            r += &bed(&format!("i{i}"), '+', &[(5000 + i, 5500)], "PA0#PR10"); // genomic A's
        }
        r += &bed("km", '+', &[(7100, 7900)], t); // inside single-exon reference T3
        let (x1, x2) = ((20000, 20100), (20200, 20300));
        for i in 0..3 {
            r += &bed(&format!("f{i}"), '-', &[x1, x2, (20400, 20590 + 5 * i)], t);
        }
        r += &bed("t0", '-', &[x1, (20200, 20280)], t); // 5′ cut on minus (5′ = end)
        r += &bed("u0", '-', &[(20250, 20300), (20400, 20600)], t); // 3′ cut with a tail
        let reference = [
            bed("T1", '+', &[(1000, 1100), e2, e3, e4], ""),
            bed("T2", '+', &[(1000, 1100), e3, e4], ""),
            bed("T3", '+', &[(7000, 8000)], ""),
            bed("T4", '+', &[(1, 2)], "").replace("chr1", "chr2"),
        ]
        .concat();

        // balanced (tail switch on: 39/51 reads have PA >= 20) with an explicit read fraction
        let o = call(&r, &reference, &["--min-read-fraction", "0.8"]).unwrap();
        assert_eq!(
            o[4],
            "reads_in\t51\nreads_multi\t38\nreads_mono\t13\nsites_corrected\t2\n\
             reads_corrected\t3\nreads_correction_skipped\t0\nchains\t10\nchains_absorbed\t2\n\
             models_known\t3\nmodels_novel\t2\nmodels_mono\t0\nexcluded_support\t3\n\
             excluded_fraction\t3\nexcluded_tail\t2\nexcluded_intronic\t1\n\
             excluded_intraprimed\t10\nmodels_out\t5\ntail_rate\t0.7647\n"
        );
        let fate: HashMap<&str, &str> = o[2]
            .lines()
            .skip(1)
            .map(|l| l.split_once('\t').unwrap())
            .map(|(read, m)| (read.split("__").next().unwrap(), m))
            .collect();
        assert_eq!((o[2].lines().count(), fate.len()), (52, 51)); // every read exactly once
        let a = fate["a0"];
        for (read, want) in [
            ("a1", a),
            ("a11", a),
            ("b0", a),
            ("m1", a),
            ("t0", fate["f0"]),
            ("d0", "excluded:support"),
            ("u0", "excluded:support"),
            ("e0", "excluded:tail"),
            ("n0", "excluded:fraction"),
            ("m2", "excluded:intronic"),
            ("i0", "excluded:intraprimed"),
        ] {
            assert_eq!(fate[read], want, "{read}");
        }
        let row = |file: &str, col: usize, read: &str| -> String {
            let hit = |l: &&str| {
                l.split('\t')
                    .nth(col)
                    .unwrap()
                    .starts_with(&format!("{read}__"))
            };
            file.lines().find(hit).unwrap().into()
        };
        assert_eq!(
            row(&o[0], 3, "a0"),
            "chr1\t995\t1800\ta0__FC0#TC0#PA30#PR30#IY990#CN16\t0\t+\t995\t1800\t0\t4\t105,100,100,200\t0,205,405,605"
        );
        assert_eq!(
            row(&o[0], 3, "k0"),
            "chr1\t1000\t1800\tk0__FC0#TC0#PA30#PR30#IY990#CN1#SG\t0\t+\t1000\t1800\t0\t3\t100,100,200\t0,400,600"
        );
        let support = |read: &str| {
            row(&o[1], 1, read)
                .split('\t')
                .skip(2)
                .collect::<Vec<_>>()
                .join(" ")
        };
        assert_eq!(
            support("a0"),
            "+ 1100-1200,1300-1400,1500-1600 known 16 15 4 16 995 1800 ."
        );
        assert_eq!(
            support("c1"),
            "+ 1500-1600 novel 10 10 0 10 1421 1800 protected"
        );
        assert_eq!(
            support("k0"),
            "+ 1100-1400,1500-1600 known 1 1 0 1 1000 1800 corrected"
        );
        assert_eq!(support("km"), "+ . known 1 1 0 1 7100 7900 .");
        assert_eq!(
            support("f2"),
            "- 20100-20200,20300-20400 novel 4 4 1 4 20600 20000 ."
        );
        let excluded: Vec<(&str, &str)> = o[3]
            .lines()
            .map(|l| {
                (
                    l.split('\t').nth(3).unwrap(),
                    l.split('\t').nth(12).unwrap(),
                )
            })
            .map(|(name, why)| (name.split("__").next().unwrap(), why))
            .collect();
        assert_eq!(
            excluded,
            [
                ("e0", "tail"),
                ("d0", "support"),
                ("n0", "fraction"),
                ("m2", "intronic"),
                ("i1", "intraprimed"),
                ("u0", "support")
            ]
        );

        // sensitive keeps every distinct chain and the tail-less mono cluster; intronic stays out
        let s = call(&r, &reference, &["--preset", "sensitive"]).unwrap();
        let tail = s[4].lines().skip(8).collect::<Vec<_>>().join(" ");
        assert_eq!(
            tail,
            "models_known\t3 models_novel\t6 models_mono\t1 excluded_support\t0 excluded_fraction\t0 \
             excluded_tail\t0 excluded_intronic\t1 excluded_intraprimed\t0 models_out\t10 tail_rate\t0.7647"
        );
    }

    #[test]
    fn chain_edges() {
        // empty input still writes all five outputs
        let o = call("", "", &[]).unwrap();
        assert_eq!(
            (o[0].as_str(), o[2].as_str(), o[3].as_str()),
            ("", "read\tmodel\n", "")
        );
        assert_eq!(o[1].lines().count(), 1);
        assert!(o[4].ends_with("models_out\t0\ntail_rate\t0.0000\n"));
        // a strand other than + or - is a hard error naming file and line
        let bad = bed("x", '+', &[(1, 10)], "") + &bed("y", '.', &[(1, 10)], "");
        let e = call(&bad, "", &[]).unwrap_err().to_string();
        assert_eq!(e, "a.bed:2: strand must be + or -");
        // no repeated 3' end: the median read, not the most extreme one, sets the model end
        let ends: String = [598, 600, 603]
            .iter()
            .enumerate()
            .map(|(i, &end)| {
                bed(
                    &format!("m{i}"),
                    '+',
                    &[(100, 200), (300, end)],
                    "PA30#PR30#IY990",
                )
            })
            .collect();
        let o = call(&ends, "", &[]).unwrap();
        assert!(o[0].starts_with("chr1\t100\t600\tm1__"), "{}", o[0]);
    }
}
