// Copyright (c) 2026 Alejandro Gonzalez-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! `collapse chain`: transcript models from segmented reads by intron chain (ARK PLAN.md §5).
//!
//! Junction wobble is corrected per strand, reads are grouped by intron chain, truncated
//! subchains are absorbed into their best-supported parent, mono-exonic reads are absorbed,
//! excluded as intronic or clustered, and models pass support, tail and read-fraction
//! thresholds. A 3′ cut whose end fits a reference transcript end (`reaches`) and that the
//! parent's own reads do not outnumber is not a truncation. A 5′ cut stays apart with a cluster
//! of its own read starts or at the annotated 5′ end of a reference; a bare CDS start (CDS-only
//! reference) is no TSS and keeps nothing apart (`tss`). With `--subchain-merge`, a kept model whose intron chain is a contiguous run of a
//! kept longer model's chain then joins it (step 11). Every model is a real read (its
//! representative) tagged `#CN<n_reads>`; where reads hold a reference CDS that the consensus ends
//! cut, the representative is one of them (`cds_rep`). Each `--bed` file is one sample:
//! `samples.tsv` has every model's support per sample (with the model's BED name, the start of its
//! name in every later table), `introns.tsv` every distinct intron once.

use std::cmp::Reverse;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::cli::{ChainArgs, Preset, SubchainMerge};
use crate::utils::Dsu;

/// An intron with this many reads never moves in junction correction.
const PROTECT_READS: u32 = 20;

/// A model end this far inside a reference CDS end still reaches it (alignment end jitter).
const CDS_SLACK: u32 = 5;

/// Half-open genomic interval (intron or exon).
type Iv = (u32, u32);

/// Chrom, name, strand (true = +), start, end and introns of a BED12 line.
type Bed<'a> = (&'a str, &'a str, bool, u32, u32, Vec<Iv>);

const OUTPUTS: [&str; 8] = [
    "models.bed",
    "support.tsv",
    "members.tsv",
    "excluded.bed",
    "counts.tsv",
    "rejected.bed",
    "samples.tsv",
    "introns.tsv",
];

/// A segmented read: its BED12 line and the iso-segment tags the caller uses.
struct Read<'a> {
    line: &'a str,
    /// index of the `--bed` file it came from (one file is one sample)
    src: usize,
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
    /// cluster2 singleton (`#SG`)
    sg: bool,
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
    /// annotated introns, sorted
    introns: Vec<Iv>,
    /// first intron in transcript orientation -> 5′ ends, each with whether it is a bare CDS start
    first: HashMap<Iv, Vec<(u32, bool)>>,
    /// last intron in transcript orientation -> 3′ ends, each with whether it is a bare CDS end
    last: HashMap<Iv, Vec<(u32, bool)>>,
    /// single-exon transcripts, sorted
    mono: Vec<(u32, u32, usize)>,
    /// coding rows with >= 1 intron inside the CDS: thick start, thick end, those introns
    cds: Vec<(u32, u32, Vec<Iv>)>,
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
    /// its 3′ end fits the 3′ end of a reference transcript with the same last intron
    ref_three: bool,
    /// kept apart from a longer container by `--subchain-merge ends`
    sub_kept: bool,
    /// kept by `--cds-pool`
    cds_pooled: bool,
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

/// Runs `collapse chain`, writing
/// `<prefix>.{models.bed,support.tsv,members.tsv,excluded.bed,counts.tsv,rejected.bed,samples.tsv,introns.tsv}`.
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
    let [models_w, support_w, members_w, excluded_w, counts_w, rejected_w, samples_w, introns_w] =
        out
    else {
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
    for (src, (path, text)) in beds.iter().enumerate() {
        for (n, line) in text.lines().enumerate() {
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let (chrom, name, plus, start, end, introns) =
                parse_bed(line).map_err(|e| format!("{path}:{}: {e}", n + 1))?;
            let (mut pa, mut pr, mut tc, mut iy, mut fg, mut sg) = (0, 0, 0, 0, false, false);
            for tag in name.split_once("__").map_or("", |t| t.1).split('#') {
                let v = tag.get(2..).and_then(|v| v.parse().ok()).unwrap_or(0);
                match tag.get(..2) {
                    Some("PA") => pa = v,
                    Some("PR") => pr = v,
                    Some("TC") => tc = v,
                    Some("IY") => iy = v,
                    Some("FG") => fg = true,
                    Some("SG") => sg = true,
                    _ => {}
                }
            }
            reads.push(Read {
                line,
                src,
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
                sg,
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
        // a transcript end that is also a thick bound has no UTR: a CDS-only row (TOGA2 projection)
        let thick = |i: usize| line.split('\t').nth(i).and_then(|v| v.parse::<u32>().ok());
        let (ts, te) = (thick(6).unwrap_or(0), thick(7).unwrap_or(0));
        let (left, right) = ((start, ts < te && ts == start), (end, ts < te && te == end));
        let r = refs.entry((chrom, plus)).or_default();
        let coding: Vec<Iv> = introns
            .iter()
            .filter(|iv| ts < te && iv.0 >= ts && iv.1 <= te)
            .copied()
            .collect();
        if !coding.is_empty() {
            r.cds.push((ts, te, coding));
        }
        match (introns.first(), introns.last()) {
            (Some(&a), Some(&b)) => {
                r.starts.extend(introns.iter().map(|i| i.0));
                r.ends.extend(introns.iter().map(|i| i.1));
                let ((first, five), (last, three)) = if plus {
                    ((a, left), (b, right))
                } else {
                    ((b, right), (a, left))
                };
                r.first.entry(first).or_default().push(five);
                r.last.entry(last).or_default().push(three);
                r.introns.extend(&introns);
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
        r.introns.sort_unstable();
        r.introns.dedup();
        r.mono.sort_unstable();
    }
    let none = Ref::default();

    // cluster2 input marks its singletons; without them every read is one molecule
    let clustered = reads.iter().any(|r| r.sg);

    // tail switch (step 7): on when enough of the input reads carry a polyA tail
    let tail_rate = reads.iter().filter(|r| r.tail()).count() as f64 / reads.len().max(1) as f64;
    let tail_on = tail_rate >= tail_min;

    let mut buckets: BTreeMap<(&str, bool), Vec<usize>> = BTreeMap::new();
    for (i, r) in reads.iter().enumerate() {
        buckets.entry((r.chrom, r.plus)).or_default().push(i);
    }

    // 2. junction correction per strand, by whole intron: a rare intron moves onto a nearby one with
    // 5x its support, so a minor alternative site (NAGNAG) is judged against its own intron's
    // neighbours, not against every read sharing the other site, and a corrected chain only ever
    // carries introns that some read or the reference has
    let (mut sites_moved, mut reads_moved, mut reads_kept) = (0, 0, 0);
    for (key, idx) in &buckets {
        let rf = refs.get(key).unwrap_or(&none);
        let (mut support, mut sites) = (HashMap::new(), (HashMap::new(), HashMap::new()));
        for &i in idx {
            for &iv in &reads[i].introns {
                *support.entry(iv).or_insert(0) += 1;
                *sites.0.entry(iv.0).or_insert(0) += 1;
                *sites.1.entry(iv.1).or_insert(0) += 1;
            }
        }
        let mv = moves(&support, &sites, rf, args.junction_wobble);
        sites_moved += mv.len();
        if mv.is_empty() {
            continue;
        }
        for &i in idx {
            let r = &mut reads[i];
            let new: Vec<Iv> = r
                .introns
                .iter()
                .map(|iv| *mv.get(iv).unwrap_or(iv))
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

    let (mut n_chains, mut n_absorbed, mut n_sub, mut n_sub_kept) = (0, 0, 0, 0);
    let mut all: Vec<Model> = Vec::new();
    let mut introns: Vec<(&str, Iv, bool, usize, bool)> = Vec::new();
    for (key, idx) in &buckets {
        let (rf, plus) = (refs.get(key).unwrap_or(&none), key.1);

        // every distinct intron of the corrected reads, once, with its reads (introns.tsv)
        let mut seen: HashMap<Iv, usize> = HashMap::new();
        for &i in idx {
            for &iv in &reads[i].introns {
                *seen.entry(iv).or_insert(0) += 1;
            }
        }
        introns.extend(
            seen.into_iter()
                .map(|(iv, n)| (key.0, iv, plus, n, rf.introns.binary_search(&iv).is_ok())),
        );

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
            .map(|(chain, full)| {
                let m = model(&reads, chain, full);
                let last = if plus {
                    chain[chain.len() - 1]
                } else {
                    chain[0]
                };
                let ends = rf.last.get(&last);
                Model {
                    known: rf.chains.contains(chain),
                    ref_three: ends.is_some_and(|v| v.iter().any(|&x| reaches(m.three, x, !plus))),
                    rep: cds_rep(&reads, &m, rf, plus),
                    ..m
                }
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
                // a 3′ cut stays apart with a tail-supported read, or (FLNC reads carry no tail tags)
                // when a reference ends where it ends and B's own reads do not outnumber it
                if three_inner
                    && (tail || (ms[a].ref_three && ms[a].full.len() >= ms[b].full.len()))
                {
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

        // 7b. --junction-support: a novel chain resting on one singleton molecule stays only if
        // every novel intron is also carried by a read outside it (singleton audit: one-read
        // novel junctions are 18% flagged and 7% non-canonical, independently seen ones are not)
        if args.junction_support {
            let mut carried: HashMap<Iv, usize> = HashMap::new();
            for &i in idx {
                for &iv in &reads[i].introns {
                    *carried.entry(iv).or_insert(0) += 1;
                }
            }
            for m in ms.iter_mut() {
                let own: Vec<&Read> = m
                    .full
                    .iter()
                    .chain(&m.partial)
                    .map(|&i| &reads[i])
                    .collect();
                if m.reason.is_some()
                    || m.known
                    || m.mol != 1
                    || (clustered && !own.iter().all(|r| r.sg))
                {
                    continue;
                }
                let lone = m.chain.iter().any(|iv| {
                    rf.introns.binary_search(iv).is_err()
                        && carried[iv] <= own.iter().filter(|r| r.introns.contains(iv)).count()
                });
                if lone {
                    m.reason = Some("junction");
                }
            }
        }

        // 7c. --cds-pool: chains short of support whose representatives hold the same reference
        // CDS intact, held by no kept model, are one coding transcript with different UTR
        // structures: their molecules count together and the best chain is kept
        if args.cds_pool && !rf.cds.is_empty() {
            let mut by_first: HashMap<Iv, Vec<usize>> = HashMap::new();
            for (c, (_, _, coding)) in rf.cds.iter().enumerate() {
                by_first.entry(coding[0]).or_default().push(c);
            }
            let holds = |r: &Read| -> Vec<usize> {
                let ex = exons(r.start, &r.introns, r.end);
                let mut out = Vec::new();
                for iv in &r.introns {
                    for &c in by_first.get(iv).map_or(&[][..], |v| v) {
                        let (ts, te, coding) = &rf.cds[c];
                        if r.introns
                            .iter()
                            .filter(|x| x.0 >= *ts && x.1 <= *te)
                            .eq(coding)
                            && ex.iter().any(|e| e.0 <= *ts && *ts < e.1)
                            && ex.iter().any(|e| e.0 < *te && *te <= e.1)
                        {
                            out.push(c);
                        }
                    }
                }
                out
            };
            let mut held: HashSet<usize> = HashSet::new();
            let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
            for (j, m) in ms.iter().enumerate().filter(|(_, m)| !m.chain.is_empty()) {
                match m.reason {
                    None => held.extend(holds(&reads[m.rep])),
                    Some("support") => {
                        for c in holds(&reads[m.rep]) {
                            groups.entry(c).or_default().push(j);
                        }
                    }
                    _ => {}
                }
            }
            for (c, js) in groups {
                // a chain restored for another CDS of the group already carries this one
                if held.contains(&c) || js.iter().any(|&j| ms[j].reason.is_none()) {
                    continue;
                }
                let mut mol: Vec<Iv> = js
                    .iter()
                    .flat_map(|&j| &ms[j].full)
                    .map(|&i| (reads[i].start, reads[i].end))
                    .collect();
                mol.sort_unstable();
                mol.dedup();
                if mol.len() < novel {
                    continue;
                }
                let unannotated = |j: usize| {
                    let novel = |iv: &&Iv| rf.introns.binary_search(iv).is_err();
                    ms[j].chain.iter().filter(novel).count()
                };
                let best = *js
                    .iter()
                    .max_by_key(|&&j| (Reverse(unannotated(j)), ms[j].n_reads(), Reverse(j)))
                    .expect("a group has a chain");
                // it stands for a reference CDS: exempt from the read fraction like a known chain
                ms[best].reason = None;
                ms[best].known = true;
                ms[best].cds_pooled = true;
                held.insert(c);
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
            // models with equal support share one decision: the reads of the better-supported ones
            let (mut cum, mut before, mut tie, mut dropped) = (0, 0, usize::MAX, false);
            for &m in &locus {
                let n = ms[m].n_reads();
                if n != tie {
                    (before, tie) = (cum, n);
                }
                if !ms[m].known && before as f64 >= fraction * total {
                    ms[m].reason = Some("fraction");
                    dropped = true;
                }
                cum += n;
            }
            // --sample-read-fraction: in a pool the well-expressed isoforms of every sample add up
            // and push a model that one sample expresses below the fraction. The same rule over
            // the locus reads of each sample (one `--bed` file) restores what a sample keeps alone
            if !(args.sample_read_fraction && dropped && beds.len() > 1) {
                continue;
            }
            let per: Vec<Vec<usize>> = locus
                .iter()
                .map(|&m| {
                    let mut n = vec![0; beds.len()];
                    for &i in ms[m].full.iter().chain(&ms[m].partial) {
                        n[reads[i].src] += 1;
                    }
                    n
                })
                .collect();
            for src in 0..beds.len() {
                let total = per.iter().map(|n| n[src]).sum::<usize>() as f64;
                let mut order: Vec<usize> = (0..locus.len()).filter(|&x| per[x][src] > 0).collect();
                order.sort_by_key(|&x| Reverse(per[x][src]));
                let (mut cum, mut before, mut tie) = (0, 0, usize::MAX);
                for x in order {
                    let n = per[x][src];
                    if n != tie {
                        (before, tie) = (cum, n);
                    }
                    if (before as f64) < fraction * total {
                        ms[locus[x]].reason = None; // only this step's drops are in the locus
                    }
                    cum += n;
                }
            }
        }

        // 11. --subchain-merge: a kept chain that is a contiguous run of a kept longer one joins it
        if args.subchain_merge != SubchainMerge::Off {
            let (merged, kept) = subchain_merge(&reads, &mut ms, args);
            (n_sub, n_sub_kept) = (n_sub + merged, n_sub_kept + kept);
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
        let flags: Vec<&str> = [
            (r.corrected, "corrected"),
            (m.protected, "protected"),
            (m.sub_kept, "subchain_kept"),
            (m.cds_pooled, "cds_pooled"),
        ]
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
        if why == "junction" {
            for &i in m.full.iter().chain(&m.partial) {
                writeln!(rejected_w, "{}", bed_line(&reads[i], 1))?;
            }
        }
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

    // support of every model per sample (one `--bed` file each, named by its file stem)
    let label = |p: &str| {
        let name = p.rsplit('/').next().unwrap_or(p);
        name.strip_suffix(".bed").unwrap_or(name).to_string()
    };
    let labels: Vec<String> = beds.iter().map(|(p, _)| label(p)).collect();
    writeln!(
        samples_w,
        "model\tname\tsample\tn_reads\tn_full\tn_full_molecules"
    )?;
    for (n, m) in kept.iter().enumerate() {
        let line = bed_line(&reads[m.rep], m.n_reads());
        let name = line.split('\t').nth(3).unwrap_or_default();
        let mut per: BTreeMap<usize, (usize, Vec<Iv>)> = BTreeMap::new();
        for &i in &m.full {
            per.entry(reads[i].src)
                .or_default()
                .1
                .push((reads[i].start, reads[i].end));
        }
        for &i in &m.partial {
            per.entry(reads[i].src).or_default().0 += 1;
        }
        for (src, (partial, mut full)) in per {
            let n_full = full.len();
            full.sort_unstable();
            full.dedup();
            writeln!(
                samples_w,
                "{prefix}.{}\t{name}\t{}\t{}\t{n_full}\t{}",
                n + 1,
                labels[src],
                partial + n_full,
                full.len()
            )?;
        }
    }
    introns.sort_unstable();
    writeln!(introns_w, "chrom\tstart\tend\tstrand\tn_reads\tannotated")?;
    for (chrom, (s, e), plus, n, ann) in introns {
        let strand = if plus { '+' } else { '-' };
        writeln!(introns_w, "{chrom}\t{s}\t{e}\t{strand}\t{n}\t{}", ann as u8)?;
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
        ("chains_subchain_merged", n_sub),
        ("chains_subchain_kept", n_sub_kept),
        ("models_known", models("known")),
        ("models_novel", models("novel")),
        ("models_mono", models("mono")),
        ("excluded_support", excluded("support")),
        ("excluded_fraction", excluded("fraction")),
        ("excluded_tail", excluded("tail")),
        ("excluded_intronic", excluded("intronic")),
        ("excluded_intraprimed", excluded("intraprimed")),
        ("excluded_junction", excluded("junction")),
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
fn moves(
    support: &HashMap<Iv, u32>,
    (starts, ends): &(HashMap<u32, u32>, HashMap<u32, u32>),
    rf: &Ref,
    wobble: u32,
) -> HashMap<Iv, Iv> {
    let site = |v: &[u32], p: u32| v.binary_search(&p).is_ok();
    let ann = |iv: Iv| rf.introns.binary_search(&iv).is_ok();
    // an annotated intron weighs at least as much as the reads at its weaker site, so a read a few
    // bp off an annotated intron that no read spans still snaps onto it
    let weight = |iv: Iv| {
        let own = support.get(&iv).copied().unwrap_or(0);
        if !ann(iv) {
            return own;
        }
        let at = |m: &HashMap<u32, u32>, p| m.get(&p).copied().unwrap_or(0);
        own.max(at(starts, iv.0).min(at(ends, iv.1)))
    };
    let mut introns: Vec<(Iv, u32)> = support.keys().map(|&iv| (iv, weight(iv))).collect();
    introns.extend(
        rf.introns
            .iter()
            .filter(|iv| !support.contains_key(iv))
            .map(|&iv| (iv, weight(iv))),
    );
    introns.sort_unstable();
    let near = |a: Iv, b: Iv| a.0.abs_diff(b.0) <= wobble && a.1.abs_diff(b.1) <= wobble;
    let mut mv = HashMap::new();
    for &(iv, _) in &introns {
        let Some(&n) = support.get(&iv) else { continue };
        // annotated sites never move
        let (fs, fe) = (site(&rf.starts, iv.0), site(&rf.ends, iv.1));
        if ann(iv) || fs && fe {
            continue;
        }
        let lo = iv.0.saturating_sub(wobble);
        let window = introns[introns.partition_point(|x| x.0 .0 < lo)..]
            .iter()
            .take_while(|x| x.0 .0 <= iv.0 + wobble)
            .filter(|x| near(x.0, iv));
        // an intron with >= 20% of its window's reads, or with PROTECT_READS reads, is protected: at
        // depth a real minor site (SIRV604, 164 reads) can sit under 20% of a site shared by
        // several isoforms, while only 1 of 85 SIRV shift artifacts reaches 20 reads
        if n >= PROTECT_READS || 5 * n >= window.clone().map(|x| x.1).sum::<u32>() {
            continue;
        }
        // target: annotated, then highest support, then nearest, then lowest coordinate; it must
        // keep any annotated site of the moving intron
        let target = window
            .copied()
            .filter(|&(q, m)| q != iv && (ann(q) || m >= 5 * n))
            .filter(|&(q, _)| (!fs || q.0 == iv.0) && (!fe || q.1 == iv.1))
            .min_by_key(|&(q, m)| {
                (
                    !ann(q),
                    Reverse(m),
                    q.0.abs_diff(iv.0) + q.1.abs_diff(iv.1),
                    q,
                )
            });
        if let Some((q, _)) = target {
            mv.insert(iv, q);
        }
    }
    // a target that itself moved: follow it (support strictly grows, so this ends)
    let keys: Vec<Iv> = mv.keys().copied().collect();
    for p in keys {
        let mut q = mv[&p];
        while let Some(&r) = mv.get(&q) {
            q = r;
        }
        mv.insert(p, q);
    }

    mv
}

/// Step 11 (`--subchain-merge`) on one (chrom, strand): a kept multi-exon model whose intron chain
/// is a contiguous proper sub-chain of a kept longer model's chain is merged into the
/// best-supported such container (most reads, then most introns; containers are decided first,
/// longest chains first). Its reads become the container's partial reads and its molecules are
/// added, so `#CN`, support.tsv and members.tsv follow. Plain merges every sub-chain, whatever
/// its ends. `ends` skips a container when the shorter model is known (unless
/// `--subchain-known-needs-ends`), is a 3′ cut that step 5 keeps apart for its reference end, or has
/// `--subchain-min-end-reads` reads ending past the container on one side (`distinct`); with no
/// other container it stays a model, flagged `subchain_kept`.
/// Returns (merged, kept by `ends`).
fn subchain_merge(reads: &[Read], ms: &mut Vec<Model>, args: &ChainArgs) -> (usize, usize) {
    let live = |m: &Model| m.reason.is_none() && !m.chain.is_empty();
    let mut by_intron: HashMap<Iv, Vec<usize>> = HashMap::new();
    for (c, m) in ms.iter().enumerate().filter(|(_, m)| live(m)) {
        for &iv in m.chain {
            by_intron.entry(iv).or_default().push(c);
        }
    }
    let mut order: Vec<usize> = (0..ms.len()).filter(|&c| live(&ms[c])).collect();
    order.sort_by_key(|&c| Reverse(ms[c].chain.len()));
    let ends = args.subchain_merge == SubchainMerge::Ends;
    let mut merged = vec![false; ms.len()];
    let (mut n_merged, mut n_kept) = (0, 0);
    for a in order {
        let (ac, k) = (ms[a].chain, ms[a].chain.len());
        let cands = ac
            .iter()
            .map(|iv| &by_intron[iv])
            .min_by_key(|v| v.len())
            .expect("multi-exon chain");
        let (mut best, mut exempt) = (None::<usize>, false);
        for &b in cands {
            let (bc, m) = (ms[b].chain, ms[b].chain.len());
            if m <= k || merged[b] {
                continue;
            }
            // introns of a chain are distinct and sorted, so A can only start at one place in B
            let Ok(p) = bc.binary_search(&ac[0]) else {
                continue;
            };
            if p + k > m || bc[p..p + k] != *ac {
                continue;
            }
            let known_stays = ms[a].known && !args.subchain_known_needs_ends;
            let three_cut = if reads[ms[a].rep].plus {
                p + k < m
            } else {
                p > 0
            };
            let end_stays = three_cut && ms[a].ref_three && ms[a].full.len() >= ms[b].full.len();
            if ends
                && (known_stays
                    || end_stays
                    || distinct(reads, &ms[a], &reads[ms[b].rep], bc, p, args))
            {
                exempt = true;
                continue;
            }
            let rank = |x: usize| (ms[x].n_reads(), ms[x].chain.len());
            if best.is_none_or(|x| rank(b) > rank(x)) {
                best = Some(b);
            }
        }
        match best {
            Some(b) => {
                merged[a] = true;
                n_merged += 1;
                let mut moved = std::mem::take(&mut ms[a].full);
                moved.append(&mut ms[a].partial);
                ms[b].mol += ms[a].mol;
                ms[b].partial.extend(moved);
            }
            None if exempt => {
                ms[a].sub_kept = true;
                n_kept += 1;
            }
            None => {}
        }
    }
    let mut gone = merged.into_iter();
    ms.retain(|_| !gone.next().unwrap_or(false));

    (n_merged, n_kept)
}

/// `--subchain-merge ends`: whether >= `--subchain-min-end-reads` reads of `a` end past its
/// container on one side. `a`'s introns are `bc[p..]` of the container, whose representative is
/// `r`. Where `a`'s terminal exon is an internal exon of the container, a read must reach into the
/// container's intron by more than `--subchain-intron-tolerance` (an alternative first or last
/// exon); where both are terminal, beyond the container's end by more than
/// `--subchain-end-tolerance` (an alternative TSS or TES). Ends inside the container are
/// truncations and never count.
fn distinct(reads: &[Read], a: &Model, r: &Read, bc: &[Iv], p: usize, args: &ChainArgs) -> bool {
    let k = a.chain.len();
    let (it, et) = (args.subchain_intron_tolerance, args.subchain_end_tolerance);
    let (lo, lt) = if p > 0 {
        (bc[p - 1].1, it)
    } else {
        (r.start, et)
    };
    let (hi, ht) = if p + k < bc.len() {
        (bc[p + k].0, it)
    } else {
        (r.end, et)
    };
    let members = || a.full.iter().chain(&a.partial).map(|&i| &reads[i]);
    let left = members()
        .filter(|x| x.start.saturating_add(lt) < lo)
        .count();
    let right = members().filter(|x| x.end > hi.saturating_add(ht)).count();

    left.max(right) >= args.subchain_min_end_reads
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
    let rep = closest(reads, &full, five, three);
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
        ref_three: false,
        sub_kept: false,
        cds_pooled: false,
        reason: None,
    }
}

/// The read closest to both ends (ties: tail, higher IY, no FG, name).
fn closest(reads: &[Read], idx: &[usize], five: u32, three: u32) -> usize {
    *idx.iter()
        .min_by_key(|&&i| {
            let r = &reads[i];
            let d = r.five().abs_diff(five) as u64 + r.three().abs_diff(three) as u64;
            (d, !r.tail(), Reverse(r.iy), r.fg, r.name)
        })
        .expect("non-empty read set")
}

/// The representative of a multi-exon chain where the reference has a CDS-only row (no UTR, as
/// TOGA2 projections) with the chain's first or last intron. The consensus ends can fall inside
/// that CDS (the 3′ mode at an internal A-run, the trimmed 5′ end) although reads hold its first
/// and last coding base. When the consensus representative misses one of them and >= 2 full-chain
/// reads hold both, the closest of those to the consensus ends is the representative. Only the
/// representative changes: the ends that steps 5 and 11 judge stay the consensus ends.
fn cds_rep(reads: &[Read], m: &Model, rf: &Ref, plus: bool) -> usize {
    let (a, b) = (m.chain[0], m.chain[m.chain.len() - 1]);
    let (first, last) = if plus { (a, b) } else { (b, a) };
    let bare = |ends: Option<&Vec<(u32, bool)>>| -> Vec<u32> {
        let ends = ends.into_iter().flatten();
        ends.filter(|x| x.1).map(|x| x.0).collect()
    };
    let (start, stop) = (bare(rf.first.get(&first)), bare(rf.last.get(&last)));
    // at or upstream of a CDS start, at or downstream of a CDS end (no slack: every coding base)
    let holds = |i: usize| {
        let (five, three) = (reads[i].five(), reads[i].three());
        let up = |x: &u32| if plus { five <= *x } else { five >= *x };
        let down = |x: &u32| if plus { three >= *x } else { three <= *x };
        (start.is_empty() || start.iter().any(up)) && (stop.is_empty() || stop.iter().any(down))
    };
    if (start.is_empty() && stop.is_empty()) || holds(m.rep) {
        return m.rep;
    }
    let good: Vec<usize> = m.full.iter().copied().filter(|&i| holds(i)).collect();
    if good.len() < 2 {
        return m.rep;
    }

    closest(reads, &good, m.five, m.three)
}

/// Whether a model's 3′ end fits a reference 3′ end `x`: within 50 nt of an annotated end, or,
/// where the reference has no UTR (`bare`: the end is its CDS stop), at or beyond it, so that the
/// model holds that end of the CDS. `low`: beyond means a lower coordinate.
fn reaches(end: u32, (x, bare): (u32, bool), low: bool) -> bool {
    if !bare {
        x.abs_diff(end) <= 50
    } else if low {
        end <= x.saturating_add(CDS_SLACK)
    } else {
        end.saturating_add(CDS_SLACK) >= x
    }
}

/// A 5′ truncation is kept separate when >= 10 of its 5′ ends fall within a 50-nt window, or
/// its consensus 5′ end is within 50 nt of the annotated 5′ end of a reference with the same
/// first intron. A bare CDS start (a reference row without a 5′ UTR) is a start codon, not a
/// TSS, and keeps nothing apart: on S2 that protection made 613 false pass models for 24 true.
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
            .is_some_and(|v| v.iter().any(|&(x, bare)| !bare && x.abs_diff(m.five) <= 50))
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
        call_beds(&[("a.bed", reads)], reference, flags)
    }

    /// `call` over several input files (one sample each).
    fn call_beds(
        beds: &[(&str, &str)],
        reference: &str,
        flags: &[&str],
    ) -> Result<Vec<String>, Box<dyn Error>> {
        let cli = [
            "chain", "--bed", "a.bed", "--ref", "ref.bed", "--prefix", "t",
        ];
        let args = ChainArgs::parse_from(cli.iter().chain(flags));
        let mut out = vec![Vec::new(); OUTPUTS.len()];
        let beds: Vec<(String, &str)> = beds.iter().map(|(p, t)| (p.to_string(), *t)).collect();
        chain(&beds, ("ref.bed", reference), &args, &mut out)?;
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
             chains_subchain_merged\t0\nchains_subchain_kept\t0\nmodels_known\t3\nmodels_novel\t2\nmodels_mono\t0\nexcluded_support\t3\n\
             excluded_fraction\t3\nexcluded_tail\t2\nexcluded_intronic\t1\n\
             excluded_intraprimed\t10\nexcluded_junction\t0\nmodels_out\t5\ntail_rate\t0.7647\n"
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
        let tail = s[4].lines().skip(10).collect::<Vec<_>>().join(" ");
        assert_eq!(
            tail,
            "models_known\t3 models_novel\t6 models_mono\t1 excluded_support\t0 excluded_fraction\t0 \
             excluded_tail\t0 excluded_intronic\t1 excluded_intraprimed\t0 excluded_junction\t0 models_out\t10 tail_rate\t0.7647"
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
        // NAGNAG (SIRV604/612): acceptor 403 has 2 of 3 reads at its own intron, but acceptor 400 is
        // also used by 10 reads of another donor; site-level wobble moved 403, intron-level keeps it
        let mut ng = String::new();
        for (i, (a, n)) in [(400, 1), (403, 2)].iter().enumerate() {
            for j in 0..*n {
                ng += &bed(
                    &format!("g{i}{j}"),
                    '+',
                    &[(100 + j, 200), (*a, 600)],
                    "PA30#PR30#IY990",
                );
            }
        }
        for j in 0..10 {
            ng += &bed(
                &format!("h{j}"),
                '+',
                &[(100 + j, 300), (400, 600)],
                "PA30#PR30#IY990",
            );
        }
        let o = call(&ng, "", &["--junction-wobble", "5"]).unwrap();
        assert!(o[4].contains("sites_corrected\t0\n"), "{}", o[4]);
        assert!(o[0].contains("\t2\t100,197\t0,303\n"), "{}", o[0]); // its own 2-read model
    }

    #[test]
    fn junction_support() {
        // cluster2 input: s0 is a singleton with novel introns no other read has, s1 a singleton
        // whose intron cluster h1 also carries, h0 a cluster (no #SG) with a unique intron
        let (t, sg) = ("PA30#PR30#IY990", "PA30#PR30#IY990#SG");
        let r = bed("s0", '+', &[(100, 200), (1000, 1100), (2000, 2100)], sg)
            + &bed("s1", '+', &[(5050, 5100), (5300, 5500)], sg)
            + &bed("h1", '+', &[(5000, 5100), (5300, 5400), (5600, 5700)], t)
            + &bed("h0", '+', &[(8000, 8100), (8300, 8400)], t);
        let names = |o: &str| {
            o.lines()
                .map(|l| l.split('\t').nth(3).unwrap()[..2].to_string())
                .collect::<Vec<_>>()
        };
        let o = call(&r, "", &["--preset", "sensitive"]).unwrap();
        assert_eq!(
            (names(&o[0]), o[5].as_str()),
            (vec!["s0".into(), "h1".into(), "s1".into(), "h0".into()], "")
        );
        let o = call(&r, "", &["--preset", "sensitive", "--junction-support"]).unwrap();
        assert_eq!(
            (names(&o[0]), names(&o[5])),
            (
                vec!["h1".into(), "s1".into(), "h0".into()],
                vec!["s0".into()]
            )
        );
        assert!(o[4].contains("excluded_junction\t1\n"), "{}", o[4]);
    }

    #[test]
    fn subchain_merge() {
        // container B (5 reads, 5 exons); F a protected 5′ cut inside B's second exon (10 reads);
        // G starts 48-50 bp into B's last intron (alternative first exon, 3 reads); H starts ~200 bp
        // before B (alternative TSS, 10 reads); K is a known chain. H and K carry tails, which keep
        // step 5 from absorbing their inner 3′ ends
        let (t, no) = ("PA30#PR30#IY990", "PA0#PR0#IY990");
        let (e2, e3, e4, e5) = ((1200, 1300), (1400, 1500), (1600, 1700), (1800, 2000));
        let mut r = String::new();
        for i in 0..5 {
            r += &bed(
                &format!("b{i}"),
                '+',
                &[(1000 + i, 1100), e2, e3, e4, e5],
                no,
            );
        }
        for i in 0..10 {
            r += &bed(&format!("f{i}"), '+', &[(1250 + i, 1300), e3, e4, e5], no);
            r += &bed(&format!("h{i}"), '+', &[(800 + i, 1100), e2], t);
        }
        for i in 0..3 {
            r += &bed(&format!("g{i}"), '+', &[(1550 + i, 1700), (1800, 1900)], no);
        }
        for i in 0..2 {
            r += &bed(&format!("k{i}"), '+', &[(1210 + i, 1300), (1400, 1450)], t);
        }
        let reference = bed("T1", '+', &[e2, e3], "");
        let run = |flags: &[&str]| {
            let flags: Vec<&str> = ["--preset", "sensitive"]
                .iter()
                .chain(flags)
                .copied()
                .collect();
            call(&r, &reference, &flags).unwrap()
        };
        let fate = |o: &[String]| -> HashMap<String, String> {
            o[2].lines()
                .skip(1)
                .map(|l| l.split_once('\t').unwrap())
                .map(|(read, m)| (read.split("__").next().unwrap().into(), m.into()))
                .collect()
        };
        // each model as <first letter of its representative>:<#CN>
        let models = |o: &[String]| -> Vec<String> {
            o[0].lines()
                .map(|l| l.split('\t').nth(3).unwrap())
                .map(|n| format!("{}:{}", &n[..1], n.rsplit('#').next().unwrap()))
                .collect()
        };
        let count = |o: &[String], k: &str| -> String {
            o[4].lines()
                .find_map(|l| l.strip_prefix(&format!("{k}\t")))
                .unwrap()
                .into()
        };

        // off: every chain is its own model
        let o = run(&[]);
        assert_eq!(models(&o), ["h:CN10", "b:CN5", "k:CN2", "f:CN10", "g:CN3"]);
        assert_eq!(count(&o, "chains_subchain_merged"), "0");

        // plain: every sub-chain joins B; reads, molecules and members follow
        let o = run(&["--subchain-merge", "plain"]);
        assert_eq!(models(&o), ["b:CN30"]);
        assert_eq!(
            (
                count(&o, "chains_subchain_merged"),
                count(&o, "chains_subchain_kept")
            ),
            ("4".into(), "0".into())
        );
        let row: Vec<&str> = o[1].lines().nth(1).unwrap().split('\t').collect();
        assert_eq!(&row[5..8], ["30", "30", "25"]); // n_reads, n_molecules, n_partial
        let f = fate(&o);
        assert_eq!((o[2].lines().count(), f.len()), (31, 30));
        assert!(f.values().all(|m| m == "t.1"), "{f:?}");

        // ends: F is a truncation and joins B; G (intron start), H (before B's start), K (known) stay
        let o = run(&["--subchain-merge", "ends"]);
        assert_eq!(models(&o), ["h:CN10", "b:CN15", "k:CN2", "g:CN3"]);
        assert_eq!(
            (
                count(&o, "chains_subchain_merged"),
                count(&o, "chains_subchain_kept")
            ),
            ("1".into(), "3".into())
        );
        let f = fate(&o);
        assert_eq!(
            ["b0", "f3", "h0", "k1", "g2"].map(|r| f[r].as_str()),
            ["t.2", "t.2", "t.1", "t.3", "t.4"]
        );
        let flags: Vec<&str> = o[1]
            .lines()
            .skip(1)
            .map(|l| l.rsplit('\t').next().unwrap())
            .collect();
        assert_eq!(
            flags,
            ["subchain_kept", ".", "subchain_kept", "subchain_kept"]
        );

        // known needs ends: K's reads end inside B's exons, so K joins B like F; G and H keep their ends
        let o = run(&["--subchain-merge", "ends", "--subchain-known-needs-ends"]);
        assert_eq!(models(&o), ["h:CN10", "b:CN17", "g:CN3"]);
        assert_eq!(
            (
                count(&o, "chains_subchain_merged"),
                count(&o, "chains_subchain_kept")
            ),
            ("2".into(), "2".into())
        );

        // thresholds: G's 3 reads are too few at 4; H's ~200 bp fit a 300 bp end tolerance
        let o = run(&["--subchain-merge", "ends", "--subchain-min-end-reads", "4"]);
        assert_eq!(models(&o), ["h:CN10", "b:CN18", "k:CN2"]);
        let o = run(&[
            "--subchain-merge",
            "ends",
            "--subchain-end-tolerance",
            "300",
        ]);
        assert_eq!(models(&o), ["b:CN25", "k:CN2", "g:CN3"]);
        // G's reads reach 50, 49 and 48 bp into B's intron: all 3 count at 47, only 2 at 48
        let o = run(&[
            "--subchain-merge",
            "ends",
            "--subchain-intron-tolerance",
            "47",
        ]);
        assert_eq!(models(&o).last().unwrap(), "g:CN3");
        let o = run(&[
            "--subchain-merge",
            "ends",
            "--subchain-intron-tolerance",
            "48",
        ]);
        assert_eq!(models(&o), ["h:CN10", "b:CN18", "k:CN2"]);
    }

    #[test]
    fn reference_ends() {
        // FLNC-like reads (no tail tags) against CDS-only references (thick = bounds, as `bed` writes)
        let no = "PA0#PR0#IY990";
        let mut r = String::new();
        let mut reference = String::new();
        // 3′ cut S (5 reads, a 5′ UTR intron before reference R1's chain) ends past R1's stop inside
        // the fourth exon of L (2 reads): S outnumbers L's own reads and stays apart
        reference += &bed("R1", '+', &[(1000, 1100), (1200, 1300), (1400, 1500)], "");
        let s_ex = |i: u32| [(800 + i, 900), (1000, 1100), (1200, 1300), (1400, 1560)];
        let l_ex = |i: u32| {
            [
                (800 + i, 900),
                (1000, 1100),
                (1200, 1300),
                (1400, 1600),
                (1700, 1800 + 10 * i),
            ]
        };
        for i in 0..5 {
            r += &bed(&format!("s{i}"), '+', &s_ex(i), no);
        }
        for i in 0..2 {
            r += &bed(&format!("l{i}"), '+', &l_ex(i), no);
        }
        // the same cut with fewer reads (T, 2) than the longer chain's own (M, 3) is absorbed
        reference += &bed(
            "R2",
            '+',
            &[(11000, 11100), (11200, 11300), (11400, 11500)],
            "",
        );
        for i in 0..2 {
            let t = [(11000 + i, 11100), (11200, 11300), (11400, 11560)];
            r += &bed(&format!("t{i}"), '+', &t, no);
        }
        for i in 0..3 {
            let m = [
                (11000 + i, 11100),
                (11200, 11300),
                (11400, 11600),
                (11700, 11800),
            ];
            r += &bed(&format!("m{i}"), '+', &m, no);
        }
        // 5′ cuts sharing a reference's first intron. A bare CDS start keeps nothing apart: P
        // starts 90 bp before the CDS start of R3, Q 20 bp inside the CDS of R4, both are absorbed.
        // R6 has a 5′ UTR: H starts 5 bp before its annotated 5′ end and stays apart
        for (base, name, cut, five) in [
            (21000, 'c', 'p', 21210),
            (31000, 'd', 'q', 31320),
            (61000, 'g', 'h', 61210),
        ] {
            let (e3, e4) = ((base + 500, base + 600), (base + 700, base + 900));
            let row = bed("R", '+', &[(base + 300, base + 400), e3, e4], "");
            reference += &if cut == 'h' {
                bed("R", '+', &[(base + 215, base + 400), e3, e4], "")
                    .replace("+\t61215", "+\t61300")
            } else {
                row
            };
            for i in 0..4 {
                let c = [(base + i, base + 100), (base + 200, base + 400), e3, e4];
                r += &bed(&format!("{name}{i}"), '+', &c, no);
            }
            for i in 0..2 {
                r += &bed(
                    &format!("{cut}{i}"),
                    '+',
                    &[(five + i, base + 400), e3, e4],
                    no,
                );
            }
        }
        // the 3′ cut again on the minus strand (3′ = start): Y (5 reads) inside the second exon of Z (2)
        reference += &bed(
            "R5",
            '-',
            &[(51400, 51500), (51600, 51700), (51800, 51900)],
            "",
        );
        let (e3, e4) = ((51600, 51700), (51800, 51900));
        for i in 0..5 {
            let y = [(51340, 51500), e3, e4, (52000, 52100 + i)];
            r += &bed(&format!("y{i}"), '-', &y, no);
        }
        for i in 0..2 {
            let z = [(51100, 51200), (51300, 51500), e3, e4, (52000, 52100 + i)];
            r += &bed(&format!("z{i}"), '-', &z, no);
        }
        // each model as <first letter of its representative>:<#CN>[:<flags>]
        let models = |o: &[String]| -> Vec<String> {
            let flags: HashMap<&str, &str> = o[1]
                .lines()
                .skip(1)
                .map(|l| l.split('\t').collect::<Vec<_>>())
                .map(|f| (f[1], f[11]))
                .collect();
            o[0].lines()
                .map(|l| l.split('\t').nth(3).unwrap())
                .map(|n| {
                    let rep = n.rsplit_once("#CN").unwrap().0;
                    let cn = n.rsplit('#').next().unwrap();
                    format!("{}:{cn}:{}", &n[..1], flags[rep])
                })
                .collect()
        };
        let o = call(&r, &reference, &[]).unwrap();
        assert_eq!(
            models(&o),
            [
                "s:CN5:.",
                "l:CN2:.",
                "m:CN5:.",
                "c:CN6:.",
                "d:CN6:.",
                "z:CN2:.",
                "y:CN5:.",
                "g:CN4:.",
                "h:CN2:protected"
            ]
        );
        // step 11 keeps S for the same reason (it is not a known chain); plain merges everything
        let o = call(&r, &reference, &["--subchain-merge", "ends"]).unwrap();
        assert_eq!(
            models(&o),
            [
                "s:CN5:subchain_kept",
                "l:CN2:.",
                "m:CN5:.",
                "c:CN6:.",
                "d:CN6:.",
                "z:CN2:.",
                "y:CN5:subchain_kept",
                "g:CN4:.",
                "h:CN2:protected,subchain_kept"
            ]
        );
        let o = call(&r, &reference, &["--subchain-merge", "plain"]).unwrap();
        assert_eq!(
            models(&o),
            ["l:CN7:.", "m:CN5:.", "c:CN6:.", "d:CN6:.", "z:CN7:.", "g:CN6:."]
        );

        // 3′: an annotated end (UTR) fits within 50 nt; a bare CDS stop fits when the model reaches it
        assert!(reaches(1050, (1000, false), true) && !reaches(1051, (1000, false), true));
        assert!(reaches(900, (1000, true), true) && reaches(1005, (1000, true), true));
        assert!(!reaches(1006, (1000, true), true));
        assert!(reaches(1100, (1000, true), false) && !reaches(994, (1000, true), false));

        // read fraction: U, V and W tie at 2 reads behind X (16 of 22 locus reads); they are kept
        // or dropped together, whatever their order
        let mut r = String::new();
        for i in 0..16 {
            let x = [(41000 + i, 41100), (41200, 41300), (41400, 41500)];
            r += &bed(&format!("x{i}"), '+', &x, no);
        }
        for (name, acceptor) in [('u', 41210), ('v', 41220), ('w', 41230)] {
            for i in 0..2 {
                let e = [(41000 + i, 41100), (acceptor, 41300), (41400, 41500)];
                r += &bed(&format!("{name}{i}"), '+', &e, no);
            }
        }
        let o = call(&r, "", &["--min-read-fraction", "0.9"]).unwrap();
        assert_eq!(models(&o), ["u:CN2:.", "v:CN2:.", "w:CN2:.", "x:CN16:."]);
        let o = call(&r, "", &["--min-read-fraction", "0.7"]).unwrap();
        assert_eq!(models(&o), ["x:CN16:."]);
        assert!(o[4].contains("excluded_fraction\t6\n"), "{}", o[4]);

        // --sample-read-fraction: sample b expresses U (3 of its 6 locus reads); the 20 reads of X
        // in the pool put U below the fraction, the rule within sample b keeps it. V (one read in
        // each sample) is below the fraction in the pool and in both samples.
        let x = |i: u32| [(51000 + i, 51100), (51200, 51300), (51400, 51500)];
        let alt = |i: u32, acceptor: u32| [(51000 + i, 51100), (acceptor, 51300), (51400, 51500)];
        let (mut a, mut b) = (String::new(), String::new());
        for i in 0..20 {
            *(if i < 18 { &mut a } else { &mut b }) += &bed(&format!("x{i}"), '+', &x(i), no);
        }
        for i in 0..3 {
            b += &bed(&format!("u{i}"), '+', &alt(i, 51210), no);
        }
        a += &bed("v0", '+', &alt(0, 51220), no);
        b += &bed("v1", '+', &alt(1, 51220), no);
        let pool = [("a.bed", a.as_str()), ("b.bed", b.as_str())];
        let o = call_beds(&pool, "", &["--min-read-fraction", "0.8"]).unwrap();
        assert_eq!(models(&o), ["x:CN20:."]);
        let flags = ["--min-read-fraction", "0.8", "--sample-read-fraction"];
        let o = call_beds(&pool, "", &flags).unwrap();
        assert_eq!(models(&o), ["u:CN3:.", "x:CN20:."]);
        assert!(o[4].contains("excluded_fraction\t2\n"), "{}", o[4]);
        // one input file: the flag changes nothing
        let one = a.clone() + &b;
        let o = call_beds(&[("a.bed", one.as_str())], "", &flags).unwrap();
        assert_eq!(models(&o), ["x:CN20:."]);
    }

    #[test]
    fn cds_pool() {
        // U (one molecule, a PCR twin) and W (one read) hold the CDS of reference R behind
        // different 5′ UTR introns; each is short of support, together they are two molecules,
        // and U has fewer unannotated introns
        let no = "PA0#PR0#IY990";
        let cds = [(70300, 70400), (70500, 70600), (70700, 70900)];
        let reference = bed("R", '+', &cds, "");
        let tail = [(70200, 70400), (70500, 70600), (70700, 70950)];
        let mut r = String::new();
        for i in 0..2 {
            let u = [&[(70000, 70100)][..], &tail[..]].concat();
            r += &bed(&format!("u{i}"), '+', &u, no);
        }
        let w = [&[(69800, 69900), (70010, 70100)][..], &tail[..]].concat();
        r += &bed("w0", '+', &w, no);
        // each model as <first letter of its representative>:CN<reads>:<flags>
        let models = |o: &[String]| -> Vec<String> {
            let support = o[1]
                .lines()
                .skip(1)
                .map(|l| l.split('\t').collect::<Vec<_>>());
            let names = o[0].lines().map(|l| l.split('\t').nth(3).unwrap());
            names
                .zip(support)
                .map(|(n, f)| format!("{}:CN{}:{}", &n[..1], f[5], f[11]))
                .collect()
        };
        let o = call(&r, &reference, &[]).unwrap();
        assert_eq!(models(&o), [""; 0]);
        let o = call(&r, &reference, &["--cds-pool"]).unwrap();
        assert_eq!(models(&o), ["u:CN2:cds_pooled"]);
        assert!(o[4].contains("models_known\t1\n") && o[4].contains("excluded_support\t1\n"));
        // one molecule in all is not enough
        let one: String = r.lines().take(2).map(|l| format!("{l}\n")).collect();
        let o = call(&one, &reference, &["--cds-pool"]).unwrap();
        assert_eq!(models(&o), [""; 0]);
    }

    #[test]
    fn reference_cds_representative() {
        // CDS-only references (thick = bounds); FLNC-like reads in two input files (two samples)
        let no = "PA0#PR0#IY990";
        let (mut a, mut b) = (String::new(), String::new());
        let mut reference = String::new();
        // A (+): six reads end at an A-run 20 bp before the stop of RA, two reach past it. The 3′
        // mode stays the consensus end, the representative is the closer of the two that hold the CDS
        reference += &bed("RA", '+', &[(1000, 1100), (1200, 1300), (1400, 1500)], "");
        let ex = |s: u32, e: u32| [(s, 1100), (1200, 1300), (1400, e)];
        for i in 0..5 {
            a += &bed(&format!("a{i}"), '+', &ex(990 + i, 1480), no);
        }
        b += &bed("a5", '+', &ex(990, 1480), no);
        b += &bed("h0", '+', &ex(990, 1520), no);
        b += &bed("h1", '+', &ex(991, 1530), no);
        // B (-, 5′ = end): the trimmed 5′ end falls 20 bp inside the CDS of RB, two reads reach its start
        reference += &bed("RB", '-', &[(2000, 2100), (2200, 2300), (2400, 2500)], "");
        let ex = |e: u32| [(1990, 2100), (2200, 2300), (2400, e)];
        a += &bed("g0", '-', &ex(2510), no);
        a += &bed("g1", '-', &ex(2505), no);
        for i in 0..18 {
            a += &bed(&format!("n{i}"), '-', &ex(2480), no);
        }
        // C (+): as A with one read past the stop of RC: one read is no evidence, nothing changes
        reference += &bed(
            "RC",
            '+',
            &[(11000, 11100), (11200, 11300), (11400, 11500)],
            "",
        );
        let ex = |s: u32, e: u32| [(s, 11100), (11200, 11300), (11400, e)];
        for i in 0..5 {
            a += &bed(&format!("c{i}"), '+', &ex(10990 + i, 11480), no);
        }
        a += &bed("k0", '+', &ex(10990, 11520), no);
        // D: a novel chain with one read in each sample
        a += &bed("d0", '+', &[(30000, 30100), (30200, 30300)], no);
        b += &bed("d1", '+', &[(30001, 30100), (30200, 30300)], no);

        let o = call_beds(&[("x/a.bed", &a), ("b.bed", &b)], &reference, &[]).unwrap();
        let col = |text: &str, cols: &[usize]| -> Vec<String> {
            text.lines()
                .map(|l| l.split('\t').collect::<Vec<_>>())
                .map(|f| cols.iter().map(|&c| f[c]).collect::<Vec<_>>().join(" "))
                .collect()
        };
        assert_eq!(
            col(&o[0], &[1, 2, 3, 5]),
            [
                "990 1520 h0__FC0#TC0#PA0#PR0#IY990#CN8 +",
                "1990 2505 g1__FC0#TC0#PA0#PR0#IY990#CN20 -",
                "10990 11480 c0__FC0#TC0#PA0#PR0#IY990#CN6 +",
                "30000 30300 d0__FC0#TC0#PA0#PR0#IY990#CN2 +",
            ]
        );
        // five_prime / three_prime are still the consensus ends
        assert_eq!(
            col(&o[1], &[0, 4, 9, 10])[1..],
            [
                "t.1 known 990 1480",
                "t.2 known 2480 1990",
                "t.3 known 10990 11480",
                "t.4 novel 30000 30300",
            ]
        );
        assert_eq!(
            o[6],
            "model\tname\tsample\tn_reads\tn_full\tn_full_molecules\n\
             t.1\th0__FC0#TC0#PA0#PR0#IY990#CN8\ta\t5\t5\t5\n\
             t.1\th0__FC0#TC0#PA0#PR0#IY990#CN8\tb\t3\t3\t3\n\
             t.2\tg1__FC0#TC0#PA0#PR0#IY990#CN20\ta\t20\t20\t3\n\
             t.3\tc0__FC0#TC0#PA0#PR0#IY990#CN6\ta\t6\t6\t6\n\
             t.4\td0__FC0#TC0#PA0#PR0#IY990#CN2\ta\t1\t1\t1\n\
             t.4\td0__FC0#TC0#PA0#PR0#IY990#CN2\tb\t1\t1\t1\n"
        );
        assert_eq!(
            o[7],
            "chrom\tstart\tend\tstrand\tn_reads\tannotated\n\
             chr1\t1100\t1200\t+\t8\t1\nchr1\t1300\t1400\t+\t8\t1\n\
             chr1\t2100\t2200\t-\t20\t1\nchr1\t2300\t2400\t-\t20\t1\n\
             chr1\t11100\t11200\t+\t6\t1\nchr1\t11300\t11400\t+\t6\t1\n\
             chr1\t30100\t30200\t+\t2\t0\n"
        );
    }
}
