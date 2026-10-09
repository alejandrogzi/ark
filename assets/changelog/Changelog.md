<p align="center">
  <p align="center">
    <img width=100 align="center" src="../figures/logo.png" >
  </p>

<p align="center">
  <picture>
    <source
      media="(prefers-color-scheme: dark)"
      srcset="../figures/hillerlab-dark.png"
    >
    <source
      media="(prefers-color-scheme: light)"
      srcset="../figures/hillerlab-light.png"
    >
    <img
      width="200"
      alt="Hiller Lab"
      src="../figures/hillerlab-light.png"
    >
  </picture>
</p>


  <span>
    <h1 align="center">
        ark
    </h1>
  </span>

  <span>
    <h2 align="center">
        CHANGELOG
    </h2>
  </span>

  <p align="center">
    <a href="https://github.com/hillerlab/ark" reference="_blank">
      <img alt="GitHub License" src="https://img.shields.io/github/license/hillerlab/ark?color=blue">
    </a>
  </p>

  <p align="center">
    <samp>
        <span> The Hiller Lab at the Senckenberg Research Institute </span>
        <br>
        <br>
        <a href="https://github.com/alejandrogzi/ark/blob/master/assets/docs/usage.md">usage</a> .
        <a href="https://github.com/hillerlab/ark/blob/main/assets/pipeline/ark.mermaid">pipeline</a> .
        <a href="https://hillerlab.com/">us</a> 
    </samp>
  </p>

</p>

# Changelog

## [v2.1.3] - unreleased

- **minimap2 now gets `--cs=<tag>` instead of `-cs <tag>`.** `-cs long` is parsed as `-c -s long` (peak DP score 0, no `cs:Z` tag), so pass 1 and the fragment pass reported low-scoring split alignments that should have been dropped. Both `MINIMAP2_ALIGN` blocks now pass `--cs=${params.minimap2_align_cigar_tag}`.
- **Veredict buckets are exclusive (veredict 0.0.9).** RT routing takes precedence over every other bucket and artifact routing over flaw-based bucketing, so each row lands in exactly one bucket. Before, the artifact mask was recomputed from the unfiltered schema and discarded the RT exclusion, so RT models leaked into the flaw buckets.
- **`artifacts.bed` is published end to end.** `VEREDICT` now emits `veredict/*.artifacts.bed`, and `POLISH` joins, detaches duplicates, converts (`BEDTOBIGBED_ARTIFACTS`) and mixes it into the published bigBeds (`12_POLISH/BED`, `BB`); the detach/join publish prefix covers `JOIN_VEREDICT_ARTIFACTS`.
- **`-resume` no longer skips POLISH after APARENT.** `APARENT_PREDICT` bedGraphs are mandatory outputs now, and a chunk with no signal on one strand still emits an (empty) file. `BEDGRAPHTOBIGWIG` deletes the bedGraph at its real path when `bigtools_keep_bedgraph` is false, so on `-resume` Nextflow saw the output missing and reran APARENT instead of emitting nothing; empty bedGraphs now produce no bigWig (`bigwig` output is optional).
- **`BEDTOBIGBED` passes `-p no`.** bigtools 0.5.6 binary-searches chromosome boundaries once a BED is >= 200 MB, skips a short chromosome between two long ones, then panics with "File is not sorted". Reading in order avoids it.

## [v2.1.2] - unreleased

- **`PBCCS` uses its task directory as `TMPDIR`.** ccs names its temp files `thread.<i>_<j>.<chunk>.bam`, so two runs sharing a `TMPDIR` (e.g. a lab-wide scratch) overwrote each other's files and failed. Seen with two giraffe runs on the same subreads.
- **CI gold regenerated with the published v2.1.1 images.** The v2.1.1 gold came from collapse 0.1.0 and the old intron tables, so `master` CI failed on `flnc`.
  - The clustered `flnc` input (`sensitive`) now drops its 5 single-read novel chains (`excluded_junction`; 8 → 3 models). xORF and NMD outputs shrink to match.
  - Both entrypoints gain an `excluded_junction` counter and a `rejected.bed`, which is empty for `subreads`.
  - The `subreads` intron table keeps its 46 introns but loses 5 duplicate rows.

## [v2.1.1] - 2026-10-07

These fixes were found by running the containerized CI and real SRA data against the published v2.1.0 images.

- **Regression in v2.1.0: POLISH wrote no `pass` outputs.** `PREPOLISH` fed the uncollapsed reads to APARENT too. APARENT scans 3′ UTRs (ORF thick end to transcript end), and reads have no ORF yet, so no chunks were made. The PAS caller and the verdict step then never ran.
  - `PREPOLISH` now takes two inputs: the ORF-annotated models (`reads`, for APARENT) and the uncollapsed reads (`evidence`, for intron frequencies).
  - With `reconstruct_engine none` both inputs are the same channel, as in v2.0.x.
- **lima hung on SRA-renamed reads.** lima 26.2.1 stalls on FASTQ whose names are not `movie/zmw/ccs` (SRA's `SRRxxx.N`), even with `--per-read`.
  - The new `FASTX_RENAME` gives raw-CCS reads CCS names (`<sample>/<n>/ccs`) before `LIMA_FASTX`. The sample as movie keeps names unique across pools.
  - On SRR27664179 (raw CCS), 3,400 of 3,783 reads pass, and the output inspects as `fl` (primer rate 0, polyA3 0.996, polyT5 0).
  - `LIMA_FASTX` now publishes only its `*.lima.*` reports.
- **CI gold regenerated** from a containerized run on the published images:
  - only files whose line counts changed are updated;
  - new outputs added: `05_FASTX` and `07B_RECONSTRUCT`;
  - outputs of the cluster2 path removed (subreads now uses `cluster_engine none`).
- **Intron tables now count each read once.** v2.0.x counted a read once per predicted ORF, inflating `seen / spanned`.
- **Harness:** the lima stub asserts CCS-style names, so the renaming cannot silently disappear.
- **`collapse chain` 0.1.1: junction correction moves whole introns and protects well-supported ones.**
  - **Finding:** on LRGASP WTC11 SIRV-Set 4, 0.1.0 merged two real NAGNAG-type isoforms. SIRV604 vs SIRV612 is a 3 bp acceptor shift with 164 reads, and SIRV307 has 58. The guard was evaluated per site, so a minor acceptor was compared with every read at a major acceptor shared by several isoforms.
  - **Introns move as a unit.** An intron moves only onto an existing read or reference intron within 5 bp at both ends that has at least 5× its support, so a corrected model only carries introns that exist.
  - **Protection.** An intron with ≥ 20 reads, or ≥ 20% of its window, never moves. Annotated sites never move, and an annotated intron weighs at least as much as the reads at its weaker site.
  - **SIRV (60 multi-exon isoforms), recall / precision / F1:**
    - `balanced`: 0.933 / 0.824 / 0.875 before, 0.967 / 0.795 / 0.872 now;
    - `strict`: 0.817 / 0.925 / 0.867 before, 0.850 / 0.944 / 0.895 now.
  - **Other data:** mouse chr19 is unchanged (607 models). On giraffe NC_137277.1 (1.74M reads), models go from 25,585 to 25,732 (+0.6%).
  - **Images:** `isox-rs` (collapse 0.1.1, aparent 0.0.4) and `isotools` (v0.0.46) modules now use `:latest`, matching `xloci:latest`.
- **New option `reconstruct_junction_support` (`collapse chain --junction-support`), on by default; `--reconstruct_junction_support false` turns it off.** A novel chain resting on one molecule (a cluster2 singleton, `#SG`; or any single read when the input has no singletons) is dropped unless each of its novel introns is carried by another read. The dropped reads are written to `<prefix>.rejected.bed` (BED12, published next to the models) for the browser, and counted as `excluded_junction`.
  - Why: with clustered input, `sensitive` keeps every chain. On the giraffe muscle slice, 2,799 of 2,905 cluster2-only models rested on one singleton; their novel junctions were flagged by PREPOLISH in 18% of models and non-canonical in 7%. A plain "≥ 2 molecules" rule would also drop 905 annotated chains.
  - A record without `#SG` in a clustered input is a cluster of ≥ 2 reads and is never dropped.
  - Benchmarks (all tools on the same 20,465 SIRV reads; giraffe `bam_isoseq` with and without the option):
    - SIRV `sensitive`: precision 0.311 → 0.541 at the same recall (0.983); after cluster2, 0.365 → 0.667 at recall 0.967. `balanced` and `strict` are unchanged: they already need 2 reads per novel chain.
    - Giraffe: 5,250 → 3,578 models, POLISH pass 3,192 → 2,309, xORF input 4,257 → 2,148 (`none` path: 2,166). All 1,628 known models are kept, as are 829/831 annotated-intron extras and 373/375 independently supported extras. 98.6% of flagged and 94.6% of non-canonical extras are dropped.
    - Cost: unconfirmed single-molecule canonical chains go too (624, of which 402 pass POLISH; 231 use new splice sites). Turn the option off when rare novel isoforms matter more than precision.
  - Needs collapse 0.1.1 (`isox-rs:latest` once published); the v2.1.0 image has no `--junction-support`.
- **`ISOTOOLS_INTRON_RETENTION` panicked on corrected models** (`Intron not found`; giraffe genome).
  - Wobble correction moves a model's donor and acceptor independently, so a model can carry an intron that no uncollapsed read has, and the read-built intron table lacks it.
  - With a reconstruction engine on, iso-intron now runs with `--allow-missing`: it warns and skips such introns. Their sites carry at least 5× more read support, so they are not RT-artifact candidates.
  - Root cause, fixed by collapse 0.1.1. On giraffe NC_137275.1, `1049593-1051084(+)` is in no read: 0.1.0 moved a 4-read donor (1049588) onto 1049593, which 538 reads use with another acceptor. 0.1.0 made 2 such introns on that chromosome (of 19,164), 0.1.1 none. `--allow-missing` stays as a guard.
- **`XLOCI_EXTRACT_INTRONS` no longer fails on introns near a contig end.** xloci 0.0.6 panicked with "Feature coordinate 54 is underflowing by 100 bases" when an intron sat closer than the 100 bp flank to a contig start; this was seen on an ERCC spike-in contig and can also happen on short scaffolds. The step now runs with `--ignore-errors`, so only that intron is skipped.
- **`-resume` no longer reruns alignment and everything after it** (present since v2.0.0). `SAMTOOLS_BAM` deleted the SAM it got from `ARK_ALIGN`, so the cached alignment task had a missing output on every resume.
  - `ARK_ALIGN`, `ARK_ALIGN_FRAGMENTS`, `MINIMAP2_ALIGN` and `FLAIR_ALIGN` now pipe minimap2 into `samtools sort` and index in the same task (biocontainers mulled image: minimap2 2.31, samtools 1.23.1). No SAM is written and the `SAMTOOLS_BAM` step is gone for them; deSALT keeps it.
  - BAMs are published where they were (`06_ARK_ALIGN/BAM`, `06_ARK_ALIGN/FRAGMENTS`, …); the `SAM` symlink folders are gone. `minimap2_align_keep_sam` now applies to deSALT only.
  - The flnc e2e matches gold, and a second `-resume` run caches all 75 tasks.
- **Intron tables hold each intron once.** xloci wrote one row per read and intron; on the giraffe pilot, NC_137281.1 had 17,056,867 rows for 59,054 introns (32 GB), and intronIC was killed at 72 GB. `XLOCI_EXTRACT_INTRONS` now keeps the first row per intron name (the name is the coordinates). intronIC only needs each intron once, and iso-classify takes frequencies from the reads, so results are unchanged; the CI gold already had unique rows.
- **APARENT no longer fails on transcripts that end near a contig end** (`aparent` 0.0.4, `assets/rust/aparent`). The 100 bp downstream flank ran past the end of short contigs ("Interval 944-1151 exceeds chromosome length 1051", ERCC-00130) and the whole chunking task exited. The flank is now clipped at the contig end. Ships with the next `isox-rs` image.
  - `XISO_APARENT_CHUNK` now uses its own `aparent_chunker_{upstream,downstream}_flank` params instead of the xloci ones (both default to 100).
- **Genome FASTA headers are cut to their first word** (new `FASTA_CLEAN`, which replaces `GUNZIP_FASTA` and also runs on plain FASTA).
  - Before, `XLOCI_EXTRACT_INTRONS` failed with "Chromosome chr12 ... not found in genome" on any FASTA whose headers carry descriptions, as NCBI, Ensembl and LRGASP headers do: xloci keyed sequences by the whole header line, while minimap2 and the other tools key them by the first word.
  - IUPAC ambiguity codes become N/n, case kept: xloci 0.0.6 panicked on them (`ERROR: Invalid base`, GRCh38 chr21, which has 94 such bases genome-wide).
  - Output names are unchanged.
- **xORF `PREDICT` no longer fails on chunks without DIAMOND hits** (fixed in the `modules/xorf` working tree, `modules/predict/predict.py`). Synthetic ERCC spike-ins, and any small rescued chunk without protein hits, made pandas raise `EmptyDataError`. An empty DIAMOND or RNAsamba table now yields empty prediction files. Rebuild `ghcr.io/hillerlab/orf-predict` to ship it.
- **CI references labelled correctly.** The CI genome is human chr19 (hg38). Its annotation is a mouse-to-human TOGA projection, so it carries mouse gene names in hg38 coordinates.
  - The CI passed `--global_species_name mm39` and mm39 selenocysteine codons, so xORF could never match a Sec codon on chr19 (e.g. GPX4).
  - `assets/test_data/data/hg38.selenocysteine.bed` (47 Sec codons from GENCODE v38) replaces `mm39.selenocysteine.bed`, and the species is now `hg38`.
  - No CI read overlaps a selenoprotein, so the gold is unchanged; both e2e runs still match it.

## [v2.1.0] - 2026-10-07

BREAKING release. Transcript models, not individual reads, now reach ORF prediction. Refined reads are no longer clustered before alignment. Every FASTA/FASTQ input is classified and cleaned before alignment. Design, evidence and pilot plan: `PLAN.md` and `reports/2026-10-07_isoseq_fastq_reconstruction_review.md`.

### BREAKING: `reconstruct_engine chain` (default) collapses reads into transcript models before ORF prediction

- New `RECONSTRUCT` subworkflow (`src/subworkflows/reconstruct/main.nf`) runs after polyA segmentation and fusion detection.
  - The fusion detector still reads reads: its recover rule is a read ratio.
  - Its free and fusion outputs are each collapsed per `[sample_id, chr]` by the new `CHAIN_COLLAPSE` module (`collapse chain`, collapse crate 0.1.0). Only the models go to xORF.
- **Models:** each model is a real read line (iso-segment tags intact) with its support appended as `#CN<n>`; `#SG` now means a one-read model.
- **`collapse chain`** (`assets/rust/collapse/src/chain.rs`) works per strand.
  - **Splice sites:** an unannotated site within 5 bp of an annotated site, or of one with 5× its reads, moves to it. Annotated sites, and sites holding ≥ 20% of their window, never move.
  - **Chains and molecules:** reads are grouped by intron chain. Reads with identical ends count as one molecule (PCR twins).
  - **Truncated copies:** 5′- and tailless 3′-truncated copies are absorbed into the best-supported compatible parent. Tight TSS clusters (≥ 10 starts within 50 nt) and annotated starts are protected.
  - **Mono-exonic reads:** absorbed when inside a model's exon, excluded when intronic, clustered otherwise.
  - **Support:** a known chain needs 1 molecule, a novel chain 2 (`balanced`). When ≥ 70% of reads carry a polyA tail, a novel chain also needs one tailed read. Novel models must fall within the first 99% of their locus reads, ranked by support (loci = models sharing exonic bp).
  - **Ends:** the 3′ end is the mode of tailed 3′ ends (the median when no end repeats). The 5′ end is the most upstream after trimming the outer 10%.
  - **Outputs:** `models.bed`, `support.tsv`, `members.tsv.gz`, `excluded.bed` (reason in column 13) and `counts.tsv`.
  - **Speed:** about 3.9 M reads in 5 s and 1.2 GB on a laptop.
- **Outputs:** every read is accounted for, either in a model's member list or in `07B_RECONSTRUCT/chain/<sample>/*.excluded.bed` with its reason.
- **Presets:** `reconstruct_preset` picks `sensitive` (only removes redundancy, keeps every distinct chain), `balanced` (default) or `strict`. `reconstruct_min_support_novel`, `reconstruct_min_support_mono`, `reconstruct_min_read_fraction` and `reconstruct_junction_wobble` override a preset. Clustered inputs always run `sensitive`.
- **Intron classification** counts reads (`seen / spanned`), so `PREPOLISH` now takes the uncollapsed free reads (new `evidence` emit of `SPLIT_ALIGN_CLEAN_CHUNKS`), re-keyed to the xORF ids. The POLISH join keys are plain strings now; xORF builds GString ids.
- **`collapse_shrink_twins`** is ignored with a reconstruction engine (warning). `reconstruct_engine none` restores per-read ORF calling and warns.
- **IsoQuant/isocall:** `isoquant` and `isocall` are reserved engine names. Validation rejects them ("not implemented in v2.1.0"). `RECONSTRUCT` documents the contract a future engine has to meet.

### BREAKING: `cluster_engine none` is the default; `cdhit` and `rattle` are removed

- **`none`:** refined (FLNC) BAMs become FASTA once (`SAMTOOLS_FASTA`) and pool through `POOL_READS` like `flnc` inputs. This is the order of PacBio's current Kinnex workflow: FLNC → align → call.
- **`isoseq`** (cluster2) remains selectable but is no longer used by default. Its output is marked `clustered`.
- **Removed:**
  - the `CDHIT_EST` and `RATTLE` modules, their config and `cdhit_identity`;
  - their CI cases.

  Passing `cluster_engine cdhit|rattle` or `cdhit_identity` fails validation.
- **`--entrypoint cluster`** reads `*_flnc.bam` only again. FASTA/FASTQ go to `flnc`.

### FASTA/FASTQ inputs are classified and normalized per file (`FASTX_PREPARE`)

- **New `FASTX_PREPARE` subworkflow:** `iso-fastx inspect` (isotools v0.0.45) classifies every `flnc` file from its first reads, and the file is routed by state:

  | State | What happens |
  |---|---|
  | `ccs` (primers, both orientations) | `LIMA_FASTX` (`lima --isoseq --peek-guess`, FASTQ in and out); the per-pair outputs are joined per input |
  | `mixed` (primer-free, partly reversed) | `FASTX_ORIENT` (`iso-fastx orient`) flips reads that start with polyT |
  | `fl` / `flnc` | used as is |
  | `clustered` | used as is, with the `sensitive` preset |
  | `empty` | skipped with a warning |
  | `subreads` / `ambiguous` | the run stops with the measured numbers and a hint |

- **Primers:** `global_primers` when given, else the detected kit picks a bundled set: `assets/primers/isoseqx.fasta` (official `IsoSeq_v2_primers_12.fasta`) or `assets/primers/express.fasta` (NEB/Clontech).
- **`flnc_input_state`** forces one state. `ccs` then needs `--global_primers`.
- **Why:** ARK aligns with `minimap2 -uf` and takes the strand from the alignment flag, which assumes oriented reads. Raw-CCS SRA FASTQ is about 50/50 forward/reverse.
- **`isotools_adapter_remove_adapters` now defaults to `false`.** Inputs are cleaned before alignment; the module stays available.

### Read-count-aware polishing

- iso-orphan, iso-utr and iso-pas (isotools v0.0.45) weight records by `#CN`, so their support thresholds and ratios mean the same on models as on reads. Untagged reads weigh 1.

### Images

- `isotools` modules use `ghcr.io/alejandrogzi/isotools:v0.0.45`. `isox-rs` and `isox-py` modules use `:v2.1.0` instead of `:latest`. Publish both tags before running v2.1.0.

### Chores

- Collapse crate 0.1.0. `test_binkey_to_bytes` expected 31 bytes and already failed on v2.0.28; it now expects 39.
- isotools v0.0.45 brings `iso-fastx`, the NEB/IsoSeqX 5′ primers in iso-adapter's database, and `#CN` weights; see its changelog.

### CI

- `assets/ci/test_isoseq.py`:
  - drops the cd-hit/RATTLE cases;
  - adds `cluster_engine none` cases and asserts one `collapse chain` call per result group with the right preset;
  - adds a `flnc` routing case (ccs → lima, mixed → orient, clustered → sensitive) and a subreads rejection;
  - adds validator checks for the removed engines and the reserved reconstruction engines.
- The containerized gold has to be regenerated from a real run once the images are published.

## [v2.0.28] - 2026-10-07

- New `cluster_engine` (`isoseq` default | `cdhit` | `rattle`) de-duplicates reads without PacBio tags (e.g. SRA IsoSeq FASTQs), which `isoseq cluster2` cannot use. The `cluster` entrypoint now also reads `*.fasta[.gz]` / `*.fastq[.gz]` (rejected with a hint under `isoseq`), and every Iso-Seq entrypoint can use the new engines. Neither tool reads BAM, so BAMs pass once through the new `SAMTOOLS_FASTA` (FASTA, because neither tool uses qualities: RATTLE's consensi are identical from FASTA and FASTQ). FASTA/FASTQ inputs are never converted. `multi_sample` pools without a separate merge: RATTLE takes every file in one `-i a,b,c`, and `CDHIT_EST` streams them into the one plain FASTA that `cd-hit-est` needs anyway. Both emit the `BAM_TO_FA` contract (`<id>.{hq,singletons}.fasta.gz` + `meta.singleton`), so alignment onward is untouched. cd-hit: one-member clusters are singletons and the other representatives are hq (`-c cdhit_identity -n 10 -d 0`, `.clstr.gz` published). RATTLE: `cluster --iso` → `correct -r 1` → `polish`; `transcriptome.fq` is hq and `uncorrected.fq` (clusters of one read) is singletons. Empty inputs are dropped before RATTLE, which segfaults on them. PBINDEX is skipped for `cluster` BAMs bound for the new engines. Default `isoseq` behaviour is unchanged. Covered by new `test_isoseq.py` cases and verified against the real containers on synthetic reads with planted cluster sizes.
- `isoseq_cluster2_mode` is renamed `cluster_mode`, since it now drives every engine and `flnc` pooling. A leftover `isoseq_cluster2_mode` fails validation instead of silently falling back to the `both` default.
- `FXSPLIT` no longer swallows every failure: only its benign empty-input exit code (1, "No FASTA records found", no chunks — probed against the real image) is ignored. Anything else follows the global ladder — retry resource codes (including 137 OOM-kill, with memory scaling per attempt), fail fast otherwise — so a killed split aborts the run instead of silently starving alignment, segmentation, and collapse of chunks. Covered by a new `test_fxsplit.py` harness (exit 1 completes chunkless, exit 137 fails loudly), run in CI.
- `flnc` entrypoint honors `cluster_mode` (renamed from `isoseq_cluster2_mode` below): `multi_sample`/`both` concatenate all input reads per hq/singleton class into one pooled sample named `global_prefix` (new `FASTX_CONCAT` module + `POOL_READS` subworkflow), so adapter removal, polyA segmentation, and twin-collapsing see all samples together per chromosome instead of running per sample. `multi_sample` flnc outputs move from per-sample names to `global_prefix`; do not name an input sample `global_prefix`.

## [v2.0.27] - 2026-10-02

This release rolls up everything since v2.0.26 (PRs #43-#47 plus direct fixes) alongside new work on this branch: it fixes a dead second pass (fragment detection never ran for any current entrypoint), adds `refine` and `cluster` restart checkpoints, wires the iso-classify intron track end to end, and implements xORF database merging through `custom_database`. The trackDb template is now generated inline instead of read from disk, the genome-browser upload wiring is corrected, resume after alignment works again, and a full end-to-end CI suite guards the pipeline. CI gold has to be regenerated (see below).

### New checkpoint: `--entrypoint refine` (PR #47)

- Runs restart from LIMA output (`02_LIMA`) with `--entrypoint refine`: every `*.<5p>--<3p>.bam` is loaded with its `.pbi` (missing indexes rebuilt via `PBINDEX`), and CCS, Skera, and LIMA are skipped, including for Kinnex libraries.
- Sample identities keep the full primer pair end to end (`movie.IsoSeqX_bc01_5p--IsoSeqX_3p`, `movie.NEB_5p--NEB_Clontech_3p`): the complete `<5p>--<3p>` suffix is parsed without assuming IsoSeqX names, a BAM missing the suffix fails before refinement, and a stable `sample_id` survives alignment into per-chromosome grouping (dots and suffixes preserved). `flnc` FASTA sample ids are normalized so `.hq` + `.singletons` files of one sample group together.
- Clustering reuses the same refined reads for `per_sample` / `multi_sample` / `both`; pooling goes through cluster2's native FOFN input, so the intermediate pooled-BAM merge and its unused imports are gone. `cluster2` handles a single BAM path as well as a list.
- Added the `refine` value to entrypoint validation, preprocessing routing, config comments, and `params.json`, plus README documentation (checkpoint, BAM naming, clustering modes, barcode assignments) and a standalone `test_isoseq.py` Nextflow regression harness (eight scenarios + negative checks) run in CI.

### BREAKING CHANGE: xORF database merging through `custom_database` (PR #46)

- `xorf_custom_database` now accepts staged FASTA (`.fa` / `.fasta`, optionally `.gz`), which is merged with the default SwissProt sequences and reindexed: new `FASTA_MERGE` + `DIAMOND_MAKEDB` modules, new `xorf_raw_database` parameter (default UniProt SwissProt FASTA), and a reworked preprocessing branch (`.dmnd` / `.dmnd.gz` replace the default database outright; anything else errors with the accepted formats). `xorf` submodule bumped accordingly.
- `workflows/ark.nf` passes the new `xorf_raw_database` through, along with the xORF v0.0.41 inputs (`run_only`, `database_versions`).

### End-to-end CI test suite (PR #45)

- Added a containerised e2e suite: `params.subreads.json` / `params.flnc.json` fixtures, `chr19` fixtures and `flnc` + `subreads` golden trees, and `compare.py` (line-count comparison since `R`-number ids are order-dependent, diffs shown in CI, PBCCS QC noise and `.sam` files ignored). Follow-ups minimized storage (dropped weights, gzipped chr), bumped the Nextflow version, and corrected fixture paths.

### Tool-flag compatibility and resource fixes

- BREAKING: `iso-classify` v0.0.13 renamed `--isoseq` / `--toga` to `--input` / `--reference` (old flags removed in isotools v0.0.42); the intron module passes the new flags, so the container must be rebuilt from isotools >= v0.0.42.
- `bed2gtf` v2 renamed `-i` / `--input` to `-b` / `--input` (`-i` is now `--isoforms`): the module passes `-b`, fixing `invalid isoforms row` failures on 12-column BED input.
- Skera/Kinnex: MAS adapter primer link updated to the Kinnex-full-length-RNA path; `PBSKERA_SPLIT` moved to the new `process_extreme` label (16 CPUs, 32 GB, 72 h) and the never-produced `non_passing.pbi` output was dropped.
- `ISOSEQ_CLUSTER2` moved to `process_extreme`; `ISOTOOLS_TRUNCATION_DETECTOR` gains `-O cds`.

### Docs and assets (PRs #43-#44)

- Updated Hiller Lab logo, then theme-dependent dark/light logos plus the illustrated changelog header.

### Fragment detection fix and post-alignment merge removal

- `ISOTOOLS_FIND_FRAGMENTS` was never submitted: it consumed `ch_pooled_reads`, which was only filled for the legacy `isoseq` / `map` entrypoints removed in v2.0.25, so the process and the three processes behind it (`ARK_ALIGN_FRAGMENTS`, `SAMTOOLS_BAM_FRAGMENTS`, `ISOTOOLS_SEGMENT_POLYA_FRAGMENTS`) ran zero times with `ark` + second pass on, for every entrypoint.
- Each chunk BAM now gets back the chunk FASTA it was aligned from. `FXSPLIT` already emits those, and `sample_id`, `singleton`, and `chunk` survive every `meta.clone()`, so they form the key: chunk FASTAs are grouped by `[ sample_id, singleton, chunk ]` (`groupTuple`) and combined with the aligned BAMs (`combine(by: 0)`, `multiMap`), so colliding inputs (e.g. `X.fasta.gz` and `X.hq.fasta.gz` under `flnc`) both reach `--reads` instead of being mispaired. One code path covers every entrypoint with cigar extension on or off, and versions are mixed on both paths.
- The post-alignment `map` + `multi_sample` / `both` merge branch is deleted, not revived: `flnc` stays per-sample, as since v2.0.25. With it go the seven `SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_*` includes, their `withName` blocks, the deleted `src/modules/custom/samtools/merge/main.nf`, and the now-readerless `prefix`, `cluster_mode`, and `entrypoint` inputs of `SPLIT_ALIGN_CLEAN_CHUNKS`. `isoseq_cluster2_mode` therefore only affects the Iso-Seq entrypoints, where pooling happens before clustering. Reasons: the merged BAM was segmented with a single `--singleton` setting (singleton reads lost their `SG` tag), `both` / `multi_sample` are the defaults so ordinary `flnc` runs would change sample names, the merged BAM went through cigar extension and segmentation unchunked, and separately clustered files share read names (`transcript/N`) that `iso-align` groups by.
- `ISOTOOLS_FIND_FRAGMENTS` now publishes to `06_ARK_ALIGN/FRAGMENTS` next to the other ark fragment outputs (it never published anything, so nothing moves) with `ext.prefix = <id>.<chunk>[.singleton]`, since `meta.id` is shared by every chunk of a sample and outputs would otherwise overwrite each other.

### New checkpoint: `--entrypoint cluster`

- Restarts from `03_ISOSEQ_REFINE`: reads every `*_flnc.bam`, skips CCS, Skera, LIMA, and refine, and goes straight to `ISOSEQ_CLUSTER2` under the same `isoseq_cluster2_mode` rules (`per_sample` / `multi_sample` / `both`).
- The sample id is the BAM name without `_flnc.bam`, i.e. the id the full pipeline gives that sample (`<movie>.<5p>--<3p>`); a BAM not ending in `_flnc.bam` is an error.
- `--global_primers` is not required and never opened for this entrypoint. Added to the validator in `src/main.nf`, the preprocessing Iso-Seq branch, `nextflow.config`, `params.json`, and the README alongside the `refine` checkpoint.

### Intron BED4 track end to end

- `iso-classify intron` now runs with `--intron-track`. The module emits the track as a second output, renames the upstream comma typo (`<prefix>.introns_track,bed` -> `<prefix>.introns_track.bed`, fixed upstream in isotools), and drops the file when empty.
- `PREPOLISH` emits `intron_track`; `workflows/ark.nf` groups it by `meta.name` through `JOIN_INTRONS` (sorted BED in `12_POLISH/BED/<name>.introns.bed`) and `BEDTOBIGBED_INTRONS` without autosql (BED4), mixed into `ch_additional_bbs` so it lands in `12_POLISH/BB` with the fusion and NMD tracks. `POLISH` additionally emits `rt`, and `JOIN_INTRONS` was added to the `12_POLISH/BED` publish pattern in `nextflow.config`.

### TrackDb generated, upload wiring fixed

- `TRACKDB` no longer reads `assets/as/track.as` (the path was wrong twice: single-quoted `${projectDir}` plus one `..` too many). The stanza is now a heredoc inside the module like the autosql schemas, so the `sed` step and the `schema` input are gone. `assets/as/` is deleted: `track.as` moved into the module, `base.as` / `schema.as` were stale copies that still said `isopipe`. Version reporting switches from `sed` to `bash`.
- Fixed while moving the template, because each `bigDataUrl` must equal the file `RSYNC_SSH` uploads: `retention.bb` -> `retentions.bb`, `truncation.bb` -> `truncations.bb`, `orphans.bb` -> `pass.scraps.bb`, `duplicates.bb` -> `pass.duplicates.bb`. The `{BIGWIG_TRACK}` placeholder was never substituted and reached the output literally; it is now `<track>_bigwig`. The `spliceai.acceptor.reverse` subtrack pointed at the `aparent.forward` file.
- `workflows/ark.nf` now loads `scraps` instead of the never-existing `ISOTOOLS_POLISH.out.orphans` (which broke `load_track = true`), feeds `LOAD_NMD_TRACK` the per-sample bigBed instead of the per-chromosome BEDs, and adds `LOAD_RT_TRACK` for the `rt` subtrack the trackDb always listed but never uploaded.

### Resume and correctness fixes

- Resume after alignment produced nothing with the default `minimap2_align_keep_sam = false`: `SAMTOOLS_BAM` deletes the SAM, `MINIMAP2_ALIGN` declared it `optional: true`, and on `-resume` the cached aligner emitted nothing while everything downstream was silently skipped (exit 0, no results). The `sam` output is required again, so a cached task with a missing SAM is re-executed. Output declarations are not part of the task hash, so nothing is invalidated.
- `SAMTOOLS_BAM` ran `samtools index -@ {task.cpus}` (missing `$`), i.e. single-threaded; fixed to `-@ ${task.cpus}`. This changes the task hash of every BAM conversion: runs resumed across this change redo BAM conversion and everything after it.
- `COLLAPSE` built its prefix with `".${meta.chr}" ?: ''` (always truthy), producing `*.pass.null.collapsed.bed`; fixed to `<sample>.pass.collapsed.bed`. Gold path changes.
- `splicing/main.nf` called `GUNZIP_MINISPLICE` without including it (included now); the `def x = f(<take input>)` declarations that Nextflow 24.10 refuses to compile are gone, so the pipeline compiles on 24.10.5 and 25.10.2; `SPLICEAI_RUN`'s bare directory is wrapped as `[ meta, dir ]` for `SPLICEAI_DERIVE` and emitted as `bigwigs`, so intron classification also sees pipeline-computed bigWigs.
- `SPLIT_ALIGN_CLEAN_CHUNKS` and `PREPROCESSING` ignored their `aligner` take-input and read `params.aligner`; they use the input now, and the unknown-aligner messages list `ark`.
- `PREPROCESSING` overwrote `ch_versions` on the custom-FASTA database path instead of mixing into it.
- Removed with no caller left: `ISOTOOLS_FIND_FRAGMENTS_ULTRA` / `_DESALT` `withName` blocks, the `samtools/merge` module, the two `GAWK_JOIN_BEDGRAPH_*` includes in `prepolish`, and the `ultra_do_second_pass` parameter in `params.json`.

### CI and test harness

- `ci.yml` now triggers on `master` (the default branch), not `main`.
- Extended `assets/ci/test_isoseq.py` with stub `iso-cigar` / `iso-align` (`iso-align` asserts `--bam` and every `--reads` exist) and distinct genome/annotation names: `refine` per/multi/both, `flnc` two-sample, and `cluster` (no `lima`/`refine`, same `cluster2` ids as `refine`, non-`_flnc` rejected) with the second pass on and cigar extension on/off, asserting `iso-align` runs once per chunk BAM with exactly that chunk's FASTA and that `samtools merge` is never called; malformed LIMA/refined names fail before any tool runs; the CLI validator covers `refine` needing primers, `cluster` not needing them, and unknown entrypoints.

### Chores

- Bumped the pipeline version to 2.0.27 in the Nextflow manifest.
- Entry checkpoints are now `subreads, ccs, refine, cluster, flnc` in validation, config comments, and `params.json` (`refine` via PR #47, `cluster` on this branch); the full aligner set (`ark, mm2, ultra, desalt, pbmm2, flair`) is listed in the unknown-aligner errors.
- Added `INFO`/`WARN` channel-shape comments across `genome`, `spliceai`, `splicing`, `prepolish`, `polish`, `split_align`, `track`, and `workflows/ark.nf` with no logic change.
- Gold has to be regenerated from a real containerised run: new `06_ARK_ALIGN/FRAGMENTS/*.report.tsv` (plus fragment outputs where reads qualify), `11_PREPOLISH/CLASSIFY/*.introns_track.bed`, `12_POLISH/BED/*.introns.bed`, possible downstream line-count shifts from fragment realignment, and the `*.pass.null.collapsed.bed` -> `*.pass.collapsed.bed` rename. The `subreads` gold already predates the primer-pair sample ids. CI will be red until gold is regenerated.

## [v2.0.26] - 2026-07-31

This release adds native support for Kinnex/MAS-Seq libraries through skera-based demultiplexing, and reworks the primer-removal and refinement segment of the ISOSEQ subworkflow so that per-primer-pair reads are tracked individually instead of being merged back into a single pool. The new `skera_is_kinnex_library` parameter activates a `PBSKERA_SPLIT` step between CCS generation and LIMA, with the MAS-Seq adapter primer set downloaded automatically at runtime when the option is enabled. On the refinement side, the LIMA output channel was restructured from the ground up: the multi-sample merge that previously collapsed per-barcode BAMs back into one file has been removed entirely, and each barcode-specific BAM is now paired with its `.pbi` index and handed to ISOSEQ_REFINE as an independent tuple. Several entrypoint dispatch and versioning-channel corrections carried over from the v2.0.25 maintenance cycle are included as well, along with a LIMA container update, a dedicated high-core resource label, and a number of output-emission fixes in the new skera module.

### Kinnex library demultiplexing

- Added the `skera_is_kinnex_library` parameter (default: `false`). When enabled, CCS BAMs are demultiplexed with skera through the new custom `PBSKERA_SPLIT` module (biocontainers/pbskera 1.4.0) before reaching LIMA for primer removal. The parameter is declared in `nextflow.config` and `params.json` and propagated through the preprocessing and ISOSEQ subworkflows.
- Added the `skera_kinnex_primers` parameter, which defaults to the MAS-Seq Adapter v3 8-primer FASTA hosted on the PacBio Cloud downloads server. When demultiplexing is active, the primer set is fetched at runtime through a dedicated WGET step instead of requiring a manual download and a local path.
- `PBSKERA_SPLIT` now emits the full set of skera outputs with correct filename prefixes: the demultiplexed `*.skera.bam` together with its `.pbi` index, the `*.non_passing.bam` file and index, `*.found_adapters.csv.gz`, `*.summary.csv`, `*.summary.json`, `*.ligations.csv`, and `*.read_lengths.csv`. Outputs are published to a dedicated `02_PBSKERA_SPLIT` directory.
- When `skera_is_kinnex_library` is enabled, LIMA is invoked with `--overwrite-biosample-names` so the biosample names produced during demultiplexing are carried through primer removal instead of being reset.
- Fixed a stray whitespace character in the pbskera container string that caused the container identifier to resolve incorrectly.

### Primer removal and refinement rework

- LIMA now expects the CCS BAM paired with its `.pbi` index as a three-element input tuple, and its container was bumped from lima 2.9.0 to 26.2.1. Both the merged CCS output from PBMERGE and the demultiplexed output from PBSKERA_SPLIT are joined with their corresponding index files before being passed to LIMA.
- Removed the merge step that previously recombined multi-barcode LIMA outputs: the single/merge branching logic and the `PBMERGE_MULTI_LIMA` invocation in the ISOSEQ subworkflow are gone. LIMA output is instead flattened into one channel item per BAM, keyed as `sample::stem`.
- Each per-primer-pair BAM is now identified by its barcode: the filename stem is parsed for the `IsoSeqX_bcNN_5p--IsoSeqX_3p` pattern and the extracted barcode is stored in the metadata alongside a new `parent_id` field. The sample id becomes `sample.barcode`, preserving per-barcode identity through the downstream steps.
- Each barcode-specific BAM is joined with its matching `.pbi` and passed to `ISOSEQ_REFINE` as an individual tuple; the refine process input was updated to accept the index alongside the BAM. Reads are therefore refined per primer pair rather than as a merged pool, and the polyA filtering operates with per-barcode context.
- LIMA's resource label was changed from `process_low` to the new `process_high_core` label, and `--log-level INFO` was added to its arguments. The label was introduced at 32 CPUs, 32 GB, and 16 hours per attempt, then adjusted down to 16 CPUs with the same memory and time limits to better match actual usage.

### Entrypoint dispatch and channel wiring fixes

- Corrected the entrypoint dispatch inside the ISOSEQ subworkflow to match the v2.0.25 entrypoint model: the subworkflow now receives `subreads` or `ccs` (with `flnc` unreachable). `subreads` runs the PBCCS chunking and merge path, while `ccs` passes BAMs straight to primer removal. Previously the subworkflow expected `ccs`/`flnc` values it could no longer receive from preprocessing.
- Moved the PBCCS, PBINDEX, and PBMERGE version-channel mixing inside the `subreads` branch. Those processes only execute for that entrypoint, so mixing their `versions.yml` outputs unconditionally left the version channel in an inconsistent state for `ccs` runs.
- Gave the LIMA step a dedicated output channel (`ch_lima_out_bams`) instead of reusing the input channel, preventing a channel reassignment that could break the workflow DAG depending on entrypoint or demultiplexing configuration. As part of this cleanup the index element was temporarily stripped from the LIMA input tuple, a detail later revised when the `.pbi` became part of LIMA's input contract.

### Chores

- Bumped the pipeline version to 2.0.26 in the Nextflow manifest.
- Corrected the `ark` project name rendering in the README header.

## [v2.0.25] - 2026-07-29

This release replaces the previous two-entrypoint model (`isoseq` / `map`) with three granular modes — `subreads`, `ccs`, and `flnc` — that map directly to the real PacBio processing stages users are working with. The new model gives you control over where your data enters the pipeline: raw subreads, already-called CCS reads, or full-length non-chimeric reads that have already gone through primer removal and clustering. The internal ISOSEQ subworkflow was refactored to branch on these entrypoints, running PBCCS and chunk merging only when needed (i.e., for `subreads` and `ccs`), while `flnc` reads skip straight to LIMA and refinement. The default entrypoint was changed from `isoseq` to `subreads`, and the workflow banner now prints the active entrypoint at launch to make it immediately visible which mode is in use. Additional structural comments were added throughout the main workflow and the ISOSEQ subworkflow to clarify the stage boundaries in the code.

### New entrypoint model

- Replaced the previous `isoseq` and `map` entrypoints with three new options: `subreads`, `ccs`, and `flnc`. Validation in `main.nf` was updated to reject unknown values with a clear error message listing the valid options.
- Changed the default `entrypoint` parameter from `"isoseq"` to `"subreads"` in `nextflow.config` and `params.json`.
- The workflow banner now displays `Entrypoint: ${params.entrypoint}` alongside the input, output, and genome paths, making it easier to confirm which mode is active at a glance.

### ISOSEQ subworkflow refactoring

- The `ISOSEQ` subworkflow now accepts an `entrypoint` input parameter (expected values: `ccs` or `flnc`; `subreads` is handled upstream). A `switch` statement dispatches execution into two branches:
  - **`ccs`**: runs `PBCCS` on chunked BAMs to generate circular consensus sequences, groups chunks by parent sample via `groupTuple`, merges them with `PBMERGE`, and passes the merged BAMs to `LIMA` for primer removal.
  - **`flnc`**: passes the raw input BAMs directly to `LIMA`, skipping CCS generation and merging entirely.
- An explicit error is thrown if an unrecognized entrypoint reaches the subworkflow.

### Preprocessing dispatch update

- The preprocessing subworkflow now checks for `subreads` or `ccs` to route into the `ISOSEQ` branch (the `ISOSEQ` call now also passes the `entrypoint` value), and checks for `flnc` to enter the direct FASTQ-reading branch (formerly the `map` branch).

### Codebase organization

- Added ASCII-delimited section headers throughout `workflows/ark.nf`, `subworkflows/isoseq/main.nf`, and `subworkflows/preprocessing/main.nf` to visually separate stages such as Autosql, Preprocessing, Alignment, ORF calling, NMD calling, Polishing, and Tracking. These are purely cosmetic and have no effect on pipeline logic.

## [v2.0.24] - 2026-07-24

This release consolidates the pipeline's alignment architecture under a single unified subworkflow, renames the project from isopipe to ARK, and introduces FLAIR as a sixth aligner option. The previously separate desalt, ultra, and pbmm2 alignment subworkflows have been merged into a single `split_align` subworkflow that dispatches to the correct backend via a switch statement, significantly reducing code duplication. Additional work includes a new ORF renaming engine for xORF, programmatic AutoSQL schema generation, and several channel-wiring corrections across downstream subworkflows.

### Pipeline rename to ARK

- Renamed the primary workflow from `ISOPIPE` to `ARK` and the corresponding entry file from `workflows/isopipe.nf` to `workflows/ark.nf`. The pipeline description was updated to "A Reference pipeline to annotate euKaryotes at high resolution" and all manifest metadata (name, homePage) now point to the `alejandrogzi/ark` repository.
- Updated `main.nf` to include `ARK` from the new workflow path and to log the Ark banner (version, description, authors, lab) instead of the former isopipe banner.

### Unified alignment architecture

- Removed three dedicated alignment subworkflows — `desalt_align`, `ultra_align`, and `pbmm2_align` — and consolidated all six supported aligners into the single `split_align` subworkflow. The subworkflow now accepts an `aligner` string parameter and uses a switch statement to route chunked reads through the correct alignment module, BAM conversion, and multi-sample merge path.
- Added dedicated index modules for each aligner: `ARK_INDEX` and `FLAIR_INDEX` (both minimap2-based with appropriate `-x` presets), alongside the existing `MINIMAP2_INDEX`, `ULTRA_INDEX`, `DESALT_INDEX`, and `PBMM2_INDEX`. Each aligner now has its own publish directory under `06_<ALIGNER>_INDEX`.
- Each aligner's BAM output now routes through a uniquely named samtools conversion process (`SAMTOOLS_BAM_ARK_ALIGN`, `SAMTOOLS_BAM_MINIMAP2_ALIGN`, `SAMTOOLS_BAM_DESALT_ALIGN`, `SAMTOOLS_BAM_PBMM2_ALIGN`, `SAMTOOLS_BAM_FLAIR_ALIGN`) to prevent channel collisions in multi-sample merging.
- Multi-sample BAM merging was similarly split into per-aligner `SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_*` processes, ensuring correct version tracking when multiple aligners are exercised across pipeline invocations.

### New aligner: FLAIR

- Added **FLAIR** as a selectable aligner (`aligner = "flair"`). The integration includes a `FLAIR_INDEX` module that builds the minimap2 index with the `-x splice` preset and a `FLAIR_ALIGN` process aliased from `MINIMAP2_ALIGN` that runs with the standard splice preset. The `split_align` subworkflow routes FLAIR-aligned output through a dedicated `SAMTOOLS_BAM_FLAIR_ALIGN` conversion step.

### Ark as default aligner

- Changed the default `aligner` parameter from `"minimap2"` to `"ark"`. The `"mm2"` aligner runs minimap2 with the PacBio Kinnex standard preset (`-uf -ax splice:hq`), replacing the previous verbose flag set (`-a -c --eqx -uf -C5 -G ... -ax splice:hq --secondary=... -cs ... --junc-bonus ... --junc-pen ...`) that now runs through `"ark"`.
- The second-pass re-alignment (cigar extension + fragment detection) is now gated exclusively behind `params.aligner == 'ark'`, since other aligners handle fragment resolution differently or do not benefit from the additional pass.

### ORF renaming engine

- Added `rename_predictions.py` (v0.0.4), a standalone Python script that rewrites ORF prediction identifiers in BED12/TSV files with human-readable, information-rich names. Supported features include hash-based rebasing of root IDs (`--rebase`), ORF score appending (`--append-orf-score`), protein name appending (`--append-protein-name`), custom prefix injection (`--custom-prefix`), and full deactivation (`--deactivate`).
- Added five new xORF parameters: `xorf_rename_deactivate` (default: `false`), `xorf_rename_rebase` (default: `false`), `xorf_rename_append_orf_score` (default: `true`), `xorf_rename_append_protein_name` (default: `true`), and `xorf_rename_custom_prefix` (default: `null`). These are propagated through both `XORF_PREDICT_ORFS` and `XORF_PREDICT_FUSION_ORFS` process invocations.

### AutoSQL schema generation

- Added `AUTOSQL_BASE` and `AUTOSQL_SCHEMA` Nextflow processes that programmatically generate AutoSQL `.as` schema files at runtime, replacing the previous reliance on static files from `assets/as/`. `AUTOSQL_BASE` emits a minimal BED12-compatible schema, while `AUTOSQL_SCHEMA` emits the extended schema including read-status, metadata, collapsed-reads, and ORF metadata columns. Both are imported in the main ARK workflow and their outputs are wired to the `autosql` and `schema` channels used by downstream bed-to-bigbed conversion.

### Parameter and validation changes

- Renamed the `minimap2` aligner option to `mm2` throughout the parameter schema, documentation, and validation logic. The valid aligner set is now: `mm2`, `ultra`, `desalt`, `pbmm2`, `flair`, `ark`.
- Removed the `ultra_do_second_pass` parameter entirely. The Ultra two-pass logic from v2.0.22 is superseded by the unified second-pass gating in `split_align`, which only activates for the `ark` aligner.
- Improved validation error messages with explicit `ERROR:` prefixes and structured multi-line formatting via `stripIndent()`.

### Channel wiring fixes

- Corrected the annotation channel format in `PREPOLISH` and `POLISH` subworkflows: the annotation is now expected as `[ val(meta), [ annotation ] ]` from preprocessing rather than being re-mapped with `[id:annotation.baseName]` at the subworkflow boundary. This prevents inconsistent channel shapes when the annotation flows through multi-sample paths.
- The `POLISH` subworkflow now receives the `autosql` channel as an explicit input parameter rather than constructing it from a static asset path, ensuring the AutoSQL schema stays in sync with the dynamically generated output.

### Chores

- Added `.gitattributes` to the repository root.
- Updated `params.json` schema documentation to reflect the new aligner list and the `xorf_rename_deactivate` parameter.
- Removed dedicated `ext.args` blocks for `MINIMAP2_ALIGN` fragment-detection processes, replacing them with the simplified Ark preset configuration.

## [v2.0.23] - 2026-07-18

This release expands the pipeline's alignment options with two new backends — deSALT and pbmm2 — giving users four aligners to choose from depending on their sequencing platform and sensitivity requirements. The minimap2 branch also receives several correctness fixes, particularly around second-pass gating and index preset configuration, along with minor resource tuning and a new xORF parameter.

### New aligner backends

- Added **deSALT** as a selectable aligner (`aligner = "desalt"`), a splice-aware long-read aligner well-suited for isoform discovery. The integration includes a new `DESALT_INDEX` module for building the deSALT genome index, a `DESALT_ALIGN` module for alignment with configurable annotation guidance via `desalt_use_annotation`, and a dedicated `desalt_align` subworkflow that handles chunking, alignment, BAM conversion, multi-sample merging, adapter removal, cigar extension, polyA segmentation, twin collapse, and fusion detection.
- Added **pbmm2** as a selectable aligner (`aligner = "pbmm2"`), PacBio's official minimap2 wrapper optimized for CCS reads. The integration includes a `PBMM2_INDEX` module that builds the index using the `splice` preset, a `PBMM2_ALIGN` module using the `ISOSEQ` preset, and a dedicated `pbmm2_align` subworkflow covering the full post-alignment processing pipeline.
- Updated the `aligner` parameter documentation and schema to reflect all four supported options: `minimap2`, `ultra`, `desalt`, and `pbmm2`.
- Extended the preprocessing subworkflow with conditional index-building branches for deSALT and pbmm2, accepting new `desalt_index` and `pbmm2_index` parameters to support pre-built index paths.

### Minimap2 indexing and second-pass fixes

- Corrected the `MINIMAP2_INDEX` process configuration by adding `ext.args` with the appropriate `-x` splice alignment preset, ensuring the minimap2 index is built with the correct parameters for transcript-aware alignment.
- Fixed `ext.use_junc_bed` assignment for both `MINIMAP2_ALIGN` and `MINIMAP2_ALIGN_FRAGMENTS` to properly reflect `!params.minimap2_align_use_splice_scores`, preventing the junction BED from being applied in scenarios where splice scores should take precedence.
- Gated the entire fragment detection and second-pass re-alignment logic in the `split_align` subworkflow behind the new `minimap2_align_do_second_pass` parameter (default: `true`). When disabled, the pipeline skips fragment detection and minimap2 re-alignment entirely, proceeding directly to segmentation with the primary BAM. This mirrors the behavior already present in the Ultra backend and avoids unnecessary computation when the second pass is not needed.

### xORF improvements

- Added the `xorf_skip_netstart` parameter (default: `true`) to allow users to bypass NetStart-based start codon predictions during ORF calling, streamlining results when only canonical start codons are of interest.

### Resource tuning

- Reduced the `process_medium_high_memory` resource label from 6 CPUs and 96 GB to 4 CPUs and 36 GB per task attempt, providing a more cost-effective baseline for medium-memory workloads without compromising stability.

### Chores

- Bumped `isotools` and `xorf` submodules to their latest versions.

## [v2.0.22] - 2026-06-19

This release introduces a two-pass alignment strategy for the Ultra backend that significantly improves sensitivity for fragmented reads, along with finer-grained control over minimap2 junction scoring during splice site annotation. Several channel wiring bugs introduced in v2.0.21 have also been corrected.

### Ultra two-pass alignment

- Implemented a second pass in the Ultra alignment branch: after the initial Ultra alignment, reads undergo optional cigar extension and fragment detection via isotools, and the resulting fragments are re-aligned using minimap2 for higher-resolution junction placement. This is controlled by the new `ultra_do_second_pass` parameter (default: `true`).
- Added the `ultra_max_intron_size` parameter (default: 300 kb) to configure the maximum intron length passed to the Ultra aligner.
- Added new process modules for the ultra second pass: `MINIMAP2_ALIGN_FRAGMENTS_ULTRA`, `SAMTOOLS_BAM_FRAGMENTS_ULTRA`, and `ISOTOOLS_FIND_FRAGMENTS_ULTRA`, each with dedicated process configuration in `nextflow.config`.
- Updated the preprocessing subworkflow to conditionally build a minimap2 index when Ultra is selected and `ultra_do_second_pass` is enabled, ensuring the second pass has the required index available.
- Changed the Ultra merged BAM output directory from `06_ULTRA_ALIGN/MERGED` to `06_ULTRA_ALIGN/ULTRA_MERGED` to avoid ambiguity with other merge paths, and added `FRAGMENTS` publish directories for fragment-level SAM and BAM outputs.
- The Ultra second pass fragments now flow through `ISOTOOLS_SEGMENT_POLYA_FRAGMENTS` and are mixed with the primary segmented output for downstream processing.

### Minimap2 scoring granularity

- The `--junc-bed` flag in `MINIMAP2_ALIGN` is now conditionally applied based on the new `ext.use_junc_bed` configuration value, allowing the pipeline to omit the junction bed when splice scores are available and should take precedence.
- The `SPLICEAI_DERIVE` process arguments (`--include-ss-from-regions`, `--position-for-ss-regions`, `--bonus`) are now only applied when both `minimap2_align_use_junc_bed` and `minimap2_align_use_splice_scores` are enabled, preventing redundant or conflicting scoring modes.

### Bug fixes

- Corrected the annotation channel format in the `SPLICING` subworkflow to prevent the annotation path from being incorrectly wrapped in a meta tuple when passed to `SPLICEAI_DERIVE`.
- Fixed `BED2GTF` invocation in the preprocessing subworkflow to supply a proper meta map instead of relying on the bed file's base name directly.
- Ensured `global_output_dir` is properly referenced in the Slurm scheduler script (`assets/sh/do_isopipe.sh`).

## [v2.0.21] - 2026-06-18

This release introduces a major architectural shift with the integration of the Ultra aligner alongside the existing minimap2 pipeline, giving users the flexibility to choose between alignment backends. It also brings significant improvements to splicing score annotation, Veredict classification, and HPC job orchestration.

### Ultra aligner integration
- Added a dedicated ultra aligner branch with new subworkflow `ultra_align/main.nf` and nf-core modules for `ultra/index` and `ultra/align`. The pipeline now branches at alignment-time: if `aligner = "minimap2"` the existing minimap2 path is used; if `aligner = "ultra"` the new Ultra path is taken, including Ultra-specific index preparation, alignment, BAM merging, and adapter removal.
- Updated `main.nf`, `nextflow.config`, and preprocessing subworkflows to conditionally dispatch channels based on the selected aligner, and added `ultra_index` and `ultra_use_annotation` parameters.
- Fixed merging step to be Ultra-specific and corrected input Ultra index channel wiring.
- Updated collateral modules to align with the new aligner branching logic.

### Splicing score annotation
- Bumped `splicing-rs` to v0.0.6, adding `--bonus` and `--bonus-score` flags to reward annotation-derived spliced sites during spliceai derive, improving the quality of junction-level scoring for annotated transcripts.
- Propagated the `--bonus` flag through the spliceai derive module.

### Veredict classification
- Veredict.py updated to v0.0.8, implementing a new `ARTIFACT` output category that allows the pipeline to explicitly flag and separate technical artifacts from genuine biological signals in the polishing step.

### xloci improvements
- Added `--unmask` option to xloci intron extraction, providing finer control over repeat masking during intron boundary delineation.

### Infrastructure and HPC
- Added a new HPC matrix scheduler script (`assets/sh/do_isopipe.sh`) to streamline job array submission and resource allocation across Slurm-based clusters.
- Fixed annotation channel tuple propagation across downstream processes to prevent channel misalignment in multi-sample runs.

## [v2.0.20] - 2026-05-28

### Breaking Changes
- v2.0.20 -> conditional cigar extension step!
- v2.0.20 -> splicing v0.0.5, include SS from --regions in final derived score tsv
- fix collapse impl
- static -> add additional collapsing step to remove extremely deep duplicates after segmenting
- V2.0.19 static -> iso-align impl + pooled reads to find fragments
- v2.0.19 -> re-align implementation of fragments

### Features
- SS from regions + SS CDS on default in spliceAi derive

### Fixes
- remove --collapse-mode from collapse mod
- ensure collapses fills up new channel
- avoid grep code 1 on non-duplicates vs duplicates
- update schema + include RT in output

### Chores
- bump minimap2 image to latest version
- pre-release iso-align + fragment logic
- bump submodules

## [v2.0.19] - 2026-04-24

### Breaking Changes
- v2.0.19 -> impl collpase for passes
- v2.0.19 -> fix veredict RT masking
- v2.0.19 -> include collapse in schemas

### Features
- v2.0.18 static -> replace gawk join begraph with bigwig merge and single bg to bw
- static -> add global_prefix to avoid colliding files in scratch

### Chores
- drop ucsc binaries
- bg baseName for bw + bigwigmerge publish

## [v2.0.18] - 2026-04-22

### Breaking Changes
- v2.0.18 -> add bigWig to track + increase classify core count + add process_high_low_memory_fast
- v2.0.18 -> drop validate_id return error, keep as warn

### Features
- add genepred lint in preprocessing

### Chores
- modify process for adapter

## [v2.0.17] - 2026-04-22

### Breaking Changes
- v2.0.17 -> relax validation id allowing asymmetry + warnings
- v2.0.17 -> add isotools adapter

## [v2.0.16] - 2026-04-20

### Breaking Changes
- v2.0.16 -> fix APARENT padding lower limit to 10bp

### Chores
- bump isotools intron process
- rm scrap

## [v2.0.15] - 2026-04-19

### Breaking Changes
- v2.0.15 -> staging options in config + split align merging on file name stripped of cigar extension
- v2.0.15 -> aparent interval logic re-impl
- isotools cigar impl -> duplicate isotools segment [add bai to input channel]
- static v2.0.14 -> add splicing scores for orphan + port logic into preprocess + splicing
- static v2.0.14 -> strip WGET duplicated processes + allow local path avoiding download

### Fixes
- clone meta for extended in isotools-cigar
- update WGET calls for output file

### Chores
- bump submodules
- bump xorf
- copyright + drop extract
- drop minimap correction -> replaced by iso-cigar
- bump isotools
- add FILL marker to xorf_selenocysteine
- Change condition to clean up temporary BAM files

## [v2.0.14] - 2026-04-10

### Breaking Changes
- v2.0.14 -> splicing v0.0.4; automate bigwig naming tokens [plus/minus + donor/acceptor]
- v2.0.14 -> isotools orphan impl + bump updates + loading track statements
- v2.0.14 -> trackDb impl + nf mod
- BREAKING CHANGE: v2.0.14 -> spliceAi impl [chunk, prediction, publish] + side tools bigwigmerge, wigtobigwig + rust/py impl + docker img + sync subworkflow inputs
- v2.0.14 -> spliceai pipeline impl
- v2.0.14 -> spliceai impl

### Features
- add gxf2bed mod
- feat/fix: add detach + join fusions/nmd -> bb + publish

### Fixes
- add procps to spliceai img

### Chores
- bump xorf
- add spliceai img
- fmt + correct output channel names
- fmt + img version

## [v2.0.13] - 2026-04-09

### Notes
- WARN: revert versioning, typo -> v2.0.13

## [v2.0.12] - 2026-04-06

### Breaking Changes
- v2.0.12 -> track upload impl

### Features
- add ssh impl + include rsync/ssh in rs img
- --autoseql + bed12 base schema
- publish joined bigbeds

### Fixes
- remove empty outputs + add sorting step on bed extension

## [v2.0.11] - 2026-04-05

### Breaking Changes
- v2.0.11 veredict nf mod
- v2.0.11 -> isotools + generic join + bedtobigbed impl
- v2.0.11 -> veredict impl v0.0.3 + new bb schema

## [v2.0.10] - 2026-04-04

### Breaking Changes
- v2.0.10 -> v0.0.6 aparent, fix threshold type to float
- v2.0.10 -> merge bed + introns as single input for polishing step

## [v2.0.9] - 2026-04-04

### Breaking Changes
- v2.0.9 -> isotools polish + re-factoring xloci, bigtools, aparent (--threshold)

## [v2.0.8] - 2026-04-02

### Breaking Changes
- v2.0.8 -> aparent v0.0.4, fix ghost outputs

### Fixes
- force single channel for aparent joins
- cover lima multi bam output

### Chores
- upgrade isoseq cluster2 reqs

## [v2.0.7] - 2026-03-26

### Breaking Changes
- v2.0.7 -> bump aparent, fix printing statement

## [v2.0.6] - 2026-03-26

### Breaking Changes
- v2.0.6 -> aparent reverse sorted output + bg join + bigwig conversion impl
- v2.0.6 -> bg/bw + ghost outputs fix + conversion impl

### Chores
- emit chromsizes

## [v2.0.5] - 2026-03-26

### Breaking Changes
- v2.0.5 -> aparent predict impl!
- v2.0.5 -> prepolish placeholders from other subworkflows

### Features
- aparent mod inclusion [chunk + predict]
- wget mod

### Fixes
- include pandas

### Chores
- fix CLI tool name
- intronIC fixed args
- update args in xloci intron
- update entry command
- bump submodule
- bump xorf

## [v2.0.4] - 2026-03-25

### Breaking Changes
- v2.0.4 -> prepolish subworkflow impl
- update img; include aparent, publish py img

### Features
- aparent ini + docker py
- add intronIC mod
- add isotools classify intron mod
- add xloci intron mod
- multi-sample to BAM before merging -> change SAMTOOLS_MERGE input

### Chores
- modify align process specs
- fmt
- samtools merge as process_low_long
- watermark

## [v2.0.3] - 2026-03-16

### Breaking Changes
- v2.0.3 -> splicing v0.0.3; adjusting coords to match minisplice coords, avoids spurious alignment shift

### Fixes
- derive junc-bed as value Channel instead of direct config path + config publishes minimap2 sam/bam

## [v2.0.2] - 2026-03-15

### Breaking Changes
- v2.0.2 -> samtools merge in mapping branch

### Features
- open port for minisplice preloaded scores

### Fixes
- match minimap2 index emitted channel preloaded/raw; explicit emptiness with null splice_scores

### Chores
- fix entrypoint docs + update latest img in spliceai derive

## [v2.0.1] - 2026-03-14

### Breaking Changes
- bump v2.0.1 -> isox-rs/splicing v0.0.2

### Fixes
- correct building path
- update bigwig arg

### Chores
- add docs
- change version to first release of isox-rs

## [v2.0.0] - 2026-03-13

### Breaking Changes
- v2.0.0 -> naming, re-factoring + CI publishing workflow
- v2.0.0 -> splicing scores prediction through spliceai/minisplice + compatibility
- v2.0.0 -> split align/bam process; bumps to latest minimap2 version + --spsc
- v2.0.0 codebase upload
- v2.0.0 init!

### Features
- updated assets
- isolate all params in main.nf

### Fixes
- drop tmps

### Chores
- bump time, bytes, slab dependencies

## [v1.0.0] - 2025-12-20

### Fixes
- missing comma leading to reading column names out of frame
