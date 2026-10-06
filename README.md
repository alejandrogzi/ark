<p align="center">
  <p align="center">
    <img width=100 align="center" src="./assets/figures/logo.png" >
  </p>

<p align="center">
  <picture>
    <source
      media="(prefers-color-scheme: dark)"
      srcset="./assets/figures/hillerlab-dark.png"
    >
    <source
      media="(prefers-color-scheme: light)"
      srcset="./assets/figures/hillerlab-light.png"
    >
    <img
      width="200"
      alt="Hiller Lab"
      src="./assets/figures/hillerlab-light.png"
    >
  </picture>
</p>

  <span>
    <h1 align="center">
        ark
    </h1>
  </span>

  <p align="center">
    <a href="https://github.com/hillerlab/ark" reference="_blank">
      <img alt="GitHub License" src="https://img.shields.io/github/license/hillerlab/ark?color=blue">
    </a>
  </p>

  <p align="center">
    <samp>
        <span> A Reference pipeline to annotate euKaryotes at high resolution </span>
        <br>
        <span> The Hiller Lab at the Senckenberg Research Institute </span>
        <br>
        <br>
        <a href="https://github.com/alejandrogzi/ark/blob/master/assets/docs/usage.md">usage</a> .
        <a href="https://github.com/hillerlab/ark/blob/main/assets/pipeline/ark.mermaid">pipeline</a> .
        <a href="https://hillerlab.com/">us</a> 
    </samp>
  </p>

</p>

---

## Usage

> [!NOTE]
> Requirements: Nextflow ≥ 25.04.6, Docker or Apptainer, Java.

```bash
git clone https://github.com/hillerlab/ark.git
cd ark
```

Edit `params.json` (set relevant options/filters), then:
```bash
# Docker
nextflow run main.nf -params-file params.json -profile docker

# Apptainer / Singularity
nextflow run main.nf -params-file params.json -profile apptainer
```

Smoke test:
```bash
nextflow run main.nf -profile test,apptainer
```

> [!NOTE]
> You can also specify these options directly in `params.json`.

To restart from LIMA outputs, point the input at `02_LIMA` and provide the same
primer FASTA used for demultiplexing:

```bash
nextflow run src/main.nf -params-file src/params.json -profile docker \
  --entrypoint refine --global_input_dir /path/to/02_LIMA \
  --global_primers /path/to/primers.fasta
```

The `refine` checkpoint reads every `*.bam` and its matching `*.bam.pbi`, creating
missing indexes, and skips CCS, Skera, and LIMA. Each BAM is refined independently,
preserving its primer-pair group. BAM names must end in
`.<5p_primer>--<3p_primer>.bam`.

To restart from refined reads, point the input at `03_ISOSEQ_REFINE`; no primer
FASTA is needed:

```bash
nextflow run src/main.nf -params-file src/params.json -profile docker \
  --entrypoint cluster --global_input_dir /path/to/03_ISOSEQ_REFINE
```

The `cluster` checkpoint reads every `*.bam`, creating missing indexes, and skips
CCS, Skera, LIMA, and refine. BAM names must end in `_flnc.bam`; the rest of the
name is kept as the sample ID, so a full run's
`movie.IsoSeqX_bc01_5p--IsoSeqX_3p_flnc.bam` stays
`movie.IsoSeqX_bc01_5p--IsoSeqX_3p`.

`cluster_mode` controls clustering of those refined reads: `per_sample`
keeps each input primer-pair group separate; `multi_sample` clusters all groups
together under `global_prefix`; `both` produces both sets using the same refined
reads. Sample IDs retain the full primer pair, for example
`movie.IsoSeqX_bc01_5p--IsoSeqX_3p` or `movie.NEB_5p--NEB_Clontech_3p`.
Their tissue names depend on your experimental barcode assignments; the pipeline
does not infer them. Each input BAM/primer pair is treated as a separate sample,
including when files come from different sequencing runs.

`cluster_engine` picks the clustering tool: `isoseq` (default, `isoseq cluster2`;
needs tagged PacBio BAMs), `cdhit` (`cd-hit-est`, identity `cdhit_identity`,
default 0.99) or `rattle` (RATTLE cluster/correct/polish at isoform level). Both
alternatives work from any entrypoint. With them, `cluster` also reads
`*.fasta[.gz]` / `*.fastq[.gz]`, for example SRA reads that lost their PacBio
tags. Each file is one sample, named after the file without `.fast[aq](.gz)`.
BAMs are converted to FASTA once; FASTA/FASTQ inputs are used as they are.
Either way, hq holds the cluster representatives (cd-hit) or the polished
consensi (RATTLE), and singletons holds the one-read clusters. Outputs go to
`04_CDHIT_EST` / `04_RATTLE`; cd-hit also writes the read-to-cluster `.clstr`.
`cluster_mode` applies unchanged; the pooled run takes every sample's file at
once. An old `isoseq_cluster2_mode` setting is rejected: it is now
`cluster_mode`.

For `flnc`, the same three modes pool input FASTA/FASTQ files instead of
clustering: `per_sample` aligns each file separately; `multi_sample`
concatenates all files per hq/singleton class into one pooled sample named
`global_prefix`, so adapter removal, segmentation, and twin-collapsing see
all samples together per chromosome; `both` produces both sets. Do not name
an input sample `global_prefix`: its reads would merge with the pool.

A helper sh script is provided to run the pipeline on a SLURM cluster. See details below.

<details>
<summary>Click to expand</summary>


Edit the path variables at the top of `assets/hpc/ark.sh` (cache dir, container image, manifest path), then submit:

```bash
sbatch --array=1-<N> ark.sh
```

Each array task spawns one Nextflow head job that submits all compute as child SLURM jobs.

PREDICT_ORFS run as SLURM job arrays. Partition routing, array sizes, and resource tiers are documented inline in `nextflow.config` — edit there to match your cluster.

</details>

---

## Output

```
results/
├── 00_?/       *bed
└── pipeline_info/    timeline, trace, DAG
```

---

## Where to edit

| File | What |
|------|------|
| `params.json` | Genome paths, alignment settings, checkpoints — per run |
| `nextflow.config` | Compute resources, profiles, container, SLURM — rarely |
