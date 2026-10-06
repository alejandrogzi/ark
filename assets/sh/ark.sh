#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# SLURM submission script for isopipe
#
# Each array task runs one soft-masking job. Nextflow itself is the
# "main job" — it submits all compute work as child SLURM jobs and only needs
# a small memory footprint.
#
# MANIFEST FILE FORMAT (species_list)
# ─────────────────────────────────────
# A tab-separated file with one run per line, no header:
#
#   <prefix>  <species_name>  <input_dir>  <genome> <annotation> <primers> <repeats> <spliceai> <database> <selenocysteine> <entrypoint>
#
# Example (hg38):
#   hg38  Homo sapiens  /path/to/input_dir  hg38.2bit  hg38.bed  primers.fa  repeats.bed  splceai_dir prot.dmdb  selenocysteines.bed isoseq
#
# Paths must be absolute. 
#
# USAGE
# ─────
# Edit the four path variables below, then submit with:
#   sbatch --array=1-<N> do_isopipe.sh
# where <N> is the number of lines in your manifest file.
# ─────────────────────────────────────────────────────────────────────────────

#SBATCH --job-name=ISOPIPE
#SBATCH --array=1-10        # set upper bound to number of lines in species_list
#SBATCH -t 2-0
#SBATCH --output=/path/to/logs/%A.%a.out  # MODIFY THIS!
#SBATCH --error=/path/to/logs/%A.%a.err   # MODIFY THIS!
#SBATCH --mem=20G          # memory for the Nextflow process itself (not compute jobs)
#SBATCH -p long            # partition name
#SBATCH -q long            # queue name

# ── Load required modules (adjust to your cluster's module system) ────────────
module load nextflow
module load openjdk

# ── Environment ───────────────────────────────────────────────────────────────
export SLURM_SKIP_EPILOG=1

# Directory where Apptainer caches pulled container images
export NXF_APPTAINER_CACHEDIR=/scratch/$USER/isopipe/apptainer

# Give Nextflow's JVM enough heap for large runs (thousands of jobs)
export NXF_OPTS="-Xms4g -Xmx16g"

# ── Paths — edit these ────────────────────────────────────────────────────────
species_list="/path/to/manifest.tsv"   # tab-separated manifest (see format above)
working_dir="/path/to/output"          # one subdirectory per species will be created here
pipeline_dir="/path/to/isopipe"        # cloned pipeline repo

# ── Parse manifest line for this array task ───────────────────────────────────
row=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$species_list")
prefix=$(echo "$row" | cut -f1)
species_name=$(echo "$row" | cut -f2)
input_dir=$(echo "$row" | cut -f3)
genome=$(echo "$row"   | cut -f4)
annotation=$(echo "$row" | cut -f5)
primers=$(echo "$row" | cut -f6)
repeats=$(echo "$row" | cut -f7)
spliceai=$(echo "$row" | cut -f8)
database=$(echo "$row" | cut -f9)
selenocysteine=$(echo "$row" | cut -f10)
entrypoint=$(echo "$row" | cut -f11)


if [[ -z "$prefix" || -z "$species_name" || -z "$input_dir" || -z "$genome" || -z "$annotation" || -z "$primers" || -z "$repeats" || -z "$spliceai" || -z "$database" || -z "$selenocysteine" ]]; then
    echo "ERROR: could not parse line ${SLURM_ARRAY_TASK_ID} of ${species_list}" >&2
    exit 1
fi

# ── Per-pair working directory ─────────────────────────────────────────────────
run_dir="${working_dir}/${prefix}_isopipe"
mkdir -p "${run_dir}/logs"

# ── Write params.json for this pair ───────────────────────────────────────────
# Scientific parameters go here; infrastructure stays in nextflow.config.
cat > "${run_dir}/params.json" <<EOF
{
    "// 1)": "Entrypoint / Global input options [ entrypoint: isoseq, map ] ───────────────────────────────────",
    "entrypoint": "${entrypoint}",
    "global_input_dir": "${input_dir}",
    "global_output_dir": "${run_dir}",
    "global_genome": "${genome}",
    "global_annotation": "${annotation}",
    "global_primers": ${primers}",
    "global_repeats": "${repeats}",
    "global_prefix": "${prefix}",
    "global_species_name": "${species_name}",

    "// 2)": "CCS / Isoseq cluster options [ clustering: per_sample, multi_sample, both ] ─────",
    "ccs_chunk": 1000,
    "cluster_mode": "multi_sample",

    "// 2b)": "Chunking options ────────────────────────────────────────",
    "fxsplit_chunk_size": 50000,

    "// 3)": "Aligner options [minimap2, ultra] ────────────────────────────────────────",
    "aligner": "minimap2",

    "// 4)": "Minimap2 options [ scoring options: spliceai, minisplice ] ────────────────────────────────────────",
    "minimap2_index_path": null,
    "minimap2_align_use_junc_bed": true,
    "minimap2_align_use_splice_scores": true,
    "minimap2_align_splicing_algorithm": "spliceai",

    "// 5)": "Ultra options ────────────────────────────────────────",
    "ultra_index": null,
    "ultra_use_annotation": true,

    "// 6)": "SpliceAI/minisplice options ────────────────────────────────────────",
    "spliceai_bigwigs_dir": ${spliceai},
    "spliceai_scores_path": null,
    "minisplice_scores_path": null,

    "// 7)": "Isotools adapter / collapse / cigar options ────────────────────────────────────────",
    "isotools_adapter_remove_adapters": false,
    "isotools_cigar_extension_extend": true,
    "collapse_shrink_twins": false,

    "// 8)": "Isotools segment options ────────────────────────────────────────",
    "isotools_segment_identity_threshold": 98,
    "isotools_segment_max_five_prime_clip": 20,
    "isotools_segment_max_three_prime_clip": 20,

    "// 9)": "xORF options ────────────────────────────────────────",
    "xorf_protein_database": ${database},
    "xorf_chunk_size": 10,
    "xorf_predict_keep_raw": true,
    "xorf_predict_min_score_max_predictions": 0.50,
    "xorf_predict_max_predictions": 3,
    "xorf_predict_threshold": 0.01,
    "xorf_selenocysteine_codons": ${selenocysteine},

    "// 10)": "Track upload options ────────────────────────────────────────",
    "load_track": false
}
EOF

# cd into run_dir so each run's .nextflow.log is saved there
cd "$run_dir"

nextflow run "${pipeline_dir}/main.nf" \
    -params-file "${run_dir}/params.json" \
    -profile     apptainer,slurm \
    -w           "${run_dir}/work"
