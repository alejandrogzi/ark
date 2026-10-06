process RATTLE {
    tag "$meta.id"
    label 'process_extreme'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rattle:1.0--h5ca1c30_0' :
        'biocontainers/rattle:1.0--h5ca1c30_0' }"

    input:
    tuple val(meta), path(reads) // 1..N FASTA/FASTQ, plain or .gz; N > 1 is a pooled sample

    output:
    tuple val(meta1), path("*.singletons.fasta.gz") , optional: true, emit: singletons
    tuple val(meta2), path("*.hq.fasta.gz")         , optional: true, emit: hq
    tuple val(meta) , path("*.transcriptome.fq.gz") , optional: true, emit: transcriptome
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args  ?: '' // rattle cluster
    def args2   = task.ext.args2 ?: '' // rattle correct
    def args3   = task.ext.args3 ?: '' // rattle polish
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0' // WARN: rattle has no --version; keep in sync with environment.yml and the container

    meta1 = meta + [ singleton: true ]
    meta2 = meta + [ singleton: false ]
    """
    # INFO: RATTLE reads every file itself (-i a,b,c; FASTA/FASTQ, .gz), so pooled samples are not merged.
    # INFO: It segfaults when any of them holds no reads, so only files with at least one byte are passed.
    inputs=""
    for f in ${reads}; do
        if [ -n "\$({ case "\$f" in *.gz) zcat "\$f" ;; *) cat "\$f" ;; esac || true; } | head -c 1)" ]; then
            inputs="\${inputs:+\$inputs,}\$f"
        fi
    done

    # INFO: RATTLE always writes 4-line records, but each header keeps its source's > or @ (mixed in pooled
    # INFO: runs) and gains ",gene_cluster_N,transcript_cluster_M": emit FASTA with the bare read name
    to_fasta() {
        awk 'NR % 4 == 1 { h = substr(\$1, 2); sub(/,gene_cluster_[0-9]+,transcript_cluster_[0-9]+\$/, "", h); print ">" h }
             NR % 4 == 2' "\$1" | gzip > "\$2"
    }

    if [ -n "\$inputs" ]; then
        rattle cluster -i "\$inputs" -o . -t ${task.cpus} ${args}
        rattle correct -i "\$inputs" -c clusters.out -o . -t ${task.cpus} ${args2}
        # INFO: corrected reads are unused; RATTLE also leaves a decompressed copy of every .gz input here
        rm -f corrected.fq \$(for f in ${reads}; do case "\$f" in *.gz) echo "\${f%.gz}" ;; esac; done)

        # INFO: polished consensi -> hq; uncorrected (clusters at or below correct -r reads) -> singletons
        if [ -s consensi.fq ]; then
            rattle polish -i consensi.fq -o . -t ${task.cpus} ${args3}
            to_fasta transcriptome.fq ${prefix}.hq.fasta.gz
            gzip -c transcriptome.fq > ${prefix}.transcriptome.fq.gz
        fi
        if [ -s uncorrected.fq ]; then
            to_fasta uncorrected.fq ${prefix}.singletons.fasta.gz
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rattle: ${VERSION}
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0'

    meta1 = meta + [ singleton: true ]
    meta2 = meta + [ singleton: false ]
    """
    touch ${prefix}.hq.fasta.gz
    touch ${prefix}.singletons.fasta.gz
    touch ${prefix}.transcriptome.fq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rattle: ${VERSION}
    END_VERSIONS
    """
}
