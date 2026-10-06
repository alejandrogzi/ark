process CDHIT_EST {
    tag "$meta.id"
    label 'process_extreme'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cd-hit:4.8.1--h5ca1c30_13' :
        'biocontainers/cd-hit:4.8.1--h5ca1c30_13' }"

    input:
    tuple val(meta), path(reads) // 1..N FASTA/FASTQ, plain or .gz; N > 1 is a pooled sample

    output:
    tuple val(meta1), path("*.singletons.fasta.gz"), optional: true, emit: singletons
    tuple val(meta2), path("*.hq.fasta.gz")        , optional: true, emit: hq
    tuple val(meta) , path("*.clstr.gz")           , emit: clstr
    path  "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    meta1 = meta + [ singleton: true ]
    meta2 = meta + [ singleton: false ]
    """
    # INFO: cd-hit-est takes one -i and seeks it to write representatives, so every staged file is
    # INFO: streamed into one plain FASTA (this is the multi_sample merge). Format comes from the first
    # INFO: byte (@ = 4-line FASTQ); headers are cut at the first space, which is what -d 0 IDs keep.
    for f in ${reads}; do
        case "\$f" in *.gz) zcat "\$f" ;; *) cat "\$f" ;; esac \\
        | awk 'NR == 1 { fq = /^@/ }
               fq && NR % 4 == 1 { print ">" substr(\$1, 2); next }
               fq && NR % 4 != 2 { next }
               /^>/ { print \$1; next }
               { print }'
    done > reads.fa

    cd-hit-est \\
        -i reads.fa \\
        -o ${prefix}.cdhit.fa \\
        -d 0 \\
        -M 0 \\
        -T ${task.cpus} \\
        ${args}

    # INFO: a .clstr member line is "0	2457nt, >read... *"; one-member clusters are singletons, every
    # INFO: other representative is hq (the same disjoint split BAM_TO_FA makes on is:i:1).
    # ponytail: split by read ID, so pooled inputs need unique IDs (PacBio movie/zmw and SRA ids are)
    awk '/^>Cluster/ { if (n == 1) print id; n = 0; next }
         { n++; id = \$3; sub(/^>/, "", id); sub(/[.][.][.]\$/, "", id) }
         END { if (n == 1) print id }' ${prefix}.cdhit.fa.clstr > singletons.ids

    awk 'FILENAME == "singletons.ids" { single[\$1]; next }
         /^>/ { out = (substr(\$1, 2) in single) ? "singletons.fa" : "hq.fa" }
         { print > out }' singletons.ids ${prefix}.cdhit.fa

    for c in hq singletons; do
        if [ -s \$c.fa ]; then gzip -c \$c.fa > ${prefix}.\$c.fasta.gz; fi
    done
    gzip ${prefix}.cdhit.fa.clstr
    rm -f reads.fa hq.fa singletons.fa singletons.ids ${prefix}.cdhit.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cd-hit: \$( cd-hit-est -h | head -n 1 | sed 's/^.*====== CD-HIT version //;s/ (built on .*) ======\$//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    meta1 = meta + [ singleton: true ]
    meta2 = meta + [ singleton: false ]
    """
    touch ${prefix}.hq.fasta.gz
    touch ${prefix}.singletons.fasta.gz
    touch ${prefix}.cdhit.fa.clstr.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cd-hit: \$( cd-hit-est -h | head -n 1 | sed 's/^.*====== CD-HIT version //;s/ (built on .*) ======\$//' )
    END_VERSIONS
    """
}
