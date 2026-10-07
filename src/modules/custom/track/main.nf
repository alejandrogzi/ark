process TRACKDB {
    tag "${prefix} trackDb"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '':
        'ghcr.io/alejandrogzi/isox-rs:v2.1.0' }"

    input:
    val browser
    val species
    val track
    val additional_columns
    val prefix

    output:
    path "*.ra",          emit: schema
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // INFO: the trackDb stanza is written inline, like the autosql schemas: no template file to locate.
    // INFO: bigDataUrl names must match the <id>.bb files RSYNC_SSH uploads.
    """
    cat <<-'EOF' > ${prefix}.schema.ra
    ##########################################################
    # Isopipe track description
    ##########################################################
    track ${track}
    compositeTrack on
    shortLabel HL ISOPIPE annotations
    longLabel HL Isopipe annotation track
    group genes
    priority 2
    visibility pack
    itemRgb on
    type bigBed 37
    searchPriority 2.07207

            track ${prefix}.pass
            parent ${track}
            subtrack ${track}
            shortLabel HL pass
            longLabel HL_Pass
            group genes
            priority 1
            visibility pack
            itemRgb On
            type bigBed ${additional_columns}
            bigDataUrl ${browser}/${species}/${prefix}.pass.bb

            track ${prefix}.retention
            parent ${track}
            subtrack ${track}
            shortLabel HL retention
            longLabel HL_Retention
            group genes
            priority 1
            visibility pack
            itemRgb On
            type bigBed ${additional_columns}
            bigDataUrl ${browser}/${species}/${prefix}.retentions.bb

            track ${prefix}.intraprimming
            parent ${track}
            subtrack ${track}
            shortLabel HL intraprimming
            longLabel HL_Intraprimming
            group genes
            priority 1
            visibility pack
            itemRgb On
            type bigBed ${additional_columns}
            bigDataUrl ${browser}/${species}/${prefix}.intraprimming.bb

            track ${prefix}.truncation
            parent ${track}
            subtrack ${track}
            shortLabel HL truncation
            longLabel HL_Truncation
            group genes
            priority 1
            visibility pack
            itemRgb On
            type bigBed ${additional_columns}
            bigDataUrl ${browser}/${species}/${prefix}.truncations.bb

            track ${prefix}.rt
            parent ${track}
            subtrack ${track}
            shortLabel HL rt
            longLabel HL_rt
            group genes
            priority 1
            visibility pack
            itemRgb On
            type bigBed ${additional_columns}
            bigDataUrl ${browser}/${species}/${prefix}.rt.bb

            track ${prefix}.trash
            parent ${track}
            subtrack ${track}
            shortLabel HL trash
            longLabel HL_trash
            group genes
            priority 1
            visibility pack
            itemRgb On
            type bigBed ${additional_columns}
            bigDataUrl ${browser}/${species}/${prefix}.trash.bb

            track ${prefix}.orphans
            parent ${track}
            subtrack ${track}
            shortLabel HL orphans
            longLabel HL_orphans
            group genes
            priority 1
            visibility pack
            itemRgb On
            type bigBed ${additional_columns}
            bigDataUrl ${browser}/${species}/${prefix}.pass.scraps.bb

            track ${prefix}.duplicates
            parent ${track}
            subtrack ${track}
            shortLabel HL duplicates
            longLabel HL_duplicates
            group genes
            priority 1
            visibility pack
            itemRgb On
            type bigBed ${additional_columns}
            bigDataUrl ${browser}/${species}/${prefix}.pass.duplicates.bb

            track ${prefix}.fusions
            parent ${track}
            subtrack ${track}
            shortLabel HL fusions
            longLabel HL_fusions
            group genes
            priority 1
            visibility pack
            itemRgb On
            type bigBed 12
            bigDataUrl ${browser}/${species}/${prefix}.fusions.bb

            track ${prefix}.nmd
            parent ${track}
            subtrack ${track}
            shortLabel HL nmd
            longLabel HL_nmd
            group genes
            priority 1
            visibility pack
            itemRgb On
            type bigBed 12
            bigDataUrl ${browser}/${species}/${prefix}.nmd.bb

    ##########################################################
    # Isopipe bigWig track description [ spliceAi + aparent ]
    ##########################################################
    track ${track}_bigwig
    compositeTrack on
    shortLabel HL ISOPIPE bigWig track
    longLabel HL Isopipe bigWig track
    group genes
    priority 2
    visibility full
    itemRgb on
    type bigWig
    searchPriority 2.07207

            track ${prefix}.aparent.forward
            parent ${track}_bigwig
            subtrack ${track}_bigwig
            shortLabel HL aparent forward
            longLabel HL_aparent_forward
            group genes
            priority 1
            visibility full
            itemRgb On
            type bigWig
            bigDataUrl ${browser}/${species}/${prefix}.aparent.forward.bw
            color 58,130,27

            track ${prefix}.aparent.reverse
            parent ${track}_bigwig
            subtrack ${track}_bigwig
            shortLabel HL aparent reverse
            longLabel HL_aparent_reverse
            group genes
            priority 1
            visibility full
            itemRgb On
            type bigWig
            bigDataUrl ${browser}/${species}/${prefix}.aparent.reverse.bw
            color 212,154,25

            track ${prefix}.spliceai.donor.forward
            parent ${track}_bigwig
            subtrack ${track}_bigwig
            shortLabel HL spliceai donor forward
            longLabel HL_spliceai_donor_forward
            group genes
            priority 1
            visibility full
            itemRgb On
            type bigWig
            bigDataUrl ${browser}/${species}/${prefix}.spliceai.donor.forward.bw
            color 38,53,201

            track ${prefix}.spliceai.donor.reverse
            parent ${track}_bigwig
            subtrack ${track}_bigwig
            shortLabel HL spliceai donor reverse
            longLabel HL_spliceai_donor_reverse
            group genes
            priority 1
            visibility full
            itemRgb On
            type bigWig
            bigDataUrl ${browser}/${species}/${prefix}.spliceai.donor.reverse.bw
            color 114,33,148

            track ${prefix}.spliceai.acceptor.forward
            parent ${track}_bigwig
            subtrack ${track}_bigwig
            shortLabel HL spliceai acceptor forward
            longLabel HL_spliceai_acceptor_forward
            group genes
            priority 1
            visibility full
            itemRgb On
            type bigWig
            bigDataUrl ${browser}/${species}/${prefix}.spliceai.acceptor.forward.bw
            color 38,53,201

            track ${prefix}.spliceai.acceptor.reverse
            parent ${track}_bigwig
            subtrack ${track}_bigwig
            shortLabel HL aparent forward
            longLabel HL_aparent_forward
            group genes
            priority 1
            visibility full
            itemRgb On
            type bigWig
            bigDataUrl ${browser}/${species}/${prefix}.spliceai.acceptor.reverse.bw
            color 114,33,148
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version 2>&1 | head -n1 | sed 's/Bash (GNU bash) //g; s/  .*//')
    END_VERSIONS
    """

    stub:
    """
    touch ${prefix}.schema.ra

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version 2>&1 | head -n1 | sed 's/Bash (GNU bash) //g; s/  .*//')
    END_VERSIONS
    """
}
