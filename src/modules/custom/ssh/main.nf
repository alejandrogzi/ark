process RSYNC_SSH {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '':
        'ghcr.io/alejandrogzi/isox-rs:latest' }"

    input:
    tuple val(meta), path(bigbed)
    val user
    val server
    val target_dir
    val web
    val species

    output:
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bb = bigbed.baseName + '.bb'
    """
    ssh \\
    ${user}@${server} \\
    mkdir -p ${target_dir}/${species}/isopipe

    rsync \\
    -av \\
    ${bigbed} \\
    ${user}@${server}:${target_dir}/${species}/isopipe/${bb}

    ssh \\
    ${user}@${server} \\
    ln -sf ${target_dir}/${species}/isopipe/${bb} ${web}/${species}/${bb}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rsync: \$(rsync --version | head -n1 | sed 's/rsync version //g; s/  .*//')
        ssh: \$(ssh --version | sed 's/OpenSSH_//g')
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rsync: \$(rsync --version | head -n1 | sed 's/rsync version //g; s/  .*//')
        ssh: \$(ssh --version | sed 's/OpenSSH_//g')
    END_VERSIONS
    """
}
