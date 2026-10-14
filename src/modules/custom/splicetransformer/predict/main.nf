// INFO: ported from hillerlab/core modules/nextflow/splicetransformer/predict (container pinned, weights input, stub fixed)
// INFO: writes <chunk>.{donor,acceptor}_{plus,minus}.wig per chunk; probabilities < --min-prob (default 0.001) are written as 0
// INFO: the image (>= 0.0.5) bundles readable weights (/models/SpTransformer_pytorch.ckpt) and a CUDA build of torch;
// INFO: `weights` (params.splicetransformer_model_path) overrides the bundled checkpoint
process SPLICETRANSFORMER_PREDICT {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'ghcr.io/hillerlab/splicetransformer:0.0.5' }"

    input:
    tuple val(meta), path(fasta)
    path(weights) // INFO: SpTransformer_pytorch.ckpt, or [] for the bundled one

    output:
    tuple val(meta), path("splicetransformer/*.donor_plus.wig"),     emit: donor_plus
    tuple val(meta), path("splicetransformer/*.donor_minus.wig"),    emit: donor_minus
    tuple val(meta), path("splicetransformer/*.acceptor_plus.wig"),  emit: acceptor_plus
    tuple val(meta), path("splicetransformer/*.acceptor_minus.wig"), emit: acceptor_minus
    tuple val(meta), path("splicetransformer/*.wig"),                emit: all
    path "versions.yml",                                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def model = weights ? "--model ${weights}" : ''
    """
    splicetransformer predict \\
        $args \\
        $model \\
        --outdir splicetransformer \\
        --sequence ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        splicetransformer predict: \$( splicetransformer predict --version | sed 's/splicetransformer-predict //g' )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p splicetransformer
    for track in donor_plus donor_minus acceptor_plus acceptor_minus; do
        touch splicetransformer/${meta.id}.\${track}.wig
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        splicetransformer predict: stub
    END_VERSIONS
    """
}
