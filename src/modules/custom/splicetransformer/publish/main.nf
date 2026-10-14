// INFO: ported from hillerlab/core modules/nextflow/splicetransformer/publish
// INFO: copies the four merged bigWigs ({donor,acceptor}_{plus,minus}.bw) into one splicetransformer/ dir;
// INFO: the splicing crate classifies them by the donor/acceptor and plus/minus tokens in the name
process SPLICETRANSFORMER_PUBLISH {
    tag "publish"
    label 'process_single'

    input:
    tuple val(meta), path(donor_plus)
    tuple val(meta1), path(donor_minus)
    tuple val(meta2), path(acceptor_plus)
    tuple val(meta3), path(acceptor_minus)

    output:
    path("splicetransformer"), emit: splicetransformer

    script:
    """
    mkdir -p splicetransformer
    cp ${donor_plus} splicetransformer/
    cp ${donor_minus} splicetransformer/
    cp ${acceptor_plus} splicetransformer/
    cp ${acceptor_minus} splicetransformer/
    """

    stub:
    """
    mkdir -p splicetransformer
    touch splicetransformer/${donor_plus.name} splicetransformer/${donor_minus.name} \\
        splicetransformer/${acceptor_plus.name} splicetransformer/${acceptor_minus.name}
    """
}
