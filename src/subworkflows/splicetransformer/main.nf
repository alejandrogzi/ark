/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SPLICETRANSFORMER_CHUNK } from '../../modules/custom/splicetransformer/chunk/main.nf'
include { SPLICETRANSFORMER_PREDICT } from '../../modules/custom/splicetransformer/predict/main.nf'
include { SPLICETRANSFORMER_PUBLISH } from '../../modules/custom/splicetransformer/publish/main.nf'

include { UCSC_WIGTOBIGWIG as WIGTOBIGWIG_DONOR_PLUS } from '../../modules/nf-core/ucsc/wigtobigwig/main.nf'
include { UCSC_WIGTOBIGWIG as WIGTOBIGWIG_DONOR_MINUS } from '../../modules/nf-core/ucsc/wigtobigwig/main.nf'
include { UCSC_WIGTOBIGWIG as WIGTOBIGWIG_ACCEPTOR_PLUS } from '../../modules/nf-core/ucsc/wigtobigwig/main.nf'
include { UCSC_WIGTOBIGWIG as WIGTOBIGWIG_ACCEPTOR_MINUS } from '../../modules/nf-core/ucsc/wigtobigwig/main.nf'

include { BIGWIGMERGE as BIGWIGMERGE_DONOR_PLUS } from '../../modules/custom/bigtools/bigwigmerge/main.nf'
include { BIGWIGMERGE as BIGWIGMERGE_DONOR_MINUS } from '../../modules/custom/bigtools/bigwigmerge/main.nf'
include { BIGWIGMERGE as BIGWIGMERGE_ACCEPTOR_PLUS } from '../../modules/custom/bigtools/bigwigmerge/main.nf'
include { BIGWIGMERGE as BIGWIGMERGE_ACCEPTOR_MINUS } from '../../modules/custom/bigtools/bigwigmerge/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// INFO: chunks are named tmp.<contig>.chunk.<n>.fa(.gz) (same chunker layout as SpliceAI); chunk numbers
// INFO: restart per contig, so the id keeps the contig (two contigs' chunk 0 would otherwise collide)
def chunkId(chunk) {
    def m = (chunk.name =~ /^(?:tmp\.)?(.+)\.chunk\.(\d+)\./)
    return m.find() ? "${m.group(1)}.${m.group(2)}" : chunk.baseName
}

// INFO: ported from hillerlab/core subworkflows/splicetransformer; the hillerlab `splicetransformer derive`
// INFO: step is not used: SPLICING calibrates these bigWigs with the same splicing crate as SpliceAI's
workflow SPLICETRANSFORMER {
    take:
      genome         // channel: [ val(meta), [ genome ] ]
      chromsizes     // channel: [  chromsizes ]
      compression    // bool
      weights        // path [ SpTransformer_pytorch.ckpt; null = bundled in the image ]
      ch_versions    // channel: [ path(version) ]

    main:
      // INFO: runs SpliceTransformer genome-wide; only reached from SPLICING (splice_model = splicetransformer)
      // INFO: when neither precomputed scores nor a bigwig dir were given

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          CHUNKING [ GENOME -> FASTA CHUNKS ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ch_chunks = Channel.empty()
      SPLICETRANSFORMER_CHUNK(genome)

      // INFO: one item per chunk -> [ [ id: '<genome>.<contig>.<chunk>' ], chunk.fa(.gz) ]
      (compression ? SPLICETRANSFORMER_CHUNK.out.fasta_gz : SPLICETRANSFORMER_CHUNK.out.fasta)
          .flatMap {
              meta, fa ->
              def fas = fa instanceof List ? fa : [fa]
              fas.collect { it ->
                  [ [ id: meta.id + '.' + chunkId(it) ], it ]
              }
          }
          .set { ch_chunks }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          PREDICTION [ DONOR / ACCEPTOR x PLUS / MINUS ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: weights as a value channel so every chunk task gets them ([] -> bundled model)
      ch_weights = weights ? Channel.value(file(weights, checkIfExists: true)) : Channel.value([])
      SPLICETRANSFORMER_PREDICT(ch_chunks, ch_weights) // INFO: four wiggles per chunk

      // INFO: wig -> bigWig per chunk
      WIGTOBIGWIG_DONOR_PLUS(SPLICETRANSFORMER_PREDICT.out.donor_plus, chromsizes)
      WIGTOBIGWIG_DONOR_MINUS(SPLICETRANSFORMER_PREDICT.out.donor_minus, chromsizes)
      WIGTOBIGWIG_ACCEPTOR_PLUS(SPLICETRANSFORMER_PREDICT.out.acceptor_plus, chromsizes)
      WIGTOBIGWIG_ACCEPTOR_MINUS(SPLICETRANSFORMER_PREDICT.out.acceptor_minus, chromsizes)

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          MERGE + PUBLISH
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: collect every chunk of a track and merge into one genome-wide bigWig per track
      BIGWIGMERGE_DONOR_PLUS(
          WIGTOBIGWIG_DONOR_PLUS.out.bigwig
            .map { meta, bigwig -> bigwig }
            .collect()
            .map { bws -> [ [ id : 'donor_plus' ], bws ] },
      )
      BIGWIGMERGE_DONOR_MINUS(
          WIGTOBIGWIG_DONOR_MINUS.out.bigwig
            .map { meta, bigwig -> bigwig }
            .collect()
            .map { bws -> [ [ id : 'donor_minus' ], bws ] },
      )
      BIGWIGMERGE_ACCEPTOR_PLUS(
          WIGTOBIGWIG_ACCEPTOR_PLUS.out.bigwig
            .map { meta, bigwig -> bigwig }
            .collect()
            .map { bws -> [ [ id : 'acceptor_plus' ], bws ] },
      )
      BIGWIGMERGE_ACCEPTOR_MINUS(
          WIGTOBIGWIG_ACCEPTOR_MINUS.out.bigwig
            .map { meta, bigwig -> bigwig }
            .collect()
            .map { bws -> [ [ id : 'acceptor_minus' ], bws ] },
      )

      // INFO: copies the four bigWigs into one splicetransformer/ dir (published to 00_SPLICETRANSFORMER)
      // WARN: emits a bare path("splicetransformer"), no meta
      SPLICETRANSFORMER_PUBLISH(
          BIGWIGMERGE_DONOR_PLUS.out.bigwig,
          BIGWIGMERGE_DONOR_MINUS.out.bigwig,
          BIGWIGMERGE_ACCEPTOR_PLUS.out.bigwig,
          BIGWIGMERGE_ACCEPTOR_MINUS.out.bigwig,
      )

      ch_versions = ch_versions.mix(SPLICETRANSFORMER_CHUNK.out.versions)
      ch_versions = ch_versions.mix(SPLICETRANSFORMER_PREDICT.out.versions.first())

    emit:
      splicetransformer = SPLICETRANSFORMER_PUBLISH.out.splicetransformer
      versions = ch_versions
}
