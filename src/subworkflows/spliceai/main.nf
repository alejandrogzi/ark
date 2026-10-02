/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SPLICEAI_CHUNK } from '../../modules/custom/spliceai/chunk/main.nf'
include { SPLICEAI_PREDICT } from '../../modules/custom/spliceai/predict/main.nf'
include { SPLICEAI_PUBLISH } from '../../modules/custom/spliceai/publish/main.nf'

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


workflow SPLICEAI {
    take:
      genome         // channel: [ val(meta), [ genome ] ]
      chromsizes     // channel: [  chromsizes ] 
      compression    // bool
      ch_versions    // channel: [ path(version) ]

    main:
      // INFO: runs the SpliceAI model genome-wide; only reached from SPLICING when neither
      // INFO: precomputed scores nor a bigwig dir were given

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          CHUNKING [ GENOME -> FASTA CHUNKS ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ch_chunks = Channel.empty()
      SPLICEAI_CHUNK(genome)

      // INFO: one item per chunk -> [ [ id: '<genome>.<chunk>' ], chunk.fa(.gz) ]
      if (compression) {
        SPLICEAI_CHUNK.out.fasta_gz
            .flatMap { 
                meta, fa ->
                def fas = fa instanceof List ? fa : [fa]
                fas.collect { it ->
                    // INFO: format of chunks is tmp.chr1.chunk.1.fasta.gz
                    // INFO: baseName strips .gz, so the chunk number is parts[-2]
                    def parts = it.baseName.split('\\.')
                    def chunk = parts[-2]
                    [ [ id: meta.id + '.' + chunk ], it ]
                }
            }
            .set { ch_chunks }
      } else {
        SPLICEAI_CHUNK.out.fasta
            .flatMap { 
                meta, fa ->
                def fas = fa instanceof List ? fa : [fa]
                fas.collect { it ->
                    def parts = it.baseName.split('\\.')
                    def chunk = parts[-1]
                    [ [ id: meta.id + '.' + chunk ], it ]
                }
            }
            .set { ch_chunks }
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          PREDICTION [ DONOR / ACCEPTOR x PLUS / MINUS ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      SPLICEAI_PREDICT(ch_chunks) // INFO: four wiggles per chunk

      // INFO: wig -> bigWig per chunk
      WIGTOBIGWIG_DONOR_PLUS(SPLICEAI_PREDICT.out.donor_plus, chromsizes)
      WIGTOBIGWIG_DONOR_MINUS(SPLICEAI_PREDICT.out.donor_minus, chromsizes)
      WIGTOBIGWIG_ACCEPTOR_PLUS(SPLICEAI_PREDICT.out.acceptor_plus, chromsizes)
      WIGTOBIGWIG_ACCEPTOR_MINUS(SPLICEAI_PREDICT.out.acceptor_minus, chromsizes)

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

      // INFO: copies the four bigWigs into one spliceai/ dir (published to 00_SPLICEAI)
      // WARN: emits a bare path("spliceai"), no meta
      SPLICEAI_PUBLISH(
          BIGWIGMERGE_DONOR_PLUS.out.bigwig,
          BIGWIGMERGE_DONOR_MINUS.out.bigwig,
          BIGWIGMERGE_ACCEPTOR_PLUS.out.bigwig,
          BIGWIGMERGE_ACCEPTOR_MINUS.out.bigwig,
      )
      
    emit:
      spliceai = SPLICEAI_PUBLISH.out.spliceai
      versions = ch_versions
}
