/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { XLOCI_INTRON as XLOCI_EXTRACT_INTRONS } from '../../modules/custom/xloci/intron/main.nf'
include { INTRONIC as IIC_PREDICT_SPLICEOSOME } from '../../modules/custom/intronic/main.nf'
include { ISOTOOLS_CLASSIFY_INTRON } from '../../modules/custom/isotools/classify/intron/main.nf'

include { APARENT_CHUNK as XISO_APARENT_CHUNK } from '../../modules/custom/aparent/chunk/main.nf'
include { APARENT_PREDICT } from '../../modules/custom/aparent/predict/main.nf'

include { BEDGRAPHTOBIGWIG as BIGTOOLS_BEDGRAPHTOBIGWIG_FORWARD } from '../../modules/custom/bigtools/bedgraphtobigwig/main.nf'
include { BEDGRAPHTOBIGWIG as BIGTOOLS_BEDGRAPHTOBIGWIG_REVERSE } from '../../modules/custom/bigtools/bedgraphtobigwig/main.nf'

include { BIGWIGMERGE as BIGTOOLS_BIGWIGMERGE_FORWARD } from '../../modules/custom/bigtools/bigwigmerge/main.nf'
include { BIGWIGMERGE as BIGTOOLS_BIGWIGMERGE_REVERSE } from '../../modules/custom/bigtools/bigwigmerge/main.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPOLISH {
    take:
      reads                  // channel: [ val(meta), [ reads ] ] ORF-annotated (APARENT reads their 3' UTRs)
      evidence               // channel: [ val(meta), [ reads ] ] reads behind the intron frequencies (same meta.id)
      genome                 // Channel.value(path)
      chrom_sizes            // Channel.value(path)
      repeats                // path
      annotation             // channel: [ val(meta), [ annotation ] ]
      bigwigs                // channel: [ val(meta), [ bigwigs ] ]
      aparent_weights        // channel: [ val(meta), [ aparent_weights ] ]
      ch_versions            // [ meta, versions.yml ]

    main:
      // INFO: reads arrive one item per sample@chromosome (NMD-filtered reads BED)
      // INFO: meta -> [ id: '<sample>@<chr>', name: '<sample>', chr, ... ]
      ch_genome = genome.map { genome -> [ [id:genome.baseName], genome ] }
      ch_reference_transcripts = annotation

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          INTRON CLASSIFICATION [ XLOCI, INTRONIC, ISO-CLASSIFY ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      XLOCI_EXTRACT_INTRONS(ch_genome, evidence) // INFO: extract the introns of every read
      IIC_PREDICT_SPLICEOSOME(XLOCI_EXTRACT_INTRONS.out.tsv) // INFO: intronIC U2/U12 scores per intron

      // INFO: repeats are optional; [[:], []] keeps the process schedulable without them
      if (repeats) {
          Channel.value([
              [ id: "repeats" ],
              file(repeats, checkIfExists: true)
          ]).set { ch_repeats }
      } else {
          ch_repeats = Channel.value([[:], []])
      }

      // INFO: two queue channels would pair by arrival order; join by meta.id so each
      // INFO: sample@chromosome BED gets its own intronIC output -> [ meta, reads, iic ]
      evidence
        .map { meta, read -> tuple(meta.id, meta, read) }
        .join(
          IIC_PREDICT_SPLICEOSOME.out.iic
            .map { meta, iic -> tuple(meta.id, meta, iic) }
        )
        .map { id, read_meta, read, iic_meta, iic ->
          tuple(read_meta, read, iic)
        }
        .set { ch_classify_inputs }

      ISOTOOLS_CLASSIFY_INTRON(
        ch_classify_inputs,
        ch_genome,
        ch_reference_transcripts,
        ch_repeats,
        bigwigs,
      ) // INFO: emits per-intron .tsv and, with --intron-track, a BED4 of artifact/unclear/RT introns

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          POLYA SIGNAL PREDICTION [ APARENT ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

     // INFO: split each sample@chromosome BED into chunk TSVs, one APARENT task per chunk
     XISO_APARENT_CHUNK(
        reads,
        ch_genome
     )

     XISO_APARENT_CHUNK.out.chunks
        .flatMap { 
            meta, chunk_tsv ->
            def chunk_tsvs = chunk_tsv instanceof List ? chunk_tsv : [chunk_tsv]
            chunk_tsvs.withIndex().collect { it, idx ->
                [ meta + [ chunk: idx ], it ]
            }
        }
        .set { ch_aparent_chunks } // INFO: [ meta + [ chunk: idx ], chunk.tsv ]

      APARENT_PREDICT(
        ch_aparent_chunks,
        aparent_weights
      ) // INFO: forward/reverse polyA-signal bedGraphs per chunk

      // INFO: bedGraph -> bigWig per chunk, then every chunk of every sample is collected
      // INFO: and merged into ONE forward and ONE reverse bigWig (POLISH PAS caller peaks)
      BIGTOOLS_BEDGRAPHTOBIGWIG_FORWARD(
        APARENT_PREDICT.out.bg_forward, chrom_sizes
      )
      BIGTOOLS_BEDGRAPHTOBIGWIG_FORWARD.out.bigwig
        .map { meta, bw -> bw }
        .collect()
        .map { bws -> [ [id:'aparent.forward', strand:'forward'], bws ] }
        .set { ch_joined_aparent_bws_forward }
      BIGTOOLS_BIGWIGMERGE_FORWARD(ch_joined_aparent_bws_forward)

      BIGTOOLS_BEDGRAPHTOBIGWIG_REVERSE(
        APARENT_PREDICT.out.bg_reverse, chrom_sizes
      )
      BIGTOOLS_BEDGRAPHTOBIGWIG_REVERSE.out.bigwig
        .map { meta, bw -> bw }
        .collect()
        .map { bws -> [ [id:'aparent.reverse', strand:'reverse'], bws ] }
        .set { ch_joined_aparent_bws_reverse }
      BIGTOOLS_BIGWIGMERGE_REVERSE(ch_joined_aparent_bws_reverse)

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          VERSIONING
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ch_versions = ch_versions.mix(XLOCI_EXTRACT_INTRONS.out.versions)
      ch_versions = ch_versions.mix(IIC_PREDICT_SPLICEOSOME.out.versions)
      ch_versions = ch_versions.mix(ISOTOOLS_CLASSIFY_INTRON.out.versions)
      ch_versions = ch_versions.mix(XISO_APARENT_CHUNK.out.versions)
      ch_versions = ch_versions.mix(APARENT_PREDICT.out.versions)
      ch_versions = ch_versions.mix(BIGTOOLS_BEDGRAPHTOBIGWIG_FORWARD.out.versions)
      ch_versions = ch_versions.mix(BIGTOOLS_BEDGRAPHTOBIGWIG_REVERSE.out.versions)
      ch_versions = ch_versions.mix(BIGTOOLS_BIGWIGMERGE_FORWARD.out.versions)
      ch_versions = ch_versions.mix(BIGTOOLS_BIGWIGMERGE_REVERSE.out.versions)

    emit:
      introns               = ISOTOOLS_CLASSIFY_INTRON.out.tsv      // [ meta, tsv ] per sample@chromosome
      intron_track          = ISOTOOLS_CLASSIFY_INTRON.out.track    // [ meta, bed ] per sample@chromosome, unsorted
      aparent_plus          = BIGTOOLS_BIGWIGMERGE_FORWARD.out.bigwig
      aparent_minus         = BIGTOOLS_BIGWIGMERGE_REVERSE.out.bigwig
      versions              = ch_versions
}
