/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CHAIN_COLLAPSE } from '../../modules/custom/collapse/chain/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// INFO: engine extension point. An engine takes the segmented reads of one [ sample_id, chr ] per
// INFO: fusion-detector class and returns transcript models whose BED12 names are real read names
// INFO: (iso-segment tags intact) plus #CN<support>, so ORF calling and polishing stay engine-agnostic.
// INFO: isoquant/isocall (reserved, PLAN.md §6.1) additionally need HQ BAMs from segmentation and a
// INFO: GTF -> BED12 adapter that picks representatives from the engine's read -> transcript table.
workflow RECONSTRUCT {
    take:
      ch_reads   // [ meta + [ class: free|fusion ], bed ] per [ sample_id, chr ]
      ch_ref     // [ meta, reference transcripts BED12 ]
      engine     // string [ chain ]

    main:
      ch_versions = Channel.empty()

      switch (engine) {
        case 'chain':
          CHAIN_COLLAPSE(ch_reads, ch_ref)
          ch_models  = CHAIN_COLLAPSE.out.models
          ch_support = CHAIN_COLLAPSE.out.support
          ch_versions = ch_versions.mix(CHAIN_COLLAPSE.out.versions)
          break

        default:
          error "ERROR: reconstruct_engine '${engine}' is not implemented (v2.1.0 ships: chain)"
      }

    emit:
      models   = ch_models  // [ meta, models.bed ]
      support  = ch_support // [ meta, support.tsv, counts.tsv ]
      versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
