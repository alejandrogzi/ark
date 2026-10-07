#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPROCESSING } from '../subworkflows/preprocessing/main.nf'

include { SPLIT_ALIGN_CLEAN_CHUNKS } from '../subworkflows/split_align/main.nf'

include { PREPOLISH as ISOTOOLS_PREPOLISH } from '../subworkflows/prepolish/main.nf'
include { POLISH as ISOTOOLS_POLISH } from '../subworkflows/polish/main.nf'

include { LOAD_TRACK as LOAD_PASS_TRACK } from '../subworkflows/track/main.nf'
include { LOAD_TRACK as LOAD_DUPLICATES_TRACK } from '../subworkflows/track/main.nf'
include { LOAD_TRACK as LOAD_ORPHANS_TRACK } from '../subworkflows/track/main.nf'
include { LOAD_TRACK as LOAD_TRASH_TRACK } from '../subworkflows/track/main.nf'
include { LOAD_TRACK as LOAD_RETENTIONS_TRACK } from '../subworkflows/track/main.nf'
include { LOAD_TRACK as LOAD_TRUNCATIONS_TRACK } from '../subworkflows/track/main.nf'
include { LOAD_TRACK as LOAD_INTRAPRIMMING_TRACK } from '../subworkflows/track/main.nf'
include { LOAD_TRACK as LOAD_RT_TRACK } from '../subworkflows/track/main.nf'
include { LOAD_TRACK as LOAD_FUSIONS_TRACK } from '../subworkflows/track/main.nf'
include { LOAD_TRACK as LOAD_NMD_TRACK } from '../subworkflows/track/main.nf'

include { XORF as XORF_PREDICT_ORFS } from '../../modules/xorf/src/subworkflows/xorf/main.nf'
include { XORF as XORF_PREDICT_FUSION_ORFS } from '../../modules/xorf/src/subworkflows/xorf/main.nf'

include { WGET as WGET_APARENT_WEIGHTS } from '../modules/nf-core/wget/main.nf'
include { WGET as WGET_SAMBA_WEIGHTS } from '../modules/nf-core/wget/main.nf'

include { ISOTOOLS_NMD as ISOTOOLS_NMD_FILTER } from '../modules/custom/isotools/nmd/main.nf'

include { GAWK_JOIN as JOIN_FUSIONS } from '../modules/custom/gawk/join/main.nf'
include { GAWK_JOIN as JOIN_NMD } from '../modules/custom/gawk/join/main.nf'
include { GAWK_JOIN as JOIN_INTRONS } from '../modules/custom/gawk/join/main.nf'

include { BEDTOBIGBED as BEDTOBIGBED_FUSIONS } from '../modules/custom/bigtools/bedtobigbed/main.nf'
include { BEDTOBIGBED as BEDTOBIGBED_NMD } from '../modules/custom/bigtools/bedtobigbed/main.nf'
include { BEDTOBIGBED as BEDTOBIGBED_INTRONS } from '../modules/custom/bigtools/bedtobigbed/main.nf'

include { PUBLISH as PUBLISH_ADDITIONAL_BIGBEDS } from '../modules/custom/publish/main.nf'
include { TRACKDB } from '../modules/custom/track/main.nf'

include { AUTOSQL_BASE } from '../modules/custom/autosql/base/main.nf'
include { AUTOSQL_SCHEMA } from '../modules/custom/autosql/schema/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ARK {
    main:
      ch_versions = Channel.empty()
      ch_reads = Channel.empty()

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          AUTOSQL
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: .as schemas for the bigBed tracks (base: fusion/NMD tracks, schema: polish tracks)
      AUTOSQL_BASE()
      AUTOSQL_SCHEMA()

      autosql = AUTOSQL_BASE.out.autosql
      schema = AUTOSQL_SCHEMA.out.autosql

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          PRE-PROCESSING [ INDEXES, SPLICE SCORES, READS ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: genome, indexes, splice scores, annotation and reads [ meta, fasta ] (one per sample and hq/singleton class)
      PREPROCESSING(
        params.entrypoint,
        params.global_input_dir,
        params.global_primers,
        params.global_genome,
        params.global_annotation,
        params.ccs_chunk,
        params.cluster_mode,
        params.cluster_engine,
        params.flnc_input_state,
        params.xorf_protein_database,
        params.xorf_custom_database,
        params.xorf_raw_database,
        params.minimap2_index_path,
        params.minimap2_align_use_splice_scores,
        params.minimap2_align_splicing_algorithm,
        params.spliceai_bigwigs_dir,
        params.minisplice_scores_path,
        params.spliceai_scores_path,
        params.spliceai_chunk_compression,
        params.global_prefix,
        params.aligner,
        params.ultra_use_annotation,
        params.ultra_index,
        params.desalt_index,
        params.pbmm2_index,
        params.skera_is_kinnex_library,
        ch_versions
      )

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ALIGNMENT [ SPLIT_ALIGN_CLEAN_CHUNKS ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: chunk, align, segment and detect fusions; both outputs are per sample and chromosome
      // INFO: [ [ id: sample_id, single_end: true, chr ], bed ]
      if (params.aligner in ['mm2', 'ultra', 'pbmm2', 'desalt', 'ark', 'flair']) {
        SPLIT_ALIGN_CLEAN_CHUNKS(
          PREPROCESSING.out.reads,
          PREPROCESSING.out.genome,
          PREPROCESSING.out.genome_index,
          PREPROCESSING.out.reference_transcripts,
          PREPROCESSING.out.splice_scores,
          params.aligner,
          params.minimap2_align_use_junc_bed,
          params.isotools_adapter_remove_adapters,
          params.collapse_shrink_twins,
          params.isotools_cigar_extension_extend,
          params.minimap2_align_do_second_pass,
          params.reconstruct_engine,
          ch_versions
        )

        ch_reads = SPLIT_ALIGN_CLEAN_CHUNKS.out.reads
        ch_fusions = SPLIT_ALIGN_CLEAN_CHUNKS.out.fusions
      } else {
        error """
        ERROR: Unsupported aligner: '${params.aligner}'.
        Valid aligners are: mm2, pbmm2, desalt, ultra, ark, flair
        """.stripIndent()
        System.exit(1)
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          WEIGHTS [ XORF ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: RNAsamba weights for xORF, local file if given, otherwise downloaded
      ch_samba_weights = Channel.empty()
      if (params.xorf_samba_local_weights) {
        ch_samba_weights = Channel.value(
          file(params.xorf_samba_local_weights, checkIfExists: true)
        ).map { path -> [ [id : path.baseName ], path ] }
      } else {
        WGET_SAMBA_WEIGHTS(
          Channel.value(
            params.xorf_samba_weights
          ).map { url -> [ [id : url.tokenize('/')[-1]], url ] }
        )
        ch_samba_weights = WGET_SAMBA_WEIGHTS.out.outfile
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ORF CALLING [ XORF ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // channel comes per chromosome previously splitted at segmentation step
      // metadata -> [ id: sampleId, single_end: true, chr: meta.chr ]
      XORF_PREDICT_ORFS(
          ch_reads,
          PREPROCESSING.out.genome,
          PREPROCESSING.out.database,
          params.global_output_dir,
          params.xorf_chunk_size,
          ch_samba_weights,
          params.xorf_predict_keep_raw,
          params.xorf_selenocysteine_codons,
          params.xorf_skip_netstart,
          params.xorf_rename_deactivate,
          false, // xorf_do_polishing
          true,  // xorf_skip_joined_concat
          false, // xorf_run_only_on
          null,  // xorf_run_only_mode
          null,  // xorf_run_only_target
          Channel.empty() // xorf_database_versions
      )
      XORF_PREDICT_ORFS.out.files
          .map { meta, bed, tsv -> [ meta, bed ] }
          .set { ch_orf_predictions_bed }

      // INFO: same ORF calling on fusion reads; per-chromosome beds are joined into one
      // INFO: <name>.fusions.bed + bigBed per sample
      XORF_PREDICT_FUSION_ORFS(
          ch_fusions,
          PREPROCESSING.out.genome,
          PREPROCESSING.out.database,
          params.global_output_dir,
          params.xorf_chunk_size,
          ch_samba_weights,
          params.xorf_predict_keep_raw,
          params.xorf_selenocysteine_codons,
          params.xorf_skip_netstart,
          params.xorf_rename_deactivate,
          false,
          true,
          false, // xorf_run_only_on
          null,  // xorf_run_only_mode
          null,  // xorf_run_only_target
          Channel.empty() // xorf_database_versions
      )
      XORF_PREDICT_FUSION_ORFS.out.files
          .map { meta, bed, tsv -> [ meta.name, meta, bed ] }
          .groupTuple()
          .map { name, metas, files ->
              [ [ id: name + '.fusions', name: name ], files ]
          }
          .set { ch_fusion_orf_predictions_bed }
      JOIN_FUSIONS(ch_fusion_orf_predictions_bed, 'bed')
      BEDTOBIGBED_FUSIONS(JOIN_FUSIONS.out.output, PREPROCESSING.out.chrom_sizes, autosql)

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          NMD CALLING [ ISOTOOLS_NMD_FILTER ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: per sample@chromosome: reads (feed pre-polishing) and NMD candidates
      // INFO: (joined into one <name>.nmd.bed + bigBed per sample)
      ISOTOOLS_NMD_FILTER(ch_orf_predictions_bed)

      ISOTOOLS_NMD_FILTER.out.nmd
          .map { meta, bed -> [ meta.name, meta, bed ] }
          .groupTuple()
          .map { name, metas, files ->
              [ [ id: name + '.nmd', name: name ], files ]
          }
          .set { ch_nmd_bed }
      JOIN_NMD(ch_nmd_bed, 'bed')
      BEDTOBIGBED_NMD(JOIN_NMD.out.output, PREPROCESSING.out.chrom_sizes, autosql)

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          APARENT WEIGHTS [ POLYA PEAKS ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: APARENT model weights, local file if given, otherwise downloaded
      ch_aparent_weights = Channel.empty()
      if (params.aparent_predict_weights_local_path) {
          ch_aparent_weights = Channel.value(
            file(params.aparent_predict_weights_local_path, checkIfExists: true)
          ).map { path -> [ [id : path.baseName ], path ] }
      } else {
          WGET_APARENT_WEIGHTS(
              Channel.value(
                params.aparent_predict_weights
              ).map { url -> [ [id : url.tokenize('/')[-1]], url ] }
          )
          ch_aparent_weights = WGET_APARENT_WEIGHTS.out.outfile
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          PRE-POLISHING [ APARENT, INTRON, RETENTION, INTRAPRIMMING, BIGWIG  ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: intron classification counts reads (seen / spanned), so with a reconstruction engine it reads
      // INFO: the uncollapsed evidence, re-keyed to the xORF meta ([ id: <sample>@<chr>, name: <sample>, chr ])
      // INFO: so the POLISH join on meta.id still pairs it with the models
      ch_prepolish_reads = params.reconstruct_engine == 'none'
          ? ISOTOOLS_NMD_FILTER.out.reads
          : SPLIT_ALIGN_CLEAN_CHUNKS.out.evidence.map { meta, bed ->
              [ [ id: "${meta.id}@${meta.chr}", name: meta.id, chr: meta.chr ], bed ]
            }

      // INFO: per sample@chromosome: intron classification (tsv + BED4 track) and polyA peaks (bigWig per strand)
      ISOTOOLS_PREPOLISH(
          ch_prepolish_reads,
          PREPROCESSING.out.genome,
          PREPROCESSING.out.chrom_sizes,
          params.global_repeats,
          PREPROCESSING.out.reference_transcripts,
          PREPROCESSING.out.bigwigs,
          ch_aparent_weights,
          ch_versions
      )

      // INFO: polish input per sample@chromosome, joined on meta.id: [ meta, reads, intron tsv, ORF tsv ]
      // INFO: keys as plain strings: xORF builds ids as GStrings, and a GString never equals a String
      ISOTOOLS_NMD_FILTER.out.reads
        .map { meta, read -> tuple(meta.id.toString(), meta, read) }
        .join(
          ISOTOOLS_PREPOLISH.out.introns
            .map { meta, introns -> tuple(meta.id.toString(), introns) }
        )
        .map { id, meta, read, introns ->
          tuple(id, meta, read, introns)
        }
        .join(
          XORF_PREDICT_ORFS.out.files
              .map { meta, bed, tsv -> tuple(meta.id.toString(), tsv) }
        )
        .map { id, meta, read, introns, tsv ->
          tuple(meta, read, introns, tsv)
        }
        .set { ch_full_length_reads }

      // INFO: artifact / unclear / RT introns, joined into one <name>.introns.bed + bigBed per sample.
      // INFO: BED4, so no autosql
      ISOTOOLS_PREPOLISH.out.intron_track
          .map { meta, bed -> [ meta.name, meta, bed ] }
          .groupTuple()
          .map { name, metas, files ->
              [ [ id: name + '.introns', name: name ], files ]
          }
          .set { ch_intron_bed }
      JOIN_INTRONS(ch_intron_bed, 'bed')
      BEDTOBIGBED_INTRONS(JOIN_INTRONS.out.output, PREPROCESSING.out.chrom_sizes, [])

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ADDITIONAL BIGBEDS [ BEDTOBIGBED_FUSIONS, BEDTOBIGBED_NMD, BEDTOBIGBED_INTRONS ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: fusion, NMD and intron bigBeds grouped per sample and published to 12_POLISH/BB
      ch_additional_bbs = Channel.empty()
      ch_additional_bbs = ch_additional_bbs.mix(BEDTOBIGBED_FUSIONS.out.bigbed)
      ch_additional_bbs = ch_additional_bbs.mix(BEDTOBIGBED_NMD.out.bigbed)
      ch_additional_bbs = ch_additional_bbs.mix(BEDTOBIGBED_INTRONS.out.bigbed)
      ch_additional_bbs.map { meta, file -> [meta.name, meta, file] }
         .groupTuple()
         .map { name, metas, files -> [ [ id: name ], files] }
         .set { ch_additional_bbs }
      PUBLISH_ADDITIONAL_BIGBEDS(ch_additional_bbs)

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          POLISHING [ ISOTOOLS_POLISH ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: sorts reads into pass / duplicates / scraps / trash / retentions / truncations / intrapriming,
      // INFO: published per sample as BED (12_POLISH/BED) and bigBed (12_POLISH/BB)
      ISOTOOLS_POLISH(
          ch_full_length_reads,
          PREPROCESSING.out.reference_transcripts,
          ISOTOOLS_PREPOLISH.out.aparent_plus,
          ISOTOOLS_PREPOLISH.out.aparent_minus,
          PREPROCESSING.out.chrom_sizes,
          PREPROCESSING.out.bigwigs,
          schema,
          ch_versions
      )

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          TRACKING [ TRACKDB ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: optional: writes trackDb and copies the bigBeds to the genome browser server
      if (params.load_track) {
          TRACKDB(
            params.load_track_browser,
            params.global_species_name,
            params.load_track_name,
            ISOTOOLS_POLISH.out.additional_columns,
            ISOTOOLS_POLISH.out.sample,
          )

          LOAD_PASS_TRACK(
            ISOTOOLS_POLISH.out.pass,
            params.load_track_user,
            params.load_track_server,
            params.load_track_target_dir,
            params.load_track_web,
            params.global_species_name,
            ch_versions
          )

          LOAD_DUPLICATES_TRACK(
            ISOTOOLS_POLISH.out.duplicates,
            params.load_track_user,
            params.load_track_server,
            params.load_track_target_dir,
            params.load_track_web,
            params.global_species_name,
            ch_versions
          )

          LOAD_ORPHANS_TRACK(
            ISOTOOLS_POLISH.out.scraps, // INFO: reads the orphan finder set aside
            params.load_track_user,
            params.load_track_server,
            params.load_track_target_dir,
            params.load_track_web,
            params.global_species_name,
            ch_versions
          )

          LOAD_TRASH_TRACK(
            ISOTOOLS_POLISH.out.trash,
            params.load_track_user,
            params.load_track_server,
            params.load_track_target_dir,
            params.load_track_web,
            params.global_species_name,
            ch_versions
          )

          LOAD_RETENTIONS_TRACK(
            ISOTOOLS_POLISH.out.retentions,
            params.load_track_user,
            params.load_track_server,
            params.load_track_target_dir,
            params.load_track_web,
            params.global_species_name,
            ch_versions
          )

          LOAD_TRUNCATIONS_TRACK(
            ISOTOOLS_POLISH.out.truncations,
            params.load_track_user,
            params.load_track_server,
            params.load_track_target_dir,
            params.load_track_web,
            params.global_species_name,
            ch_versions
          )

          LOAD_INTRAPRIMMING_TRACK(
            ISOTOOLS_POLISH.out.intraprimming,
            params.load_track_user,
            params.load_track_server,
            params.load_track_target_dir,
            params.load_track_web,
            params.global_species_name,
            ch_versions
          )

          LOAD_RT_TRACK(
            ISOTOOLS_POLISH.out.rt, // INFO: trackDb lists an rt subtrack; without this its file was never uploaded
            params.load_track_user,
            params.load_track_server,
            params.load_track_target_dir,
            params.load_track_web,
            params.global_species_name,
            ch_versions
          )

          LOAD_FUSIONS_TRACK(
            BEDTOBIGBED_FUSIONS.out.bigbed,
            params.load_track_user,
            params.load_track_server,
            params.load_track_target_dir,
            params.load_track_web,
            params.global_species_name,
            ch_versions
          )

          LOAD_NMD_TRACK(
            BEDTOBIGBED_NMD.out.bigbed, // INFO: the per-sample bigBed, not the per-chromosome BEDs
            params.load_track_user,
            params.load_track_server,
            params.load_track_target_dir,
            params.load_track_web,
            params.global_species_name,
            ch_versions
          )

          ch_versions = ch_versions.mix(TRACKDB.out.versions)
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          VERSIONING
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */


      ch_versions = ch_versions.mix(ISOTOOLS_NMD_FILTER.out.versions)
      ch_versions = ch_versions.mix(JOIN_FUSIONS.out.versions)
      ch_versions = ch_versions.mix(JOIN_NMD.out.versions)
      ch_versions = ch_versions.mix(JOIN_INTRONS.out.versions)
      ch_versions = ch_versions.mix(BEDTOBIGBED_FUSIONS.out.versions)
      ch_versions = ch_versions.mix(BEDTOBIGBED_NMD.out.versions)
      ch_versions = ch_versions.mix(BEDTOBIGBED_INTRONS.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
