/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PBCCS } from '../../modules/nf-core/pbccs/main.nf'

include { PBTK_PBINDEX as PBINDEX } from '../../modules/nf-core/pbtk/pbindex/main.nf'

include { PBTK_PBMERGE as PBMERGE } from '../../modules/nf-core/pbtk/pbmerge/main.nf'

include { PBSKERA_SPLIT } from '../../modules/custom/pbskera/split/main.nf'
include { WGET as WGET_SKERA_PRIMERS } from '../../modules/nf-core/wget/main.nf'

include { LIMA } from '../../modules/nf-core/lima/main.nf'
include { ISOSEQ_REFINE } from '../../modules/nf-core/isoseq/refine/main.nf'
include { ISOSEQ_CLUSTER2 } from '../../modules/custom/isoseq/cluster2/main.nf'
include { ISOSEQ_CLUSTER2 as ISOSEQ_CLUSTER2_MULTI_SAMPLE } from '../../modules/custom/isoseq/cluster2/main.nf'

include { BAM_TO_FA } from '../../modules/custom/bamtofa/main.nf'
include { BAM_TO_FA as BAM_TO_FA_MULTI_SAMPLE } from '../../modules/custom/bamtofa/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
 
workflow ISOSEQ {
    take:
      global_input_dir       // path
      global_primers         // path
      ccs_chunk              // int
      isoseq_cluster2_mode   // string
      prefix                 // string
      entrypoint             // [ subreads, ccs, refine ] (flnc unreachable)
      is_kinnex_library      // bool

    main:
      ch_versions = Channel.empty()
      ch_primers = Channel.value(file(global_primers, checkIfExists: true))

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          CHANNELING/INDEXING [ SUBREADS, CCS, LIMA OUTPUTS ]
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      Channel
          .fromPath("${global_input_dir}/*.bam", checkIfExists: true)
          .map { bam ->
              def pbi = bam + '.pbi'
              return [
                  [
                      id:         bam.baseName,
                      single_end: true,
                      indexed:    file(pbi).exists()
                  ],
                  bam,
                  file(pbi).exists() ? file(pbi) : []   // placeholder if missing
              ]
          }
          .set { ch_bam }


      ch_bam
          .branch {
              indexed:     it[0].indexed
              not_indexed: true
          }
          .set { ch_bam_branched }


      PBINDEX(
        ch_bam_branched.not_indexed
        .map { 
          meta, bam, _pbi -> [ meta, bam ] 
        }
      )

       PBINDEX.out.pbi
          .join(
              ch_bam_branched.not_indexed.map { meta, bam, _pbi -> [ meta, bam ] },
              by: 0   // join on meta
          )
          .map { meta, pbi, bam ->
              def meta_updated = meta + [ indexed: true ]
              [ meta_updated, bam, pbi ]
          }.set { ch_bam_reindexed }

      ch_bam = ch_bam_branched.indexed.mix(ch_bam_reindexed)
      ch_versions = ch_versions.mix(PBINDEX.out.versions)

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ENTRYPOINT BRANCHING  [ SUBREADS, CCS, REFINE ]
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */
  
      ch_ccs_bams = Channel.empty()
      switch (entrypoint) {
        case 'subreads':

          ch_bam
              .combine(Channel.of(1..ccs_chunk))   // INFO: cartesian product: N_bam × chunk combos
              .map { meta, bam, pbi, chunk_idx ->
                  [ meta + [chunk: chunk_idx], bam, pbi ]
              }
              .set { ch_chunks }

          PBCCS(ch_chunks, ccs_chunk) // INFO: generate CCS from raw reads
          PBCCS.out.bam // INFO: update meta: update id (+chunkX) and store former id
          .map {
              def chunk   = it[0].chunk
              def parent  = it[0].id
              def child   = it[0].id + "." + chunk
              return [ [id:child, parent:parent, single_end:true, chunk:chunk], it[1] ]
          }
          .set { ch_pbccs_bam_updated }

          // INFO: group all chunks belonging to the same parent sample
          ch_pbccs_bam_updated
              .map { meta, bam -> [ meta.parent, meta, bam ] }   // INFO: key by parent
              .groupTuple(by: 0, size: ccs_chunk)              // INFO: wait for all chunks
              .map { parent, metas, bams ->
                  // Reconstruct a clean meta for the merged output
                  def meta_merged = [ id: parent, single_end: true ]
                  [ meta_merged, bams ]                           // INFO: bams is now a List
              }
              .set { ch_pbccs_merged }

          PBMERGE(ch_pbccs_merged) // INFO: merge chunks
          ch_ccs_bams = PBMERGE.out.bam.join(PBMERGE.out.pbi)

          ch_versions = ch_versions.mix(PBMERGE.out.versions)
          ch_versions = ch_versions.mix(PBCCS.out.versions)

        break

        case 'ccs':
          ch_ccs_bams = ch_bam
        break

        case 'refine':
          // LIMA BAMs are already demultiplexed and primer-trimmed.
        break

        default:
          error """
          ERROR: Unknown entrypoint -> options at this step are: subreads, ccs, refine
          """.stripIndent()
          System.exit(1)
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          SKERA [ DEMULTIPLEX KINNEX ARRAY ]
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ch_skera_demux_bams = Channel.empty()
      if (entrypoint != 'refine' && is_kinnex_library) {
        WGET_SKERA_PRIMERS(
          Channel.value(
            params.skera_kinnex_primers
          ).map { url -> [ [id : url.tokenize('/')[-1]], url ] }
        )

        ch_skera_demux_primers = WGET_SKERA_PRIMERS.out.outfile

        PBSKERA_SPLIT(
          ch_ccs_bams.map{ meta, bam, pbi -> [ meta, bam ] },
          ch_skera_demux_primers.map { meta, primers -> primers }
        )
        ch_skera_demux_bams = PBSKERA_SPLIT.out.bam.join(PBSKERA_SPLIT.out.pbi)
      } else {
        ch_skera_demux_bams = ch_ccs_bams
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          LIMA [ REMOVE PRIMERS ]
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ch_lima_out_bams = Channel.empty()
      if (entrypoint == 'refine') {
          ch_lima_out_bams = ch_bam
      } else {
          LIMA(ch_skera_demux_bams, ch_primers)
          ch_lima_out_bams = LIMA.out.bam
              .flatMap { meta, bams ->
                  def bam_files = bams instanceof List ? bams : [bams]
                  bam_files.collect { bam ->
                      tuple(meta, bam, file("${bam}.pbi", checkIfExists: true))
                  }
              }
          ch_versions = ch_versions.mix(LIMA.out.versions)
      }

      // Keep one item per primer/barcode pair for both fresh runs and checkpoints.
      ch_lima_out_bams = ch_lima_out_bams.map { meta, bam, pbi ->
          def matcher = bam.baseName =~ /^(.*)\.([^.]+--[^.]+)$/
          if (!matcher.matches()) error "Unexpected LIMA BAM name: ${bam.name}"

          def pair = matcher.group(2)
          def parent = entrypoint == 'refine'
              ? matcher.group(1).replaceFirst(/_fl$/, '')
              : meta.id
          tuple(meta + [ parent_id: parent, id: "${parent}.${pair}", barcode: pair ], bam, pbi)
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ISOSEQ_REFINE [ DISCARD CCS WITHOUT POLYA TAILS ]
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ISOSEQ_REFINE(ch_lima_out_bams, ch_primers) // INFO: discard CCS without polyA tails

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ISOSEQ_CLUSTER2 [ CLUSTER READS ]
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // Refine each sample once; only clustering needs separate and pooled branches.
      ch_pbccs_merged_flnc_clustered_fa = Channel.empty()
      if (isoseq_cluster2_mode in ['per_sample', 'both']) {
        ISOSEQ_CLUSTER2(ISOSEQ_REFINE.out.bam) // INFO: cluster reads
        BAM_TO_FA(ISOSEQ_CLUSTER2.out.bam)

        ch_pbccs_merged_flnc_clustered_fa  = ch_pbccs_merged_flnc_clustered_fa.mix(BAM_TO_FA.out.singletons)
        ch_pbccs_merged_flnc_clustered_fa  = ch_pbccs_merged_flnc_clustered_fa.mix(BAM_TO_FA.out.hq)

        ch_versions = ch_versions.mix(ISOSEQ_CLUSTER2.out.versions)
        ch_versions = ch_versions.mix(BAM_TO_FA.out.versions)
      }

      if (isoseq_cluster2_mode in ['multi_sample', 'both']) {
        // cluster2 accepts a FOFN, so pooling does not need an intermediate BAM merge.
        ISOSEQ_REFINE.out.bam
          .map { meta, bam -> bam  }
          .collect()
          .map { bams -> [ [ id: prefix, single_end: true ], bams ] }
          .set { ch_pooled_bams }

        ISOSEQ_CLUSTER2_MULTI_SAMPLE(ch_pooled_bams)
        BAM_TO_FA_MULTI_SAMPLE(ISOSEQ_CLUSTER2_MULTI_SAMPLE.out.bam)

        ch_pbccs_merged_flnc_clustered_fa = ch_pbccs_merged_flnc_clustered_fa.mix(BAM_TO_FA_MULTI_SAMPLE.out.singletons)
        ch_pbccs_merged_flnc_clustered_fa = ch_pbccs_merged_flnc_clustered_fa.mix(BAM_TO_FA_MULTI_SAMPLE.out.hq)

        ch_versions = ch_versions.mix(BAM_TO_FA_MULTI_SAMPLE.out.versions)
        ch_versions = ch_versions.mix(ISOSEQ_CLUSTER2_MULTI_SAMPLE.out.versions)
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          VERSIONING
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ch_versions = ch_versions.mix(ISOSEQ_REFINE.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        OUTPUT CHANNELS
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    emit:
        reads   = ch_pbccs_merged_flnc_clustered_fa
        versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
