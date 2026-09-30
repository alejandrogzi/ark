/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MINIMAP2_ALIGN } from '../../modules/custom/minimap2/align/main.nf'
include { MINIMAP2_ALIGN as ARK_ALIGN } from '../../modules/custom/minimap2/align/main.nf'
include { MINIMAP2_ALIGN as ARK_ALIGN_FRAGMENTS } from '../../modules/custom/minimap2/align/main.nf'
include { MINIMAP2_ALIGN as FLAIR_ALIGN } from '../../modules/custom/minimap2/align/main.nf'
include { DESALT_ALIGN } from '../../modules/custom/desalt/align/main.nf'
include { ULTRA_ALIGN } from '../../modules/nf-core/ultra/align/main.nf'
include { PBMM2_ALIGN } from '../../modules/nf-core/pbmm2/align/main.nf'

include { FXSPLIT } from '../../modules/custom/fxsplit/main.nf'

include { SAMTOOLS_BAM } from '../../modules/custom/samtools/bam/main.nf'
include { SAMTOOLS_BAM as SAMTOOLS_BAM_PBMM2_ALIGN } from '../../modules/custom/samtools/bam/main.nf'
include { SAMTOOLS_BAM as SAMTOOLS_BAM_DESALT_ALIGN } from '../../modules/custom/samtools/bam/main.nf'
include { SAMTOOLS_BAM as SAMTOOLS_BAM_MINIMAP2_ALIGN } from '../../modules/custom/samtools/bam/main.nf'
include { SAMTOOLS_BAM as SAMTOOLS_BAM_ARK_ALIGN } from '../../modules/custom/samtools/bam/main.nf'
include { SAMTOOLS_BAM as SAMTOOLS_BAM_FLAIR_ALIGN } from '../../modules/custom/samtools/bam/main.nf'
include { SAMTOOLS_BAM as SAMTOOLS_BAM_FRAGMENTS } from '../../modules/custom/samtools/bam/main.nf'

include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_BAM_MULTI_SAMPLE } from '../../modules/custom/samtools/merge/main.nf'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_ARK } from '../../modules/custom/samtools/merge/main.nf'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_MINIMAP2 } from '../../modules/custom/samtools/merge/main.nf'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_ULTRA } from '../../modules/custom/samtools/merge/main.nf'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_DESALT } from '../../modules/custom/samtools/merge/main.nf'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_PBMM2 } from '../../modules/custom/samtools/merge/main.nf'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_FLAIR } from '../../modules/custom/samtools/merge/main.nf'

include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_PBMM2 } from '../../modules/custom/samtools/index/main.nf'

include { ISOTOOLS_SEGMENT as ISOTOOLS_SEGMENT_POLYA } from '../../modules/custom/isotools/segment/main.nf'
include { ISOTOOLS_SEGMENT as ISOTOOLS_SEGMENT_POLYA_FRAGMENTS } from '../../modules/custom/isotools/segment/main.nf'
include { ISOTOOLS_FUSION as ISOTOOLS_FUSION_DETECTOR } from '../../modules/custom/isotools/fusion/main.nf'
include { ISOTOOLS_CIGAR as ISOTOOLS_CIGAR_EXTENSION } from '../../modules/custom/isotools/cigar/main.nf'
include { ISOTOOLS_ADAPTER as ISOTOOLS_REMOVE_ADAPTERS } from '../../modules/custom/isotools/adapter/main.nf'
include { ISOTOOLS_ALIGN as ISOTOOLS_FIND_FRAGMENTS } from '../../modules/custom/isotools/align/main.nf'


include { COLLAPSE as COLLAPSE_TWINS } from '../../modules/custom/collapse/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SPLIT_ALIGN_CLEAN_CHUNKS {
    take:
      ch_reads                 // [ meta, reads ]
      ch_genome                // [ genome ]
      ch_genome_index          // [ meta, index ]
      ch_reference_transcripts // [ meta, bed ]
      prefix                   // string
      ch_splice_scores         // [ meta, scores ]
      cluster_mode             // string [ per_sample, multi_sample, both ]
      entrypoint               // string [ isoseq, map ]
      aligner                  // string [ minimap2, ultra, pbmm2, desalt, ark, flair ]
      aligner_use_annotation   // bool
      remove_adapters          // bool
      collapse_twins           // bool
      cigar_extension          // bool
      do_second_pass           // bool
      ch_versions              // [ meta, versions.yml ]

    main:
      ch_reads = ch_reads.map { meta, reads ->
          [ meta + [ sample_id: meta.sample_id ?: meta.id ], reads ]
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          CHUNKING [ isoseq, map ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      FXSPLIT(ch_reads)
      FXSPLIT.out.fastx_gz
          .flatMap {
              meta, fa ->
              def fas = fa instanceof List ? fa : [fa]
              fas.collect { it ->
                  def parts = it.baseName.split('_')
                  def chunk = parts.size() > 2 ? parts[2] : 0
                  [ meta + [ chunk: chunk ], it ]
              }
          }
          .set { ch_fastx_gz }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ALIGNMENT [ ark, mm2, ultra, desalt, pbmm2, flair ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ch_aligned_bam = Channel.empty()
      ch_aligned_bai = Channel.empty()

      switch (params.aligner) {
        case 'ark':
            if (aligner_use_annotation) {
                ARK_ALIGN(
                  ch_fastx_gz,
                  ch_genome_index,
                  ch_splice_scores,
                  ch_reference_transcripts
                )
            } else {
                ARK_ALIGN(
                  ch_fastx_gz,
                  ch_genome_index,
                  ch_splice_scores,
                  Channel.value([[:], []])
                )
            }

            SAMTOOLS_BAM_ARK_ALIGN(ARK_ALIGN.out.sam)
            ch_aligned_bam = SAMTOOLS_BAM_ARK_ALIGN.out.bam
            ch_aligned_bai = SAMTOOLS_BAM_ARK_ALIGN.out.bai
            ch_versions = ch_versions.mix(SAMTOOLS_BAM_ARK_ALIGN.out.versions)
            ch_versions = ch_versions.mix(ARK_ALIGN.out.versions)
          break

        case 'mm2':
            if (aligner_use_annotation) {
                MINIMAP2_ALIGN(
                  ch_fastx_gz,
                  ch_genome_index,
                  ch_splice_scores,
                  ch_reference_transcripts
                )
            } else {
                MINIMAP2_ALIGN(
                  ch_fastx_gz,
                  ch_genome_index,
                  ch_splice_scores,
                  Channel.value([[:], []])
                )
            }

            SAMTOOLS_BAM_MINIMAP2_ALIGN(MINIMAP2_ALIGN.out.sam)
            ch_aligned_bam = SAMTOOLS_BAM_MINIMAP2_ALIGN.out.bam
            ch_aligned_bai = SAMTOOLS_BAM_MINIMAP2_ALIGN.out.bai
            ch_versions = ch_versions.mix(SAMTOOLS_BAM_MINIMAP2_ALIGN.out.versions)
            ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
          break

        case 'pbmm2':
            PBMM2_ALIGN(
              ch_fastx_gz,
              ch_genome_index
            )

            SAMTOOLS_INDEX_PBMM2(PBMM2_ALIGN.out.bam)
            ch_aligned_bam = SAMTOOLS_INDEX_PBMM2.out.bam
            ch_aligned_bai = SAMTOOLS_INDEX_PBMM2.out.bai
            ch_versions = ch_versions.mix(PBMM2_ALIGN.out.versions)
          break

        case 'desalt':
            DESALT_ALIGN(
              ch_fastx_gz,
              ch_genome_index,
              Channel.value([[:], []])
            )

            SAMTOOLS_BAM_DESALT_ALIGN(DESALT_ALIGN.out.sam)
            ch_aligned_bam = SAMTOOLS_BAM_DESALT_ALIGN.out.bam
            ch_aligned_bai = SAMTOOLS_BAM_DESALT_ALIGN.out.bai
            ch_versions = ch_versions.mix(SAMTOOLS_BAM_DESALT_ALIGN.out.versions)
            ch_versions = ch_versions.mix(DESALT_ALIGN.out.versions)
          break

        case 'ultra':
            ULTRA_ALIGN(
              ch_fastx_gz,
              ch_genome.map { genome -> [ [id:genome.baseName], genome ] },
              ch_genome_index,
            )

            ch_aligned_bam = ULTRA_ALIGN.out.bam
            ch_aligned_bai = ULTRA_ALIGN.out.bai
            ch_versions = ch_versions.mix(ULTRA_ALIGN.out.versions)
          break

        case 'flair':
          FLAIR_ALIGN(
            ch_fastx_gz,
            ch_genome_index,
            Channel.value([[:], []]),
            Channel.value([[:], []])
          )

          SAMTOOLS_BAM_FLAIR_ALIGN(FLAIR_ALIGN.out.sam)
          ch_aligned_bam = SAMTOOLS_BAM_FLAIR_ALIGN.out.bam
          ch_aligned_bai = SAMTOOLS_BAM_FLAIR_ALIGN.out.bai
          ch_versions = ch_versions.mix(SAMTOOLS_BAM_FLAIR_ALIGN.out.versions)
          ch_versions = ch_versions.mix(FLAIR_ALIGN.out.versions)
        break

        default:
          error """
          ERROR: Unknown aligner: '${params.aligner}'.
          Valid aligners are: minimap2, pbmm2, desalt, ultra, flair
          """.stripIndent()
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          CLUSTERING [ isoseq, map ] + BRANCHING
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ch_pooled_reads = Channel.empty()
      if (entrypoint == "isoseq") {
        // INFO: reads were already merged in isoseq clustering
        ch_pooled_reads = ch_reads
      } else if (entrypoint == "map") {
        // INFO: fastq reads need to be merged if multi-sample
        if (cluster_mode == "per_sample") {
          ch_pooled_reads = ch_reads
        } else if (cluster_mode == "multi_sample") {
          ch_aligned_bam
              .map { meta, bam -> bam }
              .collect()
              .map { bams -> [
                [
                  id: prefix,
                  single_end: false,
                  singleton: false,
                  chunk: 0
                ], bams ] }
              .set { ch_joined_bam }

          switch (params.aligner) {
            case 'ark':
              SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_ARK(ch_joined_bam)
              ch_aligned_bam = SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_ARK.out.bam
              ch_aligned_bai = SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_ARK.out.bai
              ch_versions = ch_versions.mix(SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_ARK.out.versions)
            break

            case 'mm2':
              SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_MINIMAP2(ch_joined_bam)
              ch_aligned_bam = SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_MINIMAP2.out.bam
              ch_aligned_bai = SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_MINIMAP2.out.bai
              ch_versions = ch_versions.mix(SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_MINIMAP2.out.versions)
            break

            case 'ultra':
              SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_ULTRA(ch_joined_bam)
              ch_aligned_bam = SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_ULTRA.out.bam
              ch_aligned_bai = SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_ULTRA
              ch_versions = ch_versions.mix(SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_ULTRA.out.versions)
            break

            case 'desalt':
              SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_DESALT(ch_joined_bam)
              ch_aligned_bam = SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_DESALT.out.bam
              ch_aligned_bai = SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_DESALT.out.bai
              ch_versions = ch_versions.mix(SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_DESALT.out.versions)
            break

            case 'pbmm2':
              SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_PBMM2(ch_joined_bam)
              ch_aligned_bam = SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_PBMM2.out.bam
              ch_aligned_bai = SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_PBMM2.out.bai
              ch_versions = ch_versions.mix(SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_PBMM2.out.versions)
            break

            case 'flair':
              SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_FLAIR(ch_joined_bam)
              ch_aligned_bam = SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_FLAIR.out.bam
              ch_aligned_bai = SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_FLAIR.out.bai
              ch_versions = ch_versions.mix(SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_FLAIR.out.versions)
            break

            default:
              error """
              ERROR: Unknown aligner: '${params.aligner}'.
              Valid aligners are: ark, mm2, ultra, desalt, pbmm2, flair
              """.stripIndent()
              System.exit(1)
          }

          ch_reads
            .map { meta, reads -> reads }
            .collect()
            .map { reads -> [ [ id: 'pooled.reads' ], reads ] }
            .set { ch_pooled_reads }
        } else if (cluster_mode == "both") {
          ch_aligned_bam
              .map { meta, bam -> bam }
              .collect()
              .map { bams -> [ [ id: prefix ], bams ] }
              .set { ch_joined_bam }

          switch (params.aligner) {
            case 'ark':
              SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_ARK(ch_joined_bam)
              ch_aligned_bam = SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_ARK.out.bam
              ch_aligned_bai = SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_ARK.out.bai
              ch_versions = ch_versions.mix(SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_ARK.out.versions)
            break

            case 'mm2':
              SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_MINIMAP2(ch_joined_bam)
              ch_aligned_bam = SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_MINIMAP2.out.bam
              ch_aligned_bai = SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_MINIMAP2.out.bai
              ch_versions = ch_versions.mix(SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_MINIMAP2.out.versions)
            break

            case 'ultra':
              SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_ULTRA(ch_joined_bam)
              ch_aligned_bam = SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_ULTRA.out.bam
              ch_aligned_bai = SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_ULTRA
              ch_versions = ch_versions.mix(SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_ULTRA.out.versions)
            break

            case 'desalt':
              SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_DESALT(ch_joined_bam)
              ch_aligned_bam = SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_DESALT.out.bam
              ch_aligned_bai = SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_DESALT.out.bai
              ch_versions = ch_versions.mix(SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_DESALT.out.versions)
            break

            case 'pbmm2':
              SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_PBMM2(ch_joined_bam)
              ch_aligned_bam = SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_PBMM2.out.bam
              ch_aligned_bai = SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_PBMM2.out.bai
              ch_versions = ch_versions.mix(SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_PBMM2.out.versions)
            break

            case 'flair':
              SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_FLAIR(ch_joined_bam)
              ch_aligned_bam = SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_FLAIR.out.bam
              ch_aligned_bai = SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_FLAIR.out.bai
              ch_versions = ch_versions.mix(SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_FLAIR.out.versions)
            break

            default:
              error """
              ERROR: Unknown aligner: '${params.aligner}'.
              Valid aligners are: ark, mm2, ultra, desalt, pbmm2, flair
              """.stripIndent()
              System.exit(1)
          }

          ch_reads
            .map { meta, reads -> reads }
            .collect()
            .map { reads -> [ [ id: 'pooled.reads' ], reads ] }
            .set { ch_pooled_reads }
        }
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ADAPTER REMOVAL [ optional ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      if (remove_adapters) {
        ISOTOOLS_REMOVE_ADAPTERS(
          ch_aligned_bam,
          ch_aligned_bai,
        )
        ch_aligned_bam = ISOTOOLS_REMOVE_ADAPTERS.out.bam
        ch_aligned_bai = ISOTOOLS_REMOVE_ADAPTERS.out.bai
        ch_versions = ch_versions.mix(ISOTOOLS_REMOVE_ADAPTERS.out.versions)
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          CIGAR EXTENSION + FRAGMENT DETECTION [ ark only ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ch_fragments_bam = Channel.empty()
      if (do_second_pass && params.aligner == 'ark') {
        if (cigar_extension) {
          ISOTOOLS_CIGAR_EXTENSION(
            ch_aligned_bam,
            ch_aligned_bai,
            ch_genome.map { genome -> [ [id:genome.baseName], genome ] },
            ch_reference_transcripts
          )

          // WARN: replacing ch_aligned_bam with [ meta, bam, bai ]
          ch_aligned_bam = ISOTOOLS_CIGAR_EXTENSION.out.extended

          ISOTOOLS_FIND_FRAGMENTS(
            ISOTOOLS_CIGAR_EXTENSION.out.extended,
            ch_pooled_reads
          )

          ch_versions = ch_versions.mix(ISOTOOLS_CIGAR_EXTENSION.out.versions)
        } else {
          ISOTOOLS_FIND_FRAGMENTS(
            ch_aligned_bam.join(ch_aligned_bai),
            ch_pooled_reads
          )

          // WARN: replacing ch_aligned_bam with [ meta, bam, bai ]
          ch_aligned_bam = ch_aligned_bam.join(ch_aligned_bai)
          ch_versions = ch_versions.mix(ISOTOOLS_FIND_FRAGMENTS.out.versions)
        }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            RE-ALIGNMENT [ ark only ]
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        ARK_ALIGN_FRAGMENTS(
          ISOTOOLS_FIND_FRAGMENTS.out.fasta,
          ch_genome_index,
          ch_splice_scores,
          ch_reference_transcripts
        )

        SAMTOOLS_BAM_FRAGMENTS(ARK_ALIGN_FRAGMENTS.out.sam)
        SAMTOOLS_BAM_FRAGMENTS.out.bam
            .join(SAMTOOLS_BAM_FRAGMENTS.out.bai)
            .set { ch_fragments_bam }


      } else {
        ch_aligned_bam = ch_aligned_bam.join(ch_aligned_bai)
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          POLYA SEGMENTATION 
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ISOTOOLS_SEGMENT_POLYA(ch_aligned_bam) // INFO: polyA tails + cigar 
      ISOTOOLS_SEGMENT_POLYA_FRAGMENTS(ch_fragments_bam) // INFO: fragments + cigar

      ISOTOOLS_SEGMENT_POLYA.out.hq_bed
          .mix(ISOTOOLS_SEGMENT_POLYA_FRAGMENTS.out.hq_bed)
          .set { ch_aligned_segmented }

      ch_aligned_segmented
          .flatMap { meta, bed ->   // ← no extra braces
              def beds = bed instanceof List ? bed : [bed]
              beds.collect { it ->
                  [ meta + [ chr: it.name.split('@')[0] ], it ]
              }
          }
          .map { meta, bed ->
              [ [ meta.chr, meta.sample_id ], meta, bed ]
          }
          .groupTuple(by: 0)
          .map { key, metas, beds ->
              def meta = metas[0]
              def group_meta = [ id: meta.sample_id, single_end: true, chr: meta.chr ]
              [ group_meta, beds ]
          }
          .set { ch_aligned_segmented_hq_per_chr }
  
      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          COLLAPSE [ optional ] + FUSION DETECTION
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ch_aligned_segmented_collapsed = Channel.empty()
      if (collapse_twins) {
        COLLAPSE_TWINS(ch_aligned_segmented_hq_per_chr)
        ch_aligned_segmented_collapsed = COLLAPSE_TWINS.out.collapsed
      } else {
        ch_aligned_segmented_collapsed = ch_aligned_segmented_hq_per_chr
      }

      ISOTOOLS_FUSION_DETECTOR(
        ch_aligned_segmented_collapsed, 
        ch_reference_transcripts
      )

      ch_versions = ch_versions.mix(FXSPLIT.out.versions)
      ch_versions = ch_versions.mix(ISOTOOLS_SEGMENT_POLYA.out.versions)
      ch_versions = ch_versions.mix(ISOTOOLS_FUSION_DETECTOR.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        OUTPUT CHANNELS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    emit:
      reads    = ISOTOOLS_FUSION_DETECTOR.out.free_fusion
      fusions  = ISOTOOLS_FUSION_DETECTOR.out.fusion
      versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
