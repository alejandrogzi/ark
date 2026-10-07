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

include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_PBMM2 } from '../../modules/custom/samtools/index/main.nf'

include { ISOTOOLS_SEGMENT as ISOTOOLS_SEGMENT_POLYA } from '../../modules/custom/isotools/segment/main.nf'
include { ISOTOOLS_SEGMENT as ISOTOOLS_SEGMENT_POLYA_FRAGMENTS } from '../../modules/custom/isotools/segment/main.nf'
include { ISOTOOLS_FUSION as ISOTOOLS_FUSION_DETECTOR } from '../../modules/custom/isotools/fusion/main.nf'
include { ISOTOOLS_CIGAR as ISOTOOLS_CIGAR_EXTENSION } from '../../modules/custom/isotools/cigar/main.nf'
include { ISOTOOLS_ADAPTER as ISOTOOLS_REMOVE_ADAPTERS } from '../../modules/custom/isotools/adapter/main.nf'
include { ISOTOOLS_ALIGN as ISOTOOLS_FIND_FRAGMENTS } from '../../modules/custom/isotools/align/main.nf'


include { COLLAPSE as COLLAPSE_TWINS } from '../../modules/custom/collapse/main.nf'

include { RECONSTRUCT } from '../reconstruct/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SPLIT_ALIGN_CLEAN_CHUNKS {
    take:
      ch_reads                 // [ meta, reads ]; one hq or singleton FASTA per sample, meta.singleton set
      ch_genome                // [ genome ]
      ch_genome_index          // [ meta, index ]
      ch_reference_transcripts // [ meta, bed ]
      ch_splice_scores         // [ meta, scores ]
      aligner                  // string [ ark, mm2, ultra, desalt, pbmm2, flair ]
      aligner_use_annotation   // bool
      remove_adapters          // bool
      collapse_twins           // bool
      cigar_extension          // bool
      do_second_pass           // bool
      reconstruct_engine       // string [ chain, none ]
      ch_versions              // [ meta, versions.yml ]

    main:
      // INFO: sample_id groups hq + singleton files of one sample later on (per-chromosome
      // INFO: grouping, fragment reads). flnc sets it from the file name; Iso-Seq reads fall back to id.
      ch_reads = ch_reads.map { meta, reads ->
          [ meta + [ sample_id: meta.sample_id ?: meta.id ], reads ]
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          CHUNKING [ FXSPLIT ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: each FASTA is split into chunks named tmp_chunk_<N>_<prefix>.fasta.gz;
      // INFO: '_' field 2 is the chunk number. ch_fastx_gz: [ meta + [ chunk: N ], chunk.fasta.gz ]
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

      // INFO: one alignment per chunk; meta is the chunk meta (id, sample_id, singleton, chunk).
      // INFO: ch_aligned_bam: [ meta, bam ], ch_aligned_bai: [ meta, bai ]. They stay separate because
      // INFO: the pass-1 modules take bam and bai as two inputs; both come from one task, so order pairs them.
      ch_aligned_bam = Channel.empty()
      ch_aligned_bai = Channel.empty()

      switch (aligner) {
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
          ERROR: Unknown aligner: '${aligner}'.
          Valid aligners are: ark, mm2, ultra, desalt, pbmm2, flair
          """.stripIndent()
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ADAPTER REMOVAL [ optional ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: removes adapter sequence from soft-clipped read ends; keeps the [ meta, bam ] + [ meta, bai ] shape
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

      // INFO: second pass. iso-cigar rescues missed 3' junctions (optional); iso-align then picks reads
      // INFO: whose pass-1 split alignment hints at an intron over the pass-1 cap, to realign them below.
      ch_fragments_bam = Channel.empty()
      if (do_second_pass && aligner == 'ark') {
        // WARN: replacing ch_aligned_bam with [ meta, bam, bai ] on both paths
        if (cigar_extension) {
          ISOTOOLS_CIGAR_EXTENSION(
            ch_aligned_bam,
            ch_aligned_bai,
            ch_genome.map { genome -> [ [id:genome.baseName], genome ] },
            ch_reference_transcripts
          )

          ch_aligned_bam = ISOTOOLS_CIGAR_EXTENSION.out.extended
          ch_versions = ch_versions.mix(ISOTOOLS_CIGAR_EXTENSION.out.versions)
        } else {
          ch_aligned_bam = ch_aligned_bam.join(ch_aligned_bai)
        }

        // INFO: iso-align needs each chunk BAM plus the chunk FASTA it was aligned from.
        // INFO: key [ sample_id, singleton, chunk ] is set before alignment and survives every meta.clone().
        // INFO: groupTuple + combine instead of join: inputs that collide on the key (e.g. X.fasta.gz and
        // INFO: X.hq.fasta.gz under flnc) both reach --reads instead of being mispaired.
        ch_fastx_gz
            .map { meta, fasta -> [ [ meta.sample_id, meta.singleton, meta.chunk ], fasta ] }
            .groupTuple()
            .set { ch_chunk_reads } // [ key, [ fasta, ... ] ]

        ch_aligned_bam
            .map { meta, bam, bai -> [ [ meta.sample_id, meta.singleton, meta.chunk ], meta, bam, bai ] }
            .combine(ch_chunk_reads, by: 0)
            .multiMap { key, meta, bam, bai, reads ->
                bam:   [ meta, bam, bai ]
                reads: [ meta, reads ]
            }
            .set { ch_find_fragments }

        ISOTOOLS_FIND_FRAGMENTS(
          ch_find_fragments.bam,
          ch_find_fragments.reads
        )
        ch_versions = ch_versions.mix(ISOTOOLS_FIND_FRAGMENTS.out.versions)

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            RE-ALIGNMENT [ ark only ]
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        // INFO: fragment reads [ meta + id '.fragments', fasta ] realigned with the larger intron cap
        // INFO: (minimap2_realign_intron_size); ch_fragments_bam: [ meta, bam, bai ]
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
        // WARN: replacing ch_aligned_bam with [ meta, bam, bai ]
        ch_aligned_bam = ch_aligned_bam.join(ch_aligned_bai)
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          POLYA SEGMENTATION 
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ISOTOOLS_SEGMENT_POLYA(ch_aligned_bam) // INFO: polyA tails + cigar 
      ISOTOOLS_SEGMENT_POLYA_FRAGMENTS(ch_fragments_bam) // INFO: fragments + cigar

      // INFO: hq_bed: [ meta, [ <chr>@*.hq.bed, ... ] ], one bed per chromosome per chunk.
      // INFO: fragment beds are mixed in first so realigned reads land in the same per-chromosome groups.
      ISOTOOLS_SEGMENT_POLYA.out.hq_bed
          .mix(ISOTOOLS_SEGMENT_POLYA_FRAGMENTS.out.hq_bed)
          .set { ch_aligned_segmented }

      // INFO: regroup by [ chr, sample_id ]: all chunks, hq + singleton, of one sample on one chromosome.
      // INFO: ch_aligned_segmented_hq_per_chr: [ [ id: sample_id, single_end: true, chr, clustered ], [ beds ] ]
      ch_aligned_segmented
          .flatMap { meta, bed ->
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
              def group_meta = [ id: meta.sample_id, single_end: true, chr: meta.chr, clustered: metas.any { it.clustered } ]
              [ group_meta, beds ]
          }
          .set { ch_aligned_segmented_hq_per_chr }
  
      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          COLLAPSE [ optional ] + FUSION DETECTION
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: same per-chromosome shape in and out; fusion emits free (non-fusion) reads and fusions
      // INFO: twin collapse only applies without a reconstruction engine (collapse chain supersedes it)
      if (collapse_twins && reconstruct_engine != 'none') {
        log.warn "collapse_shrink_twins is ignored with reconstruct_engine ${reconstruct_engine}"
      }
      ch_aligned_segmented_collapsed = Channel.empty()
      if (collapse_twins && reconstruct_engine == 'none') {
        COLLAPSE_TWINS(ch_aligned_segmented_hq_per_chr)
        ch_aligned_segmented_collapsed = COLLAPSE_TWINS.out.collapsed
      } else {
        ch_aligned_segmented_collapsed = ch_aligned_segmented_hq_per_chr
      }

      ISOTOOLS_FUSION_DETECTOR(
        ch_aligned_segmented_collapsed, 
        ch_reference_transcripts
      )

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          TRANSCRIPT RECONSTRUCTION [ RECONSTRUCT: chain ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: fusion detection decides by read ratios, so it reads reads; both of its outputs are then
      // INFO: collapsed into models per [ sample_id, chr ] (meta.class keeps the two apart)
      ch_free_reads   = ISOTOOLS_FUSION_DETECTOR.out.free_fusion
      ch_fusion_reads = ISOTOOLS_FUSION_DETECTOR.out.fusion
      ch_support      = Channel.empty()
      if (reconstruct_engine != 'none') {
        RECONSTRUCT(
          ch_free_reads.map { meta, bed -> [ meta + [ class: 'free' ], bed ] }
            .mix(ch_fusion_reads.map { meta, bed -> [ meta + [ class: 'fusion' ], bed ] }),
          ch_reference_transcripts,
          reconstruct_engine
        )

        RECONSTRUCT.out.models
          .map { meta, bed -> [ meta.findAll { k, v -> k != 'class' }, bed, meta.class ] }
          .branch { meta, bed, cls ->
            free:   cls == 'free'
            fusion: true
          }
          .set { ch_models }

        ch_free_reads   = ch_models.free.map { meta, bed, cls -> [ meta, bed ] }
        ch_fusion_reads = ch_models.fusion.map { meta, bed, cls -> [ meta, bed ] }
        ch_support      = RECONSTRUCT.out.support
        ch_versions = ch_versions.mix(RECONSTRUCT.out.versions)
      }

      ch_versions = ch_versions.mix(FXSPLIT.out.versions)
      ch_versions = ch_versions.mix(ISOTOOLS_SEGMENT_POLYA.out.versions)
      ch_versions = ch_versions.mix(ISOTOOLS_FUSION_DETECTOR.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        OUTPUT CHANNELS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    emit:
      reads    = ch_free_reads                              // models (chain) or reads (none)
      fusions  = ch_fusion_reads                            // fusion models (chain) or reads (none)
      evidence = ISOTOOLS_FUSION_DETECTOR.out.free_fusion   // uncollapsed free reads (intron frequencies)
      support  = ch_support                                 // [ meta, support.tsv, counts.tsv ]
      versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
