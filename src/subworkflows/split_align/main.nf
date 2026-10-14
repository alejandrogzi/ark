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

include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_PBMM2 } from '../../modules/custom/samtools/index/main.nf'

include { ISOTOOLS_SEGMENT as ISOTOOLS_SEGMENT_POLYA } from '../../modules/custom/isotools/segment/main.nf'
include { ISOTOOLS_SEGMENT as ISOTOOLS_SEGMENT_POLYA_FRAGMENTS } from '../../modules/custom/isotools/segment/main.nf'
include { ISOTOOLS_FUSION as ISOTOOLS_FUSION_DETECTOR } from '../../modules/custom/isotools/fusion/main.nf'
include { ISOTOOLS_CIGAR as ISOTOOLS_CIGAR_EXTENSION } from '../../modules/custom/isotools/cigar/main.nf'
include { ISOTOOLS_ADAPTER as ISOTOOLS_REMOVE_ADAPTERS } from '../../modules/custom/isotools/adapter/main.nf'
include { ISOTOOLS_ALIGN as ISOTOOLS_FIND_FRAGMENTS } from '../../modules/custom/isotools/align/main.nf'
include { ISOTOOLS_ALIGN_ACCEPT as ISOTOOLS_ACCEPT_FRAGMENTS } from '../../modules/custom/isotools/accept/main.nf'
include { SAMTOOLS_SELECT_READS as SAMTOOLS_DROP_ACCEPTED } from '../../modules/custom/samtools/select/main.nf'
include { SAMTOOLS_SELECT_READS as SAMTOOLS_KEEP_ACCEPTED } from '../../modules/custom/samtools/select/main.nf'
include { FASTX_CONCAT as FASTX_CONCAT_FRAGMENTS } from '../../modules/custom/fastx/concat/main.nf'


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
      fragments_from_pass1     // bool: iso-align scans the pass-1 BAM instead of the iso-cigar BAM
      fragments_accept         // bool: a realigned read replaces its pass-1 records only if accepted
      fragments_pool           // bool: realign the fragment reads of all chunks of a sample in one task
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

            ch_aligned_bam = ARK_ALIGN.out.bam
            ch_aligned_bai = ARK_ALIGN.out.bai
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

            ch_aligned_bam = MINIMAP2_ALIGN.out.bam
            ch_aligned_bai = MINIMAP2_ALIGN.out.bai
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

          ch_aligned_bam = FLAIR_ALIGN.out.bam
          ch_aligned_bai = FLAIR_ALIGN.out.bai
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
        ch_pass1_bam = ch_aligned_bam.join(ch_aligned_bai) // [ meta, bam, bai ] before iso-cigar
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
          ch_aligned_bam = ch_pass1_bam
        }

        // INFO: iso-cigar extends split primaries across the long gap and drops their SA tag, which moves
        // INFO: them out of iso-align's gap window; scanning the pass-1 BAM keeps them selectable. The id
        // INFO: suffix is kept so fragment names and joins do not change.
        ch_find_input = (fragments_from_pass1 && cigar_extension)
            ? ch_pass1_bam.map { meta, bam, bai -> [ meta + [ id: meta.id + '.extended' ], bam, bai ] }
            : ch_aligned_bam

        // INFO: iso-align needs each chunk BAM plus the chunk FASTA it was aligned from.
        // INFO: key [ sample_id, singleton, chunk ] is set before alignment and survives every meta.clone().
        // INFO: groupTuple + combine instead of join: inputs that collide on the key (e.g. X.fasta.gz and
        // INFO: X.hq.fasta.gz under flnc) both reach --reads instead of being mispaired.
        ch_fastx_gz
            .map { meta, fasta -> [ [ meta.sample_id, meta.singleton, meta.chunk ], fasta ] }
            .groupTuple()
            .set { ch_chunk_reads } // [ key, [ fasta, ... ] ]

        ch_find_input
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
        // INFO: same hint level as pass 1; without annotation hints the fragments get no --junc-bed either.
        // INFO: pooling: a chunk has only ~100 fragment reads, so a per-chunk task spends nearly all its time loading
        // INFO: the minimap2 index and the splice-score table. With fragments_pool the FASTAs of all chunks of a
        // INFO: sample are concatenated and realigned once (meta.chunk = 'pooled'); minimap2 aligns each read
        // INFO: independently, so the alignments are the same as per chunk.
        ch_fragments_fasta = ISOTOOLS_FIND_FRAGMENTS.out.fasta
        if (fragments_pool) {
          ISOTOOLS_FIND_FRAGMENTS.out.fasta
              .map { meta, fasta -> [ [ meta.sample_id, meta.singleton, meta.id ], meta, fasta ] }
              .groupTuple(by: 0)
              .map { key, metas, fastas ->
                  def single = metas[0].singleton ? '.singleton' : ''
                  [ metas[0] + [ chunk: 'pooled', outfile: "${metas[0].id}.pooled${single}.fasta" ], fastas ]
              }
              .set { ch_fragments_pool_in }

          FASTX_CONCAT_FRAGMENTS(ch_fragments_pool_in)
          ch_fragments_fasta = FASTX_CONCAT_FRAGMENTS.out.reads.map { meta, fasta -> [ meta - [ outfile: meta.outfile ], fasta ] }
          ch_versions = ch_versions.mix(FASTX_CONCAT_FRAGMENTS.out.versions)
        }

        ARK_ALIGN_FRAGMENTS(
          ch_fragments_fasta,
          ch_genome_index,
          ch_splice_scores,
          aligner_use_annotation ? ch_reference_transcripts : Channel.value([[:], []])
        )

        ARK_ALIGN_FRAGMENTS.out.bam
            .join(ARK_ALIGN_FRAGMENTS.out.bai)
            .set { ch_fragments_bam }
        ch_versions = ch_versions.mix(ARK_ALIGN_FRAGMENTS.out.versions)

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ACCEPTANCE [ optional ]
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        // INFO: without acceptance a realigned read reaches iso-segment twice (pass-1 + realigned records).
        // INFO: with it, iso-align (accept mode) decides per read, comparing primary records only:
        // INFO:   1. AS(realigned) > AS(pass-1), strictly. AS comes from the pass-1 BAM because iso-cigar
        // INFO:      drops AS/NM on the records it edits; a missing AS on either side rejects the read.
        // INFO:   2. the realigned primary is hq, as iso-segment defines it: identity >= isotools_fragments_min_identity
        // INFO:      (98; (M - mismatches) / (M + I + D), mismatches = NM - indels; no NM = reject) and soft/hard clip
        // INFO:      <= isotools_fragments_max_clip (20 bp) on both ends. Accepted reads are therefore hq downstream too.
        // INFO: a rejected read keeps its base records (iso-cigar BAM). An accepted read loses every base record
        // INFO: (primary + supplementary, SAMTOOLS_DROP_ACCEPTED), and the fragments BAM keeps only accepted
        // INFO: reads (SAMTOOLS_KEEP_ACCEPTED), so each read reaches iso-segment exactly once.
        // INFO: rule text and threshold sweep: isotools iso-align/README.md.
        // INFO: key [ sample_id, singleton, chunk, id without the '.extended' / '.fragments' suffixes ] pairs a
        // INFO: chunk's base, pass-1 and fragments BAMs. The id is needed: under flnc, X.fasta.gz and X.hq.fasta.gz
        // INFO: share [ sample_id, singleton, chunk ]. failOnDuplicate turns any remaining collision into an error.
        if (fragments_accept) {
          def chunkKey = { meta ->
              [ meta.sample_id, meta.singleton, meta.chunk, meta.id.toString().replaceFirst(/(\.extended)?(\.fragments)?$/, '') ]
          }
          // INFO: pooled fragments: one realigned BAM per sample, shared by every chunk of that sample (chunk dropped)
          def poolKey = { meta -> chunkKey(meta).withIndex().findAll { v, i -> i != 2 }.collect { v, i -> v } }

          if (fragments_pool) {
            // INFO: each chunk's pass-1 BAM is paired with its sample's pooled BAM; iso-align decides only the reads
            // INFO: that have a pass-1 primary in that chunk (the others are missing_score there and accepted, or
            // INFO: not, in their own chunk). meta.chunk is set back to the pass-1 chunk for the joins below.
            ch_pass1_bam
                .map { meta, bam, bai -> [ poolKey(meta), meta, bam, bai ] }
                .combine(ch_fragments_bam.map { meta, bam, bai -> [ poolKey(meta), meta, bam, bai ] }, by: 0)
                .map { key, pmeta, pass1_bam, pass1_bai, meta, bam, bai -> [ meta + [ chunk: pmeta.chunk ], pass1_bam, pass1_bai, bam, bai ] }
                .set { ch_accept_input }
          } else {
            ch_pass1_bam
                .map { meta, bam, bai -> [ chunkKey(meta), bam, bai ] }
                .join(ch_fragments_bam.map { meta, bam, bai -> [ chunkKey(meta), meta, bam, bai ] }, failOnDuplicate: true)
                .map { key, pass1_bam, pass1_bai, meta, bam, bai -> [ meta, pass1_bam, pass1_bai, bam, bai ] }
                .set { ch_accept_input }
          }

          ISOTOOLS_ACCEPT_FRAGMENTS(ch_accept_input)
          ch_versions = ch_versions.mix(ISOTOOLS_ACCEPT_FRAGMENTS.out.versions)

          ch_accepted = ISOTOOLS_ACCEPT_FRAGMENTS.out.names.map { meta, names -> [ chunkKey(meta), names ] }

          // INFO: chunks without fragments (no names list) or with no accepted read (empty list) pass through
          // INFO: unchanged; the meta filter drops names lists without a base BAM (remainder fills meta with null).
          ch_aligned_bam
              .map { meta, bam, bai -> [ chunkKey(meta), meta, bam, bai ] }
              .join(ch_accepted, remainder: true, failOnDuplicate: true)
              .filter { key, meta, bam, bai, names -> meta != null }
              .branch { key, meta, bam, bai, names ->
                  filter: names != null && names.size() > 0
                      return [ meta, bam, bai, names ]
                  keep: true
                      return [ meta, bam, bai ]
              }
              .set { ch_base_split }

          SAMTOOLS_DROP_ACCEPTED(ch_base_split.filter)
          ch_aligned_bam = ch_base_split.keep.mix(SAMTOOLS_DROP_ACCEPTED.out.bam)

          // INFO: chunks with no accepted read emit no fragments BAM. Pooled: each chunk keeps its own accepted reads
          // INFO: from the sample's pooled BAM, so the fragments BAMs stay one per chunk downstream.
          if (fragments_pool) {
            ISOTOOLS_ACCEPT_FRAGMENTS.out.names
                .filter { meta, names -> names.size() > 0 }
                .map { meta, names -> [ poolKey(meta), meta, names ] }
                .combine(ch_fragments_bam.map { meta, bam, bai -> [ poolKey(meta), bam, bai ] }, by: 0)
                .map { key, meta, names, bam, bai -> [ meta, bam, bai, names ] }
                .set { ch_fragments_accept }
          } else {
            ch_fragments_bam
                .map { meta, bam, bai -> [ chunkKey(meta), meta, bam, bai ] }
                .join(ch_accepted, failOnDuplicate: true)
                .filter { key, meta, bam, bai, names -> names.size() > 0 }
                .map { key, meta, bam, bai, names -> [ meta, bam, bai, names ] }
                .set { ch_fragments_accept }
          }

          SAMTOOLS_KEEP_ACCEPTED(ch_fragments_accept)
          ch_fragments_bam = SAMTOOLS_KEEP_ACCEPTED.out.bam
          ch_versions = ch_versions.mix(SAMTOOLS_DROP_ACCEPTED.out.versions, SAMTOOLS_KEEP_ACCEPTED.out.versions)
        }


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
