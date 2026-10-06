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

include { SAMTOOLS_FASTA } from '../../modules/custom/samtools/fasta/main.nf'
include { CDHIT_EST } from '../../modules/custom/cdhit/est/main.nf'
include { RATTLE } from '../../modules/custom/rattle/main.nf'

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
      cluster_mode           // string [ per_sample, multi_sample, both ]
      cluster_engine         // string [ isoseq, cdhit, rattle ]
      prefix                 // string
      entrypoint             // [ subreads, ccs, refine, cluster ] (flnc unreachable)
      is_kinnex_library      // bool

    main:
      // INFO: stage each entrypoint enters at (input dir in brackets):
      //   subreads [raw subreads]      -> PBCCS -> SKERA (Kinnex only) -> LIMA -> ISOSEQ_REFINE -> ISOSEQ_CLUSTER2
      //   ccs      [CCS BAMs]          -> SKERA (Kinnex only) -> LIMA -> ISOSEQ_REFINE -> ISOSEQ_CLUSTER2
      //   refine   [02_LIMA]           -> ISOSEQ_REFINE -> ISOSEQ_CLUSTER2
      //   cluster  [03_ISOSEQ_REFINE]  -> ISOSEQ_CLUSTER2
      // INFO: cluster_engine cdhit/rattle swaps ISOSEQ_CLUSTER2 + BAM_TO_FA for SAMTOOLS_FASTA -> CDHIT_EST/RATTLE;
      // INFO: their cluster entrypoint also takes tag-less FASTA/FASTQ (e.g. SRA reads)
      ch_versions = Channel.empty()

      // WARN: only LIMA and ISOSEQ_REFINE read primers; cluster runs neither, so the file is not required or opened
      ch_primers = entrypoint == 'cluster'
          ? Channel.empty()
          : Channel.value(file(global_primers, checkIfExists: true))

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          CHANNELING/INDEXING [ SUBREADS, CCS, LIMA OUTPUTS, REFINE OUTPUTS ]
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: cluster also reads FASTA/FASTQ: [ [id: file name w/o .fast[aq](.gz), single_end], fastx ]
      def input_glob = entrypoint == 'cluster' ? '*.{bam,fasta,fasta.gz,fastq,fastq.gz}' : '*.bam'
      Channel
          .fromPath("${global_input_dir}/${input_glob}", checkIfExists: true)
          .branch { f ->
              bam:   f.name.endsWith('.bam')
              fastx: true
          }
          .set { ch_input }

      ch_fastx = ch_input.fastx.map { fastx ->
          if (cluster_engine == 'isoseq') {
              error "ERROR: ${fastx.name} carries no PacBio tags for isoseq cluster2; use --cluster_engine cdhit or rattle"
          }
          [ [ id: fastx.name.replaceFirst(/\.fast[aq](\.gz)?$/, ''), single_end: true ], fastx ]
      }

      // INFO: every *.bam in global_input_dir -> [ [id: file name w/o .bam, single_end, indexed], bam, pbi or [] ]
      ch_input.bam
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


      // INFO: PacBio tools need a .pbi next to each BAM; build it only where it is missing
      // INFO: (cluster BAMs bound for cdhit/rattle only meet samtools, so they need none)
      def needs_pbi = !(entrypoint == 'cluster' && cluster_engine != 'isoseq')
      ch_bam
          .branch {
              indexed:     it[0].indexed || !needs_pbi
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

      // INFO: ch_bam is now [ [id, single_end, indexed: true], bam, pbi ] for every input BAM
      ch_bam = ch_bam_branched.indexed.mix(ch_bam_reindexed)
      ch_versions = ch_versions.mix(PBINDEX.out.versions)

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ENTRYPOINT BRANCHING  [ SUBREADS, CCS, REFINE, CLUSTER ]
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: ch_ccs_bams = [ meta, bam, pbi ] CCS reads for Skera/LIMA; stays empty for refine and cluster
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
          ch_ccs_bams = PBMERGE.out.bam.join(PBMERGE.out.pbi) // INFO: [ [id: movie, single_end], bam, pbi ]

          ch_versions = ch_versions.mix(PBMERGE.out.versions)
          ch_versions = ch_versions.mix(PBCCS.out.versions)

        break

        case 'ccs':
          ch_ccs_bams = ch_bam
        break

        case 'refine':
          // LIMA BAMs are already demultiplexed and primer-trimmed.
        break

        case 'cluster':
          // INFO: refined (FLNC) BAMs are picked up at ISOSEQ_REFINE below and go straight to clustering
        break

        default:
          error """
          ERROR: Unknown entrypoint -> options at this step are: subreads, ccs, refine, cluster
          """.stripIndent()
          System.exit(1)
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          SKERA [ DEMULTIPLEX KINNEX ARRAY ]
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: Kinnex arrays hold several transcripts per read; Skera splits them before LIMA (meta unchanged)
      ch_skera_demux_bams = Channel.empty()
      if (!(entrypoint in ['refine', 'cluster']) && is_kinnex_library) {
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

      // INFO: LIMA demultiplexes by primer pair and writes one <id>_fl.<5p>--<3p>.bam per pair found
      // INFO: refine reads those BAMs back from 02_LIMA; cluster skips LIMA, so the channel stays empty
      ch_lima_out_bams = Channel.empty()
      if (entrypoint == 'refine') {
          ch_lima_out_bams = ch_bam
      } else if (entrypoint != 'cluster') {
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
      // INFO: meta becomes [ ..., id: <movie>.<5p>--<3p>, parent_id: <movie>, barcode: <5p>--<3p> ];
      // INFO: this id is the sample id from here on (03_ISOSEQ_REFINE files are named <id>_flnc.bam)
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

      // INFO: ch_refined_bams = [ meta, <id>_flnc.bam ], the single input of every clustering branch below
      ch_refined_bams = Channel.empty()
      if (entrypoint == 'cluster') {
          // INFO: restart from 03_ISOSEQ_REFINE; strip _flnc.bam to recover the sample id the full run used
          ch_refined_bams = ch_bam.map { meta, bam, _pbi ->
              def matcher = bam.name =~ /^(.+)_flnc\.bam$/
              if (!matcher.matches()) error "Unexpected refined BAM name: ${bam.name}"
              [ [ id: matcher.group(1), single_end: true ], bam ]
          }
      } else {
          ISOSEQ_REFINE(ch_lima_out_bams, ch_primers) // INFO: discard CCS without polyA tails
          ch_refined_bams = ISOSEQ_REFINE.out.bam
          ch_versions = ch_versions.mix(ISOSEQ_REFINE.out.versions)
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ISOSEQ_CLUSTER2 [ CLUSTER READS ]
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // Refine each sample once; only clustering needs separate and pooled branches.
      // INFO: per_sample clusters each sample alone; multi_sample pools all samples into one set named global_prefix;
      // INFO: both runs the two from the same refined reads. BAM_TO_FA splits each result into hq and singletons
      // INFO: and adds meta.singleton (false / true); meta.sample_id is set later in SPLIT_ALIGN (falls back to id).
      ch_pbccs_merged_flnc_clustered_fa = Channel.empty()
      if (cluster_engine == 'isoseq' && cluster_mode in ['per_sample', 'both']) {
        ISOSEQ_CLUSTER2(ch_refined_bams) // INFO: cluster reads, one task per sample id
        BAM_TO_FA(ISOSEQ_CLUSTER2.out.bam)

        ch_pbccs_merged_flnc_clustered_fa  = ch_pbccs_merged_flnc_clustered_fa.mix(BAM_TO_FA.out.singletons)
        ch_pbccs_merged_flnc_clustered_fa  = ch_pbccs_merged_flnc_clustered_fa.mix(BAM_TO_FA.out.hq)

        ch_versions = ch_versions.mix(ISOSEQ_CLUSTER2.out.versions)
        ch_versions = ch_versions.mix(BAM_TO_FA.out.versions)
      }

      if (cluster_engine == 'isoseq' && cluster_mode in ['multi_sample', 'both']) {
        // cluster2 accepts a FOFN, so pooling does not need an intermediate BAM merge.
        // INFO: per-sample metas are dropped; the pooled item is [ [id: global_prefix, single_end], [bams] ]
        ch_refined_bams
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

      if (cluster_engine != 'isoseq') {
        // INFO: neither engine reads BAM: refined BAMs become FASTA once (both modes reuse it);
        // INFO: cluster-entrypoint FASTA/FASTQ go in as they are. Items: [ meta, reads ]
        SAMTOOLS_FASTA(ch_refined_bams)
        ch_engine_reads = SAMTOOLS_FASTA.out.fasta.mix(ch_fastx)
        ch_versions = ch_versions.mix(SAMTOOLS_FASTA.out.versions)

        ch_engine_in = Channel.empty()
        if (cluster_mode in ['per_sample', 'both']) {
          ch_engine_in = ch_engine_in.mix(ch_engine_reads)
        }
        if (cluster_mode in ['multi_sample', 'both']) {
          // INFO: one task takes every file of the pool (rattle reads them all, cd-hit merges in-task);
          // INFO: sorted by name since cd-hit is order-dependent. [ [id: global_prefix, single_end], [reads] ]
          ch_engine_in = ch_engine_in.mix(
            ch_engine_reads
              .map { meta, reads -> reads }
              .collect(sort: { a, b -> a.name <=> b.name })
              .map { reads -> [ [ id: prefix, single_end: true ], reads ] }
          )
        }

        // INFO: same contract as BAM_TO_FA: [ meta + [singleton], <id>.{hq,singletons}.fasta.gz ]
        if (cluster_engine == 'cdhit') {
          CDHIT_EST(ch_engine_in)
          ch_pbccs_merged_flnc_clustered_fa = CDHIT_EST.out.singletons.mix(CDHIT_EST.out.hq)
          ch_versions = ch_versions.mix(CDHIT_EST.out.versions)
        } else {
          RATTLE(ch_engine_in)
          ch_pbccs_merged_flnc_clustered_fa = RATTLE.out.singletons.mix(RATTLE.out.hq)
          ch_versions = ch_versions.mix(RATTLE.out.versions)
        }
      }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        OUTPUT CHANNELS
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    emit:
        reads   = ch_pbccs_merged_flnc_clustered_fa // [ meta + [singleton], <id>.{hq,singletons}.fasta.gz ]
        versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
