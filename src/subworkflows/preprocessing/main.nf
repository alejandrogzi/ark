/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ISOSEQ } from '../isoseq/main.nf'
include { POOL_READS } from '../pool_reads/main.nf'
include { GENOME } from '../genome/main.nf'

include { GXF2BED } from '../../modules/custom/gxf2bed/main.nf'
include { BED2GTF } from '../../modules/custom/bed2gtf/main.nf'

include { MINIMAP2_INDEX } from '../../modules/nf-core/minimap2/index/main.nf'
include { MINIMAP2_INDEX as PBMM2_INDEX } from '../../modules/nf-core/minimap2/index/main.nf'
include { MINIMAP2_INDEX as ARK_INDEX } from '../../modules/nf-core/minimap2/index/main.nf'
include { MINIMAP2_INDEX as FLAIR_INDEX } from '../../modules/nf-core/minimap2/index/main.nf'
include { ULTRA_INDEX } from '../../modules/nf-core/ultra/index/main.nf'
include { DESALT_INDEX } from '../../modules/custom/desalt/index/main.nf'

include { GENEPRED_LINT } from '../../modules/custom/genepred/lint/main.nf'

include { SPLICING as MINISPLICE_GENOMIC_SPLICE_SCORES } from '../splicing/main.nf'
include { SPLICING as SPLICEAI_GENOMIC_SPLICE_SCORES } from '../splicing/main.nf'

include { WGET as WGET_PROTEIN_DATABASE } from '../../modules/nf-core/wget/main.nf'
include { UNTAR as UNTAR_PROTEIN_DATABASE } from '../../modules/nf-core/untar/main.nf'
include { GUNZIP as GUNZIP_DATABASE } from '../../modules/custom/gunzip/main.nf'
include { FASTA_MERGE } from '../../modules/custom/diamond/merge/main.nf'
include { DIAMOND_MAKEDB } from '../../modules/custom/diamond/makedb/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPROCESSING {
    take:
      entrypoint             // [ subreads, ccs, refine, cluster, flnc ]
      global_input_dir       // path
      global_primers         // path
      genome                 // path
      annotation             // path
      ccs_chunk              // int
      isoseq_cluster2_mode   // string
      protein_database       // path
      custom_database        // path
      raw_database           // path
      minimap2_index         // path
      use_splice_scores      // bool
      splice_score_algorithm // string
      bigwigs                // path
      minisplice             // path
      spliceai               // path
      compression            // bool
      global_prefix          // string
      aligner                // string [ ark, mm2, ultra, desalt, pbmm2, flair ]
      ultra_use_annotation   // bool
      ultra_index            // path
      desalt_index           // path
      pbmm2_index            // path
      is_kinnex_library      // bool
      ch_versions            // [ meta, versions.yml ]

    main:
      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          GENOME PROCESSING
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ch_genome = GENOME(genome)

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          PROTEIN DATABASE PROCESSING
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

    ch_database = Channel.empty()
    if (custom_database) {
      // INFO: .dmnd / .dmnd.gz -> replaces the default database entirely
      if (custom_database.endsWith('.dmnd')) {
          ch_database = Channel.fromPath(custom_database, checkIfExists: true)
      } else if (custom_database.endsWith('.dmnd.gz')) {
          GUNZIP_DATABASE(
              Channel.value(
                  [ [id: custom_database.tokenize('/')[-1]], custom_database ]
              )
          )
          GUNZIP_DATABASE.out.gunzip
            .map { meta, it -> it }
            .set { ch_database }
      // INFO: .fa / .fasta / .fa.gz / .fasta.gz -> downloaded raw database + custom
      // sequences are merged, then reindexed with DIAMOND_MAKEDB
      } else if (custom_database.endsWith('.fa') || custom_database.endsWith('.fasta')
                  || custom_database.endsWith('.fa.gz') || custom_database.endsWith('.fasta.gz')) {
          WGET_PROTEIN_DATABASE(
              Channel.value(raw_database)
              .map { it -> [ [ id: 'uniprot_sprot.fasta.gz' ], it ] }
          )
          FASTA_MERGE(
              Channel.fromPath(custom_database, checkIfExists: true)
                .map { it -> [ [ id: it.baseName ], it ] }
                .combine(WGET_PROTEIN_DATABASE.out.outfile.map { meta, raw -> raw })
          )
          DIAMOND_MAKEDB(
              FASTA_MERGE.out.fasta.map { meta, fasta -> [ [ id: 'merged_database' ], fasta ] },
              [],
              [],
              []
          )
          DIAMOND_MAKEDB.out.db
            .map { meta, it -> it }
            .set { ch_database }
          ch_versions = ch_versions.mix(FASTA_MERGE.out.versions).mix(DIAMOND_MAKEDB.out.versions)
      } else {
          error """
          ERROR: custom_database extension not recognized.
          Please provide a custom database in one of these formats:
            .dmnd/.dmnd.gz   -> replaces the default database entirely
            .fa/.fasta/.fa.gz/.fasta.gz -> appended to the default SwissProt database and reindexed
          """.stripIndent()
      }
    } else {
      WGET_PROTEIN_DATABASE(
        Channel.value(protein_database)
        .map { it -> [ [ id: 'uniprot_sprot.tar.gz' ], it ] }
      )
      UNTAR_PROTEIN_DATABASE(WGET_PROTEIN_DATABASE.out.outfile)
      ch_database = UNTAR_PROTEIN_DATABASE.out.contents.map { meta, it -> it }
    }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ANNOTATION PROCESSING -> GTF/BED CHANNELS
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ch_reference_transcripts = Channel.value(file(annotation, checkIfExists: true))
      ch_reference_transcripts_gtf = Channel.empty()
      GENEPRED_LINT(ch_reference_transcripts.map { file -> [ [ id:file.baseName ], file ] })

      if (annotation.endsWith(".gtf.gz")
        | annotation.endsWith(".gtf") 
        | annotation.endsWith(".gff.gz") 
        | annotation.endsWith(".gff")) {
          ch_reference_transcripts_gtf = ch_reference_transcripts
            .map { gtf -> [ [id:gtf.baseName], gtf ] }

          GXF2BED(ch_reference_transcripts_gtf)
          ch_reference_transcripts = GXF2BED.out.bed
      } else {
          BED2GTF(
            ch_reference_transcripts.map { bed -> [ [id:bed.baseName], bed ] },
            Channel.of([[:], []])
          )

          ch_reference_transcripts_gtf = BED2GTF.out.gtf
          ch_reference_transcripts = ch_reference_transcripts.map { bed -> [ [id:bed.baseName], bed ] }
      }


      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          INDEX BUILDING [ mm2, ultra, pbmm2, desalt, ark, flair ]
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ch_genome_index = Channel.empty()

      switch (aligner) {
        case 'ark':
          if (minimap2_index) {
              ch_genome_index = Channel.value(file(minimap2_index, checkIfExists: true))
                  .map { path -> [ [id:path.name], path ] }
          } else {
              ARK_INDEX(
                  ch_genome.genome.map { genome -> [ [id:genome.baseName], genome ] }
              )
              ch_genome_index = ARK_INDEX.out.index
              ch_versions = ch_versions.mix(ARK_INDEX.out.versions)
          }
        break
        
        case 'mm2':
          if (minimap2_index) {
              ch_genome_index = Channel.value(file(minimap2_index, checkIfExists: true))
                  .map { path -> [ [id:path.name], path ] }
          } else {
              MINIMAP2_INDEX(
                  ch_genome.genome.map { genome -> [ [id:genome.baseName], genome ] }
              )
              ch_genome_index = MINIMAP2_INDEX.out.index
              ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
          }
        break

        case 'ultra':
          if (ultra_index) {
              ch_genome_index = Channel.value(file(ultra_index, checkIfExists: true))
                  .map { path -> [ [id:path.name], path ] }
          } else {
              if (ultra_use_annotation) {
                ULTRA_INDEX(
                    ch_genome.genome.map { genome -> [ [id:genome.baseName], genome ] },
                    ch_reference_transcripts_gtf
                )
              } else {
                ULTRA_INDEX(
                    ch_genome.genome.map { genome -> [ [id:genome.baseName], genome ] },
                    Channel.of([[:], []])
                )
              }

              ch_genome_index = ULTRA_INDEX.out.index
              ch_versions = ch_versions.mix(ULTRA_INDEX.out.versions)
          }
        break

        case 'desalt':
          if (desalt_index) {
              ch_genome_index = Channel.value(file(desalt_index, checkIfExists: true))
                  .map { path -> [ [id:path.name], path ] }
          } else {
              DESALT_INDEX(
                  ch_genome.genome.map { genome -> [ [id:genome.baseName], genome ] }
              )

              ch_genome_index = DESALT_INDEX.out.index
              ch_versions = ch_versions.mix(DESALT_INDEX.out.versions)
          }
        break

        case 'pbmm2':
          if (pbmm2_index) {
              ch_genome_index = Channel.value(file(pbmm2_index, checkIfExists: true))
                  .map { path -> [ [id:path.name], path ] }
          } else {
              PBMM2_INDEX(
                  ch_genome.genome.map { genome -> [ [id:genome.baseName], genome ] }
              )

              ch_genome_index = PBMM2_INDEX.out.index
              ch_versions = ch_versions.mix(PBMM2_INDEX.out.versions)
          } 
        break

        case 'flair':
          if (minimap2_index) {
              ch_genome_index = Channel.value(file(minimap2_index, checkIfExists: true))
                  .map { path -> [ [id:path.name], path ] }
          } else {
              FLAIR_INDEX(
                  ch_genome.genome.map { genome -> [ [id:genome.baseName], genome ] }
              )
              ch_genome_index = FLAIR_INDEX.out.index
              ch_versions = ch_versions.mix(FLAIR_INDEX.out.versions)
          }
        break

        default:
          error """
          ERROR: Unknown aligner: '${aligner}'.
          Valid aligners are: ark, mm2, ultra, desalt, pbmm2, flair
          """.stripIndent()
          System.exit(1)
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          GENOMIC SPLICE SCORE BUILDERS [ ark only ]
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ch_splice_scores = Channel.empty()
      ch_spliceai_bigwigs = Channel.value([[:], []])
      if (use_splice_scores && aligner == "ark") {
          if (splice_score_algorithm == "spliceai") {
              SPLICEAI_GENOMIC_SPLICE_SCORES(
                  ch_genome.genome,
                  ch_reference_transcripts,
                  bigwigs,
                  splice_score_algorithm,
                  minisplice,
                  spliceai,
                  ch_genome.chrom_sizes,
                  compression,
                  ch_versions
              )
              ch_splice_scores = SPLICEAI_GENOMIC_SPLICE_SCORES.out.scores
              ch_spliceai_bigwigs = SPLICEAI_GENOMIC_SPLICE_SCORES.out.bigwigs
          } else if (splice_score_algorithm == "minisplice") {
              MINISPLICE_GENOMIC_SPLICE_SCORES(
                  ch_genome.genome,
                  ch_reference_transcripts,
                  bigwigs,
                  splice_score_algorithm,
                  minisplice,
                  spliceai,
                  ch_genome.chrom_sizes,
                  compression,
                  ch_versions
              )             
              ch_splice_scores = MINISPLICE_GENOMIC_SPLICE_SCORES.out.scores
          }
      } else {
          // INFO: keep ARK_ALIGN schedulable when splice scores are disabled
          ch_splice_scores = Channel.value([[:], []])
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ENTRYPOINTS [ subreads, ccs, refine, cluster, flnc ]
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: every entrypoint ends here with ch_reads = [ meta, reads FASTA/FASTQ ] for SPLIT_ALIGN
      // INFO: Iso-Seq entrypoints: BAMs -> ISOSEQ (CCS/LIMA/refine/cluster2 from the stage the entrypoint names)
      ch_reads = Channel.empty()
      if (entrypoint in ['subreads', 'ccs', 'refine', 'cluster']) {
          ISOSEQ(
              global_input_dir,
              global_primers,
              ccs_chunk,
              isoseq_cluster2_mode,
              global_prefix,
              entrypoint,
              is_kinnex_library
          )

          ch_reads = ch_reads.mix(ISOSEQ.out.reads)
          ch_versions = ch_versions.mix(ISOSEQ.out.versions)
      } else if (entrypoint == "flnc") {
          // INFO: flnc: FASTA/FASTQ files from global_input_dir, no Iso-Seq steps. meta.id = file baseName;
          // INFO: sample_id = name minus .hq/.singletons and .fast[aq](.gz), so X.hq + X.singletons are one sample X
          // INFO: singleton = name contains "singleton"
          // INFO: pooling mirrors the Iso-Seq entrypoints (pool at cluster2): multi_sample/both concatenate
          // INFO: every input file per hq/singleton class into one pooled sample named global_prefix, so
          // INFO: adapters/segmentation/collapse see all samples together per chromosome downstream
          Channel
              .fromPath("${global_input_dir}/*.fast*", checkIfExists: true)
              .map { fastx ->
                  def singleton = fastx.baseName.contains("singleton")
                  return [
                      [
                          id:         fastx.baseName,
                          sample_id:  fastx.name.replaceFirst(/(?:\.(?:hq|singletons))?\.fast[aq](?:\.gz)?$/, ''),
                          single_end: true,
                          singleton:  singleton
                      ],
                      fastx,
                  ]
              }
              .set { ch_flnc_reads }

          POOL_READS(
              ch_flnc_reads,
              isoseq_cluster2_mode,
              global_prefix
          )
          ch_reads = ch_reads.mix(POOL_READS.out.reads)
          ch_versions = ch_versions.mix(POOL_READS.out.versions)
      }
      
      ch_versions = ch_versions.mix(ch_genome.versions)

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          OUTPUT CHANNELS
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

    emit:
      genome                    = ch_genome.genome  // Channel.value(path)
      database                  = ch_database       // Channel.value(path)
      reads                     = ch_reads          // channel: [ val(meta), [ reads ] ]
      genome_index              = ch_genome_index   // channel: [ val(meta), [ index ] ]
      reference_transcripts     = ch_reference_transcripts     // channel: [ val(meta), [ annotation ] ]
      reference_transcripts_gtf = ch_reference_transcripts_gtf // channel: [ val(meta), [ annotation ] ]
      splice_scores             = ch_splice_scores              // channel: [ val(meta), [ scores ] ]
      bigwigs                   = ch_spliceai_bigwigs           // channel: [ val(meta), [ bigwigs ] ]
      chrom_sizes               = ch_genome.chrom_sizes         // channel: [ val(meta), [ chrom_sizes ] ]
      versions                  = ch_versions                   // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
