/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTX_INSPECT } from '../../modules/custom/isotools/fastx/inspect/main.nf'
include { FASTX_ORIENT } from '../../modules/custom/isotools/fastx/orient/main.nf'
include { LIMA as LIMA_FASTX } from '../../modules/nf-core/lima/main.nf'
include { FASTX_CONCAT as FASTX_CONCAT_LIMA } from '../../modules/custom/fastx/concat/main.nf'

include { POOL_READS } from '../pool_reads/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FASTX_PREPARE {
    take:
      input_dir    // path: directory with *.fast[aq](.gz) (flnc entrypoint)
      primers      // path or null: user primer FASTA (global_primers)
      input_state  // string [ auto, ccs, fl, flnc, mixed, clustered ]
      mode         // string [ per_sample, multi_sample, both ]
      prefix       // string (global_prefix)

    main:
      ch_versions = Channel.empty()

      // INFO: sample_id = name minus .hq/.singletons and .fast[aq](.gz), so X.hq + X.singletons are one sample X;
      // INFO: singleton = name contains "singleton". [ [id, sample_id, single_end, singleton], fastx ]
      Channel
          .fromPath("${input_dir}/*.fast*", checkIfExists: true)
          .map { fastx ->
              [
                  [
                      id:         fastx.baseName,
                      sample_id:  fastx.name.replaceFirst(/(?:\.(?:hq|singletons))?\.fast[aq](?:\.gz)?$/, ''),
                      single_end: true,
                      singleton:  fastx.baseName.contains("singleton")
                  ],
                  fastx
              ]
          }
          .set { ch_files }

      def user_primers = primers ? file(primers, checkIfExists: true) : []

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          STATE [ iso-fastx inspect ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: [ meta + [state, kit], fastx, row ] where row is the inspect TSV row (empty when forced)
      if (input_state == 'auto') {
          FASTX_INSPECT(ch_files, user_primers)
          ch_state = FASTX_INSPECT.out.inspected.map { meta, fastx, tsv ->
              def row = tsv.splitCsv(header: true, sep: '\t')[0]
              [ meta + [ state: row.state, kit: row.kit ], fastx, row ]
          }
          ch_versions = ch_versions.mix(FASTX_INSPECT.out.versions)
      } else {
          ch_state = ch_files.map { meta, fastx -> [ meta + [ state: input_state, kit: 'custom' ], fastx, [:] ] }
      }

      ch_state
          .branch { meta, fastx, row ->
              fail:  meta.state in ['subreads', 'ambiguous']
              empty: meta.state == 'empty'
              ccs:   meta.state == 'ccs'
              mixed: meta.state == 'mixed'
              ready: true // INFO: fl, flnc and clustered go in as they are
          }
          .set { ch_routed }

      ch_routed.empty.subscribe { meta, fastx, row -> log.warn "FASTX: ${fastx.name} holds no reads, skipped" }

      // INFO: never emits; mixed into the output below so the error is always evaluated
      ch_failed = ch_routed.fail.map { meta, fastx, row ->
          def hint = meta.state == 'subreads'
              ? 'raw subreads cannot become HiFi from FASTQ: fetch the original BAM (ENA submitted_ftp or NCBI Cloud Data Delivery) and use --entrypoint subreads'
              : 'set --flnc_input_state if you know the state'
          error "ERROR: ${fastx.name} looks like '${meta.state}' (${row.findAll { k, v -> k != 'file' }}); ${hint}"
      }

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          NORMALIZATION [ lima on raw CCS, orientation of mixed files ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: user primers win; otherwise the kit iso-fastx detected picks a bundled primer set
      ch_routed.ccs
          .map { meta, fastx, row ->
              [ meta, fastx, user_primers ?: file("${moduleDir}/../../../assets/primers/${meta.kit}.fasta", checkIfExists: true) ]
          }
          .multiMap { meta, fastx, kit_primers ->
              reads:   [ meta, fastx, [] ]
              primers: kit_primers
          }
          .set { ch_lima_in }

      LIMA_FASTX(ch_lima_in.reads, ch_lima_in.primers)
      ch_versions = ch_versions.mix(LIMA_FASTX.out.versions)

      // INFO: lima --isoseq writes one file per primer pair; join them back into one file per input
      LIMA_FASTX.out.fastqgz
          .mix(LIMA_FASTX.out.fastq, LIMA_FASTX.out.fastagz, LIMA_FASTX.out.fasta)
          .map { meta, files ->
              def reads = files instanceof List ? files : [files]
              def ext = (reads[0].name =~ /\.fast[aq](\.gz)?$/)[0][0]
              [ meta + [ outfile: "${meta.id}.fl${ext}" ], reads ]
          }
          .set { ch_lima_out }

      FASTX_CONCAT_LIMA(ch_lima_out)
      ch_versions = ch_versions.mix(FASTX_CONCAT_LIMA.out.versions)

      FASTX_ORIENT(ch_routed.mixed.map { meta, fastx, row -> [ meta, fastx ] })
      ch_versions = ch_versions.mix(FASTX_ORIENT.out.versions)

      // INFO: clustered records stand for many reads each: collapse chain then keeps every distinct chain
      ch_routed.ready
          .map { meta, fastx, row -> [ meta, fastx ] }
          .mix(
              FASTX_CONCAT_LIMA.out.reads.map { meta, fastx -> [ meta.findAll { k, v -> k != 'outfile' }, fastx ] },
              FASTX_ORIENT.out.reads,
              ch_failed
          )
          .map { meta, fastx -> [ meta + [ clustered: meta.state == 'clustered' ], fastx ] }
          .set { ch_prepared }

      POOL_READS(ch_prepared, mode, prefix)
      ch_versions = ch_versions.mix(POOL_READS.out.versions)

    emit:
      reads    = POOL_READS.out.reads // [ meta + [state, kit, clustered], fastx ], per sample and/or pooled
      versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
