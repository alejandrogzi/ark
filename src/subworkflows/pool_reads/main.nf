/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTX_CONCAT } from '../../modules/custom/fastx/concat/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow POOL_READS {
    take:
      ch_reads  // [ meta, fastx ]; meta carries sample_id + singleton (flnc file channeling)
      mode      // string [ per_sample, multi_sample, both ]
      prefix    // string (global_prefix: sample_id of the pooled branch)

    main:
      ch_versions = Channel.empty()
      ch_out = Channel.empty()

      if (mode in ['per_sample', 'both']) {
          ch_out = ch_out.mix(ch_reads)
      }

      if (mode in ['multi_sample', 'both']) {
          if (!prefix) {
              error """
              ERROR: pooling flnc reads needs --global_prefix (empty with cluster_mode '${mode}')
              """.stripIndent()
          }

          // INFO: one pooled file per hq/singleton class; singleton flag survives so SEGMENT
          // INFO: still gets --singleton on the right reads. Sorted here for deterministic
          // INFO: meta/error output; FASTX_CONCAT re-sorts in-task (channel order is not a guarantee).
          ch_reads
              .map { meta, fastx -> [ meta.singleton, fastx ] }
              .groupTuple(by: 0)
              .map { singleton, files ->
                  def sorted = files.sort { a, b -> a.name <=> b.name }
                  def exts = sorted.collect { f ->
                      def m = (f.name =~ /(\.fast[aq])(\.gz)?$/)
                      m.find() ? m.group(0) : null
                  }
                  def bad = []
                  sorted.eachWithIndex { f, i -> if (!exts[i]) { bad << f.name } }
                  if (bad) {
                      error "ERROR: cannot pool flnc reads with unrecognized read-file suffix (need .fasta/.fastq[.gz]): ${bad.join(', ')}"
                  }
                  if (exts.unique().size() > 1) {
                      error "ERROR: cannot pool flnc reads with mixed formats into one file: ${sorted.collect { it.name }.join(', ')}"
                  }
                  // INFO: sample_id == prefix makes the pool flow as one sample downstream;
                  // INFO: id stays distinct from every per-sample id so 'both' never mixes branches.
                  def cls = singleton ? '.singletons' : ''
                  def id = "${prefix}.pooled${cls}"
                  def meta = [
                      id:         id,
                      sample_id:  prefix,
                      single_end: true,
                      singleton:  singleton,
                      outfile:    "${id}${exts[0]}"
                  ]
                  [ meta, sorted ]
              }
              .set { ch_pool_in }

          FASTX_CONCAT(ch_pool_in)
          FASTX_CONCAT.out.reads
              .map { meta, pooled -> [ meta.findAll { k, v -> k != 'outfile' }, pooled ] }
              .set { ch_pooled }
          ch_out = ch_out.mix(ch_pooled)
          ch_versions = ch_versions.mix(FASTX_CONCAT.out.versions)
      }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        OUTPUT CHANNELS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    emit:
      reads    = ch_out // per-sample and/or pooled [ meta, fastx ], gated by mode
      versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
