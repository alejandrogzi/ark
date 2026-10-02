/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MINISPLICE_DOWNLOAD } from '../../modules/custom/minisplice/download/main.nf'
include { MINISPLICE_PREDICT } from '../../modules/custom/minisplice/predict/main.nf'
include { SPLICEAI_DERIVE } from '../../modules/custom/spliceai/derive/main.nf'
include { GUNZIP as GUNZIP_SPLICEAI } from '../../modules/custom/gunzip/main.nf'
include { GUNZIP as GUNZIP_MINISPLICE } from '../../modules/custom/gunzip/main.nf'
include { SPLICEAI as SPLICEAI_RUN } from '../spliceai/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow SPLICING {
    take:
      genome         // channel: [ val(meta), [ genome ] ]
      annotation     // channel: [ val(meta), [ annotation ] ]
      bigwigs        // path
      algorithm      // string
      minisplice     // path
      spliceai       // path
      chromsizes     // channel: [  chromsizes ]
      compression    // bool
      ch_versions    // channel: [ path(version) ]

    main:
      // INFO: emits genome-wide splice scores as [ meta, scores ], used by ARK_ALIGN
      // INFO: decision tree per algorithm:
      // INFO:   spliceai   -> precomputed scores > derive from given bigwigs > run SpliceAI + derive
      // INFO:   minisplice -> precomputed scores > download model + predict
      ch_spliceai_bigwigs = Channel.value([[:], []])

      if (algorithm == "spliceai") {
          /*
          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              SPLICEAI [ SCORES, BIGWIGS, MODEL ]
          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          */

          // INFO: user bigwig dir -> [ [ id: spliceai ], dir ], emitted as `bigwigs`
          // INFO: (downstream: iso-classify --bigwig and the POLISH orphan finder)
          if (bigwigs) {
            Channel.value([
                    [ id: "spliceai" ],
                    file(bigwigs, checkIfExists: true)
            ]).set { ch_spliceai_bigwigs }
          }

          // INFO: precomputed scores win over everything else (gunzipped if needed)
          if (spliceai) {
            // WARN: no `def` here: Nextflow 24.10 refuses `def x = f(<take input>)` in a workflow body
            spliceai_scores = file(spliceai, checkIfExists: true)

            if (spliceai_scores.toString().endsWith(".gz")) {
                GUNZIP_SPLICEAI([ [ id: "spliceai" ], spliceai_scores ])
                ch_scores = GUNZIP_SPLICEAI.out.gunzip
            } else {
                Channel.value([
                        [ id: "spliceai" ],
                        spliceai_scores
                ]).set { ch_scores }
            }
          } else {
            // INFO: if bigwig dir is given, derive scores from it (no model run)
            if (bigwigs) {
              SPLICEAI_DERIVE(
                  genome.map { genome -> [ [id:genome.baseName], genome ] },
                  annotation,
                  ch_spliceai_bigwigs
              )
              ch_scores = SPLICEAI_DERIVE.out.scores
              ch_versions = ch_versions.mix(SPLICEAI_DERIVE.out.versions)
            } else {
                // INFO: nothing given -> run SpliceAI genome-wide, then derive scores from its bigwigs
                SPLICEAI_RUN(
                    genome.map { genome -> [ [id:genome.baseName], genome ] },
                    chromsizes,
                    compression,
                    ch_versions
                )

                // INFO: SPLICEAI_RUN emits a bare dir; wrap it as [ meta, dir ] like a user-given
                // INFO: bigwig dir, so SPLICEAI_DERIVE and the `bigwigs` emit see the same shape
                ch_spliceai_bigwigs = SPLICEAI_RUN.out.spliceai.map { dir -> [ [ id: "spliceai" ], dir ] }

                SPLICEAI_DERIVE(
                    genome.map { genome -> [ [id:genome.baseName], genome ] },
                    annotation,
                    ch_spliceai_bigwigs
                )

                ch_scores = SPLICEAI_DERIVE.out.scores
                ch_versions = ch_versions.mix(SPLICEAI_RUN.out.versions)
                ch_versions = ch_versions.mix(SPLICEAI_DERIVE.out.versions)
            }
          }
      } else if (algorithm == "minisplice") {
          /*
          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              MINISPLICE [ SCORES, MODEL ]
          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          */

          // INFO: precomputed scores win (gunzipped if needed), else download model + predict
          if (minisplice) {
              minisplice_scores = file(minisplice, checkIfExists: true)

              if (minisplice_scores.toString().endsWith(".gz")) {
                  GUNZIP_MINISPLICE([ [ id: "minisplice" ], minisplice_scores ])
                  ch_scores = GUNZIP_MINISPLICE.out.gunzip
              } else {
                  Channel.value([
                          [ id: "minisplice" ],
                          minisplice_scores
                  ]).set { ch_scores }
              }
          } else {
            MINISPLICE_DOWNLOAD()
            MINISPLICE_PREDICT(
                genome.map { genome -> [ [id:genome.baseName], genome ] },
                MINISPLICE_DOWNLOAD.out.model,
                MINISPLICE_DOWNLOAD.out.calibration
            )

            ch_scores = MINISPLICE_PREDICT.out.scores
            ch_versions = ch_versions.mix(MINISPLICE_DOWNLOAD.out.versions)
            ch_versions = ch_versions.mix(MINISPLICE_PREDICT.out.versions)
          }
      }

    emit:
      scores = ch_scores
      bigwigs = ch_spliceai_bigwigs
      versions = ch_versions
}
