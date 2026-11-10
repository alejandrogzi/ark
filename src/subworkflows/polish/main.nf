/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ISOTOOLS_INTRON_RETENTION } from '../../modules/custom/isotools/intron/main.nf'
include { ISOTOOLS_PAS_CALLER } from '../../modules/custom/isotools/pas/main.nf'
include { ISOTOOLS_ORPHAN as ISOTOOLS_ORPHAN_FINDER } from '../../modules/custom/isotools/orphan/main.nf'
include { ISOTOOLS_TRUNCATION_DETECTOR } from '../../modules/custom/isotools/utr/main.nf'
include { VEREDICT as ISOTOOLS_PLUGIN_VEREDICT } from '../../modules/custom/veredict/main.nf'
include { RESCUE_FLNC_FROM_UNPREDICTED_ORFS } from '../../modules/custom/rescue/unpredicted/main.nf'
include { RESCUE_FLNC_FROM_ARTIFACTS } from '../../modules/custom/rescue/artifacts/main.nf'
include { RESCUE_FLNC_FROM_ARTIFACTS as RESCUE_FLNC_FROM_ORPHANS } from '../../modules/custom/rescue/artifacts/main.nf'

include { GAWK_JOIN as JOIN_VEREDICT_PASSES } from '../../modules/custom/gawk/join/main.nf'
include { GAWK_JOIN as JOIN_VEREDICT_TRASH } from '../../modules/custom/gawk/join/main.nf'
include { GAWK_JOIN as JOIN_VEREDICT_RETENTIONS } from '../../modules/custom/gawk/join/main.nf'
include { GAWK_JOIN as JOIN_VEREDICT_TRUNCATIONS } from '../../modules/custom/gawk/join/main.nf'
include { GAWK_JOIN as JOIN_VEREDICT_INTRAPRIMMING } from '../../modules/custom/gawk/join/main.nf'
include { GAWK_JOIN as JOIN_VEREDICT_RT } from '../../modules/custom/gawk/join/main.nf'
include { GAWK_JOIN as JOIN_VEREDICT_ARTIFACTS } from '../../modules/custom/gawk/join/main.nf'
include { GAWK_JOIN as JOIN_UNPREDICTED } from '../../modules/custom/gawk/join/main.nf'

include { BEDTOBIGBED as BEDTOBIGBED_PASSES } from '../../modules/custom/bigtools/bedtobigbed/main.nf'
include { BEDTOBIGBED as BEDTOBIGBED_TRASH } from '../../modules/custom/bigtools/bedtobigbed/main.nf'
include { BEDTOBIGBED as BEDTOBIGBED_RETENTIONS } from '../../modules/custom/bigtools/bedtobigbed/main.nf'
include { BEDTOBIGBED as BEDTOBIGBED_TRUNCATIONS } from '../../modules/custom/bigtools/bedtobigbed/main.nf'
include { BEDTOBIGBED as BEDTOBIGBED_INTRAPRIMMING } from '../../modules/custom/bigtools/bedtobigbed/main.nf'
include { BEDTOBIGBED as BEDTOBIGBED_DUPLICATES } from '../../modules/custom/bigtools/bedtobigbed/main.nf'
include { BEDTOBIGBED as BEDTOBIGBED_SCRAPS } from '../../modules/custom/bigtools/bedtobigbed/main.nf'
include { BEDTOBIGBED as BEDTOBIGBED_RT } from '../../modules/custom/bigtools/bedtobigbed/main.nf'
include { BEDTOBIGBED as BEDTOBIGBED_ARTIFACTS } from '../../modules/custom/bigtools/bedtobigbed/main.nf'
include { BEDTOBIGBED as BEDTOBIGBED_UNPREDICTED } from '../../modules/custom/bigtools/bedtobigbed/main.nf'

include { PUBLISH as PUBLISH_BIGBEDS } from '../../modules/custom/publish/main.nf'

include { DETACH_DUPLICATES as DETACH_PASS_DUPLICATES } from '../../modules/custom/detach/main.nf'
include { DETACH_DUPLICATES as DETACH_RETENTION_DUPLICATES } from '../../modules/custom/detach/main.nf'
include { DETACH_DUPLICATES as DETACH_TRUNCATION_DUPLICATES } from '../../modules/custom/detach/main.nf'
include { DETACH_DUPLICATES as DETACH_INTRAPRIMMING_DUPLICATES } from '../../modules/custom/detach/main.nf'
include { DETACH_DUPLICATES as DETACH_RT_DUPLICATES } from '../../modules/custom/detach/main.nf'
include { DETACH_DUPLICATES as DETACH_ARTIFACT_DUPLICATES } from '../../modules/custom/detach/main.nf'
include { DETACH_DUPLICATES as DETACH_ORPHANS_DUPLICATES } from '../../modules/custom/detach/main.nf'
include { DETACH_DUPLICATES as DETACH_TRASH_DUPLICATES } from '../../modules/custom/detach/main.nf'

include { COLLAPSE as COLLAPSE_PASSES } from '../../modules/custom/collapse/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow POLISH {
    take:
      reads                  // channel: [ val(meta), [ reads ], [ introns ], [ orfs ] ]
      unpredicted            // channel: [ val(meta), models, xorf bed | [], nmd bed | [] ] every sample@chr that went into ORF calling
      annotation             // channel: [ val(meta), [ annotation ] ]
      forward_peaks          // channel: [ val(meta), [ forward_peaks ] ]
      reverse_peaks          // channel: [ val(meta), [ reverse_peaks ] ]
      chrom_sizes            // path
      splice_scores          // channel: [ val(meta), [ splice_scores ] ]
      autosql                // channel: [ autosql ]
      rescue_artifacts       // boolean: run RESCUE_FLNC_FROM_ARTIFACTS on the VEREDICT buckets
      ch_versions            // [ meta, versions.yml ]

    main:
      // INFO: reads -> [ meta, reads.bed, introns.tsv, orfs.tsv ], one item per sample@chromosome
      // INFO: meta.id = '<sample>@<chr>', meta.name = '<sample>'
      ch_reference_transcripts = annotation

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          DESCRIPTORS [ RETENTION, INTRAPRIMING, TRUNCATION ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      ISOTOOLS_INTRON_RETENTION(
        reads
      ) // INFO: intron retention per read, uses the classified introns

      ISOTOOLS_PAS_CALLER(
        reads,
        ch_reference_transcripts,
        forward_peaks,
        reverse_peaks
      ) // INFO: polyA sites vs. annotation + APARENT peaks -> intrapriming descriptor

      ISOTOOLS_TRUNCATION_DETECTOR(
        reads
      ) // INFO: truncation descriptor (iso-utr, reads against themselves)

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          VEREDICT [ PASS, TRASH, RETENTIONS, TRUNCATIONS, INTRAPRIMING, RT, ARTIFACTS ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: join by meta -> [ meta, reads, introns, orfs, retentions, intrapriming, truncations ]
      // INFO: VEREDICT reads the four descriptors (orfs + the three above); introns are ignored here
      reads
        .join(ISOTOOLS_INTRON_RETENTION.out.descriptor)
        .join(ISOTOOLS_PAS_CALLER.out.descriptor)
        .join(ISOTOOLS_TRUNCATION_DETECTOR.out.descriptor)
        .set { ch_polished_reads }

      ISOTOOLS_PLUGIN_VEREDICT(
        ch_polished_reads
      )

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          UNPREDICTED + RESCUE [ RESCUE_FLNC_FROM_UNPREDICTED_ORFS, RESCUE_FLNC_FROM_ARTIFACTS ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: unpredicted = models that went into ORF calling minus models with a row in any
      // INFO: output (VEREDICT buckets, NMD). Every sample@chr that went in is kept, even one with
      // INFO: no VEREDICT run (xORF kept nothing on it). Payloads are one element each so a
      // INFO: remainder join keeps the positions. Keys as strings: xORF ids are GStrings
      ISOTOOLS_PLUGIN_VEREDICT.out.buckets
          .map { meta, beds -> [ meta.id.toString(), beds instanceof List ? beds : [ beds ] ] }
          .set { ch_veredict_buckets }

      unpredicted
          .map { meta, models, predicted, nmd -> [ meta.id.toString(), [ meta, models, predicted, nmd ] ] }
          .join(ch_veredict_buckets, remainder: true)
          .filter { id, inputs, beds -> inputs != null }
          .map { id, inputs, beds ->
              def (meta, models, predicted, nmd) = inputs
              [ meta, models, predicted ?: [], (beds ?: []) + (nmd ? [ nmd ] : []) ]
          }
          .set { ch_unpredicted_inputs }

      // INFO: rescued models take a reference CDS (#RSU) and join pass below
      RESCUE_FLNC_FROM_UNPREDICTED_ORFS(
        ch_unpredicted_inputs,
        ch_reference_transcripts
      )

      // INFO: NMD calls per sample@chr, as they left xORF
      unpredicted
          .filter { meta, models, predicted, nmd -> nmd }
          .map { meta, models, predicted, nmd -> [ meta, nmd ] }
          .set { ch_nmd_calls }

      // INFO: buckets + NMD -> rt/artifact rows on annotated introns re-bucketed, reference-CDS rows
      // INFO: to pass; every input file comes back without the rows that moved
      if (rescue_artifacts) {
          // INFO: one task per sample@chr with VEREDICT buckets and/or NMD calls: a chromosome whose
          // INFO: ORF-called models all went to NMD has no VEREDICT run but still has NMD rows
          ISOTOOLS_PLUGIN_VEREDICT.out.buckets
              .map { meta, beds -> [ meta.id.toString(), meta, beds ] }
              .join(
                reads.map { meta, read, introns, orfs -> [ meta.id.toString(), introns ] }
              )
              .map { id, meta, beds, introns -> [ id, [ meta, beds, introns ] ] }
              .join(
                ch_nmd_calls.map { meta, nmd -> [ meta.id.toString(), [ meta, nmd ] ] },
                remainder: true
              )
              .filter { id, veredict, nmd -> veredict != null || nmd != null }
              .map { id, veredict, nmd ->
                  def meta = veredict ? veredict[0] : nmd[0]
                  [ meta, veredict ? veredict[1] : [], veredict ? veredict[2] : [], nmd ? nmd[1] : [] ]
              }
              .set { ch_rescue_artifacts }

          RESCUE_FLNC_FROM_ARTIFACTS(
            ch_rescue_artifacts,
            ch_reference_transcripts
          )

          ch_bucket_pass          = RESCUE_FLNC_FROM_ARTIFACTS.out.pass
          ch_bucket_trash         = RESCUE_FLNC_FROM_ARTIFACTS.out.trash
          ch_bucket_retentions    = RESCUE_FLNC_FROM_ARTIFACTS.out.retentions
          ch_bucket_truncations   = RESCUE_FLNC_FROM_ARTIFACTS.out.truncations
          ch_bucket_intraprimming = RESCUE_FLNC_FROM_ARTIFACTS.out.intraprimming
          ch_bucket_rt            = RESCUE_FLNC_FROM_ARTIFACTS.out.rt
          ch_bucket_artifacts     = RESCUE_FLNC_FROM_ARTIFACTS.out.artifacts
          ch_bucket_nmd           = RESCUE_FLNC_FROM_ARTIFACTS.out.nmd
          ch_rescue_audit         = RESCUE_FLNC_FROM_UNPREDICTED_ORFS.out.audit.mix(RESCUE_FLNC_FROM_ARTIFACTS.out.audit)
          ch_versions = ch_versions.mix(RESCUE_FLNC_FROM_ARTIFACTS.out.versions)
      } else {
          ch_bucket_pass          = ISOTOOLS_PLUGIN_VEREDICT.out.pass
          ch_bucket_trash         = ISOTOOLS_PLUGIN_VEREDICT.out.trash
          ch_bucket_retentions    = ISOTOOLS_PLUGIN_VEREDICT.out.retentions
          ch_bucket_truncations   = ISOTOOLS_PLUGIN_VEREDICT.out.truncations
          ch_bucket_intraprimming = ISOTOOLS_PLUGIN_VEREDICT.out.intraprimming
          ch_bucket_rt            = ISOTOOLS_PLUGIN_VEREDICT.out.rt
          ch_bucket_artifacts     = ISOTOOLS_PLUGIN_VEREDICT.out.artifacts
          ch_bucket_nmd           = ch_nmd_calls
          ch_rescue_audit         = RESCUE_FLNC_FROM_UNPREDICTED_ORFS.out.audit
      }
      ch_bucket_pass = ch_bucket_pass.mix(RESCUE_FLNC_FROM_UNPREDICTED_ORFS.out.rescued)

      // INFO: unpredicted models (BED12, no CDS, #UPNO / #UPVD) -> one BED + bigBed per sample
      RESCUE_FLNC_FROM_UNPREDICTED_ORFS.out.unpredicted
          .map { meta, bed -> [ meta.name, meta, bed ] }
          .groupTuple()
          .map { name, metas, files ->
              [ [ id: name + '.unpredicted', name: name ], files ]
          }
          .set { ch_unpredicted }
      JOIN_UNPREDICTED(ch_unpredicted, 'bed')
      BEDTOBIGBED_UNPREDICTED(JOIN_UNPREDICTED.out.output, chrom_sizes, [])

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          PER-SAMPLE TRACKS [ JOIN -> DETACH -> BEDTOBIGBED ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: VEREDICT runs per sample@chromosome; each class is regrouped by meta.name
      // INFO: so JOIN concatenates + sorts all chromosomes into one BED per sample,
      // INFO: DETACH splits '#DU'-tagged duplicates out and BEDTOBIGBED builds one track
      // INFO: meta from here on -> [ id: '<sample>.<class>', name: '<sample>' ]

      // INFO: pass -> detach -> orphan finder (hq / scraps) -> collapse -> bigBed
      // INFO: duplicates and scraps get their own bigBed
      ch_bucket_pass
          .map { meta, pass -> [ meta.name, meta, pass ] }    
          .groupTuple()                                      
          .map { name, metas, files ->
              [ [ id: name + '.pass', name: name ], files ]
          }
          .set { ch_pass }
      JOIN_VEREDICT_PASSES(ch_pass, 'bed')
      DETACH_PASS_DUPLICATES(JOIN_VEREDICT_PASSES.out.output)
      ISOTOOLS_ORPHAN_FINDER(DETACH_PASS_DUPLICATES.out.pass, ch_reference_transcripts, splice_scores)

      // INFO: scraps with a reference CDS go back to hq (#RSS); hq and scraps are rewritten under the
      // INFO: orphan finder's file names (12_POLISH/BED/orphans); empty files are dropped as before
      if (rescue_artifacts) {
          ISOTOOLS_ORPHAN_FINDER.out.hq
              .map { meta, hq -> [ meta.name, hq ] }
              .join(
                ISOTOOLS_ORPHAN_FINDER.out.scraps.map { meta, scraps -> [ meta.name, scraps ] },
                remainder: true
              )
              .map { name, hq, scraps -> [ [ id: name + '.pass', name: name ], [ hq, scraps ].findAll { it }, [], [] ] }
              .set { ch_rescue_orphans }

          RESCUE_FLNC_FROM_ORPHANS(
            ch_rescue_orphans,
            ch_reference_transcripts
          )

          ch_hq     = RESCUE_FLNC_FROM_ORPHANS.out.pass
          ch_scraps = RESCUE_FLNC_FROM_ORPHANS.out.scraps.map { meta, bed -> [ meta + [ id: meta.id + '.scraps' ], bed ] }
          ch_rescue_audit = ch_rescue_audit.mix(RESCUE_FLNC_FROM_ORPHANS.out.audit)
          ch_versions = ch_versions.mix(RESCUE_FLNC_FROM_ORPHANS.out.versions)
      } else {
          ch_hq     = ISOTOOLS_ORPHAN_FINDER.out.hq
          ch_scraps = ISOTOOLS_ORPHAN_FINDER.out.scraps
      }

      COLLAPSE_PASSES(ch_hq.filter { meta, bed -> bed.size() > 0 })
      BEDTOBIGBED_PASSES(COLLAPSE_PASSES.out.collapsed, chrom_sizes, autosql)
      BEDTOBIGBED_DUPLICATES(DETACH_PASS_DUPLICATES.out.duplicates, chrom_sizes, autosql)
      BEDTOBIGBED_SCRAPS(ch_scraps.filter { meta, bed -> bed.size() > 0 }, chrom_sizes, autosql)

      // INFO: one <sample>.rescued.tsv per sample (header once): id, original_id, origin, destination, rule, reference
      ch_rescue_audit
          .collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${params.global_output_dir}/12_POLISH/BED") { meta, tsv ->
              [ "${meta.name}.rescued.tsv", tsv ]
          }

      ch_bucket_trash
          .map { meta, trash -> [ meta.name, meta, trash ] }    
          .groupTuple()                                      
          .map { name, metas, files ->
               [ [ id: name + '.trash', name: name ], files ]
          }
          .set { ch_trash }
      JOIN_VEREDICT_TRASH(ch_trash, 'bed')
      DETACH_TRASH_DUPLICATES(JOIN_VEREDICT_TRASH.out.output)
      BEDTOBIGBED_TRASH(DETACH_TRASH_DUPLICATES.out.pass, chrom_sizes, autosql)

      ch_bucket_retentions
          .map { meta, retentions -> [ meta.name, meta, retentions ] }    
          .groupTuple()                                      
          .map { name, metas, files ->
                [ [ id: name + '.retentions', name: name ], files ]
          }
          .set { ch_retentions }
      JOIN_VEREDICT_RETENTIONS(ch_retentions, 'bed')
      DETACH_RETENTION_DUPLICATES(JOIN_VEREDICT_RETENTIONS.out.output)
      BEDTOBIGBED_RETENTIONS(DETACH_RETENTION_DUPLICATES.out.pass, chrom_sizes, autosql)

      ch_bucket_truncations
          .map { meta, truncations -> [ meta.name, meta, truncations ] }    
          .groupTuple()                                      
          .map { name, metas, files ->
               [ [ id: name + '.truncations', name: name ], files ]
          }
          .set { ch_truncations }
      JOIN_VEREDICT_TRUNCATIONS(ch_truncations, 'bed')
      DETACH_TRUNCATION_DUPLICATES(JOIN_VEREDICT_TRUNCATIONS.out.output)
      BEDTOBIGBED_TRUNCATIONS(DETACH_TRUNCATION_DUPLICATES.out.pass, chrom_sizes, autosql)

      ch_bucket_intraprimming
          .map { meta, intraprimming -> [ meta.name, meta, intraprimming ] }    
          .groupTuple()                                      
          .map { name, metas, files ->
               [ [ id: name + '.intraprimming', name: name ], files ]
          }
          .set { ch_intraprimming }
      JOIN_VEREDICT_INTRAPRIMMING(ch_intraprimming, 'bed')
      DETACH_INTRAPRIMMING_DUPLICATES(JOIN_VEREDICT_INTRAPRIMMING.out.output)
      BEDTOBIGBED_INTRAPRIMMING(DETACH_INTRAPRIMMING_DUPLICATES.out.pass, chrom_sizes, autosql)

      ch_bucket_rt
          .map { meta, rt -> [ meta.name, meta, rt ] }    
          .groupTuple()                                      
          .map { name, metas, files ->
               [ [ id: name + '.rt', name: name ], files ]
          }
          .set { ch_rt }
      JOIN_VEREDICT_RT(ch_rt, 'bed')
      DETACH_RT_DUPLICATES(JOIN_VEREDICT_RT.out.output)
      BEDTOBIGBED_RT(DETACH_RT_DUPLICATES.out.pass, chrom_sizes, autosql)

      ch_bucket_artifacts
          .map { meta, artifacts -> [ meta.name, meta, artifacts ] }
          .groupTuple()
          .map { name, metas, files ->
               [ [ id: name + '.artifacts', name: name ], files ]
          }
          .set { ch_artifacts }
      JOIN_VEREDICT_ARTIFACTS(ch_artifacts, 'bed')
      DETACH_ARTIFACT_DUPLICATES(JOIN_VEREDICT_ARTIFACTS.out.output)
      BEDTOBIGBED_ARTIFACTS(DETACH_ARTIFACT_DUPLICATES.out.pass, chrom_sizes, autosql)

      /*
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          PUBLISH BIGBEDS
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: all bigBeds of one sample -> [ [ id: '<sample>' ], [ *.bb ] ] for 12_POLISH/BB
      ch_bbs = Channel.empty()
      ch_bbs = ch_bbs.mix(BEDTOBIGBED_PASSES.out.bigbed)
      ch_bbs = ch_bbs.mix(BEDTOBIGBED_DUPLICATES.out.bigbed)
      ch_bbs = ch_bbs.mix(BEDTOBIGBED_SCRAPS.out.bigbed)
      ch_bbs = ch_bbs.mix(BEDTOBIGBED_TRASH.out.bigbed)
      ch_bbs = ch_bbs.mix(BEDTOBIGBED_RETENTIONS.out.bigbed)
      ch_bbs = ch_bbs.mix(BEDTOBIGBED_TRUNCATIONS.out.bigbed)
      ch_bbs = ch_bbs.mix(BEDTOBIGBED_INTRAPRIMMING.out.bigbed)
      ch_bbs = ch_bbs.mix(BEDTOBIGBED_RT.out.bigbed)
      ch_bbs = ch_bbs.mix(BEDTOBIGBED_ARTIFACTS.out.bigbed)
      ch_bbs = ch_bbs.mix(BEDTOBIGBED_UNPREDICTED.out.bigbed)

      ch_bbs.map { meta, file -> [meta.name, meta, file] }
         .groupTuple()
         .map { name, metas, files -> [ [ id: name ], files] }
         .set { ch_bbs }
      PUBLISH_BIGBEDS(ch_bbs)

    emit:
      pass                  = BEDTOBIGBED_PASSES.out.bigbed
      duplicates            = BEDTOBIGBED_DUPLICATES.out.bigbed
      scraps                = BEDTOBIGBED_SCRAPS.out.bigbed
      trash                 = BEDTOBIGBED_TRASH.out.bigbed
      retentions            = BEDTOBIGBED_RETENTIONS.out.bigbed
      truncations           = BEDTOBIGBED_TRUNCATIONS.out.bigbed
      intraprimming         = BEDTOBIGBED_INTRAPRIMMING.out.bigbed
      rt                    = BEDTOBIGBED_RT.out.bigbed
      artifacts             = BEDTOBIGBED_ARTIFACTS.out.bigbed
      unpredicted           = BEDTOBIGBED_UNPREDICTED.out.bigbed
      nmd                   = ch_bucket_nmd // [ meta, nmd.bed ] per sample@chr, without rescued rows
      bigbeds               = ch_bbs
      additional_columns    = ISOTOOLS_PLUGIN_VEREDICT.out.additional_bed_columns
      sample                = ch_bbs.map { meta, file -> meta.id }
      versions              = ch_versions.mix(RESCUE_FLNC_FROM_UNPREDICTED_ORFS.out.versions)
}
