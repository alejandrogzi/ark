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

include { GAWK_JOIN as JOIN_VEREDICT_PASSES } from '../../modules/custom/gawk/join/main.nf'
include { GAWK_JOIN as JOIN_VEREDICT_TRASH } from '../../modules/custom/gawk/join/main.nf'
include { GAWK_JOIN as JOIN_VEREDICT_RETENTIONS } from '../../modules/custom/gawk/join/main.nf'
include { GAWK_JOIN as JOIN_VEREDICT_TRUNCATIONS } from '../../modules/custom/gawk/join/main.nf'
include { GAWK_JOIN as JOIN_VEREDICT_INTRAPRIMMING } from '../../modules/custom/gawk/join/main.nf'
include { GAWK_JOIN as JOIN_VEREDICT_RT } from '../../modules/custom/gawk/join/main.nf'

include { BEDTOBIGBED as BEDTOBIGBED_PASSES } from '../../modules/custom/bigtools/bedtobigbed/main.nf'
include { BEDTOBIGBED as BEDTOBIGBED_TRASH } from '../../modules/custom/bigtools/bedtobigbed/main.nf'
include { BEDTOBIGBED as BEDTOBIGBED_RETENTIONS } from '../../modules/custom/bigtools/bedtobigbed/main.nf'
include { BEDTOBIGBED as BEDTOBIGBED_TRUNCATIONS } from '../../modules/custom/bigtools/bedtobigbed/main.nf'
include { BEDTOBIGBED as BEDTOBIGBED_INTRAPRIMMING } from '../../modules/custom/bigtools/bedtobigbed/main.nf'
include { BEDTOBIGBED as BEDTOBIGBED_DUPLICATES } from '../../modules/custom/bigtools/bedtobigbed/main.nf'
include { BEDTOBIGBED as BEDTOBIGBED_SCRAPS } from '../../modules/custom/bigtools/bedtobigbed/main.nf'
include { BEDTOBIGBED as BEDTOBIGBED_RT } from '../../modules/custom/bigtools/bedtobigbed/main.nf'

include { PUBLISH as PUBLISH_BIGBEDS } from '../../modules/custom/publish/main.nf'

include { DETACH_DUPLICATES as DETACH_PASS_DUPLICATES } from '../../modules/custom/detach/main.nf'
include { DETACH_DUPLICATES as DETACH_RETENTION_DUPLICATES } from '../../modules/custom/detach/main.nf'
include { DETACH_DUPLICATES as DETACH_TRUNCATION_DUPLICATES } from '../../modules/custom/detach/main.nf'
include { DETACH_DUPLICATES as DETACH_INTRAPRIMMING_DUPLICATES } from '../../modules/custom/detach/main.nf'
include { DETACH_DUPLICATES as DETACH_RT_DUPLICATES } from '../../modules/custom/detach/main.nf'
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
      annotation             // channel: [ val(meta), [ annotation ] ]
      forward_peaks          // channel: [ val(meta), [ forward_peaks ] ]
      reverse_peaks          // channel: [ val(meta), [ reverse_peaks ] ]
      chrom_sizes            // path
      splice_scores          // channel: [ val(meta), [ splice_scores ] ]
      autosql                // channel: [ autosql ]
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
          VEREDICT [ PASS, TRASH, RETENTIONS, TRUNCATIONS, INTRAPRIMING, RT ]
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
          PER-SAMPLE TRACKS [ JOIN -> DETACH -> BEDTOBIGBED ]
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      */

      // INFO: VEREDICT runs per sample@chromosome; each class is regrouped by meta.name
      // INFO: so JOIN concatenates + sorts all chromosomes into one BED per sample,
      // INFO: DETACH splits '#DU'-tagged duplicates out and BEDTOBIGBED builds one track
      // INFO: meta from here on -> [ id: '<sample>.<class>', name: '<sample>' ]

      // INFO: pass -> detach -> orphan finder (hq / scraps) -> collapse -> bigBed
      // INFO: duplicates and scraps get their own bigBed
      ISOTOOLS_PLUGIN_VEREDICT.out.pass
          .map { meta, pass -> [ meta.name, meta, pass ] }    
          .groupTuple()                                      
          .map { name, metas, files ->
              [ [ id: name + '.pass', name: name ], files ]
          }
          .set { ch_pass }
      JOIN_VEREDICT_PASSES(ch_pass, 'bed')
      DETACH_PASS_DUPLICATES(JOIN_VEREDICT_PASSES.out.output)
      ISOTOOLS_ORPHAN_FINDER(DETACH_PASS_DUPLICATES.out.pass, ch_reference_transcripts, splice_scores)
      COLLAPSE_PASSES(ISOTOOLS_ORPHAN_FINDER.out.hq)
      BEDTOBIGBED_PASSES(COLLAPSE_PASSES.out.collapsed, chrom_sizes, autosql)
      BEDTOBIGBED_DUPLICATES(DETACH_PASS_DUPLICATES.out.duplicates, chrom_sizes, autosql)
      BEDTOBIGBED_SCRAPS(ISOTOOLS_ORPHAN_FINDER.out.scraps, chrom_sizes, autosql)

      ISOTOOLS_PLUGIN_VEREDICT.out.trash
          .map { meta, trash -> [ meta.name, meta, trash ] }    
          .groupTuple()                                      
          .map { name, metas, files ->
               [ [ id: name + '.trash', name: name ], files ]
          }
          .set { ch_trash }
      JOIN_VEREDICT_TRASH(ch_trash, 'bed')
      DETACH_TRASH_DUPLICATES(JOIN_VEREDICT_TRASH.out.output)
      BEDTOBIGBED_TRASH(DETACH_TRASH_DUPLICATES.out.pass, chrom_sizes, autosql)

      ISOTOOLS_PLUGIN_VEREDICT.out.retentions
          .map { meta, retentions -> [ meta.name, meta, retentions ] }    
          .groupTuple()                                      
          .map { name, metas, files ->
                [ [ id: name + '.retentions', name: name ], files ]
          }
          .set { ch_retentions }
      JOIN_VEREDICT_RETENTIONS(ch_retentions, 'bed')
      DETACH_RETENTION_DUPLICATES(JOIN_VEREDICT_RETENTIONS.out.output)
      BEDTOBIGBED_RETENTIONS(DETACH_RETENTION_DUPLICATES.out.pass, chrom_sizes, autosql)

      ISOTOOLS_PLUGIN_VEREDICT.out.truncations
          .map { meta, truncations -> [ meta.name, meta, truncations ] }    
          .groupTuple()                                      
          .map { name, metas, files ->
               [ [ id: name + '.truncations', name: name ], files ]
          }
          .set { ch_truncations }
      JOIN_VEREDICT_TRUNCATIONS(ch_truncations, 'bed')
      DETACH_TRUNCATION_DUPLICATES(JOIN_VEREDICT_TRUNCATIONS.out.output)
      BEDTOBIGBED_TRUNCATIONS(DETACH_TRUNCATION_DUPLICATES.out.pass, chrom_sizes, autosql)

      ISOTOOLS_PLUGIN_VEREDICT.out.intraprimming
          .map { meta, intraprimming -> [ meta.name, meta, intraprimming ] }    
          .groupTuple()                                      
          .map { name, metas, files ->
               [ [ id: name + '.intraprimming', name: name ], files ]
          }
          .set { ch_intraprimming }
      JOIN_VEREDICT_INTRAPRIMMING(ch_intraprimming, 'bed')
      DETACH_INTRAPRIMMING_DUPLICATES(JOIN_VEREDICT_INTRAPRIMMING.out.output)
      BEDTOBIGBED_INTRAPRIMMING(DETACH_INTRAPRIMMING_DUPLICATES.out.pass, chrom_sizes, autosql)

      ISOTOOLS_PLUGIN_VEREDICT.out.rt 
          .map { meta, rt -> [ meta.name, meta, rt ] }    
          .groupTuple()                                      
          .map { name, metas, files ->
               [ [ id: name + '.rt', name: name ], files ]
          }
          .set { ch_rt }
      JOIN_VEREDICT_RT(ch_rt, 'bed')
      DETACH_RT_DUPLICATES(JOIN_VEREDICT_RT.out.output)
      BEDTOBIGBED_RT(DETACH_RT_DUPLICATES.out.pass, chrom_sizes, autosql)

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
      bigbeds               = ch_bbs
      additional_columns    = ISOTOOLS_PLUGIN_VEREDICT.out.additional_bed_columns
      sample                = ch_bbs.map { meta, file -> meta.id }
      versions              = ch_versions
}
