#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS/WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ARK as MAIN } from './workflows/ark.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATION FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def validateFullRun() {
    def problems = []

    if (!params.global_input_dir) { problems << 'missing required --global_input_dir' }
    if (!params.global_output_dir) { problems << 'missing required --global_output_dir' }
    if (!params.global_genome) { problems << 'missing required --global_genome' }
    if (!params.global_annotation) { problems << 'missing required --global_annotation' }
    if (!params.global_repeats) { problems << 'missing required --global_repeats' }
    
    if (!(params.entrypoint in ['subreads', 'ccs', 'refine', 'cluster', 'flnc'])) {
      problems << 'ERROR: Unknown entrypoint option -> options are: subreads, ccs, refine, cluster, flnc'
    }

    if (params.entrypoint in ['subreads', 'ccs', 'refine'] && !params.global_primers) {
      problems << 'missing required --global_primers'
    }

    // INFO: an old params file would otherwise fall back silently to the cluster_mode default
    if (params.containsKey('isoseq_cluster2_mode')) {
      problems << 'isoseq_cluster2_mode was renamed to cluster_mode'
    }

    if (!(params.cluster_mode in ['per_sample', 'multi_sample', 'both'])) {
      problems << 'ERROR: Unknown cluster_mode option -> options are: per_sample, multi_sample, both'
    }

    if (params.cluster_engine in ['cdhit', 'rattle'] || params.containsKey('cdhit_identity')) {
      problems << 'cluster_engine cdhit/rattle (and cdhit_identity) were removed in v2.1.0 -> use none (default) or isoseq'
    } else if (!(params.cluster_engine in ['none', 'isoseq'])) {
      problems << 'ERROR: Unknown cluster_engine option -> options are: none, isoseq'
    }
    if (params.reconstruct_engine in ['isoquant', 'isocall']) {
      problems << "reconstruct_engine ${params.reconstruct_engine} is reserved but not implemented in v2.1.0 -> use chain"
    } else if (!(params.reconstruct_engine in ['none', 'chain'])) {
      problems << 'ERROR: Unknown reconstruct_engine option -> options are: chain, none'
    }
    if (!(params.reconstruct_preset in ['sensitive', 'balanced', 'strict'])) {
      problems << 'ERROR: Unknown reconstruct_preset option -> options are: sensitive, balanced, strict'
    }
    if (!(params.flnc_input_state in ['auto', 'ccs', 'fl', 'flnc', 'mixed', 'clustered'])) {
      problems << 'ERROR: Unknown flnc_input_state option -> options are: auto, ccs, fl, flnc, mixed, clustered'
    }
    // INFO: auto picks a bundled primer set from the detected kit; a forced ccs state skips that detection
    if (params.entrypoint == 'flnc' && params.flnc_input_state == 'ccs' && !params.global_primers) {
      problems << 'flnc_input_state ccs needs --global_primers (or leave it on auto to detect the kit)'
    }

    if (!(params.aligner in ['mm2', 'ultra', 'desalt', 'pbmm2', 'flair', 'ark'])) { 
      problems << 'ERROR: Unknown aligner option -> options are: mm2, ultra, desalt, pbmm2, flair, ark'
    }

    if (problems) {
        error "Parameter validation failed:\n  - " + problems.join('\n  - ')
        System.exit(1)
    }
    if (params.reconstruct_engine == 'none') {
      log.warn 'reconstruct_engine none: every segmented read reaches ORF prediction'
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ARK {  
    validateFullRun()

    log.info """
    > ark v${workflow.manifest.version}
    > A Reference pipeline to annotate euKaryotes at high resolution
    > The Hiller Lab at the Senckenberg Research Institute
  
    Authors : ${workflow.manifest.author}
    Github  :  ${workflow.manifest.homePage}

      Entrypoint: ${params.entrypoint}
      Input     : ${params.global_input_dir}
      Output    : ${params.global_output_dir}
      Genome    : ${params.global_genome}
      Annotation: ${params.global_annotation}
      Repeats   : ${params.global_repeats}
      Aligner   : ${params.aligner}
      Clustering: ${params.cluster_engine}
      Reconstruct: ${params.reconstruct_engine} (${params.reconstruct_preset})
      Database  : ${params.xorf_protein_database} (custom: ${params.xorf_custom_database})
      Profile   : ${workflow.profile}
    """.stripIndent()

    MAIN () 
}

workflow { ARK () }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
