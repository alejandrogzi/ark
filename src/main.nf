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

    if (!(params.isoseq_cluster2_mode in ['per_sample', 'multi_sample', 'both'])) { 
      problems << 'ERROR: Unknown isoseq_cluster2_mode option -> options are: per_sample, multi_sample, both'
    }

    if (!(params.aligner in ['mm2', 'ultra', 'desalt', 'pbmm2', 'flair', 'ark'])) { 
      problems << 'ERROR: Unknown aligner option -> options are: mm2, ultra, desalt, pbmm2, flair, ark'
    }

    if (problems) {
        error "Parameter validation failed:\n  - " + problems.join('\n  - ')
        System.exit(1)
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
