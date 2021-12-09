#!/usr/bin/env nextflow

/*
 * pipeline input parameters
 */

params.NCBI_tmap = '/home/wslab/test1_annotate/nextflow_scripts/NCBI_compare.stringtie_chr33.gtf.tmap'     
params.NCBI_transcripts = '/home/wslab/test1_annotate/nextflow_scripts/NCBI_transcripts.fa'
params.genome_name = 'galGal6'
params.conda = '/home/wslab/test1_annotate/nextflow_scripts/environment.yml'
params.outdir ='./'

println """\
	 I S O F O R M - I D E N T I F I C A T I O N   P I P E L I N E
         =============================================================
         NCBI_tmap         : ${params.NCBI_tmap}
         NCBI_transcripts  : ${params.NCBI_transcripts}
         genome_name       : ${params.genome_name}    
         environment       : ${params.conda}
         outdir            : ${params.outdir}
         """
         .stripIndent()

process isoform_identification {
  echo true
  stageInMode 'copy'
  conda "${params.conda}"

  publishDir "${params.outdir}", mode: 'copy'

  shell:
  '''
  isoform-identification -m "!{params.NCBI_tmap}" -t "!{params.NCBI_transcripts}" -g "!{params.genome_name}"
  '''
}

