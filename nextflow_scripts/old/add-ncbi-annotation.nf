#!/usr/bin/env nextflow

/*
 * pipeline input parameters
 */

params.stringtie = '/home/wslab/test1_annotate/nextflow_scripts/stringtie_chr33.gtf'     
params.NCBI_annotation = '/home/wslab/test1_annotate/nextflow_scripts/galGal6_ncbiRefSeq.gtf'
params.ref_annotation = '/home/wslab/test1_annotate/nextflow_scripts/galGal6.gtf'
params.genome = '/home/wslab/test1_annotate/nextflow_scripts/galGal6.fa'
params.config = '/home/wslab/test1_annotate/nextflow_scripts/gawn_config.sh'
params.threads = '10'
params.outdir = '/home/wslab/test1_annotate/nextflow_scripts/'
params.conda = '/home/wslab/test1_annotate/nextflow_scripts/environment.yml'

println """\
	 A D D - N C B I - A N N O T A T I O N   P I P E L I N E
         =======================================================
         stringtie         : ${params.stringtie}
         NCBI_annotation   : ${params.NCBI_annotation}
         ref_annotation    : ${params.ref_annotation}    
         genome            : ${params.genome}
         config_file       : ${params.config}
         threads           : ${params.threads}
         outdir            : ${params.outdir}
         environment       : ${params.conda}
         """
         .stripIndent()

process add_ncbi_annotation {
  echo true
  stageInMode 'copy'
  conda "${params.conda}"

  shell:
  '''
  add-ncbi-annotation -a "!{params.stringtie}" -n "!{params.NCBI_annotation}" -r "!{params.ref_annotation}" -g "!{params.genome}" -c "!{params.config}" -t "!{params.threads}" -o "!{params.outdir}"
  '''
}
