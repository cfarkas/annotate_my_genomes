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
         A N N O T A T E - M Y - G E N O M E S   P I P E L I N E
         =======================================================
         stringtie         : ${params.stringtie}
         ref_annotation    : ${params.ref_annotation}    
         genome            : ${params.genome}
         config_file       : ${params.config}
         threads           : ${params.threads}
         outdir            : ${params.outdir}
         environment       : ${params.conda}
         """
         .stripIndent()

process annotate_my_genomes {
  echo true
  stageInMode 'copy'
  conda "${params.conda}"

  shell:
  '''
  annotate-my-genomes -a "!{params.stringtie}" -r "!{params.ref_annotation}" -g "!{params.genome}" -c "!{params.config}" -t "!{params.threads}" -o "!{params.outdir}"
  '''
}

