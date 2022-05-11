#!/usr/bin/env nextflow

/* 
 * pipeline input parameters 
 */


params.genome = 'galGal6'
params.outdir = './'
params.conda = '/home/wslab/test1_annotate/environment.yml'

println """\
         G E N O M E - D O W N L O A D   P I P E L I N E    
         ===============================================
         genome        : ${params.genome}
         outdir        : ${params.outdir}
         environment   : ${params.conda}
         """
         .stripIndent()

process genome_download {
  echo true
  conda "${params.conda}"
  
  publishDir "${params.outdir}", mode: 'copy'

  output:
  file "${params.genome}*"
  
  shell:
  '''
  #!/usr/bin/env bash
  
  genome="!{params.genome}"
  wget http://hgdownload.cse.ucsc.edu/goldenpath/${genome}/bigZips/${genome}.2bit
  wget http://hgdownload.soe.ucsc.edu/goldenPath/${genome}/database/refGene.txt.gz
  wget http://hgdownload.soe.ucsc.edu/goldenPath/${genome}/database/ncbiRefSeq.txt.gz
  if [ -f twoBitToFa ]; then
      echo "twoBitToFa script found. Continue:"
      echo ""
      : 
  else
      echo "Downloading twoBitToFa script"
      wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa 
  fi
  chmod 755 twoBitToFa
  ./twoBitToFa ${genome}.2bit ${genome}.fa
  samtools faidx ${genome}.fa

  if [ -f genePredToGtf ]; then
      echo "genePredToGtf script found. Continue:"
      echo ""
      : 
  else
      echo "Downloading genePredToGtf script"
      wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf 
  fi
  chmod 755 genePredToGtf
  gunzip refGene.txt.gz
  gunzip ncbiRefSeq.txt.gz
  cut -f 2- refGene.txt | ./genePredToGtf file stdin -source=${genome}_Ref  ${genome}.gtf
  cut -f 2- ncbiRefSeq.txt | ./genePredToGtf file stdin -source=${genome}_Ref  ${genome}_ncbiRefSeq.gtf
  echo ""
  echo "All done. ${genome} FASTA and GTF files are located in the current working directory (or specified directory with --outdir)"
  '''
}

