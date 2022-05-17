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


process check_paths {
  echo true
  stageInMode 'copy'
  conda "${params.conda}"

  output:
  val '${m_DIR}' into records1
  val '${t_DIR}' into records2
  val '${c_DIR}' into records3
  val '${o_DIR}' into records4

  shell:
  '''
  # mandatory arguments
  if [ ! "!{params.NCBI_tmap}" ] || [ ! "!{params.NCBI_transcripts}" ] || [ ! "!{params.genome_name}" ] || [ ! "!{params.conda}" ] || [ ! "!{params.outdir}" ]; then
    echo ""
    echo "arguments -m, -t, -g, -c and -o must be provided"
    echo ""
    echo "$usage" >&2; exit 1
  fi

  # Conditions : output folder
  if [ ! -d "!{params.outdir}" ]; then
    echo ""
    echo "Output directory: !{params.outdir} not found. Please create the output directory first, before running the pipeline."
    echo ""
    exit 9999 # die with error code 9999
  fi

  # Conditions : Input existance
  if [ ! -e "!{params.NCBI_tmap}" ]; then
    echo ""
    echo "!{params.NCBI_tmap} does not exist. Check your -m input"
    echo ""
    exit 9999 # die with error code 9999
  fi

  if [ ! -e "!{params.NCBI_transcripts}" ]; then
    echo ""
    echo "!{params.NCBI_transcripts} does not exist. Check your -t input"
    echo ""
    exit 9999 # die with error code 9999
  fi

  if [ ! -e "!{params.conda}" ]; then
    echo ""
    echo "!{params.conda} does not exist. Check your -c input"
    echo ""
    exit 9999 # die with error code 9999
  fi

  # Conditions : Getting absolute path of inputs
  echo ""
  m_DIR="$( cd "$( dirname "!{params.NCBI_tmap}" )" && pwd )"
  echo ""
  echo "::: The absolute path of -m is $m_DIR"
  echo ""
  t_DIR="$( cd "$( dirname "!{params.NCBI_transcripts}" )" && pwd )"
  echo ""
  echo "::: The absolute path of -t is $t_DIR"
  echo ""
  c_DIR="$( cd "$( dirname "!{params.conda}" )" && pwd )"
  echo ""
  echo "::: The absolute path of -c is $c_DIR"
  echo ""
  o_DIR="$( cd "$( dirname "!{params.outdir}" )" && pwd )"
  echo ""
  echo "::: The absolute path of -o is $o_DIR"
  echo ""
  '''
}


process check_inputs {
  echo true
  stageInMode 'copy'
  conda "${params.conda}"

  output:
  val '${NCBI_tmap}' into records5
  val '${NCBI_transcripts}' into records6
  val '${genome_name}' into records7
  val '${anaconda_env}' into records8
  val '${outdir}' into records9
  file 'ncbiRefSeqLink.txt' into ncbiRefSeqLink

  shell:
  '''
  printf "::: Defining Variables :::\n"
  echo""
  FILE1="!{params.NCBI_tmap}"
  basename "$FILE1"
  NCBI_tmap="$(basename -- $FILE1)"
  echo "The NCBI tmap file used as input is the following: $NCBI_tmap"
  echo ""
  FILE2="!{params.NCBI_transcripts}"
  basename "$FILE2"
  NCBI_transcripts="$(basename -- $FILE2)"
  echo "The NCBI transcripts used as input is the following: $NCBI_transcripts"
  echo ""
  FILE3="!{params.conda}"
  basename "$FILE3"
  anaconda_env="$(basename -- $FILE3)"
  echo "The anaconda environment file is the following: $anaconda_env"
  echo ""
  FILE4="!{params.outdir}"
  basename "$FILE4"
  outdir="$(basename -- $FILE4)"
  echo "The outdir folder name is the following: $outdir"
  echo ""
  genome_name="!{params.genome_name}"

  if [ -f ncbiRefSeqLink.txt ]; then
      echo "::: ncbiRefSeqLink.txt file found. Continue:"
      echo ""
      :
  else
      echo "::: Downloading ncbiRefSeqLink.txt file"
      wget http://hgdownload.cse.ucsc.edu/goldenpath/${genome_name}/database/ncbiRefSeqLink.txt.gz
      gunzip ncbiRefSeqLink.txt.gz
      echo ""
      echo "Number of lines in ncbiRefSeqLink.txt:"
      cat ncbiRefSeqLink.txt | wc -l
      echo "Continue with python processing steps:"
      echo ""
  fi
  '''
}


process python_inputs {
  echo true
  stageInMode 'copy'
  conda "${params.conda}"

  input:
  val '${m_DIR}' from records1
  val '${t_DIR}' from records2
  val '${NCBI_tmap}' from records5
  val '${NCBI_transcripts}' from records6

  output:
  file 'stringtie_for_script.tmap' into records10
  file 'transcripts_Isoform2.tab' into records11

  shell:
  '''
  # Inputs for python
  cp "!{params.NCBI_tmap}" ./stringtie_for_script.tmap
  seqkit fx2tab "!{params.NCBI_transcripts}" > transcripts_Isoform.tab
  # Formatting transcripts_Isoform.tab if gene= is present in file
  sed -i 's/gene=/\t/'g transcripts_Isoform.tab
  awk '{print $1"\t"$NF}' transcripts_Isoform.tab > transcripts_Isoform2.tab
  '''
}


process gffcompare_parser {
  echo true
  stageInMode 'copy'
  conda "${params.conda}"

  input:
  file 'ncbiRefSeqLink.txt' from ncbiRefSeqLink
  file 'stringtie_for_script.tmap' from records10
  file 'transcripts_Isoform2.tab' from records11

  output:
  file 'Ref_Transcript_Annotation.csv' into records12
  file 'Novel_Transcript_Annotation.csv' into records13

  shell:
  '''
  python << END

  import sys
  import pandas as pd
  df = pd.read_csv('stringtie_for_script.tmap', sep = '\t')
  print(df.sample(10))
  print("Total number of transcripts:", df.shape[0])
  print("")
  df2 = df[~df.ref_id.astype(str).str.contains('-')]
  novel_transcripts = df[df.ref_id.astype(str).str.contains('-')]
  df3 = df2[["ref_gene_id", "ref_id", "class_code", "qry_gene_id", "qry_id", "num_exons", "FPKM", "TPM"]]
  df_novel_transcripts = novel_transcripts[["ref_gene_id", "ref_id", "class_code", "qry_gene_id", "qry_id", "num_exons", "FPKM", "TPM"]]
  print("Reference transcripts:")
  print(df3.sample(10))
  print("")
  print("Novel transcripts:")
  print(df_novel_transcripts.sample(10))
  print("")
  colnames=['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18']
  dfA1 = pd.read_csv('ncbiRefSeqLink.txt', sep = '\t', low_memory=False, names=colnames, header=None)
  print(dfA1.head(10))
  dfA2 = dfA1[['0', '1', '2', '3', '5', '14', '16']]
  dfA2 = dfA2.rename(columns={'0': 'ref_id', '1': 'Annotation Status', '2' : 'NCBI RefSeq Gene ID', '3' : 'Transcript Description', '5' : 'NCBI RefSeq Protein ID', '14' : 'Alternative Gene Name', '16' : 'RefSeq Transcript Info'})
  print("ncbiRefSeqLink annotation:")
  print(dfA2.sample(10))
  print("")
  colnames = ['qry_id', 'cds_seq', 'none']
  cds = pd.read_csv('transcripts_Isoform2.tab', sep = '\t', names=colnames)
  cds2 = cds[["qry_id", "cds_seq"]]
  print("transcripts file:")
  print(cds2.sample(10))
  print("")
  result1 = pd.merge(df3, dfA2, on='ref_id', how='inner')
  result1.sample(10)
  result2 = pd.merge(result1, cds2, on='qry_id', how='inner')
  result2.sample(10)
  result3 = pd.merge(df_novel_transcripts, cds2, on='qry_id', how='inner')
  result3.sample(10)
  print("Number of Joined Transcripts (reference):", result2.shape[0])
  print("")
  print("Number of Joined Transcripts (novel):", result3.shape[0])
  print("")
  result2.to_csv('Ref_Transcript_Annotation.csv', index=False)
  result3.to_csv('Novel_Transcript_Annotation.csv', index=False)
  print("::: Done. Ref_Transcript_Annotation.csv and Novel_Transcript_Annotation.csv were succesfully produced")
  print("")
  END
  '''
}


process output_pipeline {
  echo true
  stageInMode 'copy'
  conda "${params.conda}"

  input:
  file 'Ref_Transcript_Annotation.csv' from records12
  file 'Novel_Transcript_Annotation.csv' from records13

  shell:
  '''
  echo ""
  printf "::: Moving results to the output directory :::\n"
  cp Ref_Transcript_Annotation.csv "!{params.outdir}"
  cp Novel_Transcript_Annotation.csv "!{params.outdir}"
  echo ""
  echo "------------------------------------------------------------"
  echo "------------------------------------------------------------"
  echo "::: INFO: all done" :::
  echo ""
  echo "The following files are available in the output directory : "
  echo ""
  echo "Ref_Transcript_Annotation.csv contains annotation and coordinates of known transcripts"
  echo ""
  echo "Novel_Transcript_Annotation.csv contains collection of novel transcripts"
  echo ""
  echo "------------------------------------------------------------"
  echo "------------------------------------------------------------"
  echo ""
  '''
}
