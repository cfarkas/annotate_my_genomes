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


process check_paths {
  echo true
  stageInMode 'copy'
  conda "${params.conda}"

  output:
  val "${a_DIR}"
  val "${n_DIR}"
  val "${r_DIR}"
  val "${g_DIR}"
  val "${c_DIR}"
  val "${o_DIR}"

  shell:
  '''
  # mandatory arguments
  if [ ! "!{params.stringtie}" ] || [ ! "!{params.NCBI_annotation}" ] || [ ! "!{params.ref_annotation}" ] || [ ! "!{params.genome}" ] || [ ! "!{params.config}" ] || [ ! "!{params.threads}" ] || [ ! "!{params.conda}" ]|| [ ! "!{params.outdir}" ]; then
    echo ""
    echo "arguments -a, -n, -r, -g, -c, -t and -o must be provided"
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
  if [ ! -e "!{params.stringtie}" ]; then
    echo ""
    echo "!{params.stringtie} does not exist. Check your -a input"
    echo ""
    exit 9999 # die with error code 9999
  fi

  if [ ! -e "!{params.NCBI_annotation}" ]; then
    echo ""
    echo "!{params.NCBI_annotation} does not exist. Check your -n input"
    echo ""
    exit 9999 # die with error code 9999
  fi

  if [ ! -e "!{params.ref_annotation}" ]; then
    echo ""
    echo "!{params.ref_annotation} does not exist. Check your -r input"
    echo ""
    exit 9999 # die with error code 9999
  fi

  if [ ! -e "!{params.genome}" ]; then
    echo ""
    echo "!{params.genome} does not exist. Check your -g input"
    echo ""
    exit 9999 # die with error code 9999
  fi
  if [ ! -e "!{params.config}" ]; then
    echo ""
    echo "!{params.config} file does not exist. Check your -c input"
    echo ""
    exit 9999 # die with error code 9999
  fi

  # Conditions : Getting absolute path of inputs
  echo ""
  a_DIR="$( cd "$( dirname "!{params.stringtie}" )" && pwd )"
  echo ""
  echo "::: The absolute path of -a is $a_DIR"
  echo ""
  n_DIR="$( cd "$( dirname "!{params.NCBI_annotation}" )" && pwd )"
  echo ""
  echo "::: The absolute path of -n is $n_DIR"
  echo ""
  r_DIR="$( cd "$( dirname "!{params.ref_annotation}" )" && pwd )"
  echo ""
  echo "::: The absolute path of -r is $r_DIR"
  echo ""
  g_DIR="$( cd "$( dirname "!{params.genome}" )" && pwd )"
  echo ""
  echo "::: The absolute path of -g is $g_DIR"
  echo ""
  c_DIR="$( cd "$( dirname "!{params.config}" )" && pwd )"
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
  val "${stringtie_input}"
  val "${ncbi_reference_gtf}"
  val "${reference_gtf}"
  val "${reference_genome}"
  val "${gawn_config}"
  val "${threads}"
  val "${output_folder}"

  shell:
  '''
  begin=`date +%s`

  printf "::: Defining Variables :::\n"
  echo""
  FILE1="!{params.stringtie}"
  basename "$FILE1"
  stringtie_input="$(basename -- $FILE1)"
  echo "The stringtie file used as input is the following: $stringtie_input"
  echo ""
  FILE2="!{params.NCBI_annotation}"
  basename "$FILE2"
  ncbi_reference_gtf="$(basename -- $FILE2)"
  echo "The NCBI reference GTF used as input is the following: $ncbi_reference_gtf"
  echo ""
  FILE3="!{params.ref_annotation}"
  basename "$FILE3"
  reference_gtf="$(basename -- $FILE3)"
  echo "The reference GTF used as input is the following: $reference_gtf"
  echo ""
  FILE4="!{params.genome}"
  basename "$FILE4"
  reference_genome="$(basename -- $FILE4)"
  echo "The reference genome used as input is the following: $reference_genome"
  echo ""
  FILE5="!{params.config}"
  basename "$FILE5"
  gawn_config="$(basename -- $FILE5)"
  echo "The gawn_config file used as input is the following: $gawn_config"
  echo ""
  FILE6="!{params.threads}"
  basename "$FILE6"
  threads="$(basename -- $FILE6)"
  echo "The number of threads for calculation are the following: $threads"
  echo ""
  FILE7="!{params.outdir}"
  basename "$FILE7"
  output_folder="$(basename -- $FILE7)"
  echo "The output folder is the following: $output_folder"
  echo ""
  '''
}


process gffcompare_stringtie_fix {
  echo true
  stageInMode 'copy'
  conda "${params.conda}"

  output:
  file 'sed.script' into records0
  file 'final_annotated.gtf' into records1
  file 'NCBI_compare.stringtie_input_file.gtf.tmap' into tmap_records
  file 'gffcompare_outputs_NCBI' into gffcompare_outputs_NCBI_dir
  file 'Stats.txt' into Stats_records

  shell:
  '''
  echo ""
  echo "====  1. Overlapping StringTie transcripts with UCSC GTF:"
  echo ""
  cp "!{params.stringtie}" ./stringtie_input_file.gtf
  gffcompare -R -r "!{params.NCBI_annotation}" -s "!{params.genome}" -o NCBI_compare stringtie_input_file.gtf
  echo "::: INFO: Overlap was done."
  echo ""
  echo "::: INFO: Generating Stats"
  # Stats
  exec 3<> Stats.txt
  echo "Number of assembled genes:" >> Stats.txt
  cat NCBI_compare.stringtie_input_file.gtf.tmap | sed "1d" | cut -f4 | sort | uniq | wc -l >> Stats.txt
  echo "" >> Stats.txt
  echo "Number of novel genes:" >> Stats.txt
  cat NCBI_compare.stringtie_input_file.gtf.tmap | awk '$3=="u"{print $0}' | cut -f4 | sort | uniq | wc -l >> Stats.txt
  echo "" >> Stats.txt
  echo "Number of novel transcripts:" >> Stats.txt
  cat NCBI_compare.stringtie_input_file.gtf.tmap | awk '$3=="u"{print $0}' | cut -f5 | sort | uniq | wc -l >> Stats.txt
  echo "" >> Stats.txt
  echo "Number of transcripts matching annotation:" >> Stats.txt
  cat NCBI_compare.stringtie_input_file.gtf.tmap | awk '$3=="="{print $0}' | cut -f5 | sort | uniq | wc -l >> Stats.txt
  exec 3>&-
  echo ""
  printf "::: INFO: Stats were Done\n"
  echo ""
  echo ""
  echo "====  2. Replacing gene_id field in final_annotated.gtf file with UCSC gene_id's:"
  echo ""
  awk '{print $4"\t"$1}' NCBI_compare.stringtie_input_file.gtf.tmap > NCBI_compare.stringtie_input_file.gtf.tmap.1
  tail -n +2 NCBI_compare.stringtie_input_file.gtf.tmap.1 > NCBI_compare.stringtie_input_file.gtf.tmap.2
  awk '$2 != "-"' NCBI_compare.stringtie_input_file.gtf.tmap.2 > namelist
  awk '!a[$0]++' namelist > namelist_unique
  tac namelist_unique > namelist_unique_sorted
  rm namelist namelist_unique
  awk '{print $1}' namelist_unique_sorted  > A
  awk '{print $2}' namelist_unique_sorted  > B
  sed 's/^/"/' A > A.1
  sed 's/$/"/' A.1 > A.2
  sed 's/^/"/' B > B.1
  sed 's/$/"/' B.1 > B.2
  paste -d'\t' A.2 B.2 > namelist
  rm A A.1 A.2 B B.1 B.2
  ###############################
  # Getting gene names replaced #
  ###############################
  awk '{print $1}' namelist > fileA
  awk '{print $2}' namelist > fileB
  paste -d % fileA fileB > sed.script
  sed -i -e 's/^/s%/' sed.script
  sed -i -e 's/$/%/' sed.script
  cat stringtie_input_file.gtf | parallel --pipe -j "!{params.threads}" sed -f sed.script > final_annotated.gtf
  rm -f fileA fileB *tmap.1 *tmap.2
  # sorting GTF file
  rm -r -f gff3sort
  git clone https://github.com/cfarkas/gff3sort.git
  perl ./gff3sort/gff3sort.pl final_annotated.gtf > final_annotated.sorted.gtf
  rm final_annotated.gtf
  mv final_annotated.sorted.gtf final_annotated.gtf
  printf "::: INFO: Gene_id field was replaced in the StringTie.gtf file and final_annotated.gtf was generated with these changes\n"
  echo ""
  printf "::: INFO: Moving gffcompare results to gffcompare_outputs folder ...\n"
  rm -r -f gffcompare_outputs_NCBI
  mkdir gffcompare_outputs_NCBI
  mv *.loci *.stats *.refmap *.tmap *.tracking ./gffcompare_outputs_NCBI
  cp ./gffcompare_outputs_NCBI/NCBI_compare.stringtie_input_file.gtf.tmap ./
  echo ""
  printf "::: INFO: gffcompare_stringtie_fix process done\n"
  echo ""
  '''
}


process gffread_stringtie {
  echo true
  stageInMode 'copy'
  conda "${params.conda}"

  input:
  file 'final_annotated.gtf' from records1

  output:
  file 'NCBI_transcripts.fa' into records2

  shell:
  '''
  echo ""
  echo "====  3. Obtaining Transcripts in FASTA format with gffread"
  echo ""
  gffread -w NCBI_transcripts.fa -g "!{params.genome}" final_annotated.gtf
  printf "::: INFO: Number of transcripts:\n"
  grep ">" NCBI_transcripts.fa -c
  echo ""
  printf "::: INFO: gffread_stringtie process done\n"
  echo ""
  echo "====  4. We will perform gene annotation using GAWN pipeline:"
  echo ""
  '''
}


process gawn_annotation {
  echo true
  stageInMode 'copy'
  conda "${params.conda}"

  input:
  file 'NCBI_transcripts.fa' from records2

  output:
  file 'transcriptome.swissprot' into records3
  file 'transcriptome.hits' into records4
  file 'gawn' into gawn_dir

  shell:
  '''
  rm -r -f gawn
  git clone https://github.com/cfarkas/gawn.git
  cd gawn/02_infos/
  dir2=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
  echo ""
  cd ..
  cd ..
  dir1=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
  cd ${dir1}
  cp "!{params.genome}" ${dir1}/gawn/03_data/genome.fasta
  cp NCBI_transcripts.fa ${dir1}/gawn/03_data/transcriptome.fasta
  rm -rf ${dir2}/gawn_config.sh
  echo ""
  printf "::: INFO: Starting GAWN transcript annotation\n"
  echo ""
  cd ${dir1}/gawn/
  ./gawn "!{params.config}"
  echo ""
  printf "::: INFO: GAWN annotation is done. The novel transcripts were annotated in ./gawn/04_annotation/ \n"
  echo ""
  cd ${dir1}
  cp ${dir1}/gawn/04_annotation/transcriptome.swissprot ${dir1}
  cp ${dir1}/gawn/04_annotation/transcriptome.hits ${dir1}
  printf "::: INFO: Number of swissprot hits:\n"
  cat transcriptome.swissprot | wc -l
  echo ""
  printf "::: INFO: Number of transcriptome hits:\n"
  cat transcriptome.hits | wc -l
  echo ""
  printf "::: INFO: gawn_annotation process done. Transcriptome hits were succesfully extracted :::\n"
  echo ""
  echo "==== 5. We will identify long non-coding RNAs by using FEELnc: "
  echo ""
  '''
}


process feelnc_annotation {
  echo true
  stageInMode 'copy'
  conda "${params.conda}"

  input:
  file 'final_annotated.gtf' from records1
  file 'NCBI_transcripts.fa' from records2
  file 'transcriptome.hits' from records4

  output:
  file 'candidate_lncRNA_classes.txt' into records5
  file 'coding-genes.gtf' into records6
  file 'final_annotated.gtf' into records7
  file 'feelnc_codpot_out' into feelnc_codpot_out_dir

  shell:
  '''
  dir1=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
  cd ${dir1}
  grep "NM_" "!{params.ref_annotation}" > NM_coding.gtf
  echo ""
  printf "::: INFO: 1/3) Filtering transcripts :::\n"
  # Filter
  FEELnc_filter.pl -i final_annotated.gtf -a NM_coding.gtf -b transcript_biotype=protein_coding > candidate_lncRNA.gtf
  # rm -r -f ${g_DIR}/${reference_genome}.index
  printf "::: INFO: 2/3) Evaluating coding potential :::\n"
  # Coding_Potential
  FEELnc_codpot.pl -i candidate_lncRNA.gtf -a NM_coding.gtf -b transcript_biotype=protein_coding -g "!{params.genome}" --mode=shuffle
  printf "::: INFO: 3/3) Classifiyng lncRNA transcripts :::\n"
  # Classifier
  FEELnc_classifier.pl -i feelnc_codpot_out/candidate_lncRNA.gtf.lncRNA.gtf -a NM_coding.gtf > candidate_lncRNA_classes.txt
  echo ""
  printf "::: INFO: FEELnc calculations were done. The output is called candidate_lncRNA_classes.txt :::\n"
  echo ""
  printf "::: INFO: Parsing GAWN and FEELnc outputs :::\n"
  cd ${dir1}
  awk '{print $3}' candidate_lncRNA_classes.txt > lncRNA_genes
  tail -n +2 lncRNA_genes > lncRNA_transcripts
  rm lncRNA_genes
  grep -w -F -f lncRNA_transcripts final_annotated.gtf > merged.fixed.lncRNAs.gtf
  grep --invert-match -F -f lncRNA_transcripts final_annotated.gtf > merged.fixed.coding.gtf
  rm final_annotated.gtf
  sed -i 's/StringTie/lncRNA/' merged.fixed.lncRNAs.gtf
  awk '{print $1"\t"$2}' transcriptome.hits > coding_list
  awk -F'\t' '$2!=""' coding_list > coding_transcripts
  awk '{print $1}' coding_transcripts > coding_transcripts.tab
  rm coding_lis* coding_transcripts lncRNA_transcripts
  grep -w -F -f coding_transcripts.tab merged.fixed.coding.gtf > coding-genes.gtf
  grep --invert-match -F -f coding_transcripts.tab merged.fixed.coding.gtf > other-genes.gtf
  cat coding-genes.gtf merged.fixed.lncRNAs.gtf other-genes.gtf > final_annotated.gtf
  rm coding_transcripts.tab
  # sorting GTF file
  git clone https://github.com/cfarkas/gff3sort.git
  perl ./gff3sort/gff3sort.pl coding-genes.gtf > coding-genes.sorted.gtf
  rm coding-genes.gtf
  mv coding-genes.sorted.gtf coding-genes.gtf
  echo ""
  echo "::: INFO: Number of lines in coding-genes.gtf file:"
  cat coding-genes.gtf | wc -l
  echo ""
  printf "::: INFO: feelnc_annotation process done. Long-non-coding RNAs were succesfully extracted :::\n"
  echo ""
  printf "==== 6. We will predict coding regions from transcripts with coding potential by using TransDecoder:"
  '''
}


process TransDecoder_annotation {
  echo true
  stageInMode 'copy'
  conda "${params.conda}"

  input:
  file 'coding-genes.gtf' from records6

  output:
  file 'coding.hits' into records8
  file 'coding-transcripts.fa.transdecoder.gff3' into records9
  file 'coding-transcripts.fa.transdecoder_dir' into coding_transcripts_fa_transdecoder_dir
  file 'transdecoder' into transdecoder_dir


  shell:
  '''
  gffread -w coding-transcripts.fa -g "!{params.genome}" coding-genes.gtf
  TransDecoder.LongOrfs -m 60 -t coding-transcripts.fa
  TransDecoder.Predict -t coding-transcripts.fa --single_best_only
  awk '{print $1}' coding-transcripts.fa.transdecoder.bed > coding.sequences
  tail -n +2 coding.sequences > coding.hits && rm coding.sequences
  mkdir transdecoder
  mv coding-transcripts.fa.transdecoder.bed coding-transcripts.fa.transdecoder.cds coding-transcripts.fa.transdecoder.gff3 coding-transcripts.fa.transdecoder.pep ./transdecoder
  cp ./transdecoder/coding-transcripts.fa.transdecoder.gff3 ./
  echo ""
  printf "::: INFO: TransDecoder_annotation process done. ORFs were succesfully extracted :::\n"
  echo ""
  printf "==== 7. We will convert gff3 format to GTF and we will format coding sequences and proteins:"
  '''
}


process format_sequences_and_proteins {
  echo true
  stageInMode 'copy'
  conda "${params.conda}"

  input:
  file 'sed.script' from records0
  file 'NCBI_transcripts.fa' from records2
  file 'transcriptome.swissprot' from records3
  file 'final_annotated.gtf' from records7
  file 'coding.hits' from records8
  file 'coding-transcripts.fa.transdecoder.gff3' from records9

  output:
  file 'coding_transcripts.gtf' into records10
  file 'cds.fa' into records11
  file 'prot.fa' into records12
  file 'final_annotated.gtf' into records13
  file 'sed.script' into records14
  file 'transcriptome.swissprot' into records15

  shell:
  '''
  sed 's/Name=.*$//' coding-transcripts.fa.transdecoder.gff3 > coding-transcripts.fa.test.gff3
  sed -i 's/ID=GENE[.]/ID=/'g coding-transcripts.fa.test.gff3
  sed -i 's/Parent=GENE[.]/Parent=/'g coding-transcripts.fa.test.gff3
  sed -i 's/~~/;protein_id=/'g coding-transcripts.fa.test.gff3
  gffread coding-transcripts.fa.test.gff3 -T -P -g NCBI_transcripts.fa -o coding_transcripts.gtf
  rm coding-transcripts.fa.test.gff3
  # removing protein id by expansion
  sed -i 's/[.]p[0-9]//'g coding_transcripts.gtf
  sed -i 's/[.]p[0-9][0-9]//'g coding_transcripts.gtf
  sed -i 's/[.]p[0-9][0-9][0-9]//'g coding_transcripts.gtf
  sed -i 's/[.]p[0-9][0-9][0-9][0-9]//'g coding_transcripts.gtf
  sed -i 's/[.]p[0-9][0-9][0-9][0-9][0-9]//'g coding_transcripts.gtf
  # obtaining cds.fa and prot.fa from coding_transcripts.gtf
  echo ""
  echo "::: INFO: Obtaining cds.fa and prot.fa from coding_transcripts.gtf"
  echo ""
  gffread -x cds.fa -g NCBI_transcripts.fa coding_transcripts.gtf
  gffread -y prot.fa -g NCBI_transcripts.fa coding_transcripts.gtf
  echo "::: INFO: done"
  rm -rf coding-transcripts.fa coding-genes.gtf merged.fixed.lncRNAs.gtf other-genes.gtf
  grep "StringTie" final_annotated.gtf > genes.gtf
  grep "lncRNA" final_annotated.gtf > lncRNAs.gtf
  grep -w -F -f coding.hits genes.gtf > coding-genes.gtf
  grep --invert-match -F -f coding.hits genes.gtf > other-genes.gtf
  sed -i 's/StringTie/coding/' coding-genes.gtf
  cat coding-genes.gtf lncRNAs.gtf other-genes.gtf > final_annotated.gtf
  echo ""
  echo "::: INFO: Parsing transcriptome hits"
  echo ""
  grep -w -F -f coding.hits transcriptome.swissprot > coding.annotation
  rm transcriptome.swissprot
  mv coding.annotation transcriptome.swissprot
  echo "::: INFO: done"
  # sorting GTF file
  echo ""
  echo "::: INFO: Sorting final_annotated.gtf"
  echo ""
  git clone https://github.com/cfarkas/gff3sort.git
  perl ./gff3sort/gff3sort.pl final_annotated.gtf > final_annotated.sorted.gtf
  echo "::: INFO: done"
  rm final_annotated.gtf
  mv final_annotated.sorted.gtf final_annotated.gtf
  rm -rf coding-genes.gtf lncRNAs.gtf other-genes.gtf transcriptome.hits
  ### Novel coding genes and correspondent proteins
  echo ""
  echo "::: INFO: Now, we will obtain novel coding transcripts (cds) and correspondent proteins: "
  echo ""
  '''
}


process transcriptome_metrics {
  echo true
  stageInMode 'copy'
  conda "${params.conda}"

  input:
  file 'final_annotated.gtf' from records13
  
  output:
  file 'known-genes-coding.gtf' into records16
  file 'novel-genes-coding.gtf' into records17
  file 'novel-transcripts-lncRNA.fa' into records18
  file 'known-transcripts-lncRNA.fa' into records19
  file ''
  
  shell:
  '''
  wget https://raw.githubusercontent.com/cfarkas/annotate_my_genomes/master/additional_scripts/transcriptome_metrics.sh
  bash transcriptome_metrics.sh -f final_annotated.gtf -g "!{params.genome}"
  cp ./transcriptome_metrics/known-genes-coding.gtf ./
  cp ./transcriptome_metrics/novel-genes-coding.gtf ./
  cp ./transcriptome_metrics/novel-transcripts-lncRNA.fa ./
  cp ./transcriptome_metrics/known-transcripts-lncRNA.fa ./
  '''
}


process awk_extracted_hits {
  echo true
  stageInMode 'copy'
  conda "${params.conda}"

  input:
  file 'novel-genes-coding.gtf' from records17

  output:
  file 'novel_annotated.tab' into records20

  shell:
  '''
  awk '{print $9"\t"$10"\t"$11"\t"$12}' novel-genes-coding.gtf > novel_annotated.tab
  '''
}


process novel_hits {
  echo true
  stageInMode 'copy'
  conda "${params.conda}"

  input:
  file 'cds.fa' from records11
  file 'prot.fa' from records12
  file 'final_annotated.gtf' from records13

  output:
  file 'novel-cds.fa' into records21
  file 'novel-prot.fa' into records22
  file 'final_annotated.gff' into records23

  shell:
  '''
  awk '{print $(NF)}' novel_annotated.tab > novel-coding-transcripts.matches
  sed -i 's/;//g' novel-coding-transcripts.matches
  sed -i 's/"//g' novel-coding-transcripts.matches
  awk '!a[$0]++' novel-coding-transcripts.matches > novel-coding-transcripts.tab && rm novel-coding-transcripts.matches
  mv novel-coding-transcripts.tab novel-coding-transcripts.matches
  seqkit fx2tab cds.fa > cds.tab
  seqkit fx2tab prot.fa > prot.tab
  grep -w -F -f novel-coding-transcripts.matches cds.tab > novel-coding-cds.tab
  grep -w -F -f novel-coding-transcripts.matches prot.tab > novel-coding-prot.tab
  seqkit tab2fx novel-coding-cds.tab > novel-cds.fa && seqkit tab2fx novel-coding-prot.tab > novel-prot.fa
  rm -r -f novel-coding-cds.tab novel-coding-prot.tab novel-coding-transcripts.matches cds.tab prot.tab
  # obtaining final gff file
  echo ""
  echo "::: Obtaining final gff file"
  echo ""
  gffread -E -F --merge final_annotated.gtf -o final_annotated.gff
  rm -r -f gff3sort
  echo "done"
  echo ""
  rm -r -f merged.fixed.coding.gtf namelist namelist_unique_sorted coding.hits
  echo "------------------------------------------------------"
  echo "::: INFO: annotate_my_genomes pipeline ended correctly"
  echo "------------------------------------------------------"
  echo ""
  '''
}


process output_pipeline {
  echo true
  stageInMode 'copy'
  conda "${params.conda}"

  input:
  file 'NCBI_compare.stringtie_input_file.gtf.tmap' from tmap_records
  file 'Stats.txt' from Stats_records
  file 'NCBI_transcripts.fa' from records2
  file 'candidate_lncRNA_classes.txt' from records5
  file 'coding_transcripts.gtf' from records10
  file 'cds.fa' from records11
  file 'prot.fa' from records12
  file 'final_annotated.gtf' from records13
  file 'sed.script' from records14
  file 'transcriptome.swissprot' from records15
  file 'novel-cds.fa' from records21
  file 'novel-prot.fa' from records22
  file 'final_annotated.gff' from records23
  file 'known-genes-coding.gtf' from records16
  file 'novel-genes-coding.gtf' from records17
  file 'novel-transcripts-lncRNA.fa' from records18
  file 'known-transcripts-lncRNA.fa' from records19
  file 'gffcompare_outputs_NCBI' from gffcompare_outputs_NCBI_dir
  file 'gawn' from gawn_dir
  file 'feelnc_codpot_out' from feelnc_codpot_out_dir
  file 'coding-transcripts.fa.transdecoder_dir' from coding_transcripts_fa_transdecoder_dir
  file 'transdecoder' from transdecoder_dir

  shell:
  '''
  echo ""
  printf "::: Moving results to the output directory :::\n"
  rm -r -f output_files
  mkdir output_files
  mv NCBI_compare.stringtie_input_file.gtf.tmap gffcompare.tmap
  mv gffcompare.tmap candidate_lncRNA_classes.txt final_annotated.gtf final_annotated.gff NCBI_transcripts.fa cds.fa prot.fa coding_transcripts.gtf Stats.txt transcriptome.swissprot novel-cds.fa novel-prot.fa sed.script novel-transcripts-lncRNA.fa known-transcripts-lncRNA.fa known-genes-coding.gtf novel-genes-coding.gtf ./output_files

  if [ -z "$(ls -A "!{params.outdir}")" ]; then
     echo ""
     echo "Output folder is empty. We will work inside the provided output folder: "
     mv output_files "!{params.outdir}"
     mv gffcompare_outputs_NCBI "!{params.outdir}"
     mv gawn "!{params.outdir}"
     mv feelnc_codpot_out "!{params.outdir}"
     mv coding-transcripts.fa.transdecoder_dir "!{params.outdir}"
     mv transdecoder "!{params.outdir}"
     echo ""
  else
     echo ""
     echo "Output folder is not empty. Creating temporary folder:"
     sec=$(date "+%Y%m%d_%H%M%S")
     mkdir add_ncbi_annotation_$sec
     mv output_files ./add_ncbi_annotation_$sec
     mv gffcompare_outputs_NCBI ./add_ncbi_annotation_$sec
     mv gawn ./add_ncbi_annotation_$sec
     mv feelnc_codpot_out ./add_ncbi_annotation_$sec
     mv coding-transcripts.fa.transdecoder_dir ./add_ncbi_annotation_$sec
     mv transdecoder ./add_ncbi_annotation_$sec
     mv add_ncbi_annotation_$sec "!{params.outdir}"
  fi

  echo ""
  echo "------------------------------------------------------------"
  echo "------------------------------------------------------------"
  echo "::: INFO: all done" :::
  echo ""
  echo "The following files are available in the output directory : "
  echo ""
  echo "Transcript discoveries are summarized in Stats.txt file. GAWN protein annotation is named transcriptome.hits"
  echo ""
  echo "gffcompare.tmap file contains Best Reference Transcript for each assembled transcript"
  echo ""
  echo "GTF file named final_annotated.gtf (and correspondent gff file) contain novel genes and lncRNA classification (second field in GTF file)"
  echo ""
  echo "candidate_lncRNA_classes.txt contained detailed long non-coding classification of transcripts".
  echo ""
  echo "Associated FASTA file to this GTF correspond to NCBI_transcripts.fa file"
  echo ""
  echo "TransDecoder GTF file suitable to parse transcripts.fa (coding_transcripts.gtf), contains all coding transcripts resolved by TransDecoder"
  echo ""
  echo "Predicted coding sequences and correspondent protein sequences were named cds.fa and prot.fa, respectively"
  echo ""
  echo "Novel predicted coding sequences and correspondent protein sequences were named novel-cds.fa and novel-prot.fa, respectively"
  echo ""
  echo "Novel and Known predicted lncRNAs were named novel-transcripts-lncRNA.fa and known-transcripts-lncRNA.fa, respectively"
  echo ""
  echo "Novel and Known coding genes were named novel-genes-coding.gtf and known-genes-coding.gtf, respectively"
  echo ""
  echo "------------------------------------------------------------"
  echo "------------------------------------------------------------"
  echo ""
  '''
}
