#!/bin/bash

set -e

{

usage="$(basename "$0") [-h] [-a <stringtie.gtf>] [-n <NCBI_reference.gtf>] [-r <reference_genome.gtf>] [-g <reference_genome.fasta>] [-c <gawn_config>] [-t <threads>] [-o <output>]
This pipeline will Overlap StringTie transcripts (GTF format) with current NCBI annotation and will annotate novel transcripts.
Arguments:
    -h  show this help text
    -a  StringTie GTF
    -n  NCBI gene annotation (in GTF format)
    -r  UCSC gene annotation (in GTF format)
    -g  Reference genome (in fasta format)
    -c  GAWN config file (path to gawn_config.sh in annotate_my_genomes folder)
    -t  Number of threads for processing (integer)
    -o  output folder (must exist)"
options=':ha:n:r:g:c:t:o:'
while getopts $options option; do
  case "$option" in
    h) echo "$usage"; exit;;
    a) a=$OPTARG;;
    n) n=$OPTARG;;
    r) r=$OPTARG;;
    g) g=$OPTARG;;
    c) c=$OPTARG;;
    t) t=$OPTARG;;
    o) o=$OPTARG;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
  esac
done

# mandatory arguments
if [ ! "$a" ] || [ ! "$n" ] || [ ! "$r" ] || [ ! "$g" ] || [ ! "$c" ] || [ ! "$t" ] || [ ! "$o" ]; then
  echo "arguments -a, -n, -r, -g, -c, -t and -o must be provided"
  echo "$usage" >&2; exit 1
fi

# Conditions : output folder
if [ ! -d "$o" ]; then
  echo "Output directory: $o not found. Please create the output directory first, before running the pipeline."
  exit 9999 # die with error code 9999
fi

# Conditions : Input existance if [ ! "$a" ] || [ ! "$n" ] || [ ! "$r" ] || [ ! "$g" ] || [ ! "$c" ] || [ ! "$t" ] || [ ! "$o" ]; then

if [ ! -e "$a" ]; then
    echo "$a does not exist. Check your -a input"
    exit 9999 # die with error code 9999
fi

if [ ! -e "$n" ]; then
    echo "$n does not exist. Check your -n input"
    exit 9999 # die with error code 9999
fi

if [ ! -e "$r" ]; then
    echo "$r does not exist. Check your -r input"
    exit 9999 # die with error code 9999
fi

if [ ! -e "$g" ]; then
    echo "$g does not exist. Check your -g input"
    exit 9999 # die with error code 9999
fi

if [ ! -e "$c" ]; then
    echo "$c does not exist. Check your -c input"
    exit 9999 # die with error code 9999
fi

# Conditions : Getting absolute path of inputs
echo ""
a_DIR="$( cd "$( dirname "$a" )" && pwd )"
echo ""
echo "::: The absolute path of -a is $a_DIR"
echo ""
n_DIR="$( cd "$( dirname "$n" )" && pwd )"
echo ""
echo "::: The absolute path of -n is $n_DIR"
echo ""
r_DIR="$( cd "$( dirname "$r" )" && pwd )"
echo ""
echo "::: The absolute path of -r is $r_DIR"
echo ""
g_DIR="$( cd "$( dirname "$g" )" && pwd )"
echo ""
echo "::: The absolute path of -g is $g_DIR"
echo ""
c_DIR="$( cd "$( dirname "$c" )" && pwd )"
echo ""
echo "::: The absolute path of -c is $c_DIR"
echo ""
o_DIR="$( cd "$( dirname "$o" )" && pwd )"
echo ""
echo "::: The absolute path of -o is $o_DIR"
echo ""


begin=`date +%s`
#    .---------- constant part!
#    vvvv vvvv-- the code from above
YELLOW='\033[1;33m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color


printf "${YELLOW}::: Defining Variables :::\n"
echo ""
echo "Defining variables:"
echo""
FILE1="$a"
basename "$FILE1"
stringtie_input="$(basename -- $FILE1)"
echo "The stringtie file used as input is the following: $stringtie_input"
echo ""
FILE2="$n"
basename "$FILE2"
ncbi_reference_gtf="$(basename -- $FILE2)"
echo "The NCBI reference GTF used as input is the following: $ncbi_reference_gtf"
echo ""
FILE3="$r"
basename "$FILE3"
reference_gtf="$(basename -- $FILE3)"
echo "The reference GTF used as input is the following: $reference_gtf"
echo ""
FILE4="$g"
basename "$FILE4"
reference_genome="$(basename -- $FILE4)"
echo "The reference genome used as input is the following: $reference_genome"
echo ""
FILE5="$c"
basename "$FILE5"
gawn_config="$(basename -- $FILE5)"
echo "The gawn_config file used as input is the following: $gawn_config"
echo ""
FILE6="$t"
basename "$FILE6"
threads="$(basename -- $FILE6)"
echo "The number of threads for calculation are the following: $threads"
echo ""
FILE7="$o"
basename "$FILE7"
output_folder="$(basename -- $FILE7)"
echo "The output folder is the following: $output_folder"
echo ""

printf "${YELLOW}:::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 0. Defining directories :::\n"
printf "${YELLOW}:::::::::::::::::::::::::::::::${CYAN}\n"
echo ""

dir0=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

sec=$(date "+%Y%m%d_%H%M%S")
# mkdir add_ncbi_annotation_$sec

if [ -z "$(ls -A ${o_DIR}/${output_folder})" ]; then
   echo ""
   echo "Output folder is empty. We will work inside the provided output folder: "
   cd ${o_DIR}/${output_folder}
   echo ""
else
   echo ""
   echo "Output folder is not empty. Creating temporary folder:"
   sec=$(date "+%Y%m%d_%H%M%S")
   cd ${o_DIR}/${output_folder}
   mkdir add_ncbi_annotation_$sec && cd add_ncbi_annotation_$sec
fi

# cd ${o_DIR}/${output_folder}

if [ -f $stringtie_input ]; then
    echo ""
    echo "$stringtie_input file found in output directory. Continue."
    echo ""
    : 
else
    echo ""
    echo "Copying $stringtie_input file into the output directory:"
    cp ${a_DIR}/${stringtie_input} ./
    echo ""
fi

# cp ${a_DIR}/${stringtie_input} ${o_DIR}/${output_folder}
# cd add_ncbi_annotation_$sec

dir1=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
echo ""
echo "Current Working Directory:"
echo ""
echo $dir1
echo ""
printf "${YELLOW}::: Done :::\n"
echo ""

printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 1. Overlapping StringTie transcripts with NCBI annotation :::\n"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
echo ""

gffcompare -R -r ${n_DIR}/${ncbi_reference_gtf} -s ${g_DIR}/${reference_genome} -o NCBI_compare ${stringtie_input}
echo "Done."
printf "${PURPLE}::: Done :::\n"
echo ""

printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 2. Writting novel discoveries to Stats.txt :::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"

echo ""
# Stats
exec 3<> Stats.txt
echo "Number of assembled genes:" >> Stats.txt
cat NCBI_compare.${stringtie_input}.tmap | sed "1d" | cut -f4 | sort | uniq | wc -l >> Stats.txt
echo "" >> Stats.txt
echo "Number of novel genes:" >> Stats.txt
cat NCBI_compare.${stringtie_input}.tmap | awk '$3=="u"{print $0}' | cut -f4 | sort | uniq | wc -l >> Stats.txt
echo "" >> Stats.txt
echo "Number of novel transcripts:" >> Stats.txt
cat NCBI_compare.${stringtie_input}.tmap | awk '$3=="u"{print $0}' | cut -f5 | sort | uniq | wc -l >> Stats.txt
echo "" >> Stats.txt
echo "Number of transcripts matching annotation:" >> Stats.txt
cat NCBI_compare.${stringtie_input}.tmap | awk '$3=="="{print $0}' | cut -f5 | sort | uniq | wc -l >> Stats.txt
exec 3>&-
printf "${PURPLE}Done\n"
echo ""

printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 3. Replacing gene_id field in final_annotated.gtf file with NCBI gene_id's :::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"

echo ""
#######################################
# Merging novel transcripts with ref. #
#######################################
awk '{print $4"\t"$1}' NCBI_compare.${stringtie_input}.tmap > NCBI_compare.${stringtie_input}.tmap.1
tail -n +2 NCBI_compare.${stringtie_input}.tmap.1 > NCBI_compare.${stringtie_input}.tmap.2
awk '$2 != "-"' NCBI_compare.${stringtie_input}.tmap.2 > namelist
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
cat ${a_DIR}/${stringtie_input} | parallel --pipe -j ${t} sed -f sed.script > final_annotated.gtf
rm -f fileA fileB *tmap.1 *tmap.2
# sorting GTF file
rm -r -f gff3sort
git clone https://github.com/cfarkas/gff3sort.git
perl ./gff3sort/gff3sort.pl final_annotated.gtf > final_annotated.sorted.gtf
rm final_annotated.gtf
mv final_annotated.sorted.gtf final_annotated.gtf
printf "${PURPLE}::: Done. Gene_id field was replaced in the StringTie.gtf file and final_annotated.gtf was generated with these changes\n"
echo ""
printf "${PURPLE}::: Moving gffcompare results to gffcompare_outputs folder ...\n"
echo ""
rm -r -f gffcompare_outputs_NCBI
mkdir gffcompare_outputs_NCBI
mv *.loci *.stats *.refmap *.tmap *.tracking ./gffcompare_outputs_NCBI
echo ""
printf "${PURPLE}::: Done\n"
echo ""

printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 4. Obtaining Transcripts in FASTA format with gffread :::\n"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"

echo ""
gffread -w NCBI_transcripts.fa -g ${g_DIR}/${reference_genome} final_annotated.gtf
echo ""
printf "${PURPLE}::: Done. NCBI_transcripts.fa are located in current directory\n"
echo ""

printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 5. Performing gene annotation by using GAWN pipeline :::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"

################################################################
# Configuring Gawn Inputs, config file and running GAWN pipeline
################################################################

echo ""
printf "${PURPLE}::: Downloading GAWN annotation folder. See https://github.com/enormandeau/gawn.git${CYAN}\n"
echo ""
rm -r -f gawn
git clone https://github.com/cfarkas/gawn.git
cd gawn/02_infos/
dir2=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
echo "Done"
echo ""
cd ${dir1}
cp ${g_DIR}/${reference_genome} ${dir1}/gawn/03_data/genome.fasta
cp NCBI_transcripts.fa ${dir1}/gawn/03_data/transcriptome.fasta
rm ${dir2}/gawn_config.sh
cp ${c_DIR}/${gawn_config} ${dir2}/gawn_config.sh
echo ""
printf "${PURPLE}::: Starting GAWN transcript annotation${CYAN}\n"
echo ""
cd ${dir1}/gawn/
./gawn 02_infos/gawn_config.sh
echo ""
printf "${PURPLE}::: Done. The novel transcripts were annotated in ./gawn/04_annotation/ :::${CYAN}\n"
echo ""

#################################
# Extracting transcriptome hits #
#################################

printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 6. Extracting transcriptome hits :::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::${CYAN}\n"

echo ""
cd ${dir1}/
cp ${dir1}/gawn/04_annotation/transcriptome.swissprot ${dir1}
cp ${dir1}/gawn/04_annotation/transcriptome.hits ${dir1}
printf "${PURPLE}::: Done. transcriptome hits were succesfully extracted :::${CYAN}\n"
echo ""

############################################
# FEELnc long noncoding RNA identification #
############################################

cd ${dir1}

printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 7. Classifying protein-coding and long non-coding transcripts with FEELnc :::\n"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"

grep "NM_" ${r_DIR}/${reference_gtf} > NM_coding.gtf
echo ""
printf "${PURPLE}::: 1/3) Filtering transcripts :::${CYAN}\n"
# Filter
FEELnc_filter.pl -i final_annotated.gtf -a NM_coding.gtf -b transcript_biotype=protein_coding > candidate_lncRNA.gtf
rm -r -f ${g_DIR}/${reference_genome}.index
printf "${PURPLE}::: 2/3) Evaluating coding potential :::${CYAN}\n"
# Coding_Potential
FEELnc_codpot.pl -i candidate_lncRNA.gtf -a NM_coding.gtf -b transcript_biotype=protein_coding -g ${g_DIR}/${reference_genome} --mode=shuffle
printf "${PURPLE}::: 3/3) Classifiyng lncRNA transcripts :::${CYAN}\n"
# Classifier
FEELnc_classifier.pl -i feelnc_codpot_out/candidate_lncRNA.gtf.lncRNA.gtf -a NM_coding.gtf > candidate_lncRNA_classes.txt
echo ""
printf "${PURPLE}::: FEELnc calculations were done. The output is called candidate_lncRNA_classes.txt :::\n"
echo ""

printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 8. Parsing GAWN and FEELnc outputs :::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"

echo ""
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
perl ./gff3sort/gff3sort.pl coding-genes.gtf > coding-genes.sorted.gtf
rm coding-genes.gtf
mv coding-genes.sorted.gtf coding-genes.gtf
echo "All done"
echo ""
##########################################
# Gene Prediction Step with TransDecoder #
##########################################
cd ${dir1}
echo ""

printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 9. Predicting coding regions from transcripts with coding potential using TransDecoder :::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"

echo ""
gffread -w coding-transcripts.fa -g ${g_DIR}/${reference_genome} coding-genes.gtf
TransDecoder.LongOrfs -m 60 -t coding-transcripts.fa
TransDecoder.Predict -t coding-transcripts.fa --single_best_only
awk '{print $1}' coding-transcripts.fa.transdecoder.bed > coding.sequences
tail -n +2 coding.sequences > coding.hits && rm coding.sequences
echo ""
printf "${PURPLE}::: Done. coding-transcripts.fa.transdecoder.gff3 file is present in current directory...${CYAN}\n"
echo ""

printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 10. Converting gff3 to GTF format and formatting coding sequences and proteins :::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"

echo ""
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
#
# obtaining cds.fa and prot.fa from coding_transcripts.gtf
echo ""
echo "::: Obtaining cds.fa and prot.fa from coding_transcripts.gtf"
echo ""
gffread -x cds.fa -g NCBI_transcripts.fa coding_transcripts.gtf
gffread -y prot.fa -g NCBI_transcripts.fa coding_transcripts.gtf
echo "done"
rm coding-transcripts.fa coding-genes.gtf merged.fixed.lncRNAs.gtf other-genes.gtf
grep "StringTie" final_annotated.gtf > genes.gtf
grep "lncRNA" final_annotated.gtf > lncRNAs.gtf
grep -w -F -f coding.hits genes.gtf > coding-genes.gtf
grep --invert-match -F -f coding.hits genes.gtf > other-genes.gtf
sed -i 's/StringTie/coding/' coding-genes.gtf
cat coding-genes.gtf lncRNAs.gtf other-genes.gtf > final_annotated.gtf
echo ""
echo "::: Parsing transcriptome hits"
echo ""
grep -w -F -f coding.hits transcriptome.swissprot > coding.annotation
rm transcriptome.swissprot
mv coding.annotation transcriptome.swissprot
echo "done"
# sorting GTF file
echo ""
echo "::: Sorting final_annotated.gtf"
echo ""
perl ./gff3sort/gff3sort.pl final_annotated.gtf > final_annotated.sorted.gtf
echo "done"
rm final_annotated.gtf
mv final_annotated.sorted.gtf final_annotated.gtf
rm coding-genes.gtf lncRNAs.gtf other-genes.gtf transcriptome.hits
### Novel coding genes and correspondent proteins
echo ""
echo "::: Obtaining novel coding transcripts (cds) and correspondent proteins"
echo ""
#
wget https://raw.githubusercontent.com/cfarkas/annotate_my_genomes/master/additional_scripts/transcriptome_metrics.sh
bash transcriptome_metrics.sh -f final_annotated.gtf -g ${g_DIR}/${reference_genome}
cp ./transcriptome_metrics/known-genes-coding.gtf ./
cp ./transcriptome_metrics/novel-genes-coding.gtf ./
cp ./transcriptome_metrics/novel-transcripts-lncRNA.fa ./
cp ./transcriptome_metrics/known-transcripts-lncRNA.fa ./
#
perl -lne 'print "@m" if @m=(/((?:transcript_id|gene_id)\s+\S+)/g);' novel-genes-coding.gtf > novel_annotated.tab
awk '{print $(NF)}' novel_annotated.tab > novel-coding-transcripts.matches
sed -i 's/;//g' novel-coding-transcripts.matches
sed -i 's/"//g' novel-coding-transcripts.matches
awk '!a[$0]++' novel-coding-transcripts.matches > novel-coding-transcripts.tab && rm novel-coding-transcripts.matches
mv novel-coding-transcripts.tab novel-coding-transcripts.matches
#
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
rm merged.fixed.coding.gtf namelist namelist_unique_sorted coding.hits

###############################
# Configuring Summary Results #
###############################

printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 11. Moving results to the specified directory :::\n"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"

echo ""
printf "${PURPLE}::: Moving results to the specified directory :::${CYAN}\n"
rm -r -f output_files
mkdir output_files
mv candidate_lncRNA_classes.txt final_annotated.gtf final_annotated.gff NCBI_transcripts.fa cds.fa prot.fa Stats.txt coding_transcripts.gtf transcriptome.swissprot novel-cds.fa novel-prot.fa sed.script novel-transcripts-lncRNA.fa known-transcripts-lncRNA.fa known-genes-coding.gtf novel-genes-coding.gtf ./output_files
rm -r -f *feelncfilter.log genes.gtf pipeliner* NM_coding.gtf candidate_lncRNA.gtf* coding-transcripts.fa.transdecoder_dir.__* NCBI_transcripts.fa.fai
rm -r -f transdecoder
mkdir transdecoder
mv coding-transcripts.fa.transdecoder.* ./transdecoder
mv NCBI_compare.annotated.gtf ./gffcompare_outputs_NCBI
cp ${dir1}/gffcompare_outputs_NCBI/NCBI_compare.${stringtie_input}.tmap  ./
mv NCBI_compare.${stringtie_input}.tmap gffcompare.tmap
mv gffcompare.tmap ./output_files/

# cd ${dir0}
# mv add_ncbi_annotation_$sec ${o_DIR}/${output_folder}

cd ${dir0}

echo "Done"
echo ""
printf "${YELLOW}::: Done:::\n"
echo ""
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
echo ""
echo "The following files are available in ${dir1}/output_files : "
echo ""
echo "Transcript discoveries are summarized in Stats.txt file. GAWN protein annotation is called transcriptome.hits"
echo ""
echo "gffcompare.tmap file contains Best Reference Transcript for each assembled transcript"
echo ""
echo "GTF file named final_annotated.gtf (and correspondent gff file) contain novel genes and lncRNA classification (second field in GTF file)"
echo ""
echo "candidate_lncRNA_classes.txt contained detailed long non-coding classification of transcripts"
echo ""
echo "Associated FASTA file to this GTF correspond to NCBI_transcripts.fa file"
echo ""
echo "TransDecoder GTF file suitable to parse NCBI_transcripts.fa (coding_transcripts.gtf), contains all coding transcripts resolved by TransDecoder"
echo ""
echo "Predicted coding sequences and correspondent protein sequences were named cds.fa and prot.fa, respectively"
echo ""
echo "Novel predicted coding sequences and correspondent protein sequences were named novel-cds.fa and novel-prot.fa, respectively"
echo ""
echo "Novel and Known predicted lncRNAs were named novel-transcripts-lncRNA.fa and known-transcripts-lncRNA.fa, respectively"
echo ""
echo "Novel and Known coding genes were named novel-genes-coding.gtf and known-genes-coding.gtf, respectively"
echo ""
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${NC}\n"
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed
#
} | tee logfile_add_ncbi_annotation_$seconds
#
