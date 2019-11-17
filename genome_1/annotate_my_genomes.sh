#!/bin/bash
{

dir1=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
stringtie_gtf=${1}
reference_genome_name=${2}
threads=${3}

if [ "$1" == "-h" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {stringtie_gtf} {reference_genome_name} {threads}"
  echo ""
  echo "This script will Overlap StringTie transcripts (GTF format) with UCSC reference genome GTF and annotate novel transcripts"
  echo ""
  echo "{stringtie_gtf}: Name of the StringTie GTF file"
  echo ""
  echo "{reference_genome_name}: Name of the current assembly genome (use UCSC genome names)"
  echo ""
  echo "{threads}: Number of threads for parallel text processing (Integer)"
  echo ""
  exit 0
fi

if [ "$1" == "-help" ]; then
    echo ""
  echo "Usage: bash ./`basename $0` {stringtie_gtf} {reference_genome_name} {threads}"
  echo ""
  echo "This script will Overlap StringTie transcripts (GTF format) with UCSC reference genome GTF and annotate novel transcripts"
  echo ""
  echo "{stringtie_gtf}: Name of the StringTie GTF file"
  echo ""
  echo "{reference_genome_name}: Name of the current assembly genome (use UCSC genome names)"
  echo ""
  echo "{threads}: Number of threads for parallel text processing (Integer)"
  echo ""
  exit 0
fi
if [ "$1" == "--h" ]; then
   echo ""
  echo "Usage: bash ./`basename $0` {stringtie_gtf} {reference_genome_name} {threads}"
  echo ""
  echo "This script will Overlap StringTie transcripts (GTF format) with UCSC reference genome GTF and annotate novel transcripts"
  echo ""
  echo "{stringtie_gtf}: Name of the StringTie GTF file"
  echo ""
  echo "{reference_genome_name}: Name of the current assembly genome (use UCSC genome names)"
  echo ""
  echo "{threads}: Number of threads for parallel text processing (Integer)"
  echo ""
  exit 0
fi
if [ "$1" == "--help" ]; then
    echo ""
  echo "Usage: bash ./`basename $0` {stringtie_gtf} {reference_genome_name} {threads}"
  echo ""
  echo "This script will Overlap StringTie transcripts (GTF format) with UCSC reference genome GTF and annotate novel transcripts"
  echo ""
  echo "{stringtie_gtf}: Name of the StringTie GTF file"
  echo ""
  echo "{reference_genome_name}: Name of the current assembly genome (use UCSC genome names)"
  echo ""
  echo "{threads}: Number of threads for parallel text processing (Integer)"
  echo ""
  exit 0
fi

[ $# -eq 0 ] && { echo "Usage: bash ./`basename $0` {stringtie_gtf} {reference_genome_name} {threads}"; exit 1; }

if [ $# -ne 3 ]; then
  echo 1>&2 "Usage: bash ./`basename $0` {stringtie_gtf} {reference_genome_name} {threads}"
  exit 3
fi

begin=`date +%s`
echo ""

echo "::: Downloading Reference genome and current GTF annotation"
echo ""
bash genome_download.sh ${2}
echo ""
echo "Done"
echo ""

echo "::: Overlapping StringTie transcripts with Reference"
echo ""
gffcompare -R -r ${2}.gtf -o merged ${1}
echo ""
echo "Done"

########################################
# Stats
########################################
exec 3<> Stats.txt
echo "Number of assembled genes:" >> Stats.txt
cat merged.${1}.tmap | sed "1d" | cut -f4 | sort | uniq | wc -l >> Stats.txt
echo "" >> Stats.txt
echo "Number of novel genes:" >> Stats.txt
cat merged.${1}.tmap | awk '$3=="u"{print $0}' | cut -f4 | sort | uniq | wc -l >> Stats.txt
echo "" >> Stats.txt
echo "Number of novel transcripts:" >> Stats.txt
cat merged.${1}.tmap | awk '$3=="u"{print $0}' | cut -f5 | sort | uniq | wc -l >> Stats.txt
echo "" >> Stats.txt
echo "Number of transcripts matching annotation:" >> Stats.txt
cat merged.${1}.tmap | awk '$3=="="{print $0}' | cut -f5 | sort | uniq | wc -l >> Stats.txt
exec 3>&-
echo ""

echo "::: Writting novel discoveries to Stats.txt"
echo ""
echo "Done"
echo ""

echo "::: Replacing gene_id field in merged.annotated.gtf file with reference gene_id's"
echo ""

########################################
# Merging novel transcripts with ref. 
########################################
awk '{print $4"\t"$1}' merged.${1}.tmap > merged.${1}.tmap.1
tail -n +2 merged.${1}.tmap.1 > merged.${1}.tmap.2
awk '!/-/' merged.${1}.tmap.2 > namelist
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
rm A A.1 A.2 B B.1 B.2 refGene.txt namelist_unique_sorted

#######################################
# Getting gene names replaced in a loop
#######################################
awk '{print $1}' namelist > fileA
awk '{print $2}' namelist > fileB
tmpfile=/tmp/Asasuser.$$
exec 3< fileA
exec 4< fileB
while read -r astring <&3
do
        read -r bstring <&4
        echo "s/"$astring"/"$bstring"/" >> "$tmpfile"
done
exec 3<&- 4<&-
        cat merged.annotated.gtf | parallel --pipe -j ${3} sed -f "$tmpfile" > merged_with_reference.gtf
rm -f "$tmpfile" fileA fileB
echo "Done"
echo ""
echo "A GTF is called merged_with_reference.gtf with fixed gene names is located in the current directory."
echo ""
echo "::: Obtaining Transcripts in FASTA format with gffread"
echo ""
gffread -w transcripts.fa -g ${2}.fa merged_with_reference.gtf
echo ""
echo "Moving gffcompare results to gffcompare_outputs folder ..."
mkdir gffcompare_outputs
mv *.annotated.gtf *.loci *.stats *.refmap *.tmap *.tracking ./gffcompare_outputs
echo ""
echo "Done. Continue with GAWN annotation..."

################################################################
# Configuring Gawn Inputs, config file and running GAWN pipeline
################################################################
echo ""
echo "Downloading GAWN annotation folder. See https://github.com/enormandeau/gawn.git"
echo ""
git clone https://github.com/enormandeau/gawn.git
cd gawn/02_infos/
dir2=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
echo "Done"
echo ""
cd /${dir1}/
cp ${2}.fa /${dir1}/gawn/03_data/genome.fasta
cp transcripts.fa /${dir1}/gawn/03_data/transcriptome.fasta
rm /${dir2}/gawn_config.sh
cp gawn_config.sh /${dir2}/gawn_config.sh

echo ""

echo "::: Starting GAWN transcript annotation"
echo ""

cd /${dir1}/gawn/
./gawn 02_infos/gawn_config.sh

echo ""
echo "Done. The novel transcripts are annotated in ./gawn/05_results/"
echo ""

##########################################
# CNIT long noncoding identification
##########################################
echo "::: Classifying protein-coding and long non-coding transcripts with CNIT"
echo ""
cd /${dir1}/
cd ..
cd CNIT/
python CNIT.py -f /${dir1}/transcripts.fa -o /${dir1}/CNIT_output -m 've'
cd /${dir1}/
echo ""
echo "Done. The transcripts were classified in ./CNIT_output/"
echo ""
echo "::: Configuring Summary results"
cp /${dir1}/gawn/05_results/transcriptome_annotation_table.tsv /${dir1}/
cp /${dir1}/CNIT_output/CNIT.index /${dir1}/
awk '{print $1"\t"$2}' transcriptome_annotation_table.tsv > transcriptome_annotation
awk '{print $2"\t"$3"\t"$4"\t"$7}' CNIT.index > CNIT
tail -n +2 transcriptome_annotation > transcriptome_annotation_file
tail -n +2 CNIT > CNIT_file

##########################################
# Extracting GO terms for each transcript
##########################################
cut -d$'\t' -f 1,6 transcriptome_annotation_table.tsv > transcripts_GO
tr ';' '\t' < transcripts_GO > transcripts_GO_sep
column -t transcripts_GO_sep > transcripts_GO.tab
tail -n +2 transcripts_GO.tab > transcriptsGO.tab
rm transcripts_GO*

##########################################
# Gene Prediction Step with Augustus
##########################################
cd /${dir1}/
echo ""
echo "::: Predicting gene models from transcripts with AUGUSTUS (gff3 format). Progress will be printed for each transcript."
echo ""
wget http://augustus.gobics.de/binaries/augustus.2.5.5.tar.gz
gunzip augustus.2.5.5.tar.gz
tar -xvf augustus.2.5.5.tar
cd augustus.2.5.5/src/
make
cd ..
cd ..
export AUGUSTUS_CONFIG_PATH=./augustus.2.5.5/config/
./augustus.2.5.5/src/augustus --species=human --progress=true --UTR=on --uniqueGeneId=true --gff3=on transcripts.fa > augustus.gff3
echo "Done. augustus.gff3 file is present in current directory..."
echo ""
echo "Converting gff3 to GTF format, collecting coding sequences and proteins with gffread..."
gffread augustus.gff3 -T -o coding_transcript.gtf
gffread -x cds.fa -g transcripts.fa coding_transcript.gtf
gffread -y prot.fa -g transcripts.fa coding_transcript.gtf

# Re-formatting
cat cds.fa |rev|cut -d"." -f1 --complement|rev > transcripts_CDS.fa
cat prot.fa |rev|cut -d"." -f1 --complement|rev > transcripts_proteins.fa
rm cds.fa prot.fa

##########################################
# Re-formatting coding_transcripts.gtf
##########################################
sed 's/.t1"/"/' coding_transcript.gtf > coding_transcripts.gtf
echo ""
echo "Done. AUGUSTUS predicted transcripts were summarized in coding_transcripts.gtf file located in current directory"
rm A B coding_transcript.gtf 
##########################################
# Clean-up enviroment and building summary
##########################################
rm transcriptome_annotation_table.tsv transcriptome_annotation CNIT.index CNIT
join --nocheck-order transcriptome_annotation_file CNIT_file > file1
awk '!a[$0]++' genes1.tab > genes.tab
sed -e 's/^/s%/' -e 's/\t/%/' -e 's/$/%g/' genes.tab |
sed -f - file1 > summary.txt
sed 's/ /\t/g' summary.txt > summary.tab
rm transcriptome_annotation_file CNIT_file file1

perl -lne 'print "@m" if @m=(/((?:transcript_id|gene_id)\s+\S+)/g);' merged_with_reference.gtf > transcript_gene_names.txt
sed -i 's/transcript_id //g' transcript_gene_names.txt
sed -i 's/;/\t/g' transcript_gene_names.txt
sed -i 's/gene_id//g' transcript_gene_names.txt
sed -i 's/"//g' transcript_gene_names.txt
awk '{print $1"\t"$2}' transcript_gene_names.txt > transcript_gene_names.tab
join --nocheck-order transcript_gene_names.tab summary.tab > file1
awk '!a[$0]++' file1 > final_summary.tab

##########################################
# Extracting GO terms for each Gene
##########################################
join --nocheck-order final_summary.tab transcriptsGO.tab > summary_with_GO.tab
awk '{$1=$3=$4=$5=$6=""; print $0}' summary_with_GO.tab > genesGO.tab
sed -i 's/^ //g' genesGO.tab
awk '!a[$0]++' genesGO.tab > genesGO.tab1
mv genesGO.tab1 genesGO.tab
rm summary_with_GO.tab summary.tab merged.annotated.gtf

############################################
# Moving results to merged_annotation folder
############################################
echo ""
echo "Moving results to merged_annotation folder"
mkdir merged_annotation
mv merged_with_reference.gtf final_summary.tab Stats.txt transcripts.fa transcriptsGO.tab genesGO.tab transcripts_CDS.fa transcripts_proteins.fa coding_transcripts.gtf logfile ./merged_annotation
cp /${dir1}/gawn/05_results/transcriptome_annotation_table.tsv /${dir1}/merged_annotation/
rm summary.txt genes.tab genes1.tab genes2.tab transcripts.fa.fai namelist
echo ""
echo "All Done. The transcripts were classified in final_summary.tab file located in ./merged_annotation"
echo ""
echo "The above Classification includes:"
echo ""
echo "Field #1: StringTie transcript id (STRG)"
echo "Field #2: StringTie gene id (STRG) or gene symbol from UCSC (if transcript matched with a reference gene)"
echo "Field #3: Uniprot closest BLASTp match from GAWN annotation"
echo "Field #4: Coding or non coding classfication by CNIT"
echo "Field #5: Coding probability score (CNIT Score)"
echo "Field #6: Transcript length (bp)"
echo ""
echo ""
echo "Transcript discoveries are summarized in Stats.txt file located in ./merged_annotation . GAWN annotation is named transcriptome_annotation_table.tsv"
echo ""
echo "A new GTF file suitable for gene count quantification is named merged_with_reference.gtf and is located in ./merged_annotation"
echo ""
echo "Associated FASTA file to this GTF, named transcripts.fa is located in ./merged_annotation"
echo ""
echo "AUGUSTUS GTF file suitable for transcript count quantification is named coding_transcripts.gtf. This GTF file contains all coding transcripts resolved by AUGUSTUS and is located in ./merged_annotation"
echo ""
echo "Associated Transcript coding sequences (transcripts_CDS.fa) and correspondent protein sequences (transcripts_proteins.fa) with coding_transcripts.gtf are located in ./merged_annotation"
echo ""
echo "GO terms associated to each transcript (and gene), named transcriptsGO.tab and genesGO.tab are located in ./merged_annotation"
echo ""
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed
#
} | tee logfile
#
