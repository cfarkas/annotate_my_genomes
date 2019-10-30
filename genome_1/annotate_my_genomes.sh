# annotate_my_genomes 10/10/19
#!/bin/bash

# Preeliminars 

## Obtaining StringTie
#git clone https://github.com/gpertea/stringtie
#cd stringtie
#make release
#sudo cp stringtie /usr/local/bin/

## Obtaining gffcompare and gffread
#git clone https://github.com/gpertea/gclib
#git clone https://github.com/gpertea/gffcompare
#git clone https://github.com/gpertea/gffread
#cd gffcompare
#make release
#sudo cp gffcompare trmap /usr/local/bin/
#cd ..
#cd gffread
#make release
#sudo cp gffread /usr/local/bin/
#cd ..

## Installing up-to-date ncbi-blast+ version
# sudo apt-get remove ncbi-blast+
# wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.9.0+-x64-linux.tar.gz
# tar -xzvf ncbi-blast-2.9.0+-x64-linux.tar.gz
# rm ncbi-blast-2.9.0+-x64-linux.tar.gz
# cd ncbi-blast-2.9.0+/bin/
# sudo cp * /usr/local/bin/ 

dir1=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
stringtie_gtf=${1}
reference_genome_name=${2}

if [ "$1" == "-h" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {stringtie_gtf} {reference_genome_name}"
  echo ""
  echo "This script will Overlap StringTie transcripts (GTF format) with UCSC reference genome GTF and annotate novel transcripts with GAWN"
  echo ""
  echo "{stringtie_gtf}: Name of the StringTie GTF file"
  echo ""
  echo "{reference_genome_name}: Name of the current assembly genome (use UCSC genome names)"
  echo ""
  exit 0
fi

if [ "$1" == "-help" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {stringtie_gtf} {reference_genome_name}"
  echo ""
  echo "This script will Overlap StringTie transcripts (GTF format) with UCSC reference genome GTF and annotate novel transcripts with GAWN"
  echo ""
  echo "{stringtie_gtf}: Name of the StringTie GTF file"
  echo ""
  echo "{reference_genome_name}: Name of the current assembly genome (use UCSC genome names)"
  echo ""
  exit 0
fi
if [ "$1" == "--h" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {stringtie_gtf} {reference_genome_name}"
  echo ""
  echo "This script will Overlap StringTie transcripts (GTF format) with UCSC reference genome GTF and annotate novel transcripts with GAWN"
  echo ""
  echo "{stringtie_gtf}: Name of the StringTie GTF file"
  echo ""
  echo "{reference_genome_name}: Name of the current assembly genome (use UCSC genome names)"
  echo ""
  exit 0
fi
if [ "$1" == "--h" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {stringtie_gtf} {reference_genome_name}"
  echo ""
  echo "This script will Overlap StringTie transcripts (GTF format) with UCSC reference genome GTF and annotate novel transcripts with GAWN"
  echo ""
  echo "{stringtie_gtf}: Name of the StringTie GTF file"
  echo ""
  echo "{reference_genome_name}: Name of the current assembly genome (use UCSC genome names)"
  echo ""
  exit 0
fi
if [ "$1" == "--help" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {stringtie_gtf} {reference_genome_name}"
  echo ""
  echo "This script will Overlap StringTie transcripts (GTF format) with UCSC reference genome GTF and annotate novel transcripts with GAWN"
  echo ""
  echo "{stringtie_gtf}: Name of the StringTie GTF file"
  echo ""
  echo "{reference_genome_name}: Name of the current assembly genome (use UCSC genome names)"
  echo ""
  exit 0
fi

[ $# -eq 0 ] && { echo "Usage: bash ./`basename $0` {stringtie_gtf} {reference_genome_name}"; exit 1; }

if [ $# -ne 2 ]; then
  echo 1>&2 "Usage: bash ./`basename $0` {stringtie_gtf} {reference_genome_name}"
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
echo ""
echo "::: Replacing gene_id field in merged.annotated.gtf file with reference gene_id's"
echo ""
########################################
# Stats
########################################
exec 3<> Stats.txt
echo "Number of assembled genes"
cat merged.${1}.tmap | sed "1d" | cut -f4 | sort | uniq | wc -l
echo ""
echo "Number of novel genes"
cat merged.${1}.tmap | awk '$3=="u"{print $0}' | cut -f4 | sort | uniq | wc -l
echo ""
echo "Number of novel transcripts"
cat merged.${1}.tmap | awk '$3=="u"{print $0}' | cut -f5 | sort | uniq | wc -l
echo ""
echo "Number of transcripts matching annotation"
cat merged.${1}.tmap | awk '$3=="="{print $0}' | cut -f5 | sort | uniq | wc -l
exec 3>&-
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
sed -e 's/^/s%/' -e 's/\t/%/' -e 's/$/%g/' namelist |
sed -f - merged.annotated.gtf > merged.annotated_with_reference.gtf
########################################
# Collecting Annotated Transcripts names
########################################
perl -lne 'print "@m" if @m=(/((?:transcript_id|gene_id)\s+\S+)/g);' merged.annotated_with_reference.gtf > transcript_gene_names.txt
sed -i 's/transcript_id //g' transcript_gene_names.txt
sed -i 's/;/\t/g' transcript_gene_names.txt
sed -i 's/gene_id//g' transcript_gene_names.txt
sed -i 's/"//g' transcript_gene_names.txt
awk '{print $1"\t"$2}' transcript_gene_names.txt > transcript_gene_names.tab
cat transcript_gene_names.tab | tr [.] '\t' > gene_names.tab
awk '{print $4"."$3}' gene_names.tab > gene_names2.tab
paste -d'\t' transcript_gene_names.tab gene_names2.tab > filec
awk '{print $1"\t"$3}' filec > gene_names.tab
awk '$2 !~ /STRG./' gene_names.tab > genes.tab
##########################################
# Adding "" to each name in genes.tab file
##########################################
awk '{print $1}' genes.tab > STRG.tab
awk '{print $2}' genes.tab > GENES.tab
sed 's/^/"/' STRG.tab > STRG-left.tab
sed 's/^/"/' GENES.tab > GENES-left.tab
sed 's/$/"/' STRG-left.tab > STRG-left-right.tab
sed 's/$/"/' GENES-left.tab > GENES-left-right.tab
paste -d'\t' STRG-left-right.tab GENES-left-right.tab > filed
uniq filed > genes.tab
rm filed GENES* STRG*
rm gene_names2.tab filec transcript_gene_names.tab transcript_gene_names.txt gene_names.tab
echo "Done"
echo ""
echo "The new annotated GTF is called merged.annotated_with_reference.gtf and is located in the current directory"
echo ""
#rm merged.${1}.tmap.1 merged.${1}.tmap.2
echo ""
echo "::: Obtaining Transcripts in FASTA format with gffread"
echo ""
gffread -w transcripts.fa -g galGal6.fa merged.annotated_with_reference.gtf
echo ""
echo "Done"
echo ""
echo "Downloading GAWN annotation folder. See https://github.com/enormandeau/gawn.git"
echo ""
git clone https://github.com/enormandeau/gawn.git
cd gawn/02_infos/
dir2=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
echo "Done"
echo ""

# Configuring Gawn Inputs and config file

cd /${dir1}/
cp galGal6.fa /${dir1}/gawn/03_data/genome.fasta
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
cd ..
mkdir gffcompare_outputs
mv *.annotated.gtf *.loci *.stats *.refmap *.tmap *.tracking ./gffcompare_outputs
echo ""
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
cp /${dir1}/gawn/05_results/transcriptome_annotation_table.tsv /${dir1}/
cp /${dir1}/CNIT_output/CNIT.index /${dir1}/
sed -i 's/"//g' namelist
awk '{print $1"\t"$2}' transcriptome_annotation_table.tsv > transcriptome_annotation
awk '{print $2"\t"$3"\t"$4"\t"$7}' CNIT.index > CNIT
tail -n +2 transcriptome_annotation > transcriptome_annotation_file
tail -n +2 CNIT > CNIT_file
rm transcriptome_annotation_table.tsv transcriptome_annotation CNIT.index CNIT
join --nocheck-order transcriptome_annotation_file CNIT_file > file1
sed -e 's/^/s%/' -e 's/\t/%/' -e 's/$/%g/' genes.tab |
sed -f - file1 > summary.txt
sed 's/ /\t/g' summary.txt > summary.tab
sed -e 's/^/s%/' -e 's/\t/%/' -e 's/$/%g/' genes.tab |
sed -f - merged.annotated_with_reference.gtf > merged_with_reference_final.gtf
rm transcriptome_annotation_file CNIT_file namelist file1
mkdir merged_annotation
mv merged_with_reference_final.gtf summary.tab ./merged_annotation
rm summary.txt
echo ""
echo "All Done. The transcripts were classified in summary.tab file located in ./merged_annotation"
echo ""
echo "The Classification includes:"
echo ""
echo "Field #1: StringTie transcript id (STRG) including Reference id's"
echo "Field #2: Uniprot closest BLASTp match from GAWN annotation"
echo "Field #3: Coding or non coding classfication by CNIT"
echo "Field #4: Coding probability score (CNIT Score)"
echo "Field #5: Transcript length (bp)"
echo ""
echo "A new GTF file named merged_with_reference_final.gtf is located in ./merged_annotation"
echo ""
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed
#
