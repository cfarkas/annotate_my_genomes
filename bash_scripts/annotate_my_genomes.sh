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
echo "#######################################################"
echo ""
echo "Downloading Reference genome and current GTF annotation"
echo ""
echo "#######################################################"
echo ""
bash genome_download.sh ${2}
echo "################################################"
echo ""
echo "Overlapping StringTie transcripts with Reference"
echo ""
echo "################################################"
echo ""
gffcompare -R -r ${2}.gtf -o merged ${1}
echo "#############################################################################"
echo ""
echo "Replacing gene_id field in merged.annotated.gtf file with reference gene_id's"
echo ""
echo "#############################################################################"
echo ""
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
echo "############################################################"
echo ""
echo "All done. the new annotated GTF is called merged.annotated_with_reference.gtf and is located in the current directory"
echo ""
echo "############################################################"
echo ""
rm merged.${1}.tmap.1 merged.${1}.tmap.2 refGene.txt namelist
echo "##################################################"
echo ""
echo "Obtaining Transcripts in FASTA format with gffread"
echo ""
echo "##################################################"
echo ""
gffread -w transcripts.fa -g galGal6.fa merged.annotated_with_reference.gtf
echo "########################################################################"
echo ""
echo "Downloading GAWN annotation. See https://github.com/enormandeau/gawn.git"
echo ""
echo "########################################################################"
echo ""
git clone https://github.com/enormandeau/gawn.git
cd gawn/02_infos/
dir2=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

# Configuring Gawn Inputs and config file

cd /${dir1}/
cp galGal6.fa /${dir2}/03_data/genome.fasta
cp transcripts.fa /${dir2}/03_data/transcriptome.fasta
rm /${dir2}/gawn_config.sh
cp gawn_config.sh /${dir2}/gawn_config.sh

echo "###################################"
echo ""
echo "Starting GAWN transcript annotation"
echo ""
echo "###################################"
echo ""
cd /${dir1}/gawn/
./gawn 02_infos/gawn_config.sh
echo "###################################################################"
echo ""
echo "All Done. The novel transcripts are annotated in ./gawn/05_results/"
echo ""
echo "###################################################################"
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed
###
