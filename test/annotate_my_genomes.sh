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
rm A A.1 A.2 B B.1 B.2

###############################
# Getting gene names replaced #
###############################
awk '{print $1}' namelist > fileA
awk '{print $2}' namelist > fileB
paste -d : fileA fileB | sed 's/\([^:]*\):\([^:]*\)/s%\1%\2%/' > sed.script
cat merged.annotated.gtf | parallel --pipe -j ${3} sed -f sed.script > merged_with_reference.gtf
rm -f sed.script fileA fileB
echo "Done"
echo ""
################################
# Formatting Transcripts names #
################################
grep -v "gene_id \"STRG." merged_with_reference.gtf > annotated_genes.gtf
sed 's/\ /\t/g' annotated_genes.gtf > annotated_genes.tab
perl -lne 'print "@m" if @m=(/((?:transcript_id|gene_id)\s+\S+)/g);' annotated_genes.gtf > transcript_gene_names.txt
sed -i 's/transcript_id //g' transcript_gene_names.txt
sed -i 's/;/\t/g' transcript_gene_names.txt
sed -i 's/gene_id//g' transcript_gene_names.txt
sed 's/"//g' transcript_gene_names.txt > outfile
sed -i 's/"//g' transcript_gene_names.txt
awk '{print $1"\t"$2}' transcript_gene_names.txt > transcript_gene_names.tab
tr '.' '\t' < transcript_gene_names.tab > transcripts_sep
awk '{print $4"."$3}' < transcripts_sep > filea
sed 's/^/"/' filea > filea-left
sed 's/$/";/' filea-left > filea-left-right
awk 'FNR==NR{a[NR]=$1;next}{$10=a[FNR]}1' filea-left-right annotated_genes.tab > merged
sed 's/\ /\t/g' merged > annotated_genes.tab
grep "gene_id \"STRG." merged_with_reference.gtf > STRG_genes.gtf
cat annotated_genes.tab STRG_genes.gtf > merged.gtf

##################
# Validating GTF #
##################
cd /${dir1}/
cd ..
perl validate_gtf.pl -f /${dir1}/merged.gtf
cd /${dir1}/
echo ""
echo "The merged.gtf file was succesfully validated"
rm outfile merged transcript_gene_names.txt transcripts_sep transcript_gene_names.tab filea filea-left filea-left-right annotated_genes.tab annotated_genes.gtf STRG_genes.gtf merged.gtf
echo ""
echo "A new annotated GTF is called merged.fixed.gtf and is located in the current directory ..."
echo ""
echo "::: Obtaining Transcripts in FASTA format with gffread"
echo ""
gffread -w transcripts.fa -g ${2}.fa merged.fixed.gtf
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

############################################
# FEELnc long noncoding RNA identification #
############################################
cd /${dir1}/
echo "::: Classifying protein-coding and long non-coding transcripts with FEELnc"
git clone https://github.com/tderrien/FEELnc.git
echo ""
cp {reference_genome_name}.fa {reference_genome_name}.gtf merged.fixed.gtf /${dir1}/FEELnc/
### Cloning FEELnc in current directory
git clone https://github.com/tderrien/FEELnc.git
cd FEELnc
export FEELNCPATH=${PWD}
export PERL5LIB=$PERL5LIB:${FEELNCPATH}/lib/ #order is important to avoid &Bio::DB::IndexedBase::_strip_crnl error with bioperl >=v1.7
export PATH=$PATH:${FEELNCPATH}/scripts/
export PATH=$PATH:${FEELNCPATH}/utils/
echo ""
### Testing FEELnc first
echo "Testing FEELnc is works ..."
cd test/
# Filter
FEELnc_filter.pl -i transcript_chr38.gtf -a annotation_chr38.gtf -b transcript_biotype=protein_coding > candidate_lncRNA.gtf
# Coding_Potential
FEELnc_codpot.pl -i candidate_lncRNA.gtf -a annotation_chr38.gtf -b transcript_biotype=protein_coding -g genome_chr38.fa --mode=shuffle
# Classifier
FEELnc_classifier.pl -i feelnc_codpot_out/candidate_lncRNA.gtf.lncRNA.gtf -a annotation_chr38.gtf > candidate_lncRNA_classes.txt
echo "FEELnc Test done"
cd ..
echo ""
### Running FEELnc
echo "Running FEELnc on merged.fixed.gtf file ..."
# Filter
FEELnc_filter.pl -i merged.fixed.gtf -a {reference_genome_name}.gtf -b transcript_biotype=protein_coding > candidate_lncRNA.gtf
# Coding_Potential
FEELnc_codpot.pl -i candidate_lncRNA.gtf -a {reference_genome_name}.gtf -b transcript_biotype=protein_coding -g {reference_genome_name}.fa --mode=shuffle
# Classifier
FEELnc_classifier.pl -i feelnc_codpot_out/candidate_lncRNA.gtf.lncRNA.gtf -a {reference_genome_name}.gtf > candidate_lncRNA_classes.txt
echo "FEELnc calculations were done"
echo ""
echo "::: Parsing FEELnc output"
cp candidate_lncRNA_classes.txt /${dir1}/
cd /${dir1}/
awk '{print $3}' candidate_lncRNA_classes.txt > lncRNA_genes
tail -n +2 lncRNA_genes > lncRNA_transcripts
rm lncRNA_genes
grep -w -F -f lncRNA_transcripts merged.fixed.gtf > merged.fixed.lncRNAs.gtf
grep --invert-match -F -f lncRNA_transcripts merged.fixed.gtf > merged.fixed.coding.gtf
sed -i 's/StringTie/lncRNA/' merged.fixed.lncRNAs.gtf
sed -i 's/StringTie/coding/' merged.fixed.coding.gtf
cat merged.fixed.coding.gtf merged.fixed.lncRNAs.gtf > final_annotated.gtf
echo ""
echo "Done. The transcripts were classified and added to final_annotated.gtf file..."
echo ""

##########################################
# Extracting GO terms for each transcript
##########################################
cut -d$'\t' -f 1,6 transcriptome_annotation_table.tsv > transcripts_GO
tr ';' '\t' < transcripts_GO > transcripts_GO_sep
column -t transcripts_GO_sep > transcripts_GO.tab
tail -n +2 transcripts_GO.tab > transcriptsGO.tab
rm transcripts_GO*

##########################################
# Extracting GO terms for each Gene
##########################################

grep "STRG." transcriptsGO.tab > STRG_transcriptsGO.tab
grep -v "STRG." transcriptsGO.tab > Annotated_transcriptsGO.tab

# Working with "STRG" genes
tr '.' '\t' < STRG_transcriptsGO.tab > STRG_transcripts_GO_splitname
awk '{$3=""; print $0}' STRG_transcripts_GO_splitname > STRG_transcripts_GO_splitname1
awk '{print $1"."$2}' < STRG_transcripts_GO_splitname1 > STRG_transcripts_GO_splitname2
paste -d'\t' STRG_transcripts_GO_splitname2 STRG_transcripts_GO_splitname1 > STRG_transcripts_GO_merged
awk '{$2=$3=""; print $0}' STRG_transcripts_GO_merged > STRG_genes_GO
awk '!a[$0]++' STRG_genes_GO > STRG_genesGO
awk '{print NF}' STRG_genesGO > STRG_numbers
paste -d'\t' STRG_numbers STRG_genesGO > STRG_genesGO_with_numbers
awk '{ if ($1>1) { print } }' STRG_genesGO_with_numbers > STRG_GO
awk '{$1=""; print $0}' STRG_GO > STRG_genes_with_GO
column -t STRG_genes_with_GO > STRG_genes_withGO
rm STRG_genes_GO STRG_transcripts_GO_merged STRG_transcripts_GO_splitname* STRG_transcriptsGO.tab STRG_numbers STRG_genesGO STRG_genesGO_with_numbers STRG_genes_with_GO STRG_GO

# Working with Annotated genes
tr '.' '\t' < Annotated_transcriptsGO.tab > Annotated_transcripts_GO_splitname
awk '{$2=""; print $0}' Annotated_transcripts_GO_splitname > Annotated_genes_GO
awk '!a[$0]++' Annotated_genes_GO > Annotated_genesGO
awk '{print NF}' Annotated_genesGO > Annotated_numbers
paste -d'\t' Annotated_numbers Annotated_genesGO > Annotated_genesGO_with_numbers
awk '{ if ($1>1) { print } }' Annotated_genesGO_with_numbers > Annotated_GO
awk '{$1=""; print $0}' Annotated_GO > Annotated_genes_with_GO
column -t Annotated_genes_with_GO > genes_withGO.tab
rm Annotated_transcriptsGO.tab Annotated_transcripts_GO_splitname Annotated_genes_GO Annotated_genesGO Annotated_numbers Annotated_genesGO_with_numbers Annotated_GO Annotated_genes_with_GO

# Joining Files in order to create "genesGO.tab" file 
cat STRG_genes_withGO >> genes_withGO.tab
sed "s/^ *//;s/ *$//;s/ \{1,\}/ /g" genes_withGO.tab > genes_with_GO.tab
sed 's/ /\t/' genes_with_GO.tab > genesGO.tab
awk '{$2=""; print $0}' genesGO.tab > genesGO1.tab
rm genesGO.tab
awk '!a[$0]++' genesGO1.tab > genesGO.tab
rm genesGO1.tab
rm STRG_genes_withGO genes_withGO.tab genes_with_GO.tab
sed 's/ /\t/' genesGO.tab > genesGO1.tab
mv genesGO1.tab genesGO.tab

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
./augustus.2.5.5/src/augustus --species=human --progress=true --UTR=off --uniqueGeneId=true --gff3=on transcripts.fa > augustus.gff3
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
rm coding_transcript.gtf 

#############################
# Configuring Summary Results
#############################

############################################
# Moving results to merged_annotation folder
############################################
echo ""
echo "Moving results to merged_annotation folder"
mkdir output_files
mv candidate_lncRNA_classes.txt final_annotated.gtf Stats.txt transcripts.fa transcriptsGO.tab genesGO.tab transcripts_CDS.fa transcripts_proteins.fa coding_transcripts.gtf logfile ./output_files
cp /${dir1}/gawn/05_results/transcriptome_annotation_table.tsv /${dir1}/output_files/
rm genes1.tab genes2.tab transcripts.fa.fai namelist namelist_unique_sorted transcripts_conc.tab
echo ""
echo "All Done. The transcripts were classified in ./output_files"
echo ""
echo "Transcript discoveries are summarized in Stats.txt file located in ./output_files . GAWN annotation is named transcriptome_annotation_table.tsv"
echo ""
echo "A new GTF file suitable for gene count quantification is named final_annotated.gtf and is located in ./output_files. This GTF contains the annotated lncRNA/coding GTF in the second field".
echo ""
echo "candidate_lncRNA_classes.txt contained detailed long non-coding classification of transcripts".
echo ""
echo "Associated FASTA file to this GTF, named transcripts.fa is located in ./output_files"
echo ""
echo "AUGUSTUS GTF file suitable for transcript count quantification is named coding_transcripts.gtf. This GTF file contains all coding transcripts resolved by AUGUSTUS and is located in ./output_files"
echo ""
echo "Associated Transcript coding sequences (transcripts_CDS.fa) and correspondent protein sequences (transcripts_proteins.fa) with coding_transcripts.gtf are located in ./output_files"
echo ""
echo "GO terms associated to each transcript (and gene), named transcriptsGO.tab and genesGO.tab are located in ./output_files"
echo ""
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed
#
} | tee logfile
#
