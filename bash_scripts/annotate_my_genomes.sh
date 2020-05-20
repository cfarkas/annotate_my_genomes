#!/bin/bash
{

dir1=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
stringtie_gtf=${1}
reference_genome_gtf=${2}
reference_genome_fasta=${3}
threads=${4}

if [ "$1" == "-h" ]; then
  echo ""
  echo "Usage: ./`basename $0` [stringtie_gtf] [reference_genome_gtf] [reference_genome_fasta] [threads]"
  echo ""
  echo "This pipeline will Overlap StringTie transcripts (GTF format) with current UCSC annotation and will annotate novel transcripts"
  echo ""
  echo "[stringtie_gtf]: StringTie GTF file"
  echo ""
  echo "[reference_genome_gtf]: UCSC reference GTF file."
  echo ""
  echo "[reference_genome_fasta]: Current UCSC assembly genome in fasta format."
  echo ""
  echo "[threads]: Number of threads for parallel text processing (Integer)"
  echo ""
  exit 0
fi

if [ "$1" == "-help" ]; then
  echo ""
  echo "Usage: ./`basename $0` [stringtie_gtf] [reference_genome_gtf] [reference_genome_fasta] [threads]"
  echo ""
  echo "This pipeline will Overlap StringTie transcripts (GTF format) with current UCSC annotation and will annotate novel transcripts"
  echo ""
  echo "[stringtie_gtf]: StringTie GTF file"
  echo ""
  echo "[reference_genome_gtf]: UCSC reference GTF file."
  echo ""
  echo "[reference_genome_fasta]: Current UCSC assembly genome in fasta format."
  echo ""
  echo "[threads]: Number of threads for parallel text processing (Integer)"
  echo ""
  exit 0
fi

if [ "$1" == "--h" ]; then
  echo ""
  echo "Usage: ./`basename $0` [stringtie_gtf] [reference_genome_gtf] [reference_genome_fasta] [threads]"
  echo ""
  echo "This pipeline will Overlap StringTie transcripts (GTF format) with current UCSC annotation and will annotate novel transcripts"
  echo ""
  echo "[stringtie_gtf]: StringTie GTF file"
  echo ""
  echo "[reference_genome_gtf]: UCSC reference GTF file."
  echo ""
  echo "[reference_genome_fasta]: Current UCSC assembly genome in fasta format."
  echo ""
  echo "[threads]: Number of threads for parallel text processing (Integer)"
  echo ""
  exit 0
fi

if [ "$1" == "--help" ]; then
  echo ""
  echo "Usage: ./`basename $0` [stringtie_gtf] [reference_genome_gtf] [reference_genome_fasta] [threads]"
  echo ""
  echo "This pipeline will Overlap StringTie transcripts (GTF format) with current UCSC annotation and will annotate novel transcripts"
  echo ""
  echo "[stringtie_gtf]: StringTie GTF file"
  echo ""
  echo "[reference_genome_gtf]: UCSC reference GTF file."
  echo ""
  echo "[reference_genome_fasta]: Current UCSC assembly genome in fasta format."
  echo ""
  echo "[threads]: Number of threads for parallel text processing (Integer)"
  echo ""
  exit 0
fi

[ $# -eq 0 ] && { echo "Usage: ./`basename $0` [stringtie_gtf] [reference_genome_gtf] [reference_genome_fasta] [threads]"; exit 1; }

if [ $# -ne 4 ]; then
  echo 1>&2 "Usage: ./`basename $0` [stringtie_gtf] [reference_genome_gtf] [reference_genome_fasta] [threads]"
  exit 3
fi

begin=`date +%s`
#    .---------- constant part!
#    vvvv vvvv-- the code from above
YELLOW='\033[1;33m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color
echo ""
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 1. Overlapping StringTie transcripts with Reference :::\n"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
echo ""
gffcompare -R -r ${2} -s ${3} -o UCSC_compare ${1}
printf "${PURPLE}Done\n"
echo ""
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 2. Writting novel discoveries to Stats.txt :::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
echo ""
# Stats
exec 3<> Stats.txt
echo "Number of assembled genes:" >> Stats.txt
cat UCSC_compare.${1}.tmap | sed "1d" | cut -f4 | sort | uniq | wc -l >> Stats.txt
echo "" >> Stats.txt
echo "Number of novel genes:" >> Stats.txt
cat UCSC_compare.${1}.tmap | awk '$3=="u"{print $0}' | cut -f4 | sort | uniq | wc -l >> Stats.txt
echo "" >> Stats.txt
echo "Number of novel transcripts:" >> Stats.txt
cat UCSC_compare.${1}.tmap | awk '$3=="u"{print $0}' | cut -f5 | sort | uniq | wc -l >> Stats.txt
echo "" >> Stats.txt
echo "Number of transcripts matching annotation:" >> Stats.txt
cat UCSC_compare.${1}.tmap | awk '$3=="="{print $0}' | cut -f5 | sort | uniq | wc -l >> Stats.txt
exec 3>&-
printf "${PURPLE}Done\n"
echo ""
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 3. Replacing gene_id field in merged.annotated.gtf file with reference gene_id's :::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
echo ""
########################################
# Merging novel transcripts with ref. 
########################################
awk '{print $4"\t"$1}' UCSC_compare.${1}.tmap > UCSC_compare.${1}.tmap.1
tail -n +2 UCSC_compare.${1}.tmap.1 > UCSC_compare.${1}.tmap.2
awk '!/-/' UCSC_compare.${1}.tmap.2 > namelist
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
cat ${1} | parallel --pipe -j ${4} sed -f sed.script > merged.gtf
rm -f sed.script fileA fileB
printf "${PURPLE}::: Done. Gene_id field was replaced in the stringtie GTF file and merged.gtf was generated with these changes\n"
echo ""
#######################################
# Re-formatting merged.gtf file #
#######################################
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 4. Re-formatting merged.gtf using standard gff/gtf specifications with gffread:::\n"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
gffread -E -F --merge merged.gtf -o merged.gff3
gffread merged.gff3 -T -o merged_fixed.gtf
rm merged.gtf merged.gff3 
echo ""
echo ""
printf "${PURPLE}::: Re-formatting was done. The new GTF file is called merged_fixed.gtf\n"
printf "${PURPLE}::: Continue with FEELnc long non-coding classification...\n"
echo ""
############################################
# FEELnc long noncoding RNA identification #
############################################
cd /${dir1}/
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 5. Classifying protein-coding and long non-coding transcripts with FEELnc :::\n"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
### Cloning FEELnc in current directory
git clone https://github.com/tderrien/FEELnc.git
echo ""
cp ${3} ${2} merged_fixed.gtf /${dir1}/FEELnc/
cd FEELnc
export FEELNCPATH=${PWD}
export PERL5LIB=$PERL5LIB:${FEELNCPATH}/lib/ #order is important to avoid &Bio::DB::IndexedBase::_strip_crnl error with bioperl >=v1.7
export PATH=$PATH:${FEELNCPATH}/scripts/
export PATH=$PATH:${FEELNCPATH}/utils/
echo ""
### Testing FEELnc first
echo ""
printf "${PURPLE}::: Testing if FEELnc works ...\n"
echo ""
cd test/
# Filter
FEELnc_filter.pl -i transcript_chr38.gtf -a annotation_chr38.gtf -b transcript_biotype=protein_coding > candidate_lncRNA.gtf
# Coding_Potential
FEELnc_codpot.pl -i candidate_lncRNA.gtf -a annotation_chr38.gtf -b transcript_biotype=protein_coding -g genome_chr38.fa --mode=shuffle
# Classifier
FEELnc_classifier.pl -i feelnc_codpot_out/candidate_lncRNA.gtf.lncRNA.gtf -a annotation_chr38.gtf > candidate_lncRNA_classes.txt
echo ""
printf "${PURPLE}::: FEELnc Test done. Continue with merged_fixed.gtf file :::\n"
echo ""
cd ..
echo ""
### Running FEELnc
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 6.  Running FEELnc on merged_fixed.gtf file ...\n"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
# Filter
FEELnc_filter.pl -i merged_fixed.gtf -a ${2} -b transcript_biotype=protein_coding > candidate_lncRNA.gtf
# Coding_Potential
FEELnc_codpot.pl -i candidate_lncRNA.gtf -a ${2} -b transcript_biotype=protein_coding -g ${3} --mode=shuffle
# Classifier
FEELnc_classifier.pl -i feelnc_codpot_out/candidate_lncRNA.gtf.lncRNA.gtf -a ${2} > candidate_lncRNA_classes.txt
echo ""
printf "${PURPLE}::: FEELnc calculations were done :::\n"
echo ""
printf "${YELLOW}::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 7. Parsing FEELnc output :::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::${CYAN}\n"
cp candidate_lncRNA_classes.txt /${dir1}/
cd /${dir1}/
awk '{print $3}' candidate_lncRNA_classes.txt > lncRNA_genes
tail -n +2 lncRNA_genes > lncRNA_transcripts
rm lncRNA_genes
grep -w -F -f lncRNA_transcripts merged_fixed.gtf > merged.fixed.lncRNAs.gtf
grep --invert-match -F -f lncRNA_transcripts merged_fixed.gtf > merged.fixed.coding.gtf
sed -i 's/StringTie/lncRNA/' merged.fixed.lncRNAs.gtf
sed -i 's/StringTie/coding/' merged.fixed.coding.gtf
cat merged.fixed.coding.gtf merged.fixed.lncRNAs.gtf > final_annotated.gtf
gffread -E -F --merge final_annotated.gtf -o final_annotated.gff
rm merged_fixed.gtf
echo ""
printf "${PURPLE}::: Parsing is done. The transcripts were classified and added to final_annotated.gtf file...\n"
echo ""
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 8. Obtaining Transcripts in FASTA format with gffread :::\n"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
echo ""
gffread -w transcripts.fa -g ${3} final_annotated.gtf
echo ""
printf "${PURPLE}::: Done. transcripts.fa are located in current directory\n"
echo ""
echo ""
printf "${PURPLE}::: Moving gffcompare results to gffcompare_outputs folder ...\n"
echo ""
mkdir gffcompare_outputs
mv *.loci *.stats *.refmap *.tmap *.tracking ./gffcompare_outputs
echo ""
echo "Done"
echo ""
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 9. Performing gene annotation by using GAWN pipeline :::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
################################################################
# Configuring Gawn Inputs, config file and running GAWN pipeline
################################################################
echo ""
printf "${PURPLE}::: Downloading GAWN annotation folder. See https://github.com/enormandeau/gawn.git${CYAN}\n"
echo ""
git clone https://github.com/enormandeau/gawn.git
cd gawn/02_infos/
dir2=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
echo "Done"
echo ""
cd /${dir1}/
cp ${3} /${dir1}/gawn/03_data/genome.fasta
cp transcripts.fa /${dir1}/gawn/03_data/transcriptome.fasta
rm /${dir2}/gawn_config.sh
cp gawn_config.sh /${dir2}/gawn_config.sh
echo ""
printf "${PURPLE}::: Starting GAWN transcript annotation${CYAN}\n"
echo ""
cd /${dir1}/gawn/
./gawn 02_infos/gawn_config.sh
echo ""
echo ""
printf "${PURPLE}::: Done. The novel transcripts were annotated in ./gawn/05_results/ :::${CYAN}\n"
echo ""
###########################################
# Extracting GO terms for each transcript #
###########################################
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 10. Extracting GO terms for each transcript :::\n"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
echo ""
cd /${dir1}/
cp /${dir1}/gawn/05_results/transcriptome_annotation_table.tsv /${dir1}/
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
echo ""
printf "${PURPLE}::: Done. GO terms were succesfully extracted :::${CYAN}\n"
echo ""
######################################
# Gene Prediction Step with Augustus #
######################################
cd /${dir1}/
echo ""
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 11. Predicting gene models from transcripts with AUGUSTUS (gff3 format) :::\n"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
echo ""
printf "${PURPLE}::: Progress will be printed for each transcript :::\n"
echo ""
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
echo ""
printf "${PURPLE}::: Done. augustus.gff3 file is present in current directory...${CYAN}\n"
echo ""
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 12. Converting gff3 to GTF format, collecting coding sequences and proteins with gffread and AGAT :::\n"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
gffread augustus.gff3 -T -o coding_transcripts.gtf
agat_sp_extract_sequences.pl -g augustus.gff3 -f transcripts.fa -o cds.fa
agat_sp_extract_sequences.pl -g augustus.gff3 -f transcripts.fa -o prot.fa --protein
#############################
# Configuring Summary Results
#############################
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 13. Configuring Summary Results :::\n"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::${CYAN}\n"
############################################
# Moving results to merged_annotation folder
############################################
echo ""
printf "${PURPLE}::: Moving results to output_files folder :::${CYAN}\n"
mkdir output_files
mv candidate_lncRNA_classes.txt final_annotated.gtf final_annotated.gff Stats.txt transcripts.fa transcriptsGO.tab genesGO.tab cds.fa prot.fa coding_transcripts.gtf logfile augustus.gff3 ./output_files
cp /${dir1}/gawn/05_results/transcriptome_annotation_table.tsv /${dir1}/output_files/
rm transcripts.fa.fai namelist* isoforms_per_gene_concatenated.tab lncRNA_transcripts merged.fixed.coding.gtf merged.fixed.lncRNAs.gtf merged.gtf UCSC_compare* transcriptome_annotation_table.tsv
rm refGene.tx*
echo ""
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
echo "All Done. The transcripts were classified in ./output_files"
echo ""
echo "Transcript discoveries are summarized in Stats.txt file located in ./output_files . GAWN annotation is named transcriptome_annotation_table.tsv"
echo ""
echo "GTF file final_annotated.gtf (standard GTF) and correspondent gff file (final_annotated.gff) are located in ./output_files. These files contains the annotated lncRNA/coding GTF in the second field".
echo ""
echo "candidate_lncRNA_classes.txt contained detailed long non-coding classification of transcripts".
echo ""
echo "Associated FASTA file to this GTF, named transcripts.fa is located in ./output_files"
echo ""
echo "AUGUSTUS GTF file suitable for transcript count quantification is named coding_transcripts.gtf. This GTF file contains all coding transcripts resolved by AUGUSTUS and is located in ./output_files"
echo ""
echo "Associated Transcript coding sequences (cds.fa) and correspondent protein sequences (prot.fa) with coding_transcripts.gtf are located in ./output_files"
echo ""
echo "GO terms associated to each transcript (and gene), named transcriptsGO.tab and genesGO.tab are located in ./output_files"
echo ""
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${NC}\n"
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed
#
} | tee logfile
#
