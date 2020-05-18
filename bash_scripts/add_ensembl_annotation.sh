#!/bin/bash
{

dir1=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
final_annotated=${1}
Ensembl_reference_genome_gtf=${2}
Ensembl_reference_genome_fasta=${3}
UCSC_genome_prefix=${4}
threads=${5}

if [ "$1" == "-h" ]; then
  echo ""
  echo "Usage: ./`basename $0` [final_annotated] [Ensembl_reference_genome_gtf] [Ensembl_reference_genome_fasta] [UCSC_genome_prefix] [threads]"
  echo ""
  echo "This pipeline will Overlap final_annotated.gtf transcripts (GTF format) with current Ensembl annotation"
  echo ""
  echo "[final_annotated]: annotate_my_genomes gtf output (final_annotated.gtf)"
  echo ""
  echo "[Ensembl_reference_genome_gtf]: Ensembl reference GTF file, available here: ftp://ftp.ensembl.org/pub/release-100/gtf/ (check XXX.100.gtf.gz)"
  echo ""
  echo "[Ensembl_reference_genome_fasta]: Ensembl genome assembly in fasta format, available here: ftp://ftp.ensembl.org/pub/release-100/fasta/ (check dna, XXX.dna.toplevel.fa.gz file)"
  echo ""
  echo "[UCSC_genome_prefix]: UCSC genome assembly in fasta format, check here: https://genome.ucsc.edu/cgi-bin/hgGateway"
  echo ""
  echo "[threads]: Number of threads for parallel text processing (Integer)"
  echo ""
  exit 0
fi

if [ "$1" == "-help" ]; then
  echo ""
  echo "Usage: ./`basename $0` [final_annotated] [Ensembl_reference_genome_gtf] [Ensembl_reference_genome_fasta] [UCSC_genome_prefix] [threads]"
  echo ""
  echo "This pipeline will Overlap final_annotated.gtf transcripts (GTF format) with current Ensembl annotation"
  echo ""
  echo "[final_annotated]: annotate_my_genomes gtf output (final_annotated.gtf)"
  echo ""
  echo "[Ensembl_reference_genome_gtf]: Ensembl reference GTF file, available here: ftp://ftp.ensembl.org/pub/release-100/gtf/ (check XXX.100.gtf.gz)"
  echo ""
  echo "[Ensembl_reference_genome_fasta]: Ensembl genome assembly in fasta format, available here: ftp://ftp.ensembl.org/pub/release-100/fasta/ (check dna, XXX.dna.toplevel.fa.gz file)"
  echo ""
  echo "[UCSC_genome_prefix]: UCSC genome assembly in fasta format, check here: https://genome.ucsc.edu/cgi-bin/hgGateway"
  echo ""
  echo "[threads]: Number of threads for parallel text processing (Integer)"
  echo ""
  exit 0
fi

if [ "$1" == "--h" ]; then
  echo ""
  echo "Usage: ./`basename $0` [final_annotated] [Ensembl_reference_genome_gtf] [Ensembl_reference_genome_fasta] [UCSC_genome_prefix] [threads]"
  echo ""
  echo "This pipeline will Overlap final_annotated.gtf transcripts (GTF format) with current Ensembl annotation"
  echo ""
  echo "[final_annotated]: annotate_my_genomes gtf output (final_annotated.gtf)"
  echo ""
  echo "[Ensembl_reference_genome_gtf]: Ensembl reference GTF file, available here: ftp://ftp.ensembl.org/pub/release-100/gtf/ (check XXX.100.gtf.gz)"
  echo ""
  echo "[Ensembl_reference_genome_fasta]: Ensembl genome assembly in fasta format, available here: ftp://ftp.ensembl.org/pub/release-100/fasta/ (check dna, XXX.dna.toplevel.fa.gz file)"
  echo ""
  echo "[UCSC_genome_prefix]: UCSC genome assembly in fasta format, check here: https://genome.ucsc.edu/cgi-bin/hgGateway"
  echo ""
  echo "[threads]: Number of threads for parallel text processing (Integer)"
  echo ""
  exit 0
fi

if [ "$1" == "--help" ]; then
  echo ""
  echo "Usage: ./`basename $0` [final_annotated] [Ensembl_reference_genome_gtf] [Ensembl_reference_genome_fasta] [UCSC_genome_prefix] [threads]"
  echo ""
  echo "This pipeline will Overlap final_annotated.gtf transcripts (GTF format) with current Ensembl annotation"
  echo ""
  echo "[final_annotated]: annotate_my_genomes gtf output (final_annotated.gtf)"
  echo ""
  echo "[Ensembl_reference_genome_gtf]: Ensembl reference GTF file, available here: ftp://ftp.ensembl.org/pub/release-100/gtf/ (check XXX.100.gtf.gz)"
  echo ""
  echo "[Ensembl_reference_genome_fasta]: Ensembl genome assembly in fasta format, available here: ftp://ftp.ensembl.org/pub/release-100/fasta/ (check dna, XXX.dna.toplevel.fa.gz file)"
  echo ""
  echo "[UCSC_genome_prefix]: UCSC genome assembly in fasta format, check here: https://genome.ucsc.edu/cgi-bin/hgGateway"
  echo ""
  echo "[threads]: Number of threads for parallel text processing (Integer)"
  echo ""
  exit 0
fi

[ $# -eq 0 ] && { echo "Usage: ./`basename $0` [final_annotated] [Ensembl_reference_genome_gtf] [Ensembl_reference_genome_fasta] [UCSC_genome_prefix] [threads]"; exit 1; }

if [ $# -ne 5 ]; then
  echo 1>&2 "Usage: ./`basename $0` [final_annotated] [Ensembl_reference_genome_gtf] [Ensembl_reference_genome_fasta] [UCSC_genome_prefix] [threads]"
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
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 1. Donwload USCS reference genome in fasta format :::\n"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
echo ""
rm -r -f augustus.* FEELnc gawn 
./genome-download ${4}
printf "${PURPLE}Done. ${4}.fa is located in current directory\n"
echo ""
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 2. Obtaining Ensembl transcripts in fasta format :::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
echo ""
gffread -w ensembl_transcripts.fa -g ${3} ${2}
printf "${PURPLE}Done\n"
echo ""
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 3.  Mapping ensembl transcripts to UCSC genome, using ${5} threads :::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
echo ""
minimap2 -ax splice ${4}.fa ensembl_transcripts.fa > ensembl_aligned.sam -t ${5}
samtools view -S -b ensembl_aligned.sam -@ ${5} > ensembl_aligned.bam
samtools sort ensembl_aligned.bam -@ ${5} > ensembl_aligned.sorted.bam
samtools index ensembl_aligned.sorted.bam -@ ${5}
printf "${PURPLE}Done. Mapped transcripts are called ensembl_aligned.sorted.bam\n"
echo ""
rm ensembl_aligned.sam ensembl_aligned.bam
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 4. Obtaining bed file from alignments by using bedtools bamtobed :::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
echo ""
bamToBed -i ensembl_aligned.sorted.bam > ensembl_aligned.bed
printf "${PURPLE}Done\n"
echo ""
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 5. Converting bed file to GTF using AGAT and gffread :::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
agat_convert_bed2gff.pl --bed ensembl_aligned.bed -o ensembl_aligned.gff
sed -i 's/ID=//'g ensembl_aligned.gff
sed -i 's/Name=/ID=/'g ensembl_aligned.gff
gffread ensembl_aligned.gff --gene2exon -o ensembl_aligned.gff3
gffread ensembl_aligned.gff3 -T -o ensembl_aligned.gtf
echo ""
printf "${PURPLE}Done. ensembl_aligned.gtf contain Ensembl transcripts mapped to UCSC genome coordinates\n"
echo ""
printf "${PURPLE}::: Removing intermediate files\n"
rm ensembl_aligned.bed ensembl_aligned.sorted.bam* ensembl_aligned.sam ensembl_aligned.bam ensembl_aligned.gf*
printf "${PURPLE}Done\n"
echo ""
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 6. Overlapping final_annotated.gtf transcripts with Ensembl GTF :::\n"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
echo ""
grep "STRG." ${1} > STRG.gtf
grep -v "STRG." ${1} > non_STRG.gtf
gffcompare -R -r ensembl_aligned.gtf -s ${4}.fa -o UCSC_compare STRG.gtf
printf "${PURPLE}Done\n"
echo ""
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 7. Writting novel discoveries to Stats.txt :::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
echo ""
# Stats
exec 3<> Stats.txt
echo "Number of assembled genes:" >> Stats.txt
cat UCSC_compare.STRG.gtf.tmap | sed "1d" | cut -f4 | sort | uniq | wc -l >> Stats.txt
echo "" >> Stats.txt
echo "Number of novel genes:" >> Stats.txt
cat UCSC_compare.STRG.gtf.tmap | awk '$3=="u"{print $0}' | cut -f4 | sort | uniq | wc -l >> Stats.txt
echo "" >> Stats.txt
echo "Number of novel transcripts:" >> Stats.txt
cat UCSC_compare.STRG.gtf.tmap | awk '$3=="u"{print $0}' | cut -f5 | sort | uniq | wc -l >> Stats.txt
echo "" >> Stats.txt
echo "Number of transcripts matching annotation:" >> Stats.txt
cat UCSC_compare.STRG.gtf.tmap | awk '$3=="="{print $0}' | cut -f5 | sort | uniq | wc -l >> Stats.txt
exec 3>&-
printf "${PURPLE}Done\n"
echo ""
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 8. Replacing gene_id field in final_annotated.gtf file with Ensenbl gene_id's :::\n"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
echo ""
########################################
# Merging novel transcripts with ref. 
########################################
awk '{print $4"\t"$1}' UCSC_compare.STRG.gtf.tmap > UCSC_compare.STRG.gtf.tmap.1
tail -n +2 UCSC_compare.STRG.gtf.tmap.1 > UCSC_compare.STRG.gtf.tmap.2
awk '!/-/' UCSC_compare.STRG.gtf.tmap.2 > namelist
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
cat STRG.gtf | parallel --pipe -j ${5} sed -f sed.script > merged_with_reference.gtf
rm -f sed.script fileA fileB
printf "${PURPLE}::: Done. Gene_id field was replaced in the final_annotated.gtf file and merged_with_reference.gtf was generated with these changes\n"
echo ""
printf "${YELLOW}::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 9. Formatting Isoforms :::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::${CYAN}\n"
################################
# Formatting Transcripts names #
################################
# Extracting replaced genes names from merged_with_reference.gtf file
grep -v "gene_id \"STRG." merged_with_reference.gtf > annotated_genes.gtf
sed 's/\ /\t/g' annotated_genes.gtf > annotated_genes.tab
perl -lne 'print "@m" if @m=(/((?:transcript_id|gene_id)\s+\S+)/g);' annotated_genes.gtf > transcript_gene_names.txt
sed -i 's/transcript_id //g' transcript_gene_names.txt
sed -i 's/;/\t/g' transcript_gene_names.txt
sed -i 's/gene_id//g' transcript_gene_names.txt
sed -i 's/"//g' transcript_gene_names.txt
sed -i 's/"//g' transcript_gene_names.txt
# generating replaced gene names with matched original stringtie isoforms
awk '{print $1"\t"$2}' transcript_gene_names.txt > transcript_gene_names.tab
# removing duplicates
awk '!a[$0]++' transcript_gene_names.tab > transcript_gene_names.unique.tab
# selecting column of replaced genes names and iterate numbers to obtain fixed isoforms numbers
awk '{print $1}' < transcript_gene_names.unique.tab > replaced_gene_names.tab
# iterate numbers in each unique gene_id
awk '{ printf "%06d.%d\t%s\n",(!a[$1]++? ++c:c),a[$1],$0 }' replaced_gene_names.tab > replaced_gene_names_iterate.tab
# generating isoforms IDs
tr '.' '\t' < replaced_gene_names_iterate.tab > replaced_gene_names_iterate_sep.tab
awk '{print $5}' replaced_gene_names_iterate_sep.tab > isoforms_per_gene
# align transcript_gene_names.unique.tab with new isoforms
paste -d'\t' transcript_gene_names.unique.tab isoforms_per_gene > isoforms_per_gene_concatenated
awk '{print $1"\t"$2"."$3}' isoforms_per_gene_concatenated > isoforms_per_gene_concatenated.tab
# generate file for sed script, as "namelist"
awk '{print $1}' isoforms_per_gene_concatenated.tab  > A
awk '{print $2}' isoforms_per_gene_concatenated.tab  > B
sed 's/^/"/' A > A.1
sed 's/$/"/' A.1 > A.2
sed 's/^/"/' B > B.1
sed 's/$/"/' B.1 > B.2
paste -d'\t' A.2 B.2 > namelist_isoforms
rm A A.1 A.2 B B.1 B.2 transcript_gene* isoforms_per_gene isoforms_per_gene_concatenated replaced_*
##################################
# Getting isoform names replaced #
##################################
awk '{print $1}' namelist_isoforms > fileA
awk '{print $2}' namelist_isoforms > fileB
paste -d : fileA fileB | sed 's/\([^:]*\):\([^:]*\)/s%\1%\2%/' > sed.script
cat merged_with_reference.gtf | parallel --pipe -j ${5} sed -f sed.script > merged.gtf
rm -f sed.script fileA fileB annotated_genes*
echo ""
printf "${PURPLE}::: Done. Gene_id field was replaced and intermediate merged.gtf file was generated with these changes. Continue with GTF validation\n"
echo ""
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 10. Re-generating final GTF file to annotate :::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
echo ""
cat non_STRG.gtf merged.gtf > annotated_ensembl.gtf
sed -i 's/lncRNA/StringTie/' annotated_ensembl.gtf
sed -i 's/coding/StringTie/' annotated_ensembl.gtf
printf "${PURPLE}::: All done. annotated_ensembl.gtf contained full classification of transcripts. Continue with lncRNA classification ...\n"
echo ""
############################################
# FEELnc long noncoding RNA identification #
############################################
cd /${dir1}/
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 11. Classifying protein-coding and long non-coding transcripts with FEELnc :::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
rm -r -f FEELnc
git clone https://github.com/tderrien/FEELnc.git
echo ""
cp ensembl_aligned.gtf ${4}.fa annotated_ensembl.gtf /${dir1}/FEELnc/
### Cloning FEELnc in current directory
git clone https://github.com/tderrien/FEELnc.git
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
printf "${PURPLE}::: FEELnc Test done. Continue with final_annotated_ensembl.gtf file :::\n"
echo ""
cd ..
echo ""
### Running FEELnc
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 12.  Running FEELnc on final_annotated_ensembl.gtf file :::\n"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
# Filter
FEELnc_filter.pl -i annotated_ensembl.gtf -a ensembl_aligned.gtf -b transcript_biotype=protein_coding > candidate_lncRNA.gtf
# Coding_Potential
FEELnc_codpot.pl -i candidate_lncRNA.gtf -a ensembl_aligned.gtf -b transcript_biotype=protein_coding -g ${4}.fa --mode=shuffle
# Classifier
FEELnc_classifier.pl -i feelnc_codpot_out/candidate_lncRNA.gtf.lncRNA.gtf -a ensembl_aligned.gtf > candidate_lncRNA_classes.txt
echo ""
printf "${PURPLE}::: FEELnc calculations were done :::\n"
echo ""
printf "${YELLOW}::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 13. Parsing FEELnc output :::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::${CYAN}\n"
cp candidate_lncRNA_classes.txt /${dir1}/
cd /${dir1}/
awk '{print $3}' candidate_lncRNA_classes.txt > lncRNA_genes
tail -n +2 lncRNA_genes > lncRNA_transcripts
rm lncRNA_genes
grep -w -F -f lncRNA_transcripts annotated_ensembl.gtf > merged.fixed.lncRNAs.gtf
grep --invert-match -F -f lncRNA_transcripts annotated_ensembl.gtf > merged.fixed.coding.gtf
rm annotated_ensembl.g*
sed -i 's/StringTie/lncRNA/' merged.fixed.lncRNAs.gtf
sed -i 's/StringTie/coding/' merged.fixed.coding.gtf
cat merged.fixed.coding.gtf merged.fixed.lncRNAs.gtf > final_annotated_ensembl.gtf
agat_sp_ensembl_output_style.pl -g final_annotated_ensembl.gtf -o final_annotated_ensembl.gff
echo ""
printf "${PURPLE}::: Parsing is done. The transcripts were classified and added to final_annotated_ensembl.gtf file...\n"
echo ""
rm merged.fixed.gff merged.fixed.gtf merged.fixed.lncRNAs.gtf merged.fixed.coding.gtf
echo ""
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 14. Obtaining Transcripts in FASTA format with gffread :::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
echo ""
gffread -w ensembl_transcripts.fa -g ${4}.fa final_annotated_ensembl.gtf
echo ""
printf "${PURPLE}::: Done. transcripts.fa are located in current directory\n"
echo ""
echo ""
printf "${PURPLE}::: Moving gffcompare results to gffcompare_outputs_ensembl folder ...\n"
echo ""
mkdir gffcompare_outputs_ensembl
mv *.loci *.stats *.refmap *.tmap *.tracking ./gffcompare_outputs_ensembl
echo ""
echo "Done"
echo ""
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 14. Performing gene annotation by using GAWN pipeline :::\n"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
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
cp ${4}.fa /${dir1}/gawn/03_data/genome.fasta
cp ensembl_transcripts.fa /${dir1}/gawn/03_data/transcriptome.fasta
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
printf "${YELLOW}::: 15. Extracting GO terms for each transcript :::\n"
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
printf "${YELLOW}::: 16. Predicting gene models from transcripts with AUGUSTUS (gff3 format) :::\n"
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
./augustus.2.5.5/src/augustus --species=human --progress=true --UTR=off --uniqueGeneId=true --gff3=on ensembl_transcripts.fa > augustus.gff3
echo ""
printf "${PURPLE}::: Done. augustus.gff3 file is present in current directory...${CYAN}\n"
echo ""
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 17. Converting gff3 to GTF format, collecting coding sequences and proteins with gffread...\n"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
gffread augustus.gff3 -T -o coding_transcript.gtf
gffread -x cds.fa -g ensembl_transcripts.fa coding_transcript.gtf
gffread -y prot.fa -g ensembl_transcripts.fa coding_transcript.gtf
# Re-formatting
cat cds.fa |rev|cut -d"." -f1 --complement|rev > transcripts_CDS.fa
cat prot.fa |rev|cut -d"." -f1 --complement|rev > transcripts_proteins.fa
rm cds.fa prot.fa
##########################################
# Re-formatting coding_transcripts.gtf
##########################################
sed 's/.t1"/"/' coding_transcript.gtf > coding_transcripts.gtf
echo ""
printf "${PURPLE}::: Done. AUGUSTUS predicted transcripts were summarized in coding_transcripts.gtf file located in current directory :::${CYAN}\n"
echo ""
rm coding_transcript.gtf 
#############################
# Configuring Summary Results
#############################
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 18. Configuring Summary Results :::\n"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::${CYAN}\n"
############################################
# Moving results to merged_annotation folder
############################################
echo ""
printf "${PURPLE}::: Moving results to output_files_ensembl folder :::${CYAN}\n"
mkdir output_files_ensembl
mv candidate_lncRNA_classes.txt final_annotated_ensembl.gtf final_annotated_ensembl.gff Stats.txt ensembl_transcripts.fa ensembl_aligned.gtf transcriptsGO.tab genesGO.tab transcripts_CDS.fa transcripts_proteins.fa coding_transcripts.gtf augustus.gff3 logfile_ensembl ./output_files_ensembl
cp /${dir1}/gawn/05_results/transcriptome_annotation_table.tsv /${dir1}/output_files_ensembl/
rm transcripts.fa.fai namelist* isoforms_per_gene_concatenated.tab lncRNA_transcripts merged.gtf merged_with_reference.gtf UCSC_compare* non_STRG.gtf STRG.gtf refGene.txt transcriptome_annotation_table.tsv 
echo ""
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n" 
echo "All Done. The transcripts were classified in ./output_files_ensembl"
echo ""
echo "Transcript discoveries are summarized in Stats.txt file located in ./output_files_ensembl . GAWN annotation is named transcriptome_annotation_table.tsv"
echo ""
echo "GTF file final_annotated_ensembl.gtf (standard GTF) and correspondent gff file (final_annotated_ensembl.gff) are located in ./output_files_ensembl. These files contains the annotated lncRNA/coding GTF in the second field".
echo ""
echo "candidate_lncRNA_classes.txt contained detailed long non-coding classification of transcripts".
echo ""
echo "Associated FASTA file to this GTF, named ensembl_transcripts.fa is located in ./output_files_ensembl"
echo ""
echo "AUGUSTUS GTF file suitable for transcript count quantification is named coding_transcripts.gtf. This GTF file contains all coding transcripts resolved by AUGUSTUS and is located in ./output_files_ensembl"
echo ""
echo "Associated Transcript coding sequences (transcripts_CDS.fa) and correspondent protein sequences (transcripts_proteins.fa) with coding_transcripts.gtf are located in ./output_files_ensembl"
echo ""
echo "GO terms associated to each transcript (and gene), named transcriptsGO.tab and genesGO.tab are located in ./output_files_ensembl"
echo ""
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${NC}\n"
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed
#
} | tee logfile_ensembl
#
