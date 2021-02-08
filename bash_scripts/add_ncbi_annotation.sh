#!/bin/bash

set -e

{

dir1=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
stringtie_gtf=${1}
NCBI_reference_genome_gtf=${2}
reference_genome_gtf=${3}
reference_genome_fasta=${4}
threads=${5}

if [ "$1" == "-h" ]; then
  echo ""
  echo "Usage: ./`basename $0` [stringtie_gtf] [NCBI_reference_genome_gtf] [reference_genome_gtf] [reference_genome_fasta] [threads]"
  echo ""
  echo "This pipeline will Overlap StringTie transcripts (GTF format) with current NCBI annotation and will annotate novel transcripts"
  echo ""
  echo "[stringtie_gtf]: StringTie GTF file"
  echo ""
  echo "[NCBI_reference_genome_gtf]: NCBI reference GTF file."
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
  echo "Usage: ./`basename $0` [stringtie_gtf] [NCBI_reference_genome_gtf] [reference_genome_gtf] [reference_genome_fasta] [threads]"
  echo ""
  echo "This pipeline will Overlap StringTie transcripts (GTF format) with current NCBI annotation and will annotate novel transcripts"
  echo ""
  echo "[stringtie_gtf]: StringTie GTF file"
  echo ""
  echo "[NCBI_reference_genome_gtf]: NCBI reference GTF file."
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
  echo "Usage: ./`basename $0` [stringtie_gtf] [NCBI_reference_genome_gtf] [reference_genome_gtf] [reference_genome_fasta] [threads]"
  echo ""
  echo "This pipeline will Overlap StringTie transcripts (GTF format) with current NCBI annotation and will annotate novel transcripts"
  echo ""
  echo "[stringtie_gtf]: StringTie GTF file"
  echo ""
  echo "[NCBI_reference_genome_gtf]: NCBI reference GTF file."
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
  echo "Usage: ./`basename $0` [stringtie_gtf] [NCBI_reference_genome_gtf] [reference_genome_gtf] [reference_genome_fasta] [threads]"
  echo ""
  echo "This pipeline will Overlap StringTie transcripts (GTF format) with current NCBI annotation and will annotate novel transcripts"
  echo ""
  echo "[stringtie_gtf]: StringTie GTF file"
  echo ""
  echo "[NCBI_reference_genome_gtf]: NCBI reference GTF file."
  echo ""
  echo "[reference_genome_gtf]: UCSC reference GTF file."
  echo ""
  echo "[reference_genome_fasta]: Current UCSC assembly genome in fasta format."
  echo ""
  echo "[threads]: Number of threads for parallel text processing (Integer)"
  echo ""
  exit 0
fi

[ $# -eq 0 ] && { echo "Usage: ./`basename $0` [stringtie_gtf] [NCBI_reference_genome_gtf] [reference_genome_gtf] [reference_genome_fasta] [threads]"; exit 1; }

if [ $# -ne 5 ]; then
  echo 1>&2 "Usage: ./`basename $0` [stringtie_gtf] [NCBI_reference_genome_gtf] [reference_genome_gtf] [reference_genome_fasta] [threads]"
  exit 3
fi

begin=`date +%s`
#    .---------- constant part!
#    vvvv vvvv-- the code from above
YELLOW='\033[1;33m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color
echo "Cleaning directory..."
rm -r -f augustus.* FEELnc gawn 
echo ""
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 1. Overlapping StringTie transcripts with NCBI annotation :::\n"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
echo ""
gffcompare -R -r ${2} -s ${4} -o UCSC_compare ${1}
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
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 3. Replacing gene_id field in final_annotated.gtf file with NCBI gene_id's :::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
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
cat ${1} | parallel --pipe -j ${5} sed -f sed.script > final_annotated.gtf
rm -f fileA fileB *tmap.1 *tmap.2
printf "${PURPLE}::: Done. Gene_id field was replaced in the StringTie.gtf file and final_annotated.gtf was generated with these changes\n"
echo ""
printf "${PURPLE}::: Moving gffcompare results to gffcompare_outputs folder ...\n"
echo ""
mkdir gffcompare_outputs_NCBI
mv *.loci *.stats *.refmap *.tmap *.tracking ./gffcompare_outputs_NCBI
echo ""
printf "${PURPLE}::: Continue with protein-coding annotation\n" 
echo ""
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 4. Obtaining Transcripts in FASTA format with gffread :::\n"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
echo ""
gffread -w NCBI_transcripts.fa -g ${4} final_annotated.gtf
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
git clone https://github.com/cfarkas/gawn.git
cd gawn/02_infos/
dir2=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
echo "Done"
echo ""
cd /${dir1}/
cp ${4} /${dir1}/gawn/03_data/genome.fasta
cp NCBI_transcripts.fa /${dir1}/gawn/03_data/transcriptome.fasta
rm /${dir2}/gawn_config.sh
cp gawn_config.sh /${dir2}/gawn_config.sh
echo ""
printf "${PURPLE}::: Starting GAWN transcript annotation${CYAN}\n"
echo ""
cd /${dir1}/gawn/
./gawn 02_infos/gawn_config.sh
echo ""
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
cd /${dir1}/
cp /${dir1}/gawn/04_annotation/transcriptome.hits /${dir1}/
printf "${PURPLE}::: Done. transcriptome hits were succesfully extracted :::${CYAN}\n"
echo ""
############################################
# FEELnc long noncoding RNA identification #
############################################
cd /${dir1}/
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 7. Classifying protein-coding and long non-coding transcripts with FEELnc :::\n"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
### Cloning FEELnc in current directory
git clone https://github.com/cfarkas/FEELnc.git
echo ""
cp ${4} ${3} final_annotated.gtf /${dir1}/FEELnc/
cd FEELnc
grep "NM_" ${3} > NM_coding.gtf
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
printf "${PURPLE}::: FEELnc Test done. Continue with final_annotated.gtf file :::\n"
echo ""
cd ..
### Running FEELnc
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 8. Running FEELnc on final_annotated.gtf file :::\n"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
echo ""
# Filter
FEELnc_filter.pl -i final_annotated.gtf -a NM_coding.gtf -b transcript_biotype=protein_coding > candidate_lncRNA.gtf
# Coding_Potential
FEELnc_codpot.pl -i candidate_lncRNA.gtf -a NM_coding.gtf -b transcript_biotype=protein_coding -g ${4} --mode=shuffle
# Classifier
FEELnc_classifier.pl -i feelnc_codpot_out/candidate_lncRNA.gtf.lncRNA.gtf -a NM_coding.gtf > candidate_lncRNA_classes.txt
echo ""
printf "${PURPLE}::: FEELnc calculations were done. The output is called candidate_lncRNA_classes.txt:::\n"
echo ""
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 9. Parsing GAWN and FEELnc outputs :::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
echo ""
cp candidate_lncRNA_classes.txt /${dir1}/
cd /${dir1}/
awk '{print $3}' candidate_lncRNA_classes.txt > lncRNA_genes
tail -n +2 lncRNA_genes > lncRNA_transcripts
rm lncRNA_genes
grep -w -F -f lncRNA_transcripts final_annotated.gtf > merged.fixed.lncRNAs.gtf
grep --invert-match -F -f lncRNA_transcripts final_annotated.gtf > merged.fixed.coding.gtf
rm final_annotated.gtf
sed -i 's/StringTie/lncRNA/' merged.fixed.lncRNAs.gtf
awk '{print $1"\t"$2}' transcriptome.hits > coding_list
awk -F'\t' '$2!=""' coding_list > coding_list.tab
tail -n +2 "coding_list.tab" > coding_transcripts
awk '{print $1}' coding_transcripts > coding_transcripts.tab
rm coding_lis* coding_transcripts lncRNA_transcripts
grep -w -F -f coding_transcripts.tab merged.fixed.coding.gtf > coding-genes.gtf
grep --invert-match -F -f coding_transcripts.tab merged.fixed.coding.gtf > other-genes.gtf
sed -i 's/StringTie/coding/' coding-genes.gtf
cat coding-genes.gtf merged.fixed.lncRNAs.gtf other-genes.gtf > final_annotated.gtf
gffread -E -F --merge final_annotated.gtf -o final_annotated.gff
rm coding_transcripts.tab
echo "All done"
echo ""
##########################################
# Gene Prediction Step with TransDecoder #
##########################################
cd /${dir1}/
echo ""
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 10. Predicting coding regions from transcripts with coding potential using TransDecoder :::\n"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
echo ""
echo ""
gffread -w coding-transcripts.fa -g ${4} coding-genes.gtf
TransDecoder.LongOrfs -m 60 -t coding-transcripts.fa
TransDecoder.Predict -t coding-transcripts.fa --single_best_only
awk '{print $1}' coding-transcripts.fa.transdecoder.bed > coding.sequences
grep "STRG." coding.sequences > coding.hits && rm coding.sequences
echo ""
printf "${PURPLE}::: Done. coding-transcripts.fa.transdecoder.gff3 file is present in current directory...${CYAN}\n"
echo ""
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 11. Converting gff3 to GTF format and formatting coding sequences and proteins :::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
gffread coding-transcripts.fa.transdecoder.gff3 -T -o coding_transcripts.gtf
sed -i 's/GENE.//'g coding_transcripts.gtf
sed 's/.[0-9]~/~/'g coding-transcripts.fa.transdecoder.cds > cds.fa
sed -i 's/GENE./gene_id="/'g cds.fa
sed -i 's/~~/" cds_id=/'g cds.fa
sed 's/.[0-9]~/~/'g coding-transcripts.fa.transdecoder.pep > prot.fa
sed -i 's/GENE./gene_id="/'g prot.fa
sed -i 's/~~/" cds_id=/'g prot.fa
sed 's/.[0-9]~/~/'g coding-transcripts.fa.transdecoder.bed > cds.bed
sed -i 's/GENE./gene_id="/'g cds.bed
sed -i 's/~~/";cds_id=/'g cds.bed
cat cds.fa | parallel --pipe -j ${5} sed -f sed.script > cds.fixed.fa
cat cds.bed | parallel --pipe -j ${5} sed -f sed.script > cds.fixed.bed
cat prot.fa | parallel --pipe -j ${5} sed -f sed.script > prot.fixed.fa
rm cds.fa cds.bed prot.fa
mv cds.fixed.fa cds.fa
mv cds.fixed.bed cds.bed 
mv prot.fixed.fa prot.fa
echo ""
###############################
# Configuring Summary Results #
###############################
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 12. Moving results to output_files_NCBI folder :::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
echo ""
printf "${PURPLE}::: Moving results to output_files_NCBI folder :::${CYAN}\n"
mkdir output_files_NCBI
mv candidate_lncRNA_classes.txt final_annotated.gtf final_annotated.gff NCBI_transcripts.fa cds.fa prot.fa cds.bed Stats.txt coding_transcripts.gtf logfile ./output_files_NCBI/
cp /${dir1}/gawn/04_annotation/transcriptome.hits /${dir1}/output_files_NCBI/
echo ""
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
echo "All Done. The transcripts were classified in ./output_files_NCBI"
echo ""
echo "Transcript discoveries are summarized in Stats.txt file located in ./output_files_NCBI. GAWN protein annotation is named transcriptome.hits"
echo ""
echo "GTF file named final_annotated.gtf (and correspondent gff file) are located in ./output_files_NCBI, containing novel genes and lncRNA classification (second field in GTF file)"
echo ""
echo "candidate_lncRNA_classes.txt contained detailed long non-coding classification of transcripts"
echo ""
echo "Associated FASTA file to this GTF (NCBI_transcripts.fa) is located in ./output_files_NCBI"
echo ""
echo "TransDecoder GTF file suitable for transcript count quantification (coding_transcripts.gtf) contains all coding transcripts resolved by TransDecoder and is located in ./output_files_NCBI"
echo ""
echo "Associated Transcript coding sequences (cds.fa) and correspondent protein sequences (prot.fa) are located in ./output_files_NCBI"
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
