#!/bin/bash

set -e

{

dir1=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
usage="$(basename "$0") [-h] [-a <stringtie.gtf>] [-n <NCBI_reference.gtf>] [-r <reference_genome.gtf>] [-g <reference_genome.fasta>] [-t <threads>]
This pipeline will Overlap StringTie transcripts (GTF format) with current NCBI annotation and will annotate novel transcripts.
Arguments:
    -h  show this help text
    -a  StringTie GTF
    -n  NCBI gene annotation (in GTF format)
    -r  UCSC gene annotation (in GTF format)
    -g  Reference genome (in fasta format)
    -t  Number of threads for processing (integer)"
options=':ha:n:r:g:t:'
while getopts $options option; do
  case "$option" in
    h) echo "$usage"; exit;;
    a) a=$OPTARG;;
    n) n=$OPTARG;;
    r) r=$OPTARG;;
    g) g=$OPTARG;;
    t) t=$OPTARG;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
  esac
done

# mandatory arguments
if [ ! "$a" ] || [ ! "$n" ] || [ ! "$r" ] || [ ! "$g" ] || [ ! "$t" ]; then
  echo "arguments -a, -n, -r, -g and -t must be provided"
  echo "$usage" >&2; exit 1
fi

begin=`date +%s`
#    .---------- constant part!
#    vvvv vvvv-- the code from above
YELLOW='\033[1;33m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color
echo "Cleaning directory..."
rm -r -f FEELnc gawn gff3sort gffcompare_outputs_NCBI 
echo ""

printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 1. Overlapping StringTie transcripts with NCBI annotation :::\n"
printf "${YELLOW}:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"
echo ""

gffcompare -R -r ${n} -s ${g} -o UCSC_compare ${a}
printf "${PURPLE}Done\n"
echo ""

printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 2. Writting novel discoveries to Stats.txt :::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"

echo ""
# Stats
exec 3<> Stats.txt
echo "Number of assembled genes:" >> Stats.txt
cat UCSC_compare.${a}.tmap | sed "1d" | cut -f4 | sort | uniq | wc -l >> Stats.txt
echo "" >> Stats.txt
echo "Number of novel genes:" >> Stats.txt
cat UCSC_compare.${a}.tmap | awk '$3=="u"{print $0}' | cut -f4 | sort | uniq | wc -l >> Stats.txt
echo "" >> Stats.txt
echo "Number of novel transcripts:" >> Stats.txt
cat UCSC_compare.${a}.tmap | awk '$3=="u"{print $0}' | cut -f5 | sort | uniq | wc -l >> Stats.txt
echo "" >> Stats.txt
echo "Number of transcripts matching annotation:" >> Stats.txt
cat UCSC_compare.${a}.tmap | awk '$3=="="{print $0}' | cut -f5 | sort | uniq | wc -l >> Stats.txt
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
awk '{print $4"\t"$1}' UCSC_compare.${a}.tmap > UCSC_compare.${a}.tmap.1
tail -n +2 UCSC_compare.${a}.tmap.1 > UCSC_compare.${a}.tmap.2
awk '$2 != "-"' UCSC_compare.${a}.tmap.2 > namelist
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
cat ${a} | parallel --pipe -j ${t} sed -f sed.script > final_annotated.gtf
rm -f fileA fileB *tmap.1 *tmap.2
# sorting GTF file
git clone https://github.com/cfarkas/gff3sort.git
perl ./gff3sort/gff3sort.pl final_annotated.gtf > final_annotated.sorted.gtf
rm final_annotated.gtf
mv final_annotated.sorted.gtf final_annotated.gtf
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
gffread -w NCBI_transcripts.fa -g ${g} final_annotated.gtf
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
cp ${g} /${dir1}/gawn/03_data/genome.fasta
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
cp /${dir1}/gawn/04_annotation/transcriptome.swissprot /${dir1}/
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

grep "NM_" ${r} > NM_coding.gtf
echo ""
printf "${PURPLE}::: 1/3) Filtering transcripts :::${CYAN}\n"
# Filter
FEELnc_filter.pl -i final_annotated.gtf -a NM_coding.gtf -b transcript_biotype=protein_coding > candidate_lncRNA.gtf
rm -r -f ${g}.index
printf "${PURPLE}::: 2/3) Evaluating coding potential :::${CYAN}\n"
# Coding_Potential
FEELnc_codpot.pl -i candidate_lncRNA.gtf -a NM_coding.gtf -b transcript_biotype=protein_coding -g ${g} --mode=shuffle
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
cd /${dir1}/
echo ""

printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 9. Predicting coding regions from transcripts with coding potential using TransDecoder :::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"

echo ""
echo ""
gffread -w coding-transcripts.fa -g ${g} coding-genes.gtf
TransDecoder.LongOrfs -m 60 -t coding-transcripts.fa
TransDecoder.Predict -t coding-transcripts.fa --single_best_only
awk '{print $1}' coding-transcripts.fa.transdecoder.bed > coding.sequences
grep "STRG." coding.sequences > coding.hits && rm coding.sequences
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
# removing transcript id by expansion
sed -i 's/[.][0-9]"/"/'g coding_transcripts.gtf
sed -i 's/[.][0-9][0-9]"/"/'g coding_transcripts.gtf
sed -i 's/[.][0-9][0-9][0-9]"/"/'g coding_transcripts.gtf
# removing protein id by expansion
sed -i 's/[.]p[0-9]//'g coding_transcripts.gtf
sed -i 's/[.]p[0-9][0-9]//'g coding_transcripts.gtf
sed -i 's/[.]p[0-9][0-9][0-9]//'g coding_transcripts.gtf
sed -i 's/[.]p[0-9][0-9][0-9][0-9]//'g coding_transcripts.gtf
sed -i 's/[.]p[0-9][0-9][0-9][0-9][0-9]//'g coding_transcripts.gtf
cat coding_transcripts.gtf | parallel --pipe -j ${t} sed -f sed.script > coding_transcripts.fixed.gtf
rm coding_transcripts.gtf
mv coding_transcripts.fixed.gtf coding_transcripts.gtf
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
rm coding-genes.gtf lncRNAs.gtf other-genes.gtf sed.script transcriptome.hits
### Novel coding genes and correspondent proteins
echo ""
echo "::: Obtaining novel coding transcripts (cds) and correspondent proteins"
echo ""
perl -lne 'print "@m" if @m=(/((?:transcript_id|gene_id)\s+\S+)/g);' coding_transcripts.gtf > final_annotated.tab
sed -i 's/transcript_id //g' final_annotated.tab
sed -i 's/;/\t/g' final_annotated.tab
sed -i 's/gene_id//g' final_annotated.tab
sed -i 's/"//g' final_annotated.tab
awk '!a[$0]++' final_annotated.tab > transcripts_and_genes.tab && rm final_annotated.tab
awk '{print $2"\t"$1}' transcripts_and_genes.tab > coding-genes-and-transcripts.tab && rm transcripts_and_genes.tab
awk '$1 ~ /STRG./' coding-genes-and-transcripts.tab > novel-coding-genes.matches
awk '{print $2}' novel-coding-genes.matches > novel-coding-transcripts.matches
seqkit fx2tab cds.fa > cds.tab
seqkit fx2tab prot.fa > prot.tab
grep -w -F -f novel-coding-transcripts.matches cds.tab > novel-coding-cds.tab
grep -w -F -f novel-coding-transcripts.matches prot.tab > novel-coding-prot.tab
seqkit tab2fx novel-coding-cds.tab > novel-cds.fa && seqkit tab2fx novel-coding-prot.tab > novel-prot.fa
rm novel-coding-cds.tab novel-coding-prot.tab novel-coding-transcripts.matches novel-coding-genes.matches coding-genes-and-transcripts.tab cds.tab prot.tab
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

printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
printf "${YELLOW}::: 11. Moving results to output_files_NCBI folder :::\n"
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::${CYAN}\n"

echo ""
printf "${PURPLE}::: Moving results to output_files_NCBI folder :::${CYAN}\n"
mkdir output_files_NCBI
mv candidate_lncRNA_classes.txt final_annotated.gtf final_annotated.gff NCBI_transcripts.fa cds.fa prot.fa Stats.txt coding_transcripts.gtf transcriptome.swissprot logfile novel-cds.fa novel-prot.fa ./output_files_NCBI/
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
echo "TransDecoder GTF file suitable to parse NCBI_transcripts.fa (coding_transcripts.gtf), contains all coding transcripts resolved by TransDecoder and is located in ./output_files_NCBI"
echo ""
echo "Predicted coding sequences (cds.fa) and correspondent protein sequences (prot.fa) are located in ./output_files_NCBI"
echo ""
echo "Novel predicted coding sequences (novel-cds.fa) and correspondent protein sequences (novel-prot.fa) are located in ./output_files_NCBI"
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
