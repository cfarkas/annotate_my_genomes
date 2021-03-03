
#!/bin/bash

set -e

usage="$(basename "$0") [-h] [-f <final_annotated.gtf>] [-g <reference_genome.fasta>]
This script will obtain metrics from the annotated StringTie transcripts (final_annotated.gtf) and output them into -transcriptome_metrics- sudirectory.
Arguments:
    -h  show this help text
    -f  Name of the StringTie annotated GTF from the pipeline
    -g  Reference genome (in fasta format)
options=':hf:g:'
while getopts $options option; do
  case "$option" in
    h) echo "$usage"; exit;;
    f) f=$OPTARG;;
    g) g=$OPTARG;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
  esac
done

# mandatory arguments
if [ ! "$f" ] || [ ! "$g" ]; then
  echo "arguments -f and -g must be provided"
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
rm -r -f transcriptome_metrics
echo ""
echo "done"
echo ""
echo "===> Working on transcriptome metrics"
echo ""
perl -lne 'print "@m" if @m=(/((?:transcript_id|gene_id)\s+\S+)/g);' ${f} > final_annotated.tab
sed -i 's/transcript_id //g' final_annotated.tab
sed -i 's/;/\t/g' final_annotated.tab
sed -i 's/gene_id//g' final_annotated.tab
sed -i 's/"//g' final_annotated.tab
awk '!a[$0]++' final_annotated.tab > genes_and_transcripts.tab && rm final_annotated.tab
awk '{print $1"\t"$2}' genes_and_transcripts.tab > genes-and-transcripts.tab && rm genes_and_transcripts.tab
awk '{print $1}' genes-and-transcripts.tab > genes.tab
# Novel genes list
grep "STRG." genes.tab > novel-genes.tab
# Known genes list
grep -v "STRG." genes.tab > known-genes.tab
echo "::: Parsing final_annotated.gtf file to obtain novel/known and coding/lncRNA transcripts, respectively."
echo ""
grep -w -F -f novel-genes.tab final_annotated.gtf > novel-genes.gtf
grep -w -F -f known-genes.tab final_annotated.gtf > known-genes.gtf
grep "coding" known-genes.gtf > known-genes-coding.gtf
grep "lncRNA" known-genes.gtf > known-genes-lncRNA.gtf
grep "StringTie" known-genes.gtf > known-genes-other.gtf  # other = no lncRNA and no protein-coding
grep "coding" novel-genes.gtf > novel-genes-coding.gtf
grep "lncRNA" novel-genes.gtf > novel-genes-lncRNA.gtf
grep "StringTie" novel-genes.gtf > novel-genes-other.gtf  # other = no lncRNA and no protein-coding
echo "::: We will use gffread to obtain reconciled and novel transcripts in the parsed GTF file"
echo ""
gffread -w known-transcripts-coding.fa -g ${g} known-genes-coding.gtf
gffread -w known-transcripts-lncRNA.fa -g ${g} known-genes-lncRNA.gtf
gffread -w known-transcripts-other.fa -g ${g} known-genes-other.gtf
gffread -w novel-transcripts-coding.fa -g ${g} novel-genes-coding.gtf
gffread -w novel-transcripts-lncRNA.fa -g ${g} novel-genes-lncRNA.gtf
gffread -w novel-transcripts-other.fa -g ${g} novel-genes-other.gtf
exec 3<> transcriptome_metrics.txt
echo "Number of reconciled coding transcripts:" >> transcriptome_metrics.txt
grep ">" known-transcripts-coding.fa -c >> transcriptome_metrics.txt
echo "" >> transcriptome_metrics.txt
echo "Number of reconciled non-coding transcripts:" >> transcriptome_metrics.txt
grep ">" known-transcripts-lncRNA.fa -c >> transcriptome_metrics.txt
echo "" >> transcriptome_metrics.txt
echo "Number of other expressed features, annotated:" >> transcriptome_metrics.txt
grep ">" known-transcripts-other.fa -c >> transcriptome_metrics.txt
echo "" >> transcriptome_metrics.txt
echo "Number of non-annotated coding transcripts:" >> transcriptome_metrics.txt
grep ">" novel-transcripts-coding.fa -c >> transcriptome_metrics.txt
echo "" >> transcriptome_metrics.txt
echo "Number of non-annotated non-coding novel transcripts:" >> transcriptome_metrics.txt
grep ">" novel-transcripts-lncRNA.fa -c >> transcriptome_metrics.txt
echo "" >> transcriptome_metrics.txt
echo "Number of novel expressed features:" >> transcriptome_metrics.txt
grep ">" novel-transcripts-other.fa -c >> transcriptome_metrics.txt
exec 3>&-
echo "::: Done. transcriptome_metrics.txt contains metrics of classified transcripts. Continue with gene metrics.."
echo ""
echo ""
echo "===> Working on gene metrics"
echo ""
exec 3<> gene_metrics.txt
# known coding genes counts
perl -lne 'print "@m" if @m=(/((?:transcript_id|gene_id)\s+\S+)/g);' known-genes-coding.gtf > known-genes-coding.tab
sed -i 's/transcript_id //g' known-genes-coding.tab
sed -i 's/;/\t/g' known-genes-coding.tab
sed -i 's/gene_id//g' known-genes-coding.tab
sed -i 's/"//g' known-genes-coding.tab
awk '{print $1}' known-genes-coding.tab > known-genes-coding.tabular && rm known-genes-coding.tab
awk '!a[$0]++' known-genes-coding.tabular > known-genes-coding.tab && rm known-genes-coding.tabular
echo "Number of reconciled coding genes:" >> gene_metrics.txt
cat known-genes-coding.tab | wc -l >> gene_metrics.txt
echo "" >> gene_metrics.txt
echo "Number of reconciled non-coding genes:" >> gene_metrics.txt
perl -lne 'print "@m" if @m=(/((?:transcript_id|gene_id)\s+\S+)/g);' known-genes-lncRNA.gtf > known-genes-lncRNA.tab
sed -i 's/transcript_id //g' known-genes-lncRNA.tab
sed -i 's/;/\t/g' known-genes-lncRNA.tab
sed -i 's/gene_id//g' known-genes-lncRNA.tab
sed -i 's/"//g' known-genes-lncRNA.tab
awk '{print $1}' known-genes-lncRNA.tab > known-genes-lncRNA.tabular && rm known-genes-lncRNA.tab
awk '!a[$0]++' known-genes-lncRNA.tabular > known-genes-lncRNA.tab && rm known-genes-lncRNA.tabular
cat known-genes-lncRNA.tab | wc -l >> gene_metrics.txt
echo "" >> gene_metrics.txt
echo "Number of other reconciled genes:" >> gene_metrics.txt
perl -lne 'print "@m" if @m=(/((?:transcript_id|gene_id)\s+\S+)/g);' known-genes-other.gtf > known-genes-other.tab
sed -i 's/transcript_id //g' known-genes-other.tab
sed -i 's/;/\t/g' known-genes-other.tab
sed -i 's/gene_id//g' known-genes-other.tab
sed -i 's/"//g' known-genes-other.tab
awk '{print $1}' known-genes-other.tab > known-genes-other.tabular && rm known-genes-other.tab
awk '!a[$0]++' known-genes-other.tabular > known-genes-other.tab && rm known-genes-other.tabular
cat known-genes-other.tab | wc -l >> gene_metrics.txt
echo "" >> gene_metrics.txt
echo "Number of non-annotated coding genes:" >> gene_metrics.txt
perl -lne 'print "@m" if @m=(/((?:transcript_id|gene_id)\s+\S+)/g);' novel-genes-coding.gtf > novel-genes-coding.tab
sed -i 's/transcript_id //g' novel-genes-coding.tab
sed -i 's/;/\t/g' novel-genes-coding.tab
sed -i 's/gene_id//g' novel-genes-coding.tab
sed -i 's/"//g' novel-genes-coding.tab
awk '{print $1}' novel-genes-coding.tab > novel-genes-coding.tabular && rm novel-genes-coding.tab
awk '!a[$0]++' novel-genes-coding.tabular > novel-genes-coding.tab && rm novel-genes-coding.tabular
cat novel-genes-coding.tab | wc -l >> gene_metrics.txt
echo "" >> gene_metrics.txt
echo "Number of non-annotated non-coding genes:" >> gene_metrics.txt
perl -lne 'print "@m" if @m=(/((?:transcript_id|gene_id)\s+\S+)/g);' novel-genes-lncRNA.gtf > novel-genes-lncRNA.tab
sed -i 's/transcript_id //g' novel-genes-lncRNA.tab
sed -i 's/;/\t/g' novel-genes-lncRNA.tab
sed -i 's/gene_id//g' novel-genes-lncRNA.tab
sed -i 's/"//g' novel-genes-lncRNA.tab
awk '{print $1}' novel-genes-lncRNA.tab > novel-genes-lncRNA.tabular && rm novel-genes-lncRNA.tab
awk '!a[$0]++' novel-genes-lncRNA.tabular > novel-genes-lncRNA.tab && rm novel-genes-lncRNA.tabular
cat novel-genes-lncRNA.tab | wc -l >> gene_metrics.txt
echo "" >> gene_metrics.txt
echo "Number of non-annotated other genes:" >> gene_metrics.txt
perl -lne 'print "@m" if @m=(/((?:transcript_id|gene_id)\s+\S+)/g);' novel-genes-other.gtf > novel-genes-other.tab
sed -i 's/transcript_id //g' novel-genes-other.tab
sed -i 's/;/\t/g' novel-genes-other.tab
sed -i 's/gene_id//g' novel-genes-other.tab
sed -i 's/"//g' novel-genes-other.tab
awk '{print $1}' novel-genes-other.tab > novel-genes-other.tabular && rm novel-genes-other.tab
awk '!a[$0]++' novel-genes-other.tabular > novel-genes-other.tab && rm novel-genes-other.tabular
cat novel-genes-other.tab | wc -l  >> gene_metrics.txt
exec 3>&-
echo "::: gene_metrics.txt were succesfully generated."
echo ""
mkdir transcriptome_metrics
mv genes-and-transcripts.tab genes.tab novel-genes.tab known-genes.tab novel-genes.gtf known-genes.gtf known-genes-coding.gtf known-genes-lncRNA.gtf known-genes-other.gtf novel-genes-coding.gtf novel-genes-lncRNA.gtf novel-genes-other.gtf known-transcripts-coding.fa known-transcripts-lncRNA.fa known-transcripts-other.fa novel-transcripts-coding.fa novel-transcripts-lncRNA.fa novel-transcripts-other.fa known-genes-coding.tab known-genes-lncRNA.tab known-genes-other.tab novel-genes-coding.tab novel-genes-lncRNA.tab novel-genes-other.tab ./transcriptome_metrics/
mv transcriptome_metrics.txt gene_metrics.txt ./transcriptome_metrics/
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${NC}\n"
echo "All Done. The transcripts were classified in the ./transcriptome_metrics subdirectory."
echo ""
echo "Transcript discoveries are summarized in transcriptome_metrics.txt file."
echo ""
echo "Gene discoveries are summarized in gene_metrics.txt file."
echo ""
echo "known-genes-coding.gtf, known-genes-lncRNA.gtf and known-genes-other.gtf contains reconciled annotation with reference, in GTF format"
echo ""
echo "novel-genes-coding.gtf, novel-genes-lncRNA.gtf and novel-genes-other.gtf contains novel annotation with reference, in GTF format"
echo ""
echo "known-transcripts-coding.fa, known-transcripts-lncRNA.fa and known-transcripts-other.fa contains reconciled classified transcripts, in FASTA format"
echo ""
echo "novel-transcripts-coding.fa, novel-transcripts-lncRNA.fa and novel-transcripts-other.fa contains novel classified transcripts, in FASTA format"
echo ""
printf "${YELLOW}::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::${NC}\n"
