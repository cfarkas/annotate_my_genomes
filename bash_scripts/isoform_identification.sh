#!/bin/bash

set -e

{

usage="$(basename "$0") [-h] [-m <gffcompare tmap file>] [-t <transcripts.fa file>] [-g <genome name>]
This script will produce an annotated csv table of transcripts, by using the tmap output file from the pipeline.
Arguments:
    -h  show this help text
    -m  gffcompare tmap output file. In example: UCSC_compare.stringtie.gtf.tmap or NCBI_compare.stringtie.gtf.tmap
    -t  transcripts file, output of annotate-my-genomes (transcripts.fa) or add-ncbi-annotation (NCBI_transcripts.fa) programs
    -g  UCSC genome name. In example: mm10, galGal6, hg38, rn6."
options=':hm:t:g:'
while getopts $options option; do
  case "$option" in
    h) echo "$usage"; exit;;
    m) m=$OPTARG;;
    t) t=$OPTARG;;
    g) g=$OPTARG;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
  esac
done

# mandatory arguments
if [ ! "$m" ] || [ ! "$t" ] || [ ! "$g" ]; then
  echo ""
  echo "arguments -m, -t and -g must be provided"
  echo ""
  echo "$usage" >&2; exit 1
fi

# Conditions : Input existance

if [ ! -e "$m" ]; then
    echo ""
    echo "$m does not exist. Check your -m input"
    echo ""
    exit 9999 # die with error code 9999
fi

if [ ! -e "$t" ]; then
    echo ""
    echo "$t does not exist. Check your -t input"
    echo ""
    exit 9999 # die with error code 9999
fi

# Conditions : Getting absolute path of inputs
echo ""
m_DIR="$( cd "$( dirname "$m" )" && pwd )"
echo ""
echo "::: The absolute path of -m is $m_DIR"
echo ""
t_DIR="$( cd "$( dirname "$t" )" && pwd )"
echo ""
echo "::: The absolute path of -t is $t_DIR"
echo ""

printf "${YELLOW}::: Defining Inputs :::\n"
echo ""
echo "Defining variables:"
echo""
FILE1="$m"
basename "$FILE1"
tmap_input="$(basename -- $FILE1)"
echo "The tmap file used as input is the following: $tmap_input"
echo ""
echo""
FILE2="$t"
basename "$FILE2"
transcripts_input="$(basename -- $FILE2)"
echo "The transcript file used as input is the following: $transcripts_input"

if [ -f ncbiRefSeqLink.txt ]; then
    echo "ncbiRefSeqLink.txt file found. Continue:"
    echo ""
    :
else
    echo "Downloading ncbiRefSeqLink.txt file"
    wget http://hgdownload.cse.ucsc.edu/goldenpath/${g}/database/ncbiRefSeqLink.txt.gz
    gunzip ncbiRefSeqLink.txt.gz
fi


if [ -f gffcompare_parser.py ]; then
    echo "gffcompare_parser.py script found. Continue:"
    echo ""
    :
else
    echo "Downloading gffcompare_parser.py script"
    wget https://raw.githubusercontent.com/cfarkas/annotate_my_genomes/master/additional_scripts/gffcompare_parser.py
    chmod 755 gffcompare_parser.py
fi

# Execute seqkit to convert fasta file to tabular file
seqkit fx2tab ${t_DIR}/${transcripts_input} > transcripts_Isoform.tab
# Execute gffcompare_parser.py
python gffcompare_parser.py --tmap_file ${m_DIR}/${tmap_input} --transcripts_file transcripts_Isoform.tab --genome ${g}
rm transcripts_Isoform.tab
echo "All done."
}
