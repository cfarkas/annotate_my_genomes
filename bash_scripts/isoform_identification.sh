#!/bin/bash

set -e

{

usage="$(basename "$0") [-h] [-m <gffcompare tmap file>] [-t <transcripts.fa file>] [-g <genome name>]
This script will produce an annotated csv table of transcripts, by using the tmap output file from add-ncbi-annotation pipeline.
Arguments:
    -h  show this help text
    -m  NCBI gffcompare tmap output file. As example: NCBI_compare.stringtie.gtf.tmap
    -t  transcripts file, output of add-ncbi-annotation program. As example: NCBI_transcripts.fa
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
echo ""
printf "${YELLOW}::: Defining Inputs :::\n"
echo  ""
FILE1="$m"
basename "$FILE1"
tmap_input="$(basename -- $FILE1)"
echo "The tmap file used as input is the following: $tmap_input"
echo""
FILE2="$t"
basename "$FILE2"
transcripts_input="$(basename -- $FILE2)"
echo "The transcript file used as input is the following: $transcripts_input"
echo ""
if [ -f ncbiRefSeqLink.txt ]; then
    echo "::: ncbiRefSeqLink.txt file found. Continue:"
    echo ""
    :
else
    echo "::: Downloading ncbiRefSeqLink.txt file"
    wget http://hgdownload.cse.ucsc.edu/goldenpath/${g}/database/ncbiRefSeqLink.txt.gz
    gunzip ncbiRefSeqLink.txt.gz
fi

# Inputs for python
cp ${m_DIR}/${tmap_input} ./stringtie_for_script.tmap
seqkit fx2tab ${t_DIR}/${transcripts_input} > transcripts_Isoform.tab

# Execute gffcompare_parser.py
python << END

import sys
import pandas as pd

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    OKRED = '\033[91m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

df = pd.read_csv('stringtie_for_script.tmap', sep = '\t')
print(df.sample(10))
print("Total number of transcripts:", df.shape[0])
print("")

df2 = df[~df.ref_id.astype(str).str.contains('-')]
novel_transcripts = df[df.ref_id.astype(str).str.contains('-')]

df3 = df2[["ref_gene_id", "ref_id", "class_code", "qry_gene_id", "qry_id", "num_exons", "FPKM", "TPM"]]
df_novel_transcripts = novel_transcripts[["ref_gene_id", "ref_id", "class_code", "qry_gene_id", "qry_id", "num_exons", "FPKM", "TPM"]]

print("Reference transcripts:")
print(df3.sample(10))
print("")

print("Novel transcripts:")
print(df_novel_transcripts.sample(10))
print("")

colnames=['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18']
dfA1 = pd.read_csv('ncbiRefSeqLink.txt', sep = '\t', low_memory=False, names=colnames, header=None)
print(dfA1.head(10))

dfA2 = dfA1[['0', '1', '2', '3', '5', '14', '16']]

dfA2 = dfA2.rename(columns={'0': 'ref_id', '1': 'Annotation Status', '2' : 'NCBI RefSeq Gene ID', '3' : 'Transcript Description', '5' : 'NCBI RefSeq Protein ID', '14' : 'Alternative Gene Name', '16' : 'RefSeq Transcript Info'})
print("ncbiRefSeqLink annotation:")
print(dfA2.sample(10))
print("")

colnames = ['qry_id', 'cds_seq', 'none']
cds = pd.read_csv('transcripts_Isoform.tab', sep = '\t', names=colnames)
cds2 = cds[["qry_id", "cds_seq"]]
print("transcripts file:")
print(cds2.sample(10))
print("")

result1 = pd.merge(df3, dfA2, on='ref_id', how='inner')
result1.sample(10)
result2 = pd.merge(result1, cds2, on='qry_id', how='inner')
result2.sample(10)
result3 = pd.merge(df_novel_transcripts, cds2, on='qry_id', how='inner')
result3.sample(10)
print("Number of Joined Transcripts (reference):", result2.shape[0])
print("")
print("Number of Joined Transcripts (novel):", result3.shape[0])
print("")
result2.to_csv('Ref_Transcript_Annotation.csv', index=False)
result3.to_csv('Novel_Transcript_Annotation.csv', index=False)
print(bcolors.OKGREEN + "::: Done. Ref_Transcript_Annotation.csv and Novel_Transcript_Annotation.csv were succesfully produced" + bcolors.ENDC)
print("")
END

rm -r -f transcripts_Isoform.tab stringtie_for_script.tmap ncbiRefSeqLink.txt
echo "::: All done. :::"
}
