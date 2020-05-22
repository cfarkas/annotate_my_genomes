#!/bin/bash

dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
merged_with_reference=${1}
reference_genome=${2}
gene_name=${3}

if [ "$1" == "-h" ]; then
  echo ""
  echo "Usage: ./`basename $0` [final_annotated] [reference_genome] [gene_name]"
  echo ""
  echo "This program will obtain and align all transcripts coming from a given gene, in order to obtain a consensus"
  echo ""
  echo "[final_annotated]: Name of the StringTie annotated GTF from the pipeline"
  echo ""
  echo "[reference_genome]: Masked Ensenbl genome Assembly (fasta format)"
  echo ""
  echo "[gene_name]: gene_name to process"
  echo ""
  exit 0
fi

if [ "$1" == "-help" ]; then
  echo ""
  echo "Usage: ./`basename $0` [merged_with_reference] [reference_genome] [gene_name]"
  echo ""
  echo "This script will obtain and align all transcripts coming from a given gene, in order to obtain a consensus"
  echo ""
  echo "[final_annotated]: Name of the StringTie annotated GTF from the pipeline"
  echo ""
  echo "[reference_genome]: Masked Ensenbl genome Assembly (fasta format)"
  echo ""
  echo "[gene_name]: gene_name to process"
  echo ""
  exit 0
fi

if [ "$1" == "--h" ]; then
  echo ""
  echo "Usage: ./`basename $0` [merged_with_reference] [reference_genome] [gene_name]"
  echo ""
  echo "This script will obtain and align all transcripts coming from a given gene, in order to obtain a consensus"
  echo ""
  echo "[final_annotated]: Name of the StringTie annotated GTF from the pipeline"
  echo ""
  echo "[reference_genome]: Masked Ensenbl genome Assembly (fasta format)"
  echo ""
  echo "[gene_name]: gene_name to process"
  echo ""
  exit 0
fi

if [ "$1" == "--help" ]; then
  echo ""
  echo "Usage: ./`basename $0` [merged_with_reference] [reference_genome] [gene_name]"
  echo ""
  echo "This script will obtain and align all transcripts coming from a given gene, in order to obtain a consensus"
  echo ""
  echo "[final_annotated]: Name of the StringTie annotated GTF from the pipeline"
  echo ""
  echo "[reference_genome]: Masked Ensenbl genome Assembly (fasta format)"
  echo ""
  echo "[gene_name]: gene_name to process"
  echo ""
  exit 0
fi

[ $# -eq 0 ] && { echo "Usage: ./`basename $0` [final_annotated] [reference_genome] [gene_name]"; exit 1; }

if [ $# -ne 3 ]; then
  echo 1>&2 "Usage: ./`basename $0` [final_annotated] [reference_genome] [gene_name]"
  exit 3
fi
echo "Working in $dir"
echo ""
echo "Obtaining GTF for the given gene_name using merged_with_reference.gtf file"
grep "\<${3}\>" ${1} > ${3}.gtf
echo "Done."
echo ""
echo "Obtaining gene-associated transcripts in fasta format"
gffread -w ${3}.fa -g ${2} ${3}.gtf
echo "Done."
echo ""
echo "Aligning transcript sequences with Clustal Omega"
clustalo -i ${3}.fa -o ${3}.aln
echo "Done"
echo ""
echo "Obtaining consensus sequence from alignment with EMBOSS consensus"
em_cons -sequence ${3}.aln -outseq ${3}.cons
echo ""
echo "All done. ${3}.cons file contain suitable sequence to be validated for qPCR."
echo ""
echo "${3}.fa contain transcript sequences in fasta format associated with ${3} gene"
