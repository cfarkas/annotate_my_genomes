#!/bin/bash

dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
merged_with_reference=${1}
reference_genome=${2}
gene_name=${3}

if [ "$1" == "-h" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {merged_with_reference} {reference_genome} {gene_name}"
  echo ""
  echo "This script will align all transcripts coming from a given gene"
  echo ""
  echo "{merged_with_reference}: Name of the StringTie annotated GTF from the pipeline"
  echo ""
  echo "{reference_genome}: Genome Assembly (fasta format)"
  echo ""
  echo "{gene_name}: gene_name to process"
  echo ""
  exit 0
fi

if [ "$1" == "-help" ]; then
    echo ""
  echo "Usage: bash ./`basename $0` {merged_with_reference} {reference_genome} {gene_name}"
  echo ""
  echo "This script will align all transcripts coming from a given gene"
  echo ""
  echo "{merged_with_reference}: Name of the StringTie annotated GTF from the pipeline"
  echo ""
  echo "{reference_genome}: Genome Assembly (fasta format)"
  echo ""
  echo "{gene_name}: gene_name to process"
  echo ""
  exit 0
fi
if [ "$1" == "--h" ]; then
   echo ""
  echo "Usage: bash ./`basename $0` {merged_with_reference} {reference_genome} {gene_name}"
  echo ""
  echo "This script will align all transcripts coming from a given gene"
  echo ""
  echo "{merged_with_reference}: Name of the StringTie annotated GTF from the pipeline"
  echo ""
  echo "{reference_genome}: Genome Assembly (fasta format)"
  echo ""
  echo "{gene_name}: gene_name to process"
  echo ""
  exit 0
fi
if [ "$1" == "--help" ]; then
    echo ""
  echo "Usage: bash ./`basename $0` {merged_with_reference} {reference_genome} {gene_name}"
  echo ""
  echo "This script will align all transcripts coming from a given gene"
  echo ""
  echo "{merged_with_reference}: Name of the StringTie annotated GTF from the pipeline"
  echo ""
  echo "{reference_genome}: Genome Assembly (fasta format)"
  echo ""
  echo "{gene_name}: gene_name to process"
  echo ""
  exit 0
fi

[ $# -eq 0 ] && { echo "Usage: bash ./`basename $0` {merged_with_reference} {reference_genome} {gene_name}"; exit 1; }

if [ $# -ne 3 ]; then
  echo 1>&2 "Usage: bash ./`basename $0` {merged_with_reference} {reference_genome} {gene_name}"
  exit 3
fi
echo ""
echo "Obtaining GTF for the given gene_name using merged_with_reference.gtf file"
grep "${3}" merged_with_reference.gtf > ${3}.gtf
echo "Done."
echo ""
echo "Obtaining gene-associated transcripts in fasta format"
gffread -w ${3}.fa -g ${2} ${3}.gtf
echo "Done."
echo ""