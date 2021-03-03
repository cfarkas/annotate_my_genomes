#!/bin/bash

set -e

dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
usage="$(basename "$0") [-h] [-f <final_annotated.gtf>] [-g <reference_genome.fasta>] [-i <gene_name>]
This program will obtain and align all transcripts coming from a given gene, in order to obtain a consensus.
Arguments:
    -h  show this help text
    -f  Name of the StringTie annotated GTF from the pipeline
    -g  Reference genome (in fasta format)
    -i  Gene Symbol"
options=':hf:g:i:'
while getopts $options option; do
  case "$option" in
    h) echo "$usage"; exit;;
    f) f=$OPTARG;;
    g) g=$OPTARG;;
    i) i=$OPTARG;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
  esac
done

# mandatory arguments
if [ ! "$f" ] || [ ! "$g" ] || [ ! "$i" ]; then
  echo "arguments -f, -g and -i must be provided"
  echo "$usage" >&2; exit 1
fi

echo "Working in $dir"
echo ""
echo "Obtaining GTF for the given gene_name using final_annotated.gtf file"
grep "\<${i}\>" ${f} > ${i}.gtf
echo "Done."
echo ""
echo "Obtaining gene-associated transcripts in fasta format"
gffread -w ${i}.fa -g ${g} ${i}.gtf
echo "Done."
echo ""
echo "Aligning transcript sequences with Clustal Omega"
clustalo -i ${i}.fa -o ${i}.aln
echo "Done"
echo ""
echo "Obtaining consensus sequence from alignment with EMBOSS consensus"
em_cons -sequence ${i}.aln -outseq ${i}.cons
echo ""
echo "All done. ${i}.cons file contain suitable sequence to be validated for qPCR."
echo ""
echo "${i}.fa contain transcript sequences in fasta format associated with ${i} gene"
