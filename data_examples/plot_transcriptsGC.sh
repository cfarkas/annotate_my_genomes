#!/bin/bash

GTF=${1}
Reference_genome=${2}

if [ "$1" == "-h" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {GTF} {Reference_genome}"
  echo ""
  echo "This script will plot GC content of every transcript annotated in a given GTF file using reference genome in fasta format"
  echo ""
  echo "GTF: GTF file containing transcripts to plot (such as final_annotated.gtf)"
  echo ""
  echo "Reference_genome: Reference genome in fasta format (e.g.: galGal6.fa)."
  exit 0
fi

if [ "$1" == "--help" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {GTF} {Reference_genome}"
  echo ""
  echo "This script will plot GC content of every transcript annotated in a given GTF file using reference genome in fasta format"
  echo ""
  echo "GTF: GTF file containing transcripts to plot (such as final_annotated.gtf)"
  echo ""
  echo "The GTF genome must be also indexed. e.g.: samtools faidx hg19.fa"
  echo ""
  echo "Reference_genome: Reference genome in fasta format (e.g.: galGal6.fa)."
  exit 0
fi
if [ "$1" == "--h" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {GTF} {Reference_genome}"
  echo ""
  echo "This script will plot GC content of every transcript annotated in a given GTF file using reference genome in fasta format"
  echo ""
  echo "GTF: GTF file containing transcripts to plot (such as final_annotated.gtf)"
  echo ""
  echo "Reference_genome: Reference genome in fasta format (e.g.: galGal6.fa)."
  exit 0
fi

if [ "$1" == "-help" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {GTF} {Reference_genome}"
  echo ""
  echo "This script will plot GC content of every transcript annotated in a given GTF file using reference genome in fasta format"
  echo ""
  echo "GTF: GTF file containing transcripts to plot (such as final_annotated.gtf)"
  echo ""
  echo "The GTF genome must be also indexed. e.g.: samtools faidx hg19.fa"
  echo ""
  echo "Reference_genome: Reference genome in fasta format (e.g.: galGal6.fa)."
  exit 0
fi

[ $# -eq 0 ] && { echo "Usage: bash ./`basename $0` {GTF} {Reference_genome}"; exit 1; }

if [ $# -ne 2 ]; then
  echo 1>&2 "Usage: bash ./`basename $0` {GTF} {Reference_genome}"
  exit 3
fi

grep "\<transcript\>" ${1} > transcript_seqs.gtf
cat transcript_seqs.gtf | awk '{print $1"\t"$4"\t"$5}' > gtf.bed
cat gtf.bed | tr '[ ]' '[\t]' > gtf.tab.bed
rm transcript_seqs.gtf gtf.bed
bedtools nuc -fi ${2} -bed gtf.tab.bed > GC_nuc
Rscript plot_transcriptGC.R
rm gtf.tab.bed
echo ""
echo "All done. transcript_GC_content.pdf is located in current directory. GC_nuc file contain all AT/GC content per transcript"
echo ""
