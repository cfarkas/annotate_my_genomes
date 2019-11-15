#!/bin/bash

genome=${1}

if [ "$1" == "-h" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {genome}"
  echo ""
  echo "This script will download and index the specified genome from the uscs server (goldenpath) including annotation"
  echo ""
  echo "{genome}: Name of the genome assembly"
  echo ""
  exit 0
fi

if [ "$1" == "-help" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {genome}"
  echo ""
  echo "This script will download and index the specified genome from the uscs server (goldenpath) including annotation"
  echo ""
  echo "{genome}: Name of the genome assembly"
  echo ""
  exit 0
fi
if [ "$1" == "--h" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {genome}"
  echo ""
  echo "This script will download and index the specified genome from the uscs server (goldenpath) including annotation"
  echo ""
  echo "{genome}: Name of the genome assembly"
  echo ""
  exit 0
fi

if [ "$1" == "--help" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {genome}"
  echo ""
  echo "This script will download and index the specified genome from the uscs server (goldenpath) including annotation"
  echo ""
  echo "{genome}: Name of the genome assembly"
  echo ""
  exit 0
fi

[ $# -eq 0 ] && { echo "Usage: bash ./`basename $0` {genome}"; exit 1; }

if [ $# -ne 1 ]; then
  echo 1>&2 "Usage: bash ./`basename $0`  {genome}"
  exit 3
fi

if [ -f {genome}.fa ]; then
    echo "{genome}.fa is already in current folder. Stop."
    exit 3
    : 
else
    echo "Continue..."
    # Obtaining {genome}.fa genome, and indexing
    wget http://hgdownload.cse.ucsc.edu/goldenpath/${genome}/bigZips/${genome}.2bit
fi
if [ -f twoBitToFa ]; then
    echo "twoBitToFa script found. Continue:"
    echo ""
    : 
else
    echo "Downloading twoBitToFa script"
    wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa 
fi

chmod 755 twoBitToFa
./twoBitToFa ${genome}.2bit ${genome}.fa
samtools faidx ${genome}.fa
echo ""
echo "All Done. ${genome}.fa is located in the current directory"
