#!/bin/bash

mkdir swissprot
cd swissprot
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz
gunzip swissprot.tar.gz
tar -xvf swissprot.tar
cd ..
gunzip CNIT.tar.gz
tar -xvf CNIT.tar
###
echo "make Done"

