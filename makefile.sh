#!/bin/bash

dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
mkdir swissprot
cd swissprot
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz
gunzip swissprot.tar.gz
tar -xvf swissprot.tar
SWISSPROT_PATH=$PWD
echo "$SWISSPROT_PATH"
cd ..
cd test
echo 'SWISSPROT_DB="'$SWISSPROT_PATH'/swissprot"' >> gawn_config.sh
echo '#' >> gawn_config.sh
cd ..
mkdir bin
mkdir genome_1
mkdir get_transcripts
cp ./test/gawn_config.sh ./genome_1/
cp ./test/gawn_config.sh ./bash_scripts/
git clone https://github.com/neurobin/shc.git
cd shc/
./autogen.sh
./configure
make
cd ..
echo ""
echo "make done. Continue with install"
#
