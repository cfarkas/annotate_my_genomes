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
cp gawn_config.sh ${dir}/genome_1
cp gawn_config.sh ${dir}/bash_scripts
cd ..
git clone https://github.com/neurobin/shc.git
cd shc/
./autogen.sh
./configure
make
cd /${dir1}/
./shc/src/shc -f ./bash_scripts/annotate_my_genomes.sh -o annotate_my_genomes
./shc/src/shc -f ./bash_scripts/get_transcripts.sh -o get_transcripts
mkdir bin
mkdir genome_1
mv annotate_my_genomes get_transcripts ./bin/
cp ./bin/annotate_my_genomes ./test/
cp ./bin/annotate_my_genomes ./genome_1/
echo "bin folder containing executable binaries are made"
echo ""
echo "make done"
#
