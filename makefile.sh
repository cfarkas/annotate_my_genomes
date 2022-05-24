#!/bin/bash

dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
rm -r -f swissprot
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
cp ./test/gawn_config.sh ./
git clone https://github.com/cfarkas/shc.git
cd shc/
./autogen.sh
./configure
make
cd ..
echo ""
echo "make done. Continue with install"
# Install
./shc/src/shc -f ./bash_scripts/annotate_my_genomes.sh -o ./annotate-my-genomes
./shc/src/shc -f ./bash_scripts/get_transcripts.sh -o ./get-transcripts
./shc/src/shc -f ./bash_scripts/genome_download.sh -o ./genome-download
./shc/src/shc -f ./bash_scripts/genome_download_macOSX.sh -o ./genome-download-macOSX
./shc/src/shc -f ./bash_scripts/add_ncbi_annotation.sh -o ./add-ncbi-annotation
./shc/src/shc -f ./bash_scripts/isoform_identification.sh -o ./isoform-identification
mv annotate-my-genomes get-transcripts genome-download genome-download-macOSX add-ncbi-annotation isoform-identification ./bin/
cp ./bin/annotate-my-genomes ./test/
cp ./bin/annotate-my-genomes ./genome_1/
cp ./bin/genome-download ./test/
cp ./bin/genome-download ./genome_1/
cp ./bin/genome-download-macOSX ./test/
cp ./bin/genome-download-macOSX ./genome_1/
cp ./bin/add-ncbi-annotation ./test/
cp ./bin/add-ncbi-annotation ./genome_1/
cp ./bin/isoform-identification ./test/
cp ./bin/isoform-identification ./genome_1/
cp ./bin/get-transcripts ./get_transcripts/
cp ./bin/genome-download ./get_transcripts/
echo "::: All done. Binaries are located in ./bin/ folder. :::" 
echo ""
echo "::: With sudo privileges, users can do : sudo cp ./bin/* /usr/local/bin/ :::"
echo ""
#
