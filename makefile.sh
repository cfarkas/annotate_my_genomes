#!/bin/bash

# make
dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
mkdir bin
mkdir genome_1
mkdir get_transcripts
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
./shc/src/shc -f ./bash_scripts/add_ncbi_annotation.sh -o ./add-ncbi-annotation
./shc/src/shc -f ./bash_scripts/IsoSeq_annotation.sh -o ./annotate-IsoSeq
mv annotate-my-genomes get-transcripts genome-download add-ncbi-annotation annotate-IsoSeq ./bin/
cp ./bin/annotate-my-genomes ./test/
cp ./bin/annotate-my-genomes ./genome_1/
cp ./bin/genome-download ./test/
cp ./bin/genome-download ./genome_1/
cp ./bin/add-ncbi-annotation ./test/
cp ./bin/add-ncbi-annotation ./genome_1/
cp ./bin/annotate-IsoSeq ./test/
cp ./bin/annotate-IsoSeq ./genome_1/
cp ./bin/get-transcripts ./get_transcripts/
cp ./bin/genome-download ./get_transcripts/
echo ""
echo "All done. Binaries are located in ./bin/ ./genome_1/ and ./test/ folders"
#
