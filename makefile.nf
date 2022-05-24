#!/usr/bin/env nextflow

/*
 * pipeline input parameters
 */

params.workdir = './'
params.conda = './environment.yml'

println """\
                     M A K E F I L E
         =======================================
         working_directory  :  ${params.workdir}
         environment        :  ${params.conda}
         """
         .stripIndent()


process make_and_install {
  echo true
  stageInMode 'copy'
  conda "${params.conda}"

  shell:
  '''
  cd "!{params.workdir}"
  dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
  rm -r -f swissprot
  mkdir swissprot
  cd swissprot
  wget --ignore-length ftp://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz
  gunzip swissprot.tar.gz
  tar -xvf swissprot.tar
  SWISSPROT_PATH=$PWD
  echo "$SWISSPROT_PATH"
  cd "!{params.workdir}"/test/
  exec 3<> gawn_config.sh
  #!/bin/bash >> gawn_config.sh
  echo "" >> gawn_config.sh
  # Modify the following parameter values according to your experiment >> gawn_config.sh
  # Do not modify the parameter names or remove parameters >> gawn_config.sh
  # Do not add spaces around the equal (=) sign >> gawn_config.sh
  echo "" >> gawn_config.sh
  # Global parameters >> gawn_config.sh
  NCPUS=10                    # Number of CPUs to use for analyses (int, 1+) >> gawn_config.sh
  echo "" >> gawn_config.sh
  # Genome indexing >> gawn_config.sh
  SKIP_GENOME_INDEXING=1      # 1 to skip genome indexing, 0 to index it >> gawn_config.sh
  echo "" >> gawn_config.sh
  # Genome annotation with transcriptome >> gawn_config.sh
  # NOTE: do not use compressed fasta files >> gawn_config.sh
  GENOME_NAME="genome.fasta"                  # Name of genome fasta file found in 03_data >> gawn_config.sh
  TRANSCRIPTOME_NAME="transcriptome.fasta"    # Name of transcriptome fasta file found in 03_data >> gawn_config.sh
  echo "" >> gawn_config.sh
  # Path to swissprot database >> gawn_config.sh
  echo 'SWISSPROT_DB="'$SWISSPROT_PATH'/swissprot"' >> gawn_config.sh
  echo '#' >> gawn_config.sh
  exec 3>&-
  cd "!{params.workdir}"
  mkdir bin
  mkdir genome_1
  mkdir get_transcripts
  cp "!{params.workdir}"/test/gawn_config.sh "!{params.workdir}"/genome_1/
  cp "!{params.workdir}"/test/gawn_config.sh "!{params.workdir}"/bash_scripts/
  cp "!{params.workdir}"/test/gawn_config.sh "!{params.workdir}"
  cd "!{params.workdir}"
  git clone https://github.com/cfarkas/shc.git
  cd shc/
  ./autogen.sh
  ./configure
  make
  # Install
  cd "!{params.workdir}"
  "!{params.workdir}"/shc/src/shc -f "!{params.workdir}"/bash_scripts/annotate_my_genomes.sh -o ./annotate-my-genomes
  "!{params.workdir}"/shc/src/shc -f "!{params.workdir}"/bash_scripts/annotate_my_genomes.sh -o ./annotate-my-genomes
  "!{params.workdir}"/shc/src/shc -f "!{params.workdir}"/bash_scripts/get_transcripts.sh -o ./get-transcripts
  "!{params.workdir}"/shc/src/shc -f "!{params.workdir}"/bash_scripts/genome_download.sh -o ./genome-download
  "!{params.workdir}"/shc/src/shc -f "!{params.workdir}"/bash_scripts/genome_download_macOSX.sh -o ./genome-download-macOSX
  "!{params.workdir}"/shc/src/shc -f "!{params.workdir}"/bash_scripts/add_ncbi_annotation.sh -o ./add-ncbi-annotation
  "!{params.workdir}"/shc/src/shc -f "!{params.workdir}"/bash_scripts/isoform_identification.sh -o ./isoform-identification
  mv "!{params.workdir}"/annotate-my-genomes "!{params.workdir}"/get-transcripts "!{params.workdir}"/genome-download "!{params.workdir}"/genome-download-macOSX "!{params.workdir}"/add-ncbi-annotation "!{params.workdir}"/isoform-identification "!{params.workdir}"/bin/
  cp "!{params.workdir}"/bin/annotate-my-genomes "!{params.workdir}"/test/
  cp "!{params.workdir}"/bin/annotate-my-genomes "!{params.workdir}"/genome_1/
  cp "!{params.workdir}"/bin/genome-download "!{params.workdir}"/test/
  cp "!{params.workdir}"/bin/genome-download "!{params.workdir}"/genome_1/
  cp "!{params.workdir}"/bin/genome-download-macOSX "!{params.workdir}"/test/
  cp "!{params.workdir}"/bin/genome-download-macOSX "!{params.workdir}"/genome_1/
  cp "!{params.workdir}"/bin/add-ncbi-annotation "!{params.workdir}"/test/
  cp "!{params.workdir}"/bin/add-ncbi-annotation "!{params.workdir}"/genome_1/
  cp "!{params.workdir}"/bin/isoform-identification "!{params.workdir}"/test/
  cp "!{params.workdir}"/bin/isoform-identification "!{params.workdir}"/genome_1/
  cp "!{params.workdir}"/bin/get-transcripts "!{params.workdir}"/get_transcripts/
  cp "!{params.workdir}"/bin/genome-download "!{params.workdir}"/get_transcripts/
  echo ""
  echo "::: All done. Binaries are located in "!{params.workdir}"/bin/ folder. :::"
  echo ""
  echo "::: With sudo privileges, users can do : sudo cp ./bin/* /usr/local/bin/ :::"
  echo ""
  echo ""
  '''
}
