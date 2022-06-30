ARG UBUNTU_VER=20.04
ARG CONDA_VER=latest
ARG OS_TYPE=x86_64
ARG PY_VER=3.8.11
ARG TF_VER=2.5.0

FROM ubuntu:${UBUNTU_VER}

# System packages 
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc
RUN apt-get update && apt-get install -yq build-essential g++ python-dev autotools-dev libicu-dev libbz2-dev libboost-all-dev zlib1g-dev curl wget unzip sed jq vim nano libidn11 libnet-perl perl-doc liblmdb-dev && apt-get install -y git && apt install -y make && apt install -y autoconf
RUN apt install -y parallel


# Install make
RUN apt update && apt install -y make && apt install -y autoconf

# cmake
RUN wget https://github.com/Kitware/CMake/releases/download/v3.17.3/cmake-3.17.3.tar.gz --no-check-certificate && tar -zxvf cmake-3.17.3.tar.gz && cd cmake-3.17.3 && apt-get install libssl-dev && ./bootstrap && make && make install && cd /

# R and dependences
RUN apt install -y dirmngr gnupg apt-transport-https ca-certificates software-properties-common && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' && apt install -y r-base
RUN R -e "install.packages('ROCR',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('randomForest',dependencies=TRUE, repos='http://cran.rstudio.com/')"

# stringtie
RUN git clone https://github.com/gpertea/stringtie && cd stringtie && make release && make test && ./run_tests.sh && cp stringtie /usr/local/bin/ && cp stringtie /usr/bin/ && cd SuperReads_RNA && ./install.sh && cd /

# FEELnc
RUN apt-get install -y libcurl4 libcurl4-openssl-dev && apt-get install -y libxml-dom-xpath-perl && apt-get install -y cpanminus
RUN cpanm Parallel::ForkManager Bio::DB::SeqFeature

# Installing KmerInShort 
RUN git clone --recursive https://github.com/rizkg/KmerInShort && cd KmerInShort && mkdir build;  cd build;  cmake ..;  make -j 8 && cp KmerInShort /usr/local/bin/ && cp KmerInShort /usr/bin/ && cd /

# Installing fasta_ushuffle
RUN wget -O fasta_ushuffle.zip https://github.com/agordon/fasta_ushuffle/archive/refs/heads/master.zip --no-check-certificate && unzip fasta_ushuffle.zip && cd fasta_ushuffle-master/ && make  && cp fasta_ushuffle ushuffle /usr/local/bin/ && cp fasta_ushuffle ushuffle /usr/bin/ && cd / 

# Installing FEELnc
RUN git clone https://github.com/tderrien/FEELnc.git && cd /FEELnc && export FEELNCPATH=$(pwd) && export PERL5LIB=$PERL5LIB:${FEELNCPATH}/lib/ && export PATH=$PATH:${FEELNCPATH}/scripts/ && export PATH=$PATH:${FEELNCPATH}/utils/ && export PATH=$PATH:${FEELNCPATH}/bin/LINUX/ && cp -r ${FEELNCPATH}/bin/LINUX/ ~/bin/ 
ENV PATH=/FEELnc/bin/LINUX:${PATH}
ENV FEELNCPATH=/FEELnc
ENV PERL5LIB=:/FEELnc/lib/
ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/FEELnc/bin/LINUX:/FEELnc/utils:/FEELnc/scripts/

# FEELnc Test
RUN cd /FEELnc/test/ && FEELnc_filter.pl -i transcript_chr38.gtf -a annotation_chr38.gtf -b transcript_biotype=protein_coding > candidate_lncRNA.gtf && FEELnc_codpot.pl -i candidate_lncRNA.gtf -a annotation_chr38.gtf -b transcript_biotype=protein_coding -g genome_chr38.fa --mode=shuffle && FEELnc_classifier.pl -i feelnc_codpot_out/candidate_lncRNA.gtf.lncRNA.gtf -a annotation_chr38.gtf > candidate_lncRNA_classes.txt && cd /

# Installing gffcompare and gclib
RUN git clone https://github.com/gpertea/gclib && git clone https://github.com/gpertea/gffcompare && git clone https://github.com/gpertea/gffread
RUN cd /gffcompare && make release && cp gffcompare trmap /usr/local/bin/ && cp gffcompare trmap /usr/bin/ && cd /
RUN cd /gffread && make release && cp gffread /usr/local/bin/ && cd /

# Installing ncbi-blast+
RUN apt-get remove -y ncbi-blast+
RUN apt-get install -y ncbi-blast+

# gmap 
RUN apt-get install -y gmap

# bedtools
RUN apt-get install -y bedtools

# samtools
RUN apt-get install -y samtools && apt-get install -y bcftools

# transdecoder (TransDecoder.LongOrfs 5.5.0)
RUN wget https://github.com/TransDecoder/TransDecoder/archive/refs/tags/TransDecoder-v5.5.0.zip --no-check-certificate && unzip TransDecoder-v5.5.0.zip && mv TransDecoder-TransDecoder-v5.5.0 TransDecoder-v5.5.0 && apt-get install -y hmmer
RUN cd /TransDecoder-v5.5.0 && ln -s /TransDecoder-v5.5.0/TransDecoder.LongOrfs /usr/local/bin/ && ln -s /TransDecoder-v5.5.0/TransDecoder.Predict /usr/local/bin/ && ln -s /TransDecoder-v5.5.0/TransDecoder.LongOrfs /usr/bin/ && ln -s /TransDecoder-v5.5.0/TransDecoder.Predict /usr/bin/ 

# seqkit
RUN wget https://github.com/shenwei356/seqkit/releases/download/v0.12.1/seqkit_linux_386.tar.gz --no-check-certificate && gunzip seqkit_linux_386.tar.gz && tar -xvf seqkit_linux_386.tar && cp seqkit /usr/local/bin/ && cp seqkit /usr/bin/ && cd /

# emboss
RUN apt-get install -y emboss

# Clustalo
RUN apt-get install -y clustalo

# Cufflinks
RUN apt-get install -y cufflinks

# gawk 
RUN apt install -y gawk

# minimap2
RUN apt-get -y install minimap2

# pandas
RUN apt-get -y install python3-pip
RUN pip install pandas
RUN pip install numpy

# annotate_my_genomes
RUN git clone https://github.com/cfarkas/annotate_my_genomes.git && cd annotate_my_genomes && chmod 755 ./makefile.sh && ./makefile.sh
RUN ln -s /annotate_my_genomes/bash_scripts/ /usr/local/bin/add_ncbi_annotation.sh
RUN ln -s /annotate_my_genomes/bash_scripts/ /usr/local/bin/annotate_my_genomes.sh 
RUN ln -s /annotate_my_genomes/bash_scripts/ /usr/local/bin/genome_download.sh
RUN ln -s /annotate_my_genomes/bash_scripts/ /usr/local/bin/get_transcripts.sh
RUN ln -s /annotate_my_genomes/bash_scripts/ /usr/bin/add_ncbi_annotation.sh
RUN ln -s /annotate_my_genomes/bash_scripts/ /usr/bin/annotate_my_genomes.sh
RUN ln -s /annotate_my_genomes/bash_scripts/ /usr/bin/genome_download.sh
RUN ln -s /annotate_my_genomes/bash_scripts/ /usr/bin/get_transcripts.sh 
ENV PATH=/annotate_my_genomes/bin:${PATH}
