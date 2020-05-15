# annotate_my_genomes 
### :microscope: :hatching_chick: :hatched_chick: 
Genome annotation pipeline using long sequencing reads from non-model (and model) vertebrate organisms.

# Pipeline Outline
  Often, genomes from non-model organisms (and even from model organisms) contain reference genome annotation available in GTF format (Gene Transfer Format), but these annotations may fail to capture all genome features. Novel genes and novel transcripts can be absent from reference genome annotations due tissue or stage-specific gene expression when using RNA-seq data for transcript characterization.
  
  annotate_my_genomes are a set of bash scripts that aim to annotate transfrags obtained by genome-guided transcriptome assembly strategies (StringTie) coming from long read RNA-Seq alignments in vertebrate genomes (i.e. PacBio/Oxford Nanopore technologies). Transcripts are classified by its coding potential, probable gene function and identified as novel or reconciliated with the current reference annotation from Ensembl. Also, coding sequences in nucleotides and correspondent proteins sequences can be reconstructed from these procedures. 
  
  The pipeline is designed for:
  
- Use as input GTF file coming from StringTie tool using long sequencing reads settings (for documentation, please see http://ccb.jhu.edu/software/stringtie/ and the documentation in this repository).
- Conciliate current gene annotation from an organism with this GTF and expand this annotation by annotating novel transcripts with GAWN (Genome Annotation Without Nightmares, please see https://github.com/enormandeau/gawn). A concilliated GTF file is generated with annotated gene names and corresponding StringTie assembled transfrags (transcripts). As example:
```
gene_id "YF5"; transcript_id "YF5.4"
(Genes denoted with "STRG" prefix are novel).
```
- The resulting GTF file is validated by using "Validate GTF" tool from Brent Lab: http://mblab.wustl.edu/software.html#validategtf
- Perform gene prediction on reconstructed transcripts with Augustus software. Please see (http://augustus.gobics.de/)
- Assess coding potential of each assembled transcript with FEELnc tool (https://github.com/tderrien/FEELnc).
- Assign to each transcripts and genes gene ontology terms (GO) and output formatted tables compatibles with WEGO annotation server: (http://wego.genomics.org.cn/). 

This pipeline requieres to run:
- StringTie assembled transcripts (in GTF format)
- USCS reference genome annotation (in GTF format)
- UCSC genome assembly (masked, fasta format)

The two last requirements can be downloaded by using genome-download script provided in this repository

# Dependences:

### gcc and g++ compilers, version >= 6 
Ubuntu/linux may come with GCC/G++ compilers <= version 6. To complile augustus gene prediction tool properly, users must upgrade old gcc/g++ compilers as follows: 
```
sudo apt-get update 
sudo apt-get install build-essential software-properties-common -y
sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y
## Important: if you have problems with this steps, users may do the following in order to fix it:
# sudo apt-get remove python3-apt
# sudo apt-get install python3-apt
sudo apt-get update
sudo apt-get install gcc-6 g++-6 -y
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-6 60 --slave /usr/bin/g++ g++ /usr/bin/g++-6

# To check new GCC and G++ versions, do the following: 
g++ --version
gcc --version
```

### Obtaining and installing StringTie (v2.0 release needed)

```
git clone https://github.com/gpertea/stringtie
cd stringtie
make release
sudo cp stringtie /usr/local/bin/
```
### Obtaining and installing gffcompare and gffread

```
git clone https://github.com/gpertea/gclib
git clone https://github.com/gpertea/gffcompare
git clone https://github.com/gpertea/gffread
cd gclib
make release
cd ..
cd gffcompare
make release
sudo cp gffcompare trmap /usr/local/bin/
cd ..
cd gffread
make release
sudo cp gffread /usr/local/bin/
cd ..
```

### Installing up-to-date ncbi-blast+ version (v2.9.0)
#### If you have ncbi-blast+ version > 2.7 is OK, older binaries also work, but GAWN pipeline strongly recommends > v2.7.1+

```
# To remove older ncbi-blast+ binaries from your system 
sudo apt-get remove ncbi-blast+

# Installing ncbi-blast+ version 2.9.0

wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.9.0+-x64-linux.tar.gz
tar -xzvf ncbi-blast-2.9.0+-x64-linux.tar.gz
rm ncbi-blast-2.9.0+-x64-linux.tar.gz
cd ncbi-blast-2.9.0+/bin/
sudo cp * /usr/local/bin/ 
```
### Installing GMAP genomic aligner program 
#### (for documentation, please see http://research-pub.gene.com/gmap/)

```
sudo apt-get install gmap
```

### Installing dependences for FEELnc : FlExible Extraction of LncRNA 
#### (for documentation, please see https://github.com/tderrien/FEELnc)

```
## Perl Requirements: Parallel and BioPerl (with sudo privileges)

sudo perl -MCPAN -e shell
install Parallel::ForkManager
install BioPerl
quit

# Installing KmerInShort 

git clone --recursive https://github.com/rizkg/KmerInShort
cd KmerInShort
mkdir build;  cd build;  cmake ..;  make -j 8
sudo cp KmerInShort /usr/local/bin/

# Installing fasta_ushuffle

git clone git://github.com/agordon/fasta_ushuffle.git
cd fasta_ushuffle
make
sudo cp ushuffle /usr/local/bin/
sudo cp fasta_ushuffle /usr/local/bin/
```

### Obtaining and installing up-to-date SAMtools with htslib (version >= 1.9)
(Old samtools version can also work). Users needs to install version up to date of these three packages. Users can first install htslib v1.9 and then samtools with bcftools v1.9, respectively. For downloading these packages, see http://www.htslib.org/download/). The latter can be accomplish by downloading the three packages, decompressing it, and doing the following:
```
wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2
bzip2 -d htslib-1.10.2.tar.bz2
tar -xvf htslib-1.10.2.tar
rm htslib-1.10.2.tar
cd htslib-1.10.2    # and similarly for samtools
sudo ./configure --prefix=/usr/local/bin
sudo make
sudo make install
# this step is only for samtools
sudo cp samtools /usr/local/bin/

# Similarly as htslib, samtools and bcftools can be downloaded as follows:

wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2
```

Then in a terminal type
>samtools

to check 1.10 version (using htslib v1.10)

### Obtaining and installing EMBOSS toolkit (Open Source software for molecular biology)
Complete instructions can be found at the webpage: https://ssbio.readthedocs.io/en/latest/instructions/emboss.html. A way to install it can be the following:
```
sudo apt-get install emboss
```
### Obtaining and installing Clustal Omega (DNA/Protein alignment program)
Description can be found at the webpage: http://www.clustal.org/omega/. A way to install it can be the following:
```
sudo apt-get install clustalo
```

### Installing AGAT: Another Gff Analysis Toolkit (AGAT). Suite of tools to handle gene annotations in any GTF/GFF format.
For more information, please visit: https://github.com/NBISweden/AGAT
```
# Install through bioconda recipe:
conda install agat

# Standard install
git clone https://github.com/NBISweden/AGAT.git # Clone AGAT
cd AGAT                                         # move into AGAT folder
perl Makefile.PL                                # Check all the dependencies*
make                                            # Compile
make test                                       # Test
make install                                    # Install
```

# Installation: 

Clone this repository every time you need to work with a different assembly and/or genome. In the folder you want to annotate the GTF file, do the following: 
  
```
git clone https://github.com/cfarkas/annotate_my_genomes.git
cd annotate_my_genomes
# make
bash makefile.sh
```
Binaries are located in bin, genome_1 and test folders, respectively.

### Obtaining StringTie GTF file for annotation

#### 1) Alignment of long sequencing reads using minimap aligner (e.g.: against galGal6 genome from UCSC, using 30 threads). You can use gmap as well. 
```
# Installing minimap2
git clone https://github.com/lh3/minimap2
cd minimap2 && make
sudo cp minimap2 /usr/local/bin/

# Convert PacBio subreads (bam files) to fastq
samtools bam2fq m54027_190807_082031.subreads.bam > reads.fastq

# Download Gallus gallus v6 fasta file (use genome-download program from this repository, located in ./bin folder)
./genome-download galGal6

# Aligning with minimap2
minimap2 -ax splice galGal6.fa reads.fastq > aln_galGal6.sam -t 30
samtools view -S -b aln_galGal6.sam -@ 30 > aln_galGal6.bam
samtools sort aln_galGal6.bam -@ 30 > aln_galGal6.sorted.bam
samtools index aln_galGal6.sorted.bam -@ 30
```

If users also sequenced with Illumina, short reads can be aligned by using hisat2 (https://ccb.jhu.edu/software/hisat2/manual.shtml) and merged with minimap alignments as follows:
```
# Install hisat2
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.0.4-Linux_x86_64.zip
unzip hisat2-2.0.4-Linux_x86_64.zip
sudo cp hisat2-2.0.4/hisat2* /usr/local/bin/

# build hisat2 index
hisat2-build galGal6.fa -p 40 ./galGal6

# align short illumina reads agains galGal6 genome from USCS, using 30 threads
hisat2 -x ./galGal6 -p 30 -1 41-A3_S1_R1_001.fastq.gz -2 41-A3_S1_R2_001.fastq.gz | samtools view -bS - > 41.bam

# merge with PacBio alignments (aln_galGal6.bam) with illumina alignments (41.bam), using 30 threads:
samtools merge PacBio_Illumina_merged.bam aln_galGal6.bam 41.bam -@ 30 -f

# Sort and Index
samtools sort -o PacBio_Illumina_merged.sorted.bam PacBio_Illumina_merged.bam -@ 30
samtools index PacBio_Illumina_merged.sorted.bam
```
PacBio_Illumina_merged.sorted.bam can be used as input for StringTie. 

#### 2) Obtaining GTF (transcripts.gtf) from the above alignment using StringTie (e.g.: using -p: 30 threads, -L: long read settings)
```
# Alignments from long reads (PacBio)
stringtie -p 1 -L -v -a 4 -o transcripts.gtf aln_galGal6.sorted.bam

# Alignments from long and short reads (PacBio + Illumina)
stringtie -p 1 -v -a 4 -o transcripts.gtf PacBio_Illumina_merged.sorted.bam

# If the above fails, users can increase -j and -c parameters (useful for large BAM file processing)
stringtie -p 1 -j 2 -c 2 -v -a 4 -o transcripts.gtf PacBio_Illumina_merged.sorted.bam
```

#

# Quickstart (Running the test)

1) Optionally, edit number of cpus in /test/gawn_config.sh:

- NCPUS=10
  - Increase this value to speed-up things :rocket:

2) Run the pipeline with a set of transcripts from chromosome 33, Gallus gallus genome version "6". Users need to specify the stringtie output (GTF format), UCSC reference genome (GTF annotation and fasta file) and the number of threads for text processing (5 for this example). Go to /annotate_my_genomes/test and do the following:

```
# Download Gallus gallus v6 fasta assembly with matched GTF file (Masked fasta file, from UCSC repository)
./genome-download galGal6

# Execute in folder
./annotate-my-genomes stringtie_chr33.gtf galGal6.gtf galGal6.fa 5
```
#
## Usage examples

#### To download reference genome files (Ensembl), please visit: https://uswest.ensembl.org/info/data/ftp/index.html
#### To download reference genome sequences (UCSC), use genome-download program from this repository. To check genome names, please visit: https://genome.ucsc.edu/cgi-bin/hgGateway

(Optional) Edit NCPUS value in gawn_config.sh file in "genome_1" folder. Default is 10

- For mouse assembly using "target.gtf" in genome_1 folder, using 30 threads for text processing:
```
./genome-download mm10
./annotate-my-genomes target.gtf mm10.gtf mm10.fa 30
```
- For rabbit assembly using "target.gtf" in genome_1 folder, using 30 threads for text processing:
```
./genome-download oryCun2
./annotate-my-genomes target.gtf oryCun2.gtf oryCun2.fa 30
```

#
## Downstream analysis using outputs:

### (1) Gene quantification procedure examples using output GTF file (merged.fixed.gtf):
- Install HTSeq-count: (please see https://htseq.readthedocs.io/en/release_0.11.1/index.html)
```
sudo apt-get install build-essential python2.7-dev python-numpy python-matplotlib python-pysam python-htseq
```

- Gene-level quantification using "merged.fixed.gtf" GTF file
```
htseq-count --stranded=no --format bam condition1.bam condition2.bam merged.fixed.gtf > gene_counts
```

### (2) "genesGO.tab" and "transcriptsGO.tab" output files usage:
- A gene list in tabular format (e.g.: coming from differential expression analysis from genes or transcripts) can be intersected with these output files by using the galaxy framework (https://usegalaxy.org/). A gene list can be the following:

```
less gene_list.tab

STRG.14047
GMPPB
STRG.26267
MSX2
STRG.12490
HMGA1
...
```

can be intersected with "genesGO.tab" output file, containing annotated GO terms. 
```
less genesGO.tab

ARF5    GO:0005737 GO:0005794 GO:0005886 GO:0005525 GO:0006886 GO:0016192
GCC1    GO:0005829 GO:0005794 GO:0000139 GO:0005886
MIR1651 GO:0005829 GO:0043231 GO:0019903 GO:0043666
LMF2    GO:0005737 GO:0016604 GO:0048471 GO:0001691 GO:0017112 GO:0001558 GO:0043087 GO:0007283
...
```
Following these simple steps in galaxy: https://usegalaxy.org/u/carlosfarkas/h/joining-user-genelist-with-genegotab-file

The file number 4, "Cut on data 3" file can be inputted in WEGO 2.0 server (http://wego.genomics.org.cn/) in native format, to explore and plot GO terms. 

- IMPORTANT: if you copied the gene list from excel and you paste it in a txt file, you must convert this file to tabular and remove non ASCII characters. Galaxy (https://usegalaxy.org/) can easily do this by uploading the .txt file and executing Text Manipulation --> Convert delimiters to TAB in a single step.   

### (3) I need the transcript sequences matching each gene. Also validate conserved regions with qPCR. What can I do?:

The above gene list in tabular format can also be used to extract: 
- Transcripts sequences associated to each gene. 
- Align transcript sequences in order to obtain consensus sequences

This can be accomplished by copying merged.fixed.gtf file and an user-provided gene list in tabular format (such as gene_list.tab) to get_transcripts folder/ . Execute the following (e.g.: for chicken genome):

```
# Downloading masked Gallus gallus v6 genome
./genome-download galGal6

# generate "commands" file using provided list of genes 
awk '{print "./get-transcripts merged_with_reference.gtf galGal6.fa " $0}' gene_list.tab > commands

# Execute "commands" file
bash commands
``` 

- {gene_name}.cons files contain conserved regions within transcripts and could suitable for PCR primer picking. Users can go to https://www.ncbi.nlm.nih.gov/tools/primer-blast/ , paste this sequences and pick appropiate primers, specifying the genome to discard off-targets. Aditionally, users can compare a precomputed primer list for each gene here: https://gecftools.epfl.ch/getprime

### (4) I need to annotate and characterize the different types of long-noncoding RNAs in the transcriptome:

- For a detailed characterization of lncRNAs, users can use FEELnc. See here for install requirements: https://github.com/tderrien/FEELnc
- final_annotated.gtf already contains the annotated lncRNA/coding GTF in the second field. 

To obtain coding and non-coding transcripts use final_annotated.gtf file as follows:

```
grep "coding" final_annotated.gtf > final_annotated.coding.gtf
grep "lncRNA" final_annotated.gtf > final_annotated.lncRNAs.gtf

gffread -w noncoding.fa -g galGal6.fa final_annotated.coding.gtf
gffread -w coding.fa -g galGal6.fa final_annotated.lncRNAs.gtf

# Number of coding transcripts
grep ">" -c coding.fa

# Number of noncoding transcripts
grep ">" -c noncoding.fa
```

### More Scenarios?

To see more examples, please visit and clone https://github.com/cfarkas/annotate_my_genomes_examples.git

### Notes
Compiling automatically uses Shell script compiler shc to make binaries, please check: https://github.com/neurobin/shc. 
