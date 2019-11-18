# annotate_my_genomes 
### :microscope: :hatching_chick: :hatched_chick: 
Genome annotation pipeline using long sequencing reads from non-model (and model) vertebrate organisms.

## Pipeline Outline
  Often, genomes from non-model organisms (and even from model organisms) contain reference genome annotation available in GTF format (Gene Transfer Format), but these annotations may fail to capture all genome features. Novel genes and novel transcripts can be absent from reference genome annotations due tissue or stage-specific gene expression when using RNA-seq data for transcript characterization
  
  annotate_my_genomes are a set of bash scripts that aim to annotate transfrags obtained by genome-guided transcriptome assembly strategies (StringTie) coming from long read RNA-Seq alignments in vertebrate genomes (i.e. PacBio/Oxford Nanopore technologies). Transcripts are classified by its coding potential, probable gene function and identified as novel or reconciliated with the current reference annotation. Also, coding sequences in nucleotides and correspondent proteins sequences can be reconstructed from these procedures. 
  
  The pipeline is designed for:
  
- Use as input GTF file coming from StringTie tool using long sequencing reads settings (for documentation, please see http://ccb.jhu.edu/software/stringtie/ and the documentation in this repository).
- Conciliate current gene annotation from an organism with this GTF and expand this annotation by annotating novel transcripts with GAWN (Genome Annotation Without Nightmares, please see https://github.com/enormandeau/gawn). A concilliated GTF file is generated with annotated gene names and corresponding StringTie transcripts. 
- Perform gene prediction on reconstructed transcripts with Augustus software. Please see (http://augustus.gobics.de/)
- Assess coding potential of each assembled transcript with CNIT tool (http://cnit.noncode.org/CNIT/).
- Assign to each transcripts and genes gene ontology terms (GO) and output formatted tables compatibles with WEGO annotation server: (http://wego.genomics.org.cn/). 

This pipeline requieres to run:
- StringTie assembled transcripts (in GTF format)
- genome assembly name (check UCSC format)

This repository can be cloned every time you need to work with a different genome/assembly

## Preeliminars:

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

### Installing python dependences for CNIT: A tool for identifying protein-coding and long non-coding transcripts based on intrinsic sequence composition 
#### (for documentation, please see http://cnit.noncode.org/CNIT/)

```
pip install numpy
pip install sklearn
pip install xgboost
```

### Obtaining and installing up-to-date SAMtools with htslib (version 1.9)
(Old samtools version can also work). Users needs to install version up to date of these three packages. Users can first install htslib v1.9 and then samtools with bcftools v1.9, respectively. For downloading these packages, see http://www.htslib.org/download/). The latter can be accomplish by downloading the three packages, decompressing it, and doing the following:
```
cd htslib-1.9    # and similarly for samtools
sudo ./configure --prefix=/usr/local/bin
sudo make
sudo make install
# this step is only for samtools
sudo cp samtools /usr/local/bin/
```
Then in a terminal type
>samtools

to check 1.9 version (using htslib v1.9)

## Obtaining and installing EMBOSS toolkit (Open Source software for molecular biology)
Complete instructions can be found at the webpage: https://ssbio.readthedocs.io/en/latest/instructions/emboss.html. A way to install it can be the following:
```
sudo apt-get install emboss
```
## Obtaining and installing Clustal Omega (DNA/Protein alignment program)
Description can be found at the webpage: http://www.clustal.org/omega/. A way to install it can be the following:
```
sudo apt-get install clustalo
```

### Obtaining StringTie transfrags for annotation

#### 1) Alignment of long sequencing reads using minimap aligner (e.g.: against galGal6 genome, using 30 threads). You can use gmap as well. 
```
git clone https://github.com/lh3/minimap2
cd minimap2 && make
./minimap2 -ax splice galGal6.fa long_reads.fastq > aln_galGal6.sam -t 30
samtools view -S -b aln_galGal6.sam > aln_galGal6.bam --threads 30
samtools sort aln_galGal6.bam > aln_galGal6.sorted.bam --threads 30
samtools index aln_galGal6.sorted.bam -@ 30
```

#### 2) Obtaining transfrags from the above alignment using StringTie (e.g.: using -p: 30 threads, -L: long read settings)
```
stringtie -p 30 -L -v -a 4 -o transcripts.gtf aln_4_galGal6.sorted.bam
```

## Installation
  
```
git clone https://github.com/cfarkas/annotate_my_genomes.git
cd annotate_my_genomes
bash makefile.sh
```

## Quickstart (Running the test)

1) Optionally, edit number of cpus in /test/gawn_config.sh:

- NCPUS=10
  - Increase this value to speed-up things :rocket:

2) Run the pipeline with a set of transcripts from chromosome 33, Gallus gallus genome version "6" specifying the reference genome assembly (galGal6) and the number of threads for text processing (5 for this example). Go to /annotate_my_genomes/test and type in a terminal:

```
bash annotate_my_genomes.sh stringtie_chr33.gtf galGal6 5
```

## Usage


#### Assuming you runned the test ...

1) Place your StringTie GTF assembly to annotate in "genome_1" folder (you can duplicate this folder every time you need to annotate a new GTF assembly, changing the folder name)

2) (Optional) Edit NCPUS value in gawn_config.sh file in "genome_1" folder. Default is 10. 

3) Run the pipeline in genome_1 with a GTF named "target.gtf" (as an example) with 30 threads for text processing:
```
bash annotate_my_genomes.sh target.gtf genome_assembly_name 30
```
#### To check genome_assembly_names (UCSC Genome Browser assembly ID), please visit: https://genome.ucsc.edu/cgi-bin/hgGateway


## Usage examples

- For mouse assembly using "target.gtf" in genome_1 folder, using 30 threads for text processing:
```
bash annotate_my_genomes.sh target.gtf mm10 30
```
- For rabbit assembly using "target.gtf" in genome_1 folder, using 30 threads for text processing:
```
bash annotate_my_genomes.sh target.gtf oryCun2 30
```

## Downstream analysis using outputs:

### (1) Gene quantification procedure examples using output GTF file (merged_with_reference.gtf):
- Install HTSeq-count: (please see https://htseq.readthedocs.io/en/release_0.11.1/index.html)

```
sudo apt-get install build-essential python2.7-dev python-numpy python-matplotlib python-pysam python-htseq
```

- Gene-level quantification using "merged_with_reference.gtf" GTF file

```
htseq-count -t transcript --stranded=no --format bam condition1.bam condition2.bam merged_with_reference.gtf > gene_counts
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

This can be accomplished by copying merged_with_reference.gtf file and a user-provided gene list in tabular format (such as gene_list.tab) to get_transcripts folder/ and execute the following (e.g.: for chicken genome):

```
# Downloading galGal6 genome
bash genome_download_for_transcripts.sh galGal6

# generate "commands" file using provided list of genes 
awk '{print "bash get_transcripts.sh merged_with_reference.gtf galGal6.fa " $0}' gene_list.tab > commands

# Execute "commands" file
bash commands
```

- Users can duplicate get_transcripts folder every time you need to work with a different gene list). 

- {gene_name}.cons files contain common sequences within transcripts and could suitable for PCR primer picking in conserved regions. Users can go to https://www.ncbi.nlm.nih.gov/tools/primer-blast/ , paste this sequences and pick appropiate primers, specifying the genome to discard off-targets. Aditionally, users can compare a precomputed primer list for each gene here: https://gecftools.epfl.ch/getprime
