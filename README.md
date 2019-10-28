# annotate_my_genomes 
### :microscope: :hatching_chick: :hatched_chick: 
Genome annotation pipeline using long sequencing reads from non-model (and model) organisms.

## Pipeline Outline

  annotate_my_genomes are a set of bash scripts that aim to annotate transfrags obtained by genome-guided transcriptome assembly strategies coming from long sequencing reads (i.e. PacBio/Oxford Nanopore technologies). By using this sequencing, a set of high-quality transfrags can be obtained by tools such as StringTie or related (for documentation, please see http://ccb.jhu.edu/software/stringtie/) and further characterized. 
  
  Often, genomes from non-model organisms (and even from model organisms) contain reference genome annotation available in GTF format (Gene Transfer Format), but these annotations may fail to capture all genome features. Novel genes and novel transcripts can be absent from reference genome annotations due tissue or stage-specific gene expression when using RNA-seq data for transcript characterization. The pipeline aims to conciliate current gene annotation from an organism and expand this annotation by annotating transcripts with GAWN (Genome Annotation Without Nightmares, for documentation of this pipeline, please see https://github.com/enormandeau/gawn) and classify long-noncoding RNAs with CNIT tool (http://cnit.noncode.org/CNIT/). 

The pipeline is implemented in BASH enviroment.


## Preeliminars:

## Files

### Alignment of long sequencing reads using minimap aligner (against galGal6 genome as an example, using 30 threads)
```
./minimap2 -ax splice galGal6.fa long_reads.fastq > aln_galGal6.sam -t 30
samtools view -S -b aln_galGal6.sam > aln_4_galGal6.bam --threads 30
samtools sort aln_4_galGal6.bam > aln_4_galGal6.sorted.bam --threads 30
samtools index aln_4_galGal6.sorted.bam -@ 30
```

### Obtaining transfrags from the above alignment using StringTie (using -p: 30 threads, -L: long read settings)
```
stringtie -p 30 -L -v -a 4 -o transcripts.gtf aln_4_galGal6.sorted.bam
```

## Dependences

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
#### (please see http://research-pub.gene.com/gmap/)

```
sudo apt-get install gmap
```

### Installing python dependences for CNIT: A tool for identifying protein-coding and long non-coding transcripts based on intrinsic sequence composition 
#### (please see http://cnit.noncode.org/CNIT/)

```
pip install numpy
pip install sklearn
pip install xgboost
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

2) Run the pipeline with a tiny set of transcripts (Chromosome 3 from Gallus gallus genome version "6") specifying the reference genome assembly. Go to /annotate_my_genomes/test and type in a terminal:

```
bash annotate_my_genomes.sh stringtie_chr3.gtf galGal6
```

## Usage


#### Assuming you runned the test ...

1) Place your StringTie GTF assembly to annotate in "genome_1" folder (you can duplicate this folder every time you need to annotate a new GTF assembly, changing the folder name)

2) Copy gawn_config.sh file from test and place it in "genome_1" folder.

3) Run the pipeline in genome_1 with a GTF named "target.gtf" (as an example):
```
bash annotate_my_genomes.sh target.gtf genome_assembly_name
```
#### To check genome_assembly_names (UCSC Genome Browser assembly ID), please visit: https://genome.ucsc.edu/cgi-bin/hgGateway


## Usage examples

- For mouse assembly using "target.gtf" in genome_1 folder:
```
bash annotate_my_genomes.sh target.gtf mm10
```
- For rabbit assembly using "target.gtf" in genome_1 folder:
```
bash annotate_my_genomes.sh target.gtf oryCun2
```
