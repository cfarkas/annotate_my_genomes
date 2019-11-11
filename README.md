# annotate_my_genomes 
### :microscope: :hatching_chick: :hatched_chick: 
Genome annotation pipeline using long sequencing reads from non-model (and model) vertebrate organisms.

## Pipeline Outline

  annotate_my_genomes are a set of bash scripts that aim to annotate transfrags obtained by genome-guided transcriptome assembly strategies (StringTie) coming from long read RNA-Seq alignments in vertebrate genomes (i.e. PacBio/Oxford Nanopore technologies). Transcripts are classified by its coding potential and identified as novel or reconciliated with the current reference annotation.
  
  Often, genomes from non-model organisms (and even from model organisms) contain reference genome annotation available in GTF format (Gene Transfer Format), but these annotations may fail to capture all genome features. Novel genes and novel transcripts can be absent from reference genome annotations due tissue or stage-specific gene expression when using RNA-seq data for transcript characterization. The pipeline aims to:
- Conciliate current gene annotation from an organism and expand this annotation by annotating transcripts with GAWN (Genome Annotation Without Nightmares, please see https://github.com/enormandeau/gawn). Transcripts in fasta format an concilliated GTF files are generated. 
- Predict exons with ChopStitch tool (https://github.com/bcgsc/ChopStitch.git), by implementing a Bloom filter that represents the k-mer spectrum of assembled transcripts.
- Train exons sequences along with the assembled transcripts in order to discover gene models with GlimmerHMM (https://ccb.jhu.edu/software/glimmerhmm/)
- Assess coding potential of each assembled transcript with CNIT tool (http://cnit.noncode.org/CNIT/).
- Assign to each transcripts and genes gene ontology terms (GO) and output formatted tables compatibles with WEGO annotation server: (http://wego.genomics.org.cn/). 

The input for this pipeline are GTF files coming from StringTie tool using long sequencing reads settings (for documentation, please see http://ccb.jhu.edu/software/stringtie/ and the documentation in this repository). 

The pipeline is implemented in BASH enviroment.


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
### Installing ChopStitch: exon annotation and splice graph construction using transcriptome assembly
#### (for documentation, please see https://github.com/bcgsc/ChopStitch.git)

```
git clone https://github.com/bcgsc/ChopStitch.git
cd ChopStitch/
pip install -r requirements.txt
./autogen.sh
./configure
make
sudo make install
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

2) Run the pipeline with a tiny set of transcripts from chromosome 33, Gallus gallus genome version "6" specifying the reference genome assembly (galGal6). Go to /annotate_my_genomes/test and type in a terminal:

```
bash annotate_my_genomes.sh stringtie_chr33.gtf galGal6
```

## Usage


#### Assuming you runned the test ...

1) Place your StringTie GTF assembly to annotate in "genome_1" folder (you can duplicate this folder every time you need to annotate a new GTF assembly, changing the folder name)

2) (Optional) Edit NCPUS value in gawn_config.sh file in "genome_1" folder. Default is 10. 

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

## "genesGO.tab" and "transcriptsGO.tab" output files usage:
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
