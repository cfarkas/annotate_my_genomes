# annotate_my_genomes 
### :microscope: :hatching_chick: :hatched_chick: 
Genome annotation pipeline using long sequencing reads from non-model (and model) organisms.

## Pipeline Outline

  annotate_my_genomes are a set of bash scripts that aim to annotate transfrags obtained by genome-guided transcriptome assembly strategies coming from long sequencing reads (i.e. PacBio/Oxford Nanopore technologies). By using this sequencing, a set of high-quality transfrags can be obtained by tools such as StringTie or related (for documentation, please see http://ccb.jhu.edu/software/stringtie/) and further characterized. 
  
  Often, genomes from non-model organisms (and even from model organisms) contain reference genome annotation available in GTF format (Gene Transfer Format), but these annotations often failed to capture all genome features. Novel genes and novel transcripts can be absent from reference genome annotations due tissue or stage-specific gene expression when using RNA-seq data for transcript characterization. The pipeline aims to conciliate current gene annotation from an organism and expand this annotation by annotating transcripts with GAWN (Genome Annotation Without Nightmares, for documentation of this pipeline, please see https://github.com/enormandeau/gawn). 

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
```
#1: Just edit two parameters in /test/gawn_config.sh. You can edit:

NCPUS=10 # Increase this value to speed-up things

SWISSPROT_DB="paste_pwd_path_here/swissprot" # You need to specify full path to swissprot database. You can type pwd inside swissprot folder (found in annotate_my_genomes after make), copy this path and replace paste_pwd_path_here.

# i.e.: SWISSPROT_DB="/home/lrt/brain_chicken/PacBio/annotate_my_genomes/swissprot/swissprot"

#2 Run the pipeline with a tiny set of transcripts (Chromosome 3 from Gallus gallus genome) specifying the reference genome assembly. Go to /annotate_my_genomes/test and type in a terminal: 

bash annotate_my_genomes.sh stringtie_chr3.gtf galGal6
```

## Usage

```
# Assuming you runned the test ...

#1: Rename "genome_1" folder as you wish. Put your target GTF assembly to annotate in this folder (you can duplicate this folder every time you need to annotate a new GTF assembly)

#2: Copy gawn_config.sh file from test and put it in this folder.

#3 Run the pipeline 

cd annotate_my_genomes/your_folder/
bash annotate_my_genomes.sh target.gtf genome_assembly_name

# To check genome_assembly_names (UCSC Genome Browser assembly ID), please check: https://genome.ucsc.edu/cgi-bin/hgGateway
```
