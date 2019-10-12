# annotate_my_genomes
Genome annotation pipeline using long sequencing reads.

## Pipeline Outline

annotate_my_genomes are a set of bash scripts that aim to annotate transfrags obtained by genome-guided transcriptome assembly strategies. A set of transfrags can be obtained by tools such as StringTie or related (for documentation, please see http://ccb.jhu.edu/software/stringtie/). Often, genomes from model and Non-model organisms contain reference genome annotation available in GTF format (Gene Transfer Format) but these annotations often failed to capture all genome features including novel transcripts and genes due missing, tissue and/or stage-specific gene expression. The pipeline aims to conciliate current gene annotation from an organism and expand this annotation by annotating transcripts with GAWN (Genome Annotation Without Nightmares, for documentation of this pipeline, please see https://github.com/enormandeau/gawn). 
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
git clone https://github.com/gpertea/gffcompare
git clone https://github.com/gpertea/gffread
cd gffcompare
make release
sudo cp gffcompare trmap /usr/local/bin/
cd ..
cd gffread
make release
sudo cp gffread /usr/local/bin/
cd ..
```
### Installing up-to-date ncbi-blast+ version (October 2019 version)
### If you have ncbi-blast+ version > 2.7 is OK, older binaries also work, but GAWN pipeline strongly recommends > v2.7.1+

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
### Installing GMAP genomic aligner program (please see http://research-pub.gene.com/gmap/)

```
sudo apt-get install gmap
```

## Quickstart


