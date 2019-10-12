# annotate_my_genomes 
### :microscope: :hatching_chick: :hatched_chick: 
Genome annotation pipeline using long sequencing reads from non-model (and model) organisms.

## Pipeline Outline

  annotate_my_genomes are a set of bash scripts that aim to annotate transfrags obtained by genome-guided transcriptome assembly strategies from long sequencing reads (i.e. PacBio/Oxford Nanopore technologies). By using these sequencing, a set of high-quality transfrags can be obtained by tools such as StringTie or related (for documentation, please see http://ccb.jhu.edu/software/stringtie/).  
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

## Quickstart

```
#1 Clone the repo
git clone https://github.com/cfarkas/annotate_my_genomes.git

#2 Specify Number of threads in gawn_config.sh file located in test folder. Default are NCPUS=10 but you can increase this number to speed-up things.

#3 Run the pipeline with a tiny set of transcripts (Chromosome 3 StringTie assembled transcripts in Gallus gallus genome)
cd annotate_my_genomes/test/
bash annotate_my_genomes.sh stringtie_chr3.gtf galGal6
```
