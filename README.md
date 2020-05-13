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
- Ensembl reference genome annotation (in GTF format)
- Ensembl genome assembly (masked, fasta format)

The two last requirements can be downloaded from Ensembl ftp webpage: https://uswest.ensembl.org/info/data/ftp/index.html

# Dependences:

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

# Testing if FEELnc works 

git clone https://github.com/tderrien/FEELnc
cd FEELnc/test/
# Filter
FEELnc_filter.pl -i transcript_chr38.gtf -a annotation_chr38.gtf -b transcript_biotype=protein_coding > candidate_lncRNA.gtf
# Coding_Potential
FEELnc_codpot.pl -i candidate_lncRNA.gtf -a annotation_chr38.gtf -b transcript_biotype=protein_coding -g genome_chr38.fa --mode=shuffle
# Classifier
FEELnc_classifier.pl -i feelnc_codpot_out/candidate_lncRNA.gtf.lncRNA.gtf -a annotation_chr38.gtf > candidate_lncRNA_classes.txt
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

# Installation: 

Clone this repository every time you need to work with a different assembly and/or genome. In the folder you want to annotate the GTF file, do the following: 
  
```
git clone https://github.com/cfarkas/annotate_my_genomes.git
cd annotate_my_genomes
bash makefile.sh
```

### Obtaining StringTie GTF file for annotation

#### 1) Alignment of long sequencing reads using minimap aligner (e.g.: against galGal6 genome, using 30 threads). You can use gmap as well. 
```
# Installing minimap2
git clone https://github.com/lh3/minimap2
cd minimap2 && make
sudo cp minimap2 /usr/local/bin/

# Convert PacBio subreads (bam files) to fastq
samtools bam2fq m54027_190807_082031.subreads.bam > reads.fastq

# Aligning with minimap2
minimap2 -ax splice galGal6.fa reads.fastq > aln_galGal6.sam -t 30
samtools view -S -b aln_galGal6.sam -@ 30 > aln_galGal6.bam
samtools sort aln_galGal6.bam -@ 30 > aln_galGal6.sorted.bam
samtools index aln_galGal6.sorted.bam -@ 30
```

#### 2) Obtaining GTF (transcripts.gtf) from the above alignment using StringTie (e.g.: using -p: 30 threads, -L: long read settings)
```
# Alignments from long reads (PacBio)
stringtie -p 30 -L -v -a 4 -o transcripts.gtf aln_galGal6.sorted.bam

# Alignments from long and short reads (PacBio + Illumina)
stringtie -p 30 -v -a 4 -o transcripts.gtf aln_galGal6.sorted.bam

# If the above fails, users can increase -j and -c parameters (useful for large BAM file processing)

stringtie -j 2 -c 2 -p 30 -v -a 4 -o transcripts.gtf aln_galGal6.sorted.bam
```

#

# Quickstart (Running the test)

1) Optionally, edit number of cpus in /test/gawn_config.sh:

- NCPUS=10
  - Increase this value to speed-up things :rocket:

2) Run the pipeline with a set of transcripts from chromosome 33, Gallus gallus genome version "6". Users need to specify the stringtie output (GTF format), reference genome assembly annotation (GTF format: ftp://ftp.ensembl.org/pub/release-100/gtf/gallus_gallus/Gallus_gallus.GRCg6a.100.gtf.gz), sequence (fasta format: ftp://ftp.ensembl.org/pub/release-100/fasta/gallus_gallus/dna/Gallus_gallus.GRCg6a.dna_rm.toplevel.fa.gz) and the number of threads for text processing (5 for this example). Go to /annotate_my_genomes/test and do the following:

```
# Donwload Gallus gallus v6 GTF file and decompress
wget ftp://ftp.ensembl.org/pub/release-100/gtf/gallus_gallus/Gallus_gallus.GRCg6a.100.gtf.gz
gunzip Gallus_gallus.GRCg6a.100.gtf.gz

# Download Gallus gallus v6 fasta file (Masked fasta file, with "rm" prefix) and decompress  
ftp://ftp.ensembl.org/pub/release-100/fasta/gallus_gallus/dna/Gallus_gallus.GRCg6a.dna_rm.toplevel.fa.gz
gunzip Gallus_gallus.GRCg6a.dna_rm.toplevel.fa.gz

# Execute
bash annotate_my_genomes.sh Gallus_gallus.GRCg6a.100.gtf Gallus_gallus.GRCg6a.dna_rm.toplevel.fa 5
```
#
# Usage


#### Assuming you runned the test ...

1) Place your StringTie GTF assembly to annotate in "genome_1" folder along with Ensembl reference genome (fasta and GTF files)

2) (Optional) Edit NCPUS value in gawn_config.sh file in "genome_1" folder. Default is 10. 

3) Run the pipeline in genome_1 with a GTF named "target.gtf" (as an example) with 30 threads (Important: each thread will use 1 GB of memory, check your machine):
```
bash annotate_my_genomes.sh target.gtf genome_assembly.gtf genome_assembly.fa 30
```
#### To download reference genome files (Ensembl), please visit: https://uswest.ensembl.org/info/data/ftp/index.html

#
## Usage examples

- For mouse assembly using "target.gtf" in genome_1 folder, using 30 threads for text processing:
```
bash annotate_my_genomes.sh target.gtf Mus_musculus.GRCm38.100.gtf.gz Mus_musculus.GRCm38.dna_rm.alt.fa 30
```
- For rabbit assembly using "target.gtf" in genome_1 folder, using 30 threads for text processing:
```
bash annotate_my_genomes.sh target.gtf Oryctolagus_cuniculus.OryCun2.0.100.gtf Oryctolagus_cuniculus.OryCun2.0.dna_rm.toplevel.fa 30
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
wget ftp://ftp.ensembl.org/pub/release-100/fasta/gallus_gallus/dna/Gallus_gallus.GRCg6a.dna_rm.toplevel.fa.gz
gunzip Gallus_gallus.GRCg6a.dna_rm.toplevel.fa.gz

# generate "commands" file using provided list of genes 
awk '{print "bash get_transcripts.sh merged_with_reference.gtf galGal6.fa " $0}' gene_list.tab > commands

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
