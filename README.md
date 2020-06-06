# annotate_my_genomes 
### :microscope: :hatching_chick: :hatched_chick: 
Genome annotation pipeline using long sequencing reads from non-model (and model) animal organisms.

# Pipeline Outline
  Often, genomes from non-model organisms (and even from model organisms) contain reference genome annotation available in GTF format (Gene Transfer Format), but these annotations may fail to capture all genome features. Novel genes and novel transcripts can be absent from reference genome annotations due tissue or stage-specific gene expression when using RNA-seq data for transcript characterization.
  
  annotate_my_genomes is a pipeline that aim to annotate transfrags obtained by genome-guided transcriptome assembly strategies (StringTie) coming from long read RNA-Seq alignments in vertebrate genomes (i.e. PacBio technology). Transcripts are classified by its coding potential, probable gene function and identified as novel or reconciliated with the current reference annotation from refSeq/NCBI. Also, coding sequences in nucleotides and correspondent proteins sequences can be reconstructed from these procedures, with annotation of gene_id/transcript_id. 
  
  The pipeline is designed for:
  
- Use as input GTF file coming from StringTie using long sequencing reads settings (for documentation, please see http://ccb.jhu.edu/software/stringtie/ and the documentation in this repository).
- Conciliate current gene annotation from an organism with this GTF and expand this annotation by annotating novel transcripts with GAWN (Genome Annotation Without Nightmares, please see https://github.com/enormandeau/gawn). A concilliated GTF file is generated with annotated gene names and corresponding StringTie assembled transfrags (transcripts). As example:
```
gene_id "YF5"; transcript_id "STRG.40.1" ---> 40 = assembled transfrag; 1 = isoform number.
```
- The resulting GTF file will be validated by using several tools such as AGAT: Another Gff Analysis Toolkit (https://github.com/NBISweden/AGAT), and gff utilities (http://ccb.jhu.edu/software/stringtie/gff.shtml). 
- Perform gene prediction on reconstructed transcripts with Augustus software. Please see (http://augustus.gobics.de/)
- Assess coding potential of each assembled transcript with FEELnc tool (https://github.com/tderrien/FEELnc).
- Assign to each transcripts and genes gene ontology terms (GO) and output formatted tables compatibles with WEGO annotation server: (http://wego.genomics.org.cn/). 
- Complement annotation using UCSC/NCBI annotation, using genome coordinates from UCSC genomes. 

This pipeline requieres to run:

1) StringTie assembled transcripts (in GTF format)

2) USCS/NCBI reference genome annotations (in GTF format) and genome assembly (non-masked, fasta format from UCSC). All these requirements can be downloaded by using the genome-download program provided in this repository plus genome prefix as follows: 
```
./genome-download [genome]
```
Check UCSC genome prefixes here: https://genome.ucsc.edu/cgi-bin/hgGateway. As example for latest mouse assembly (mm10)
```
./genome-download mm10
```
will download UCSC mouse genome assembly (mm10.fa), UCSC gtf (mm10.gtf) and NCBI GTF (mm10_ncbiRefSeq.gtf).

# Dependences:

### gcc and g++ compilers, version >= 6 
Ubuntu/linux may come with GCC/G++ compilers <= version 6. To complile augustus gene prediction tool properly, users must upgrade old gcc/g++ compilers as follows: 
```
sudo apt-get update 
sudo apt-get install build-essential software-properties-common -y
sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y

## Important: if you have problems with these steps, users may do the following in order to fix it:
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

### Obtaining and Installing BEDTools
Complete instructions can be found in https://bedtools.readthedocs.io/en/latest/content/installation.html. Users with privileges can accomplish with sudo: 
```
sudo apt-get install bedtools
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

2) Run the pipeline with a set of transcripts from chromosome 33, Gallus gallus genome version "6". Users need to specify the stringtie output (GTF format), UCSC reference genome (GTF annotation and fasta file) and the number of threads for text processing (20 for this example). Go to /annotate_my_genomes/test and do the following:

```
# Download Gallus gallus v6 fasta assembly (non masked) with matched GTF files (UCSC/Ensembl)
./genome-download galGal6

# Execute in folder with 20 threads as example
./annotate-my-genomes stringtie_chr33.gtf galGal6.gtf galGal6.fa 20
```
#
## Usage examples

#### To download reference genome sequences (UCSC), use genome-download program from this repository. To check UCSC genome names, please visit: https://genome.ucsc.edu/cgi-bin/hgGateway

(Optional) Edit NCPUS value in gawn_config.sh file in "genome_1" folder. Default is 10

- For mouse assembly using "target.gtf" in genome_1 folder, using 30 threads for cpu processing:
```
./genome-download mm10
./annotate-my-genomes target.gtf mm10.gtf mm10.fa 30
```
- For rabbit assembly using "target.gtf" in genome_1 folder, using 30 threads for cpu processing:
```
./genome-download oryCun2
./annotate-my-genomes target.gtf oryCun2.gtf oryCun2.fa 30
```
## Adding NCBI annotations
Users can add annotations from NCBI by using the three outputs from ./genome-download program by using ./add-ncbi-annotation. 
As example, the pipeline will work as follows (chicken assembly, inside test folder):
```
# Downloading galGal6 genome and correspondent UCSC/NCBI GTF annotations
./genome-download galGal6

# Running the pipeline on StringTie.gtf, using NCBI GTF (galGal6_ncbiRefSeq.gtf), UCSC GTF (galGal6.gtf), genome (galGal6.fa) and 30 threads for processing:
./add-ncbi-annotation StringTie.gtf galGal6_ncbiRefSeq.gtf galGal6.gtf galGal6.fa 30
```
final_annotated.gtf (located in output_files_NCBI) will contained the merged NCBI-updated annotation (in UCSC coordinates)

## Overlap NCBI annotations to IsoSeq transcripts
- IsoSeq protocol is a very-well established computational framwork to obtain high-quality assembled transcripts from PacBio IsoSeq chemistry as depicted here: https://github.com/PacificBiosciences/IsoSeq
- A suitable protocol to obtain polished transcripts from PacBio raw reads can be found here: https://github.com/PacificBiosciences/IsoSeq_SA3nUP/wiki/Tutorial:-Installing-and-Running-Iso-Seq-3-using-Conda
- To install all dependences and anaCogent5.2 enviroment, follow these instructions:
https://github.com/PacificBiosciences/IsoSeq_SA3nUP/wiki/Tutorial:-Installing-and-Running-Iso-Seq-3-using-Conda#installing-isoseq-using-anaconda

As example, we will obtain high-quality IsoSeq transcripts (polished.hq.fasta) from "subreads.bam" PacBio IsoSeq raw reads, as follows: 
```
# Activate conda enviroment
conda activate anaCogent5.2

# Consensus sequences from IsoSeq raw reads using 30 threads
ccs subreads.bam ccs.bam --noPolish --minPasses 1 -j 30

# bam to fastq conversion of ccs reads
samtools bam2fq ccs.bam > ccs.fastq

# Cluster and Polishing ccs reads to obtain high-quality transcripts
isoseq3 cluster ccs.bam unpolished.bam -j 30
isoseq3 polish unpolished.transcriptset.xml subreads.bam polished.bam -j 30
mkdir final_sequences
cp polished.hq.fasta.gz polished.lq.fasta.gz polished.hq.fastq.gz polished.lq.fastq.gz ./final_sequences/

# Decompress polished.hq.fasta.gz to downstream applications
gunzip polished.hq.fasta.gz
```

Users can overlap IsoSeq transcripts (polished.hq.fasta) or consensus sequences (ccs.fastq) reads with gene annotation from NCBI by using the three outputs from ./genome-download program by using ./annotate-IsoSeq. 
As example, the pipeline will work as with ccs.fastq reads (chicken assembly, inside genome_1 folder):
```
# Downloading galGal6 genome and correspondent UCSC/NCBI GTF annotations
./genome-download galGal6

# Running the pipeline on polished.hq.fasta, using NCBI GTF (galGal6_ncbiRefSeq.gtf), UCSC GTF (galGal6.gtf), genome (galGal6.fa) and 30 threads for processing:
./annotate-IsoSeq ccs.fastq galGal6_ncbiRefSeq.gtf galGal6.gtf galGal6.fa 30
```
final_annotated.gtf (located in output_files_IsoSeq) will contained the merged NCBI-updated annotation (in UCSC coordinates)

#
## Downstream analysis using outputs:

### (1) Gene quantification procedure examples using output GTF file (final_annotated.gtf):
- Install HTSeq-count: (please see https://htseq.readthedocs.io/en/release_0.11.1/index.html)
```
sudo apt-get install build-essential python2.7-dev python-numpy python-matplotlib python-pysam python-htseq
```

- Gene-level quantification using "final_annotated.gtf" GTF file
```
htseq-count --stranded=no --format bam condition1.bam condition2.bam final_annotated.gtf > gene_counts
```
To obtain a suitable count table for R, do
```
tac gene_counts | sed "1,5{d}" | tac > count_table
```
### (2) Transcriptome metrics

To accomplish all of these examples, users will need final_annotated.gtf file (output from this pipeline) and the reference genome in fasta format. Obtaining list of genes and matched transcripts (genes-and-transcripts.tab) from final_annotated.gtf
``` 
perl -lne 'print "@m" if @m=(/((?:transcript_id|gene_id)\s+\S+)/g);' final_annotated.gtf > final_annotated.tab
sed -i 's/transcript_id //g' final_annotated.tab
sed -i 's/;/\t/g' final_annotated.tab
sed -i 's/gene_id//g' final_annotated.tab
sed -i 's/"//g' final_annotated.tab
awk '!a[$0]++' final_annotated.tab > genes_and_transcripts.tab && rm final_annotated.tab
awk '{print $1"\t"$2}' genes_and_transcripts.tab > genes-and-transcripts.tab && rm genes_and_transcripts.tab
``` 
genes-and-transcripts.tab contains the list of assembled genes and corresponding transcripts, in tabular format. This file can be used to obtain novel/known/other transcripts and further coding/lncRNA classification of these transcripts as follows:

``` 
awk '{print $1}' genes-and-transcripts.tab > genes.tab 

# Novel genes list
grep "STRG." genes.tab > novel-genes.tab

# Known genes list
grep -v "STRG." genes.tab > known-genes.tab

# Parsing final_annotated.gtf file to obtain novel/known and coding/lncRNA transcripts, respectively. 

grep -w -F -f novel-genes.tab final_annotated.gtf > novel-genes.gtf
grep -w -F -f known-genes.tab final_annotated.gtf > known-genes.gtf
grep "coding" known-genes.gtf > known-genes-coding.gtf
grep "lncRNA" known-genes.gtf > known-genes-lncRNA.gtf
grep "StringTie" known-genes.gtf > known-genes-other.gtf  # other = no lncRNA and no protein-coding
grep "coding" novel-genes.gtf > novel-genes-coding.gtf
grep "lncRNA" novel-genes.gtf > novel-genes-lncRNA.gtf
grep "StringTie" novel-genes.gtf > novel-genes-other.gtf  # other = no lncRNA and no protein-coding

# gffread can be used to obtain transcripts in each GTF file (in example, by using galGal6.fa genome)

gffread -w known-transcripts-coding.fa -g galGal6.fa known-genes-coding.gtf
gffread -w known-transcripts-lncRNA.fa -g galGal6.fa known-genes-lncRNA.gtf
gffread -w known-transcripts-other.fa -g galGal6.fa known-genes-other.gtf
gffread -w novel-transcripts-coding.fa -g galGal6.fa novel-genes-coding.gtf
gffread -w novel-transcripts-lncRNA.fa -g galGal6.fa novel-genes-lncRNA.gtf
gffread -w novel-transcripts-other.fa -g galGal6.fa novel-genes-other.gtf

# Counting coding known transcripts
grep ">" known-transcripts-coding.fa -c 
# Counting non-coding known genes
grep ">" known-transcripts-lncRNA.fa -c
# Counting other expressed features
grep ">" known-transcripts-other.fa -c

# Counting coding novel transcripts
grep ">" novel-transcripts-coding.fa -c
# Counting non-coding novel transcripts
grep ">" novel-transcripts-lncRNA.fa -c 
# Counting other expressed features
grep ">" novel-transcripts-other.fa -c
``` 
As expected, coding transcripts will surpass the number of non-coding genes known genes, but not in novel genes.
Looking the same as above, at the gene level:
```
# known coding genes counts
perl -lne 'print "@m" if @m=(/((?:transcript_id|gene_id)\s+\S+)/g);' known-genes-coding.gtf > known-genes-coding.tab
sed -i 's/transcript_id //g' known-genes-coding.tab
sed -i 's/;/\t/g' known-genes-coding.tab
sed -i 's/gene_id//g' known-genes-coding.tab
sed -i 's/"//g' known-genes-coding.tab
awk '{print $1}' known-genes-coding.tab > known-genes-coding.tabular && rm known-genes-coding.tab
awk '!a[$0]++' known-genes-coding.tabular > known-genes-coding.tab && rm known-genes-coding.tabular
cat known-genes-coding.tab | wc -l

# known lncRNA genes counts
perl -lne 'print "@m" if @m=(/((?:transcript_id|gene_id)\s+\S+)/g);' known-genes-lncRNA.gtf > known-genes-lncRNA.tab
sed -i 's/transcript_id //g' known-genes-lncRNA.tab
sed -i 's/;/\t/g' known-genes-lncRNA.tab
sed -i 's/gene_id//g' known-genes-lncRNA.tab
sed -i 's/"//g' known-genes-lncRNA.tab
awk '{print $1}' known-genes-lncRNA.tab > known-genes-lncRNA.tabular && rm known-genes-lncRNA.tab
awk '!a[$0]++' known-genes-lncRNA.tabular > known-genes-lncRNA.tab && rm known-genes-lncRNA.tabular
cat known-genes-lncRNA.tab | wc -l

# known other-features counts
perl -lne 'print "@m" if @m=(/((?:transcript_id|gene_id)\s+\S+)/g);' known-genes-other.gtf > known-genes-other.tab
sed -i 's/transcript_id //g' known-genes-other.tab
sed -i 's/;/\t/g' known-genes-other.tab
sed -i 's/gene_id//g' known-genes-other.tab
sed -i 's/"//g' known-genes-other.tab
awk '{print $1}' known-genes-other.tab > known-genes-other.tabular && rm known-genes-other.tab
awk '!a[$0]++' known-genes-other.tabular > known-genes-other.tab && rm known-genes-other.tabular
cat known-genes-other.tab | wc -l

# novel coding gene counts
perl -lne 'print "@m" if @m=(/((?:transcript_id|gene_id)\s+\S+)/g);' novel-genes-coding.gtf > novel-genes-coding.tab
sed -i 's/transcript_id //g' novel-genes-coding.tab
sed -i 's/;/\t/g' novel-genes-coding.tab
sed -i 's/gene_id//g' novel-genes-coding.tab
sed -i 's/"//g' novel-genes-coding.tab
awk '{print $1}' novel-genes-coding.tab > novel-genes-coding.tabular && rm novel-genes-coding.tab
awk '!a[$0]++' novel-genes-coding.tabular > novel-genes-coding.tab && rm novel-genes-coding.tabular
cat novel-genes-coding.tab | wc -l

# novel lncRNA genes counts
perl -lne 'print "@m" if @m=(/((?:transcript_id|gene_id)\s+\S+)/g);' novel-genes-lncRNA.gtf > novel-genes-lncRNA.tab
sed -i 's/transcript_id //g' novel-genes-lncRNA.tab
sed -i 's/;/\t/g' novel-genes-lncRNA.tab
sed -i 's/gene_id//g' novel-genes-lncRNA.tab
sed -i 's/"//g' novel-genes-lncRNA.tab
awk '{print $1}' novel-genes-lncRNA.tab > novel-genes-lncRNA.tabular && rm novel-genes-lncRNA.tab
awk '!a[$0]++' novel-genes-lncRNA.tabular > novel-genes-lncRNA.tab && rm novel-genes-lncRNA.tabular
cat novel-genes-lncRNA.tab | wc -l  

# novel other-features counts
perl -lne 'print "@m" if @m=(/((?:transcript_id|gene_id)\s+\S+)/g);' novel-genes-other.gtf > novel-genes-other.tab
sed -i 's/transcript_id //g' novel-genes-other.tab
sed -i 's/;/\t/g' novel-genes-other.tab
sed -i 's/gene_id//g' novel-genes-other.tab
sed -i 's/"//g' novel-genes-other.tab
awk '{print $1}' novel-genes-other.tab > novel-genes-other.tabular && rm novel-genes-other.tab
awk '!a[$0]++' novel-genes-other.tabular > novel-genes-other.tab && rm novel-genes-other.tabular
cat novel-genes-other.tab | wc -l  
```

### (3) I need the transcript sequences matching each gene. Also validate conserved regions with qPCR. What can I do?:

A gene list in tabular format can also be used to extract: 
- Transcripts sequences associated to each gene. 
- Align transcript sequences in order to obtain consensus sequences

This can be accomplished by copying final_annotated.gtf file and an user-provided gene list in tabular format (such as gene_list.tab) to get_transcripts folder/ . Execute the following (e.g.: for chicken genome):

```
# Get gene list in tabular format (gene_id.tab) from final_annotated.gtf. In this example we will obtain all genes in the transcriptome:

perl -lne 'print "@m" if @m=(/((?:gene_id)\s+\S+)/g);' final_annotated.gtf > gene_id.tab
sed -i 's/gene_id//g' gene_id.tab && sed -i 's/"//g' gene_id.tab && sed -i 's/;//g' gene_id.tab
awk '{print $1}' gene_id.tab > gene_id.1.tab
awk 'seen[$0]++ == 1' gene_id.1.tab > gene_id.tab
rm gene_id.1.tab

# Downloading masked Gallus gallus v6 genome
./genome-download galGal6

# generate "commands" file using provided list of genes 
awk '{print "./get-transcripts final_annotated.gtf galGal6.fa " $0}' gene_id.tab > commands

# Execute "commands" file
bash commands
``` 

- {gene_id}.fa contains transcripts per gene in fasta format. {gene_id}.gtf contains transcripts per gene in GTF format.

- {gene_id}.cons files contain conserved regions within transcripts and could suitable for PCR primer picking. Users can go to https://www.ncbi.nlm.nih.gov/tools/primer-blast/ , paste this sequences and pick appropiate primers, specifying the genome to discard off-targets. Aditionally, users can compare a precomputed primer list for each gene here: https://gecftools.epfl.ch/getprime

### (4) Identifying paralogs in novel proteins from transcriptome

Since transcriptome_annotation_table.tsv contain most closed matched to each protein in the transcriptome (using uniprot databases from different species), users can employ donwload_proteome_uniprot.pl script to download up-to-date proteome of the species of interest (providing taxid to the script) to refine these results. Then, by using novel-transcripts-coding.fa file (generated in the previous step), the predicted transcriptome can be blasted against the reference proteome to find paralogs/missed genes as follows:
```
# Donwload gallus gallus proteome (taxid 9031)
perl download_proteome_uniprot.pl 9031

# Augustus prediction of novel-transcripts-coding.fa sequences

wget http://augustus.gobics.de/binaries/augustus.2.5.5.tar.gz
gunzip augustus.2.5.5.tar.gz
tar -xvf augustus.2.5.5.tar
cd augustus.2.5.5/src/
make
cd ..
cd ..
export AUGUSTUS_CONFIG_PATH=./augustus.2.5.5/config/
./augustus.2.5.5/src/augustus --species=human --progress=true --UTR=off --uniqueGeneId=true --gff3=on novel-transcripts-coding.fa > augustus.gff3

# Converting gff3 to GTF format, collecting coding sequences and proteins with gffread and AGAT

gffread augustus.gff3 -T -o coding_transcripts.gtf
agat_sp_extract_sequences.pl -g augustus.gff3 -f novel-transcripts-coding.fa -o novel-transcripts-coding-cds.fa
agat_sp_extract_sequences.pl -g augustus.gff3 -f novel-transcripts-coding.fa -o novel-transcripts-coding-prot.fa --protein

# making blast database
makeblastdb -in 9031.fasta -dbtype 'prot' -out 9031

# Blasting protein sequences (using 20 threads, max targets= 1, output format=xml)
blastp -db 9031 -max_hsps 1 -max_target_seqs 1 -out blast_results -query novel-transcripts-coding-prot.fa -num_threads 20 -outfmt 5

# Parsing blast results (min_identity 90, min_alignment_length 100). Thanks to: https://github.com/sandyjmacdonald
git clone https://github.com/sandyjmacdonald/blast_parser.git  # requires BioPython installed. 
python ./blast_parser/blast_parser.py -i blast_results -e 1e-20 -p 90 -a 100 > parsed_results.txt
# Removing queries with no hits
awk -F'\t' '$2!=""' parsed_results.txt > parsed_results.tab
```
parsed_results.tab will contain matched proteins to the provided reference. Proteins with identity closed to 100% are good paralog candidates (if not missed genes in the reference).

### (5) get-homologues-est analysis (find orthologs in other species)

Installing get-homologues: for documentation see https://github.com/eead-csic-compbio/get_homologues
``` 
git clone https://github.com/eead-csic-compbio/get_homologues.git
cd get_homologues/
perl install.pl   # press Y in every step
``` 
In this case we will find orthologs chicken transcripts using human annotated transcripts. Download hg38 human genome using genome-download program and obtain cds sequences with AGAT
``` 
./genome-download hg38
agat_sp_extract_sequences.pl -g hg38.gtf -f hg38.fa -o hg38_cds.fa
```  
Copy cds.fa from annotate_my_genomes output_files_NCBI as galGal6_cds.fa and perform the analysis
``` 
cp ./path/to/annotate_my_genomes/genome_1/output_files_NCBI/cds.fa ./galGal6_cds.fa
mkdir chicken_vs_human
mv galGal6_cds.fa hg38_cds.fa ./chicken_vs_human
# Execute get_homologues using 20 threads, -min% sequence identity 70 and min% coverage of shortest sequence 50
perl ./get_homologues/get_homologues-est.pl -d chicken_vs_human -M -c -z -A -n 20 -S 70 -C 50
``` 
Folder inside "chicken_vs_human_est_homologues" new folder will contain the results. 

### More Scenarios?

To see more examples, please visit and clone https://github.com/cfarkas/annotate_my_genomes_examples.git

### Notes
Compiling automatically uses Shell script compiler shc to make binaries, please check: https://github.com/neurobin/shc.
