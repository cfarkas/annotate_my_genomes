# annotate_my_genomes 
### :microscope: :hatching_chick: :hatched_chick: 
Genome annotation pipeline using long sequencing reads from non-model (and model) animal organisms.

### For documentation and detailed examples, please visit our wiki page : https://github.com/cfarkas/annotate_my_genomes.wiki.git

# Pipeline Outline
  annotate_my_genomes is a pipeline that aim to annotate transfrags obtained by genome-guided transcriptome assembly strategies (StringTie) coming from long read RNA-Seq alignments in vertebrate genomes (i.e. PacBio technology). Transcripts are classified by its coding potential, probable gene function and identified as novel or reconciliated with the current reference annotation from refSeq/NCBI, without loosing novel isoforms and exon information. Also, coding sequences in nucleotides and correspondent proteins can be inferred.   

This pipeline requieres to run:

1) StringTie assembled transcripts (in GTF format)

2) USCS/NCBI reference genome annotations (in GTF format) and genome assembly (non-masked, fasta format from UCSC). All these requirements can be downloaded once by using the genome-download program provided in this repository plus genome prefix as follows: 
```
./genome-download [genome]  # mm10 for mouse, hg38 for human, galGal6 for chicken, etc. 
```
Check UCSC genome prefixes here: https://genome.ucsc.edu/cgi-bin/hgGateway. As example using prefix mm10, users will download UCSC mouse genome assembly (mm10.fa), UCSC gtf (mm10.gtf) and NCBI GTF (mm10_ncbiRefSeq.gtf).

Running the basic pipeline as follows (as example for mouse, using 20 threads):
```
./annotate-my-genomes stringtie.gtf path/to/mm10.gtf path/to/mm10.fa 20
```
will output:
```
- final_annotated.gtf: an annotated GTF file in the "gene_id" field, containing novel genes and lncRNA classification (second field in GTF file). 
- transcripts.fa : associated transcripts from final_annotated.gtf 
- cds. fa: associated coding sequences to final_annotated.gtf
- prot.fa  associated protein sequences to final_annotated.gtf
- coding_transcripts.gtf: GTF file containing cds sequences. 
- transcriptsGO.tab: annotated Gene Ontology (GO terms) for each transcript. 
```
* Users can also employ mm10_ncbiRefSeq.gtf by using "add-ncbi-annotation" instead of "annotate-my-genomes". See an example here: https://github.com/cfarkas/annotate_my_genomes/blob/master/README.md#adding-ncbi-annotations-to-increase-annotation-of-transcripts  

# Dependences 
### (installation procedures of every dependence is detailed in our wiki page)

#### Mandatory
- gcc and g++ compilers, version >= 6 
- StringTie (v2.0 release needed)
- gffcompare and gffread
- ncbi-blast+ version (v2.9.0)
- GMAP genomic aligner program 
- BEDTools
- FEELnc : FlExible Extraction of LncRNA : https://github.com/tderrien/FEELnc
- SAMtools with htslib (version >= 1.9)  : http://www.htslib.org/download/ 
- AGAT: Another Gff Analysis Toolkit (AGAT). Suite of tools to handle gene annotations in any GTF/GFF format : https://github.com/NBISweden/AGAT

#### Optional (for downstream analysis)
- EMBOSS toolkit (Open Source software for molecular biology)
- Clustal Omega (DNA/Protein alignment program)

These packages except the last three (FEELnc, SAMtools and AGAT) can be achieved via conda in a single command-line as follows:
```
conda install -c bioconda stringtie
conda install -c bioconda gffcompare
conda install -c bioconda gffread
conda install -c bioconda blast
conda install -c bioconda gmap
conda install -c bioconda bedtools
# Optional
conda install -c bioconda emboss
conda install -c bioconda clustalo
```

To install anaconda in ubuntu/linux:
```
apt-get install libgl1-mesa-glx libegl1-mesa libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2 libxi6 libxtst6
wget https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh
chmod 755 Anaconda3-2020.02-Linux-x86_64.sh
./Anaconda3-2020.02-Linux-x86_64.sh
```

# Installation: 
```
git clone https://github.com/cfarkas/annotate_my_genomes.git
cd annotate_my_genomes
# make
bash makefile.sh
```
Binaries are located in bin, genome_1 and test folders, respectively.

# Quickstart (Running the test)

1) Optionally, edit and increase the number of cpus in /test/gawn_config.sh: NCPUS=10

2) Run the pipeline with a set of transcripts from chromosome 33, Gallus gallus genome version "6". Users need to specify the stringtie output (GTF format), UCSC reference genome (GTF annotation and fasta file) and the number of threads for text processing (20 for this example). 

Go to /annotate_my_genomes/test and do the following:

```
# Download Gallus gallus v6 fasta assembly (non masked) with matched GTF files (UCSC/Ensembl)
./genome-download galGal6

# Execute in folder with 20 threads as example
./annotate-my-genomes stringtie_chr33.gtf galGal6.gtf galGal6.fa 20
```

# Usage
(Optional) Edit NCPUS value in gawn_config.sh file in "genome_1" folder. Default is 10

- For mouse assembly using "target.gtf" in genome_1 folder, using 30 threads for cpu processing:
```
./genome-download mm10
./annotate-my-genomes target.gtf mm10.gtf mm10.fa 30
```

## Adding NCBI annotations to increase annotation of transcripts
Users can add annotations from NCBI by using the three outputs from ./genome-download program by using ./add-ncbi-annotation. 
As example, the pipeline will work as follows (chicken assembly, inside test folder):
```
# Downloading galGal6 genome and correspondent UCSC/NCBI GTF annotations
./genome-download galGal6

# Running the pipeline on StringTie.gtf, using NCBI GTF (galGal6_ncbiRefSeq.gtf), UCSC GTF (galGal6.gtf), genome (galGal6.fa) and 30 threads for processing:
./add-ncbi-annotation StringTie.gtf galGal6_ncbiRefSeq.gtf galGal6.gtf galGal6.fa 30
```
final_annotated.gtf (located in output_files_NCBI) will contained the merged NCBI-updated annotation (in UCSC coordinates)


### More Scenarios?

- For downstream analysis, please visit our wiki page : https://github.com/cfarkas/annotate_my_genomes.wiki.git
- For additional examples, please visit https://github.com/cfarkas/annotate_my_genomes_examples

### Notes
Compiling automatically uses Shell script compiler shc to make binaries, please check: https://github.com/neurobin/shc.
