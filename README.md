# annotate_my_genomes

Transcriptome annotation pipeline using long sequencing reads from non-model (and model) animal organisms.

![image](https://user-images.githubusercontent.com/7016350/108611599-a6319f00-73a5-11eb-89b7-3cfd44b00cc5.png)

## Pipeline Outline
  annotate_my_genomes is a pipeline that aims to annotate genome-guided transcriptome assemblies from StringTie, coming from long read RNA-Seq alignments in vertebrate genomes (i.e. PacBio technology). Transcripts are classified by its coding potential, probable gene function and identified as novel or reconciliated with the current reference annotation from refSeq/NCBI, without loosing novel isoforms and exon information. Also, known/novel coding sequences in nucleotides and correspondent proteins will be resolved.  

This pipeline requieres to run:

1) StringTie assembled transcripts (in GTF format)

2) USCS/NCBI reference genome annotations (in GTF format) and genome assembly (non-masked, fasta format from UCSC). All these requirements can be downloaded once by using the genome-download program provided in this repository and inputtimg genome prefix as follows: 
```
./genome-download [genome]  # mm10 for mouse, hg38 for human, galGal6 for chicken, etc. Use genome-download-macOSX instead in macOSX 
```
Check UCSC genome prefixes here: https://genome.ucsc.edu/cgi-bin/hgGateway. As example using prefix mm10, users will download UCSC mouse genome assembly (mm10.fa), UCSC gtf (mm10.gtf) and NCBI GTF (mm10_ncbiRefSeq.gtf).

Running the basic pipeline as follows (as example for mouse, using 20 threads):
```
./annotate-my-genomes -a stringtie.gtf -r /path/to/mm10.gtf -g /path/to/mm10.fa -t 20
```
will output:
```
- final_annotated.gtf: an annotated GTF file in the "gene_id" field, containing novel genes and lncRNA classification (second field in GTF file). 
- transcripts.fa : associated transcripts from final_annotated.gtf 
- cds. fa: associated coding sequences to final_annotated.gtf
- prot.fa  associated protein sequences to final_annotated.gtf
- coding_transcripts.gtf: GTF file containing cds sequences.
- novel coding sequences (novel-cds.fa) and correspondent novel protein sequences (novel-prot.fa).
```
* Users can also employ mm10_ncbiRefSeq.gtf by using "add-ncbi-annotation" instead of "annotate-my-genomes". See an example here: https://github.com/cfarkas/annotate_my_genomes/blob/master/README.md#adding-ncbi-annotations-to-increase-annotation-of-transcripts  

## Installation:  

### Option 1: Installing dependences via anaconda (recommended)
- requires miniconda, python2.7 and/or python>=3. To install miniconda, see: https://docs.conda.io/en/latest/miniconda.html
```
git clone https://github.com/cfarkas/annotate_my_genomes.git   # clone repository
cd annotate_my_genomes                                         # enter repository
conda config --add channels bioconda                           # add bioconda channel (if you haven't already done so)
conda config --add channels conda-forge                        # add conda-forge channel (if you haven't already done so)
conda create --name annotate_my_genomes python=3.6.13          # create environment
conda activate annotate_my_genomes                             # load environment

# Install packages
conda install -c conda-forge -y parallel 
conda install -c bioconda -y cufflinks
conda install -c bioconda -y stringtie
conda install -c bioconda -y gffcompare
conda install -c bioconda -y gffread
conda install -c bioconda -y gmap
conda install -c bioconda -y bedtools
conda install -c bioconda -y emboss
conda install -c bioconda -y clustalo
conda install -c bioconda -y samtools
conda install -c bioconda -y minimap2
conda install -c bioconda -y transdecoder
conda install -c bioconda -y seqkit
conda install -c conda-forge -y coreutils
conda install -c anaconda -y gawk
conda install -c conda-forge -y sed
conda install -c bioconda/label/cf201901 -y feelnc
conda install -c bioconda -y perl-local-lib

bash makefile.sh                                               # make  & install
```
- Optionally (requires sudo privileges)
```
sudo cp ./bin/* /usr/local/bin/
```
- Also install (not through conda):
ncbi-blast+ version equal or higher than v2.7.1. To install it, see here: https://github.com/cfarkas/annotate_my_genomes/wiki#5-installing-up-to-date-ncbi-blast-version-v271

After these steps, a conda enviroment called annotate_my_genomes can be managed as follows:
```
# To activate this environment, use
#
#     $ conda activate annotate_my_genomes
#
# To deactivate an active environment, use
#
#     $ conda deactivate
```
- By activating annotate_my_genomes enviroment, all binaries in the annotate_my_genomes repository can be executed.
- To install optional programs for downstream analysis, please see here: https://github.com/cfarkas/annotate_my_genomes/wiki#optional-dependences-to-run-all-the-downstream-analysis
- Uninstall as follows: 
```
conda remove --name annotate_my_genomes --all
```

### Option 2: Without using conda, program by program:

- see detailed installation steps in our wiki here: https://github.com/cfarkas/annotate_my_genomes/wiki

## Quickstart (Running the test)

1) Inside test folder, run the pipeline with a provided set of transcripts from chromosome 33, Gallus gallus genome version "6", in GTF format. Users need to specify the stringtie output (GTF format), UCSC reference genome (GTF annotation and fasta file) and the number of threads for text processing (20 for this example). 

2) Optionally, edit and increase the number of cpus in /test/gawn_config.sh: NCPUS=10

Enter /annotate_my_genomes/test/ directory and execute the following:

```
# Download Gallus gallus v6 fasta assembly (non masked) with matched GTF files (UCSC/Ensembl)
./genome-download galGal6        

# Execute pipeline on stringtie_chr33.gtf (provided file) with 20 threads:
./annotate-my-genomes -a stringtie_chr33.gtf -r galGal6.gtf -g galGal6.fa -t 20

# Include NCBI annptations on stringtie_chr33.gtf (provided file) with 20 threads:
./add-ncbi-annotation -a stringtie_chr33.gtf -n galGal6_ncbiRefSeq.gtf -r galGal6.gtf -g galGal6.fa -t 20
```

## Simplest usage
(Optional) Edit NCPUS value in gawn_config.sh file in "genome_1" folder. Default is 10
- As example, to annotate a chicken GTF file (i.e: "target.gtf") in genome_1 folder, using 20 threads for cpu processing:
```
./genome-download galGal6          
./annotate-my-genomes -a target.gtf -r galGal6.gtf -g galGal6.fa -t 20
```
- final_annotated.gtf (located in output_files_UCSC) will contained the merged NCBI-updated annotation (in UCSC coordinates)
- To produce target.gtf assembly, check stringtie parameteres here: https://github.com/cfarkas/annotate_my_genomes/wiki#ii-obtaining-stringtie-gtf-file-for-annotation

## Adding NCBI annotations to increase annotation of transcripts
Users can add annotations from NCBI by using the three outputs from ./genome-download program and input into ./add-ncbi-annotation. 
- Resuming the previous example, using add-ncbi-annotation instead of annotate-my-genomes:
```
./genome-download galGal6         
./add-ncbi-annotation -a target.gtf -n galGal6_ncbiRefSeq.gtf -r galGal6.gtf -g galGal6.fa -t 20
```
- final_annotated.gtf (located in output_files_NCBI) will contained the merged NCBI-updated annotation (in UCSC coordinates).

As example for mouse genome, change galGal6 prefix to mm10. Using 30 threads for processing "mouse.gtf" assembly:
```
./genome-download mm10            
./add-ncbi-annotation -a mouse.gtf -n mm10_ncbiRefSeq.gtf -r mm10.gtf -g mm10.fa -t 20
```
### Identifying homologs in novel proteins from transcriptome

- See this example: https://github.com/cfarkas/annotate_my_genomes/wiki#5-identifying-homologs-in-novel-proteins-from-transcriptome

### More Scenarios?

- For downstream analysis and examples, please visit our wiki page : https://github.com/cfarkas/annotate_my_genomes.wiki.git

### Notes
Compiling automatically uses Shell script compiler shc to make binaries, please check: https://github.com/neurobin/shc.
