# annotate_my_genomes

### :microscope: :hatching_chick: :hatched_chick: 
Transcriptome annotation pipeline using long sequencing reads from non-model (and model) animal organisms.

![Pipeline_Scheme](https://user-images.githubusercontent.com/7016350/85973038-d9469e80-b98e-11ea-864e-03a803368143.jpg)


## Pipeline Outline
  annotate_my_genomes is a pipeline that aims to annotate genome-guided transcriptome assemblies from StringTie, coming from long read RNA-Seq alignments in vertebrate genomes (i.e. PacBio technology). Transcripts are classified by its coding potential, probable gene function and identified as novel or reconciliated with the current reference annotation from refSeq/NCBI, without loosing novel isoforms and exon information. Also, known/novel coding sequences in nucleotides and correspondent proteins will be resolved.  

This pipeline requieres to run:

1) StringTie assembled transcripts (in GTF format)

2) USCS/NCBI reference genome annotations (in GTF format) and genome assembly (non-masked, fasta format from UCSC). All these requirements can be downloaded once by using the genome-download program provided in this repository plus genome prefix as follows: 
```
./genome-download [genome]  # mm10 for mouse, hg38 for human, galGal6 for chicken, etc. 
```
Check UCSC genome prefixes here: https://genome.ucsc.edu/cgi-bin/hgGateway. As example using prefix mm10, users will download UCSC mouse genome assembly (mm10.fa), UCSC gtf (mm10.gtf) and NCBI GTF (mm10_ncbiRefSeq.gtf).

Running the basic pipeline as follows (as example for mouse, using 20 threads):
```
./annotate-my-genomes stringtie.gtf /path/to/mm10.gtf /path/to/mm10.fa 20
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
conda env update --file environment.yml                        # install required programs
conda activate annotate_my_genomes                             # load environment
bash makefile.sh                                               # make  & install
```
Also install (not through conda):
- ncbi-blast+ version equal or higher than v2.7.1. To install it, see here: https://github.com/cfarkas/annotate_my_genomes/wiki#6-installing-up-to-date-ncbi-blast-version-v271

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

### Option 2: Without using conda, program by program:

- see detailed installation steps in our wiki here: https://github.com/cfarkas/annotate_my_genomes.wiki.git

## Quickstart (Running the test)

1) Inside test folder, run the pipeline with a provided set of transcripts from chromosome 33, Gallus gallus genome version "6", in GTF format. Users need to specify the stringtie output (GTF format), UCSC reference genome (GTF annotation and fasta file) and the number of threads for text processing (20 for this example). 

2) Optionally, edit and increase the number of cpus in /test/gawn_config.sh: NCPUS=10

Enter /annotate_my_genomes/test/ directory and execute the following:

```
# Download Gallus gallus v6 fasta assembly (non masked) with matched GTF files (UCSC/Ensembl)
./genome-download galGal6

# Execute pipeline on stringtie_chr33.gtf (provided file) with 20 threads as example
./annotate-my-genomes stringtie_chr33.gtf galGal6.gtf galGal6.fa 20
```

## Simplest usage
(Optional) Edit NCPUS value in gawn_config.sh file in "genome_1" folder. Default is 10

- As example, to annotate an user-provided mouse GTF file (named "target.gtf" as example) in genome_1 folder, using 30 threads for cpu processing:
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

- For downstream analysis and examples, please visit our wiki page : https://github.com/cfarkas/annotate_my_genomes.wiki.git

### Notes
Compiling automatically uses Shell script compiler shc to make binaries, please check: https://github.com/neurobin/shc.
