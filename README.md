# annotate_my_genomes

Transcriptome annotation pipeline using short and long sequencing reads from non-model (and model) animal organisms.

![image](https://user-images.githubusercontent.com/7016350/108611599-a6319f00-73a5-11eb-89b7-3cfd44b00cc5.png)

#### See publication here: https://doi.org/10.1093/gigascience/giac099

## I) Pipeline Outline
  ```annotate_my_genomes``` is a pipeline that aims to annotate genome-guided transcriptome assemblies from StringTie, coming from long read RNA-Seq alignments in vertebrate genomes (i.e. PacBio technology). Transcripts are classified by its coding potential, probable gene function and identified as novel or reconciliated with the current reference annotation from RefSeq/NCBI, without loosing isoform and exon information. Also, known/novel coding sequences in nucleotides and correspondent proteins will be resolved.  

This pipeline requieres to run:

1) StringTie assembled transcripts (in GTF format). Check here: https://github.com/cfarkas/annotate_my_genomes/wiki#ii-obtaining-stringtie-gtf-file-for-annotation

2) At minimum, coding UCSC/NCBI reference genome annotations (in GTF format) and genome assembly (non-masked fasta from UCSC). All these requirements can be downloaded once by using the ```genome-download``` program provided in this repository and inputting a genome prefix as follows: 
```
./genome-download [genome]  # mm10 for mouse, hg38 for human, galGal6 for chicken, etc. Use genome-download-macOSX instead in macOSX 
```
- In example, ```./genome-download mm10 ``` , will output: ```mm10.fa```, ```mm10.gtf``` and ```mm10_ncbiRefSeq.gtf``` files.
- ```mm10.gtf``` contains coding genes and ```mm10_ncbiRefSeq.gtf``` contains all NCBI annotations.

- For genomes, check UCSC genome prefixes here: http://hgdownload.soe.ucsc.edu/downloads.html

3) Finally, the basic pipeline can be runned using a mouse transcriptome as example (stringtie.gtf) and 20 threads, as follows:
```
mkdir output1
./annotate-my-genomes -a /path/to/stringtie.gtf -r /path/to/mm10.gtf -g /path/to/mm10.fa -c /path/to/annotate_my_genomes/gawn_config.sh -t 20 -o /path/to/output1
```
The latter will output inside output1 folder:
```
- final_annotated.gtf: an annotated GTF file in the "gene_id" field, containing novel genes and lncRNA classification (second field in GTF file). 
- transcripts.fa : associated transcripts from final_annotated.gtf 
- cds. fa: associated coding sequences to final_annotated.gtf
- prot.fa  associated protein sequences to final_annotated.gtf
- coding_transcripts.gtf: GTF file containing cds sequences.
- novel coding sequences (novel-cds.fa) and correspondent novel protein sequences (novel-prot.fa).
```
* Users can also employ ```mm10_ncbiRefSeq.gtf``` by using ```add-ncbi-annotation``` instead of ```annotate-my-genomes``` binary. See an example here: https://github.com/cfarkas/annotate_my_genomes/blob/master/README.md#v-adding-ncbi-annotations-to-increase-annotation-of-transcripts

## II) Installation:  

### Option 1: Via Nextflow (recommended)

- Nextflow (https://www.nextflow.io/) is a great workflow framework and a programming DSL that eases the writing of data-intensive computational pipelines. We encourage and support the usage of this framework across different platforms for reproducibility. 

### Requirements: 

- Nextflow can be installed as depicted here (https://www.nextflow.io/) or via anaconda as follows:

```
conda install -c bioconda nextflow
```
Also install (not through conda):

- ```wget``` Comes by default with Linux/Ubuntu distros
- ```sed``` editor. Comes by default with Linux/Ubuntu distros
- ```ncbi-blast+``` version equal or higher than v2.7.1. To install it, see here: https://github.com/cfarkas/annotate_my_genomes/wiki#5-installing-up-to-date-ncbi-blast-version-v271
- ```SAMtools``` . To install it, see here: https://github.com/cfarkas/annotate_my_genomes/wiki#9-obtaining-and-installing-up-to-date-samtools-with-htslib-version--19


### Installation:

In a given directory:
```
git clone https://github.com/cfarkas/annotate_my_genomes.git                        # clone repository
cd annotate_my_genomes                                                              # enter repository
current_dir=$(pwd)                                                                  # set working directory
echo $current_dir                                                                   # check working directory
nextflow run makefile.nf --workdir $current_dir --conda ./22.04_environment.yml     # make & install; use environment.yml for Ubuntu < 22.04
```

### Option 2: Installing dependences via anaconda (tested in Ubuntu 16.05, 18.04, 20.04 and 22.04 LTS)

### Requirements: 
- requires miniconda, python2.7 and/or python>=3. To install miniconda, see: https://docs.conda.io/en/latest/miniconda.html

Also install (not through conda):

- ```wget``` Comes by default with Linux/Ubuntu distros
- ```sed``` editor. Comes by default with Linux/Ubuntu distros
- ```ncbi-blast+``` version equal or higher than v2.7.1. To install it, see here: https://github.com/cfarkas/annotate_my_genomes/wiki#5-installing-up-to-date-ncbi-blast-version-v271
- ```SAMtools``` . To install it, see here: https://github.com/cfarkas/annotate_my_genomes/wiki#9-obtaining-and-installing-up-to-date-samtools-with-htslib-version--19

### Installation:

In a given directory:
```
git clone https://github.com/cfarkas/annotate_my_genomes.git   # clone repository
cd annotate_my_genomes                                         # enter repository
conda config --add channels bioconda                           # add bioconda channel (if you haven't already done so)
conda config --add channels conda-forge                        # add conda-forge channel (if you haven't already done so)
conda env create -f 22.04_environment.yml                      # create and install environment; use environment.yml for Ubuntu < 22.04
conda activate annotate_my_genomes                             # activate environment
bash makefile.sh                                               # make  & install
```
- Copy binaries to ```/usr/local/bin```
```
sudo cp ./bin/* /usr/local/bin/
```

After these steps, a conda enviroment called ```annotate_my_genomes``` can be managed as follows:
```
# To activate this environment, use
#
#     $ conda activate annotate_my_genomes
#
# To deactivate an active environment, use
#
#     $ conda deactivate
```

#### Notes: 

- By activating annotate_my_genomes enviroment, all binaries in the annotate_my_genomes repository can be executed.
- To install optional programs for downstream analysis, please see here: https://github.com/cfarkas/annotate_my_genomes/wiki#optional-dependences-to-run-all-the-downstream-analysis

- Uninstall environment as follows: 
```
conda remove --name annotate_my_genomes --all
```

- Inside the repository, there is a file called ```gawn_config.sh```. Optionally, edit and increase/decrease the number of cpus for blast processing:
```
NCPUS=10
```
To a value according to the computational capacity of your machine.

### Option 3: Run through docker:
- See installation and pipeline run here: https://hub.docker.com/r/carlosfarkas/annotate_my_genomes
```
# Run docker without sudo privileges as follows:
sudo chmod 666 /var/run/docker.sock

# Downloading the docker image
docker pull carlosfarkas/annotate_my_genomes:latest

# Downloading repository
git clone https://github.com/cfarkas/annotate_my_genomes.git && cd annotate_my_genomes

# make & install using workdir
chmod 755 makefile.sh
docker run --volume $HOME:$HOME --workdir $(pwd) carlosfarkas/annotate_my_genomes ./makefile.sh         # make & install
                    
                    OR
                    
# make & install using -it (interactively)
docker run -v $(pwd):/annotate_my_genomes -it carlosfarkas/annotate_my_genomes:latest
cd annotate_my_genomes/
bash makefile.sh     
```
### Option 4: Without using conda, program by program:

- See detailed installation steps in our wiki here: https://github.com/cfarkas/annotate_my_genomes/wiki

## III) Running the whole pipeline via nextflow (recommended)

- Inside ```annotate_my_genomes``` folder, enter into ```nextflow_scripts``` subdirectory and run the full pipeline using ```--flags``` parameters. 
- NOTE 1: Users **must provide full paths to inputs in the command line**.  We recommed to split the flags with backslashes and run the pipeline exactly as follows: 
- NOTE 2: Use environment.yml for Ubuntu < 22.04
```
cd nextflow_scripts/
```
2.1) Run ```genome-download.nf``` (i.e : output galGal6 genome)
```
nextflow run genome-download.nf \
--genome galGal6 \
--conda /path/to/22.04_environment.yml --outdir /path/to/output_folder/
```
2.2) Run ```annotate-my-genomes.nf``` . Details here: https://github.com/cfarkas/annotate_my_genomes/blob/master/README.md#b-simplest-usage
```
nextflow run annotate-my-genomes.nf \
--stringtie /path/to/stringtie.gtf \
--ref_annotation /path/to/galGal6.gtf \ 
--genome /path/to/galGal6.fa \
--config /path/to/annotate_my_genomes/gawn_config.sh \
--threads 20 \
--conda /path/to/22.04_environment.yml --outdir /path/to/output_folder/
```
2.3) Run ```add-ncbi-annotation.nf``` . Details here: https://github.com/cfarkas/annotate_my_genomes/blob/master/README.md#c-adding-ncbi-annotations-to-increase-annotation-of-transcripts
```
nextflow run add-ncbi-annotation.nf \
--stringtie /path/to/stringtie.gtf \
--NCBI_annotation /path/to/galGal6_ncbiRefSeq.gtf \
--ref_annotation /path/to/galGal6.gtf \
--genome /path/to/galGal6.fa \
--config /path/to/annotate_my_genomes/gawn_config.sh \
--threads 20  \
--conda /path/to/22.04_environment.yml --outdir /path/to/output_folder/
```
2.4) Run ```isoform-identification.nf``` . Details here: https://github.com/cfarkas/annotate_my_genomes/blob/master/README.md#d-post-processing-add-ncbi-annotation-outputs
```
nextflow run isoform-identification.nf \
--NCBI_tmap /path/to/gffcompare.tmap \
--NCBI_transcripts /path/to/NCBI_transcripts.fa \
--genome_name galGal6 \
--conda /path/to/22.04_environment.yml --outdir /path/to/output_folder/
```

#### Notes: 

- Users must provide full paths to files when running nextflow scripts.

- Inside the repository, there is a file called gawn_config.sh. Optionally, edit and increase/decrease the number of cpus for blast processing:
```
NCPUS=10
```
To a value according to the computational capacity of your machine. 


## IV) Running the whole pipeline via anaconda + binaries: 

### A) Quickstart (Running the test)

- Inside ```test``` folder, run the pipeline with a provided set of transcripts from chromosome 33, Gallus gallus genome version "6", in GTF format. 
- Users need to specify the stringtie output (GTF format), UCSC reference genome (GTF annotation and fasta file), gawn_config.sh file (check NCPUS for blast, default = 10), number of threads for text processing (20 for this example) and the output folder. 

Go to ```annotate_my_genomes/test``` directory and execute the following:

```
# Download Gallus gallus v6 fasta assembly (non masked) with matched GTF files (UCSC/Ensembl)
./genome-download galGal6        

# Execute pipeline on stringtie_chr33.gtf (provided file) with 20 threads:
mkdir output1
./annotate-my-genomes -a stringtie_chr33.gtf -r galGal6.gtf -g galGal6.fa -c gawn_config.sh -t 20 -o output1

# Include NCBI annptations on stringtie_chr33.gtf (provided file) with 20 threads:
mkdir output2
./add-ncbi-annotation -a stringtie_chr33.gtf -n galGal6_ncbiRefSeq.gtf -r galGal6.gtf -g galGal6.fa -c gawn_config.sh -t 20 -o output2
```

### B) Simplest usage
(Optional) Edit NCPUS value in ```gawn_config.sh``` file inside the repository. Default is 10
- As example, to annotate a chicken GTF file (i.e: "target.gtf") using 20 threads for cpu processing:

```
mkdir output1
./genome-download galGal6          
./annotate-my-genomes -a /path/to/target.gtf -r /path/to/galGal6.gtf -g /path/to/galGal6.fa -c /path/to/gawn_config.sh -t 20 -o /path/to/output1
```
- ```final_annotated.gtf``` (located in output1/) will contained the merged NCBI-updated annotation (in UCSC coordinates)
- To produce ```target.gtf``` assembly, check stringtie parameters here: https://github.com/cfarkas/annotate_my_genomes/wiki#ii-obtaining-stringtie-gtf-file-for-annotation

### C) Adding NCBI annotations to increase annotation of transcripts
Users can add annotations from NCBI by using the three outputs from ./genome-download program as inputs into ./add-ncbi-annotation. 
- Resuming the previous example, using add-ncbi-annotation instead of annotate-my-genomes:
```
mkdir output2
./genome-download galGal6         
./add-ncbi-annotation -a /path/to/target.gtf -n /path/to/galGal6_ncbiRefSeq.gtf -r /path/to/galGal6.gtf -g /path/to/galGal6.fa -c /path/to/gawn_config.sh -t 20 -o /path/to/output2
```
- ```final_annotated.gtf``` (located in output2/) will contained the merged NCBI-updated annotation (in UCSC coordinates).

As example for mouse genome, change galGal6 prefix to mm10. Using 30 threads for processing "mouse.gtf" assembly:
```
mkdir output3
./genome-download mm10            
./add-ncbi-annotation -a /path/to/mouse.gtf -n /path/to/mm10_ncbiRefSeq.gtf -r /path/to/mm10.gtf -g /path/to/mm10.fa -c /path/to/gawn_config.sh -t 30 -o /path/to/output3
```
### D) Post processing add-ncbi-annotation outputs

If ```stringtie.gtf``` (as an example of input GTF) was annotated with ```add-ncbi-annotation```, users can produce transcripts annotation tables (csv format) using two outputs from add-ncbi-annotation pipeline as follows:

- gffcompare.tmap (inside ```output_files``` subdirectory)
- NCBI_transcripts.fa (inside ```gffcompare_outputs_NCBI``` subdirectory)

By using isoform-identification pipeline, as follows: 

```
isoform-identification -m /path/to/gffcompare.tmap -t /path/to/NCBI_transcripts.fa -g galGal6
```
In this example:
- ```gffcompare.tmap``` correspond to the transcript map output from gffcompare
- ```NCBI_transcripts.fa``` correspond to the transcripts sequences from ```stringtie.gtf```, in fasta format 
- ```galGal6``` correspond to the NCBI genome name (in this example, Gallus gallus 6 genome, galGal6). 

The outputs ```Ref_Transcript_Annotation.csv``` and ```Novel_Transcript_Annotation.csv``` files will contain detailed annotation of transcripts. Ref_Transcript_Annotation.csv should look like this:

```
ref_gene_id	ref_id	class_code	qry_gene_id	qry_id	num_exons	FPKM	TPM	Annotation Status	NCBI RefSeq Gene ID	Transcript Description	NCBI RefSeq Protein ID	Alternative Gene Name	RefSeq Transcript Info	cds_seq
OR14J1L40	XM_025145345.1	x	STRG.16902	STRG.16902.1	3	0.089321	0.347251	Model	OR14J1L40	olfactory receptor 14J1-like 40	XP_025001113.1			AATTTCATTGGAATTAAATTTATTATACGTATGACAAACTGatatgaagaagaaacagaaacaccacATAAAATCTATCAGGCTTTTCCTAAATTTTCTGTAGTCTTGAGAGCATGATGAACATCTTTCTGATAGTGAAACCGGGTATGTTGGAGTATCTTCCTGAGGGAacccttgagctcctggttcctcatgctgtagatgagggggttcaaAGCTGGAGGCACCACTGTGTATAGAAATGACACCACCAGGTCCagagatggggaggagatggagggaggcttcaggtaggcaaacatggcagtgctgacaaacagggagagcacagccaggtgagggaggcacgtggagaaggttttgtgctgtccctgctcagagggcatcctcagcacggccctgaagatctgcacataggagaagagaatgaaagcaaagcaccCAGATGCTAAAGAGGCACTGACAATAAGAAGCCAAATGTCTTTGAGATAGGAGTGTGAGCaagagagcttgaggatctgggggatttcacagaagaactgatccacagcattgccttggcacagaggcagggaaaatgtattggcagtgtgcagcagggaattaaggacccccgtgccccaggcagctgctgccatggtggcacacgctctgctgcccagcagggtccggtagtgcaggggcttgcagatggcaac
LOC100857209	XM_015272533.2	x	STRG.16904	STRG.16904.1	3	0.099526	0.386921	Model	LOC100857209	olfactory receptor 14A16-like	XP_015128019.2			catctgcagttcctgggcatggagtcctgttcagacTGCAGGAGATAATGATGAGTCGATACCATTCTCAGAGACACTCCTCCTGCAcactttgaaaatgcatttaactCCATAGCAtgagtttattttcatgagcttcAGAATCATGTAAGAAGTAGAAACTTAAGGAGCATTTAGTTTCCTATCATTTCCTAATCATATCCCAGGCTCCTGGattttttcctcataggagCTGTTTCCACATCTCTTTTCTttacccctaaccctaacttcTATGTTCTTCAACTTCTGTTAGAGAAATCTGTTTGATTGGAGGCTAAGTACATTATTCATGACTGCAGAGAATGACAATAAtttcagctggtgctgtcctttgggggaggagaggctgaaagcacatgAGGAGATTGTTCATATAACAGCAGACTGAGAAAGGTACAATTCAGGGTACTCAGAGATGTGTTCATATTTTCTGGCTCCcttcagatttctgcctccaatccttttcccttctcttagggtataaaagaaaaatccctgccctgtctctcctcttgcaaagAGGAGCAAACACCTTTGGAAACACCCTATGGTGCAGCtgtagctgtgatACCCCTGGCTCAGGCAgaagctgtggcagcagaaggccccttCCCTGCCGGGGGGCttcttccccccacacgtctccctgcagcgccctgggcagctccccgggcaggctgagtgctgagcctggcaggcggcagagtccctgccccggcacacagcccctggggcacagcagggaccctgctctgcactacagccctgggcacccggctgcacccaaacagcacagcctgcagccgtcctgggacacgcagccctcagggctgtgctctgatgctgcagcacagaagcccTCATCTGGAACAGTAGTCTTTTTCCATAGCAAGGAAACATGAAGTACTTTCAGCCAGATCTGCTATGGGATATCCCTGATTCAGTGATCCCTCCTGGAAAAACAGCTTCATTGCCTACTGCAAGAGACTTACCCTGTCAAGCGCTGTGAGCAAtgctcctccagtgagctcacatCCTACTCACACTGTACACATCCTGtaatctctttctcttttctcttctatcTTCATGTCACCTGCAGATCATGTCTatagccctgctgtgctgtacagaagagctgctcctgtgcaCAGCTGTCTCTCCGCAGCGCTGCCTGCTTTTatgagctccctgtgtcccaggagcctggcccagctcagcagc
LOC112530844	XM_025145380.1	p	STRG.16906	STRG.16906.1	1	0.192245	0.747381	Model	LOC112530844	olfactory receptor 14A16-like	XP_025001148.1			aaatcagcgggagacaagtctcatgctttcatgatcaacaagtctcagctttattgAAGCACACGCAGGCATTTATACGATAGTTAATGAGCTACTACATATGCCAAATTGGGTTCTCTTATTGGTTAGTTCTTTACGTGAGAAAGTAACCTTCAACGCTAGATACCGTGACAGTCCCGTGATGAATGCCCGATTGTTTACCGCATACCACTCAATTTTCTTAACTGCAGCATGTTcttatcacttccttgctcctgagtGAGGGCAGCACGACCTTGCCTGGTTTAATGAGCAGGGCCCTATctccttaccagctgcatcccatCATGGCCCCTCTCCCGGAGCCAGTGCTCCGGGTCCCAAAAGCTCTCCACACTTCCCCCGTTTTCTTTTGGTACGAGCCAGGTTGTATGAATCGCATCTTGAACCACCTTTTGCTAGCATTACAGTAAACAAAGCATGATTATCAGCATACCAATCACTATCTATAAGAATACACTAGATTTATgttacacacttctacaaagcattccttgtcagtaaactaacagtaaagactacacagcacaccagtattaactacagtttcaatatcccgatgaataaaataccacagtccCCACTCTGGATCAACCACTGTACCTGACCCCCACAATTAGTGCGCTTCTGAGTCTCATAACCGccaattgctcctggcagttcccagtgtCCAAGAGACCTTtctgatgagatgttttctgcaatCTGCTAAGGGAATACCAGTCGCAGCTCAGGAGTCACGGCACTGTATATGATGTCTTGCACACCATGCGGCTATCGCTCGCCGGAGTCGCCGTTGTTGTCATCGGGTTGAGATGGGTTGTTGATGTTCGGGGCTGGCTTAgtccatttactgggaacccataatgggccagatcctgtggAAACACAGCTCTCTCCTGGaagcctcccatgatgtttacaaaattccTATTGATTCCTAATTCactcaaagtttccacaaacccTTAACACCGTACagtgatattgttcagttataaacacttgggaacagatctcacagaagcttgTCCATGTTCCCTTACACGCTTCCATgcaatcagaacacagtactagATAAACAGGTtgacactcattccctgaaaggaacacatctcactcacaccacactcactctgacatttagaacaaaaaacatAGTTTATACATAACccacaatgctgacgacgtcttttAGCTTGTATCTTAATAACACTAGTGCATTAGTCAATTAGTTGCAATtcctaccccagccggcaatctaacctgtgagctcacgtatctcggggggggggggggaagcaggcacgctccttcataccctgcgtaggacgtctcctcacgccttacgggcacccccttttctatacacatacctgaTACACcaatggatggtccttgtctgtccctgcagtgatcgggtgaggaagggagaccttccaagaaatcttggggcgcgccaaaggtgtcccctctctcaatCGATCCCGCAGCCGAACAGAGCGGATCTATTCTCGTTGCAAAATTGAGTTGtagaaatcagaccctatatccggtaaggatatagagcaggcatgcGTCTATTGATGTCTATTGAtagtgcaagggggatcactccacctaacttgcacaccgtcaggagaaattgtactatagatataggtcaaactaatacataaccaatagttgacaggaattcagatacattttcattacgtccctgaaagacacattttcatgcagtataatgagacagaagaacagagggtAGTGCTGGCGCAGTTCTCATaatttgcagttgcttgcagcttgactcacagcacctggcacagcggtctctatcacagctctgcattcctttcgcctactcccatcattgttctgtgtgagacagtgatccatagcagctgttttacttgcactgacccagggggagaaaaacatgacctcgCTGGGTCAGCCGTCCATCCACAATTTCCCTGTTCTACTATTGCCTGGCCTGTGGGTGAGTTTGGGATACCCGTACTGTGTTTTACTCCCCATGTTTGCAGAAACTCCCCAAGCCTACGACTAGTGTAGGCTGGGccattgtctgtttttattcGTAGTGATATACCCATAACTGCAAAGCAACAACTGAGATGCTTTTCTACATACAtagccttttctccaggttgagcGGTGGCCCACATAAGATGACTATATGTATCTATAGACACGTGTACATATTTCAGCTGCCCGAACTCACCCACATGCATCACATCCATCTGCCTATTTTCGTTAGCTCTAAGTCCCCTGGGGTTAACTCCTAGCCCGAGACCCATACTGCCATTATGGTGGCTGCACACTGGGCACGATCTAACAATTACCTTAGCATCCTCATATGTTATCTGATATTCCCTTCTTAGCCCCTTGGCATTCTGGTGAAACATAGAGTACGCCTCTCGGGCCAGGACATGCCGGGAGACTAAAGGTCTCTGCGCCAGTGACACCAAGCGATCAGCTCTCGCATTTCCCTCTCCCAAGTCTATCTCCCATTTATGACCTCGAACATGTATTACTGCATATGAGTGCTCCCTAATtctgattgctctctgcaactgcacgAACAGCTTGTACAGCCGCcgattctgcacttcctttatgTAGGCTTCCTCTATTTGGTGGCATACTCCAGCTACATAAAGGGAGTCGGTGACCACATTAAGGGGGCCGATTAAGTTCATCATGGCCCATACAACGGCCACCAGCTCCAATGTTTGCAATAAGTCCTTATCATCGTCTGCAATGAGGTGATGTCTCCAGGAGCCgccctgctgccaggtcactgctgctgttctagacTTCTGTCCCGCATCCGTGTAAGCCGTGATTGTGTTCTGCAAGGGCGTCTCATGCTGCTTTGGTATCCGGAGCCAACTCCATTGACCAATCCAATGTAGCGGCACGTTCGGAATCTTTTCCACTGAAACCGTACTTCCAGCTCCTAAGAGAGCATCCTGTAACTCTGGACTATGCTGCACATACCATGTCAGAGTGTCCTTCTGCATTGGCAGCTGTACACACACAGGCTCCATACCTATGATCTGCAGGGTACGTTCTCGCCCTTTcttaatcacttctgccaggagttcagttttttgaagaagtgtttttgattgctgcagtgagggacagATCCACTCTAGTACCCATACctcccccgttttctttttagattgtgCCAACGCTCCTAAAAGGTACTTTGGTCCATACCATACCATAACCTGTATGGGGAGGTCAGGGTCACGTCTCCGAACACTGCCGTGTATAATGCAGTCCATAATCTGTTGTAGTAGACGTTTGTGCTGCGTTGTCACCGTTACAGGCTGGGCCGGGTCAGTGCCCTGTAACAAAGGTCGCAACGACTCTAAGAGTTCGTTTGGGATGCCCACCACAGGGCTCAACCACTTTAAGTCCCCCAGTAACCTTTGGGCATCATGTAGAGTCTCTAGTTTAGTATCcagttgcagtttctgtggggTTACTATCGTGTTAGTCAGTGTCCATCCTAAGTACTTCCGGGGCGCGGAGAGTTGTACCTTTTCAGGGGCAAACATAAGTTCTTCCCTATTTAGGGTCTTTTCTATTTGCCaaatttgttcctgtgtgaaggcCTCTGGCTGGGCAAAAAGGATGTCCTCCATGTAATGATAAAtgaccatttgtttccattctcgCCGGAGTGGTTGTAGAGCATGATCGACATATAGTTGACATCGCGTGGGGCTATTTTTCATCCTTTGAGGTAATACTGTCCATTCAAAACGTTGATCAGGGTGTTCTCGATTCAATGCAGGCAATGTGAAGGCAAATCGTTTAGTGTCCTGAGGGTGCAGGGTAATAGTaaagaaacagtcctttaaGTCACTAATTAGTAATGGCCAATTGTAAGGTAGCATGGCAGGATTAGGCAGGGCGGGTTGAAGTGCCCCCAGTTGAGAGAGCACATTGTGGCCAATTAAGCATTGAACAGTGGGGGGTAGAGGTGCCACCGAGACAGAGGTATGGACTACTTGTTCATCAAGGTGGATTTGCAGGGGAGGTGACTTTTTCGCTAAGGATAGTCCACCTGTACCCGTCACTGTGGCTATGGCCGCTTGCAGTGGCCATTGAGGCGGCCAAATTTCTGGGCTCAATATGCTGTTGTCGGCCCCTGTATCTAATAGACCttgaagtttgatttcttcctctctgtgtttAAGTGTCACTGGTTTTTTAGGTCGATCATGCAAATTTAGTGATAGCAATGCTAAGTCCCCTGAGGAGCCAAACCCTTGCTCCCCTCGGGGAGACGATTGACACGGTGTTAAGGCTTTGGTCAATTGCTCTAGGGGTACTAACTGCGCTATCCGTTGccctttctcaatttttattggAGGAAACGGGGTGTATACCATAATCTGGATCTCACCCTGAAAGTCCGCATCTATTACCCCAGGGAGGACAAAAAGTCCGAGCATCGATGCTGAAGAACGCCCCAATAAAAGGGCCCCAACAGCGGTTCCATTTATCATTACTGGTCCCCTGATCCCTGTAGACACCCGCTCAGGTTTTGTGGTCATTAAGGTCGTGGTCACTGCGGCTGCCAAGTCCAAGCCGAGGCTTCCTGGTGTGGCTgattgcagggctgctgctggctggaaacGGCTACTTGTGTCTGTGCGTGGCCGTCGTTTCTTTCTCGCGCTGGGCTGGGGGTTTCCTGACCGGCGTCGACAGGCATTGGTATTGTGGTTGTCCATACGACATGTGTGACACCATGAACCGGTGGTTTGACACTGACGACGCATATGTCCCATGCCGCCACAGCGATAGCATTTGATGCGACCAGCAACAGGCGATCTCGGGCCTAAATTTGTTATCGCAGACGCTTGTAAGGATGCAAGAGCTGCTAGCACTTGATTGTGAGAGGCCTCAGCTTGCGCCTTTAAACTTGCCCCTAACTCCTTAATAGCCTCAATCAGAAATGCTTGGGGCCCGACTGGCACGCTTGATagcttttccagtgcctcttcAATAGTCCAATTACTCCTCAAAGTACTCAGAGTACTACGTGCTGTTGAATTACAATTTTGGAGCGCGCATTGTTTTAACATTACTCCTCTCATATACTCTGGCACCCCTGCTTTTTCAATAGCCCCGGCTACCTTATCTATGAATGCCCCAAAGTCCTCATCTCTACCTTGTCGGATCCCCATATAAAATGGCAATCCATCAGGCACCTTAATCTTGTCCATGGCCTGTCTAGCTAAATACATCGTTTCTCGACATTTATCTGGCCCTAATAATGCTTGGGCTTGTGTTCTGAAAAAAGGCCCTAGCCCTAAGAGTTCTTCGATAGTTACACCATGTAGTGGGTCTCCCGGCTGCCTAGCCTTTGAGACACTCTGATGGCACAGTTCTTGCCAATATGCattaaacaacagctgttgATGTTGTGAAGAGATCAATTTTGCTATTGCCCGACAATCGGATGGCAGCAATATCTGCGTACTCCAAATATAATCCAATATCTGCTTAGCTGGCTCGCTTTTTACCCCAAACTGACTAACTGTAGATCGTAGCTGCGATAATAATTTCCAATCTAAAGCTGTGATGGTGGCCTGCATCCCTCCCGCAGGATTAGAGGCATATATCACTGGAAACGCCATGTGCCGCACGGCCTCC
```

## V) Annotate and identify homologs in novel proteins from transcriptome

- See this example: https://github.com/cfarkas/annotate_my_genomes/wiki#5-annotate-and-identify-homologs-in-novel-proteins-from-transcriptome

## VI Annotation of BRAKER2 / TSEBRA gtf output

- The output ```braker.gtf``` from BRAKER2 pipeline (https://github.com/Gaius-Augustus/BRAKER) or ```tsebra.gtf``` from TSEBRA pipeline (https://github.com/Gaius-Augustus/TSEBRA) can be annotated using a few tools before running the pipeline.

As a requirement, the AGAT toolkit (https://github.com/NBISweden/AGAT) must be installed:
```
conda activate annotate_my_genomes
conda install -c bioconda agat
```
- Suppose you recently annotated the Gallus gallus genome (galGal6) using BRAKER2 or TSEBRA. The ```braker.gtf / tsebra.gtf``` output can be pre-processed as follows:

#### BRAKER2 run
```
agat_convert_sp_gff2gtf.pl --gff braker.gtf -o braker_fixed.gtf                        # clean and fix braker.gtf with AGAT                         
stringtie --merge -G galGal6_ncbiRefSeq.gtf braker_fixed.gtf -o braker_merged.gtf      # merge braker.gtf with reference genome GTF (i.e.: galGal6_ncbiRefSeq.gtf)
sed 's/ gene_name.*//'g braker_merged.gtf > braker_fixed.gtf                           # fix additional entries
grep "StringTie" braker_fixed.gtf > braker_stringtie.gtf                               # Exclude reference transcripts not found in braker annotation
```
- Now, ``` braker_stringtie.gtf``` can annotated as follows (i.e. using 30 threads for processing):
```
mkdir braker_annotated
add-ncbi-annotation -a braker_stringtie.gtf -n galGal6_ncbiRefSeq.gtf -r galGal6.gtf -g galGal6.fa -c gawn_config.sh -t 30 -o braker_annotated/
```

#### TSEBRA run
```
agat_convert_sp_gff2gtf.pl --gff tsebra.gtf -o tsebra_fixed.gtf                        # clean and fix tsebra.gtf with AGAT                         
stringtie --merge -G galGal6_ncbiRefSeq.gtf tsebra_fixed.gtf -o tsebra_merged.gtf      # merge tsebra.gtf with reference genome GTF (i.e.: galGal6_ncbiRefSeq.gtf)
sed 's/ gene_name.*//'g tsebra_merged.gtf > tsebra_fixed.gtf                           # fix additional entries
grep "StringTie" tsebra_fixed.gtf > tsebra_stringtie.gtf                               # Exclude reference transcripts not found in braker annotation
```
- Now, ``` tsebra_stringtie.gtf``` can annotated as follows (i.e. using 30 threads for processing):
```
mkdir tsebra_annotated
add-ncbi-annotation -a tsebra_stringtie.gtf -n galGal6_ncbiRefSeq.gtf -r galGal6.gtf -g galGal6.fa -c gawn_config.sh -t 30 -o tsebra_annotated/
```
### More Scenarios?

- For downstream analysis and examples, please visit our wiki page : https://github.com/cfarkas/annotate_my_genomes/wiki

### Notes
Compiling automatically uses Shell script compiler shc to make binaries, please check: https://github.com/neurobin/shc.
