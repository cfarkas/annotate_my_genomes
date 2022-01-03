# annotate_my_genomes

Transcriptome annotation pipeline using long sequencing reads from non-model (and model) animal organisms.

![image](https://user-images.githubusercontent.com/7016350/108611599-a6319f00-73a5-11eb-89b7-3cfd44b00cc5.png)

## I) Pipeline Outline
  annotate_my_genomes is a pipeline that aims to annotate genome-guided transcriptome assemblies from StringTie, coming from long read RNA-Seq alignments in vertebrate genomes (i.e. PacBio technology). Transcripts are classified by its coding potential, probable gene function and identified as novel or reconciliated with the current reference annotation from refSeq/NCBI, without loosing novel isoforms and exon information. Also, known/novel coding sequences in nucleotides and correspondent proteins will be resolved.  

This pipeline requieres to run:

1) StringTie assembled transcripts (in GTF format). Check here: https://github.com/cfarkas/annotate_my_genomes/wiki#ii-obtaining-stringtie-gtf-file-for-annotation

2) USCS/NCBI reference genome annotations (in GTF format) and genome assembly (non-masked, fasta format from UCSC). All these requirements can be downloaded once by using the genome-download program provided in this repository and inputting genome prefix as follows: 
```
./genome-download [genome]  # mm10 for mouse, hg38 for human, galGal6 for chicken, etc. Use genome-download-macOSX instead in macOSX 
```
- In example, ```./genome-download mm10 ``` , will output: ```mm10.fa```, ```mm10.gtf``` and ```mm10_ncbiRefSeq.gtf``` files

- For genomes, check UCSC genome prefixes here: https://genome.ucsc.edu/cgi-bin/hgGateway.

3) Finally, the basic pipeline can be runned using a mouse transcriptome as example (stringtie.gtf) and 20 threads, as follows:
```
mkdir output1
./annotate-my-genomes -a /path/to/stringtie.gtf -r /path/to/mm10.gtf -g /path/to/mm10.fa -c /path/to/gawn_config.sh -t 20 -o /path/to/output1
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
* Users can also employ mm10_ncbiRefSeq.gtf by using "add-ncbi-annotation" instead of "annotate-my-genomes". See an example here: https://github.com/cfarkas/annotate_my_genomes/blob/master/README.md#adding-ncbi-annotations-to-increase-annotation-of-transcripts

## II) Installation:  

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
conda install -c anaconda -y pandas

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
- Inside the repository, there is a file called gawn_config.sh. Optionally, edit and increase/decrease the number of cpus for blast processing:
```
NCPUS=10
```
To a value according to the computational capacity of your machine. 

### Option 2: Using Nextflow: 

- https://github.com/cfarkas/annotate_my_genomes#vii-nextflow

### Option 3: Without using conda, program by program:

- see detailed installation steps in our wiki here: https://github.com/cfarkas/annotate_my_genomes/wiki

## III) Quickstart (Running the test)

- Inside test folder, run the pipeline with a provided set of transcripts from chromosome 33, Gallus gallus genome version "6", in GTF format. 
- Users need to specify the stringtie output (GTF format), UCSC reference genome (GTF annotation and fasta file), gawn_config.sh file (check NCPUS for blast, default = 10), number of threads for text processing (20 for this example) and the output folder. 

Enter /annotate_my_genomes/test/ directory and execute the following:

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

## IV) Simplest usage
(Optional) Edit NCPUS value in gawn_config.sh file inside the repository. Default is 10
- As example, to annotate a chicken GTF file (i.e: "target.gtf") using 20 threads for cpu processing:

```
mkdir output1
./genome-download galGal6          
./annotate-my-genomes -a /path/to/target.gtf -r /path/to/galGal6.gtf -g /path/to/galGal6.fa -c /path/to/gawn_config.sh -t 20 -o /path/to/output1
```
- final_annotated.gtf (located in output1/output_files/) will contained the merged NCBI-updated annotation (in UCSC coordinates)
- To produce target.gtf assembly, check stringtie parameteres here: https://github.com/cfarkas/annotate_my_genomes/wiki#ii-obtaining-stringtie-gtf-file-for-annotation

## V) Adding NCBI annotations to increase annotation of transcripts
Users can add annotations from NCBI by using the three outputs from ./genome-download program as inputs into ./add-ncbi-annotation. 
- Resuming the previous example, using add-ncbi-annotation instead of annotate-my-genomes:
```
mkdir output2
./genome-download galGal6         
./add-ncbi-annotation -a /path/to/target.gtf -n /path/to/galGal6_ncbiRefSeq.gtf -r /path/to/galGal6.gtf -g /path/to/galGal6.fa -c /path/to/gawn_config.sh -t 20 -o /path/to/output2
```
- final_annotated.gtf (located in output2/output_files/) will contained the merged NCBI-updated annotation (in UCSC coordinates).

As example for mouse genome, change galGal6 prefix to mm10. Using 30 threads for processing "mouse.gtf" assembly:
```
mkdir output3
./genome-download mm10            
./add-ncbi-annotation -a /path/to/mouse.gtf -n /path/to/mm10_ncbiRefSeq.gtf -r /path/to/mm10.gtf -g /path/to/mm10.fa -c /path/to/gawn_config.sh -t 30 -o /path/to/output3
```
## VI) Post processing add-ncbi-annotation outputs

Users can produce transcripts annotation tables (csv format) using two outputs from add-ncbi-annotation pipeline:

- tmap output (inside output_files subdirectory)
- NCBI_transcripts.fa (inside output_files subdirectory)

By using isoform-identification pipeline, as follows: 

```
isoform-identification -m /path/to/NCBI_compare.stringtie_chr33.gtf.tmap -t /path/to/NCBI_transcripts.fa -g galGal6
```
In this example:
- NCBI_compare.stringtie.gtf.tmap correspond to the transcript map output (tmap) 
- NCBI_transcripts.fa correspond to the transcripts sequences in fasta format 
- galGal6 correspond to the NCBI genome name (in this example, Gallus gallus 6 genome, galGal6). 

The outputs will be Ref_Transcript_Annotation.csv Novel_Transcript_Annotation.csv files that will contain detailed annotation of transcripts. Ref_Transcript_Annotation.csv should look as follows:

```
ref_gene_id	ref_id	class_code	qry_gene_id	qry_id	num_exons	FPKM	TPM	Annotation Status	NCBI RefSeq Gene ID	Transcript Description	NCBI RefSeq Protein ID	Alternative Gene Name	RefSeq Transcript Info	cds_seq
OR14J1L40	XM_025145345.1	x	STRG.16902	STRG.16902.1	3	0.089321	0.347251	Model	OR14J1L40	olfactory receptor 14J1-like 40	XP_025001113.1			AATTTCATTGGAATTAAATTTATTATACGTATGACAAACTGatatgaagaagaaacagaaacaccacATAAAATCTATCAGGCTTTTCCTAAATTTTCTGTAGTCTTGAGAGCATGATGAACATCTTTCTGATAGTGAAACCGGGTATGTTGGAGTATCTTCCTGAGGGAacccttgagctcctggttcctcatgctgtagatgagggggttcaaAGCTGGAGGCACCACTGTGTATAGAAATGACACCACCAGGTCCagagatggggaggagatggagggaggcttcaggtaggcaaacatggcagtgctgacaaacagggagagcacagccaggtgagggaggcacgtggagaaggttttgtgctgtccctgctcagagggcatcctcagcacggccctgaagatctgcacataggagaagagaatgaaagcaaagcaccCAGATGCTAAAGAGGCACTGACAATAAGAAGCCAAATGTCTTTGAGATAGGAGTGTGAGCaagagagcttgaggatctgggggatttcacagaagaactgatccacagcattgccttggcacagaggcagggaaaatgtattggcagtgtgcagcagggaattaaggacccccgtgccccaggcagctgctgccatggtggcacacgctctgctgcccagcagggtccggtagtgcaggggcttgcagatggcaac
LOC100857209	XM_015272533.2	x	STRG.16904	STRG.16904.1	3	0.099526	0.386921	Model	LOC100857209	olfactory receptor 14A16-like	XP_015128019.2			catctgcagttcctgggcatggagtcctgttcagacTGCAGGAGATAATGATGAGTCGATACCATTCTCAGAGACACTCCTCCTGCAcactttgaaaatgcatttaactCCATAGCAtgagtttattttcatgagcttcAGAATCATGTAAGAAGTAGAAACTTAAGGAGCATTTAGTTTCCTATCATTTCCTAATCATATCCCAGGCTCCTGGattttttcctcataggagCTGTTTCCACATCTCTTTTCTttacccctaaccctaacttcTATGTTCTTCAACTTCTGTTAGAGAAATCTGTTTGATTGGAGGCTAAGTACATTATTCATGACTGCAGAGAATGACAATAAtttcagctggtgctgtcctttgggggaggagaggctgaaagcacatgAGGAGATTGTTCATATAACAGCAGACTGAGAAAGGTACAATTCAGGGTACTCAGAGATGTGTTCATATTTTCTGGCTCCcttcagatttctgcctccaatccttttcccttctcttagggtataaaagaaaaatccctgccctgtctctcctcttgcaaagAGGAGCAAACACCTTTGGAAACACCCTATGGTGCAGCtgtagctgtgatACCCCTGGCTCAGGCAgaagctgtggcagcagaaggccccttCCCTGCCGGGGGGCttcttccccccacacgtctccctgcagcgccctgggcagctccccgggcaggctgagtgctgagcctggcaggcggcagagtccctgccccggcacacagcccctggggcacagcagggaccctgctctgcactacagccctgggcacccggctgcacccaaacagcacagcctgcagccgtcctgggacacgcagccctcagggctgtgctctgatgctgcagcacagaagcccTCATCTGGAACAGTAGTCTTTTTCCATAGCAAGGAAACATGAAGTACTTTCAGCCAGATCTGCTATGGGATATCCCTGATTCAGTGATCCCTCCTGGAAAAACAGCTTCATTGCCTACTGCAAGAGACTTACCCTGTCAAGCGCTGTGAGCAAtgctcctccagtgagctcacatCCTACTCACACTGTACACATCCTGtaatctctttctcttttctcttctatcTTCATGTCACCTGCAGATCATGTCTatagccctgctgtgctgtacagaagagctgctcctgtgcaCAGCTGTCTCTCCGCAGCGCTGCCTGCTTTTatgagctccctgtgtcccaggagcctggcccagctcagcagc
LOC112530844	XM_025145380.1	p	STRG.16906	STRG.16906.1	1	0.192245	0.747381	Model	LOC112530844	olfactory receptor 14A16-like	XP_025001148.1			aaatcagcgggagacaagtctcatgctttcatgatcaacaagtctcagctttattgAAGCACACGCAGGCATTTATACGATAGTTAATGAGCTACTACATATGCCAAATTGGGTTCTCTTATTGGTTAGTTCTTTACGTGAGAAAGTAACCTTCAACGCTAGATACCGTGACAGTCCCGTGATGAATGCCCGATTGTTTACCGCATACCACTCAATTTTCTTAACTGCAGCATGTTcttatcacttccttgctcctgagtGAGGGCAGCACGACCTTGCCTGGTTTAATGAGCAGGGCCCTATctccttaccagctgcatcccatCATGGCCCCTCTCCCGGAGCCAGTGCTCCGGGTCCCAAAAGCTCTCCACACTTCCCCCGTTTTCTTTTGGTACGAGCCAGGTTGTATGAATCGCATCTTGAACCACCTTTTGCTAGCATTACAGTAAACAAAGCATGATTATCAGCATACCAATCACTATCTATAAGAATACACTAGATTTATgttacacacttctacaaagcattccttgtcagtaaactaacagtaaagactacacagcacaccagtattaactacagtttcaatatcccgatgaataaaataccacagtccCCACTCTGGATCAACCACTGTACCTGACCCCCACAATTAGTGCGCTTCTGAGTCTCATAACCGccaattgctcctggcagttcccagtgtCCAAGAGACCTTtctgatgagatgttttctgcaatCTGCTAAGGGAATACCAGTCGCAGCTCAGGAGTCACGGCACTGTATATGATGTCTTGCACACCATGCGGCTATCGCTCGCCGGAGTCGCCGTTGTTGTCATCGGGTTGAGATGGGTTGTTGATGTTCGGGGCTGGCTTAgtccatttactgggaacccataatgggccagatcctgtggAAACACAGCTCTCTCCTGGaagcctcccatgatgtttacaaaattccTATTGATTCCTAATTCactcaaagtttccacaaacccTTAACACCGTACagtgatattgttcagttataaacacttgggaacagatctcacagaagcttgTCCATGTTCCCTTACACGCTTCCATgcaatcagaacacagtactagATAAACAGGTtgacactcattccctgaaaggaacacatctcactcacaccacactcactctgacatttagaacaaaaaacatAGTTTATACATAACccacaatgctgacgacgtcttttAGCTTGTATCTTAATAACACTAGTGCATTAGTCAATTAGTTGCAATtcctaccccagccggcaatctaacctgtgagctcacgtatctcggggggggggggggaagcaggcacgctccttcataccctgcgtaggacgtctcctcacgccttacgggcacccccttttctatacacatacctgaTACACcaatggatggtccttgtctgtccctgcagtgatcgggtgaggaagggagaccttccaagaaatcttggggcgcgccaaaggtgtcccctctctcaatCGATCCCGCAGCCGAACAGAGCGGATCTATTCTCGTTGCAAAATTGAGTTGtagaaatcagaccctatatccggtaaggatatagagcaggcatgcGTCTATTGATGTCTATTGAtagtgcaagggggatcactccacctaacttgcacaccgtcaggagaaattgtactatagatataggtcaaactaatacataaccaatagttgacaggaattcagatacattttcattacgtccctgaaagacacattttcatgcagtataatgagacagaagaacagagggtAGTGCTGGCGCAGTTCTCATaatttgcagttgcttgcagcttgactcacagcacctggcacagcggtctctatcacagctctgcattcctttcgcctactcccatcattgttctgtgtgagacagtgatccatagcagctgttttacttgcactgacccagggggagaaaaacatgacctcgCTGGGTCAGCCGTCCATCCACAATTTCCCTGTTCTACTATTGCCTGGCCTGTGGGTGAGTTTGGGATACCCGTACTGTGTTTTACTCCCCATGTTTGCAGAAACTCCCCAAGCCTACGACTAGTGTAGGCTGGGccattgtctgtttttattcGTAGTGATATACCCATAACTGCAAAGCAACAACTGAGATGCTTTTCTACATACAtagccttttctccaggttgagcGGTGGCCCACATAAGATGACTATATGTATCTATAGACACGTGTACATATTTCAGCTGCCCGAACTCACCCACATGCATCACATCCATCTGCCTATTTTCGTTAGCTCTAAGTCCCCTGGGGTTAACTCCTAGCCCGAGACCCATACTGCCATTATGGTGGCTGCACACTGGGCACGATCTAACAATTACCTTAGCATCCTCATATGTTATCTGATATTCCCTTCTTAGCCCCTTGGCATTCTGGTGAAACATAGAGTACGCCTCTCGGGCCAGGACATGCCGGGAGACTAAAGGTCTCTGCGCCAGTGACACCAAGCGATCAGCTCTCGCATTTCCCTCTCCCAAGTCTATCTCCCATTTATGACCTCGAACATGTATTACTGCATATGAGTGCTCCCTAATtctgattgctctctgcaactgcacgAACAGCTTGTACAGCCGCcgattctgcacttcctttatgTAGGCTTCCTCTATTTGGTGGCATACTCCAGCTACATAAAGGGAGTCGGTGACCACATTAAGGGGGCCGATTAAGTTCATCATGGCCCATACAACGGCCACCAGCTCCAATGTTTGCAATAAGTCCTTATCATCGTCTGCAATGAGGTGATGTCTCCAGGAGCCgccctgctgccaggtcactgctgctgttctagacTTCTGTCCCGCATCCGTGTAAGCCGTGATTGTGTTCTGCAAGGGCGTCTCATGCTGCTTTGGTATCCGGAGCCAACTCCATTGACCAATCCAATGTAGCGGCACGTTCGGAATCTTTTCCACTGAAACCGTACTTCCAGCTCCTAAGAGAGCATCCTGTAACTCTGGACTATGCTGCACATACCATGTCAGAGTGTCCTTCTGCATTGGCAGCTGTACACACACAGGCTCCATACCTATGATCTGCAGGGTACGTTCTCGCCCTTTcttaatcacttctgccaggagttcagttttttgaagaagtgtttttgattgctgcagtgagggacagATCCACTCTAGTACCCATACctcccccgttttctttttagattgtgCCAACGCTCCTAAAAGGTACTTTGGTCCATACCATACCATAACCTGTATGGGGAGGTCAGGGTCACGTCTCCGAACACTGCCGTGTATAATGCAGTCCATAATCTGTTGTAGTAGACGTTTGTGCTGCGTTGTCACCGTTACAGGCTGGGCCGGGTCAGTGCCCTGTAACAAAGGTCGCAACGACTCTAAGAGTTCGTTTGGGATGCCCACCACAGGGCTCAACCACTTTAAGTCCCCCAGTAACCTTTGGGCATCATGTAGAGTCTCTAGTTTAGTATCcagttgcagtttctgtggggTTACTATCGTGTTAGTCAGTGTCCATCCTAAGTACTTCCGGGGCGCGGAGAGTTGTACCTTTTCAGGGGCAAACATAAGTTCTTCCCTATTTAGGGTCTTTTCTATTTGCCaaatttgttcctgtgtgaaggcCTCTGGCTGGGCAAAAAGGATGTCCTCCATGTAATGATAAAtgaccatttgtttccattctcgCCGGAGTGGTTGTAGAGCATGATCGACATATAGTTGACATCGCGTGGGGCTATTTTTCATCCTTTGAGGTAATACTGTCCATTCAAAACGTTGATCAGGGTGTTCTCGATTCAATGCAGGCAATGTGAAGGCAAATCGTTTAGTGTCCTGAGGGTGCAGGGTAATAGTaaagaaacagtcctttaaGTCACTAATTAGTAATGGCCAATTGTAAGGTAGCATGGCAGGATTAGGCAGGGCGGGTTGAAGTGCCCCCAGTTGAGAGAGCACATTGTGGCCAATTAAGCATTGAACAGTGGGGGGTAGAGGTGCCACCGAGACAGAGGTATGGACTACTTGTTCATCAAGGTGGATTTGCAGGGGAGGTGACTTTTTCGCTAAGGATAGTCCACCTGTACCCGTCACTGTGGCTATGGCCGCTTGCAGTGGCCATTGAGGCGGCCAAATTTCTGGGCTCAATATGCTGTTGTCGGCCCCTGTATCTAATAGACCttgaagtttgatttcttcctctctgtgtttAAGTGTCACTGGTTTTTTAGGTCGATCATGCAAATTTAGTGATAGCAATGCTAAGTCCCCTGAGGAGCCAAACCCTTGCTCCCCTCGGGGAGACGATTGACACGGTGTTAAGGCTTTGGTCAATTGCTCTAGGGGTACTAACTGCGCTATCCGTTGccctttctcaatttttattggAGGAAACGGGGTGTATACCATAATCTGGATCTCACCCTGAAAGTCCGCATCTATTACCCCAGGGAGGACAAAAAGTCCGAGCATCGATGCTGAAGAACGCCCCAATAAAAGGGCCCCAACAGCGGTTCCATTTATCATTACTGGTCCCCTGATCCCTGTAGACACCCGCTCAGGTTTTGTGGTCATTAAGGTCGTGGTCACTGCGGCTGCCAAGTCCAAGCCGAGGCTTCCTGGTGTGGCTgattgcagggctgctgctggctggaaacGGCTACTTGTGTCTGTGCGTGGCCGTCGTTTCTTTCTCGCGCTGGGCTGGGGGTTTCCTGACCGGCGTCGACAGGCATTGGTATTGTGGTTGTCCATACGACATGTGTGACACCATGAACCGGTGGTTTGACACTGACGACGCATATGTCCCATGCCGCCACAGCGATAGCATTTGATGCGACCAGCAACAGGCGATCTCGGGCCTAAATTTGTTATCGCAGACGCTTGTAAGGATGCAAGAGCTGCTAGCACTTGATTGTGAGAGGCCTCAGCTTGCGCCTTTAAACTTGCCCCTAACTCCTTAATAGCCTCAATCAGAAATGCTTGGGGCCCGACTGGCACGCTTGATagcttttccagtgcctcttcAATAGTCCAATTACTCCTCAAAGTACTCAGAGTACTACGTGCTGTTGAATTACAATTTTGGAGCGCGCATTGTTTTAACATTACTCCTCTCATATACTCTGGCACCCCTGCTTTTTCAATAGCCCCGGCTACCTTATCTATGAATGCCCCAAAGTCCTCATCTCTACCTTGTCGGATCCCCATATAAAATGGCAATCCATCAGGCACCTTAATCTTGTCCATGGCCTGTCTAGCTAAATACATCGTTTCTCGACATTTATCTGGCCCTAATAATGCTTGGGCTTGTGTTCTGAAAAAAGGCCCTAGCCCTAAGAGTTCTTCGATAGTTACACCATGTAGTGGGTCTCCCGGCTGCCTAGCCTTTGAGACACTCTGATGGCACAGTTCTTGCCAATATGCattaaacaacagctgttgATGTTGTGAAGAGATCAATTTTGCTATTGCCCGACAATCGGATGGCAGCAATATCTGCGTACTCCAAATATAATCCAATATCTGCTTAGCTGGCTCGCTTTTTACCCCAAACTGACTAACTGTAGATCGTAGCTGCGATAATAATTTCCAATCTAAAGCTGTGATGGTGGCCTGCATCCCTCCCGCAGGATTAGAGGCATATATCACTGGAAACGCCATGTGCCGCACGGCCTCC
```

### Identifying homologs in novel proteins from transcriptome

- See this example: https://github.com/cfarkas/annotate_my_genomes/wiki#5-identifying-homologs-in-novel-proteins-from-transcriptome

## VII) Nextflow

- Nextflow (https://www.nextflow.io/) is a great workflow framework and a programming DSL that eases the writing of data-intensive computational pipelines. We encourage and support the usage of this framework across different platforms for reproducibility. 

### Requirements: 

- Nextflow can be installed as depicted here (https://www.nextflow.io/) or via anaconda:

```
conda install -c bioconda nextflow
```
- ncbi-blast+ version equal or higher than v2.7.1. To install it, see here: https://github.com/cfarkas/annotate_my_genomes/wiki#5-installing-up-to-date-ncbi-blast-version-v271

### Installation and Complete pipeline:

1) Install annotate_my_genomes (skipping conda install):
```
git clone https://github.com/cfarkas/annotate_my_genomes.git   # clone repository
cd annotate_my_genomes                                         # enter repository
bash makefile.sh                                               # make & install
sudo cp ./bin/* /usr/local/bin/                                # required to run nextflow scripts (requires sudo privileges)
```

2) Enter into nextflow_scripts subdirectory and run the pipeline using --flags parameters as follows:
```
cd nextflow_scripts/
```
2.1) genome-download (i.e : galGal6 genome, ouputting in current directory)
```
nextflow run genome-download.nf --genome galGal6 --conda /path/to/annotate_my_genomes/environment.yml --outdir ./
```
2.1) annotate-my-genomes
```
nextflow run annotate-my-genomes.nf --stringtie /path/to/stringtie.gtf --ref_annotation /path/to/galGal6.gtf --genome /path/to/galGal6.fa --config /path/to/annotate_my_genomes/gawn_config.sh --threads 50 --conda /path/to/annotate_my_genomes/environment.yml --outdir /path/to/output_folder/
```
2.2) add-ncbi-annotation
```
nextflow run add-ncbi-annotation.nf --stringtie /path/to/stringtie.gtf --NCBI_annotation /path/to/galGal6_ncbiRefSeq.gtf --ref_annotation /path/to/galGal6.gtf --genome /path/to/galGal6.fa --config /path/to/annotate_my_genomes/gawn_config.sh --threads 50 --conda /path/to/annotate_my_genomes/environment.yml --outdir /path/to/output_folder/
```
2.3) isoform-identification (i.e.: outputting in current directory)
```
nextflow run isoform-identification.nf --NCBI_tmap /path/to/NCBI_compare.stringtie.gtf.tmap --NCBI_transcripts /path/to/NCBI_transcripts.fa --genome_name galGal6 --conda /path/to/annotate_my_genomes/environment.yml --outdir ./
```
To run nextflow parameters in each case, full paths to files are mandatory

### More Scenarios?

- For downstream analysis and examples, please visit our wiki page : https://github.com/cfarkas/annotate_my_genomes.wiki.git

### Notes
Compiling automatically uses Shell script compiler shc to make binaries, please check: https://github.com/neurobin/shc.
