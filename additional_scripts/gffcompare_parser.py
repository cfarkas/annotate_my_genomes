#!/usr/bin/env python
# coding: utf-8

### Libraries
import sys
import argparse

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    OKRED = '\033[91m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

parser = argparse.ArgumentParser(description="This script will parse gffcompare tmap file, including ncbiRefSeqLink.txt file to produce an annotated transcript csv file.")
parser.add_argument('--tmap_file', help="path to gffcompare tmap file. In example: UCSC_compare.stringtie_chr33.gtf.tmap or NCBI_compare.stringtie_chr33.gtf.tmap", type=str)
parser.add_argument('--transcripts_file', help="transcripts file, output of annotate-my-genomes (transcripts.fa) or add-ncbi-annotation (NCBI_transcripts.fa) programs", type=str)
parser.add_argument('--genome_name', help="Genome name in UCSC notation. In example: mm10, galGal6, hg38, rn6", type=str)
args = parser.parse_args()

name = sys.argv[0]
tmap_file = str(sys.argv[2])
transcripts_file = str(sys.argv[4])
genome_name = str(sys.argv[6])

print(bcolors.OKGREEN + "Command line:", str(sys.argv) + bcolors.ENDC)
print("")
print(bcolors.OKGREEN + "--- Loading python libraries --" + bcolors.ENDC)
print("")
import pandas as pd

df = pd.read_csv(tmap_file, sep = '\t')
print("Total number of transcripts:", df.shape[0])
print("")

df2 = df[~df.ref_id.astype(str).str.contains('-')]
novel_transcripts = df[df.ref_id.astype(str).str.contains('-')]

df3 = df2[["ref_gene_id", "ref_id", "class_code", "qry_gene_id", "qry_id", "num_exons", "FPKM", "TPM"]]
df_novel_transcripts = novel_transcripts[["ref_gene_id", "ref_id", "class_code", "qry_gene_id", "qry_id", "num_exons", "FPKM", "TPM"]]

df3.head(10)
df_novel_transcripts.head(10)

colnames=['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18']
dfA1 = pd.read_csv(ncbiRefSeqLink.txt, sep = '\t', low_memory=False, names=colnames, header=None)
dfA1.head(10)

dfA2 = dfA1[['0', '1', '2', '3', '5', '14', '16']]
dfA2.head(10)

dfA2 = dfA2.rename(columns={'0': 'ref_id', '1': 'Annotation Status', '2' : 'NCBI RefSeq Gene ID', '3' : 'Transcript Description', '5' : 'NCBI RefSeq Protein ID', '14' : 'Alternative Gene Name', '16' : 'RefSeq Transcript Info'})
dfA2.sample(10)

colnames = ['qry_id', 'cds_seq', 'none']
cds = pd.read_csv(transcripts_Isoform.tab, sep = '\t', names=colnames)
cds.head(10)
cds2 = cds[["qry_id", "cds_seq"]]
cds2.sample(10)

result1 = pd.merge(df3, dfA2, on='ref_id', how='inner')
result1.sample(10)
result2 = pd.merge(result1, cds2, on='qry_id', how='inner')
result2.sample(10)
result3 = pd.merge(df_novel_transcripts, cds2, on='qry_id', how='inner')
result3.sample(10)
print("Number of Joined Transcripts (reference):", result2.shape[0])
print("")
print("Number of Joined Transcripts (novel):", result3.shape[0])
print("")
result2.to_csv('Ref_Transcript_Annotation.csv', index=False)
result3.to_csv('Novel_Transcript_Annotation.csv', index=False)
print(bcolors.OKGREEN + "::: All Done. Ref_Transcript_Annotation.csv and Novel_Transcript_Annotation.csv were succesfully produced" + bcolors.ENDC)
print("")
