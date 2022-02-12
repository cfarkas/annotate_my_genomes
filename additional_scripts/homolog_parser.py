#!/usr/bin/env python
# coding: utf-8

### Libraries

# conda install -c anaconda pandas
# conda install --channel conda-forge --channel bioconda pybedtools
# pip install xlrd==1.2.0

import sys
import csv
import argparse
import pandas as pd
import numpy as np
import pybedtools

np.random.seed(5)

parser = argparse.ArgumentParser(description="This script will annotate and add genomic coordinates to Novel proteins not included in reference GTF annotations.")
parser.add_argument('--Ref', help="Ref_Transcript_Annotation.csv file, obtained with isoform-identification pipeline", type=str)
parser.add_argument('--blastp', help="parsed_results.tab file, blastp result from Novel proteins against UniProt database", type=str)
parser.add_argument('--bed', help="final_annotated.format.bed, which correspond to final_annotation.gtf file as BED format", type=str)
parser.add_argument('--eggnog', help="out.emapper.annotations.xlsx file, which correspond to eggNOG-mapper annotations", type=str)
args = parser.parse_args()

name = sys.argv[0]
REFERENCE = str(sys.argv[2])
BLASTP = str(sys.argv[4])
BED = str(sys.argv[6])
EGGNOG = str(sys.argv[8])

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


print(bcolors.OKGREEN + "1 ::: Reading Ref_Transcript_Annotation.csv ::: " + bcolors.ENDC)
print("")


df1 = pd.read_csv(REFERENCE, sep = ',')
df1 = df1.rename(columns={'qry_id': 'transcript'})
print("Number of transcripts matching reference:", df1.shape[0])
print("")


df1.head(5)

print(bcolors.OKGREEN + "2 ::: Reading parsed_results.tab ::: " + bcolors.ENDC)
print("")


df2 = pd.read_csv(BLASTP, sep = '\t')
print("Number of novel proteins sharing >=90% identity with uniprot:", df2.shape[0])
print("")


df2.head(5)


print(bcolors.OKGREEN + " ::: formatting columns of blastp_results :::" + bcolors.ENDC)
print("")

df2[['hit_id','Transcript Description']] = df2['hit_id'].str.split(" ", 1 ,expand=True)
df2


df2[['Transcript Description', 'organism']] = df2['Transcript Description'].str.split('OS=',expand=True)
df2

print(df2.head(5))

blastp_results = df2.iloc[:, 0:7]
blastp_results.columns = ['transcript', 'hit_id', 'percentage_identity', 'query_length', 'alignment_length', 'e_value', 'Transcript Description']


print(blastp_results.head(5))


print(bcolors.OKGREEN + "3 ::: Reading final_annotated.format.bed ::: " + bcolors.ENDC)
print("")


df3 = pd.read_csv(BED, sep = '\t')
print("Number of trascripts in bed format:", df3.shape[0])
print("")
df3.columns = ['chr', 'start', 'end', 'transcript', 'gene']


print(df3.head(5))


print(bcolors.OKGREEN + "4 ::: Adding genomic coordinates to blastp_results ::: " + bcolors.ENDC)
print("")


result1 = pd.merge(blastp_results, df3, on='transcript', how='inner')
result1


print(bcolors.OKGREEN + "5 ::: Adding genomic coordinates to Novel proteins and Ref_Transcript_Annotation ::: " + bcolors.ENDC)
print("")


result2 = pd.merge(df1, df3, on='transcript', how='inner')
result2


Novel_protein_hits = result1.loc[:, ['chr', 'start', 'end', 'transcript', 'Transcript Description', 'percentage_identity']]
print("Novel protein hits:")
print(Novel_protein_hits.head(5))
print("")

Reference_annotation = result2.loc[:, ['chr', 'start', 'end', 'transcript', 'Transcript Description', 'NCBI RefSeq Gene ID']]
print("Reference annotations:")
print(Reference_annotation.head(5))
print("")

print(bcolors.OKGREEN + " 6 ::: Creating bed files with annotations ::: " + bcolors.ENDC)
print("")
a = pybedtools.BedTool.from_dataframe(Novel_protein_hits)
b = pybedtools.BedTool.from_dataframe(Reference_annotation)
a.saveas('Novel_protein_with_coordinates.bed')
b.saveas('Reference_annotation_with_coordinates.bed')

print(bcolors.OKGREEN + " 7 ::: Parsing eggNOG-mapper annotations and intersect with blastp results::: " + bcolors.ENDC)
print("")

df1 = pd.read_excel(EGGNOG)

df1.head(5)

df1 = df1.rename(columns={'query': 'transcript'})

intersect1 = pd.merge(df1, blastp_results, on='transcript', how='inner')
Novel_protein_hits_coords = result1.loc[:, ['chr', 'start', 'end', 'transcript']]
intersect2 = pd.merge(Novel_protein_hits_coords, intersect1, on='transcript', how='inner')
intersect2 = intersect2.rename(columns={'hit_id': 'blastp_hit_id', 'percentage_identity': 'blastp_percentage_identity', 'e_value': 'blastp_e_value', 'Transcript Description': 'blastp_transcript_description'})
intersect2.to_csv('eggNOG-mapper-blastp-intersections.csv')

print(intersect2.head(5))
print("")

print("filtering table with intersections with assigned Gene IDs: ")
intersect3 = intersect2.replace('-', np.nan)
intersect3 = intersect3.dropna(subset=['Preferred_name'])
c = pybedtools.BedTool.from_dataframe(intersect3)
c.saveas('eggNOG-mapper-blastp-with_coordinates.bed')

print("")
print(bcolors.OKGREEN + "Novel_protein_with_coordinates.bed ==> corresponds to a BED file containing novel proteins with mapped coordinates" + bcolors.ENDC)
print("")
print(bcolors.OKGREEN + "Reference_annotation_with_coordinates.bed ==> corresponds to a BED file containing annotated proteins with mapped coordinates" + bcolors.ENDC)
print("")
print(bcolors.OKGREEN + "eggNOG-mapper-blastp-intersections.csv ==> corresponds to a csv file containing eggNOG-mapper + blastp intersections" + bcolors.ENDC)
print("")
print(bcolors.OKGREEN + "eggNOG-mapper-blastp-with_coordinates.bed ==> corresponds to a BED file containing eggNOG-mapper + blastp intersections" + bcolors.ENDC)
print("")
print("All Done")
print("")
