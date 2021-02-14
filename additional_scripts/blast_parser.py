#!/usr/bin/env python

# Takes a .xml formatted blast results file as input and prints the query and hit ids
# for sequences passing the thresholds passed via the command line arguments. For sequences
# with no hits below the thresholds, the program returns "no hits below threshold" rather
# than the hit id.

import getopt, sys
from Bio import SeqIO
from Bio.Blast import NCBIXML

## Function to parse an XML format BLAST results file.

def parse_results(result_file, e_val_thresh, ident_thresh, align_thresh):
        result_handle = open(result_file, 'r')  ## The XML file to parse.
        blast_records = NCBIXML.parse(result_handle)
        print('query_id\thit_id\tpercentage_identity\tquery_length\talignment_length\te_value')

        for record in blast_records:  ## Loop through each query.
                query_id = record.query
                if len(record.alignments) > 0:  ## Check whether there are hits.
                        e_val = record.alignments[0].hsps[0].expect
                        if e_val < e_val_thresh:  ## Is hit below E-value?
                                tot_ident = sum([hsp.identities for hsp in record.alignments[0].hsps])  ## Sum of all identities for all hsps.
                                query_len = record.query_length  ## Length of query
                                align_len = sum([hsp.align_length for hsp in record.alignments[0].hsps])  ## Length of query alignment to hit.
                                pct_ident = tot_ident/float(align_len)*100  ## Calculates percentage identity.
                                top_hit = record.alignments[0].hit_id + record.alignments[0].hit_def
                                if pct_ident > ident_thresh:  ## Checks whether above percentage identity cutoff.
                                        if align_len > align_thresh:
                                                print('%s\t%s\t%f\t%i\t%i\t%s' % (query_id, top_hit, pct_ident, query_len, align_len, str(e_val)))
                                        else:
                                                print('%s\t%s\t%s\t%s\t%s\t%s' % (query_id, '', '', '', '', ''))
                                else:
                                        print('%s\t%s\t%s\t%s\t%s\t%s' % (query_id, '', '', '', '', ''))
                        else:
                                print('%s\t%s\t%s\t%s\t%s\t%s' % (query_id, '', '', '', '', ''))
                else:
                        print('%s\t%s\t%s\t%s\t%s\t%s' % (query_id, '', '', '', '', ''))

        result_handle.close()

## How to use this.

def usage():
        print("""
\nblast_parser.py.\n
Takes a .xml formatted blast results file as input and prints the query and hit ids
for sequences passing the thresholds passed via the command line arguments. For sequences
with no hits below the thresholds, the program returns "no hits below threshold" rather
than the hit id.\n
Basic usage:
\tpython blast_parser.py -i <results.xml> -e 1e-20 -p 97 -a 100 > parsed_results.txt\n
Arguments:
\t-h, --help\t\t\tPrint this information.
\t-i, --in <results.xml>\t\tXML format BLAST results file.
\t-e, --evalue <number>\t\tExpect value.
\t-p, --pct_ident <number, 0-100>\t\tPercentage identity cutoff.
\t-a, --align_len <number>\t\t Minimum alignment length.
""")

## The main program.

def main():
        try:  ## Parses the command line arguments.
                opts, args = getopt.getopt(sys.argv[1:], 'e:i:p:a:h', ['evalue=', 'in=', 'pct_ident=', 'align_len=', 'help'])
        except getopt.GetoptError:
                usage()
                sys.exit(2)

        ## Creates variables from the arguments.

        for opt, arg in opts:
                if opt in ('-e', '--evalue'):
                        e_val_thresh = float(arg)
                elif opt in ('-p', '--pct_ident'):
                        ident_thresh = float(arg)
                elif opt in ('-a', '--align_len'):
                        align_thresh = float(arg)
                elif opt in ('-i', '--in'):
                        result_file = arg
                elif opt in ('-h', '--help'):
                        usage()
                        sys.exit(0)
                else:
                        usage()
                        sys.exit(2)

        try:  ## Tries to parse the results file.
                parse_results(result_file, e_val_thresh, ident_thresh, align_thresh)
        except:  ## Otherwise, shows usage.
                sys.exit(1)

if __name__ == "__main__":
    main()
