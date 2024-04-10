import argparse

def parse_gtf(input_gtf):
    transcript_coords = {}
    with open(input_gtf, 'r') as gtf:
        for line in gtf:
            if line.startswith('#') or line.strip() == "":
                continue  # Skip header and empty lines

            parts = line.strip().split('\t')
            if parts[2] == 'transcript':  # Process only transcript lines
                chromosome = parts[0]
                gene_id = parts[8].split('gene_id "')[1].split('"')[0]
                start, end = parts[3], parts[4]

                # Store the gene_id with its coordinates in a dictionary
                transcript_coords[gene_id] = (chromosome, start, end)
    return transcript_coords

import argparse

def parse_gtf(input_gtf):
    transcript_coords = {}
    with open(input_gtf, 'r') as gtf:
        for line in gtf:
            if line.startswith('#') or line.strip() == "":
                continue  # Skip header and empty lines

            parts = line.strip().split('\t')
            if parts[2] == 'transcript':  # Process only transcript lines
                chromosome = parts[0]
                gene_id = parts[8].split('gene_id "')[1].split('"')[0]
                start, end = parts[3], parts[4]

                # Store the gene_id with its coordinates in a dictionary
                transcript_coords[gene_id] = (chromosome, start, end)
    return transcript_coords

def update_cds_coordinates(transcript_coords, cds_gtf, updated_cds_output):
    processed_count = 0  # Counter for successfully processed lines
    skipped_count = 0  # Counter for skipped lines

    with open(cds_gtf, 'r') as cds_file, open(updated_cds_output, 'w') as output_file:
        for line in cds_file:
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue  # Skip incomplete lines

            # Extract the gene_id from the Parent attribute
            attribute_part = parts[8]
            gene_id_full = None
            for attr in attribute_part.split(';'):
                if attr.startswith('Parent='):
                    gene_id_full = attr.split('Parent=')[1]

            if gene_id_full:
                # Remove potential .p<number> part to match the gene_id format in the GTF file
                gene_id = '.'.join(gene_id_full.split('.')[:-1])

                if gene_id in transcript_coords:
                    chromosome, genomic_start, _ = transcript_coords[gene_id]
                    cds_start, cds_end = int(parts[3]), int(parts[4])

                    # Update the CDS coordinates based on the transcript's genomic start position
                    new_cds_start = int(genomic_start) + cds_start - 1
                    new_cds_end = int(genomic_start) + cds_end - 1

                    # Write the updated CDS entry to the output file
                    parts[0] = chromosome
                    parts[3], parts[4] = str(new_cds_start), str(new_cds_end)
                    output_file.write('\t'.join(parts) + '\n')

                    processed_count += 1  # Increment processed lines counter
                else:
                    print(f"Warning: gene_id '{gene_id}' not found in GTF file. Skipping CDS entry.")
                    skipped_count += 1  # Increment skipped lines counter
            else:
                print(f"Warning: No Parent attribute found for line. Skipping CDS entry.")
                skipped_count += 1  # Increment skipped lines counter

    print(f"Total lines processed: {processed_count}")
    print(f"Total lines not processed: {skipped_count}")




def main():
    parser = argparse.ArgumentParser(description='Process GTF and GFF3 files to update CDS coordinates.')
    parser.add_argument('--input_gtf', required=True, help='Input GTF file path')
    parser.add_argument('--cds_gff3', required=True, help='Input CDS GFF3 file path')
    parser.add_argument('--updated_cds_output', required=True, help='Output path for updated CDS GFF3 data')
    args = parser.parse_args()

    # Parse the GTF file to get transcript coordinates
    transcript_coords = parse_gtf(args.input_gtf)

    # Update the CDS coordinates in the GFF3 file using the parsed transcript coordinates
    update_cds_coordinates(transcript_coords, args.cds_gff3, args.updated_cds_output)

    print(f"Updated CDS coordinates have been written to {args.updated_cds_output}")

if __name__ == '__main__':
    main()
