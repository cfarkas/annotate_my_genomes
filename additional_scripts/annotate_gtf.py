import sys
from tqdm import tqdm

# Define file names from command line arguments
gtf_file = sys.argv[1]  # GTF file
hits_file = sys.argv[2]  # Hits file
annotation_table_file = sys.argv[3]  # Transcriptome annotation table file
output_file = sys.argv[4]  # Output GTF file with annotations

# Function to append blastx information and gawn_name to GTF entries
def append_annotations_to_gtf(gtf_line, blastx_info, gawn_names):
    if 'transcript' in gtf_line:
        transcript_id = gtf_line.split('transcript_id "')[1].split('"')[0]
        annotations = []
        if transcript_id in blastx_info:
            annotations.append(f'blastx "{blastx_info[transcript_id]}";')  # Added semicolon here
        if transcript_id in gawn_names:
            annotations.append(f'gawn_annotation "{gawn_names[transcript_id]}";')  # Added semicolon here
        if annotations:
            gtf_line += ' ' + ' '.join(annotations)
    return gtf_line

# Read the hits file and store the blastx info in a dictionary
blastx_info = {}
with open(hits_file, 'r') as hits:
    for line in hits:
        transcript_id, blastx_id = line.strip().split(' ')
        blastx_info[transcript_id] = blastx_id

# Read the annotation table and store the gawn_names in a dictionary
gawn_names = {}
with open(annotation_table_file, 'r') as table:
    for line in table:
        parts = line.strip().split('\t')
        if len(parts) > 2:
            transcript_id = parts[0]
            gawn_name = parts[2]
            gawn_names[transcript_id] = gawn_name

# Read the GTF file, modify entries with annotations, and write to the output file
with open(gtf_file, 'r') as gtf, open(output_file, 'w') as out_gtf:
    for line in tqdm(gtf, desc="Annotating GTF"):
        modified_line = append_annotations_to_gtf(line.strip(), blastx_info, gawn_names)
        out_gtf.write(modified_line + '\n')

print("Annotation completed. Output is in", output_file)
