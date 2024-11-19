import os
from Bio import SeqIO

# Define the input and output directories
input_folder = '02_orthogroups_cds'
output_folder = '02_protein_coding'

# Create output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Iterate over all the FASTA files in the input folder
for file_name in os.listdir(input_folder):
    if file_name.endswith('.fasta') or file_name.endswith('.fa'):  # check for FASTA files
        input_file_path = os.path.join(input_folder, file_name)
        
        # Prepare the output file path
        output_file_path = os.path.join(output_folder, file_name)

        # Open the input file and the output file
        with open(input_file_path, 'r') as input_handle, open(output_file_path, 'w') as output_handle:
            # Parse the FASTA file and filter sequences whose length is a multiple of 3
            sequences = []
            for record in SeqIO.parse(input_handle, 'fasta'):
                if len(record.seq) % 3 == 0:
                    sequences.append(record)
            
            # Write the filtered sequences to the output file
            SeqIO.write(sequences, output_handle, 'fasta')

        print(f"Processed {file_name} and saved to {output_file_path}")
