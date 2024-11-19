import os
from Bio import SeqIO

# Input file path
input_file = "agaricus_bisporus_sequences.fasta"  # Replace with the path to your input FASTA file

# Output directory
output_dir = "./01_genes_to_study"

# Create the output directory if it does not exist
os.makedirs(output_dir, exist_ok=True)

# Read the input FASTA file and write each sequence to a separate FASTA file
with open(input_file, "r") as input_handle:
    sequences = SeqIO.parse(input_handle, "fasta")
    
    for seq_record in sequences:
        # Get the sequence's ID (this will be used as the filename)
        seq_id = seq_record.id
        
        # Define the output filename for each sequence in the desired directory
        output_filename = os.path.join(output_dir, f"{seq_id}.fasta")
        
        # Write the individual sequence to a new FASTA file
        with open(output_filename, "w") as output_handle:
            SeqIO.write(seq_record, output_handle, "fasta")

        print(f"Created: {output_filename}")
