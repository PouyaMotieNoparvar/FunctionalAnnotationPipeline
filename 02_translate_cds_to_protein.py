import os
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Define paths for input and output directories
input_dir = "02_protein_coding"
output_dir = "03_translated_proteins"

# Create the output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Set up logging to save warnings to a log file
logging.basicConfig(filename='processing_warnings.log', 
                    level=logging.WARNING, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Function to check for warnings and translate sequences
def translate_sequences():
    # Iterate over all files in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith(".fasta"):  # Process only FASTA files
            input_path = os.path.join(input_dir, filename)
            output_path = os.path.join(output_dir, filename)
            
            # Read the sequences from the input file
            try:
                sequences = list(SeqIO.parse(input_path, "fasta"))
            except Exception as e:
                logging.error(f"Error reading file '{filename}': {e}")
                continue

            # Check if the file is empty
            if not sequences:
                logging.warning(f"File '{filename}' is empty.")
                continue
            
            # List to store translated sequences
            translated_records = []
            
            # Process each sequence in the file
            for record in sequences:
                # Check if the sequence starts with the start codon (ATG)
                if not record.seq.startswith("ATG"):
                    logging.warning(f"Sequence '{record.id}' in file '{filename}' does not start with the start codon (ATG).")
                
                # Check if the sequence length is a multiple of 3
                if len(record.seq) % 3 != 0:
                    logging.warning(f"Sequence '{record.id}' in file '{filename}' has a length that is not a multiple of 3.")
                
                # Translate the sequence
                protein_seq = record.seq.translate(table="Standard", to_stop=True)
                
                # Check if the translated sequence length is less than 10 amino acids
                if len(protein_seq) < 10:
                    logging.warning(f"Sequence '{record.id}' in file '{filename}' results in a protein sequence less than 10 amino acids.")
                
                # Create a SeqRecord for the translated protein
                translated_record = SeqRecord(
                    protein_seq,
                    id=record.id,
                    description="Translated protein"
                )
                translated_records.append(translated_record)
            
            # Write translated sequences to the output file
            try:
                SeqIO.write(translated_records, output_path, "fasta")
                logging.info(f"Translated sequences saved to '{output_path}'")
            except Exception as e:
                logging.error(f"Error writing to file '{output_path}': {e}")

# Run the function to translate sequences
translate_sequences()
