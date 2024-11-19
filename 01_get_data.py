import os
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML
import re

# Set up your email for Entrez
Entrez.email = "p.motienoparvar1@universityofgalway.ie"  # Replace with your email

def read_cds_from_fasta(filepath):
    """Read the CDS sequence from a specified FASTA file."""
    with open(filepath, "r") as file:
        # Read the FASTA file
        records = SeqIO.parse(file, "fasta")
        # Assume the FASTA file contains a single CDS sequence per file
        cds_seq = next(records).seq  # Get the first sequence from the file
    return str(cds_seq)

def perform_blast_search(cds_seq):
    """Perform BLAST search and retrieve accessions of all sequences."""
    blast_results = NCBIWWW.qblast("blastn", "nt", cds_seq)
    blast_records = NCBIXML.read(blast_results)
    
    accessions = []
    for alignment in blast_records.alignments:
        for hsp in alignment.hsps:
            # Remove the condition that filters out "Arabidopsis"
            acc_match = re.search(r'[A-Z]{1,2}_\d+\.\d+', alignment.title)
            if acc_match:
                accessions.append(acc_match.group(0))
    return accessions

def fetch_cds_sequences(accessions):
    """Fetch CDS regions based on GenBank annotations."""
    cds_sequences = []
    for acc in accessions:
        try:
            handle = Entrez.efetch(db="nucleotide", id=acc, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            
            for feature in record.features:
                if feature.type == "CDS":
                    cds_seq = feature.location.extract(record.seq)
                    cds_sequences.append((acc, str(cds_seq)))
        except Exception as e:
            print(f"Error fetching sequence for {acc}: {e}")
    return cds_sequences

def save_to_fasta(sequences, output_file):
    """Save sequences to a FASTA file."""
    with open(output_file, "w") as file:
        for acc, seq in sequences:
            file.write(f">{acc}\n{seq}\n")

# Main function to process all files in a directory
def process_all_files(input_folder, output_cds_folder):
    """Process all FASTA files in the input_folder and save results in the output folder."""
    if not os.path.exists(output_cds_folder):
        os.makedirs(output_cds_folder)

    for filename in os.listdir(input_folder):
        if filename.endswith(".fasta") or filename.endswith(".fa"):
            filepath = os.path.join(input_folder, filename)
            cds_seq = read_cds_from_fasta(filepath)
            
            # Perform BLAST and get accessions from all organisms
            accessions = perform_blast_search(cds_seq)
            print(f"Found {len(accessions)} sequences for {filename}.")
            
            # Fetch CDS sequences based on accessions
            cds_sequences = fetch_cds_sequences(accessions)
            
            # Define the output filename for CDS sequences
            cds_output_filename = f"{os.path.splitext(filename)[0]}_ortho_cds.fasta"
            
            # Save CDS sequences to the output FASTA file
            cds_output_filepath = os.path.join(output_cds_folder, cds_output_filename)
            save_to_fasta(cds_sequences, cds_output_filepath)
            print(f"CDS sequences saved to {cds_output_filepath}")

# Define folder paths and run the processing
if __name__ == "__main__":
    input_folder = "01_genes_to_study"  # Folder with FASTA files as input
    output_cds_folder = "02_orthogroups_cds"  # Folder to save output CDS sequences
    process_all_files(input_folder, output_cds_folder)

