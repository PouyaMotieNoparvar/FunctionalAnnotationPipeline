from Bio import Entrez, SeqIO

# Configure your email to use NCBI Entrez
Entrez.email = "p.motienoparvar1@universityofgalway.ie"  # Replace with your email

# Define the organism name and search for its CDS records
organism = "Agaricus bisporus"
search_term = f"{organism}[Organism] AND cds[Feature]"

def fetch_cds_sequences(organism_name, output_file="agaricus_bisporus_cds.fasta"):
    try:
        # Search for CDS sequences in NCBI
        print("Searching for CDS records...")
        handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=10000)
        record = Entrez.read(handle)
        handle.close()

        # Get the list of IDs
        ids = record["IdList"]
        print(f"Found {len(ids)} CDS records for {organism_name}.")

        if len(ids) == 0:
            print("No records found.")
            return

        # Fetch the sequences
        print("Fetching CDS sequences...")
        handle = Entrez.efetch(db="nucleotide", id=ids, rettype="fasta", retmode="text")
        sequences = handle.read()
        handle.close()

        # Write to output file
        with open(output_file, "w") as f:
            f.write(sequences)

        print(f"CDS sequences saved to {output_file}.")
    except Exception as e:
        print(f"An error occurred: {e}")

# Call the function
fetch_cds_sequences(organism)
