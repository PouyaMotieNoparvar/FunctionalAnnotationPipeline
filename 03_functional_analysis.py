import os
import requests
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# Function to get GO annotations from UniProt using protein accession
def get_go_terms(uniprot_id):
    # UniProt API URL for retrieving GO terms for a given protein
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.xml"
    response = requests.get(url)
    
    if response.status_code == 200:
        # Parse the XML response to extract GO terms
        go_terms = []
        from xml.etree import ElementTree
        tree = ElementTree.ElementTree(ElementTree.fromstring(response.content))
        root = tree.getroot()
        
        # Extract GO terms from the XML structure
        for entry in root.findall(".//{http://uniprot.org/uniprot}entry"):
            for go_term in entry.findall(".//{http://uniprot.org/uniprot}dbReference"):
                if 'GO' in go_term.attrib.get('type', ''):
                    go_terms.append(go_term.attrib.get('id', ''))
        return go_terms
    else:
        return []

# Function to perform BLAST and parse results
def blast_sequence(sequence):
    result_handle = NCBIWWW.qblast("blastp", "nr", sequence)
    blast_records = NCBIXML.parse(result_handle)
    result = []
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                result.append({
                    "query": blast_record.query,
                    "title": alignment.title,
                    "length": alignment.length,
                    "e_value": hsp.expect,
                    "identity": hsp.identities,
                    "query_from": hsp.query_start,
                    "query_to": hsp.query_end,
                    "subject_from": hsp.sbjct_start,
                    "subject_to": hsp.sbjct_end,
                    "alignment": hsp.sbjct,
                    "uniprot_id": alignment.title.split('|')[1]  # Extract UniProt ID from the title
                })
    return result

# Function to process a FASTA file and perform BLAST on each sequence
def process_fasta_file(fasta_file_path, output_dir):
    output_lines = []
    file_name = os.path.basename(fasta_file_path)
    output_file_name = file_name.replace(".fasta", "_blast_go_results.txt")
    
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Path to save the output file
    output_file_path = os.path.join(output_dir, output_file_name)

    with open(output_file_path, "w") as output_file:
        # Write the header
        output_file.write("Query Sequence\tTitle\tLength\tE-value\tIdentity\tQuery Start\tQuery End\t"
                          "Subject Start\tSubject End\tAlignment\tUniProt ID\tGO Terms\n")

        # Read the sequences from the FASTA file
        for record in SeqIO.parse(fasta_file_path, "fasta"):
            sequence = str(record.seq)
            results = blast_sequence(sequence)
            
            # Write the BLAST results for each sequence
            for result in results:
                # Get GO terms from UniProt
                go_terms = get_go_terms(result['uniprot_id'])
                go_terms_str = ', '.join(go_terms) if go_terms else 'No GO terms found'

                # Write the results to the file
                output_file.write(f"{result['query']}\t{result['title']}\t{result['length']}\t{result['e_value']}\t"
                                  f"{result['identity']}\t{result['query_from']}\t{result['query_to']}\t"
                                  f"{result['subject_from']}\t{result['subject_to']}\t{result['alignment']}\t"
                                  f"{result['uniprot_id']}\t{go_terms_str}\n")

    print(f"BLAST and GO results written to {output_file_path}")

# Main function to loop over all FASTA files in a given directory
def process_directory(directory_path, output_dir):
    # Loop through each file in the directory
    for file_name in os.listdir(directory_path):
        if file_name.endswith(".fasta"):
            fasta_file_path = os.path.join(directory_path, file_name)
            process_fasta_file(fasta_file_path, output_dir)

# Specify the directory containing the FASTA files and the output directory
input_directory = "./03_translated_proteins"  # Modify this path to your folder location
output_directory = "./04_GO_BLAST"  # Output directory for GO and BLAST results

process_directory(input_directory, output_directory)
