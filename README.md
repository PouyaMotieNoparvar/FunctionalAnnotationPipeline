
---

# Functional Annotation Pipeline

This pipeline provides a set of scripts designed to functionally annotate genes based on homology for genes that are not already annotated.

## Prerequisites

Before using the pipeline, it is recommended to create a virtual Python environment and install the required packages: 

- **biopython**
- **requests**

You can install them using the following command:

```bash
pip install biopython requests
```

## Tutorial: How to Use the Scripts

In the main directory, follow the steps below to run the pipeline:

### Step 1: Retrieve RNA Sequences

Use the `00_retrieve_whole_genes.py` script to retrieve RNA sequences for a specific organism (e.g., *Agaricus bisporus*).

**Important Notice**:  
You will need to manually remove large files as they are typically chromosome sequences, which may cause errors during processing.

In the script, set your email address and adjust the number of output sequences (e.g., `retmax=10000`). This will fetch 10,000 sequences into a single FASTA file named `agaricus_bisporus_cds.fasta`. You can change the name of this file in the script if necessary.

### Step 2: Split Sequences into Separate Files

Use the `00_split_sequences.py` script to split the sequences into separate files for each individual sequence.  
You can modify the script to work with a single file if needed, though the current setup is optimized for working with multiple files.

### Step 3: Retrieve Orthologue Genes via BLAST

Run the `01_get_data.py` script to retrieve orthologous gene sequences by performing a BLAST search for each gene you retrieved in the previous step. The script filters the results to keep only CDS (coding sequence) sequences.

### Step 4: Filter Sequences to Keep Protein-Coding Sequences

Use the `01_filter_out.py` script to remove sequences whose lengths are not multiples of three, ensuring that only protein-coding sequences are retained.

### Step 5: Perform Functional Analysis

Finally, run the `03_functional_analysis.py` script. This will generate a text file containing the Gene Ontology (GO) terms and protein BLAST results. The output will look like this:

```
Query Sequence    Title            Length    E-value    Identity    Query Start    Query End    Subject Start    Subject End    Alignment    UniProt ID    GO Terms
ATGGCGTG...      sp|Q9XZI2|MT1A_HUMAN    120    1e-50    98%    1    120    1    120    MGRPQST...    Q9XZI2    GO:0005737, GO:0004672
```

---

