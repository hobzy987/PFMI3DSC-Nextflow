# PFMI3DSC-Nextflow

**Protein Functional Mutation Identification by 3D Structure Comparison of Protein Families (PFMI3DSC)**

A reproducible Nextflow pipeline for predicting functional mutation hotspots in proteins, leveraging structural alignments, cancer genomics databases, and AlphaFold models.

## Quick Start

- Place your required input data in the `data/` folder.{the cleand file of protein families, mutations and hotspots is uploaded}
- Run the pipeline using Nextflow.
- Outputs will appear in the `results/` folder.

## For full documentation, see the comments in each file and [under review](#).

# PFMI3DSC-Nextflow Workflow Documentation

## Overview

This Nextflow pipeline automates the PFMI3DSC functional residue prediction workflow. It takes a list of UniProt accession IDs, fetches family and structure data, runs structure alignment, parses mutations, scores functional residues, and generates an HTML summary report for each query protein.

---

## Workflow Steps

1. **get_family**  
   - **Input:** UniProt ACCID  
   - **Output:** `family_${accid}.json`  
   - **Description:** Retrieves family information for the given protein accession.

2. **download_pdbs**  
   - **Input:** `(accid, family_json)`  
   - **Output:** `pdbs_${family_json.baseName}.txt`  
   - **Description:** Downloads PDB structures for the protein family.

3. **run_foldseek**  
   - **Input:** `(accid, pdbs_txt)`  
   - **Output:** `alignments_${pdbs_txt.baseName}.tsv`  
   - **Description:** Runs Foldseek to align structures.

4. **parse_mutations**  
   - **Input:** `(accid, family_json)`  
   - **Output:** `mutations_${family_json.baseName}.json`  
   - **Description:** Parses mutation data for the family.

5. **pfmi3dsc_score**  
   - **Input:**  
     - `(accid, alignment_file)`  
     - `(accid2, mut_json)`  
   - **Output:** `results_${alignment_file.baseName}.json`  
   - **Description:** Scores functional residues using PFMI3DSC. Only runs when `accid == accid2`.

6. **summarize_results**  
   - **Input:** `(accid, result_json)`  
   - **Output:** `results/summary_${result_json.baseName}.html`  
   - **Description:** Generates an HTML summary report for each query protein.  
   - **Output Location:** The HTML file is automatically copied to the `results/` directory in the project root.

---

## Output

- **HTML summary reports** for each query protein are saved in the `results/` directory at the root of your project:
  ```
  ./results/summary_results_alignments_pdbs_family_<ACCID>.html
  ```

---

## Cleaning Up

Before rerunning the workflow, you can clean up all intermediate and Nextflow metadata files by running:
```bash
rm -rf work/ .nextflow*
```
from your project root.

---

## Usage

Run the workflow with:
```bash
nextflow run main.nf --accids <ACCID1,ACCID2,...>
```
Replace `<ACCID1,ACCID2,...>` with your comma-separated list of UniProt accession IDs.

---

## Notes

- All processes use relative paths for outputs, and Nextflow's `publishDir` ensures only the HTML summary files are copied to the root `results/` directory.
- The workflow is written in Nextflow DSL2.
- Each process is modular and can be adapted or extended as needed.
- The protein family information is parsed using the `similar.txt` file from UniProtKB. The full file can be accessed at:  
  [https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/docs/similar.txt](https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/docs/similar.txt)
- The 3D protein structures used in this analysis are obtained from the AlphaFold Protein Structure Database, and the version currently integrated is **AlphaFold v4**, which offers improved accuracy over previous versions.

---
