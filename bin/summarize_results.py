import argparse
import json
import pandas as pd

parser = argparse.ArgumentParser(description="Summarize PFMI3DSC results and generate an HTML table")
parser.add_argument('--result_json', required=True, help='Input JSON from pfmi3dsc_score.py')
parser.add_argument('--out', required=True, help='Output HTML file')
parser.add_argument('--query_accid', required=True, help='UniProt ACCID of the query protein')
args = parser.parse_args()

with open(args.result_json) as f:
    results = json.load(f)

# Load main alignment matrix (residues as index)
alignment = pd.DataFrame(results['alignment_matrix_seq'])
if "scores" not in alignment.columns:
    alignment["scores"] = 0
if "probability" not in alignment.columns:
    alignment["probability"] = 1.0

# Always use the provided query_accid
query_protein = args.query_accid

# List of proteins (columns)
protein_cols = [c for c in alignment.columns if c not in ['scores', 'probability', 'Predicted Functional']]
# Move query protein to the first column if present
if query_protein in protein_cols:
    protein_cols = [query_protein] + [c for c in protein_cols if c != query_protein]


# Calculate thresholds (same logic as used in pfmi3dsc_score.py)
num_proteins = len(protein_cols)
threshold = 0.01 / (num_proteins - 1) if num_proteins > 1 else 0
threshold1 = num_proteins / 2

# Mark predicted functional residues
alignment['Predicted Functional'] = (alignment['probability'] < threshold) & (alignment['scores'] > threshold1)
alignment['Predicted Functional'] = alignment['Predicted Functional'].map({True: "★", False: ""})

# Reorder columns
alignment = alignment.reset_index().rename(columns={'index': 'Residue'})  # Ensure residue index is a column
if 'pos' in alignment.columns:
    alignment = alignment.rename(columns={'pos': 'Residue'})
if 'Residue' not in alignment.columns:
    alignment.insert(0, 'Residue', alignment.index + 1)
columns_order = ['Residue'] + protein_cols + ['scores', 'probability', 'Predicted Functional']
alignment = alignment[columns_order]

# Build HTML
html_parts = []
html_parts.append("<html><head><title>PFMI3DSC Functional Residue Table</title></head><body>")
html_parts.append("<h1>PFMI3DSC Predicted Functional Residues</h1>")
html_parts.append(f"<h2>Query Protein: {query_protein}</h2>")
html_parts.append(f"<p><b>Threshold for probability:</b> {threshold:.4g}<br>")
html_parts.append(f"<b>Threshold for score:</b> {threshold1:.2f}</p>")
html_parts.append(alignment.to_html(index=False, border=1, justify='center'))
html_parts.append("</body></html>")

with open(args.out, "w") as f:
    f.write('\n'.join(html_parts))

print(f"HTML table report written to {args.out}")
