import argparse
import json
import ast
import pandas as pd
import os

parser = argparse.ArgumentParser(description="Find all family members for a given UniProt ACCID")
parser.add_argument('--accid', required=True, help='UniProt ACCID (e.g., P01112)')
parser.add_argument('--data', required=True, help='Path to data directory')
parser.add_argument('--out', required=True, help='Output JSON filename')
args = parser.parse_args()

csv_file = os.path.join(args.data, 'final_database.csv')
if not os.path.exists(csv_file):
    raise FileNotFoundError(f"Cannot find {csv_file}")

df = pd.read_csv(csv_file)

# Find row(s) where 'uniprot_ACCID' contains the query ACCID
mask = df['uniprot_ACCID'].str.contains(args.accid)
if not mask.any():
    raise ValueError(f"ACCID {args.accid} not found in final_database.csv.")

row = df[mask].iloc[0]
# Parse family members and gene names from their stored list form
family_name = row['families']
try:
    uniprot_ids = ast.literal_eval(row['uniprot_ACCID'])
    gene_names = ast.literal_eval(row['gene name'])
except Exception:
    uniprot_ids = [row['uniprot_ACCID']]
    gene_names = [row['gene name']]

result = {
    "input_accid": args.accid,
    "family_name": family_name,
    "family_members": uniprot_ids,
    "gene_names": gene_names
}

with open(args.out, "w") as f:
    json.dump(result, f, indent=2)

print(f"Found {len(uniprot_ids)} family members for {args.accid}.")
