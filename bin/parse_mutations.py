import argparse
import json
import ast
import pandas as pd
import os

parser = argparse.ArgumentParser(description="Extract mutation/hotspot profiles for a protein family from final_database.csv")
parser.add_argument('--family_json', required=True, help='Path to family JSON file from get_family.py')
parser.add_argument('--data', required=True, help='Path to data directory (must contain final_database.csv)')
parser.add_argument('--out', required=True, help='Output JSON filename')
args = parser.parse_args()

# Load family JSON (get family_members and gene_names)
with open(args.family_json) as f:
    fam = json.load(f)
uniprot_ids = fam["family_members"]
gene_names = fam.get("gene_names", [])

# Load database
csv_file = os.path.join(args.data, "final_database.csv")
df = pd.read_csv(csv_file)

# Find the row for this family (matching one of the family_members)
row = None
for uid in uniprot_ids:
    mask = df['uniprot_ACCID'].str.contains(uid)
    if mask.any():
        row = df[mask].iloc[0]
        break

if row is None:
    raise ValueError("No family row found for provided family members.")

# Extract mutation and hotspot profiles from pre-combined columns
try:
    biomuta_list = ast.literal_eval(row['biomuta'])
    hotspot_list = ast.literal_eval(row['hotspot'])
except Exception as e:
    raise ValueError(f"Error parsing biomuta/hotspot columns: {e}")

# Map each UniProt ID to its mutation/hotspot profile
biomuta_profile = {uid: mut for uid, mut in zip(uniprot_ids, biomuta_list)}
hotspot_profile = {uid: hot for uid, hot in zip(uniprot_ids, hotspot_list)}

# Output
outdata = {
    "uniprot_ids": uniprot_ids,
    "gene_names": gene_names,
    "biomuta_profile": biomuta_profile,
    "hotspot_profile": hotspot_profile
}
with open(args.out, "w") as f:
    json.dump(outdata, f, indent=2)

print(f"Wrote mutation/hotspot profiles for {len(uniprot_ids)} proteins to {args.out}")
