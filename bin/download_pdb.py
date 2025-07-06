import argparse
import json
import os
import urllib.request

parser = argparse.ArgumentParser(description="Download AlphaFold PDBs for UniProt IDs in a family")
parser.add_argument('--family_json', required=True, help='Path to family JSON from get_family.py')
parser.add_argument('--out', required=True, help='Output .txt file with list of downloaded PDBs')
args = parser.parse_args()

out_dir = "pdb_files"
os.makedirs(out_dir, exist_ok=True)

with open(args.family_json) as f:
    family = json.load(f)

ids = family.get("family_members") or family.get("uniprot_ids")
if not ids:
    raise ValueError("No family members found in JSON.")

downloaded = []
for acc in ids:
    pdb_url = f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v4.pdb"
    out_file = os.path.join(out_dir, f"{acc}.pdb")
    try:
        if not os.path.exists(out_file):
            print(f"Downloading {acc} ...")
            urllib.request.urlretrieve(pdb_url, out_file)
        downloaded.append(os.path.abspath(out_file))
    except Exception as e:
        print(f"WARNING: Could not download {acc}: {e}")

# Output a text file with all successfully downloaded PDB file paths
with open(args.out, "w") as out:
    for path in downloaded:
        out.write(path + "\n")

print(f"Downloaded {len(downloaded)} AlphaFold PDBs to '{out_dir}' and listed them in '{args.out}'")
