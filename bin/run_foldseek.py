import argparse
import os
import subprocess
import shutil

# ------------------ Argument Parsing ------------------
parser = argparse.ArgumentParser(description="Run Foldseek to align family PDBs")
parser.add_argument('--pdbs_txt', required=True, help='Text file listing all PDB files')
parser.add_argument('--query_accid', required=True, help='UniProt ACCID of query protein')
parser.add_argument('--out', required=True, help='Output alignment file (.tsv)')
args = parser.parse_args()

# ------------------ Load PDB Paths ------------------
with open(args.pdbs_txt) as f:
    pdbs = [os.path.abspath(line.strip()) for line in f if line.strip()]

# ------------------ Identify Query PDB ------------------
query_pdb = None
for pdb in pdbs:
    if args.query_accid in os.path.basename(pdb):
        query_pdb = pdb
        break

if not query_pdb:
    raise ValueError(f"Could not find PDB for query ACCID {args.query_accid}")

# ------------------ Remove Query from Family List ------------------
family_pdbs = [pdb for pdb in pdbs if pdb != query_pdb]
if not family_pdbs:
    raise ValueError("No family PDBs found (excluding the query) for pairwise alignment.")

# ------------------ Optional: Use database or direct pairwise ------------------
# Option A: Run pairwise alignment with each family PDB
# (This is the default behavior, uncomment Option B below if you want to use a Foldseek database)
# Prepare merged alignment output file
aln_tsv = args.out
if os.path.exists(aln_tsv):
    os.remove(aln_tsv)  # Avoid appending to an old file

os.makedirs("./tmpFoldseek", exist_ok=True)

for target_pdb in family_pdbs:
    pair_name = os.path.basename(target_pdb).split('.')[0]
    temp_out_file = f"./tmpFoldseek/{pair_name}.tsv"

    foldseek_cmd = [
        "foldseek", "easy-search",
        query_pdb,
        target_pdb,
        temp_out_file,
        "./tmpFoldseek",
        "--format-output",
        "query,target,qaln,taln,gapopen,qstart,qend,tstart,tend,qlen,tlen,tcov,qseq,tseq"
    ]

    print(f"[Foldseek] {os.path.basename(query_pdb)} vs {os.path.basename(target_pdb)}")
    try:
        subprocess.run(foldseek_cmd, check=True)
        print(f"[Output] Written to {temp_out_file}")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Foldseek failed for {target_pdb}: {e}")
        continue

    # Append to merged output
    with open(temp_out_file, "r") as fin, open(aln_tsv, "a") as fout:
        shutil.copyfileobj(fin, fout)
    print(f"[Merge] Appended {temp_out_file} → {aln_tsv}")


# ------------------ Option B: Use Foldseek database for faster search ------------------
# (Uncomment this section if you want DB-based search instead of pairwise)

# db_dir = "pdb_files_db"
# os.makedirs(db_dir, exist_ok=True)
# for pdb in pdbs:
#     dst = os.path.join(db_dir, os.path.basename(pdb))
#     if not os.path.exists(dst):
#         print(f"Copying {pdb} to {dst}")
#         shutil.copy2(pdb, dst)

# db_name = "family_db"
# if not (os.path.exists(f"{db_dir}/{db_name}.ffindex")):
#     subprocess.run(["foldseek", "createdb", db_dir, f"{db_dir}/{db_name}"], check=True)

# out_file = args.out
# foldseek_cmd = [
#     "foldseek", "easy-search",
#     query_pdb,
#     f"{db_dir}/{db_name}",
#     out_file,
#     "./tmpFoldseek",
#     "--format-output",
#     "query,target,qaln,taln,gapopen,qstart,qend,tstart,tend,qlen,tlen,tcov,qseq,tseq"
# ]
# print("Running Foldseek with query PDB:", query_pdb)
# subprocess.run(foldseek_cmd, check=True)
# print(f"Foldseek alignment done. Output written to {out_file}")
