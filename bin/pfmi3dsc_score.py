import argparse
import json
import pandas as pd
import numpy as np
import os

def pfmi3dsc_result_matrix(uniprot_input_protein, aln_tsv, biomuta_profile, hotspot_profile, primary_list):


    s = pd.read_csv(aln_tsv, sep='\t')
    s.columns = ['query', 'target', 'qaln', 'taln', 'gapopen', 'qstart', 'qend', 'tstart', 'tend', 'qlen', 'tlen',
                 'tcov', 'qseq', 'tseq']

    query_seq_list = []
    target_seq_list = []
    protein_names = []
    protein_lengths = []

    for i in s.itertuples():
        q = list(i.qaln)
        t = list(i.taln)
        query_seq_list.append(q)
        target_seq_list.append(t)
        protein_names.append(i.target.split('.')[0])
        protein_lengths.append(i.tlen)

    alignment_matrix_seq = pd.DataFrame(target_seq_list).transpose()
    alignment_matrix_seq.columns = protein_names
    alignment_matrix_seq['pos'] = np.arange(1, len(alignment_matrix_seq)+1)
    alignment_matrix_seq = alignment_matrix_seq.set_index('pos')

    # === Score Matrix Construction ===
    score_matrix = pd.DataFrame(0, index=alignment_matrix_seq.index, columns=protein_names)
    for pname in protein_names:
        bm = set(biomuta_profile.get(pname, []))
        hs = set(hotspot_profile.get(pname, []))
        for pos in score_matrix.index:
            if pos in hs:
                score_matrix.loc[pos, pname] = 2
            elif pos in bm:
                score_matrix.loc[pos, pname] = 1

    # === Residue Score (sum over family) ===
    scores_of_residues = score_matrix.sum(axis=1)
    alignment_matrix_seq['scores'] = scores_of_residues

    # === Probability Calculation ===
    mutation_probability = {}
    hotspot_probability = {}
    other_probability = {}

    for pname, plen in zip(protein_names, protein_lengths):
        m = len([p for p in biomuta_profile.get(pname, []) if p not in hotspot_profile.get(pname, [])])
        h = len(hotspot_profile.get(pname, []))
        mutation_probability[pname] = m / plen
        hotspot_probability[pname] = h / plen
        other_probability[pname] = max(0.0, 1.0 - (m + h) / plen)

    prob_matrix = pd.DataFrame(index=alignment_matrix_seq.index, columns=protein_names)

    for pos in alignment_matrix_seq.index:
        for pname in protein_names:
            val = score_matrix.loc[pos, pname]
            if pd.isna(val):
                prob_matrix.loc[pos, pname] = 1.0
            elif val == 1:
                prob_matrix.loc[pos, pname] = mutation_probability[pname]
            elif val == 2:
                prob_matrix.loc[pos, pname] = hotspot_probability[pname]
            else:
                prob_matrix.loc[pos, pname] = other_probability[pname]

    # Remove query column for product
    columns = [c for c in protein_names if c != uniprot_input_protein]
    prob_matrix_partial = prob_matrix[columns].astype(float)
    probability_of_residues = prob_matrix_partial.prod(axis=1, skipna=True)
    alignment_matrix_seq['probability'] = probability_of_residues

    # === Final Prediction ===
    threshold = 0.01 / (len(primary_list[2]) - 1)
    threshold1 = len(primary_list[2]) / 2
    alignment_matrix_seq['Predicted Functional'] = (
        (alignment_matrix_seq['probability'] < threshold) &
        (alignment_matrix_seq['scores'] > threshold1)
    ).map({True: "★", False: ""})

    # Extract result
    result = alignment_matrix_seq.loc[alignment_matrix_seq['Predicted Functional'] == "★"].reset_index()
    result_residues = result['pos'].tolist()

    return alignment_matrix_seq, result, result_residues, score_matrix, prob_matrix


def main():
    parser = argparse.ArgumentParser(description="Score PFMI3DSC functional residues")
    parser.add_argument('--align_json', required=True, help='Alignment .tsv file from Foldseek')
    parser.add_argument('--mut_json', required=True, help='Mutations/hotspot JSON from parse_mutations')
    parser.add_argument('--out', required=True, help='Output JSON filename')
    args = parser.parse_args()

    # Load mutation/hotspot data
    with open(args.mut_json) as f:
        mutdat = json.load(f)
    biomuta = mutdat['biomuta_profile']
    hotspot = mutdat['hotspot_profile']
    uniprot_ids = mutdat['uniprot_ids']
    gene_names = mutdat['gene_names']

    # Determine the query/input protein (should match user's input, e.g. P01112)
    # Here, let's take the first uniprot_id as the query
    uniprot_input_protein = uniprot_ids[0]
    primary_list = [uniprot_input_protein, '', uniprot_ids, gene_names, uniprot_input_protein, '']

    # Run main scoring function (ports your result_matrix logic)
    alignment_matrix_seq, result, result_residues, score_matrix, prob_matrix = pfmi3dsc_result_matrix(
        uniprot_input_protein,
        args.align_json,
        biomuta,
        hotspot,
        primary_list
    )

    # Output as JSON (and optionally as CSV for easy inspection)
    outdata = {
        "uniprot_input_protein": uniprot_input_protein,
        "alignment_matrix_seq": alignment_matrix_seq.to_dict(),
        "result_residues": result_residues,
        "result_matrix": result.to_dict(),
        "score_matrix": score_matrix.to_dict(),
        "probability_matrix": prob_matrix.to_dict(),
    }
    with open(args.out, "w") as f:
        json.dump(outdata, f, indent=2)
    print(f"PFMI3DSC scoring done for {uniprot_input_protein}. Output: {args.out}")

if __name__ == "__main__":
    main()
