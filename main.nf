nextflow.enable.dsl=2

process get_family {
    input:
        val accid
        file('bin/get_family.py')
    output:
        tuple val(accid), path("family_${accid}.json")
    script:
    """
    python3 bin/get_family.py --accid ${accid} --data /Users/hobzy/Desktop/PFMI3DSC-Nextflow/data/ --out family_${accid}.json
    """
}

process download_pdbs {
    input:
        tuple val(accid), path(family_json)
        file('bin/download_pdb.py')
    output:
        tuple val(accid), path("pdbs_${family_json.baseName}.txt")
    script:
    """
    python3 bin/download_pdb.py --family_json ${family_json} --out pdbs_${family_json.baseName}.txt
    """
}

process run_foldseek {
    input:
        tuple val(accid), path(pdbs_txt)
        file('bin/run_foldseek.py')
        file('foldseek/bin/foldseek')
    output:
        tuple val(accid), path("alignments_${pdbs_txt.baseName}.tsv")
    script:
    """
    chmod +x foldseek/bin/foldseek
    export PATH=\$PWD/foldseek/bin:\$PATH
    python3 bin/run_foldseek.py --pdbs_txt ${pdbs_txt} --query_accid ${accid} --out alignments_${pdbs_txt.baseName}.tsv
    """
}

process parse_mutations {
    input:
        tuple val(accid), path(family_json)
        file('bin/parse_mutations.py')
    output:
        tuple val(accid), path("mutations_${family_json.baseName}.json")
    script:
    """
    python3 bin/parse_mutations.py --family_json ${family_json} --data /Users/hobzy/Desktop/PFMI3DSC-Nextflow/data/ --out mutations_${family_json.baseName}.json
    """
}

process pfmi3dsc_score {
    input:
        tuple val(accid), path(alignment_file)
        tuple val(accid2), path(mut_json)
        file('bin/pfmi3dsc_score.py')
    output:
        tuple val(accid), path("results_${alignment_file.baseName}.json")
    when:
        accid == accid2
    script:
    """
    python3 bin/pfmi3dsc_score.py --align_json ${alignment_file} --mut_json ${mut_json} --out results_${alignment_file.baseName}.json
    """
}

process summarize_results {
    input:
        tuple val(accid), path(result_json)
        file('bin/summarize_results.py')
    output:
        path "results/summary_${result_json.baseName}.html"
    publishDir "results", mode: 'copy', pattern: "results/summary_*.html"
    script:
    """
    mkdir -p results
    python3 bin/summarize_results.py \
        --result_json ${result_json} \
        --out results/summary_${result_json.baseName}.html \
        --query_accid ${accid}
    """
}



/*
 * Utility: Clean all Nextflow temp and metadata files before rerun
 * Usage: nextflow clean or run in terminal:
 *   rm -rf work/ .nextflow*
 */

// To clean up before rerunning the workflow, run this command in your project root:
// rm -rf work/ .nextflow*


// -------- Workflow Definition --------

workflow {

    accid_channel = Channel.from(params.accids.split(','))

    family_json_ch = get_family(accid_channel, file('bin/get_family.py'))
    pdbs_ch = download_pdbs(family_json_ch, file('bin/download_pdb.py'))
    alignments_ch = run_foldseek(pdbs_ch, file('bin/run_foldseek.py'), file('foldseek/bin/foldseek'))
    mutdata_ch = parse_mutations(family_json_ch, file('bin/parse_mutations.py'))
    results_ch = pfmi3dsc_score(alignments_ch, mutdata_ch, file('bin/pfmi3dsc_score.py'))
    summarize_results(results_ch, file('bin/summarize_results.py'))
}

// -------- End of Workflow --------