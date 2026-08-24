#!/usr/bin/env nextflow
/*
 * Reproduces the full community-scale metabolic modeling pipeline described
 * in README.md: community model (42 real biological replicates, 4-taxon
 * community) -> two-way ANOVA -> host clock coupling, plus an independent
 * thermodynamic feasibility check.
 *
 * This wraps the existing, already-validated scripts in src/ as-is; it does
 * not re-implement their logic, only sequences them and collects their real
 * outputs into results/ and figures/.
 *
 * Usage:
 *   conda env create -f environment.yml
 *   nextflow run main.nf
 */
nextflow.enable.dsl = 2

params.repo_dir = "$projectDir"

process COMMUNITY_MODEL {
    tag "4-taxon community, 42 real replicates"
    publishDir "${params.repo_dir}/results", mode: 'copy'

    output:
    path "butyrate_all_replicates_results.json", emit: results_json

    script:
    """
    python ${params.repo_dir}/src/community_modeling/butyrate.py
    """
}

process ANOVA {
    tag "two-way ANOVA: condition x ZT"
    publishDir "${params.repo_dir}/results", mode: 'copy'

    input:
    path results_json

    output:
    path "anova_summary.txt", emit: summary

    script:
    """
    python ${params.repo_dir}/src/community_modeling/butyrate_anova.py | tee anova_summary.txt
    """
}

process HOST_CLOCK {
    tag "Goodwin3 diurnal-forced host clock"
    publishDir "${params.repo_dir}/figures", mode: 'copy', pattern: "*.png"

    input:
    path results_json

    output:
    path "*.png", emit: figures

    script:
    """
    python ${params.repo_dir}/src/host_clock/goodwin3_diurnal.py
    """
}

process THERMODYNAMICS {
    tag "tFBA: rate-limiting reaction feasibility"
    publishDir "${params.repo_dir}/results", mode: 'copy'

    output:
    path "thermo_summary.txt", emit: summary

    script:
    """
    python ${params.repo_dir}/src/thermodynamics/ecfba_thermo_check.py | tee thermo_summary.txt
    """
}

workflow {
    community = COMMUNITY_MODEL()
    ANOVA(community.results_json)
    HOST_CLOCK(community.results_json)
    THERMODYNAMICS()
}
