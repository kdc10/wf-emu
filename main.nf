#!/usr/bin/env nextflow

// Developer notes
// 
// This template workflow provides a basic structure to copy in order
// to create a new workflow. Current recommended practices are:
//     i) create a simple command-line interface.
//    ii) include an abstract workflow scope named "pipeline" to be used
//        in a module fashion
//   iii) a second concreate, but anonymous, workflow scope to be used
//        as an entry point when using this workflow in isolation.

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/fastqingress'
//include { start_ping; end_ping } from './lib/ping'

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

process getVersions {
    label "wfemu"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    python --version | sed 's/^/python,/' >> versions.txt
    emu --version | sed 's/^/emu,/' >> versions.txt
    minimap2 --version | sed 's/^/minimap2,/' >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    python -c "import numpy; print(f'numpy,{numpy.__version__}')" >> versions.txt
    python -c "import pandas; print(f'pandas,{pandas.__version__}')" >> versions.txt
    """
}


process getParams {
    label "wfemu"
    cpus 1
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}

process unpackDatabase {
    label "wfemu"
    cpus 1
    storeDir "${params.store_dir}/${database.simpleName}"
    //storeDir "${params.store_dir}/"
    input:
        //path fasta
        //path taxonomy
        path database
    output:
        path "database_dir/"
        //wget "${database}" ADD THIS from the intrenal DB
    """
    mkdir database_dir
    tar xvzf "${database}" -C database_dir/
    """
    //mkdir database_dir
    //cp "${taxonomy}" database_dir/
    //cp "${fasta}" database_dir/
    //"""
}

process mapReads {
    label "wfemu"
    cpus params.threads
    input:
        tuple val(meta), path(concat_seqs), path(fastcat_stats)
        path database
    output:
        tuple(
            val(meta), path("*.sam"), emit: sam)
    script:
        def sample_id = meta["alias"]
    """
    minimap2 -ax map-ont "${database}/species_taxid.fasta" "${concat_seqs}" -t "${task.cpus}" -N ${params.num_alignments} -p 0.9 -K ${params.K} -o "${sample_id}.sam"
    """
}

process runEmu {
    label "wfemu"
    cpus params.threads
    publishDir "${params.out_dir}", mode: 'copy', overwrite: true, pattern: "*", saveAs: {name -> "emu_results"}
    input:
        tuple val(meta), path(sam)
        path database
    output:
        path "emu_results/*_read-assignment-distributions.tsv", emit: read_assignment_distributions_tsv
        path "emu_results/*_rel-abundance.tsv", emit: rel_abundance_tsv
        path "emu_results/*_unclassified.fa", emit: unclassified_fa
    script:
        def sample_id = meta["alias"]
    """
    emu abundance ${sam} --output-dir emu_results --db ${database} --min-abundance ${params.min_abundance} --output-basename ${sample_id} \
    --keep-counts --keep-read-assignments --output-unclassified
    """
}


process makeReport {
    label "wfemu"
    input:
        val metadata
        path per_read_stats
        path rel_abundance_tsv
        path "versions/*"
        path "params.json"
    output:
        path "wf-emu-*.html"
    script:
        String report_name = "wf-emu-report.html"
        String metadata = new JsonBuilder(metadata).toPrettyString()
        String stats_args = \
            (per_read_stats.name == OPTIONAL_FILE.name) ? "" : "--stats $per_read_stats"
    """
    echo '${metadata}' > metadata.json
    workflow-glue report $report_name \
        --versions versions \
        $stats_args \
        --rel_abun $rel_abundance_tsv\
        --params params.json \
        --metadata metadata.json
    """
}



// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    label "wfemu"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        path fname
    output:
        path fname
    """
    echo "Writing output files."
    """
}


// workflow module
workflow pipeline {
    take:
        samples
        database
    main:
        software_versions = getVersions()
        workflow_params = getParams()
        database = unpackDatabase(database)
        // Initial reads QC
        per_read_stats = samples.map {
            it[2] ? it[2].resolve('per-read-stats.tsv') : null
        }
        | collectFile ( keepHeader: true )
        | ifEmpty ( OPTIONAL_FILE )
        metadata = samples.map { it[0] }.toList()


        // Run Minimap2

        minimap_alignments = mapReads(
                samples
                | map { [it[0], it[1], it[2] ?: OPTIONAL_FILE ] },
                database
            )

        emu = runEmu(
            minimap_alignments
            | map{[it[0], it[1]]},
            database
            )

        // output updating files as part of this pipeline
        // TODO: add Emu pictoral Emu results to report
        report = makeReport(
            metadata, per_read_stats, emu.rel_abundance_tsv.collect(), software_versions.collect(), workflow_params
        )

        output(report.mix(
            software_versions, workflow_params, emu.rel_abundance_tsv, emu.unclassified_fa, emu.read_assignment_distributions_tsv),
        )


    emit:
        report
        workflow_params

}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {


    // Check source param is valid
    sources = params.database_sets
    source_name = params.database_set
    source_data = sources.get(source_name, false)
    if (!sources.containsKey(source_name) || !source_data) {
        keys = sources.keySet()
        throw new Exception("Source $params.source is invalid, must be one of $keys")
    }

    // Grab database files
    database_repo_path = sources[source_name]["database"]
    database = file("${projectDir}/${database_repo_path}", type: "file", checkIfExists:true)
    samples = fastq_ingress([
        "input":params.fastq,
        "sample":params.sample,
        "sample_sheet":params.sample_sheet,
        "analyse_unclassified":false,
        "fastcat_stats": params.wf.fastcat_stats,
        "fastcat_extra_args": ""])

    
    results = pipeline(samples, database)
}

