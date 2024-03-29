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

include { fastq_ingress; xam_ingress } from './lib/ingress'

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

process getVersions {
    label "wf_common"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    python --version | sed 's/^/python,/' >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    python -c "import numpy; print(f'numpy,{numpy.__version__}')" >> versions.txt
    python -c "import pandas; print(f'pandas,{pandas.__version__}')" >> versions.txt
    """
}

process addEmuToVersions {
    label "wf_emu"
    cpus 1
    input:
        path "common_versions.txt"
    output:
        path "versions.txt"
    script:
    """
    emu --version | sed 's/^/emu,/' >> common_versions.txt
    minimap2 --version | sed 's/^/minimap2,/' >> common_versions.txt
    mv common_versions.txt versions.txt
    """
}


process getParams {
    label "wf_common"
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
    label "wf_emu"
    cpus 1
    storeDir "${params.store_dir}/${params.database_set}"
    input:
        val database
    output:
        path "database_dir/${params.database_set}"
    script:
        def db_name = "${database}".tokenize('/').last()
    """
    mkdir database_dir
    osf -p 56uf7 fetch "${database}"
    tar -xf "${db_name}" -C database_dir/
    """
}

process mapReads {
    label "wf_emu"
    cpus params.threads
    input:
        tuple val(meta), path(concat_seqs), path(stats)
        path database
    output:
        tuple(
            val(meta), path("*.sam"), path(stats), emit: sam)
    script:
        def sample_id = meta["alias"]
    """
    minimap2 -ax map-ont "${database}/species_taxid.fasta" "${concat_seqs}" -t "${task.cpus}" -N "${params.num_alignments}" -p 0.9 -K "${params.K}" -o "${sample_id}.sam"
    """
}

//Emu uses SAM files
process bam2sam {
    label "wf_common"
    cpus params.threads
    input:
        tuple val(meta), path(concat_seqs), path(stats)
    output:
        tuple(
            val(meta), path("*.sam"), path(stats), emit: sam)
    script:
        def sample_id = meta["alias"]
    """
    samtools view -h "${concat_seqs}" > "${sample_id}.sam"
    """
}

process runEmu {
    label "wf_emu"
    cpus params.threads
    input:
        tuple val(meta), path(sam), path(stats)
        path database
    output:
        tuple(
        val(meta),
        path("emu_results/*_read-assignment-distributions.tsv"),
        path("emu_results/*_rel-abundance.tsv"),
        path("emu_results/*_unclassified.fa"),
        emit:emu_results)

    script:
        def sample_id = meta["alias"]
        def output_unclassified = "${params.output_unclassified}"!= null ? "--output-unclassified" : ""
        def keep_read_assignments = "${params.keep_read_assignments}"!= null ? "--keep-read-assignments" : ""
        def keep_counts = "${params.keep_counts}"!= null ? "--keep-counts" : ""
    """
    emu abundance ${sam} --output-dir emu_results --db ${database} --min-abundance ${params.min_abundance} --output-basename ${sample_id} \
    ${keep_counts} ${keep_read_assignments} ${output_unclassified}
    """
}


process combineOutput {
    label "wf_emu"
    cpus params.threads
    input:
        path emu_results
    output:
        path "*.tsv" , emit: combine_output
    script:
        def split_tables = params.split_tables ? "--split-tables" : ""
        def counts = params.counts ? "--counts" : ""
    """
    emu combine-outputs . ${params.rank} ${counts} ${split_tables}
    """
}

process makeReport {
    label "wf_metagenomics"
    input:
        val metadata
        path "read_stats/per-read-stats*.tsv.gz"
        path rel_abundance_tsv
        path "versions/*"
        path "params.json"
    output:
        path "wf-emu-*.html", emit: report_html
    script:
        String report_name = "wf-emu-report.html"
        String metadata = new JsonBuilder(metadata).toPrettyString()
        def stats_args = params.wf.stats ? "--read_stats read_stats/*" : ""
    """
    echo '${metadata}' > metadata.json
    workflow-glue report $report_name \
        --versions versions \
        ${stats_args} \
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
    label "wf_common"
    publishDir (
        params.out_dir,
        mode: "copy",
        saveAs: { dirname ? "$dirname/$fname" : fname }
    )
    input:
        tuple path(fname), val(dirname)
    output:
        path fname
    """
    echo "Writing output files"
    """
}


// workflow module
workflow pipeline {
    Pinguscript.ping_start(nextflow, workflow, params)
    take:
        samples
        input_type
        database
    main:
        common_versions = getVersions()
        software_versions = addEmuToVersions(common_versions)
        workflow_params = getParams()

        database = unpackDatabase(database)
        // Initial reads QC
        // bamstats files have a different name
        per_read_stats = samples.map {it[2].resolve('*.tsv.gz')}.collect()
        metadata = samples.map { it[0] }.toList()

        samples.view()
        // Run Minimap2
        if (input_type == 'fastq') {
            minimap_alignments = mapReads(
                    samples
                    | map { [it[0], it[1], it[2] ?: OPTIONAL_FILE ] },
                    database
            )
        } else{
            minimap_alignments = bam2sam(
                samples
                | map { [it[0], it[1], it[2] ?: OPTIONAL_FILE ] },
            )
        }
    
        emu = runEmu(
            minimap_alignments
            | map{[it[0], it[1], it[2]]},
            database
            )
        
        emu_rel_abun = emu.emu_results.map { it[2] }.collect()
        combine_output = combineOutput(emu_rel_abun)

        report = makeReport(
            metadata,
            per_read_stats,
            emu_rel_abun,
            software_versions,
            workflow_params
        )


        ch_to_publish = Channel.empty()
        | mix(
            software_versions,
            workflow_params,
            report.report_html,
            combine_output.combine_output | flatten,
        )
        | map { [it, null] }

        ch_to_publish = ch_to_publish | mix (
            emu.emu_results | map { meta, read_assignment, rel_abundance, unclassified  -> [rel_abundance, "rel_abun"]},
            emu.emu_results | map { meta, read_assignment, rel_abundance, unclassified -> [unclassified, "unclassified"]},
            emu.emu_results | map { meta, read_assignment, rel_abundance, unclassified -> [read_assignment, "read_assignment"]},
        )
    
        ch_to_publish | output
    emit:
        report
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    ArrayList fastcat_extra_args = []
    if (params.min_len) { fastcat_extra_args << "-a $params.min_len" }
    if (params.max_len) { fastcat_extra_args << "-b $params.max_len" }
    if (params.min_read_qual) { fastcat_extra_args << "-q $params.min_read_qual" }

    // Check source param is valid
    sources = params.database_sets
    source_name = params.database_set
    source_data = sources.get(source_name, false)
    if (!sources.containsKey(source_name) || !source_data) {
        keys = sources.keySet()
        throw new Exception("Source $params.source is invalid, must be one of $keys")
    }

    // Grab database files
    database = sources[source_name]["database"]
    def input_type = ['fastq', 'bam'].findAll { params[it] }
    if (input_type.size() != 1) {
        error "Only provide one of '--fastq', '--bam'."
    }
    input_type = input_type[0]

    if (input_type == 'fastq') {
        samples = fastq_ingress([
            "input":params.fastq,
            "sample":params.sample,
            "sample_sheet":params.sample_sheet,
            "analyse_unclassified":params.analyse_unclassified,
            "stats": params.wf.stats,
            "required_sample_types": [],
            "fastcat_extra_args": fastcat_extra_args.join(" "),
            "watch_path": false])
    }
    results = pipeline(samples, input_type, database)
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}