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
include { start_ping; end_ping } from './lib/ping'


process summariseReads {
    // concatenate fastq and fastq.gz in a dir

    label "wfemu"
    cpus 1
    input:
        tuple path(directory), val(sample_id), val(type)
    output:
        path "${sample_id}.stats"
    shell:
    """
    fastcat -s ${sample_id} -r ${sample_id}.stats -x ${directory} > /dev/null
    """
}


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
    storeDir "store_dir/"
    output:
        path "database_dir/${params.database}"
    """
    mkdir database_dir
    mkdir database_dir/${params.database}
    wget -qO- https://gitlab.com/treangenlab/emu/-/archive/v3.4.2/emu-v3.4.2.tar.gz | tar -C database_dir/${params.database} -xvz --strip-components=2 emu-v3.4.2/${params.database}_database/
    """
}

process mapReads {
    label "wfemu"
    cpus params.threads
    publishDir "${params.out_dir}", mode: 'copy', overwrite: true
    input:
        tuple path(directory), val(sample_id), val(type)
        path database
    output:
        path "minimap-alignments/*"
    """
    mkdir minimap-alignments
    for file in ${directory}/*
    do
        minimap2 -ax ${params.type} -t ${params.threads} -N ${params.num_alignments} -p 0.9 -K ${params.K} "${database}/species_taxid.fasta" \${file} -o "minimap-alignments/${sample_id}.sam"
    done
    """
}

process runEmu {
    label "wfemu"
    cpus params.threads
    publishDir "${params.out_dir}", mode: 'move', overwrite: true
    input:
        path minimap_alignments
        path database
    output:
        path "${params.emu_output_dir}/*"
    """
    for file in ${minimap_alignments}
    do
        echo \$file
        emu abundance \$file --threads ${task.cpus} --output-dir ${params.emu_output_dir} --db ${database} --min-abundance ${params.min_abundance} --keep-counts --keep-read-assignments --output-unclassified
    done
    """
}



process makeReport {
    label "wfemu"
    input:
        path "seqs.txt"
        path "versions/*"
        path "params.json"
    output:
        path "wf-emu-*.html"
    script:
        report_name = "wf-emu-" + params.report_name + '.html'
    """
    report.py $report_name \
        --versions versions \
        seqs.txt \
        --params params.json
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
        reads
    main:
        summary = summariseReads(reads)
        software_versions = getVersions()
        workflow_params = getParams()
        report = makeReport(summary, software_versions.collect(), workflow_params)

        database = unpackDatabase()
        minimap_alignments = mapReads(reads, database)
        emu = runEmu(minimap_alignments, database)
        // TODO: add Emu pictoral Emu results to report


    emit:
        results = summary.concat(report)
        // TODO: use something more useful as telemetry
        telemetry = workflow_params
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    start_ping()
    samples = fastq_ingress(
        params.fastq, params.out_dir, params.sample, params.sample_sheet, params.sanitize_fastq)

    //call emu pipeline
    pipeline(samples)
    output(pipeline.out.results)
    end_ping(pipeline.out.telemetry)




}

