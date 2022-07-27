#!/usr/bin/env nextflow

// Developer notes
//
// This template workflow provides a basic structure to copy in order
// to create a new workflow. Current recommended pratices are:
//     i) create a simple command-line interface.
//    ii) include an abstract workflow scope named "pipeline" to be used
//        in a module fashion
//   iii) a second concreate, but anonymous, workflow scope to be used
//        as an entry point when using this workflow in isolation.

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/fastqingress'
include { start_ping; end_ping } from './lib/ping'


process combineFastq {
    // concatenate fastq and fastq.gz in a dir

    label "wfflu"
    cpus 1
    input:
        tuple path(directory), val(meta)
    output:
        tuple val(meta.sample_id), val(meta.type), path("${meta.sample_id}.fastq.gz"), emit: fastqfiles
        path "${meta.sample_id}.stats", emit: fastqstats
    shell:
    """
    fastcat -s ${meta.sample_id} -r ${meta.sample_id}.stats -x ${directory} | seqkit seq -m 200 - > ${meta.sample_id}.fastq
    gzip ${meta.sample_id}.fastq
    """
}


process alignReads {
    // align reads to reference

    label "wfflu"
    cpus 1
    input:
        tuple val(sample_id), val(type), path(sample_fastq)
        path reference
    output:
        tuple val(sample_id), val(type), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: alignments
        tuple path("${sample_id}.bamstats"), path("${sample_id}.bam.summary"), emit: bamstats
    shell:
    """
    mini_align -i ${sample_fastq} -r ${reference} -p ${sample_id}_tmp -t $task.cpus -m

    # keep only mapped reads
    samtools view --write-index -F 4 ${sample_id}_tmp.bam -o ${sample_id}.bam##idx##${sample_id}.bam.bai

    # get stats from bam
    stats_from_bam -o ${sample_id}.bamstats -s ${sample_id}.bam.summary -t $task.cpus ${sample_id}.bam
    """
}


process getVersions {
    label "wfflu"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    """
}


process getParams {
    label "wfflu"
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


process makeReport {
    label "wfflu"
    input:
        path "seqs.txt"
        path "versions/*"
        path "params.json"
    output:
        path "wf-template-*.html"
    script:
        report_name = "wf-template-" + params.report_name + '.html'
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
    label "wfflu"
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
        reference
    main:
        fastq = combineFastq(reads)

        alignment = alignReads(fastq.fastqfiles,reference)

        software_versions = getVersions()
        workflow_params = getParams()

        output_alignments = alignment.alignments.map{ it -> return tuple(it[2], it[3]) }

        report = makeReport(fastq.fastqstats, software_versions.collect(), workflow_params)

    emit:
        results = fastq.fastqstats.concat(
            report,
            output_alignments.collect()
        )
        // TODO: use something more useful as telemetry
        telemetry = workflow_params
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    start_ping()
    samples = fastq_ingress([
        "input":params.fastq,
        "sample":params.sample,
        "sample_sheet":params.sample_sheet,
        "sanitize": params.sanitize_fastq,
        "output":params.out_dir])

  //get reference
    if (params.reference == null){
      params.remove('reference')
      params._reference = projectDir.resolve("./data/primer_schemes/V1/consensus_irma.fasta").toString()
    } else {
      params._reference = file(params.reference, type: "file", checkIfExists:true).toString()
      params.remove('reference')
    }

    pipeline(samples,params._reference)
    output(pipeline.out.results)
    end_ping(pipeline.out.telemetry)
}
