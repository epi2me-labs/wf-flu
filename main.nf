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

process alignReads {
    label "wfflu"
    cpus 2
    input:
        tuple val(meta), path(sample_fastq), path(fastq_stats)
        path reference
    output:
        tuple val(meta.alias), val(meta.type), path("${meta.alias}.bam"), path("${meta.alias}.bam.bai"), emit: alignments
        tuple path("${meta.alias}.bamstats"), path("${meta.alias}.bam.summary"), emit: bamstats
    shell:
    """
    mini_align -i ${sample_fastq} -r ${reference} -p ${meta.alias}_tmp -t ${params.align_threads} -m

    # keep only mapped reads
    samtools view --write-index -F 4 ${meta.alias}_tmp.bam -o ${meta.alias}.bam##idx##${meta.alias}.bam.bai

    # get stats from bam
    stats_from_bam -o ${meta.alias}.bamstats -s ${meta.alias}.bam.summary -t ${params.align_threads} ${meta.alias}.bam
    """
}

process coverageCalc {
      label "wfflu"
      cpus 2
      input:
          tuple val(alias), val(type), path(bam), path(bai)
      output:
          tuple val(alias), val(type), path("${alias}.depth.txt")
      """
      samtools depth -aa ${bam} -Q 20 -q 1 > ${alias}.depth.txt
      """

}

process downSample {
    label 'wfflu'
    cpus 2
    input:
        tuple val(alias), val(type), path(bam), path(bai)
        path reference
    output:
        tuple val(alias), val(type), path("${alias}_all_merged.sorted.bam"), path("${alias}_all_merged.sorted.bam.bai"), emit: alignments
    """
    # split bam
    header_count=`samtools view -H ${alias}.bam | wc -l`
    lines=\$(( ${params.downsample} + \$header_count + 1 ))

    # get the regions from the fasta
    awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length(\$0)}END{print seqlen}' ${reference} | tr "\\n" "," | tr '>' '\\n' | sed 's/,\$//g' | grep ',' > regions.txt

    # for every region we're going to downsample separatley

    while read region_string;
    do

      region=`echo \${region_string} | cut -f1 -d','`
      length=`echo \${region_string} | cut -f2 -d','`

      # get upper and lower bounds of reference span

      upper=`echo \$((\${length}+(\${length}*10/100)))`
      lower=`echo \$((\${length}-(\${length}*10/100)))`

      samtools view -H ${bam} \${region} > head.sam

      # filter reads in region and covering region

      samtools view ${bam} \${region} | awk -v upper="\${upper}" -v lower="\${lower}" '{if(length(\$10) < upper && length(\$10) > lower) print \$0}' > ${alias}_\${region}.tmp.sam;

      # ignore regions with no reads
      count=`wc -l < ${alias}_\${region}.tmp.sam`

      if [ "\${count}" -eq "0" ];
      then
        echo "no reads in \${region} so continuing"
        cat head.sam | samtools view -bh > ${alias}_\${region}_all.bam
        continue;
      fi

      cat head.sam ${alias}_\${region}.tmp.sam | samtools view -bh > ${alias}_\${region}.bam

      samtools view -h -F16 ${alias}_\${region}.bam > ${alias}_\${region}_fwd.sam;
      head -\${lines} ${alias}_\${region}_fwd.sam | samtools view -bh - > ${alias}_\${region}_fwd.bam;

      samtools view -h -f16 ${alias}_\${region}.bam > ${alias}_\${region}_rev.sam;
      head -\${lines} ${alias}_\${region}_rev.sam | samtools view -bh - > ${alias}_\${region}_rev.bam;
      samtools merge ${alias}_\${region}_all.bam ${alias}_\${region}_fwd.bam ${alias}_\${region}_rev.bam;

    done < regions.txt

    samtools merge ${alias}_all_merged.bam *_all.bam
    samtools sort ${alias}_all_merged.bam > ${alias}_all_merged.sorted.bam
    samtools index ${alias}_all_merged.sorted.bam
    echo "done"

    """
}

process medakaVariants {
    label "wfflu"
    cpus 2
    input:
        tuple val(alias), val(type), path(bam), path(bai)
        path reference
    output:
        tuple val(alias), val(type), path("${alias}.annotate.filtered.vcf")
    """

    medaka consensus ${bam} ${alias}.hdf
    medaka variant --gvcf ${reference} ${alias}.hdf ${alias}.vcf --verbose
    medaka tools annotate --debug --pad 25 ${alias}.vcf ${reference} ${bam} ${alias}.annotate.vcf

    bcftools filter -e "ALT='.'" ${alias}.annotate.vcf | bcftools filter -o ${alias}.annotate.filtered.vcf -O v -e "INFO/DP<${params.min_coverage}" -
    """

}

process makeConsensus {
    label "wfflu"
    cpus 2
    input:
        tuple val(alias), val(type), path(vcf), path(depth)
        path reference
    output:
        tuple val(alias), val(type), path("${alias}.draft.consensus.fasta")
    """
    awk '{if (\$3<${params.min_coverage}) print \$1"\t"\$2+1}' ${depth} > mask.regions
    bgzip ${vcf}
    tabix ${vcf}.gz

    bcftools consensus --mask mask.regions  --mark-del '-' --mark-ins lc --fasta-ref ${reference} -o ${alias}.draft.consensus.fasta ${vcf}.gz
    """
}

process typeFlu {
    label "wfflutyping"
    cpus 2
    input:
        tuple val(alias), val(type), path(consensus)
        path(blastdb)
    output:
        tuple val(alias), val(type), path("${alias}.insaflu.typing.txt"), emit: typing
        path "abricate.version", emit: version
    """
    abricate --version | sed 's/ /,/' > abricate.version
    abricate --datadir ${blastdb} --db insaflu -minid 70 -mincov 60 --quiet ${consensus} 1> ${alias}.insaflu.typing.txt
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
    medaka --version | sed 's/ /,/' >> versions.txt
    bcftools --version | head -1 | sed 's/ /,/' >> versions.txt
    samtools --version | grep samtools | head -1 | sed 's/ /,/' >> versions.txt
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
        val metadata
        path "versions/*"
        path "coverage/*"
        path "typing/*"
        path "fastqstats/per-file-stats.tsv"
        path "params.json"
    output:
        tuple path("wf-flu-*.html"), path("wf-flu-results.csv")
    script:
        report_name = "wf-flu-report.html"
        def metadata = new JsonBuilder(metadata).toPrettyString()
    """
    echo '${metadata}' > metadata.json
    workflow-glue report $report_name \
        --versions versions \
        --coverage coverage \
        --typing typing \
        --fastqstats fastqstats/per-file-stats.tsv \
        --params params.json \
        --metadata metadata.json
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
        samples
        reference
        blastdb
    main:

        alignment = alignReads(samples, reference)
        coverage = coverageCalc(alignment.alignments)

        // do crude downsampling
        if (params.downsample != null){
            println("Downsampling!!!")
            downsample = downSample(alignment.alignments, reference)
        } else {
            println("NOT Downsampling!!!")
            downsample = alignment
        }

        variants = medakaVariants(downsample.alignments, reference)

        for_draft = variants.join(coverage.map{it -> return tuple(it[0], it[2])})

        draft = makeConsensus(for_draft, reference)
        type = typeFlu(draft, blastdb)

        software_versions = getVersions()
        software_versions = software_versions.mix(type.version.first())
        workflow_params = getParams()

        output_alignments = alignment.alignments.map{ it -> return tuple(it[2], it[3]) }

        report = makeReport(
            samples.map { it -> return it[0] }.toList(),
            software_versions.collect(),
            coverage.map{it -> it[2] }.collect(),
            type.typing.map{it -> it[2] }.collect(),
            samples | map { it[2].resolve("per-read-stats.tsv") } | collectFile(keepHeader: true),
            workflow_params
        )

    emit:
        results = output_alignments.concat(
            report.collect(),
            variants.map{it-> it[2]},
            draft.map{it -> it[2]},
            type.typing.map{it -> it[2]},
            coverage.map{it -> it[2]}
        )
        // TODO: use something more useful as telemetry
        telemetry = workflow_params
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    if (params.disable_ping == false) {
        Pinguscript.ping_post(workflow, "start", "none", params.out_dir, params)
    }

    samples = fastq_ingress([
    "input": params.fastq,
    "fastcat_stats": true,
    "sample_sheet": params.sample_sheet])


  //get reference
    if (params.reference == null){
      params.remove('reference')
      params._reference = projectDir.resolve("./data/primer_schemes/V1/consensus_irma.fasta").toString()
    } else {
      params._reference = file(params.reference, type: "file", checkIfExists:true).toString()
      params.remove('reference')
    }

    //get reference
    if (params.blastdb == null){
      params.remove('blastdb')
      params._blastdb = projectDir.resolve("./data/primer_schemes/V1/blastdb").toString()
    } else {
      params._blastdb = file(params.reference, type: "file", checkIfExists:true).toString()
      params.remove('blastdb')
    }

    pipeline(samples,params._reference,params._blastdb)
    output(pipeline.out.results)
}

if (params.disable_ping == false) {
    workflow.onComplete {
        Pinguscript.ping_post(workflow, "end", "none", params.out_dir, params)
    }

    workflow.onError {
        Pinguscript.ping_post(workflow, "error", "$workflow.errorMessage", params.out_dir, params)
    }
}
