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

process combineFastq {
    // concatenate fastq and fastq.gz in a dir

    label "wfflu"
    cpus 1
    input:
        tuple path(directory), val(meta)
    output:
        tuple val(meta.sample_id), val(meta.type), val(meta.barcode), path("${meta.sample_id}.fastq.gz"), emit: fastqfiles
        path "${meta.sample_id}.stats", emit: fastqstats
    shell:
    """
    fastcat -s ${meta.sample_id} -q ${params.min_qscore} -r ${meta.sample_id}.stats -x ${directory} | seqkit seq -m 200 - > ${meta.sample_id}.fastq
    gzip ${meta.sample_id}.fastq
    """
}


process alignReads {
    // align reads to reference

    label "wfflu"
    cpus 1
    input:
        tuple val(sample_id), val(type), val(barcode), path(sample_fastq)
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

process coverageCalc {
  depth_threads = {params.threads >= 4  ? 4 : params.threads}
      label "wfflu"
      cpus depth_threads
      input:
          tuple val(sample_id), val(type), path(bam), path(bai)
      output:
          tuple val(sample_id), val(type), path("${sample_id}.depth.txt")
      """
      samtools depth -aa ${bam} -Q 20 -q 1 > ${sample_id}.depth.txt
      """

}

process downSample {
    label 'wfflu'
    cpus params.threads
    input:
        tuple val(sample_id), val(type), path(bam), path(bai)
        path reference
    output:
        tuple val(sample_id), val(type), path("${sample_id}_all_merged.sorted.bam"), path("${sample_id}_all_merged.sorted.bam.bai"), emit: alignments
    """
    # split bam
    header_count=`samtools view -H ${sample_id}.bam | wc -l`
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

      samtools view ${bam} \${region} | awk -v upper="\${upper}" -v lower="\${lower}" '{if(length(\$10) < upper && length(\$10) > lower) print \$0}' > ${sample_id}_\${region}.tmp.sam;

      # ignore regions with no reads
      count=`wc -l < ${sample_id}_\${region}.tmp.sam`

      if [ "\${count}" -eq "0" ];
      then
        echo "no reads in \${region} so continuing"
        cat head.sam | samtools view -bh > ${sample_id}_\${region}_all.bam
        continue;
      fi

      cat head.sam ${sample_id}_\${region}.tmp.sam | samtools view -bh > ${sample_id}_\${region}.bam

      samtools view -h -F16 ${sample_id}_\${region}.bam > ${sample_id}_\${region}_fwd.sam;
      head -\${lines} ${sample_id}_\${region}_fwd.sam | samtools view -bh - > ${sample_id}_\${region}_fwd.bam;

      samtools view -h -f16 ${sample_id}_\${region}.bam > ${sample_id}_\${region}_rev.sam;
      head -\${lines} ${sample_id}_\${region}_rev.sam | samtools view -bh - > ${sample_id}_\${region}_rev.bam;
      samtools merge ${sample_id}_\${region}_all.bam ${sample_id}_\${region}_fwd.bam ${sample_id}_\${region}_rev.bam;

    done < regions.txt

    samtools merge ${sample_id}_all_merged.bam *_all.bam
    samtools sort ${sample_id}_all_merged.bam > ${sample_id}_all_merged.sorted.bam
    samtools index ${sample_id}_all_merged.sorted.bam
    echo "done"

    """
}

process medakaVariants {
    label "wfflu"
    cpus params.threads
    input:
        tuple val(sample_id), val(type), path(bam), path(bai)
        path reference
    output:
        tuple val(sample_id), val(type), path("${sample_id}.annotate.filtered.vcf")
    """

    medaka consensus ${bam} ${sample_id}.hdf
    medaka variant --gvcf ${reference} ${sample_id}.hdf ${sample_id}.vcf --verbose
    medaka tools annotate --debug --pad 25 ${sample_id}.vcf ${reference} ${bam} ${sample_id}.annotate.vcf

    bcftools filter -e "ALT='.'" ${sample_id}.annotate.vcf | bcftools filter -o ${sample_id}.annotate.filtered.vcf -O v -e "INFO/DP<${params.min_coverage}" -
    """

}

process makeConsensus {
    label "wfflu"
    cpus params.threads
    input:
        tuple val(sample_id), val(type), path(vcf), path(depth)
        path reference
    output:
        tuple val(sample_id), val(type), path("${sample_id}.draft.consensus.fasta")
    """
    awk '{if (\$3<${params.min_coverage}) print \$1"\t"\$2+1}' ${depth} > mask.regions
    bgzip ${vcf}
    tabix ${vcf}.gz

    bcftools consensus --mask mask.regions  --mark-del '-' --mark-ins lc --fasta-ref ${reference} -o ${sample_id}.draft.consensus.fasta ${vcf}.gz
    """
}

process typeFlu {
    label "wfflutyping"
    cpus params.threads
    input:
        tuple val(sample_id), val(type), path(consensus)
        path(blastdb)
    output:
        tuple val(sample_id), val(type), path("${sample_id}.insaflu.typing.txt"), emit: typing
        path "abricate.version", emit: version
    """
    abricate --version | sed 's/ /,/' > abricate.version
    abricate --datadir ${blastdb} --db insaflu -minid 70 -mincov 60 --quiet ${consensus} 1> ${sample_id}.insaflu.typing.txt
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
        path "fastqstats/*"
        path "params.json"
    output:
        tuple path("wf-flu-*.html"), path("wf-flu-results.csv")
    script:
        report_name = "wf-flu-report.html"
        def metadata = new JsonBuilder(metadata).toPrettyString()
    """
    echo '${metadata}' > metadata.json
    report.py $report_name \
        --versions versions \
        --coverage coverage \
        --typing typing \
        --fastqstats fastqstats \
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
        reads
        reference
        blastdb
    main:
        fastq = combineFastq(reads)
        alignment = alignReads(fastq.fastqfiles,reference)
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

        // sample_ids = fastq.fastqfiles.map{it -> it[0]}.collect().map{it -> it.join(' ')}
        // sample_types = fastq.fastqfiles.map{it-> it[1]}.collect().map{it -> it.join(' ')}
        // sample_barcodes = fastq.fastqfiles.map{it-> it[2]}.collect().map{it -> it.join(' ')}


        report = makeReport(
            reads.map { it -> return it[1] }.toList(),
            software_versions.collect(),
            coverage.map{it -> it[2] }.collect(),
            type.typing.map{it -> it[2] }.collect(),
            fastq.fastqstats.collect(),
            workflow_params
        )

    emit:
        results = fastq.fastqstats.concat(
            report.collect(),
            output_alignments.collect(),
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
       "input":params.fastq,
       "sample":params.sample,
       "sample_sheet":params.sample_sheet])

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
