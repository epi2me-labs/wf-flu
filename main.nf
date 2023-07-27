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

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

process alignReads {
    label "wfflu"
    cpus 2
    input:
        tuple val(meta), path(sample_fastq), path(fastq_stats)
        path reference
    output:
        tuple val(meta), path("${meta.alias}.bam"), path("${meta.alias}.bam.bai"), emit: alignments
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
          tuple val(meta), path(bam), path(bai)
      output:
          tuple val(meta.alias), val(meta), path("${meta.alias}.depth.txt")
      """
      samtools depth -aa ${bam} -Q 20 -q 1 > ${meta.alias}.depth.txt
      """

}

process downSample {
    label 'wfflu'
    cpus 2
    input:
        tuple val(meta), path(bam), path(bai)
        path reference
    output:
        tuple val(meta), path("${meta.alias}_all_merged.sorted.bam"), path("${meta.alias}_all_merged.sorted.bam.bai"), emit: alignments
    """
    # split bam
    header_count=`samtools view -H ${meta.alias}.bam | wc -l`
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

      samtools view ${bam} \${region} | awk -v upper="\${upper}" -v lower="\${lower}" '{if(length(\$10) < upper && length(\$10) > lower) print \$0}' > ${meta.alias}_\${region}.tmp.sam;

      # ignore regions with no reads
      count=`wc -l < ${meta.alias}_\${region}.tmp.sam`

      if [ "\${count}" -eq "0" ];
      then
        echo "no reads in \${region} so continuing"
        cat head.sam | samtools view -bh > ${meta.alias}_\${region}_all.bam
        continue;
      fi

      cat head.sam ${meta.alias}_\${region}.tmp.sam | samtools view -bh > ${meta.alias}_\${region}.bam

      samtools view -h -F16 ${meta.alias}_\${region}.bam > ${meta.alias}_\${region}_fwd.sam;
      head -\${lines} ${meta.alias}_\${region}_fwd.sam | samtools view -bh - > ${meta.alias}_\${region}_fwd.bam;

      samtools view -h -f16 ${meta.alias}_\${region}.bam > ${meta.alias}_\${region}_rev.sam;
      head -\${lines} ${meta.alias}_\${region}_rev.sam | samtools view -bh - > ${meta.alias}_\${region}_rev.bam;
      samtools merge ${meta.alias}_\${region}_all.bam ${meta.alias}_\${region}_fwd.bam ${meta.alias}_\${region}_rev.bam;

    done < regions.txt

    samtools merge ${meta.alias}_all_merged.bam *_all.bam
    samtools sort ${meta.alias}_all_merged.bam > ${meta.alias}_all_merged.sorted.bam
    samtools index ${meta.alias}_all_merged.sorted.bam
    echo "done"

    """
}

process lookup_medaka_consensus_model {
    label "wfflu"
    input:
        path("lookup_table")
        val basecall_model
    output:
        stdout
    shell:
    '''
    medaka_model=$(workflow-glue resolve_medaka_model lookup_table '!{basecall_model}' "medaka_consensus")
    echo -n $medaka_model
    '''
}


process medakaVariants {
    label "medaka"
    cpus 2
    input:
        tuple val(meta), path(bam), path(bai), val(medaka_model)
        path reference
    output:
        tuple val(meta.alias), val(meta), path("${meta.alias}.annotate.filtered.vcf")
    script:
    """
    medaka consensus ${bam} ${meta.alias}.hdf --model ${medaka_model}
    medaka variant --gvcf ${reference} ${meta.alias}.hdf ${meta.alias}.vcf --verbose
    medaka tools annotate --debug --pad 25 ${meta.alias}.vcf ${reference} ${bam} ${meta.alias}.annotate.vcf

    bcftools filter -e "ALT='.'" ${meta.alias}.annotate.vcf | bcftools filter -o ${meta.alias}.annotate.filtered.vcf -O v -e "INFO/DP<${params.min_coverage}" -
    """
}


process makeConsensus {
    label "wfflu"
    cpus 2
    input:
        tuple val(meta), path(vcf), path(depth)
        path reference
    output:
        tuple val(meta), path("${meta.alias}.draft.consensus.fasta")
    """
    awk '{if (\$3<${params.min_coverage}) print \$1"\t"\$2+1}' ${depth} > mask.regions
    bgzip ${vcf}
    tabix ${vcf}.gz

    bcftools consensus --mask mask.regions  --mark-del '-' --mark-ins lc --fasta-ref ${reference} -o ${meta.alias}.draft.consensus.fasta ${vcf}.gz
    """
}

process typeFlu {
    label "wfflutyping"
    cpus 2
    input:
        tuple val(meta), path(consensus)
        path(blastdb)
    output:
        tuple val(meta), path("${meta.alias}.insaflu.typing.txt"), path(consensus), emit: typing
        path "abricate.version", emit: version
    """
    abricate --version | sed 's/ /,/' > abricate.version
    abricate --datadir ${blastdb} --db insaflu -minid 70 -mincov 60 --quiet ${consensus} 1> ${meta.alias}.insaflu.typing.txt
    """
}

process processType {
    label "wfflu"
    cpus 1
    input:
        tuple val(meta), path(typing), path(consensus)
    output:
        tuple val(meta), path("${meta.alias}.typing.json"), path(consensus)
    script:
    """
    workflow-glue process_abricate --typing ${typing} --output ${meta.alias}.typing.json
    """
}


process prepNextclade {
    label "wfflu"
    input:
        tuple val(meta), path(typing_json), path(consensus)
        path nextclade_data
    output:
        path "datasets", optional: true
    script:
    """
    mkdir datasets
    workflow-glue nextclade_helper \
        --typing ${typing_json} \
        --nextclade_datasets ${nextclade_data} \
        --consensus ${consensus} \
        --sample_alias ${meta.alias}
    
    if [ -z "\$(ls -A datasets)" ]; then
        rm -r datasets
    fi
    """
}


process nextclade {
    label "nextclade"
    input:
        tuple val(dataset), path(files)
    output:
        tuple val(dataset), path("${dataset}/${dataset}.json")
    script:
    """
    mkdir ${dataset}
    cat *.fasta > ${dataset}/${dataset}.consensus.fasta
    if [[ "${dataset}" != "flu_h1n1pdm_na" ]]; then
        nextclade dataset get --name \"${dataset}\" --output-dir "nextclade_datasets/${dataset}"
    else
        nextclade dataset get --name \"${dataset}\" --output-dir "nextclade_datasets/${dataset}" --reference MW626056
    fi
    nextclade run --input-dataset nextclade_datasets/${dataset} --output-all=${dataset}/ ${dataset}/${dataset}.consensus.fasta
    mv ${dataset}/nextclade.json ${dataset}/${dataset}.json
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
        path "processed_type/*"
        path "fastqstats/per-file-stats.tsv"
        path "nextclade/*"
        path nextclade_datasets
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
        --processed_type processed_type \
        --fastqstats fastqstats/per-file-stats.tsv \
        --nextclade_files nextclade/* \
        --nextclade_datasets ${nextclade_datasets} \
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
        nextclade_data
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

        if(params.medaka_consensus_model) {
            log.warn "Overriding Medaka consensus model with ${params.medaka_consensus_model}."
            medaka_consensus_model = Channel.fromList([params.medaka_consensus_model])
        }
        else {
            lookup_table = Channel.fromPath("${projectDir}/data/medaka_models.tsv", checkIfExists: true)
            medaka_consensus_model = lookup_medaka_consensus_model(lookup_table, params.basecaller_cfg)
        }

        bams_for_calling = downsample.alignments.combine(medaka_consensus_model)

        variants = medakaVariants(bams_for_calling, reference)

        for_draft = variants.join(coverage.map{it -> return tuple(it[0], it[2])})
        
        draft = makeConsensus(for_draft.map {it -> tuple(it[1..3])}, reference)
        flu_type = typeFlu(draft, blastdb)

        processed_type = processType(flu_type.typing)

        nextclade_prep = prepNextclade(processed_type, nextclade_data)
        
        nextclade_datasets = nextclade_prep
        | map { file(it.resolve("**"), type: "file") }
        | flatten
        | map { [it.parent.name, it] }
        | groupTuple

        nextclade_result = nextclade(nextclade_datasets)

        software_versions = getVersions()
        software_versions = software_versions.mix(flu_type.version.first())
        workflow_params = getParams()

        output_alignments = alignment.alignments.map{ it -> return tuple(it[1], it[2]) }

        report = makeReport(
            samples.map{it -> it[0]}.toList(),
            software_versions.collect(),
            coverage.map{it -> it[2] }.collect(),
            flu_type.typing.map{it -> it[1] }.collect(),
            processed_type.map{it -> it[1] }.collect(),
            samples | map { it[2].resolve("per-read-stats.tsv") } | collectFile(keepHeader: true),
            nextclade_result.map{it -> it[1]}.collect(),
            nextclade_data,
            workflow_params
        )

    emit:
        results = output_alignments.concat(
            report.collect(),
            variants.map{it-> it[2]},
            draft.map{it -> it[1]},
            flu_type.typing.map{it -> it[1]},
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

    nextclade_data = projectDir.resolve("./data/nextclade.csv").toString()

    pipeline(samples, params._reference, params._blastdb, nextclade_data)
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
