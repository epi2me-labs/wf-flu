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

include { fastq_ingress } from './lib/ingress'

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

process alignReads {
    label "wfflu"
    cpus 2
    input:
        tuple val(meta), path("reads.fastq.gz")
        path "reference.fasta"
    output:
        tuple val(meta), path("align.bam"), path("align.bam.bai"), emit: alignments
        tuple path("align.bamstats"), path("align.bam.summary"), emit: bamstats
    shell:
    """
    mini_align -i reads.fastq.gz -r reference.fasta -p align_tmp -t ${params.align_threads} -m

    # keep only mapped reads
    samtools view --write-index -F 4 align_tmp.bam -o align.bam##idx##align.bam.bai

    # get stats from bam
    stats_from_bam -o align.bamstats -s align.bam.summary -t ${params.align_threads} align.bam
    """
}

process coverageCalc {
      label "wfflu"
      cpus 2
      input:
          tuple val(meta), path("align.bam"), path("align.bam.bai")
      output:
          tuple val(meta), path("depth.txt")
      """
      samtools depth -aa align.bam -Q 20 -q 1 > depth.txt
      """

}

process downSample {
    label 'wfflu'
    cpus 2
    input:
        tuple val(meta), path("align.bam"), path("align.bam.bai")
        path "reference.fasta"
    output:
        tuple val(meta), path("all_merged.sorted.bam"), path("all_merged.sorted.bam.bai"), emit: alignments
    """
    # split bam
    header_count=`samtools view -H align.bam | wc -l`
    lines=\$(( ${params.downsample} + \$header_count + 1 ))

    # get the regions from the fasta
    awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length(\$0)}END{print seqlen}' reference.fasta | tr "\\n" "," | tr '>' '\\n' | sed 's/,\$//g' | grep ',' > regions.txt

    # for every region we're going to downsample separatley

    while read region_string;
    do

      region=`echo \${region_string} | cut -f1 -d','`
      length=`echo \${region_string} | cut -f2 -d','`

      # get upper and lower bounds of reference span

      upper=`echo \$((\${length}+(\${length}*10/100)))`
      lower=`echo \$((\${length}-(\${length}*10/100)))`

      samtools view -H align.bam \${region} > head.sam

      # filter reads in region and covering region

      samtools view align.bam \${region} | awk -v upper="\${upper}" -v lower="\${lower}" '{if(length(\$10) < upper && length(\$10) > lower) print \$0}' > \${region}.tmp.sam;

      # ignore regions with no reads
      count=`wc -l < \${region}.tmp.sam`

      if [ "\${count}" -eq "0" ];
      then
        echo "no reads in \${region} so continuing"
        cat head.sam | samtools view -bh > \${region}_all.bam
        continue;
      fi

      cat head.sam \${region}.tmp.sam | samtools view -bh > \${region}.bam

      samtools view -h -F16 \${region}.bam > \${region}_fwd.sam;
      head -\${lines} \${region}_fwd.sam | samtools view -bh - > \${region}_fwd.bam;

      samtools view -h -f16 \${region}.bam > \${region}_rev.sam;
      head -\${lines} \${region}_rev.sam | samtools view -bh - > \${region}_rev.bam;
      samtools merge \${region}_all.bam \${region}_fwd.bam \${region}_rev.bam;

    done < regions.txt

    samtools merge all_merged.bam *_all.bam
    samtools sort all_merged.bam > all_merged.sorted.bam
    samtools index all_merged.sorted.bam
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
        tuple val(meta), path("downsample.bam"), path("downsample.bam.bai"), val(medaka_model)
        path("reference.fasta")
    output:
        tuple val(meta), path("variants.annotated.filtered.vcf")
    script:
    """
    medaka consensus downsample.bam consensus.hdf --model "${medaka_model}"
    medaka variant --gvcf reference.fasta consensus.hdf variants.vcf --verbose
    medaka tools annotate --debug --pad 25 variants.vcf reference.fasta downsample.bam variants.annotated.vcf

    bcftools filter -e "ALT='.'" variants.annotated.vcf | bcftools filter -o variants.annotated.filtered.vcf -O v -e "INFO/DP<${params.min_coverage}" -
    """
}


process makeConsensus {
    label "wfflu"
    cpus 2
    input:
        tuple val(meta), path("variants.annotated.filtered.vcf"), path("depth.txt")
        path "reference.fasta"
    output:
        tuple val(meta), path("draft.consensus.fasta")
    """
    awk '{if (\$3<${params.min_coverage}) print \$1"\t"\$2+1}' depth.txt > mask.regions
    bgzip variants.annotated.filtered.vcf
    tabix variants.annotated.filtered.vcf.gz

    bcftools consensus --mask mask.regions  --mark-del '-' --mark-ins lc --fasta-ref reference.fasta -o draft.consensus.fasta variants.annotated.filtered.vcf.gz
    """
}

process typeFlu {
    label "wfflutyping"
    cpus 2
    input:
        tuple val(meta), path("consensus.fasta")
        path("blast_db")
    output:
        tuple val(meta), path("insaflu.typing.txt"), emit: typing
        path "abricate.version", emit: version
    """
    abricate --version | sed 's/ /,/' > abricate.version
    abricate --datadir blast_db --db insaflu -minid 70 -mincov 60 --quiet consensus.fasta 1> insaflu.typing.txt
    """
}

process processType {
    label "wfflu"
    cpus 1
    input:
        tuple val(meta), path("insaflu.typing.txt")
    output:
        tuple val(meta), path("processed_type.json")
    script:
    """
    workflow-glue process_abricate --typing insaflu.typing.txt --output processed_type.json
    """
}


process prepNextclade {
    label "wfflu"
    input:
        tuple val(meta), path("typing.json"), path("consensus.fasta")
        path "nextclade_data"
    output:
        path "datasets", optional: true
    script:
    String alias = meta.alias
    """
    mkdir datasets
    workflow-glue nextclade_helper \
        --typing typing.json \
        --nextclade_datasets "nextclade_data" \
        --consensus "consensus.fasta" \
        --sample_alias "${alias}"
    
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

process collectFilesInDir {
    label "wfflu"
    cpus 1
    input: tuple val(meta), path("staging_dir/*"), val(dirname)
    output: tuple val(meta), path(dirname)
    script:
    """
    mv staging_dir $dirname
    """
}

process makeReport {
    label "wfflu"
    input:
        val metadata
        path "fastcat_stats/?.gz"
        path "data/*"
        path "nextclade/*"
        path nextclade_datasets
        path "versions/*"
        path "params.json"
    output:
        tuple path("wf-flu-*.html"), path("wf-flu-results.csv")
    script:
        report_name = "wf-flu-report.html"
        def metadata = new JsonBuilder(metadata).toPrettyString()
    """
    echo '${metadata}' > metadata.json
    workflow-glue report "${report_name}"\
        --data data \
        --stats fastcat_stats \
        --versions versions \
        --params params.json \
        --nextclade_files nextclade/* \
        --nextclade_datasets "${nextclade_datasets}" \
        --metadata metadata.json
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636. This is the only way to
// publish files from a workflow whilst decoupling the publish from the process steps.
// The process takes a tuple containing the filename and the name of a sub-directory to
// put the file into. If the latter is `null`, puts it into the top-level directory.
process output {
    // publish inputs to output directory
    label "wfflu"
    cpus 1
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
        alignment = alignReads(samples.map{ meta, reads, stats -> [ meta,reads ] }, reference)
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

        for_draft = variants.join(coverage)

        draft = makeConsensus(for_draft, reference)
        flu_type = typeFlu(draft, blastdb)

        processed_type = processType(flu_type.typing)

        nextclade_prep_input = processed_type.join(draft, remainder: true)

        nextclade_prep = prepNextclade(nextclade_prep_input, nextclade_data)
        
        nextclade_datasets = nextclade_prep
        | map { file(it.resolve("**"), type: "file") }
        | flatten
        | map { [it.parent.name, it] }
        | groupTuple

        nextclade_result = nextclade(nextclade_datasets)

        software_versions = getVersions()
        software_versions = software_versions.mix(flu_type.version.first()) | collectFile()
        workflow_params = getParams()

        // output_alignments = alignment.alignments.map{ it -> return tuple(it[1], it[2]) }

        // get all the per sample results together
        ch_per_sample_results = samples
        | join(coverage)
        | join(flu_type.typing) 
        | join(processed_type)

        // collect results into a directory for the sample directory to avoid collisions
        ch_results_for_report = ch_per_sample_results
        | map {
            meta = it[0]
            rest = it[1..-1]
            [meta, rest, meta.alias]
        }
        | collectFilesInDir
        | map { meta, dirname -> dirname }



        report = makeReport(
            samples.map{it -> it[0]}.toList(),
            samples.map{it -> it[2].resolve("per-read-stats.tsv.gz")}.toList(),
            ch_results_for_report | collect,
            nextclade_result.map{it -> it[1]}.collect(),
            nextclade_data,
            software_versions.collect(),
            workflow_params
        )

        // create channel with files to publish; the channel will have the shape `[file,
        // name of sub-dir to be published in]`.

        ch_to_publish = Channel.empty()
        | mix(
            software_versions | map { [it, null] },
            workflow_params | map { [it, null] },
            report | map { [it[0] , null] },
            report | map { [it[1], null] },
            alignment.alignments
            | map { meta, bam, bai -> [bam, "$meta.alias/alignments"] },
            variants
            | map { meta, vcf -> [vcf, "$meta.alias/variants"]},
            coverage
            | map { meta, depth -> [depth, "$meta.alias/coverage"]},
            draft
            | map { meta, fa -> [fa, "$meta.alias/consensus"]},
            flu_type.typing
            | map { meta, json -> [json, "$meta.alias/typing"]}
        )

    emit: 
        results = ch_to_publish
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    Pinguscript.ping_start(nextflow, workflow, params)

    samples = fastq_ingress([
    "input": params.fastq,
    "stats": true,
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
    pipeline.out.results
    | toList
    | flatMap
    | output
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
