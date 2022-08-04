## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[conda](https://docs.conda.io/en/latest/miniconda.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker of conda is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit out website](https://labs.epi2me.io/wfindex).

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-flu --help
```

to see the options for the workflow.

**Workflow outputs**

The workflow outputs several files that are useful for interpretation and analysis:

* Per run:
  * `wf-flu-report.html`: Easy to use HTML report for all samples on the run
  * `wf-flu-results.csv`: Typing results in CSV format for onward processing
* Per sample:
  * `<SAMPLE_NAME>.stats`: Read stats
  * `<SAMPLE_NAME>.bam`: Alignment of reads to reference
  * `<SAMPLE_NAME>.bam.bai`: BAM index
  * `<SAMPLE_NAME>.annotate.filtered.vcf`: medaka called variants
  * `<SAMPLE_NAME>.draft.consensus.fasta`: Consensus FASTA
  * `<SAMPLE_NAME>.insaflu.typing.txt`: abricate typing results
