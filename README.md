# Influenza Typing Workflow

Influenza A&B typing and analysis from Nanopore data.



## Introduction

Influenza is a single-stranded RNA virus and contains a 13.5-14.5kb genome which is split into 8 segments encoding 10-14 proteins (dependent on strain).

The virus is classified using two proteins found on the outer surface of the viral capsid. You’ve probably heard of H1N1 Influenza for example. The H represents hemagglutinin and the N is neuraminidase.

This analysis workflow can be used with Oxford Nanopore Technologies sequencing data from amplified segments of the Influenza Type A and Type B genomes, to determine the most likely strain of Influenza to which the sequenced sample belongs.




## Compute requirements

Recommended requirements:

+ CPUs = 32
+ Memory = 32GB

Minimum requirements:

+ CPUs = 4
+ Memory = 2GB

Approximate run time: 30 minutes when number of cores >= samples

ARM processor support: False




## Install and run

These are instructions to install and run the workflow on command line. You can also access the workflow via the [EPI2ME application](https://labs.epi2me.io/downloads/).

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and software resources, therefore nextflow will need to be installed before attempting to run the workflow.

The workflow can currently be run using either [Docker](https://www.docker.com/products/docker-desktop) or
[singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided 
either docker or singularity is installed. This is controlled by the [`-profile`](https://www.nextflow.io/docs/latest/config.html#config-profiles) parameter as exemplified below.

It is not required to clone or download the git repository in order to run the workflow.
More information on running EPI2ME workflows can be found on our [website](https://labs.epi2me.io/wfindex).

The following command can be used to obtain the workflow. This will pull the repository in to the assets folder of nextflow and provide a list of all parameters available for the workflow as well as an example command:

```
nextflow run epi2me-labs/wf-flu -help
```
A demo dataset is provided for testing of the workflow. It can be downloaded using:
```
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-flu/wf-flu-demo.tar.gz
tar -xzvf wf-flu-demo.tar.gz
```
The workflow can be run with the demo data using:
```
nextflow run epi2me-labs/wf-flu \
--fastq test_data/fastq -profile standard
```
For further information about running a workflow on the cmd line see https://labs.epi2me.io/wfquickstart/





## Related protocols

This workflow is designed to take input sequences that have been produced from [Oxford Nanopore Technologies](https://nanoporetech.com/) devices using this protocol: (https://community.nanoporetech.com/docs/prepare/library_prep_protocols/ligation-sequencing-influenza-whole-genome)
Samples not prepared with this protocol may work sub-optimally or fail to complete succesfully.




## Inputs

### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| fastq | string | FASTQ files to use in the analysis. | This accepts one of three cases: (i) the path to a single FASTQ file; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| basecaller_cfg | string | Name of the model that was used to basecall signal data, used to select an appropriate Medaka model. | The basecaller configuration is used to automatically select the appropriate Medaka model. The automatic selection can be overridden with the 'medaka_variant_model' and 'medaka_consensus_model' parameters. The model list only shows models that are compatible with this workflow. | dna_r10.4.1_e8.2_400bps_hac |
| analyse_unclassified | boolean | Analyse unclassified reads from input directory. By default the workflow will not process reads in the unclassified directory. | If selected and if the input is a multiplex directory the workflow will also process the unclassified directory. | False |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| reference | string | Enter the full path to a custom reference genome you would like to use. | The workflow defaults to the IRMA consensus reference. This option allows you to specify a path to an alternative reference. |  |
| blastdb | string | blastdb file used for typing. | The workflow provides the INSaFLU blastdb. If you would like to supply an alternative then provide the full path to the file here. |  |
| min_coverage | integer | Coverage threshold for masking bases in the consensus. | Any bases that are covered below 20x are masked (i.e. represented by 'N') by default in the consensus, this threshold can be changed using this parameter. | 20 |
| min_qscore | number | Minimum read quality score for fastcat. | Any reads which are below quality score of 9 are not used by default. This parameter allows you to customise that. For more information on quality scores please see this blog post: https://labs.epi2me.io/quality-scores | 9 |
| downsample | integer | Number of reads to downsample to in each direction, leave blank for no downsampling. | By default the workflow will use 130 reads (65 forward, 65 reverse) for typing. However, if you wish to change the number of reads to downsample please specify here. | 130 |
| medaka_consensus_model | string | The name of a Medaka model to use. By default the workflow will select an appropriate Medaka model from the basecaller configuration provided. Entering a name here will override the automated selection and use the Medaka model named here. | The workflow will attempt to map the basecalling model used to a suitable Medaka consensus model. You can override this by providing a model with this option instead. |  |
| rbk | boolean | Set when using data created with the RBK protocol. | This prevents shorter reads being filtered out and also turns off downsampling as this is not appropriate for the shorter reads generated with RBK. | False |


### Miscellaneous Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| disable_ping | boolean | Enable to prevent sending a workflow ping. |  | False |






## Outputs

Outputs files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| Workflow report | ./wf-flu-report.html | Easy-to-use HTML report for all samples in the run. | aggregated |
| Typing results | ./wf-flu-results.csv | Typing results in CSV format for onward processing. | aggregated |
| Read alignments | ./{{ alias }}/alignments/align.bam | Read allignments per sample in BAM format. | per-sample |
| Draft consensus FASTA | ./{{ alias }}/consensus/draft.consensus.fasta | Draft consensus sequence. | per-sample |
| Read depth | ./{{ alias }}/coverage/depth.txt | Read depth per base. | per-sample |
| Insaflu typing results | ./{{ alias }}/typing/insaflu.typing.txt | Insaflu abricate typing results. | per-sample |
| Variants file | ./{{ alias }}/variants/variants.annotated.filtered.vcf | Called variants in VCF format. | per-sample |




## Pipeline overview

1. Concatenate reads and filter out short reads < 200 bases long
2. Align reads to reference with [minimap2](https://github.com/lh3/minimap2)
3. Coverage calculations with [samtools](http://www.htslib.org/))
4. Call variants using [medaka](https://github.com/nanoporetech/medaka) ([medaka blog](https://labs.epi2me.io/notebooks/Introduction_to_how_ONT's_medaka_works.html))
5. Make a (coverage masked) consensus with [bcftools](https://samtools.github.io/bcftools/bcftools.html)
6. Typing using [abricate](https://github.com/tseemann/abricate) with the [insaflu](https://insaflu.insa.pt/) database, containing the following sequences:

| Database|Gene|Accession|Details|
|----|----|----|----|
|insaflu|M1|MK576795|Type_A MK576795 A/England/7821/2019 2019/01/03 7 (MP)
|insaflu|M1|AF100378|Type_B AF100378.1 Influenza B virus B/Yamagata/16/88 segment 7 M1 matrix protein (M) and BM2 protein (BM2) genes, complete cds
|insaflu|HA|FJ966974|H1 FJ966974.1 Influenza A virus (A/California/07/2009(H1N1)) segment 4 hemagglutinin (HA) gene, complete cds
|insaflu|HA|L11142|H2 L11142.1 Influenza A virus (A/Singapore/1/57 (H2N2)) hemagglutinin (HA) gene, complete cds
|insaflu|HA|MK576794|H3 MK576794 A/England/7821/2019 2019/01/03 4 (HA)
|insaflu|HA|AF285883|H4 AF285883.2 Influenza A virus (A/Swine/Ontario/01911-2/99 (H4N6)) segment 4 hemagglutinin (HA) gene, complete cds
|insaflu|HA|EF541403|H5 EF541403.1 Influenza A virus (A/Viet Nam/1203/2004(H5N1)) segment 4 hemagglutinin (HA) gene, complete cds
|insaflu|HA|AB295613|H15 AB295613.1 Influenza A virus (A/duck/Australia/341/83(H15N8)) HA gene for haemagglutinin, complete cds
|insaflu|NA|GQ377078|N1 GQ377078.1 Influenza A virus (A/California/07/2009(H1N1)) segment 6 neuraminidase (NA) gene, complete cds
|insaflu|NA|MK576796|N2 MK576796 A/England/7821/2019 2019/01/03 6 (NA)
|insaflu|NA|AB295614|N8 AB295614.1 Influenza A virus (A/duck/Australia/341/83(H15N8)) NA gene for neuraminidase, complete cds
|insaflu|HA|AY338459|H7 AY338459.1 Influenza A virus (A/Netherlands/219/2003(H7N7)) segment 4 hemagglutinin (HA) gene, complete cds
|insaflu|HA|CY014659|H8 CY014659.1 Influenza A virus (A/turkey/Ontario/6118/1968(H8N4)) segment 4, complete sequence
|insaflu|HA|CY014694|H13 CY014694.1 Influenza A virus (A/gull/Maryland/704/1977(H13N6)) segment 4, complete sequence
|insaflu|HA|CY018765|Yamagata CY018765.1 Influenza B virus (B/Yamagata/16/1988) segment 4, complete sequence
|insaflu|HA|CY103892|H17 CY103892.1 Influenza A virus (A/little yellow-shouldered bat/Guatemala/060/2010(H17N10)) hemagglutinin (HA) gene, complete cds
|insaflu|NA|CY103894|N10 CY103894.1 Influenza A virus (A/little yellow-shouldered bat/Guatemala/060/2010(H17N10)) neuraminidase (NA) gene, complete cds
|insaflu|NA|CY125730|N3v2 CY125730.1 Influenza A virus (A/Mexico/InDRE7218/2012(H7N3)) neuraminidase (NA) gene, complete cds
|insaflu|HA|CY125945|H18 CY125945.1 Influenza A virus (A/flat-faced bat/Peru/033/2010(H18N11)) hemagglutinin (HA) gene, complete cds
|insaflu|NA|CY125947|N11 CY125947.1 Influenza A virus (A/flat-faced bat/Peru/033/2010(H18N11)) neuraminidase-like protein (NA) gene, complete cds
|insaflu|HA|CY130078|H12 CY130078.1 Influenza A virus (A/duck/Alberta/60/1976(H12N5)) hemagglutinin (HA) gene, complete cds
|insaflu|HA|CY130094|H14 CY130094.1 Influenza A virus (A/mallard/Astrakhan/263/1982(H14N5)) hemagglutinin (HA) gene, complete cds
|insaflu|NA|CY130096|N5 CY130096.1 Influenza A virus (A/mallard/Astrakhan/263/1982(H14N5)) neuraminidase (NA) gene, complete cds
|insaflu|HA|DQ376624|H6 DQ376624.1 Influenza A virus (A/chicken/Taiwan/0705/99(H6N1)) hemagglutinin (HA) gene, complete cds
|insaflu|HA|EU293864|H16 EU293864.1 Influenza A virus (A/black-headed gull/Turkmenistan/13/76(H16N3)) hemagglutinin (HA) gene, complete cds
|insaflu|HA|FJ183474|H10 FJ183474.1 Influenza A virus (A/mallard/Bavaria/3/2006(H10N7)) segment 4 hemagglutinin (HA) gene, complete cds
|insaflu|NA|FJ183475|N7 FJ183475.1 Influenza A virus (A/mallard/Bavaria/3/2006(H10N7)) segment 6 neuraminidase (NA) gene, complete cds
|insaflu|NA|GQ907296|N3v1 GQ907296.1 Influenza A virus (A/black headed gull/Mongolia/1756/2006(H16N3)) segment 6 neuraminidase (NA) gene, complete cds
|insaflu|HA|GU052203|H11 GU052203.1 Influenza A virus (A/duck/England/1/1956(H11N6)) segment 4 hemagglutinin (HA) gene, complete cds
|insaflu|NA|KC853765|N9 KC853765.1 Influenza A virus (A/Hangzhou/1/2013(H7N9)) segment 6 neuraminidase (NA) gene, complete cds
|insaflu|HA|KX879589|H9 KX879589.1 Influenza A virus (A/swine/Hong Kong/9/98(H9N2)) segment 4 hemagglutinin (HA) gene, partial cds
|insaflu|HA|M58428|Victoria M58428.1 Influenza B/Victoria/2/87, hemagglutinin (seg 4), RNA
|insaflu|NA|EU429793|N4 EU429793.1 Influenza A virus (A/turkey/Ontario/6118/1968(H8N4)) segment 6 neuraminidase (NA) mRNA, complete cds
|insaflu|NA|EU429795|N6 EU429795.1 Influenza A virus (A/duck/England/1/1956(H11N6)) segment 6 neuraminidase (NA) mRNA, complete cds

7. Clade and lineage assignment using [nextclade](https://clades.nextstrain.org/)



## Troubleshooting

 + If the workflow fails please run it with the demo data set to ensure the workflow itself is working. This will help us determine if the issue is related to the environment, input parameters or a bug.
+ See how to interpret some common nextflow exit codes [here](https://labs.epi2me.io/trouble-shooting/).



## FAQ's

If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-flu/issues) page or start a discussion on the [community](https://nanoporetech.com/support).

_Why does the workflow fail, or the report shows very low coverage?_

This can happen when users use the workflow on data that has been generate using the RBK protocol instead of the recomended [Influenza whole-genome protocol](https://community.nanoporetech.com/docs/prepare/library_prep_protocols/ligation-sequencing-influenza-whole-genome), as a result of RBK's shorter read lengths. Ensure the --rbk flag has been set to prevent over-filtering of reads.

 




## Related blog posts

 + [Importing third-party workflows into EPI2ME Labs](https://labs.epi2me.io/nexflow-for-epi2melabs/)

See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts.



