### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| fastq | string | FASTQ files to use in the analysis. | This accepts one of three cases: (i) the path to a single FASTQ file; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| analyse_unclassified | boolean | Analyse unclassified reads from input directory. By default the workflow will not process reads in the unclassified directory. | If selected and if the input is a multiplex directory the workflow will also process the unclassified directory. | False |


### Sample Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_sheet | string | A CSV file used to map barcodes to sample aliases. The sample sheet can be provided when the input data is a directory containing sub-directories with FASTQ files. | The sample sheet is a CSV file with, minimally, columns named `barcode` and `alias`. Extra columns are allowed. A `type` column is required for certain workflows and should have the following values; `test_sample`, `positive_control`, `negative_control`, `no_template_control`. |  |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all workflow results. |  | output |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| override_basecaller_cfg | string | Override auto-detected basecaller model that processed the signal data; used to select an appropriate Medaka model. | Per default, the workflow tries to determine the basecall model from the input data. This parameter can be used to override the detected value (or to provide a model name if none was found in the inputs). However, users should only do this if they know for certain which model was used as selecting the wrong option might give sub-optimal results. A list of recent models can be found here: https://github.com/nanoporetech/dorado#DNA-models. |  |
| reference | string | Enter the full path to a custom reference genome you would like to use. | The workflow defaults to the IRMA consensus reference. This option allows you to specify a path to an alternative reference. |  |
| blastdb | string | blastdb file used for typing. | The workflow provides the INSaFLU blastdb. If you would like to supply an alternative then provide the full path to the file here. |  |
| min_coverage | integer | Coverage threshold for masking bases in the consensus. | Any bases that are covered below 20x are masked (i.e. represented by 'N') by default in the consensus, this threshold can be changed using this parameter. | 20 |
| min_qscore | number | Minimum read quality score for fastcat. | Any reads which are below quality score of 9 are not used by default. This parameter allows you to customise that. For more information on quality scores please see this blog post: https://labs.epi2me.io/quality-scores | 9 |
| downsample | integer | Number of reads to downsample to in each direction, leave blank for no downsampling. | The workflow for each segment will first filter reads to include only those that are ±10% of the segment length before downsampling to the specified integer (taking an even split from forward and reverse reads). This downsampled data is then used in variant calling. |  |
| rbk | boolean | Set when using data created with the RBK protocol. | This prevents shorter reads being filtered out and also turns off downsampling as this is not appropriate for the shorter reads generated with RBK. | False |


