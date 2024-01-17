Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| Workflow report | ./wf-flu-report.html | Easy-to-use HTML report for all samples in the run. | aggregated |
| Typing results | ./wf-flu-results.csv | Typing results in CSV format for onward processing. | aggregated |
| Read alignments | ./{{ alias }}/alignments/align.bam | Read allignments per sample in BAM format. | per-sample |
| Draft consensus FASTA | ./{{ alias }}/consensus/draft.consensus.fasta | Draft consensus sequence. | per-sample |
| Read depth | ./{{ alias }}/coverage/depth.txt | Read depth per base. | per-sample |
| Insaflu typing results | ./{{ alias }}/typing/insaflu.typing.txt | Insaflu abricate typing results. | per-sample |
| Variants file | ./{{ alias }}/variants/variants.annotated.filtered.vcf | Called variants in VCF format. | per-sample |
