If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-flu/issues) page or start a discussion on the [community](https://nanoporetech.com/support).

_Why does the workflow fail, or the report shows very low coverage?_

This can happen when users use the workflow on data that has been generate using the RBK protocol instead of the recomended [Influenza whole-genome protocol](https://community.nanoporetech.com/docs/prepare/library_prep_protocols/ligation-sequencing-influenza-whole-genome), as a result of RBK's shorter read lengths. Ensure the --rbk flag has been set to prevent over-filtering of reads.

 
