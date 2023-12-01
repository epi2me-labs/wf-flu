_Why does the workflow hang after only running validate_sample_sheet and fastcat processes?_

This is likely happening because the user is running the workflow on ARM processors, such as in M1/2 MACs. 
Avoid this by either using a diferent local computer or by running the workflow on cloud.

_Why does the workflow fail, or the report shows very low coverage?_

This can happen when users use the workflow on data that has been generate using the RBK protocol instead of the recomended [Influenza whole-genome protocol](https://community.nanoporetech.com/docs/prepare/library_prep_protocols/ligation-sequencing-influenza-whole-genome), as a result of RBK's shorter read lengths.

 
