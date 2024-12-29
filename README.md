## Benchmarking Code for SAVANA Manuscript

Each folder contains scripts to `load` the SVs into an RDS and (if applicable) `compare` SVs from those RDS objects to each other

### tumour_replicates

Comparison between split tumour bams (_in-silico_ replicates). VCFs have been labelled with their replicate calls using `savana evaluate` functionality, allowing a 100bp overlap (INFO column =SOMATIC/NOT_IN_COMPARISON).

* `load_tumour_replicates.R` -> saves data from all callers into `rds/SVs_tumour_splits_svs.rds`

### srWGS_lrWGS

Comparison between srWGS (GRIDSS) and lrWGS from all algorithms

* load_srWGS_SVs.R -> saves srWGS (Illumina+GRIDSS) svs into `rds/illumina_svs.rds`
* load_srWGS_CNAs.R -> saves srWGS (Illumina+PURPLE) CNAs into `rds/illumina_cnas.rds`
* load_lrWGS_SVs_all_algorithms.R -> saves lrWGS (all callers) into `rds/SVs_tumour_lrWGS_svs.rds`
* compare_srWGS_lrWGS.R -> uses `rds/illumina_svs.rds` and `rds/SVs_tumour_lrWGS_svs.rds` to label both srWGS and lrWGS SVs and CNAs into:
	* `rds/lrwgs_svs_validation.rds` -> all lrWGS algorithms SVs, with columns for "ILLUMINA_MATCH", "ILLUMINA_MATCH_SBND", and "ILLUMINA_MATCH_CNA"
	* `rds/illumina_svs_validation.rds` -> contains a "caller"_MATCH validation column for every caller, set to TRUE if lrWGS SV from that caller within 500bp
	* `rds/illumina_sbds_validation.rds` -> same as above but for srWGS single breakends
	* `rds/illumina_cna_validation.rds` -> same as above but for srWGS CNAs

### normal_replicates

Comparison between split normal bams (replicates) from all callers

* load_normal_replicates.R -> loads and saves labelled VCFs from all callers (INFO column =SOMATIC/NOT_IN_COMPARISON) into `rds/SVs_normal_splits_svs.rds`

