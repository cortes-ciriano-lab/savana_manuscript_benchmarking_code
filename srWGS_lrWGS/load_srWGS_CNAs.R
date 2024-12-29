######################################################################################
## Script to load in cna calls from Illumina
## 20/09/2024
######################################################################################

setwd("/nfs/research/icortes/helrick/projects/thesis/Analysis_Chapter/IlluminaComparison")

sample_info = read.table("/nfs/research/icortes/DATA/SAVANA_manuscript/nanopore_samples_metadata.tsv", header=T, sep="\t")
samples = unique(sample_info$donor_id)

illumina_cnas = c()

for (sample in samples){
  print(sample)
  alias_id = sample_info$alias_id[which(sample_info$donor_id == sample)]
  cohort = sample_info$cohort[which(sample_info$donor_id == sample)]
  purple_path = paste0("../../data_common/donor_purple_results/", sample)
  cna_tsv = list.files(path=paste0(purple_path), full.names=T, pattern=".purple.cnv.somatic.tsv$")
  cn = read_tsv(cna_tsv, show_col_types = FALSE)
  cn = cn[which(cn$chromosome %in% paste0("chr",c(1:22,"X"))),]
  cn.f = cn %>%
    mutate(copyNumber=ifelse((end-start) < 1000, lag(copyNumber), copyNumber),
           start=ifelse((end-start) < 1000, lag(start), start)) %>%
    dplyr::select(chromosome, start, end, copyNumber) %>%
    group_by(chromosome, start, copyNumber) %>%
    dplyr::summarise(end=max(end)) %>%
    mutate(copyNumber=round(copyNumber)) %>%
    ungroup() %>%
    arrange(chromosome, start) %>%
    mutate(new_segment = row_number() == 1 | !(chromosome == lag(chromosome) & copyNumber == lag(copyNumber))) %>%
    mutate(segment = cumsum(new_segment)) %>%
    group_by(segment) %>%
    summarize(
      chromosome = dplyr::first(chromosome),
      start = dplyr::first(start),
      end = dplyr::last(end),
      copyNumber = dplyr::first(copyNumber)
    )
  cna_data <- cn.f %>%
    dplyr::select(chromosome, start, end, copyNumber) %>%
    dplyr::mutate(sample_id=alias_id)
  illumina_cnas = rbind(illumina_cnas, cna_data)
}

illumina_cnas <- as.data.frame(illumina_cnas) # convert to dataframe
saveRDS(illumina_cnas, file="rds/illumina_cnas.rds")
