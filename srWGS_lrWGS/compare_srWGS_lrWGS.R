######################################################################################
## Script to compare Illumina and SV calls from all algorithms and save them to RDS
## Adapted from ICC scripts by HE
## 12/12/2024
######################################################################################

setwd("/nfs/research/icortes/DATA/SAVANA_manuscript/revision/benchmarking_code/srWGS_lrWGS")

# load libraries
library(StructuralVariantAnnotation)

# read in rds objects of variants
svs_illumina <- readRDS("rds/illumina_svs.rds")
cnas_illumina <- readRDS("rds/illumina_cnas.rds")
svs_lr <- readRDS("rds/SVs_tumour_lrWGS_bps.rds")
  
sample_info = read.table("../nanopore_samples_metadata_age.tsv", header=T, sep="\t")

# remove qc failures Illumina: 
# - HMF QC STATUS FAIL
# - HMF QC STATUS WARN + < 50 high quality srWGS SVs
illumina_failures <- c("GBM-004","GBM-010","GBM-012","GBM-014","SARC-022_R2","SARC-050_R1","SARC-050_R2","SARC-050_R3")

# shorten the REGION in sample_id for readability
svs_lr$sample_id <- gsub("REGION-", "R", svs_lr$sample)
svs_lr = subset(svs_lr, select = -sample)
svs_illumina$sample_id <- gsub("REGION-", "R", svs_illumina$sample_id)
cnas_illumina$sample_id <- gsub("REGION-", "R", cnas_illumina$sample_id)
sample_info$alias_id <- gsub("REGION-", "R", sample_info$alias_id)
# remove that failed Illumina QC
samples = unique(sample_info$alias_id)
samples = samples[! samples %in% illumina_failures]

# get a list of the callers
callers <- unique(svs_lr$caller)

# compare/label the SVs
svs_illumina_validation = c()
sbd_illumina_validation = c()
cna_illumina_validation = c()
svs_lr_validation = c()

buffer = 500
for (sample in samples) {
  illumina_sample = svs_illumina[which(svs_illumina$sample_id == sample), ]
  # remove the single breakends from this set - compare separately
  illumina_sample = illumina_sample[which(illumina_sample$svtype != "SBND"),]
  # get the single-breakends from Illumina
  illumina_sbnd_sample = svs_illumina[which(svs_illumina$sample_id == sample),]
  illumina_sbnd_sample = illumina_sbnd_sample[which(illumina_sbnd_sample$svtype ==
                                                      "SBND"),]
  # get the illumina cnas
  illumina_cna_sample = cnas_illumina[which(cnas_illumina$sample_id == sample),]
  print(paste0("# Illumina SVs for Donor ", sample, ": ", nrow(illumina_sample)))
  print(paste0(
    "# Illumina SBNDS for Donor ",
    sample,
    ": ",
    nrow(illumina_sbnd_sample)
  ))
  print(paste0(
    "# Illumina CNAs for Donor ",
    sample,
    ": ",
    nrow(illumina_cna_sample)
  ))
  for (caller in callers) {
    # now get the lr SVs for each caller
    lr_svs_sample = svs_lr[which(svs_lr$sample_id == sample &
                                   svs_lr$caller == caller),]
    print(paste0(
      "# LR SVs for Donor ",
      sample,
      " from ",
      caller,
      ": ",
      nrow(lr_svs_sample)
    ))
    
    if (nrow(lr_svs_sample)) {
      lr_svs_sample$strand1 = as.vector(lr_svs_sample$strand1)
      lr_svs_sample$strand2 = as.vector(lr_svs_sample$strand2)
      
      # START VS START
      # CEILING VS CEILING
      ## ILLUMINA SVs
      illumina_sample$key = paste(illumina_sample$chrom1,
                                  ceiling((illumina_sample$start1 / buffer)) *
                                    buffer,
                                  illumina_sample$sample_id)
      lr_svs_sample$key = paste(lr_svs_sample$chrom1,
                                ceiling((lr_svs_sample$start1 / buffer)) *
                                  buffer,
                                lr_svs_sample$sample_id)
      # add a column for the caller's validation status for the illumina SVs
      illumina_sample[, paste0(caller, "_MATCH")] <- FALSE
      illumina_sample[, paste0(caller, "_MATCH")] <-
        ifelse(illumina_sample$key %in% unique(lr_svs_sample$key),
               TRUE,
               illumina_sample[, paste0(caller, "_MATCH")])
      # update the lr SVs
      lr_svs_sample[, "ILLUMINA_MATCH"] <- FALSE
      lr_svs_sample[, "ILLUMINA_MATCH"] <-
        ifelse(lr_svs_sample$key %in% unique(illumina_sample$key),
               TRUE,
               lr_svs_sample[, "ILLUMINA_MATCH"])
      ## ILLUMINA SINGLE-BREAKENDS
      illumina_sbnd_sample$key = paste(illumina_sbnd_sample$chrom1,
                                       ceiling((illumina_sbnd_sample$start1 /
                                                  buffer)) * buffer,
                                       illumina_sbnd_sample$sample_id)
      # add a column for the validation status for the illumina sbnds
      illumina_sbnd_sample[, paste0(caller, "_MATCH")] <- FALSE
      illumina_sbnd_sample[, paste0(caller, "_MATCH")] <- ifelse(
        illumina_sbnd_sample$key %in% unique(lr_svs_sample$key),
        TRUE,
        illumina_sbnd_sample[, paste0(caller, "_MATCH")]
      )
      # update the lr SVs
      lr_svs_sample[, "ILLUMINA_MATCH_SBND"] <- FALSE
      lr_svs_sample[, "ILLUMINA_MATCH_SBND"] <-
        ifelse(lr_svs_sample$key %in% unique(illumina_sbnd_sample$key),
               TRUE,
               lr_svs_sample[, "ILLUMINA_MATCH_SBND"])
      ## ILLUMINA CNAS
      illumina_cna_sample$key = paste(
        illumina_cna_sample$chromosome,
        ceiling((illumina_cna_sample$end / buffer)) *
          buffer,
        illumina_cna_sample$sample_id
      )
      # add a column for the validation status of the illumina cnas
      illumina_cna_sample[, paste0(caller, "_MATCH")] <- FALSE
      illumina_cna_sample[, paste0(caller, "_MATCH")] <- ifelse(
        illumina_cna_sample$key %in% unique(lr_svs_sample$key),
        TRUE,
        illumina_cna_sample[, paste0(caller, "_MATCH")]
      )
      # update the lr SVs
      lr_svs_sample[, "ILLUMINA_MATCH_CNA"] <- FALSE
      lr_svs_sample[, "ILLUMINA_MATCH_CNA"] <-
        ifelse(lr_svs_sample$key %in% unique(illumina_cna_sample$key),
               TRUE,
               lr_svs_sample[, "ILLUMINA_MATCH_CNA"])
      # for CNAs we also need to compare the second coordinates of SV
      lr_svs_sample$key = paste(lr_svs_sample$chrom2,
                                ceiling((lr_svs_sample$start2 / buffer)) *
                                  buffer,
                                lr_svs_sample$sample_id)
      # update cnas
      illumina_cna_sample[, paste0(caller, "_MATCH")] <- ifelse(
        illumina_cna_sample$key %in% unique(lr_svs_sample$key),
        TRUE,
        illumina_cna_sample[, paste0(caller, "_MATCH")]
      )
      # update the lr SVs
      lr_svs_sample[, "ILLUMINA_MATCH_CNA"] <-
        ifelse(lr_svs_sample$key %in% unique(illumina_cna_sample$key),
               TRUE,
               lr_svs_sample[, "ILLUMINA_MATCH_CNA"])
      
      # CEILING VS FLOOR
      illumina_sample$key = paste(illumina_sample$chrom1,
                                  ceiling((illumina_sample$start1 / buffer)) *
                                    buffer,
                                  illumina_sample$sample)
      lr_svs_sample$key = paste(lr_svs_sample$chrom1,
                                floor((lr_svs_sample$start1 / buffer)) * buffer,
                                lr_svs_sample$sample)
      # add a column for the caller's validation status for the illumina SVs
      illumina_sample[, paste0(caller, "_MATCH")] <-
        ifelse(illumina_sample$key %in% unique(lr_svs_sample$key),
               TRUE,
               illumina_sample[, paste0(caller, "_MATCH")])
      # update the lr SVs
      lr_svs_sample[, "ILLUMINA_MATCH"] <-
        ifelse(lr_svs_sample$key %in% unique(illumina_sample$key),
               TRUE,
               lr_svs_sample[, "ILLUMINA_MATCH"])
      # ALSO CHECK ILLUMINA SINGLE-BREAKENDS
      illumina_sbnd_sample$key = paste(illumina_sbnd_sample$chrom1,
                                       ceiling((illumina_sbnd_sample$start1 /
                                                  buffer)) * buffer,
                                       illumina_sbnd_sample$sample)
      # add a column for the caller's validation status for the illumina SVs
      illumina_sbnd_sample[, paste0(caller, "_MATCH")] <- ifelse(
        illumina_sbnd_sample$key %in% unique(lr_svs_sample$key),
        TRUE,
        illumina_sbnd_sample[, paste0(caller, "_MATCH")]
      )
      # update the lr SVs
      lr_svs_sample[, "ILLUMINA_MATCH_SBND"] <-
        ifelse(lr_svs_sample$key %in% unique(illumina_sbnd_sample$key),
               TRUE,
               lr_svs_sample[, "ILLUMINA_MATCH_SBND"])
      ## ILLUMINA CNAS
      illumina_cna_sample$key = paste(
        illumina_cna_sample$chromosome,
        ceiling((illumina_cna_sample$end / buffer)) *
          buffer,
        illumina_cna_sample$sample_id
      )
      # add a column for the validation status of the illumina cnas
      illumina_cna_sample[, paste0(caller, "_MATCH")] <- ifelse(
        illumina_cna_sample$key %in% unique(lr_svs_sample$key),
        TRUE,
        illumina_cna_sample[, paste0(caller, "_MATCH")]
      )
      # update the lr SVs
      lr_svs_sample[, "ILLUMINA_MATCH_CNA"] <-
        ifelse(lr_svs_sample$key %in% unique(illumina_cna_sample$key),
               TRUE,
               lr_svs_sample[, "ILLUMINA_MATCH_CNA"])
      # for CNAs we also need to compare the second coordinates of SV
      lr_svs_sample$key = paste(lr_svs_sample$chrom2,
                                floor((lr_svs_sample$start2 / buffer)) * buffer,
                                lr_svs_sample$sample)
      # update cnas
      illumina_cna_sample[, paste0(caller, "_MATCH")] <- ifelse(
        illumina_cna_sample$key %in% unique(lr_svs_sample$key),
        TRUE,
        illumina_cna_sample[, paste0(caller, "_MATCH")]
      )
      # update the lr SVs
      lr_svs_sample[, "ILLUMINA_MATCH_CNA"] <-
        ifelse(lr_svs_sample$key %in% unique(illumina_cna_sample$key),
               TRUE,
               lr_svs_sample[, "ILLUMINA_MATCH_CNA"])
      # FLOOR VS CEILING
      illumina_sample$key = paste(illumina_sample$chrom1,
                                  floor((illumina_sample$start1 / buffer)) *
                                    buffer,
                                  illumina_sample$sample)
      lr_svs_sample$key = paste(lr_svs_sample$chrom1,
                                ceiling((lr_svs_sample$start1 / buffer)) *
                                  buffer,
                                lr_svs_sample$sample)
      # add a column for the caller's validation status for the illumina SVs
      illumina_sample[, paste0(caller, "_MATCH")] <-
        ifelse(illumina_sample$key %in% unique(lr_svs_sample$key),
               TRUE,
               illumina_sample[, paste0(caller, "_MATCH")])
      # update the lr SVs
      lr_svs_sample[, "ILLUMINA_MATCH"] <-
        ifelse(lr_svs_sample$key %in% unique(illumina_sample$key),
               TRUE,
               lr_svs_sample[, "ILLUMINA_MATCH"])
      # ALSO CHECK ILLUMINA SINGLE-BREAKENDS
      illumina_sbnd_sample$key = paste(illumina_sbnd_sample$chrom1,
                                       floor((illumina_sbnd_sample$start1 /
                                                buffer)) * buffer,
                                       illumina_sbnd_sample$sample)
      lr_svs_sample$key = paste(lr_svs_sample$chrom1,
                                ceiling((lr_svs_sample$start1 / buffer)) *
                                  buffer,
                                lr_svs_sample$sample)
      # add a column for the caller's validation status for the illumina SVs
      illumina_sbnd_sample[, paste0(caller, "_MATCH")] <- ifelse(
        illumina_sbnd_sample$key %in% unique(lr_svs_sample$key),
        TRUE,
        illumina_sbnd_sample[, paste0(caller, "_MATCH")]
      )
      # update the lr SVs
      lr_svs_sample[, "ILLUMINA_MATCH_SBND"] <-
        ifelse(lr_svs_sample$key %in% unique(illumina_sbnd_sample$key),
               TRUE,
               lr_svs_sample[, "ILLUMINA_MATCH_SBND"])
      ## ILLUMINA CNAS
      illumina_cna_sample$key = paste(
        illumina_cna_sample$chromosome,
        floor((illumina_cna_sample$end / buffer)) *
          buffer,
        illumina_cna_sample$sample_id
      )
      # add a column for the validation status of the illumina cnas
      illumina_cna_sample[, paste0(caller, "_MATCH")] <- ifelse(
        illumina_cna_sample$key %in% unique(lr_svs_sample$key),
        TRUE,
        illumina_cna_sample[, paste0(caller, "_MATCH")]
      )
      # update the lr SVs
      lr_svs_sample[, "ILLUMINA_MATCH_CNA"] <-
        ifelse(lr_svs_sample$key %in% unique(illumina_cna_sample$key),
               TRUE,
               lr_svs_sample[, "ILLUMINA_MATCH_CNA"])
      # for CNAs we also need to compare the second coordinates of SV
      lr_svs_sample$key = paste(lr_svs_sample$chrom2,
                                ceiling((lr_svs_sample$start2 / buffer)) *
                                  buffer,
                                lr_svs_sample$sample)
      # update cnas
      illumina_cna_sample[, paste0(caller, "_MATCH")] <- ifelse(
        illumina_cna_sample$key %in% unique(lr_svs_sample$key),
        TRUE,
        illumina_cna_sample[, paste0(caller, "_MATCH")]
      )
      # update the lr SVs
      lr_svs_sample[, "ILLUMINA_MATCH_CNA"] <-
        ifelse(lr_svs_sample$key %in% unique(illumina_cna_sample$key),
               TRUE,
               lr_svs_sample[, "ILLUMINA_MATCH_CNA"])
      # add the caller's validated svs to full df
      svs_lr_validation = rbind(svs_lr_validation, lr_svs_sample)
    }
  }
  # add the validated illumina calls to full dfs
  svs_illumina_validation = rbind(svs_illumina_validation, illumina_sample)
  sbd_illumina_validation = rbind(sbd_illumina_validation, illumina_sbnd_sample)
  cna_illumina_validation = rbind(cna_illumina_validation, illumina_cna_sample)
}

saveRDS(svs_lr_validation, file="rds/lrwgs_svs_validation.rds")
saveRDS(svs_illumina_validation, file="rds/illumina_svs_validation.rds")
saveRDS(sbd_illumina_validation, file="rds/illumina_sbds_validation.rds")
saveRDS(cna_illumina_validation, file="rds/illumina_cna_validation.rds")


