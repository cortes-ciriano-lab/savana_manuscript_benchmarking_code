######################################################################
## Script to load lrWGS SV calls
## Adapted from ICC by HE
## 20/10/2024
######################################################################
setwd(
  "/nfs/research/icortes/DATA/SAVANA_manuscript/revision/benchmarking_code/srWGS_lrWGS"
)
#--------------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(plyr)
library(StructuralVariantAnnotation)

sample_info = read.table(
  "/nfs/research/icortes/helrick/projects/thesis/data_common/nanopore_samples_metadata_age.tsv",
  header = T,
  sep = "\t"
)

#--------------------------------------------------------------------------------------------------------------------------
# Replicate analysis
#--------------------------------------------------------------------------------------------------------------------------
samples = unique(sample_info$donor_id)
length(samples)

callers = c(
  "cuteSV-2.1.0",
  "savana-1.2.6",
  "nanomonsv-0.5.0",
  "svim-1.4.2",
  "severus-1.0",
  "sniffles-2.2.0",
  "svision-pro"
)

#------------------------------------

all = c()
annotation_info_fields = c(
  "SVTYPE",
  "SVLEN",
  "AF",
  "insLen",
  "MICROSATELLITE",
  "REPEAT",
  "CENTROMERE",
  "TELOMERE",
  "GENE",
  "SVINSLEN"
)

for (sample in samples) {
  for (caller in callers) {
    cohort = sample_info$cohort[which(sample_info$donor_id == sample)]
    print(paste(caller, sample, cohort))
    bedpe = c()
    if (caller == "nanomonsv-0.5.0") {
      path_now = paste0(
        "/nfs/research/icortes/DATA/SAVANA_manuscript/tool_benchmark/",
        caller,
        "-sbnd",
        "/results/",
        cohort,
        "/",
        sample,
        "/tumour"
      )
      file_nows = list.files(path = path_now,
                             pattern = "Genes.vcf",
                             full.names = T)
      if (length(file_nows) > 1) {
        print('More than one file matched at:')
        print(file_nows)
      }
      file_now = file_nows[1]
      print(file_now)
      if (is.na(file_now) | !file.exists(file_now)) {
        print('File is NA or file does not exist with that search pattern')
        print(path_now)
      }
      else{
        vcf <- readVcf(file_now)
        gr <- breakpointRanges(vcf, info_columns = annotation_info_fields)
        bedpe <- data.frame(
          chrom1 = seqnames(gr),
          start1 = start(gr) - 1,
          end1 = end(gr),
          chrom2 = seqnames(partner(gr)),
          start2 = start(partner(gr)) - 1,
          end2 = end(partner(gr)),
          name = names(gr),
          score = gr$QUAL,
          filter = gr$FILTER,
          svtype = gr$svtype,
          svLen = gr$svLen,
          insLen = gr$SVINSLEN,
          strand1 = strand(gr),
          strand2 = strand(partner(gr)),
          REPEAT = gr$REPEAT,
          CENTROMERE = gr$CENTROMERE,
          TELOMERE = gr$TELOMERE,
          GENE = gr$GENE,
          MICROSATELLITE = gr$MICROSATELLITE
        )
        bedpe = bedpe[which(bedpe$filter == "PASS"), ]
        bedpe$caller = caller
        bedpe$sample = sample
        bedpe$cohort = cohort
        bedpe$sample = sample_info$alias_id[match(bedpe$sample, sample_info$donor_id)]
        bedpe$cohort = sample_info$cohort[match(bedpe$cohort, sample_info$cohort)]
        # add read support information from the additional `result.txt` file
        file_now = list.files(path = path_now,
                              pattern = "tumour.nanomonsv.result.txt",
                              full.names = T)
        if (!is.na(file_now)) {
          now = read.table(
            file_now,
            sep = "\t",
            header = T,
            stringsAsFactors = F
          )
          now$AF = now$Supporting_Read_Num_Tumor / now$Checked_Read_Num_Tumor
          # match on the name == SV_ID to add AF and SUPPORT
          bedpe$AF = now$AF[match(sub("^([^_]*_[^_]*)_.*", "\\1", bedpe$name),
                                  # Extract the part of bedpe$name up to the second underscore
                                  now$SV_ID)]
          bedpe$SUPPORT = now$Supporting_Read_Num_Tumor[match(sub("^([^_]*_[^_]*)_.*", "\\1", bedpe$name),
                                                              # Extract the part of bedpe$name up to the second underscore
                                                              now$SV_ID)]
        }
        #------------------------------------------------------------------------------------
        ### single breakends are not reported in the vcf file. They are in:
        ### nanomonsv.sbnd.result.txt
        ### add them here as well
        #------------------------------------------------------------------------------------
        sbnd_file_now = list.files(path = path_now,
                                   pattern = "tumour.nanomonsv.sbnd.result.txt",
                                   full.names = T)
        if (!is.na(sbnd_file_now)) {
          now = read.table(
            sbnd_file_now,
            sep = "\t",
            header = T,
            stringsAsFactors = F
          )
          now$Contig = NULL
          names(now)[1:2] = c("chrom1", "start1")
          now$end1 = now$start1
          now$chrom2 = "."
          now$start2 = "."
          now$end2 = "."
          now$name = now$SV_ID
          now$score = NA
          now$filter = "PASS"
          now$svtype = "SBND"
          now$svLen = NA
          now$insLen = NA
          now$strand1 = now$Dir_1
          now$strand2 = NA
          now$REPEAT = NA
          now$MICROSATELLITE = NA
          now$TELOMERE = NA
          now$CENTROMERE = NA
          now$GENE = NA
          now$SUPPORT = now$Supporting_Read_Num_Tumor
          now$AF = now$Supporting_Read_Num_Tumor / now$Checked_Read_Num_Tumor
          now$end1 = now$start1
          # null out un-needed fields
          now$SV_ID = NULL
          now$Dir_1 = NULL
          now$Checked_Read_Num_Tumor = NULL
          now$Supporting_Read_Num_Tumor = NULL
          now$Checked_Read_Num_Control = NULL
          now$Supporting_Read_Num_Control = NULL
          now$caller = caller
          now$sample = sample
          now$cohort = cohort
          # need to convert the sample IDs etc because the other callers have been processed and already reformatted
          now$sample = sample_info$alias_id[match(now$sample, sample_info$donor_id)]
          now$cohort = sample_info$cohort[match(now$cohort, sample_info$cohort)]
          nrow(now)
          # add the sbnds to bedpe
          bedpe = rbind(bedpe, now)
        }
        # add the bedpe to all data
        all = rbind(all, bedpe)
      }
    }
    if (caller == "sniffles-2.2.0") {
      path_now = paste0(
        "/nfs/research/icortes/DATA/SAVANA_manuscript/tool_benchmark/",
        caller,
        "/results/",
        cohort,
        "/",
        sample
      )
      file_nows = list.files(path = path_now,
                             pattern = "Genes.vcf",
                             full.names = T)
      if (length(file_nows) > 1) {
        print('More than one file matched at:')
        print(file_nows)
      }
      file_now = file_nows[1]
      print(file_now)
      if (is.na(file_now) | !file.exists(file_now)) {
        print('File is NA or file does not exist with that search pattern')
        print(path_now)
      }
      else{
        vcf <- readVcf(file_now)
        gr <- breakpointRanges(
          vcf,
          info_columns = c(annotation_info_fields, "SUPPORT", "AF"),
          inferMissingBreakends = TRUE
        )
        bedpe <- data.frame(
          chrom1 = seqnames(gr),
          start1 = start(gr) - 1,
          end1 = end(gr),
          chrom2 = seqnames(partner(gr)),
          start2 = start(partner(gr)) - 1,
          end2 = end(partner(gr)),
          name = names(gr),
          score = gr$QUAL,
          filter = gr$FILTER,
          svtype = gr$svtype   ,
          svLen = gr$svLen ,
          insLen = gr$SVLEN,
          strand1 = strand(gr),
          strand2 = strand(partner(gr)),
          REPEAT = gr$REPEAT,
          CENTROMERE = gr$CENTROMERE,
          TELOMERE = gr$TELOMERE,
          GENE = gr$GENE,
          MICROSATELLITE = gr$MICROSATELLITE,
          SUPPORT = gr$SUPPORT,
          AF = gr$AF
        )
        bedpe = bedpe[which(bedpe$filter == "PASS"), ]
        bedpe$caller = caller
        bedpe$sample = sample
        bedpe$cohort = cohort
        # need to conver the sample IDs etc because the other callers have been processed and already reformatted
        bedpe$sample = sample_info$alias_id[match(bedpe$sample, sample_info$donor_id)]
        bedpe$cohort = sample_info$cohort[match(bedpe$cohort, sample_info$cohort)]
        nrow(bedpe)
        all = rbind(all, bedpe)
      }
    }
    if (caller == "svim-1.4.2") {
      path_now = paste0(
        "/nfs/research/icortes/DATA/SAVANA_manuscript/tool_benchmark/",
        caller,
        "/results/",
        cohort,
        "/",
        sample
      )
      file_nows = list.files(path = path_now,
                             pattern = "Genes.vcf",
                             full.names = T)
      if (length(file_nows) > 1) {
        print('More than one file matched at:')
        print(file_nows)
      }
      file_now = file_nows[1]
      print(file_now)
      if (is.na(file_now) | !file.exists(file_now)) {
        print('File is NA or file does not exist with that search pattern')
        print(path_now)
      }
      else{
        vcf <- readVcf(file_now)
        gr <- breakpointRanges(
          vcf,
          info_columns = c(annotation_info_fields, "SUPPORT"),
          inferMissingBreakends = TRUE
        )
        bedpe <- data.frame(
          chrom1 = seqnames(gr),
          start1 = start(gr) - 1,
          end1 = end(gr),
          chrom2 = seqnames(partner(gr)),
          start2 = start(partner(gr)) - 1,
          end2 = end(partner(gr)),
          name = names(gr),
          score = gr$QUAL,
          filter = gr$FILTER,
          svtype = gr$svtype,
          svLen = gr$svLen ,
          insLen = gr$SVLEN,
          strand1 = strand(gr),
          strand2 = strand(partner(gr)),
          REPEAT = gr$REPEAT,
          CENTROMERE = gr$CENTROMERE,
          TELOMERE = gr$TELOMERE,
          GENE = gr$GENE,
          MICROSATELLITE = gr$MICROSATELLITE,
          SUPPORT = gr$SUPPORT,
          AF = NA
        )
        bedpe = bedpe[which(bedpe$filter == "PASS"), ]
        bedpe$caller = caller
        bedpe$sample = sample
        bedpe$cohort = cohort
        # need to convert the sample IDs etc because the other callers have been processed and already reformatted
        bedpe$sample = sample_info$alias_id[match(bedpe$sample, sample_info$donor_id)]
        bedpe$cohort = sample_info$cohort[match(bedpe$cohort, sample_info$cohort)]
        nrow(bedpe)
        all = rbind(all, bedpe)
      }
    }
    if (caller == "cuteSV-2.1.0") {
      path_now = paste0(
        "/nfs/research/icortes/DATA/SAVANA_manuscript/tool_benchmark/",
        caller,
        "/results/",
        cohort,
        "/",
        sample
      )
      file_nows = list.files(path = path_now,
                             pattern = "Genes.vcf",
                             full.names = T)
      if (length(file_nows) > 1) {
        print('More than one file matched at:')
        print(file_nows)
      }
      file_now = file_nows[1]
      print(file_now)
      if (is.na(file_now) | !file.exists(file_now)) {
        print('File is NA or file does not exist with that search pattern')
        print(path_now)
      }
      else{
        vcf <- readVcf(file_now)
        gr <- breakpointRanges(
          vcf,
          info_columns = c(annotation_info_fields, "RE"),
          inferMissingBreakends = TRUE
        )
        bedpe <- data.frame(
          chrom1 = seqnames(gr),
          start1 = start(gr) - 1,
          end1 = end(gr),
          chrom2 = seqnames(partner(gr)),
          start2 = start(partner(gr)) - 1,
          end2 = end(partner(gr)),
          name = names(gr),
          score = gr$QUAL,
          filter = gr$FILTER,
          svtype = gr$svtype   ,
          svLen = gr$svLen ,
          insLen = gr$SVLEN,
          strand1 = strand(gr),
          strand2 = strand(partner(gr)),
          REPEAT = gr$REPEAT,
          CENTROMERE = gr$CENTROMERE,
          TELOMERE = gr$TELOMERE,
          GENE = gr$GENE,
          MICROSATELLITE = gr$MICROSATELLITE,
          SUPPORT = gr$RE,
          AF = NA
        )
        bedpe = bedpe[which(bedpe$filter == "PASS"), ]
        bedpe$caller = caller
        bedpe$sample = sample
        bedpe$cohort = cohort
        # need to convert the sample IDs etc because the other callers have been processed and already reformatted
        bedpe$sample = sample_info$alias_id[match(bedpe$sample, sample_info$donor_id)]
        bedpe$cohort = sample_info$cohort[match(bedpe$cohort, sample_info$cohort)]
        all = rbind(all, bedpe)
      }
    }
    if (caller == "severus-1.0") {
      path_now = paste0(
        "/nfs/research/icortes/DATA/SAVANA_manuscript/tool_benchmark/",
        caller,
        "/results/",
        cohort,
        "/",
        sample,
        "/somatic_SVs"
      )
      file_nows = list.files(path = path_now,
                             pattern = "Genes.vcf",
                             full.names = T)
      if (length(file_nows) > 1) {
        print('More than one file matched at:')
        print(file_nows)
      }
      file_now = file_nows[1]
      print(file_now)
      if (is.na(file_now) | !file.exists(file_now)) {
        print('File is NA or file does not exist with that search pattern')
        print(path_now)
      }
      else{
        vcf <- readVcf(file_now)
        gr <- breakpointRanges(
          vcf,
          info_columns = c(annotation_info_fields,
                           "SUPPREAD",
                           "HVAF",
                           "VAF",
                           "DV"),
          inferMissingBreakends = TRUE
        )
        bedpe <- data.frame(
          chrom1 = seqnames(gr),
          start1 = start(gr) - 1,
          end1 = end(gr),
          chrom2 = seqnames(partner(gr)),
          start2 = start(partner(gr)) - 1,
          end2 = end(partner(gr)),
          name = names(gr),
          score = gr$QUAL,
          filter = gr$FILTER,
          svtype = gr$svtype   ,
          svLen = gr$svLen ,
          insLen = gr$SVLEN,
          strand1 = strand(gr),
          strand2 = strand(partner(gr)),
          REPEAT = gr$REPEAT,
          CENTROMERE = gr$CENTROMERE,
          TELOMERE = gr$TELOMERE,
          GENE = gr$GENE,
          MICROSATELLITE = gr$MICROSATELLITE
        )
	   	# mate is stored in "MATE_ID" rather than "MATEID" as in VCF spec.
	   	# leads to breakpointRanges creating mates for these "unpartnered" breakends
	   	# these "unpartnered" breakends begin with "svrecord"
	   	# remove all BNDs that start with "svrecord" since they are duplicates
	   	# the mates are present but just not linked correctly
        bedpe = bedpe[which(bedpe$svtype != "BND" |
                              !grepl("svrecord", bedpe$name)), ]
        
        # add GT info
        SUPPORT = geno(vcf)[["DV"]]
        bedpe$SUPPORT = SUPPORT[match(gsub("_bp1|_bp2|_bp3|_bp4", "", bedpe$name) ,
                                      rownames(SUPPORT)), 1]
        # add AF info
        AF = geno(vcf)[["VAF"]]
        bedpe$AF = AF[match(gsub("_bp1|_bp2|_bp3|_bp4", "", bedpe$name) , rownames(AF)), 1]
        bedpe = bedpe[which(bedpe$filter == "PASS"), ]
        # meta info
        bedpe$caller = caller
        bedpe$sample = sample
        bedpe$cohort = cohort
        # need to convert the sample IDs etc because the other callers have been processed and already reformatted
        bedpe$sample = sample_info$alias_id[match(bedpe$sample, sample_info$donor_id)]
        bedpe$cohort = sample_info$cohort[match(bedpe$cohort, sample_info$cohort)]
        all = rbind(all, bedpe)
      }
    }
    if (caller == "svision-pro") {
      path_now = paste0(
        "/nfs/research/icortes/DATA/SAVANA_manuscript/tool_benchmark/",
        caller,
        "/results/",
        cohort,
        "/",
        sample
      )
      file_nows = list.files(path = path_now,
                             pattern = "Genes.vcf",
                             full.names = T)
      if (length(file_nows) > 1) {
        print('More than one file matched at:')
        print(file_nows)
      }
      file_now = file_nows[1]
      print(file_now)
      if (is.na(file_now) | !file.exists(file_now)) {
        print('File is NA or file does not exist with that search pattern')
        print(path_now)
      }
      else{
        vcf <- readVcf(file_now)
        gr <- breakpointRanges(
          vcf,
          info_columns = c(annotation_info_fields,
                           "SUPPORT",
                           "VAF",
                           "BKPS"
                           ),
          inferMissingBreakends = TRUE
        )
        bedpe <- data.frame(
          chrom1 = seqnames(gr),
          start1 = start(gr) - 1,
          end1 = end(gr),
          chrom2 = seqnames(partner(gr)),
          start2 = start(partner(gr)) - 1,
          end2 = end(partner(gr)),
          name = names(gr),
          score = gr$QUAL,
          filter = gr$FILTER,
          svtype = gr$svtype   ,
          svLen = gr$svLen ,
          insLen = gr$SVLEN,
          strand1 = strand(gr),
          strand2 = strand(partner(gr)),
          SUPPORT = gr$SUPPORT,
          AF = gr$VAF,
          REPEAT = gr$REPEAT,
          CENTROMERE = gr$CENTROMERE,
          TELOMERE = gr$TELOMERE,
          GENE = gr$GENE,
          MICROSATELLITE = gr$MICROSATELLITE
        )
	   	# all SVision-pro BNDs have the exact same ID: "0"
	   	# replace svrecord with caller name to distinguish
	   	# otherwise leads to confusion as ID is not unique
        bedpe$name = gsub("svrecord","svision_",bedpe$name) 
        # parse "CSV" breakpoints from BKPS in INFO column - these don't have their own line in VCF
        new_rows <- list()
        bkps <- gr$BKPS
        for (i in seq_along(bkps)) {
          if (length(bkps[[i]]) > 0) { # Check if the element in bkps is non-empty
            string <- unlist(bkps[i])
            split_values <- strsplit(string, " ")
            if (length(string) > 1){
              for (value in string[-1]) {
                split_row <- strsplit(value, "_")[[1]] # first is the type, second is ID? then chrom and pos
                new_row <- bedpe[i, ]
                new_row$svtype = split_row[[1]]
                new_row$chrom1 = split_row[[3]]
                new_row$start1 = split_row[[4]]
                new_row$end1 = split_row[[4]]
                new_row$start2 = split_row[[5]]
                new_row$end2 = split_row[[5]]
                # temporarily set the insLen to the absolute difference of start1 and end2
                # this will later be used for the svLen
                new_row$insLen = abs(as.numeric(split_row[[4]]) - as.numeric(split_row[[5]]))
                new_rows <- rbind(new_rows, new_row)
              }
            }
          }
        }
        bedpe = rbind(new_rows, bedpe)
        bedpe = bedpe[which(bedpe$filter == "PASS"), ]
        # meta info
        bedpe$caller = caller
        bedpe$sample = sample
        bedpe$cohort = cohort
        # need to convert the sample IDs etc because the other callers have been processed and already reformatted
        bedpe$sample = sample_info$alias_id[match(bedpe$sample, sample_info$donor_id)]
        bedpe$cohort = sample_info$cohort[match(bedpe$cohort, sample_info$cohort)]
        all = rbind(all, bedpe)
      }
    }
    if (caller == "savana-1.2.6") {
      path_now = paste0(
        "/nfs/research/icortes/DATA/SAVANA_manuscript/tool_benchmark/",
        caller,
        "/results/cohorts/",
        cohort,
        "/",
        sample
      )
      file_nows = list.files(path = path_now,
                             pattern = glob2rx("*Genes.vcf"),
                             full.names = T)
      if (length(file_nows) > 1) {
        print('More than one file matched at:')
        print(file_nows)
      }
      file_now = file_nows[1]
      print(file_now)
      if (is.na(file_now) | !file.exists(file_now)) {
        print('File is NA or file does not exist with that search pattern')
        print(path_now)
      }
      else{
        vcf <- readVcf(file_now)
        gr <- breakpointRanges(
          vcf,
          info_columns = c(
            annotation_info_fields,
            "TUMOUR_READ_SUPPORT",
            "TUMOUR_AF"
          ),
          inferMissingBreakends = TRUE
        )
        bedpe <- data.frame(
          chrom1 = seqnames(gr),
          start1 = start(gr) - 1,
          end1 = end(gr),
          chrom2 = seqnames(partner(gr)),
          start2 = start(partner(gr)) - 1,
          end2 = end(partner(gr)),
          name = names(gr),
          score = gr$QUAL,
          filter = gr$FILTER,
          svtype = gr$SVTYPE,
          svLen = gr$SVLEN,
          insLen   = gr$insLen,
          strand1 = strand(gr),
          strand2 = strand(partner(gr)),
          CENTROMERE = gr$CENTROMERE,
          TELOMERE = gr$TELOMERE,
          GENE = gr$GENE,
          REPEAT = gr$REPEAT,
          MICROSATELLITE = gr$MICROSATELLITE,
          SUPPORT = gr$TUMOUR_READ_SUPPORT,
          AF = unlist(gr$TUMOUR_AF)[seq(1, length(unlist(gr$TUMOUR_AF)), 2)]
        )
        bedpe = bedpe[which(bedpe$filter == "PASS"), ]
        bedpe$caller = caller
        bedpe$sample = sample
        bedpe$cohort = cohort
        # need to conver the sample IDs etc because the other callers have been processed and already reformatted
        bedpe$sample = sample_info$alias_id[match(bedpe$sample, sample_info$donor_id)]
        bedpe$cohort = sample_info$cohort[match(bedpe$cohort, sample_info$cohort)]
        # now load sbnd (we cneed to do it separately as the bedpe assumes both ends are mapped)
        begr <-
          breakendRanges(
            vcf,
            info_columns = c(
              annotation_info_fields,
              "TUMOUR_READ_SUPPORT",
              "TUMOUR_AF"
            )
          )
        bedpe_sbnd = NULL
        if (nrow(data.frame(begr)) > 0) {
          # when there are no single breakends
          bedpe_sbnd <- data.frame(
            chrom1 = seqnames(begr),
            start1 = start(begr) - 1,
            end1 = end(begr) ,
            chrom2 = ".",
            start2 = ".",
            end2 = ".",
            name = names(begr),
            score = begr$QUAL,
            filter = begr$FILTER,
            svtype = "SBND",
            svLen = begr$svLen ,
            insLen = ".",
            strand1 = strand(begr),
            strand2 = ".",
            REPEAT = begr$REPEAT,
            CENTROMERE = begr$CENTROMERE,
            TELOMERE = begr$TELOMERE,
            GENE = begr$GENE,
            MICROSATELLITE = begr$MICROSATELLITE,
            SUPPORT = begr$TUMOUR_READ_SUPPORT,
            AF = unlist(begr$TUMOUR_AF)[seq(1, length(unlist(begr$TUMOUR_AF)), 2)]
          )
          table(bedpe_sbnd$svtype)
          bedpe_sbnd = bedpe_sbnd[which(bedpe_sbnd$filter == "PASS"), ]
          # need to check here that there are still SBNDs post PASS filter
          if (nrow(bedpe_sbnd) > 0){
            bedpe_sbnd$caller = caller
            bedpe_sbnd$sample = sample
            bedpe_sbnd$cohort = cohort
            # need to convert the sample IDs etc because the other callers have been processed and already reformatted
            bedpe_sbnd$sample = sample_info$alias_id[match(bedpe_sbnd$sample, sample_info$donor_id)]
            bedpe_sbnd$cohort = sample_info$cohort[match(bedpe_sbnd$cohort, sample_info$cohort)]
          }
        }
        if (!is.null(bedpe_sbnd) > 0) {
          all = rbind(all, bedpe, bedpe_sbnd)
        } else{
          all = rbind(all, bedpe)
        }
      }
    }
  }
}

head(all)
table(all$filter)
table(all$caller)

saveRDS(all, file = "rds/SVs_tumour_lrWGS_raw.rds")
#------------------------------------------------------------------------------------
# add annotations
#------------------------------------------------------------------------------------

callers_match = data.frame(
  old = c(
    "cuteSV-2.1.0",
    "savana-1.2.6",
    "nanomonsv-0.5.0",
    "svim-1.4.2",
    "sniffles-2.2.0",
    "severus-1.0",
    "svision-pro"
  ),
  new = c("cuteSV", "SAVANA", "NanomonSV", "SVIM", "Sniffles2", "Severus", "SVision-pro")
)
all$caller = callers_match$new[match(all$caller, callers_match$old)]

# swap the strands (they're reversed for _2 edge of BNDS - both standard and bp2)
all <- all %>%
  dplyr::mutate(
    strand1_tmp=ifelse(grepl("2$", name) & (caller == "Severus" | caller == "SAVANA"), as.character(strand2),as.character(strand1)),
    strand2_tmp=ifelse(grepl("2$", name) & (caller == "Severus" | caller == "SAVANA"), as.character(strand1),as.character(strand2)),
    strand1=strand1_tmp,
    strand2=strand2_tmp) %>%
  dplyr::select(!c(strand1_tmp,strand2_tmp))

all$AF = as.numeric(all$AF)
# only keep canonical chromosomes
chromosomes <- paste0("chr", c(1:22, "X", "Y"))
chromosomes <- c(chromosomes,".") # . for sbnd
all <- all %>%
  filter((chrom1 %in% chromosomes & chrom2 %in% chromosomes))
# relevel
all$chrom1 <- factor(all$chrom1, levels = chromosomes)
all$chrom2 <- factor(all$chrom2, levels = chromosomes)

# svLen is not always correct - calculate manually for non-insertions on same chrom
all$abs_svLen = abs(all$svLen)
all$start1 = as.numeric(all$start1)
all$start2 = as.numeric(all$start2)
idx = which(all$chrom1 == all$chrom2 & all$svtype != "INS")
all$abs_svLen[idx] = abs(all$start1[idx] -  all$start2[idx])

# SVision-Pro does not follow VCF spec, manually fix the location and length
idx = which(all$caller=="SVision-pro" & (all$svLen == 1 | all$svLen == 2))
all$svLen[idx] = all$insLen[idx]
# also keep only the first SVtype (since it's been translated into breakpoints)
# translate to simplified SV types
idx = which(all$caller=="SVision-pro" & grepl("\\+", all$svtype))
all$svtype[idx] = sub("\\+.*", "", all$svtype[idx])
svision_svtypes = data.frame(
  old = c(
    "idDUP",
    "itDUP",
    "tDUP",
    "dDUP"
  ),
  new = c("DUP","DUP","DUP","DUP")
)
idx <- match(all$svtype, svision_svtypes$old)
all$svtype[!is.na(idx)] <- svision_svtypes$new[na.omit(idx)]

# use svtype as basis
all$SV_type = all$svtype
# translate into types
table(all$SV_type)
idx =which(all$SV_type =="BND" & all$strand1 =="+" &  all$strand2 =="-" & all$abs_svLen <= 500)
all$SV_type[idx] = "Deletion <= 500bp"
idx =which(all$SV_type =="BND" & all$strand1 =="+" &  all$strand2 =="-" & all$abs_svLen > 500)
all$SV_type[idx] = "Deletion > 500bp"
idx =which(all$SV_type =="DEL" & all$abs_svLen <= 500)
all$SV_type[idx] = "Deletion <= 500bp"
idx =which(all$SV_type =="DEL"  & all$abs_svLen > 500)
all$SV_type[idx] = "Deletion > 500bp"
idx =which(all$SV_type =="BND" & all$strand1 =="-" &  all$strand2 =="+")
all$SV_type[idx] = "Duplication"
idx =which(all$svtype =="DUP")
all$SV_type[idx] = "Duplication"
idx =which(all$SV_type =="BND" & all$strand1 =="-" &  all$strand2 =="-")
all$SV_type[idx] = "Inversion"
idx =which(all$SV_type =="BND" & all$strand1 =="+" &  all$strand2 =="+")
all$SV_type[idx] = "Inversion"
idx =which(all$SV_type =="INV" )
all$SV_type[idx] = "Inversion"
idx =which(all$svtype =="INS" & all$insLen > 500)
all$SV_type[idx] = "Insertion > 500bp"
idx =which(all$svtype =="INS" & all$insLen <= 500)
all$SV_type[idx] = "Insertion <= 500bp"
idx = which(as.vector(all$chrom1) != as.vector(all$chrom2) )
all$SV_type[idx] = "Translocation"
## Single-breakend
idx = which(all$svtype == "SBND")
all$SV_type[idx] = "Single breakend"
table(all$SV_type, all$caller)

#----
saveRDS(all, file = "rds/SVs_tumour_lrWGS_bps.rds")

# sanity-check: confirm the number of samples is the same
all %>% group_by(caller) %>% dplyr::mutate(count=n_distinct(sample)) %>% dplyr::select(caller,count) %>% distinct()
