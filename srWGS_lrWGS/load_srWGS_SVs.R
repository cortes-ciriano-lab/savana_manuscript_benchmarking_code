######################################################################
## Script to load Illumina SV calls from cohort
## 19/09/2024
######################################################################

setwd("/nfs/research/icortes/DATA/SAVANA_manuscript/revision/benchmarking_code/srWGS_lrWGS")

#----------------------------------------------------
library(StructuralVariantAnnotation)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(stringr)
#----------------------------------------------------

sample_info = read.table("/nfs/research/icortes/DATA/SAVANA_manuscript/nanopore_samples_metadata.tsv", header=T, sep="\t")
samples = unique(sample_info$donor_id)

all_svs = c()

for (sample in samples){
  print(sample)
  alias_id = sample_info$alias_id[which(sample_info$donor_id == sample)]
  cohort = sample_info$cohort[which(sample_info$donor_id == sample)]
  purple_path = paste0("../../data_common/donor_purple_results/", sample)
  sv_vcf_gz = list.files(path=paste0(purple_path), full.names=T, pattern=".purple.sv.RepeatMasker.Blacklist.GridssPoN.microsatellites.Centromeres.Telomeres.Genes.vcf.gz$")
  vcf <- readVcf(sv_vcf_gz)
  # read in breakpoints
  info_fields = c("VF","REF", "REFPAIR","REPEAT","MICROSATELLITE","CENTROMERE","TELOMERE","GENE")
  bpgr = breakpointRanges(vcf, info_columns = info_fields)
  # read in breakends (single breakends)
  begr = breakendRanges(vcf, info_columns = info_fields)
  
  # convert both to bedpe
  bp_bedpe <- data.frame(
    chrom1=seqnames(bpgr),
    start1=start(bpgr),
    end1=end(bpgr),
    chrom2=seqnames(partner(bpgr)),
    start2=start(partner(bpgr)),
    end2=end(partner(bpgr)),
    name=names(bpgr),
    score=bpgr$QUAL,
    filter=bpgr$FILTER,
    svtype= bpgr$svtype,
    svLen = bpgr$svLen,
    insLen   = bpgr$insLen,
    strand1=strand(bpgr),
    strand2=strand(partner(bpgr)),
    AF=bpgr$VF/(bpgr$VF + bpgr$REF + bpgr$REFPAIR),
    REPEAT=bpgr$REPEAT,
    MICROSATELLITE=bpgr$MICROSATELLITE,
    CENTROMERE=bpgr$CENTROMERE,
    TELOMERE=bpgr$TELOMERE,
    GENE=bpgr$GENE,
    LABEL="GRIDSS2 - Illumina"
  )
  be_bedpe <- data.frame(
    chrom1=seqnames(begr),
    start1=start(begr),
    end1=end(begr),
    chrom2=NA,
    start2=NA,
    end2=NA,
    name=names(begr),
    score=begr$QUAL,
    filter=begr$FILTER,
    svtype= "SBND",
    svLen = begr$svLen,
    insLen = begr$insLen,
    strand1=strand(begr),
    strand2=NA,
    AF= begr$VF/(begr$VF + begr$REF + begr$REFPAIR),
    REPEAT=begr$REPEAT,
    MICROSATELLITE=begr$MICROSATELLITE,
    CENTROMERE=begr$CENTROMERE,
    TELOMERE=begr$TELOMERE,
    GENE=begr$GENE,
    LABEL="GRIDSS2 - Illumina"
  )
  # put everything in a big bedpe
  bedpe = rbind(bp_bedpe, be_bedpe)
  # only include PASS variants
  bedpe = bedpe[which(bedpe$filter=="PASS"),]
  # remove those with svLen < 32
  bedpe = bedpe[which(abs(bedpe$svLen) >= 32 | is.na(bedpe$svLen)),]
  # add anonymous id
  bedpe$sample_id = alias_id
  
  # add to meta df
  all_svs = rbind(all_svs, bedpe)
}

# swap the strands (they're reversed for _2 edge)
all_svs <- all_svs %>% dplyr::mutate(
  strand1_tmp=ifelse(grepl("h$", name), as.character(strand2),as.character(strand1)),
  strand2_tmp=ifelse(grepl("h$", name), as.character(strand1),as.character(strand2)),
  strand1=strand1_tmp,
  strand2=strand2_tmp)

# CONVERT/INFER SV TYPE
all_svs$abs_svLen = abs(as.numeric(all_svs$svLen))
all_svs$start1 = as.numeric(all_svs$start1)
all_svs$start2 = as.numeric(all_svs$start2)
all_svs$SV_type = all_svs$svtype
table(all_svs$SV_type)
idx=which(all_svs$SV_type =="BND" & all_svs$strand1 =="+" &  all_svs$strand2 =="-" & all_svs$abs_svLen <= 500)
all_svs$SV_type[idx] = "Deletion <= 500bp"
idx=which(all_svs$SV_type =="BND" & all_svs$strand1 =="+" &  all_svs$strand2 =="-" & all_svs$abs_svLen > 500)
all_svs$SV_type[idx] = "Deletion > 500bp"
idx=which(all_svs$SV_type =="DEL" & all_svs$abs_svLen <= 500)
all_svs$SV_type[idx] = "Deletion <= 500bp"
idx =which(all_svs$SV_type =="DEL"  & all_svs$abs_svLen > 500)
all_svs$SV_type[idx] = "Deletion > 500bp"
idx=which(all_svs$SV_type =="BND" & all_svs$strand1 =="-" &  all_svs$strand2 =="+")
all_svs$SV_type[idx] = "Duplication"
idx=which(all_svs$svtype =="DUP")
all_svs$SV_type[idx] = "Duplication"
idx=which(all_svs$SV_type =="BND" & all_svs$strand1 =="-" & all_svs$strand2 =="-")
all_svs$SV_type[idx] = "Inversion"
idx =which(all_svs$SV_type =="BND" & all_svs$strand1 =="+" & all_svs$strand2 =="+")
all_svs$SV_type[idx] = "Inversion"
idx =which(all_svs$SV_type =="INV" )
all_svs$SV_type[idx] = "Inversion"
idx=which(all_svs$insLen >= (abs(all_svs$svLen)*0.7) & all_svs$insLen > 500)
all_svs$SV_type[idx] = "Insertion > 500bp"
idx=which(all_svs$insLen >= (abs(all_svs$svLen)*0.7) & all_svs$insLen <= 500)
all_svs$SV_type[idx] = "Insertion <= 500bp"
idx = which(as.vector(all_svs$chrom1) != as.vector(all_svs$chrom2) )
all_svs$SV_type[idx] = "Translocation"
# fix incorrect lengths
idx = which(all_svs$SV_type =="BND" & all_svs$strand1 =="+" &  all_svs$strand2 =="-" & all_svs$chrom1==all_svs$chrom2)
all_svs$abs_svLen[idx] = abs(all_svs$start1[idx] -  all_svs$start2[idx])
# more plottable name
idx = which(all_svs$SV_type =="SBND")
all_svs$SV_type[idx] = "Single breakend"
table(all_svs$SV_type)

saveRDS(all_svs, file="rds/illumina_svs.rds")

