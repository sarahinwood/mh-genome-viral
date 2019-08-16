#!/usr/bin/env Rscript

#######
# LOG #
#######

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

#############
# LIBRARIES #
#############

library(data.table)

###########
# GLOBALS #
###########

depth <- snakemake@input[["depth"]]

########
# MAIN #
########

depth_table <- fread(depth)
scaffold_ids <- unique(depth_table$`#CHROM`)

MEAN_SCAFFOLD_DEPTH <- function(x, depth){
  scaffold_depth<-depth[depth$`#CHROM`==x, mean(`output/samtools/sorted.bam`)]
  return(data.table(scaffold_id=x, mean_depth=scaffold_depth))
}

scaffold_depths <- lapply(scaffold_ids, MEAN_SCAFFOLD_DEPTH, depth=depth_table)
scaff_depths_table <- rbindlist(scaffold_depths)
fwrite(scaff_depths_table, snakemake@output[["mean_depth_table"]])

#write log
sessionInfo()
