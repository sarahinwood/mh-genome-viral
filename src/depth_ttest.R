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

#depth file
st_depth_file <- snakemake@input[["st_depth_file"]]
##scaffold labels
scaffold_id_table <- snakemake@input[["scaffold_id_table"]]

st_depth_names <- c("Scaffold_full_id", "BP", "depth")
st_depth <- fread(st_depth_file, col.names=st_depth_names)
scaffold_table <- fread(scaffold_id_table, header=TRUE)
st_depth_boxpl <- merge(st_depth, scaffold_table, by="Scaffold_full_id", all.x=TRUE)
test_res <- t.test(depth ~ test_label, data=st_depth_boxpl)
chars <- capture.output(print(test_res))
writeLines(chars, con=file(snakemake@output[["ttest_results"]]))

##histogram of depths
viral_contigs <- st_depth_boxpl[st_depth_boxpl$test_label=="Viral"]
pdf(snakemake@output[["viral_depth_hist"]])
viral_contigs[,hist(depth, breaks=1000, xlim=c(0,400))]
dev.off()

hic_scaffolds <- st_depth_boxpl[st_depth_boxpl$test_label=="Hi-C"]
pdf(snakemake@output[["hic_depth_hist"]])
hic_scaffolds[,hist(depth, breaks=1000, xlim=c(0,400))]
dev.off()

wilcox <- wilcox.test(depth ~ test_label, data=st_depth_boxpl)
wilcox_chars <- capture.output(print(wilcox))
writeLines(wilcox_chars, con=file(snakemake@output[["wilcox_res"]]))