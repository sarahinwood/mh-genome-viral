library(data.table)

viral_scaffold_list <- fread("data/viral_scaffold_ids.txt", header=FALSE)
mh_gc <- fread("output/bb_stats/gc.txt")
mh_gc_hist <- fread("output/bb_stats/gc_hist.out")

##plot GC% in histogram
hist(mh_gc$GC, breaks=100, xlab="GC Content", main="M.hyp Genome Scaffolds GC Content")

##subset GC content into viral and non-viral
viral_GC <- subset(mh_gc, `#Name` %in% viral_scaffold_list$V1)
##80 and 8624 have non-viral hits on them
prob_not_viral <- list("Scaffold80", "Scaffold8624")
viral_minus_prob_nv_GC <- subset(viral_GC, !(`#Name` %in% prob_not_viral))
non_viral_GC <- subset(mh_gc, !(`#Name` %in% viral_scaffold_list$V1))
##mean viral GC is 0.34, non viral 0.29 - not drastically different
##Lc mean = 0.29, Lb = 0.36, Lh = 0.29

depth <- fread("output/samtools/depth.out")
##make list of scaffold IDs to feed into below function
scaffold_ids <- unique(depth$`#CHROM`)

#want mean of output/samtools.sorted.bam column when #CHROM = each Scaffold
##calculate mean depth for a scaffold e.g. scaffold 5
depth[depth$`CHROM`=="Scaffold80",mean(depth$`output/samtools/sorted.bam`)]
##adapt above into a function to loop through each scaffold as below

EXTRACT_GO_TERMS <- function(x, trinotate_report){
  my_terms<-trinotate_report[`#gene_id`==x,unique(unlist(strsplit(gene_ontology_pfam, "`")))]
  my_accessions<-unique(gsub("\\^.*", "", my_terms))
  my_accessions<-my_accessions[!is.na(my_accessions)]
  return(data.table(gene_id=x, accessions=my_accessions))
}

x <- depth[depth$`#CHROM`==x,]
x_mean <- mean(x$`output/samtools/sorted.bam`)
return(data.table)

##function to calculate mean depth for all scaffolds
MEAN_SCAFFOLD_DEPTH <- function(x, depth){
  scaffold_depth<-depth[depth$`#CHROM`==x, mean(`output/samtools/sorted.bam`)]
  return(data.table(scaffold_id=x, mean_depth=scaffold_depth))
}
##run above function to calculate mean depth for all scaffolds in list
##gives all scaffolds at 190.14???
scaffold_depths <- lapply(scaffold_ids, MEAN_SCAFFOLD_DEPTH, depth=depth)
scaff_depths_table <- rbindlist(scaffold_depths)

hist(depth$`output/samtools/sorted.bam`, breaks=1000, xlim = c(0,500))
