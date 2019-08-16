library(data.table)
library(ggplot2)

full_viral_scaffold_list <- fread("data/viral_scaffold_ids.txt", header=FALSE)
##80 and 8624 have non-viral hits on them
dedup_viral_list <- unique(full_viral_scaffold_list)
prob_not_viral <- list("Scaffold80", "Scaffold8624")
viral_scaffold_list <- subset(dedup_viral_list, !(V1 %in% prob_not_viral))

mh_gc <- fread("output/bb_stats/gc.txt")
mh_gc_hist <- fread("output/bb_stats/gc_hist.out")

##plot GC% in histogram
hist(mh_gc$GC, breaks=100, xlab="GC Content", main="M.hyp Genome Scaffolds GC Content")

##subset GC content into viral and non-viral

viral_GC <- subset(mh_gc, (`#Name` %in% viral_scaffold_list$V1))
non_viral_GC <- subset(mh_gc, !(`#Name` %in% viral_scaffold_list$V1))
##mean viral GC is 0.34, non viral 0.29 - not drastically different
##Lc mean = 0.29, Lb = 0.36, Lh = 0.29

##read in mean scaffold depths for scaffolds from snakemake script
mean_scaff_depths <- fread('output/samtools/mean_depth_table.csv')
##plot mean depths
hist(mean_scaff_depths$mean_depth, breaks=1000, xlim = c(0,500))
##make table to depth vs GC
depth_vs_gc <- merge(mean_scaff_depths, mh_gc, by.x="scaffold_id", by.y="#Name", all=TRUE)
##add column for suspected viral contig?
depth_vs_gc$label <- ifelse(depth_vs_gc$scaffold_id %in% viral_scaffold_list$V1, 'viral', 'non-viral')

viral_depth_gc <- subset(depth_vs_gc, scaffold_id %in% viral_scaffold_list$V1)
##remove 2 scaffolds probably not viral
hist(viral_depth_gc$mean_depth, breaks=100, xlim = c(0,500))

non_viral_depth_gc <- subset(depth_vs_gc, !(scaffold_id %in% viral_scaffold_list$V1))
hist(non_viral_depth_gc$mean_depth, breaks=1000, xlim = c(0,500))

colours=

ggplot(depth_vs_gc, aes(x=GC, y=mean_depth, colour=label, alpha=label))+
  geom_point()+xlab("GC content")+ylab("mean scaffold depth")+
  scale_alpha_manual(guide='none', values = list("non-viral"=0.1, "viral"=1))

ggplot(depth_vs_gc) +
  geom_point(aes(x = GC, y = mean_depth, color = label, alpha=label)) +
  geom_point(data = subset(depth_vs_gc, label == 'viral'),
             aes(x = GC, y = mean_depth, color = label, alpha=label))+
  scale_color_manual(values = c("grey40", "red"))+
  scale_alpha_manual(guide='none', values = list("non-viral"=0.2, "viral"=0.8))+
  theme_light()
