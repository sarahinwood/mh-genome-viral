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
mean_scaff_depths <- fread('output/samtools/smalt_mean_depth_table.csv')
##plot mean depths
hist(mean_scaff_depths$mean_depth, breaks=1000, xlim = c(0,500))
##make table to depth vs GC
depth_vs_gc <- merge(mean_scaff_depths, mh_gc, by.x="scaffold_id", by.y="#Name", all=TRUE)
##add column for suspected viral contig?
depth_vs_gc$label <- ifelse(depth_vs_gc$scaffold_id %in% viral_scaffold_list$V1, 'viral',
                            ifelse(depth_vs_gc$scaffold_id %in% prob_not_viral, 'probably not viral, clustered', 'non-viral'))

viral_depth_gc <- subset(depth_vs_gc, scaffold_id %in% viral_scaffold_list$V1)
hist(viral_depth_gc$mean_depth, breaks=100, xlim = c(0,500))

non_viral_depth_gc <- subset(depth_vs_gc, !(scaffold_id %in% viral_scaffold_list$V1))
hist(non_viral_depth_gc$mean_depth, breaks=1000, xlim = c(0,500))

high_depth_nudivirus <- list("Scaffold28315", "Scaffold3939")
viral_minus_nudivirus <- subset(viral_depth_gc, !(scaffold_id %in% high_depth_nudivirus))

##read in busco data
busco <- fread("output/busco/run_mh_genome/full_table_mh_genome.tsv", header=TRUE, fill=TRUE, skip = 4)
busco_scaffolds <- busco$Contig
##create column labelling whether contigs contain BUSCO gene
depth_vs_gc$busco <- ifelse(depth_vs_gc$scaffold_id %in% busco_scaffolds, 'BUSCO', 'NA')
depth_vs_gc$plot_label <- paste(depth_vs_gc$label,depth_vs_gc$busco)

##playing around for different group means - not for plot
b <- subset(depth_vs_gc, !label == "non-viral")
c <- subset(a, !(scaffold_id %in% high_depth_nudivirus))
mean(c$mean_depth)
mean(c$GC)

##plot depth vs gc
ggplot(depth_vs_gc) +
  geom_point(data=subset(depth_vs_gc, label=="non-viral"),
             aes(x=GC, y=mean_depth, colour=label, alpha=label), size=3)+
  geom_point(data=subset(depth_vs_gc, busco == "BUSCO"),
             aes(x=GC, y=mean_depth, color=busco, alpha=busco), size=3)+
  geom_point(data = subset(depth_vs_gc, label == "probably not viral, clustered"),
             aes(x = GC, y = mean_depth, color=label, alpha=label), size=3)+ 
  geom_point(data = subset(depth_vs_gc, label == 'viral'),
             aes(x = GC, y = mean_depth, color=label, alpha=label), size=3)+
  scale_color_manual(values = c("blue","grey50", "gold", "red"))+
  scale_alpha_manual(values = c("non-viral"=0.15, "BUSCO"=0.3, "viral"=0.7, "probably not viral, clustered"=0.7))+
  theme_light()


###probably need to set axis limits if they differ from plot above and save with same dimensions (7x10 for above??)
##coverage box plot
ggplot(depth_vs_gc, aes(x=plot_label, y=mean_depth, colour=plot_label))+
  geom_boxplot(outlier.alpha = 0.2)+
  scale_colour_manual(values=c("blue","grey50", "gold", "red"))+
  theme_light()

##gc content box plot
ggplot(depth_vs_gc, aes(x=plot_label, y=GC, colour=plot_label))+
  geom_boxplot(outlier.alpha = 0.2)+
  scale_colour_manual(values=c("blue","grey50", "gold", "red"))+
  coord_flip()+
  theme_light()

##make table of scaffold id to hi-c cluster id
hi_c_scaffold_clusters <- fread("data/hi-c_genome/clusters.by_name_simple.txt", skip = c(16), fill=TRUE)
scaffold_to_cluster <- hi_c_scaffold_clusters[,c(2,4)]
setnames(scaffold_to_cluster, old=c("V2", "V4"), new=c("scaffold_id", "cluster_id"))
scaffold_to_cluster$cluster_id <- tstrsplit(scaffold_to_cluster$cluster_id, "\\(", keep=c(2))
scaffold_to_cluster$cluster_id <- tstrsplit(scaffold_to_cluster$cluster_id, "\\)", keep=c(1))
##merge with depth vs gc table
depth_vs_gc<- merge(depth_vs_gc, scaffold_to_cluster, by="scaffold_id", all.x=TRUE)
depth_vs_gc$cluster_status <- ifelse(!is.na(depth_vs_gc$cluster_id), 'clustered', 'unclustered')
fwrite(depth_vs_gc, "output/depth_vs_gc/depth_vs_gc.csv")

##plot depth vs gc
ggplot(depth_vs_gc) +
  geom_point(data=subset(depth_vs_gc, cluster_status=="unclustered"),
             aes(x=GC, y=mean_depth, color=cluster_status, alpha=cluster_status), size=3)+
  geom_point(data=subset(depth_vs_gc, cluster_status == "clustered"),
             aes(x=GC, y=mean_depth, color=cluster_status, alpha=cluster_status), size=3)+
  geom_point(data = subset(depth_vs_gc, label == "probably not viral, clustered"),
             aes(x = GC, y = mean_depth, color = label, alpha=label), size=3)+ 
  geom_point(data = subset(depth_vs_gc, label == 'viral'),
             aes(x = GC, y = mean_depth, color = label, alpha=label), size=3)+
  scale_color_manual(values = c("blue", "gold", "grey50", "red"))+
  scale_alpha_manual(values = c("unclustered"=0.15, "clustered"=0.3, "viral"=0.7, "probably not viral, clustered"=0.7))+
  theme_light()



