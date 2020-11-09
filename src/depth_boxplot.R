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
library(ggplot2)
library(viridis)

###########
# GLOBALS #
###########

st_depth_file <- snakemake@input[["st_depth_file"]]
scaffold_id_table <- snakemake@input[["scaffold_id_table"]]

########
# MAIN #
########

##read in depth file
st_depth_names <- c("Scaffold_full_id", "BP", "depth")
st_depth <- fread(st_depth_file, col.names=st_depth_names)
##merge with table of scaffold ids for plotting
scaffold_table <- fread(scaffold_id_table, header=TRUE)
st_depth_boxpl <- merge(st_depth, scaffold_table, by="Scaffold_full_id", all.x=TRUE)
#order chr numerically
st_depth_boxpl$Scaffold_id <- factor(st_depth_boxpl$Scaffold_id,
    levels=c("PGA_scaffold0", "PGA_scaffold1", "PGA_scaffold2",
      "PGA_scaffold3", "PGA_scaffold4", "PGA_scaffold5",
      "PGA_scaffold6", "PGA_scaffold7", "PGA_scaffold8",
      "PGA_scaffold9", "PGA_scaffold10", "PGA_scaffold11",
      "Scaffold1791", "Scaffold5674", "Scaffold6852",
      "Scaffold13392", "Scaffold15344", "Scaffold15995",
      "Scaffold16687", "Scaffold18966", "Scaffold27879",
      "Scaffold28498", "Scaffold29164", "Scaffold30641", "Scaffold32600", 
      ##non-LbFV containing contigs
      "Scaffold3332", "Scaffold8243", "Scaffold9703",
      "Scaffold32814", "Scaffold43528", "Scaffold45974"))
##order legend labels
st_depth_boxpl$Scaffold_label <- factor(st_depth_boxpl$Scaffold_label,
	levels=c("Hi-C Scaffold", "Viral Contig, LbFV hit", "Viral Contig"))

##boxplot, no outliers, zoomed in
pdf(snakemake@output[["boxplot_no_outliers"]])

ggplot(st_depth_boxpl, aes(x=Scaffold_id, y=depth, colour=Scaffold_label))+
  ##boxplot, make outlier points not plot - still scales Y accounting for them though
  geom_boxplot(outlier.shape=NA)+
  ##white background
  theme_light()+
  ##turn scaffold labels 90 degrees, align with middle of tick up against it
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+
  ylab("Depth")+
  ##plot means
  stat_summary(fun.y=mean, geom="point", colour="grey35")+
  ##colour-blind friendly palette
  scale_colour_viridis(discrete=TRUE)+
  ##zoom in on depth 0-3000 without disrupting plotting
  coord_cartesian(ylim = c(0, 350))
dev.off()

##full boxplot with outliers
pdf(snakemake@output[["boxplot"]])
ggplot(st_depth_boxpl, aes(x=Scaffold_id, y=depth, colour=Scaffold_label))+
  ##boxplot, make outlier points not plot - still scales Y accounting for them though
  geom_boxplot(outlier.alpha=0.1)+
  ##white background
  theme_light()+
  ##turn scaffold labels 90 degrees, align with middle of tick up against it
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+
  ylab("Depth")+
  ##plot means
  stat_summary(fun.y=mean, geom="point", colour="grey35")+
  ##colour-blind friendly palette
  scale_colour_viridis(discrete=TRUE)
dev.off()


#write log
sessionInfo()

##automatically compute values for coord_cartesian - but how to make this work with multiple boxplots?
# compute lower and upper whiskers
#ylim1 = boxplot.stats(df$y)$stats[c(1, 5)]
# scale y limits based on ylim1
#p1 = p0 + coord_cartesian(ylim = ylim1*1.05)
