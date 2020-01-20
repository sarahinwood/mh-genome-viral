library(data.table)
library(ggplot2)

depth_vs_gc <- fread("output/depth_vs_gc/depth_vs_gc.csv")

##depth vs gc  - hi-c cluster status
ggplot(depth_vs_gc) +
  geom_point(data=subset(depth_vs_gc, cluster_status=="unclustered"),
             aes(x=GC, y=mean_depth, color=cluster_status), size=3)+
  geom_point(data=subset(depth_vs_gc, cluster_status == "clustered"),
             aes(x=GC, y=mean_depth, color=cluster_status), size=3)

##plot depth vs gc full graph - all scaffolds
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

##plot depth vs gc - without hi-c unclustered scaffolds
ggplot(depth_vs_gc) +
  geom_point(data=subset(depth_vs_gc, cluster_status == "clustered"),
             aes(x=GC, y=mean_depth, color=cluster_status, alpha=cluster_status), size=3)+
  geom_point(data = subset(depth_vs_gc, label == "probably not viral, clustered"),
             aes(x = GC, y = mean_depth, color = label, alpha=label), size=3)+ 
  geom_point(data = subset(depth_vs_gc, label == 'viral'),
             aes(x = GC, y = mean_depth, color = label, alpha=label), size=3)+
  scale_color_manual(values = c("blue", "gold", "red"))+
  scale_alpha_manual(values = c("clustered"=0.3, "viral"=0.7, "probably not viral, clustered"=0.7))+
  theme_light()

###probably need to set axis limits if they differ from plot above and save with same dimensions (7x10 for above??)
##coverage box plot
ggplot(depth_vs_gc, aes(x=cluster_status, y=mean_depth, colour=plot_label))+
  geom_boxplot(outlier.alpha = 0.2)+
  scale_colour_manual(values=c("blue","grey50", "gold", "red"))+
  theme_light()

##gc content box plot
ggplot(depth_vs_gc, aes(x=cluster_status, y=GC, colour=plot_label))+
  geom_boxplot(outlier.alpha = 0.2)+
  scale_colour_manual(values=c("blue","grey50", "gold", "red"))+
  coord_flip()+
  theme_light()
