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
  geom_point(data=subset(depth_vs_gc, hic_plot_label=="non-viral, unclustered"),
             aes(x=GC, y=mean_depth, color=hic_plot_label, alpha=hic_plot_label), size=3)+
  geom_point(data=subset(depth_vs_gc, hic_plot_label == "non-viral, clustered"),
             aes(x=GC, y=mean_depth, color=hic_plot_label, alpha=hic_plot_label), size=3)+
  geom_point(data = subset(depth_vs_gc, hic_plot_label == "probably not viral, clustered"),
             aes(x = GC, y = mean_depth, color = hic_plot_label, alpha=hic_plot_label), size=3)+ 
  geom_point(data = subset(depth_vs_gc, hic_plot_label == 'viral, unclustered'),
             aes(x = GC, y = mean_depth, color = hic_plot_label, alpha=hic_plot_label), size=3)+
  scale_color_manual(name = "Scaffold is:", values = c("blue", "grey50", "gold", "red"))+
  scale_alpha_manual(name = "Scaffold is:", values = c("non-viral, unclustered"=0.15, "non-viral, clustered"=0.3, "viral, unclustered"=0.7, "probably not viral, clustered"=0.7))+
  theme_light()

##plot depth vs gc - without hi-c unclustered scaffolds
ggplot(depth_vs_gc) +
  geom_point(data=subset(depth_vs_gc, hic_plot_label == "non-viral, clustered"),
             aes(x=GC, y=mean_depth, color=hic_plot_label, alpha=hic_plot_label), size=3)+
  geom_point(data = subset(depth_vs_gc, hic_plot_label == "probably not viral, clustered"),
             aes(x = GC, y = mean_depth, color = hic_plot_label, alpha=hic_plot_label), size=3)+ 
  geom_point(data = subset(depth_vs_gc, hic_plot_label == 'viral, unclustered'),
             aes(x = GC, y = mean_depth, color = hic_plot_label, alpha=hic_plot_label), size=3)+
  scale_color_manual(name = "Scaffold is:", values = c("blue", "gold", "red"))+
  scale_alpha_manual(name = "Scaffold is:", values = c("non-viral, clustered"=0.3, "viral, unclustered"=0.7, "probably not viral, clustered"=0.7))+
  theme_light()

###probably need to set axis limits if they differ from plot above and save with same dimensions (7x10 for above??)
##coverage box plot
ggplot(depth_vs_gc, aes(x=hic_plot_label, y=mean_depth, colour=hic_plot_label))+
  scale_colour_manual(values=c("blue","grey50", "gold", "red"))+
  geom_boxplot(outlier.alpha = 0.2)+
  theme_light()+
  theme(legend.position = "none")

##gc content box plot
ggplot(depth_vs_gc, aes(x=hic_plot_label, y=GC, colour=hic_plot_label))+
  scale_colour_manual(values=c("blue","grey50", "gold", "red"))+
  geom_boxplot(outlier.alpha = 0.2)+
  coord_flip()+
  theme_light()+
  theme(legend.position = "none")
