library(data.table)
library(ggplot2)

depth_vs_gc <- fread("output/depth_vs_gc/depth_vs_gc.csv")

##plot depth vs gc, all scaffolds
##plot depth vs gc full graph - all scaffolds
ggplot(depth_vs_gc) +
  geom_point(data=subset(depth_vs_gc, busco_plot_label=="non-viral, no BUSCO"),
             aes(x=GC, y=mean_depth, color=busco_plot_label, alpha=busco_plot_label), size=3)+
  geom_point(data=subset(depth_vs_gc, busco_plot_label == "non-viral, BUSCO-containing"),
             aes(x=GC, y=mean_depth, color=busco_plot_label, alpha=busco_plot_label), size=3)+
  geom_point(data = subset(depth_vs_gc, busco_plot_label == "probably not viral, BUSCO-containing"),
             aes(x = GC, y = mean_depth, color = busco_plot_label, alpha=busco_plot_label), size=3)+ 
  geom_point(data = subset(depth_vs_gc, busco_plot_label == 'viral, no BUSCO'),
             aes(x = GC, y = mean_depth, color = busco_plot_label, alpha=busco_plot_label), size=3)+
  scale_color_manual(name = "Scaffold is:", values = c("blue", "grey50", "gold", "red"))+
  scale_alpha_manual(name = "Scaffold is:", values = c("non-viral, no BUSCO"=0.15, "non-viral, BUSCO-containing"=0.3, "viral, no BUSCO"=0.7, "probably not viral, BUSCO-containing"=0.7))+
  theme_light()

##plot depth vs gc - without 'non-busco, non-viral'
ggplot(depth_vs_gc) +
  geom_point(data=subset(depth_vs_gc, busco_plot_label == "non-viral, BUSCO-containing"),
             aes(x=GC, y=mean_depth, color=busco_plot_label, alpha=busco_plot_label), size=3)+
  geom_point(data = subset(depth_vs_gc, busco_plot_label == "probably not viral, BUSCO-containing"),
             aes(x = GC, y = mean_depth, color = busco_plot_label, alpha=busco_plot_label), size=3)+ 
  geom_point(data = subset(depth_vs_gc, busco_plot_label == 'viral, no BUSCO'),
             aes(x = GC, y = mean_depth, color = busco_plot_label, alpha=busco_plot_label), size=3)+
  scale_color_manual(name = "Scaffold is:", values = c("blue", "gold", "red"))+
  scale_alpha_manual(name = "Scaffold is:", values = c("non-viral, BUSCO-containing"=0.3, "viral, no BUSCO"=0.7, "probably not viral, BUSCO-containing"=0.7))+
  theme_light()

###probably need to set axis limits if they differ from plot above and save with same dimensions (7x10 for above??)
##coverage box plot
ggplot(depth_vs_gc, aes(x=busco_plot_label, y=mean_depth, colour=busco_plot_label))+
  scale_colour_manual(values=c("blue","grey50", "gold", "red"))+
  geom_boxplot(outlier.alpha = 0.2)+
  theme_light()+
  theme(legend.position = "none")
  
##gc content box plot
ggplot(depth_vs_gc, aes(x=busco_plot_label, y=GC, colour=busco_plot_label))+
  scale_colour_manual(values=c("blue","grey50", "gold", "red"))+
  geom_boxplot(outlier.alpha = 0.2)+
  coord_flip()+
  theme_light()+
  theme(legend.position = "none")
