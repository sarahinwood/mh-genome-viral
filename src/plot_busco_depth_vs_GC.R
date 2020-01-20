library(data.table)
library(ggplot2)

depth_vs_gc <- fread("output/depth_vs_gc/depth_vs_gc.csv")

##plot depth vs gc, all scaffolds
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

##plot depth vs gc - without 'non-busco, non-viral'
ggplot(depth_vs_gc) +
  geom_point(data=subset(depth_vs_gc, busco == "BUSCO"),
             aes(x=GC, y=mean_depth, color=busco, alpha=busco), size=3)+
  geom_point(data = subset(depth_vs_gc, label == "probably not viral, clustered"),
             aes(x = GC, y = mean_depth, color=label, alpha=label), size=3)+ 
  geom_point(data = subset(depth_vs_gc, label == 'viral'),
             aes(x = GC, y = mean_depth, color=label, alpha=label), size=3)+
  scale_color_manual(values = c("blue", "gold", "red"))+
  scale_alpha_manual(values = c("BUSCO"=0.3, "viral"=0.7, "probably not viral, clustered"=0.7))+
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
