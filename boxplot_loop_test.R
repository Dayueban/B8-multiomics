# More rigorous screening of metabolites
rm(list = ls())
library(pheatmap)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(agricolae)
library(dplyr)
library(multcompView)
library(cowplot)
a <- read.csv("Different_overall_tr.csv", header = T,
              row.names = 1, check.names = F)

metabo_list <- as.character(c("Sphingosine", "Inosine", "Glutathione", "Spermine",
                              "Hippuric_acid", "Indole_3_carboxylic_acid",
                              "Gluconic_acid", "D_Mannose"))
a1 <- a[,names(a) %in% metabo_list]
a1 <- log10(a1)
a1$group <- c(rep("Control",10),rep("Model1",10),rep("Model2",10),rep("Model3",10))


Start = 1
Stop = 8
gvec <- vector("list",length = length(Start:Stop))

colors_index <- c("#E16663", "#7292AA","#90CBD3","#EAA58E")
names(colors_index) <- c("Control", "Model1", "Model2", "Model3")
colors_index <- as.data.frame(colors_index)
colors_index$group <- rownames(colors_index)

for(i in Start:Stop){
  m=names(a1)[i]
  index=a1[,c(m,"group")]
  model = aov(index[[m]] ~ group, data=index)
  # LSD test for stat label
  out = LSD.test(model,"group", p.adj="none") # alternative fdr
  
  stat = out$groups
  stat$colors = colors_index
  index$stat=stat[as.character(index$group),]$groups
  index$colors=stat[as.character(index$group),]$colors
  max=max(index[,c(m)])
  min=min(index[,c(m)])
  x = index[,c("group",m)]
  y = x %>% group_by(group) %>% summarise_(Max=paste('max(',m,')',sep=""))
  y=as.data.frame(y)
  rownames(y)=y$group
  index$y=y[as.character(index$group),]$Max + (max-min)*0.05
  
  batch_neu <- ggplot(data = index, 
                      aes_string(x="group",y=m, fill="group"))+
    scale_fill_manual(values = c("#E16663", "#7292AA","#90CBD3","#EAA58E"))+
    geom_jitter(size = 1.2)+
    geom_boxplot(size=1.5, alpha=.6)+
    xlab("")+
    ylab("metabolites abundance (log10)")+
    ggtitle(m)+
    theme_classic()+
    theme(legend.position = "none")+
    geom_text(data=index, aes(x=group, y=y, color=group, label= stat),size=4) +
    scale_y_continuous(labels=function(y) format(y,scientific=FALSE))+
    theme(axis.title.y = element_text(size = 8,face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 8),
          axis.text.x = element_text(size = 8, color = "black", angle = 60,
                                     hjust = 0.95, vjust = 0.95),
          axis.text.y = element_text(size = 8, color = "black"))
  ggsave(batch_neu, file = paste0("plot_",m,".png"),
         units="in", width=4, height=3, dpi=300)
  gvec[[i-Start+1]] <- batch_neu
}

pp <- plot_grid(plotlist = gvec, nrow = 2, labels = c('A', 'B', 'C', 'D', 'E','F','G')) # R package cowplot
ggsave(pp, file = "merge_metabolite2.png", height = 6, width = 8, dpi = 300)

