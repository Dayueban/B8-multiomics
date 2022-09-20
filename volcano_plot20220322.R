rm(list = ls())
#a <- read.csv("MSMS/MSMSfinal_DEM.csv", header = T, check.names = F, row.names = 1 )
a <- read.csv("Model3 vs Control/1-差异分析/Different_origin.csv", header = T, check.names = F, row.names = 1 )
#a$log2FC <- as.numeric(a$log2FC)
#logFC_cutoff <- with(a, mean(abs(log2FC)) + 2*sd(abs(log2FC)))
logFC_cutoff <- 1  # fold change 大于2倍和fdr小于0.05为显著
a$change <- as.factor(ifelse(a$FDR < 0.05 & 
                               abs(a$log2FC) > logFC_cutoff,
                             ifelse(a$log2FC > logFC_cutoff, 'UP',
                                    'DOWN'),'NOT'))


this_tile <- paste0('Cutoff for log2FC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(a[a$change == 'UP',]),
                    '\nThe number of down gene is ',nrow(a[a$change == 'DOWN',]))
library(ggplot2)
library(ggpubr)
library(ggrepel)
b = data.frame(rownames(a),a)
colnames(b)[1] <- "ID"
metabo_list <- as.character(c("Sphingosine","Inosine","Glutathione","Spermine"))
g <- ggplot(data=b, aes(x=log2FC, y=-log10(FDR), color=change)) + 
  geom_point(alpha = 0.9,aes(size = VIP)) +
  geom_point() +
  theme_bw() +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle(this_tile) + theme(plot.title = element_text(size = 15, hjust = 0.5)) +
  scale_color_manual(values = c("#E16663","black","#90CBD3")) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.5) +
  xlim(-6,8) +
  ylim(0,4.5)

g <- g + theme(axis.title = element_text(size = 18),
               axis.text = element_text(size = 16),
               legend.key.height=unit(1,"cm"),
               legend.key.width=unit(1,"cm"),
               legend.text=element_text(lineheight=0.8,face="bold",size=16),
               legend.title=element_text(size=20))

b1 <- subset(b, grepl(paste(metabo_list, collapse = "|"), rownames(b)))
g1 <- g + geom_label_repel(
  data = b1,
  #data$label <- metabo_list,
  #aes(x=log2FC, y=-log10(Pvalue),label = label),
  aes(label = rownames(b1)),
  stat = "identity",
  position = "identity",
  size = 5,
  box.padding = unit(0.6,'lines'),
  label.padding = unit(0.5,'lines'),
  point.padding = unit(0.4,'lines'),segment.color = "black",show.legend = FALSE
)
print(g1)
ggsave(g1,filename = "Model3 vs Control/1-差异分析/volcano_diff2.png", width = 8, height = 8)


