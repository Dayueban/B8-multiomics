# More rigorous screening of metabolites
rm(list = ls())
library(pheatmap)
library(reshape2)
a <- read.csv("Model1 vs Control/1-差异分析/Different_origin.csv", header = T,
              row.names = 1, check.names = F)
# Sort data from smallest to largest with column "FC" in a 
a <- a[order(a$FC, decreasing = T),]
b <- read.csv("Model1 vs Control/1-差异分析/Different_metabo_intensity.csv",
              header = T, row.names = 1, check.names = F)
# select rows with multiple conditions
keep_metabo <- rownames(a[abs(a$log2FC) > 1 & a$FDR <= 0.05 & a$VIP > 1, ])
# match the dataframe b with specified metabolites
b2 <- b[keep_metabo,]

c <- data.frame(t(b2))
names(c) <- keep_metabo
c2 <- cbind(rownames(c),c)
names(c2)[1] <- "ID"

c_rs <- melt(c2,id.vars = "ID",variable.name = "Metabolites", value.name = "Relative_abundance")
# adding a new cloumn as the group information
c_rs$group <- "NA"
c_rs$ID <- as.character(c_rs$ID)
c_rs$group[grep("^C",c_rs$ID)] <- "Control"
c_rs$group[grep("^M2",c_rs$ID)] <- "Model1"
# draw barplot a
library(ggplot2)
#par(mar=c(2,2,2,2))
p <- ggplot(data=c_rs,aes(x=Metabolites,y=Relative_abundance,fill=group))+geom_boxplot(notch = FALSE)
p <- p + theme_classic()
#p <- p + theme(axis.text.x=element_text(angle=60,color=c(rep("red",51),rep("blue",14)),
#                                        vjust = 0.95,hjust = 0.95,size=10, face="bold"),
#               axis.text.y=element_text(size = 10, colour = "black"),
#               axis.title.y=element_text(size = 14,face = "bold"))
p <- p + theme(axis.text.x=element_text(angle=60,color="black",
                                        vjust = 0.95,hjust = 0.95,size=10, face="bold"),
               axis.text.y=element_text(size = 10, colour = "black"),
               axis.title.y=element_text(size = 14,face = "bold"))
p <- p + theme(plot.background=element_blank(),panel.background=element_blank())
p <- p + guides(fill=guide_legend(title=NULL))
p <- p + theme(legend.position = "top")
p <- p + theme(plot.title = element_text(size = 20))
p <- p + xlab("")
p <- p + labs(y = "Relative abundance (log10)")
p <- p + scale_fill_manual(values=c("#E16663", "#7292AA"),
                           labels = c("Control", "Model1"))
p <- p + theme(axis.line.x=element_line(colour = "black"),axis.line.y=element_line(colour = "black"),
               legend.key.height=unit(0.5,"cm"),
               legend.key.width=unit(0.5,"cm"),
               legend.text=element_text(lineheight=0.8,face="bold",size=16),
               legend.title=element_text(size=20))
#set the margin of plot
p <- p + theme(plot.margin = unit(c(1,0.5,1,0.5),"cm"))
p <- p + scale_y_continuous(trans = "log10")
#p <- p + scale_y_continuous(limits=c(miny*1.1,maxy*1.1))
##p <- p + geom_jitter() ##添加数据点
ggsave("Model1 vs Control/1-差异分析/Metabolites_diff_20220324.tiff",width = 14,height = 8)



