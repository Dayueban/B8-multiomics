#' Differential analysis of metabolites upon different treatments
#' 2021/09/10
rm(list=ls())
library(ade4)
library(RColorBrewer)
library(vegan)
library(randomForest)
library(dplyr)

Taxo_G <- read.csv("Model1 vs Control/1-差异分析/Model1toControl.csv",header=T,row.names = 1, check.names = F)
group_list <- Taxo_G$group
# table(group_list)
# Female   Male 
# 6      4
Taxo_G$group <- NULL
Taxo_G <- as.data.frame(t(Taxo_G))
#' Transform into relative abundance so that each row sums to 100
Taxo_G2 <- sweep(Taxo_G, 2, colSums(Taxo_G), '/')*100
#genus_total_data <- log10(genus_total_data + 1 )

data1 <- Taxo_G2
nr <- nrow(data1)
#which(genus_total_class == "female")
#nc <- ncol(mergedMEs_tr)
#预生成2个长度与输入文件行数相同的全为0的向量，将用于存储p value和差异倍数（log2FC）
Pvalue<-c(rep(0,nr))
log2_FC<-c(rep(0,nr))
FC <- c(rep(0,nr))
# 1~26列是处理组1,27~71列是处理组2；
#将使用循环对每一行进行t检验
#如果某一行两组的标准差都等于0，将无法进行t检验，所以这些行使用NA替代
#每一行计算得到p value和log2FC将被加入原文件的后两列；
#计算log2FC时，每组均值加0.001，是为了防止分母为0导致bug；
for (i in 1:nr){
  if(sd(data1[i,1:10])==0 && sd(data1[i,11:ncol(data1)])==0){
    Pvalue <- "NA"
    log2_FC<- "NA"
    FC <- "NA"
  }else{
    Pvalue[i] <- wilcox.test(as.numeric(data1[i,1:10]),as.numeric(data1[i,11:ncol(data1)]))$p.value
    FC[i] <- (mean(as.numeric(data1[i,i,1:10]))+0.001)/ (mean(as.numeric(data1[i,11:ncol(data1)]))+0.001)
    log2_FC[i] <- log2((mean(as.numeric(data1[i,1:10]))+0.001)/ (mean(as.numeric(data1[i,11:ncol(data1)]))+0.001))
  }
}
Pvalue<-as.numeric(Pvalue)  #在操作的过程中发现秩和检验的包输出来不是列表的形式，这里变成列表
fdr.w<- p.adjust(Pvalue,method="fdr",length(Pvalue))    #p.adjust就是计算FDR的包，这个可要记得了
Taxo_G3 <- cbind(data1,Pvalue,fdr.w,FC,log2_FC)
#output data
write.csv(Taxo_G3,file="Model1 vs Control/1-差异分析/metabolites_differ_Control2Model1.csv")



# without computing sd value
for (i in 1:nr){
  Pvalue[i] <- wilcox.test(as.numeric(data1[i,1:5]),as.numeric(data1[i,6:ncol(data1)]))$p.value
  FC[i] <- (mean(as.numeric(data1[i,i,1:5]))+0.001)/ (mean(as.numeric(data1[i,6:ncol(data1)]))+0.001)
  log2_FC[i] <- log2((mean(as.numeric(data1[i,1:5]))+0.001)/ (mean(as.numeric(data1[i,6:ncol(data1)]))+0.001))
}

Pvalue<-as.numeric(Pvalue)  #在操作的过程中发现秩和检验的包输出来不是列表的形式，这里变成列表
fdr.w<- p.adjust(Pvalue,method="fdr",length(Pvalue))    #p.adjust就是计算FDR的包，这个可要记得了
Taxo_G3 <- cbind(data1,Pvalue,fdr.w,FC,log2_FC)
#output data
write.csv(Taxo_G3,file="taxonomy/Species/species_Ct_differ.csv")


#*******************second part--drawing barplot***************
rm(list=ls())
data2 <- read.csv("taxonomy/Species/species_differ_plot.csv",header = T, check.names = F)
library(reshape2)
## reshape dataframe, from long sheet to wild sheet
data2_rs <- melt(data2,id.vars = "ID",variable.name = "Taxa", value.name = "Relative_abundance")
# add a new cloumn as the group information
data2_rs$group <- "NA"
data2_rs$ID <- as.character(data2_rs$ID)
data2_rs$group[grep("^M",data2_rs$ID)] <- "Mild"
data2_rs$group[grep("^H",data2_rs$ID)] <- "Health"

# draw barplot a
library(ggplot2)
#par(mar=c(2,2,2,2))
p <- ggplot(data=data2_rs,aes(x=Taxa,y=Relative_abundance,fill=group))+geom_boxplot(notch = FALSE)
p <- p + theme_classic()
p <- p + theme(axis.text.x=element_text(angle=60,color="black",vjust = 0.95,hjust = 0.95,size=14, face="bold"),
               axis.text.y=element_text(size = 16, colour = "black"),
               axis.title.y=element_text(size = 18,face = "bold"))
p <- p + theme(plot.background=element_blank(),panel.background=element_blank())
p <- p + guides(fill=guide_legend(title=NULL))
p <- p + theme(legend.position = "top")
p <- p + theme(plot.title = element_text(size = 20))
p <- p + xlab("")
p <- p + labs(y = "Relative abundance (log10)")
p <- p + scale_fill_manual(values=c("#80B1D3", "#FB8072"),
                           labels = c("Health", "Mild"))
p <- p + theme(axis.line.x=element_line(colour = "black"),axis.line.y=element_line(colour = "black"),
               legend.key.height=unit(1,"cm"),
               legend.key.width=unit(1,"cm"),
               legend.text=element_text(lineheight=0.8,face="bold",size=16),
               legend.title=element_text(size=20))
#set the margin of plot
p <- p + theme(plot.margin = unit(c(1,0.5,1,0.5),"cm"))
p <- p + scale_y_continuous(trans = "log10")
#p <- p + scale_y_continuous(limits=c(miny*1.1,maxy*1.1))
##p <- p + geom_jitter() ##添加数据点
ggsave("taxonomy/Species/Species_diff20210915.tiff",width = 6,height = 8)
