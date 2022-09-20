# associations between different metabolites and genes
rm(list = ls())
library(pheatmap)
library(ggplot2)
#' 第一步，筛选差异的代谢物
a <- read.csv("Model3 vs Control/1-差异分析/Different_origin.csv", header = T,
              row.names = 1, check.names = F)
b <- read.csv("Model3 vs Control/1-差异分析/Different_metabo_intensity.csv",
              header = T, row.names = 1, check.names = F)
# select rows with multiple conditions
keep_metabo <- rownames(a[abs(a$log2FC) > 1 & a$FDR <= 0.05 & a$VIP > 1, ])
# match the dataframe b with specified metabolites
b2 <- b[keep_metabo,]

#' 第二步，筛选差异的基因

temp_dir <- "../../../HT20210918140009 转录组 分析/result/mRNA/2_DEG/DESeq/"
temp_dir2 <- "../../../conjoint_analysis/"

c <- read.csv(paste(temp_dir,"C_vs_M3/C_vs_M3.DESeq.csv",sep = ""), header = T, check.names = F, row.names = 1 )
# b1 <- a[,c(1:9,16:17)] # 24h vs control
# b1 <- a[,c(1:6,10:12,18:19)] # 48h vs control
# b1 <- a[,c(1:6,13:15,20:21)] # 72h vs control
d1 <- c[c$log2FoldChange != Inf ,]
d1 <- d1[d1$log2FoldChange != -Inf, ]
d1 <- d1[d1$Name != "-", ]
names(d1)[5] <- "log2FC"
names(d1)[7] <- "FDR"
#b1$FDR <- p.adjust(b1$`pValue(M1_C)`, method = "fdr", nrow(b1))
#a$log2FC <- as.numeric(a$log2FC)
#logFC_cutoff <- with(a, mean(abs(log2FC)) + 2*sd(abs(log2FC)))
logFC_cutoff <- 0.58 # fold change 大于1.5倍和fdr小于0.05为显著
d1$change <- as.factor(ifelse(d1$FDR < 0.05 & 
                                abs(d1$log2FC) > logFC_cutoff,
                              ifelse(d1$log2FC > logFC_cutoff, 'UP',
                                     'DOWN'),'NOT'))
d1 <- d1[d1$change %in% c("UP", "DOWN"), ]

# 以Apoptosis通路的基因为例
gene_list <- read.delim("clipboard", check.names = F) #从gsea_results.csv表中core_enrichemnt列选择
gene_list <- names(gene_list)
gene_list <- strsplit(gene_list, split = "/")[[1]] #进行分割
keep_genes_df <- d1[d1$Entrez_geneID %in% gene_list, ]
# 将Entrez_ID转换为基因名称

#gene_list <- as.character(c("GRPEL1","TRIM33","LAMP1","SCARB2","MMRN2",
#                              "MGME1","NCL","ATG16L1","ISG15","ERMP1","SQSTM1",
#                              "PAWR","GBP1","GBP2"))
#d2 <- subset(d1, Name %in% gene_list)

# 加载基因表达量数据
f <- read.csv(paste(temp_dir,"Expression.csv",sep = ""), header = T, check.names = F, row.names = 1 )
# 把Model3和control组的提取出来
f_M3C <- f[, c(1:3, 10:12)]
# 再把b1的行名输入f_M3C数据框中得到差异的基因及其表达量
f_M3C <- f_M3C[rownames(keep_genes_df), ]
# 将f_M3C数据框行名换成基因名
rownames(f_M3C) <- keep_genes_df$Name
for(i in 1:6){
  f_M3C[,i] <- as.numeric(f_M3C[,i])
}
f_M3C <- t(f_M3C)


#' 第三步，代谢组数据每组有10个个体，因此将代谢组和转录组的个体ID进行匹配
b2 <- b2[, rownames(f_M3C)]
b2 <- as.data.frame(b2)
# 对代谢组数据进行转换
# Z score transformation 
z_score <- function(x){
  x <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}
b2 <- apply(b2, 2, z_score)
b2 <- t(b2)


#' 第四步，进行关联分析
# spearman correlation analysis
library(plyr)
library(psych)
library(reshape2)
corr_df <- corr.test(b2, f_M3C, method = "spearman", adjust = "fdr", alpha = .05)
corr_df_cor <- corr_df$r
corr_df_p <- corr_df$p

#corr_df_cor <- corr_df_cor[abs(corr_df_cor) > 0.6 | corr_df_p < 0.05 ]

#' jiont the r and p matrix
asso_df <- rbind(corr_df_cor,corr_df_p)
#' transpose the data frame
asso_df <- as.data.frame(t(asso_df))
# rename the column variables
colnames(asso_df) <- paste(colnames(asso_df),rep(c("cor","P_adj"),c(23,23)),sep = "_")
#' re-order the asso_df according to the adjusted p value
#asso_df <- asso_df[order(asso_df$P_adj,decreasing = FALSE),]
# write.csv(asso_df, "M3C_metabo_transc_asso.csv",row.names = TRUE)

## use ggplot2
# Reset rownames
corr_df_cor <- data.frame(row=rownames(corr_df_cor),corr_df_cor,check.names = F) # create a column called "row" 
rownames(corr_df_cor) <- NULL
corr_df_p <- data.frame(row=rownames(corr_df_p),corr_df_p,check.names = F) # create a column called "row" 
rownames(corr_df_p) <- NULL
# Melt
nbar.m <- melt(corr_df_cor)
nbap.m <- melt(corr_df_p)
# Classify (you can classify differently for nbar and for nbap also)         
nbar.m$value2<-cut(nbar.m$value,breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1),include.lowest=TRUE, label=c("(-0.75,-1)","(-0.5,-0.75)","(-0.25,-0.5)","(0,-0.25)","(0,0.25)","(0.25,0.5)","(0.5,0.75)","(0.75,1)")) # the label for the legend
nbap.m$value2<-cut(nbap.m$value,breaks=c(-Inf, 0.001, 0.01, 0.05),label=c("***", "**", "*")) 
nbar.m<-cbind.data.frame(nbar.m,nbap.m$value,nbap.m$value2) # adding the p value and its cut to the first dataset of R coefficients
names(nbar.m)[5]<-paste("valuep") # change the column names of the dataframe 
names(nbar.m)[6]<-paste("signif.")
#nbar.m$row <- factor(nbar.m$row, levels=rev(unique(as.character(nbar.m$variable)))) # reorder the variable factor
# Plotting the matrix correlation heatmap
# Set options for a blank panel
nbar.m <- nbar.m[abs(nbar.m$value) > 0.6, ]
nbar.m <- nbar.m[!(nbar.m$signif. %in% NA), ]
pa <- ggplot(nbar.m, aes(row, variable)) +
  geom_tile(aes(fill=value),colour="white") +
  #scale_fill_brewer(palette = "RdYlGn",name="Correlation")# RColorBrewer package
  scale_fill_gradient2(low="blue", high="red", guide="colorbar",name="correlation") +
  theme_classic() +
  theme(axis.text.x=element_text(face="bold",angle=45,color="black",vjust = 0.95,hjust = 0.95,size=13),
        axis.text.y=element_text(face = "bold",size = 12, color = "black"),
        axis.title.y=element_text(size = 20),axis.title.x=element_text(size = 20))+
  labs(title="Genes associated with metabolites") + 
  xlab("Metabolites") +
  ylab("Genes within NF-kappa B") +
  theme(plot.title = element_text(size = 18,color = "black", hjust = 0.5)) +
  theme(axis.line.x=element_line(colour = "black"),axis.line.y=element_line(colour = "black"),
        legend.key.height=unit(0.8,"cm"),
        legend.key.width=unit(0.8,"cm"),
        legend.text=element_text(lineheight=0.8,face="bold",size=12),
        legend.title=element_text(size=15))
# Adding the significance level stars using geom_text 
pp<- pa +
  geom_text(aes(label=signif.),size=5,color = "white",na.rm=TRUE,nudge_y = -0.15)
ggsave(paste(temp_dir2, "M3C_metabo_transc2.png"), plot = pp, width = 10,height = 10)



