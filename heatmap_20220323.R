# More rigorous screening of metabolites
rm(list = ls())
library(pheatmap)
a <- read.csv("Model1 vs Control/1-差异分析/Different_origin.csv", header = T,
              row.names = 1, check.names = F)
b <- read.csv("Model1 vs Control/1-差异分析/Different_metabo_intensity.csv",
              header = T, row.names = 1, check.names = F)
# select rows with multiple conditions
keep_metabo <- rownames(a[abs(a$log2FC) > 1 & a$FDR <= 0.05 & a$VIP > 1, ])
# match the dataframe b with specified metabolites
b2 <- b[keep_metabo,]
write.csv(b2, "Model1 vs Control/Heatmap_differ_metabo.csv") # 筛选出来的代谢物可以进行富集分析

# heatmap plots
b2 <- read.csv("Model3 vs Control/1-差异分析/Heatmap_differ_metabo3.csv", header = T, row.names = 1,
               check.names = F)
# scale before heatmap clustering
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
#fec_conj <- read.csv("MSMS/msms_Intensity2.csv")
#fec_conj_num <- read.csv("MSMS/msms_Intensity2_1.csv")
#type <- c(rep("glucuronide",13),rep("sulfate",20))
# 运行一次结果后发现M015个体和SV组聚成一个cluster，因此建议删除
#fec_conj <- fec_conj[,-16]
#fec_conj_num <- fec_conj_num[,-15]
sample <- c(rep("Control",10),rep("Model3",10))
type <- b2$type
fecPhase <- b2 
#rownames(fecPhase) <- rownames(b2)
fecPhase$type <- NULL
fecPhase$change <- NULL
fecPhase_norm <- t(apply(fecPhase, 1, cal_z_score))
#fecPhase_samplerow <- data.frame(sample)
fecPhase_samplerow <- data.frame(type)
#fecPhase_typecol <- data.frame(type)
fecPhase_typecol <- data.frame(sample)
row.names(fecPhase_samplerow) <- rownames(fecPhase)
row.names(fecPhase_typecol) <- colnames(fecPhase)
ann_colors <- list(sample = c(Control = "#E16663", Model3 = "#EAA58E"), 
                   type = c(Amino_acid = "#cc340c", Carbohydrate = "#3f60aa",
                            Cofactors_and_Vitamins = "#f18800", Fatty_acids= "#f88421",
                            Lipid = "#e4ce00", Nucleotide = "#9ec417",
                            Peptide = "#006b7b", Unknow = "#13a983",
                            Xenobiotics = "#44c1f0"))
#ann_colors <- list(sample = c(Control = "#E16663", Model3 = "#EAA58E"))
fecPhase_norm_t <- t(fecPhase_norm)
fecPhase_pheatmap <- pheatmap(fecPhase_norm_t,
                              color = colorRampPalette(c("navy", "#FEF9E7", "firebrick3"))(500),
                              annotation_col = fecPhase_samplerow,
                              annotation_row = fecPhase_typecol,
                              cutree_cols = 2,
                              cutree_rows = 2,
                              #cluster_cols = FALSE,
                              #cluster_rows = FALSE,
                              border_color = "NA",
                              fontsize_row = 7,
                              fontsize_col = 8,
                              annotation_colors = ann_colors
)
save_pheatmap_png <- function(x, filename, width=7200, height=4200, res = 600) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(fecPhase_pheatmap, "Model3 vs Control/1-差异分析/fecPhase_pheatmap0810.png")
