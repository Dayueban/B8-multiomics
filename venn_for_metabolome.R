rm(list = ls())
library(stringr)
library(qpcR)
library(readxl)

metabo <- read_xlsx(path = "Con_M1_M2_M3_venn.xlsx", sheet = 1)

# intersection based on metabolites
# 24hpi
metabo1 <- metabo$Model1_Control
# 48hpi
metabo2 <- metabo$Model2_Control
metabo2 <- metabo2[!is.na(metabo2)] # 只保留不是NA的代谢物
# 72hpi
metabo3 <- metabo$Model3_Control
metabo3 <- metabo3[!is.na(metabo3)]
## second part: draw venn plot
# load library
library(stringr)
library(qpcR)
library(VennDiagram)

## Application on rap Lyrics
# protein with upregulated
library(tidyverse)
library(hrbrthemes)
library(tm)
library(proustr)
venn.diagram(
  x = list(
    metabo1, metabo2, metabo3
  ),
  category.names = c("Model1 vs Control" , "Model2 vs Control" , "Model3 vs Control"),
  filename = "venn_protein_UP.png2.png",
  output = TRUE ,
  imagetype="png" ,
  height = 1000 , 
  width = 1000 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#7292AA","#90CBD3","#EAA58E"),
  fill = c(alpha("#7292AA",0.8), alpha('#90CBD3',0.8), alpha('#EAA58E',0.8)),
  cex = 0.8,
  fontfamily = "sans",
  cat.cex = 0.8,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  #cat.col = c("#7292AA","#90CBD3","#EAA58E"), # 标签的颜色
  cat.col = c("black", "black", "black"),
  rotation = 1
)
