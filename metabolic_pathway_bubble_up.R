# dotplot or barplot of enrichment analysis for metabolomics
rm(list = ls())
library(argparser)
library(ggplot2)
# Model1 vs Control
a <- read.csv("Model1 vs Control/1-差异分析/enrichment_Model1_up.csv", 
              header = T, check.names = F)
a$Enrichment_ratio <- a$hits / a$expected
a <- a[order(a$Enrichment_ratio, decreasing = T), ]
dim(a)
a$Group <- rep("M1_vs_C", 26)

# Model2 vs Control
b <- read.csv("Model2 vs Control/1-差异分析/enrichment_Model2_up.csv",
              header = T, check.names = F)
b$Enrichment_ratio <- b$hits / b$expected
b <- b[order(b$Enrichment_ratio, decreasing = T),]
dim(b)
b$Group <- rep("M2_vs_C", 26)

# Model3 vs Control
c <- read.csv("Model3 vs Control/1-差异分析/enrichment_Model3_up.csv",
              header = T, check.names = F)
c$Enrichment_ratio <- c$hits / c$expected
c <- c[order(c$Enrichment_ratio, decreasing = T),]
dim(c)
c$Group <- rep("M3_vs_C", 19)


# 三个表进行合并
abc30 <- rbind(a[, c(1,5,8,9)], b[, c(1,5,8,9)], c[, c(1,5,8,9)])
#abc30 <- abc30 %>%
#  mutate(Count2 = Count / 10)


xx = ggplot(abc30, aes(x = Group, y = pathway)) + 
  geom_point(aes(size = Enrichment_ratio, fill = Raw_p), alpha = 0.75, shape = 21) + 
  #scale_size_continuous(limits = c(1, 130), range = c(1,13), breaks = c(1,5,10,13)) + 
  labs( x= "", y = "metabolic pathway with upregulated metabolites", size = "Enrichment ratio", fill = "P-value")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 45, vjust = 0.95, hjust = 0.95), 
        axis.text.y = element_text(colour = "black", size = 10), 
        axis.title = element_text(color = "black", face = "bold", size = 13),
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 11, face = "bold"), panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right", panel.grid.major.y = element_line(colour = "grey95")) +  
  scale_fill_continuous(type = "viridis")
xx
ggsave("metabolic_pathway_M1M2M3C_up_20220901.png", plot = xx, width = 6, height = 6, dpi = 300)
