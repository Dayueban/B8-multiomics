# dotplot or barplot of enrichment analysis for metabolomics
rm(list = ls())
library(argparser)
library(ggplot2)
a <- read.csv("Model1 vs Control/1-差异分析/enrichment_Model1_DOWN.csv", header = T, check.names = F)
a$Enrichment_ratio <- a$hits / a$expected

a <- a[order(a$Enrichment_ratio, decreasing = T), ]

#Drawing enrichment scatter plot
p <- ggplot(a, aes(x = reorder(pathway, Enrichment_ratio), y = Enrichment_ratio, colour = Raw_p, size = Enrichment_ratio))
p <- p + geom_point() +
  coord_flip()
p <- p + scale_colour_gradientn(colours = rainbow(4), guide = "colourbar") + expand_limits(color = seq(0,1,by = 0.25))
p <- p + ggtitle("Enriched Metabolites Sets") + xlab("Pathway") +ylab("Enrichment ratio")
p <- p + theme_bw() + theme(axis.text = element_text(color = "black", size = 12),
                            axis.title = element_text(color = "black", size = 14))
p <- p + theme(panel.border = element_rect(colour = "black"))
p <- p + theme(plot.title = element_text(hjust = 0.5), legend.key = element_blank())
ggsave("Model1 vs Control/1-差异分析/enrichment_pathway_DOWN.png", plot = p, width = 8, height = 6, type = 'cairo-png')
