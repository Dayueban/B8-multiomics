## Dataset one

rm(list = ls())
options(stringsAsFactors = FALSE)
library(ggpubr)
## import dataset
#Two_variable_one <- read.table("ileum_two_variable_one.csv",header = T,
#                               check.names = FALSE,sep = ",")
Two_variable_two <- read.csv("taxonomy/Genus/genus_assoc_Ct.csv",header = T,
                               check.names = FALSE)
# In order to fetch the name of the ylab
# fake_data <- read.table("grey60_module_associated_fake.csv",header = T,
#                        check.names = FALSE,sep = ",")

#Min-Max Normalization
trans_max_min <- function(x){
  (x-min(x))/(max(x)-min(x))
}
#Two_variable_one[,2:5] <- as.data.frame(apply(Two_variable_one[,2:5],2,trans_max_min))
Two_variable_two[,2:ncol(Two_variable_two)] <- as.data.frame(apply(Two_variable_two[,2:ncol(Two_variable_two)],2,trans_max_min))
myplots <- list()  # new empty list
for (i in 1:2) {
  p1 <- ggscatter(Two_variable_two,x= "Ct value",y=colnames(Two_variable_two)[i+1],
                  color = "black",fill="black", shape = 21, size = 1.5, # Points color, shape and size
                  add.params = list(color = "red", fill = "grey66",size = 2), # Customize reg. line
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "spearman",
                  cor.coeff.args = list(method = "spearman", label.x = 0.25, label.sep = "\n",size=6),
                  xlab = "Ct value", ylab = colnames(Two_variable_two)[i+1])
  p1 <- ggpar(p1,
              font.tickslab = c(14,"bold", "black"),font.x = c(18, "bold"),
              font.y = c(18, "bold"))
  print(i)
  print(p1)
  myplots[[i]] <- p1  # add each plot into plot list
}
#Arrange on one page
p_total <- ggarrange(plotlist = myplots,labels = "B",ncol = 2,nrow = 1)
p_total
ggsave("taxonomy/Genus/Ct.value.assoc.genus.tiff",width = 14,height = 8)
