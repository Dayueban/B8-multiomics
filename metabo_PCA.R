## PCA analysis for metabolomics
## produce Figure 1B
rm(list=ls())
library(ade4)
library(RColorBrewer)
library(vegan)

## group male before and after castration
df <- read.csv("../负离子模式/metabolome_neg_PCA.csv",header=T,row.names = 1)
df_tr <- as.data.frame(t(df))
# 离差标准化
max_min <- function(x){
  x <- (x-min(x,na.rm = T)) / (max(x,na.rm = T)-min(x,na.rm = T))
}
df_tr_Nor <- apply(df_tr,2,max_min)

############
Group <- rep(c("Control","Model1","Model2","Model3"),c(10,10,10,10))
Group <- as.factor(Group)
label <- factor(Group)
#legend <- unique(Group)
legend <- levels(label)
k <- length(legend)
if (k <= 4){
  colors =c("#E16663", "#7292AA","#90CBD3","#EAA58E")
} else if (k <= 9 && k >4) {
  colors = brewer.pal(k, "Set1")
} else {
  colors = rainbow(k, star=0)
}



#    pairs  F.Model        R2 p.value p.adjusted
#1 NC_bbc_male vs NC_bbc_female 6.918418 0.1003856   0.001      0.001
df_tr_Nor2 <- vegdist(df_tr_Nor,method = "bray")
df_tr_Nor2 <- as.dist(df_tr_Nor2)
df.dudi <- dudi.pco(df_tr_Nor2,scannf=F,nf=2)
ratio <- inertia.dudi(df.dudi)
pco1 <- as.numeric(sprintf("%0.2f",(ratio$tot.inertia)[1,3]))##pca1用于坐标轴名称,这里计算的是解释变量。解释下 %0.2f 的含义：% 表示起始字符0 表示空位用0填满2 表示小数点后必须占两位f 表示转换成浮点数

pco2 <- as.numeric(sprintf("%0.2f",((ratio$tot.inertia)[2,3]-(ratio$tot.inertia)[1,3])))##同上

xmax <- max(df.dudi$li[,1])##到时候用于坐标轴刻度
# xmax <- 0.4
xmin <- min(df.dudi$li[,1])
ymax <- max(df.dudi$li[,2])
ymin <- min(df.dudi$li[,2])

tiff("../负离子模式/neg_metabo_PCA_anosim_0810.tiff",width=2000,height=1600,res=300)
#par(mfrow=c(1,2))
layout(matrix(c(1,1,2), 1, 3, byrow = TRUE))
plot(df.dudi$li[,1:2],type="n",frame.plot=T,xlab=paste("PC1","(",pco1,"%",")",sep=""),
     ylab=paste("PC2","(",pco2,"%",")",sep=""),tcl=0.2,mgp=c(2.2,0.2,0),main="Negative mode",
     las=1,xlim=c(1.1*xmin,1.1*xmax),ylim=c(1.1*ymin,1.1*ymax),cex.axis=1.2,
     cex.lab=1.2,cex.main=1.4)##frame.plot 是否给图形加框
abline(h=0,v=0,col=rgb(0,0,0,0.5))
axis(3,tcl=0.2,labels=F)
axis(4,tcl=0.2,labels=F)

s.class(df.dudi$li,label,cpoint = 1.2,addaxes=T,axesell=T,origin = c(0,0),
        col=colors,grid=F,pch=16,add.plot=T)
if (k >= 4) {
  legend("topleft",legend=legend,pch=16,col=colors,bty="n",cex=1.2)
} else {}
#dev.off()

df_ano <- anosim(df_tr_Nor,grouping = Group,permutations = 999,distance = "bray")
summary(df_ano)
# ANOSIM statistic R: 0.4148 
# Significance: 0.003 
#tiff("../负离子模式/neg_anosim_metabo_0810.tiff",height = 2000,width = 1000,res=300)
par(mar=c(6,4,3,2))
plot(df_ano,col=c("#845EC2","#E16663", "#7292AA","#90CBD3","#EAA58E"),title="Anosim anslysis",cex.lab=1.3,
     xlab = "Group", ylab = "Dissimilarity ranks between and within groups",
     cex.main=0.8,cex.sub=1.6,cex.axis=1.4,xaxt="n")
labs <- c("Between","Control","Model1","Model2","Model3")
#title(xlab = "", ylab = "Dissimilarity ranks between and within groups")
axis(1,at=seq_along(labs),labels=FALSE,xpd = TRUE)
text(x=seq_along(labs), y=par("usr")[3] - (par("usr")[4]-par("usr")[3])*0.02, srt = 30, adj=1,labels=labs, xpd = TRUE,cex=1.2)
dev.off()


