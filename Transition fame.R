library("FactoMineR")
library("factoextra")
library(ade4)
library(MASS)
library(ellipse)
library(ggplot2)
library(agricolae)
library(reshape2)
library(vegan)
library(ggthemes)
library(RVAideMemoire)
library(dplyr)

# Set working directory
setwd("G:/My Drive/labs/NMSU/Transtion/FAME/")#"C:/Users/juanp/Google Drive/labs/NMSU/Transtion/Bacteria/"

Metadata <- read.table("Metadata.txt", header = TRUE)
FAME <- read.table("FAME.txt", header = TRUE)

DATA_FAME = cbind(Metadata, FAME[1:32,])

#PERMANOVA 
FAME.dist<-vegdist(FAME[1:32,], method='bray')
PERMANOVA <- adonis2(FAME.dist~Treatment*Irrigation, DATA_FAME, permutation = 999, method ='bray')
summary(PERMANOVA)

PERMANOVA2 = data.frame(PERMANOVA)
variable  = c('variable','variable','variable','variable','variable')
PERMANOVA2$Significance <- "No Significant"
pval_thres <- 0.05
PERMANOVA2$Significance[which(PERMANOVA2$Pr..F. < pval_thres)] <- "Significant"
PERMANOVA2$Significance <- PERMANOVA2$Significance %>% factor

LC = data.frame(cbind(variable, row.names(PERMANOVA2), PERMANOVA2[,c(3,6)]))
colnames(LC) = c('variable','Effect', 'R2','Significance')

#Plot Variance
dim(LC)
colnames(LC)
dim(LC)
data_melt <- LC[-(5),]

#Reorder x axis: the order of the factor will be the same as in the data.csv file
data_melt$Effect <- as.character(data_melt$Effect)#Turn your 'treatment' column into a character vector
data_melt$Effect <- factor(data_melt$Effect, levels=unique(rev(data_melt$Effect)))#Then turn it back into a factor with the levels in the correct order


mypal = c('white','#c0047b','#eaabcd','#a9d11c','#528501','#2d4800','#dcda65','#53c0cc','#1f4f58')

ggplot(data=data_melt, aes(x=variable, y=R2,fill=Effect)) +
  geom_bar(stat="identity",aes(color=Significance), size = 0.3,width = 1,linewidth=0.5) +# 
  scale_fill_manual(values = mypal) +
  theme_few() + #guides() + #color = FALSE, fill=FALSE
  labs(y="Variance explained (%)") +
  scale_color_manual(values = c('black',"red"),na.value =  "transparent",name = "Significance vs Water")+
  theme(legend.title=element_blank(), legend.margin=margin(c(0,0,0,0)),
        axis.text.x   = element_blank())


safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499") 
DATA_FAME$Treatment= factor(DATA_FAME$Treatment, levels=c("NCC","GBL_NR","GBL_R","G"))

#Total FAME
head(DATA_FAME)
aovTotal_FAME_I <- aov(TotalBiomass~Irrigation, data=DATA_FAME[DATA_FAME$Treatment!='G',])
summary(aovTotal_FAME_I)
#             Df Sum Sq Mean Sq F value Pr(>F)
#Irrigation   1    789   789.2   0.336  0.568
#Residuals   22  51663  2348.3

aovTotal_FAME_Ir <- aov(TotalBiomass~Treatment, data=DATA_FAME[DATA_FAME$Irrigation != 'Dryland',])
summary(aovTotal_FAME_Ir)
#Df Sum Sq Mean Sq F value Pr(>F)  
#Treatment    3  12419    4140   3.038 0.0706 .
#Residuals   12  16351    1363

aovTotal_FAME_Dr <- aov(TotalBiomass~Treatment, data=DATA_FAME[DATA_FAME$Irrigation == 'Dryland',])
summary(aovTotal_FAME_Dr)
#Df Sum Sq Mean Sq F value Pr(>F)  
#Treatment    3  17743    5914   3.277 0.0587 .
#Residuals   12  21659    1805 


ggplot(data=DATA_FAME, aes(x=Irrigation, y=TotalBiomass, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  #facet_grid(Depth~.,space = "fixed",scales = "free_y") +
  ggtitle("Total Biomass")+ ylab("nmol / g soil")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Year")+
  theme_bw()


#GNegSum
head(DATA_FAME)
aovGNegSum <- aov(GNegSum  ~Treatment*Irrigation, data=DATA_FAME)
aovGNegSum_I <- aov(GNegSum~Irrigation, data=DATA_FAME[DATA_FAME$Treatment!='G',])
summary(aovGNegSum_I)
#             Df Sum Sq Mean Sq F value Pr(>F)
#Irrigation   1  0.866  0.8664   1.638  0.214
#Residuals   22 11.638  0.5290

aovGNegSum_Ir <- aov(GNegSum~Treatment, data=DATA_FAME[DATA_FAME$Irrigation != 'Dryland',])
summary(aovGNegSum_Ir)
#Df Sum Sq Mean Sq F value Pr(>F)  
#Treatment    3  2.940  0.9801   2.193  0.142
#Residuals   12  5.363  0.4469 

aovGNegSum_Dr <- aov(GNegSum~Treatment, data=DATA_FAME[DATA_FAME$Irrigation == 'Dryland',])
summary(aovGNegSum_Dr)
#Df Sum Sq Mean Sq F value Pr(>F)  
#Treatment     3  3.130  1.0433   2.814 0.0844 .
#Residuals   12  4.449  0.3707  

ggplot(data=DATA_FAME, aes(x=Irrigation, y=GNegSum, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  #facet_grid(Depth~.,space = "fixed",scales = "free_y") +
  ggtitle("Gram Neg")+ ylab("nmol / g soil")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Year")+
  theme_bw()

#GPosSum
head(DATA_FAME)
aovGPosSum <- aov(GPosSum ~Treatment*Irrigation, data=DATA_FAME)
summary(aovGPosSum)
aovGPosSum_I <- aov(GPosSum~Irrigation, data=DATA_FAME[DATA_FAME$Treatment!='G',])
summary(aovGPosSum_I)
#             Df Sum Sq Mean Sq F value Pr(>F)
#Irrigation   1   1.7   1.693   0.114  0.739
#Residuals   22 326.5  14.841 

aovGPosSum_Ir <- aov(GPosSum~Treatment, data=DATA_FAME[DATA_FAME$Irrigation != 'Dryland',])
summary(aovGPosSum_Ir)
#Df Sum Sq Mean Sq F value Pr(>F)  
#Treatment    3  76.77  25.590   2.956 0.0753 .
#Residuals   12  103.88   8.657 

aovGPosSum_Dr <- aov(GPosSum~Treatment, data=DATA_FAME[DATA_FAME$Irrigation == 'Dryland',])
summary(aovGPosSum_Dr)
#Df Sum Sq Mean Sq F value Pr(>F)  
#Treatment     3  128.2   42.74   3.883 0.0376 *
#Residuals   12  132.1   11.00    

LSDoav <- HSD.test(aovGPosSum_Dr, c("Treatment"), alpha = 0.05, group=TRUE) 
LSDoav

#GPosSum groups
#G      20.62575      a
#GBL_R  19.67500      ab
#GBL_NR 15.09825      b
#NCC    14.05675      b


ggplot(data=DATA_FAME, aes(x=Irrigation, y=GPosSum, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  #facet_grid(Depth~.,space = "fixed",scales = "free_y") +
  ggtitle("Gram Pos")+ ylab("nmol / g soil")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Year")+
  theme_bw()

#ActSum 
head(DATA_FAME)
aovActSum  <- aov(ActSum~Treatment*Irrigation, data=DATA_FAME)
aovActSum_I <- aov(ActSum~Irrigation, data=DATA_FAME[DATA_FAME$Treatment!='G',])
summary(aovActSum_I)
#             Df Sum Sq Mean Sq F value Pr(>F)
#Irrigation   1    1.82   1.822   0.463  0.503
#Residuals   22 86.60   3.936

aovActSum_Ir <- aov(ActSum~Treatment, data=DATA_FAME[DATA_FAME$Irrigation != 'Dryland',])
summary(aovActSum_Ir)
#Df Sum Sq Mean Sq F value Pr(>F)  
#Treatment    3  14.38   4.794   1.564  0.249
#Residuals   12  36.78   3.065 

aovActSum_Dr <- aov(ActSum~Treatment, data=DATA_FAME[DATA_FAME$Irrigation == 'Dryland',])
summary(aovActSum_Dr)
#             Df Sum Sq Mean Sq F value Pr(>F)  
#Treatment     3  18.18   6.060   1.496  0.266
#Residuals   12  48.62   4.052      

ggplot(data=DATA_FAME, aes(x=Irrigation, y=ActSum , fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  #facet_grid(Depth~.,space = "fixed",scales = "free_y") +
  ggtitle("Actinomycetes")+ ylab("nmol / g soil")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Year")+
  theme_bw()

#BactSum
head(DATA_FAME)
aovBactSum <- aov(BactSum ~Treatment*Irrigation, data=DATA_FAME)
summary(aovBactSum)
aovBactSum_I <- aov(BactSum~Irrigation, data=DATA_FAME[DATA_FAME$Treatment!='G',])
summary(aovBactSum_I)
#             Df Sum Sq Mean Sq F value Pr(>F)
#Irrigation   1   12.8   12.83   0.308  0.585
#Residuals   22  917.4   41.70 

aovBactSum_Ir <- aov(BactSum~Treatment, data=DATA_FAME[DATA_FAME$Irrigation != 'Dryland',])
summary(aovBactSum_Ir)
#Df Sum Sq Mean Sq F value Pr(>F)  
#Treatment    3  203.1   67.72   2.441  0.115
#Residuals   12  332.9   27.74  

aovBactSum_Dr <- aov(BactSum~Treatment, data=DATA_FAME[DATA_FAME$Irrigation == 'Dryland',])
summary(aovBactSum_Dr)
#Df Sum Sq Mean Sq F value Pr(>F)  
#Treatment     3  3  296.9   98.98    2.88   0.08 .
#Residuals   12  412.4   34.37  

ggplot(data=DATA_FAME, aes(x=Irrigation, y=BactSum, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  #facet_grid(Depth~.,space = "fixed",scales = "free_y") +
  ggtitle("Total Bacteria")+ ylab("nmol / g soil")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Year")+
  theme_bw()

#FungSum
head(DATA_FAME)
aovFungSum <- aov(FungSum~Treatment*Irrigation, data=DATA_FAME)
summary(aovFungSum)
aovFungSum_I <- aov(FungSum~Irrigation, data=DATA_FAME[DATA_FAME$Treatment!='G',])
summary(aovFungSum_I)
#             Df Sum Sq Mean Sq F value Pr(>F)
#Irrigation   25.2   25.21   0.529  0.475
#Residuals   22 1047.8   47.63                

aovFungSum_Ir <- aov(FungSum~Treatment, data=DATA_FAME[DATA_FAME$Irrigation != 'Dryland',])
summary(aovFungSum_Ir)
#Df Sum Sq Mean Sq F value Pr(>F)  
#Treatment    3  263.6   87.85   3.311 0.0572 .
#Residuals   12  318.4   26.53 

aovFungSum_Dr <- aov(FungSum~Treatment, data=DATA_FAME[DATA_FAME$Irrigation == 'Dryland',])
summary(aovFungSum_Dr)
#Df Sum Sq Mean Sq F value Pr(>F)  
#Treatment    3  319.4  106.48   2.755 0.0885 .
#Residuals   12  463.7   38.65  

ggplot(data=DATA_FAME, aes(x=Irrigation, y=FungSum, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  #facet_grid(Depth~.,space = "fixed",scales = "free_y") +
  ggtitle("Total Fungi")+ ylab("nmol / g soil")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Year")+
  theme_bw()

#AMF
head(DATA_FAME)
aovAMF <- aov(AMF~Treatment*Irrigation, data=DATA_FAME)
summary(aovAMF)
aovAMF_I <- aov(AMF~Irrigation, data=DATA_FAME[DATA_FAME$Treatment!='G',])
summary(aovAMF_I)
#             Df Sum Sq Mean Sq F value Pr(>F)
#Irrigation    1  22.12  22.124   7.713  0.011 *
#Residuals   22  63.10   2.868               

aovAMF_Ir <- aov(AMF~Treatment, data=DATA_FAME[DATA_FAME$Irrigation != 'Dryland',])
summary(aovAMF_Ir)
#Df Sum Sq Mean Sq F value Pr(>F)  
#Treatment    3  34.11  11.369   3.512 0.0492 *
#Residuals   12  38.84   3.237 
LSDoav_IR <- HSD.test(aovAMF_Ir, c("Treatment"), alpha = 0.05,  group=TRUE) 
LSDoav_IR

#AMF groups
#GBL_R  7.27900      a
#G      6.39575     ab
#GBL_NR 5.86650     ab
#NCC    3.34975      b

aovAMF_Dr <- aov(AMF~Treatment, data=DATA_FAME[DATA_FAME$Irrigation == 'Dryland',])
summary(aovAMF_Dr)
#Df Sum Sq Mean Sq F value Pr(>F)  
#Treatment    3  8.918  2.9725   4.202 0.0301 *
#Residuals   12  8.490  0.7075 

LSDoav <- HSD.test(aovAMF_Dr, c("Treatment"), alpha = 0.05, group=TRUE) 
LSDoav

#AMF groups
#GBL_R  4.31675      a
#GBL_NR 4.04550     ab
#G      3.45575     ab
#NCC    2.37225      b

ggplot(data=DATA_FAME, aes(x=Irrigation, y=AMF, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  #facet_grid(Depth~.,space = "fixed",scales = "free_y") +
  ggtitle("AMF")+ ylab("nmol / g soil")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Year")+
  theme_bw()

#SPSum
head(DATA_FAME)
aovSPSum <- aov(SPSum~Treatment*Irrigation, data=DATA_FAME)
summary(aovSPSum)

aovSPSum_I <- aov(SPSum~Irrigation, data=DATA_FAME[DATA_FAME$Treatment!='G',])
summary(aovSPSum_I)
#             Df Sum Sq Mean Sq F value Pr(>F)
#Irrigation    1  0.1    0.10   0.003  0.956
#Residuals   22  698.1   31.73                 

aovSPSum_Ir <- aov(SPSum~Treatment, data=DATA_FAME[DATA_FAME$Irrigation != 'Dryland',])
summary(aovSPSum_Ir)
#Df Sum Sq Mean Sq F value Pr(>F)  
#Treatment    3  111.6   37.21   2.871 0.0806 .
#Residuals   12  155.6   12.96

aovSPSum_Dr <- aov(SPSum~Treatment, data=DATA_FAME[DATA_FAME$Irrigation == 'Dryland',])
summary(aovSPSum_Dr)
#Df Sum Sq Mean Sq F value Pr(>F)  
#Treatment    3  252.0   84.00   2.656  0.096 .
#Residuals   12  379.6   31.63 

ggplot(data=DATA_FAME, aes(x=Irrigation, y=SPSum, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  #facet_grid(Depth~.,space = "fixed",scales = "free_y") +
  ggtitle("Saprophytic Fungi")+ ylab("nmol / g soil")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Year")+
  theme_bw()

#FBRatio      
head(DATA_FAME)
aovFBRatio <- aov(FBRatio~Treatment*Irrigation, data=DATA_FAME)
summary(aovFBRatio)

aovFBRatio_I <- aov(FBRatio~Irrigation, data=DATA_FAME[DATA_FAME$Treatment!='G',])
summary(aovFBRatio_I)
#             Df Sum Sq Mean Sq F value Pr(>F)
#Irrigation    1  0.0076 0.007597   0.817  0.376
#Residuals   22 0.2045 0.009294                 

aovFBRatio_Ir <- aov(FBRatio~Treatment, data=DATA_FAME[DATA_FAME$Irrigation != 'Dryland',])
summary(aovFBRatio_Ir)
#Df Sum Sq Mean Sq F value Pr(>F)  
#Treatment    3  0.01457 0.004857   0.474  0.706
#Residuals   12 0.12297 0.01024

aovFBRatio_Dr <- aov(FBRatio~Treatment, data=DATA_FAME[DATA_FAME$Irrigation == 'Dryland',])
summary(aovFBRatio_Dr)
#Df Sum Sq Mean Sq F value Pr(>F)  
#Treatment    3  0.05572 0.018574   4.244 0.0292 *
#Residuals   12 0.05252 0.004377 

LSDoav <- HSD.test(aovFBRatio_Dr, c("Treatment"), alpha = 0.05, group=TRUE) 
LSDoav

#FBRatio groups
#GBL   0.8965      a
#GBL-Rem  0.8225     ab
#PG       0.7715     ab
#NCC     0.7405      b

ggplot(data=DATA_FAME, aes(x=Irrigation, y=FBRatio, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  #facet_grid(Depth~.,space = "fixed",scales = "free_y") +
  ggtitle("Fungal/Bacteria Ratio")+ ylab("nmol / g soil")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Year")+
  theme_bw()


##############################################
res.pca <- prcomp(FAME[1:32,], scale = TRUE)
fviz_eig(res.pca)
summary(res.pca)

load$PC1= data.frame(res.pca$rotation*5)

title ="PCA - PCA1 vs PCA2"
ev1.F <- 76.91 
ev2.F <- 13.71 
xlab_text.F <- paste("PCA1 (", round(ev1.F,2), "%)")
ylab_text.F <- paste("PCA2 (", round(ev2.F,2), "%)")

f<-data.frame(res.pca$x[,1:2])
f <- cbind(Metadata,f)
rownames(f)<-f$Row.names
f$Row.names <- NULL

PCA_Comp_mean = aggregate(cbind(mean.x=PC1,mean.y=PC2)~Treatment*Irrigation,f,mean)
PCA_Comp_meands = aggregate(cbind(ds.x=PC1,ds.y=PC2)~Treatment*Irrigation,f,se)
PCoA_Comp_mean2 =cbind(PCA_Comp_mean,PCA_Comp_meands[,c(3:4)])


paleta_alive <- c('#FF0000','#00008B',"#FFB919","#00CC1C")

ggplot(data = PCoA_Comp_mean2,aes(mean.x,mean.y)) + 
  geom_vline(xintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_hline(yintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_errorbarh(mapping = aes(xmin =mean.x - ds.x, xmax = mean.x + ds.x),size = 0.01,alpha = 0.1) +
  geom_errorbar(mapping = aes(ymin =mean.y - ds.y, ymax = mean.y + ds.y),size = 0.01,alpha = 0.1)  + 
  geom_point(size=2.5, aes(color = Treatment,  shape=Irrigation),stroke = 1) + #size = 4,
  labs(x=xlab_text.F, 
       y=ylab_text.F, 
       title=paste("FAME analysis")) +
  scale_fill_manual(values = safe_colorblind_palette) +
  scale_color_manual(values = safe_colorblind_palette)+
  theme_few()

ggplot(data = f,aes(PC1,PC2)) + #PCoA_Comp.F_mean2 mean.x mean.y
  geom_segment(data=load,aes(x=0,y=0,xend=PC1,yend=PC2),arrow=arrow(length=unit(.3,"lines")),color="black",size=.4)+
  annotate("text",x=load$PC1*1.5,y=load$PC2,label=rownames(load),size=3)+
  geom_vline(xintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_hline(yintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_point(size=5, aes(color = Treatment, shape=Irrigation),stroke = 1) + #size = 4, #,  shape=Year
  labs(x=xlab_text.F, y=ylab_text.F) +
  scale_fill_manual(values = safe_colorblind_palette) +
  scale_color_manual(values = safe_colorblind_palette)+
  theme_bw()+ theme(legend.position = 'right')

##############RDA###########
par(mfrow=c(1,1))
Bact_Fame <- rda(FAME[1:32,]~., chemistry[(33:64),-(1:5)], scale=T)
summary(Bact_Fame)
Bact_Fame$CCA$u

L<-data.frame(rownames(Bact_Fame$CCA$u),Bact_Fame$CCA$u[, 1:2])
colnames(L) =c('ID', 'RDA1','RDA2')
sam = cbind(rownames(chemistry[,1]), ps.2@sam_data)
colnames(sam) =c('ID', 'Plot','Irrigation','Treatment','Plants','Residue','Rep')
L <- merge(L,sam, by='ID')
rownames(L)<-L$ID

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499") 

L$Treatment = factor(L$Treatment, c('NCC','GBL_NR','GBL_R','G'))

load.FAME = data.frame(Bact_Fame$CCA$biplot[,1:2])


RDAFAME = ggplot(L, aes(RDA1, RDA2)) + 
  geom_vline(xintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_hline(yintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_point(aes(colour=Treatment, shape=Irrigation),size=5) +
  theme_bw() +
  theme(legend.position="right")+ 
  labs(y='RDA2 (3.48%)', x='RDA2 (27.28%)')+
  scale_colour_manual(values = safe_colorblind_palette)+
  geom_segment(data=load.FAME,aes(x=0,y=0,xend=RDA1,yend=RDA2),arrow=arrow(length=unit(.3,"lines")),color="black",size=0.8)+
  annotate("text",x=load.FAME$RDA1, y=load.FAME$RDA2,label=rownames(load.FAME),size=4)

anova (Bact_Fame, by = 'margin', parallel = 8)
#PMC       1    203.5 1.1827  0.036 *


### 
library(ggcorrplot)
ncol(chemistry)
TablaCompleta = cbind(chemistry[33:64,],FAME[1:32,])
Corr <-cor(TablaCompleta[,-c(1:5)], method = 'spearman')
p.mat <-cor_pmat(TablaCompleta[,-c(1:5)], method = 'spearman')

ncol(TablaCompleta[,-c(1:5)])

Corr1 = as.matrix(Corr[1:11,-c(1:11)])
p.mat1 = as.matrix(p.mat[1:11,-c(1:11)])
ncol(Corr1)
nrow(Corr1)

Corr2 = as.matrix(Corr[1:12,-c(1:9)])
p.mat2 = as.matrix(p.mat[2:10,-c(1:9)])


ggcorrplot(Corr1,
           p.mat = p.mat1, insig = "blank",
           method = "circle", 
           hc.order = TRUE, #type = "upper",
           #outline.color = "white",
           ggtheme = ggplot2::theme_bw,
           colors = c("#6D9EC1", "white", "#E46726"))
