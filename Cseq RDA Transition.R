##############RDA###########
setwd("G:/My Drive/labs/NMSU/Transtion/")
Table <- read.table("RDA C seq.txt", header = TRUE)

Cseq <- rda(Table[-30,3:6]~., Table[-30,-(1:6)], scale=T)
summary(Cseq)

r<-data.frame(rownames(Cseq$CCA$u),Cseq$CCA$u[, 1:2])
colnames(r) =c('ID', 'RDA1','RDA2')
r <- cbind(Table[-30,1:3],r)
rownames(r)<-r$ID

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499") 

r$Treatment = factor(r$Treatment, c('NCC','GBL-Rem','GBL','GO','PG'))

load = data.frame(Cseq$CCA$biplot[,1:2])


ggplot(r, aes(RDA1, RDA2)) + 
  geom_vline(xintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_hline(yintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_point(aes(colour=Treatment, shape=Irrigation),size=5) +
  theme_bw() +
  theme(legend.position="none")+ 
  labs(x='RDA1(48.9%)', y='RDA (25.46%)')+
  scale_colour_manual(values = safe_colorblind_palette)+
  geom_segment(data=load,aes(x=0,y=0,xend=RDA1,yend=RDA2),arrow=arrow(length=unit(.3,"lines")),color="black",linewidth=0.8)+
  annotate("text",x=load$RDA1, y=load$RDA2,label=rownames(load),size=4)
  
Cseq_aov = anova (Cseq, by = 'margin', parallel = 8)

C_importance = data.frame(cbind(rownames(Cseq_aov),Cseq_aov[,2],Cseq_aov[,4])) 
colnames(C_importance) = c('V1','Variance', 'Signification')
C_importance$Variance = as.numeric(C_importance$Variance)
C_importance$V1 = factor(C_importance$V1, rev(c("Ascomycota","Mortierellomycota","Bacteroidota","MBC","Soil_Saprotroph",
                                                "Proteobacteria","TotalBiomass","FBRatio","Kickxellomycota",
                                                "PMC","AMF","NAG","Actinobacteriota","BG","NT0",'Dung_Saprotroph',"PMN","TN","Residual")))


ggplot(data=C_importance[-19,], aes(x=V1, y=Variance, fill=Variance)) +
  geom_bar(stat="identity")+
  theme_bw()+ labs(y= 'Variance',x= 'Parameters')+
  coord_flip()+ theme(legend.position = 'left')
