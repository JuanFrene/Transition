packages <- c('ggthemes','dplyr', "ape", "ShortRead", "Biostrings", "phyloseq", 
              "corncob","DESeq2", "microbiome", "DECIPHER", "phangorn", "tibble", "lme4", 
              "lmerTest", "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire")
sapply(packages, require, character.only = TRUE)              

# Set working directory
setwd("G:/My Drive/labs/NMSU/Transtion/Fungi/")#"C:/Users/juanp/Google Drive/labs/NMSU/Transtion/Bacteria/"

# Set plotting theme
theme_set(theme_bw())
# Assign variables for imported data
sharedfile_fun = "final.asv.shared"
taxfile_fun = "final.asv.ASV.cons.taxonomy"

mothur_data_fun <- import_mothur(mothur_shared_file = sharedfile_fun,
                             mothur_constaxonomy_file = taxfile_fun)
mothur_data_fun@otu_table
mothur_data_fun@tax_table[1:17,]

#Agregar la data
Metadata <- read.table("Metadata.txt", header = TRUE)
sample_dat = sample_data(Metadata)
sample_names(sample_dat) = rownames(Metadata)

## Building the Phyloseq Object
physeq_fun = merge_phyloseq(mothur_data_fun, sample_dat)
ps.1 = prune_taxa(taxa_sums(physeq_fun) >= 1, physeq_fun)

# Now, let's make sure that the correct information is included in our phyloseq object
ps.F = subset_taxa(ps.1, Rank1 ==  "k__Fungi" )
ps.perc.F <- transform_sample_counts(ps.F, function(x) x / sum(x) * 100) 
ps.perc.F.2 = prune_taxa(taxa_sums(ps.perc.F) >= 0.01, ps.perc.F)
ps.F.H <- transform(ps.F, "hellinger")

ps.F.H_I = subset_samples(ps.F.H, Irrigation ==  "Irrigated" )
ps.F.H_D = subset_samples(ps.F.H, Irrigation !=  "Irrigated" )


###ANALYSIS###
# 1. PCoA plot####
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499") 

PCoA <- ordinate(ps.F.H_D, "PCoA", "bray")
title ="PCoA - PCoA1 vs PCoA2"
ev1 <- (PCoA$values$Relative_eig[1])*100
ev2 <- (PCoA$values$Relative_eig[2])*100
xlab_text <- paste("PCoA1 (", round(ev1,2), "%)")
ylab_text <- paste("PCoA2 (", round(ev2,2), "%)")

x<-data.frame(PCoA$vectors[, 1:2])
x <- merge(x,ps.perc.F.2@sam_data, by=0)
rownames(x)<-x$Row.names
x$Row.names <- NULL

x$Treatment = factor(x$Treatment, c('NCC', 'GBL_NR','GBL_R','G'))

ggplot(x, aes(Axis.1, Axis.2)) + 
  geom_vline(xintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_hline(yintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_point(aes(colour=Treatment, shape=Irrigation),size=5) +
  theme_bw() +
  #stat_ellipse(aes(group = Irrigation, color = Irrigation),size = 1)+
  theme(legend.position="right")+ 
  #scale_colour_paletteer_c("pals::ocean.speed",limits = c(0.2,0.8)) +
  #scale_colour_gradientn(colours=c("Blue","Red"))+
  labs(y=ylab_text, x=xlab_text)+
  scale_colour_manual(values = safe_colorblind_palette)


PCoAmean = aggregate(cbind(mean.x=Axis.1,mean.y=Axis.2)~Treatment*Irrigation,x,mean)
PCoAds = aggregate(cbind(ds.x=Axis.1,ds.y=Axis.2)~Treatment*Irrigation,x,sd)
PCoAmean <- cbind(PCoAmean, PCoAds[,3:4])

ggplot(data = PCoAmean, aes(mean.x, mean.y)) +
  geom_vline(xintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_hline(yintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_errorbarh(mapping = aes(xmin =mean.x - ds.x ,xmax = mean.x + ds.x),size = 0.1, alpha = 0.15) +
  geom_errorbar(mapping = aes(ymin =mean.y - ds.y ,ymax = mean.y + ds.y),size = 0.1, alpha = 0.15)  + 
  scale_colour_manual(values = safe_colorblind_palette)+
  geom_point(size = 5, aes(colour = Treatment, shape=Irrigation),stroke = 1) + 
  labs(y=ylab_text, x=xlab_text)+
  theme_few()


# 2. PERMANOVA
#PERMANOVA
set.seed(1)

# Calculate bray curtis distance matrix
ps.T3_bray <- phyloseq::distance(ps.F.H_I, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(ps.F.H_I))

# Adonis test
PERMANOVA = adonis2(ps.T3_bray ~ Treatment, data = sampledf)
#Df SumOfSqs      R2      F Pr(>F)    
#Treatment             3   1.3816 0.13655 1.4787  0.001 ***
#Irrigation            1   0.5732 0.05665 1.8404  0.001 ***
#Treatment:Irrigation  3   1.0001 0.09885 1.0704  0.146    
#Residual             23   7.1633 0.70796                  
#Total                30  10.1183 1.00000 

#Dryland
#Treatment  3   1.3872 0.21965 1.1259  0.001 ***
#Residual  12   4.9283 0.78035                  
#Total     15   6.3155 1.00000

#         Df SumOfSqs      R2      F Pr(>F)  
#Treatment  3   1.3327 0.21967 1.0322  0.037 *
#Residual  11   4.7343 0.78033                
#Total     14   6.0670 1.00000  

PERMANOVA2 = data.frame(PERMANOVA)
variable  = c('variable','variable','variable')
PERMANOVA2$Significance <- "No Significant"
pval_thres <- 0.05
PERMANOVA2$Significance[which(PERMANOVA2$Pr..F. < pval_thres)] <- "Significant"
PERMANOVA2$Significance <- PERMANOVA2$Significance %>% factor

LC = data.frame(cbind(variable, row.names(PERMANOVA2), PERMANOVA2[,c(3,6)]))
colnames(LC) = c('variable','Effect', 'R2','Significance')
LC$Effect[which(LC$Effect == 'Unknow')] <- "Conditions"
LC$Effect[which(LC$Effect == 'Species:Unknow')] <- "Species:Conditions"
LC$Effect[which(LC$Effect == 'Location:Unknow')] <- "Location:Conditions"
LC$Effect[which(LC$Effect == 'Species:Location:Unknow')] <- "Species:Location:Conditions"


#Plot Variance
dim(LC)
colnames(LC)
dim(LC)
data_melt <- LC[-(3),]

#Reorder x axis: the order of the factor will be the same as in the data.csv file
data_melt$Effect <- as.character(data_melt$Effect)#Turn your 'treatment' column into a character vector
data_melt$Effect <- factor(data_melt$Effect, levels=unique(rev(data_melt$Effect)))#Then turn it back into a factor with the levels in the correct order

mypal = c('white','#c0047b','#a9d11c','#1f4f58','#528501')

ggplot(data=data_melt, aes(x=variable, y=R2,fill=Effect)) +
  geom_bar(stat="identity",aes(color=Significance), size = 0.3,width = 1,linewidth=0.4) +# 
  scale_fill_manual(values = mypal) +
  theme_few() + #guides() + #color = FALSE, fill=FALSE
  labs(y="Variance explained (%)") +
  scale_color_manual(values = c('black',"red"),na.value =  "transparent",name = "Significance vs Water")+
  theme(legend.title=element_blank(), legend.margin=margin(c(0,0,0,0)),
        axis.text.x   = element_blank())

ggsave("Fungi_variance Irrigated.pdf", height=5, width=3.5, units='in')


# 3. Simper function
otudat <- t(otu_table(ps.perc.2))
otudat <- otudat[rowSums(otudat) > 0,]
ncol(otudat)
## in case you also have empty otu
otudat <- otudat[, colSums(otudat) > 0]
Simper_I = simper(otudat, sample_data(ps.perc.2)$Irrigation)
Simper_T = simper(otudat, sample_data(ps.perc.2)$Treatment)


table_I_genus = summary(Simper_I)$Dryland_Irrigated
table_T_genus_NCC_GBL_R = summary(Simper_T)$NCC_GBL_R
table_T_genus_NCC_G = summary(Simper_T)$NCC_G
table_T_genus_NCC_GBL_NR = summary(Simper_T)$NCC_GBL_NR
table_T_genus_GBL_R_G = summary(Simper_T)$GBL_R_G
table_T_genus_GBL_R_GBL_NR = summary(Simper_T)$GBL_R_GBL_NR
table_T_genus_G_GBL_NR = summary(Simper_T)$G_GBL_NR


Simplot <- read.table("Simplot.txt", header = TRUE)

ggplot(data = Simplot[Simplot$Kingdom=='Bacteria',], aes(Taxa,Indicator)) + #
  geom_raster(aes(fill = Porcentage))+
  theme_few() +
  facet_grid(.~Kingdom)+
  geom_tile(aes(color = Significance),fill = '#00000000', size = 0.2, linewidth = 0.9, width = 1,height = 1) + #
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_40_95_c42", #kovesi.diverging_bwr_55_98_c37
                         limits = c(0,4),na.value = "#D9D9D9")+
  scale_color_manual(values = c("black","black"),na.value =  "transparent",name = "Significance") + #Significance Genotype vs Col-0
  theme(axis.text.x = element_text(angle = -45, hjust=-0.05),axis.title = element_blank())

ggplot(data = Simplot[Simplot$Kingdom!='Bacteria',], aes(Taxa,Indicator)) + #
  geom_raster(aes(fill = Porcentage))+
  theme_few() +
  facet_grid(.~Kingdom)+
  geom_tile(aes(color = Significance),fill = '#00000000', size = 0.2, linewidth = 0.9, width = 1,height = 1) + #
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_40_95_c42", #kovesi.diverging_bwr_55_98_c37
                         limits = c(0,4),na.value = "#D9D9D9")+
  scale_color_manual(values = c('black',"black"),na.value =  "transparent",name = "Significance") + #Significance Genotype vs Col-0
  theme(axis.text.x = element_text(angle = -45, hjust=-0.05),axis.title = element_blank())


# 4. Plotting Relative Abundance Bar Charts####
# phylum-level
ps.compositional <- microbiome::transform(ps.perc, "compositional")
ps.phyla.perc <-taxa_level(ps.perc.F, "Rank2")
ps.clases.perc <-taxa_level(ps.compositional, "Rank3")
# identify the 10 most abundant phylum
phylum.10 <- names(sort(taxa_sums(ps.phyla.perc), TRUE)[1:12])

clases.10 <- names(sort(taxa_sums(ps.clases.perc), TRUE)[1:15])

# now we can use the list of the top phylum and subset those from the phyloseq object
ps.phylum.10 <- prune_taxa(phylum.10, ps.phyla.perc)
ps.clases.10 <- prune_taxa(clases.10, ps.clases.perc)

melt.phylum <- psmelt(ps.phylum.10)
melt.clases <- psmelt(ps.clases.10)


safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888",
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
#,"black", "black", "black", "black", "black","black", "black", "black", "black", "black","black", "black", "black", "black", "black")

melt.phylum$ID <- paste(melt.phylum$Irrigation, melt.phylum$Treatment, sep="_")
melt.clases$ID <- paste(melt.clases$Irrigation, melt.clases$Treatment, sep="_")

melt.phylum.mean <- melt.phylum%>%group_by(ID, OTU)%>%
  summarise_all(mean)
melt.clases.mean <- melt.clases%>%group_by(ID, OTU)%>%
  summarise_all(mean)

melt.phylum.mean$ID= factor(melt.phylum.mean$ID, levels=c("Dryland_NCC","Dryland_GBL_NR","Dryland_GBL_R","Dryland_G","Irrigated_NCC","Irrigated_GBL_R","Irrigated_GBL_NR","Irrigated_G"))
ggplot(melt.phylum.mean, aes(x = ID, y = Abundance, fill = OTU)) + 
  theme_few() +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = safe_colorblind_palette) + 
  labs(fill = "Phylum")+
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5))

melt.clases.mean$ID= factor(melt.clases.mean$ID, levels=c("Dryland_NCC","Dryland_GBL_NR","Dryland_GBL_R","Dryland_G","Irrigated_NCC","Irrigated_GBL_R","Irrigated_GBL_NR","Irrigated_G"))
ggplot(melt.clases.mean, aes(x = ID, y = Abundance, fill = OTU)) + 
  theme_few() +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = safe_colorblind_palette) + 
  labs(fill = "Clases")+theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5))

#5. Heatmap
ps.genus.perc.fun <-taxa_level(ps.perc, "Rank6")

# identify the 50 most abundant genus
genus.10.fun <- names(sort(taxa_sums(ps.genus.perc.fun), TRUE)[2:51])
# now we can use the list of the top phylum and subset those from the phyloseq object
ps.genus.10.fun <- prune_taxa(genus.10.fun, ps.genus.perc.fun)

Table.fun = cbind(ps.genus.10.fun@sam_data, log10(ps.genus.10.fun@otu_table+1)) 
Table_mean.fun <- Table.fun%>%group_by(Treatment, Irrigation)%>%
  summarise_all(mean)

Mean1 <- data.frame(Table_mean.fun)
Mean0_15<- data.frame(t(Mean1[,-c(1:6)]))
names (Mean0_15) = c("Dryland_G","Irrigated_G","Dryland_GBL_NR","Irrigated_GBL_NR","Dryland_GBL_R","Irrigated_GBL_R","Dryland_NCC","Irrigated_NCC")
annotation_col = data.frame(
  Genera =substr(colnames(Mean0_15),1,3))
rownames(annotation_col)=colnames(Mean0_15)
colnames(Mean0_15)
row.names(Mean1$ID)

#annotation_col$Genera = factor(annotation_col$Genera, c('DB/','BW/','RO/'))

library(pheatmap)
pheatmap(Mean0_15, show_rownames=TRUE,show_colnames=TRUE,
         annotation_col=annotation_col, annotation_legend = TRUE,
         scale = "none",clustering_method="ward.D2",
         fontsize_row = 6, #breaks = 0:4,
         clustering_distance_cols="euclidean")




##7.  Alfadiversity
ordered(sample_sums(ps.2))

# let's calcuolate Shannon diversity and Richness (Observed)
alpha.div<-estimate_richness(physeq_fun, measures=c("Shannon", "Observed",'Chao', 'Simpson','Richness'))
even <- evenness(physeq_fun, 'pielou')

# it is easiest to append the alpha-diversity estimates to the metadata file for downstream analyses

Data = data.frame(physeq_fun@sam_data)

Data$Shannon <- paste(alpha.div$Shannon)
Data$Observed <- paste(alpha.div$Observed)
Data$Simpson <- paste(alpha.div$Simpson)
Data$Chao1 <- paste(alpha.div$Chao1)
Data$Evenness <- even$pielou
Data$Observed <- as.numeric(Data$Observed)
Data$Shannon <- as.numeric(Data$Shannon)
Data$Simpson <- as.numeric(Data$Simpson)
Data$Chao1 <- as.numeric(Data$Chao1)
Data$Plant <- as.factor(Data$Plant)
Data$Treatment <- as.factor(Data$Treatment)
Data$Irrigation <- as.factor(Data$Irrigation)
Data$Rep <- as.numeric(Data$Rep)
write.table(Data, "Alfadiversity en txt.txt")

Alfa_Fun = read.table("Alfadiversity en txt.txt")

Alfa_Fun2 = Alfa_Fun[-c(2),]
##Shannon
# the below commands will perform anova's based on flavonoid production and plant genotype
# Note that I have included Block as a random factor here, thus I have to make a linear mixed effects model before running the anova
# make sure lmerTest is a loaded package or you will not see p-values
shannon.model.fun <-lm(Shannon ~ Irrigation   + (1|Rep), data = Alfa_Fun2[Alfa_Fun2$Treatment!='G',])
anova(shannon.model.fun) 
# Df Sum Sq Mean Sq F value Pr(>F)
#Irrigation  1 0.6373 0.63735  2.0263   0.17
#Residuals  20 6.2907 0.31454 

shannon.model.fun <-lm(Shannon ~ Treatment   + (1|Rep), data = Alfa_Fun2[Alfa_Fun2$Irrigation!='Dryland',])
anova(shannon.model.fun) 
#Df  Sum Sq Mean Sq F value Pr(>F)
#Treatment  3 0.47974 0.15992  0.6656 0.5904
#Residuals 11 2.64296 0.24027 

shannon.model.fun <-lm(Shannon ~ Treatment   + (1|Rep), data = Alfa_Fun2[Alfa_Fun2$Irrigation =='Dryland',])
anova(shannon.model.fun) 
#Df  Sum Sq Mean Sq F value Pr(>F)
#Treatment  3 1.2229 0.40764  1.6575 0.2331
#Residuals 11 2.7052 0.24593 


shannon.model.fun_I <-lm(Evenness ~ Treatment   + (1|Rep), data = Data[Data$Irrigation=='Irrigated',])
anova(shannon.model.fun_I) 
#NS

shannon.model.fun_D <-lm(Evenness ~ Treatment   + (1|Rep), data = Data[Data$Irrigation!='Irrigated',])
anova(shannon.model.fun_D) 
#NS

Alfa_Fun$Treatment = factor(Alfa_Fun$Treatment, c('NCC','GBL_NR','GBL_R','G'))
ggplot(Alfa_Fun, aes(Irrigation, Shannon,  fill = Treatment)) + 
  geom_boxplot() + 
  theme_bw()+
  theme(legend.position = 'none')
  

##Venn Diagram
library(gplots)
library(eulerr)

ps.perc.N <- subset_samples(ps.perc.2, Treatment == "NCC")
ps.perc.R <- subset_samples(ps.perc.2, Treatment == "GBL_R")
ps.perc.NR <- subset_samples(ps.perc.2, Treatment == "GBL_NR")
ps.perc.G <- subset_samples(ps.perc.2, Treatment == "G")

ps.perc.N1 = prune_taxa(taxa_sums(ps.perc.N) > 0, ps.perc.N)
ps.perc.R1 = prune_taxa(taxa_sums(ps.perc.R) > 0, ps.perc.R)
ps.perc.NR1 = prune_taxa(taxa_sums(ps.perc.NR) > 0, ps.perc.NR)
ps.perc.G1 = prune_taxa(taxa_sums(ps.perc.G) > 0, ps.perc.G)

N <- rownames(otu_table(ps.perc.N1))
R <- rownames(otu_table(ps.perc.R1))
NR <- rownames(otu_table(ps.perc.NR1))
G <- rownames(otu_table(ps.perc.G1))

Venn <-venn(list(A=N, B=R, C=NR, D=G))
plot(venn(list('NCC'=N, 'GBL_NR'=NR, 'GBL_R'=R, 'G'=G)))



####Phylum analysis
library(agricolae)
ps.phyla.perc
ncol(ps.phyla.perc@otu_table)
Phylum.F = data.frame(cbind(ps.phyla.perc@sam_data,ps.phyla.perc@otu_table))
colnames(Phylum.F)

Phylum.F$Treatment = factor(Phylum.F$Treatment, c('NCC','GBL_NR','GBL_R','G'))

#p__Ascomycota
summary(aov(p__Ascomycota ~ Irrigation+ (1|Rep), data = Phylum.F[Phylum.F$Treatment!='G',]))
#NS
#Irrigation
summary(aov(p__Ascomycota ~ Treatment   + (1|Rep), data = Phylum.F[Phylum.F$Irrigation!='Dryland',]))
#Df Sum Sq Mean Sq F value Pr(>F)  
#Treatment    3  507.3  169.11   5.405 0.0157 *
#Residuals   11  344.2   31.29

summary(aov(p__Ascomycota ~ Treatment   + (1|Rep), data = Phylum.F[Phylum.F$Irrigation=='Dryland',]))
#NS

ggplot(Phylum.F, aes(Irrigation, p__Ascomycota,  fill = Treatment)) + 
  geom_boxplot() + 
  theme_few()


#p__Basidiomycota
summary(aov(p__Basidiomycota ~ Irrigation+ (1|Rep), data = Phylum.F[Phylum.F$Treatment!='G',]))
#NS
summary(aov(p__Basidiomycota ~ Treatment   + (1|Rep), data = Phylum.F[Phylum.F$Irrigation!='Dryland',]))
#NS
summary(aov(p__Basidiomycota ~ Treatment   + (1|Rep), data = Phylum.F[Phylum.F$Irrigation=='Dryland',]))
#NS

ggplot(Phylum.F, aes(Irrigation, p__Basidiomycota,  fill = Treatment)) + 
  geom_boxplot() + 
  theme_few()

#p__Mortierellomycota 
summary(aov(p__Mortierellomycota ~ Irrigation+ (1|Rep), data = Phylum.F[Phylum.F$Treatment!='G',]))
#Df Sum Sq Mean Sq F value  Pr(>F)   
#Irrigation   1  56.11   56.11   13.44 0.00144 **
#Residuals   21  87.69    4.18

ggplot(Phylum.F, aes(Irrigation, p__Mortierellomycota,  fill = Treatment)) + 
  geom_boxplot() + 
  theme_few()

summary(aov(p__Mortierellomycota ~ Treatment   + (1|Rep), data = Phylum.F[Phylum.F$Irrigation!='Dryland',]))
#NS
summary(aov(p__Mortierellomycota ~ Treatment   + (1|Rep), data = Phylum.F[Phylum.F$Irrigation=='Dryland',]))
#NS


#p__Chytridiomycota 
summary(aov(p__Chytridiomycota ~ Irrigation+ (1|Rep), data = Phylum.F[Phylum.F$Treatment!='G',]))
#NS
summary(aov(p__Chytridiomycota ~ Treatment   + (1|Rep), data = Phylum.F[Phylum.F$Irrigation!='Dryland',]))
#NS
summary(aov(p__Chytridiomycota ~ Treatment   + (1|Rep), data = Phylum.F[Phylum.F$Irrigation=='Dryland',]))
#NS

ggplot(Phylum.F, aes(Irrigation, p__Chytridiomycota,  fill = Treatment)) + 
  geom_boxplot() + 
  theme_few()

colnames(Phylum)
#
#p__Glomeromycota 
summary(aov(p__Glomeromycota ~ Irrigation+ (1|Rep), data = Phylum.F[Phylum.F$Treatment!='G',]))
#  Df Sum Sq Mean Sq F value Pr(>F)  
#Irrigation   1  7.949   7.949   7.561  0.012 *
#Residuals   21 22.076   1.051  

ggplot(Phylum.F, aes(Irrigation, p__Glomeromycota,  fill = Treatment)) + 
  geom_boxplot() + 
  theme_few()

summary(aov(p__Glomeromycota ~ Treatment   + (1|Rep), data = Phylum.F[Phylum.F$Irrigation!='Dryland',]))
#NS
summary(aov(p__Glomeromycota ~ Treatment   + (1|Rep), data = Phylum.F[Phylum.F$Irrigation=='Dryland',]))
#NS

#p__Aphelidiomycota 
summary(aov(p__Aphelidiomycota ~ Irrigation+ (1|Rep), data = Phylum.F[Phylum.F$Treatment!='G',]))
#NS

summary(aov(p__Aphelidiomycota ~ Treatment   + (1|Rep), data = Phylum.F[Phylum.F$Irrigation!='Dryland',]))
#NS
summary(aov(p__Aphelidiomycota ~ Treatment   + (1|Rep), data = Phylum.F[Phylum.F$Irrigation=='Dryland',]))
#NS

ggplot(Phylum.F, aes(Irrigation, p__Aphelidiomycota,  fill = Treatment)) + 
  geom_boxplot() + 
  theme_few()

colnames(Phylum)
#p__Kickxellomycota 
summary(aov(p__Kickxellomycota ~ Irrigation+ (1|Rep), data = Phylum.F[Phylum.F$Treatment!='G',]))
#NS

summary(Kickxe_I = aov(p__Kickxellomycota ~ Treatment   + (1|Rep), data = Phylum.F[Phylum.F$Irrigation!='Dryland',]))
#Df Sum Sq Mean Sq F value Pr(>F)  
#Treatment    3  4.825  1.6084   4.571  0.026 *
#Residuals   11  3.871  0.3519 

LSDoav <- HSD.test(Kickxe_I, c("Treatment"), alpha = 0.05, group=TRUE) 
LSDoav

ggplot(Phylum.F, aes(Irrigation, p__Kickxellomycota,  fill = Treatment)) + 
  geom_boxplot() + 
  theme_few()

summary(Kickxe_D = aov(p__Kickxellomycota ~ Treatment   + (1|Rep), data = Phylum.F[Phylum.F$Irrigation=='Dryland',]))
#NS
Kickxe_D

LSDoav <- HSD.test(Kickxe_D, c("Treatment"), alpha = 0.05, group=TRUE) 
LSDoav

colnames(Phylum)
#p__Zoopagomycota 
summary(aov(p__Zoopagomycota ~ Irrigation+ (1|Rep), data = Phylum.F[Phylum.F$Treatment!='G',]))
#NS

summary(aov(p__Zoopagomycota ~ Treatment   + (1|Rep), data = Phylum.F[Phylum.F$Irrigation!='Dryland',]))
#NS

summary(aov(p__Zoopagomycota ~ Treatment   + (1|Rep), data = Phylum.F[Phylum.F$Irrigation=='Dryland',]))
#NS          p__Zoopagomycota groups

#p__Basidiobolomycota 
summary(aov(p__Basidiobolomycota ~ Irrigation+ (1|Rep), data = Phylum.F[Phylum.F$Treatment!='G',]))
#NS

summary(aov(p__Basidiobolomycota ~ Treatment   + (1|Rep), data = Phylum.F[Phylum.F$Irrigation!='Dryland',]))
#NS

summary(aov(p__Basidiobolomycota ~ Treatment   + (1|Rep), data = Phylum.F[Phylum.F$Irrigation=='Dryland',]))


colnames(Phylum)
#p__Mucoromycota 
summary(aov(p__Mucoromycota ~ Irrigation+ (1|Rep), data = Phylum.F[Phylum.F$Treatment!='G',]))
#NS

summary(aov(p__Mucoromycota ~ Treatment   + (1|Rep), data = Phylum.F[Phylum.F$Irrigation!='Dryland',]))
#NS

summary(aov(p__Mucoromycota ~ Treatment   + (1|Rep), data = Phylum.F[Phylum.F$Irrigation=='Dryland',]))

Phylum_mean <- Phylum%>%group_by(Plot)%>%
  summarise_all(mean)
Phylum_mean2= data.frame(Phylum_mean)
#p__Neocallimastigomycota 4.100428e-06
#p__Calcarisporiellomycota 8.202836e-07
#p__unclassified 0.05363228
#p__Monoblepharomycota: 4.086703e-06 
#p__Entomophthoromycota: 7.306088e-05
#p__Mucoromycota 0.0002253888
#p__Basidiobolomycota 0.0002592296
#p__Zoopagomycota 0.001770157
#p__Kickxellomycota 0.006640763
#p__Aphelidiomycota:0.000354498
#p__Glomeromycota 0.008586809
#p__Chytridiomycota:0.008977863 
#p__Mortierellomycota:0.06650007 
#k__Fungi_unclassified: 0.1309838 
#p__Basidiomycota: 0.2044166
#p__Ascomycota:0.5175705

#### log2foldchange
#Root and Front vs water
library(reshape2)
ps.perc
tableSyncom = cbind(ps.perc.F.2@sam_data, t(ps.perc.F.2@otu_table))
ncol(tableSyncom)
tableSyncom2=tableSyncom[,-c(1,4:6)]

melted_tableSyncom <- tableSyncom2 %>% melt

ttestRat <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value
}

melted_tableSyncomG = melted_tableSyncom[melted_tableSyncom$Treatment != 'G',]

###Irrigation
log2foldchange <- c()

for(plant in melted_tableSyncomG$Treatment  %>% unique){
  melted_sub <- melted_tableSyncomG %>% subset(Treatment ==  plant) %>% droplevels
  for(specie in melted_sub$variable  %>% unique){
  melted_sub2 <- melted_sub %>% subset(variable ==  specie) %>% droplevels
  
  Irrigated = data.frame(t(melted_sub2[melted_sub2$Irrigation=='Irrigated',][,4]))
  Dryland = data.frame(t(melted_sub2[melted_sub2$Irrigation=='Dryland',][,4]))
  
  Irrigated_mean = apply(Irrigated, 1, mean) 
  Drylnad_mean = apply(Dryland, 1, mean) 
  
  Irri_Dry <- log2(Irrigated_mean+1) - log2(Drylnad_mean+1) 
  
  Irri_Dry_statistic = t.test(Irrigated,Dryland)
  pvalueIrri_Dry = Irri_Dry_statistic$p.value#Water_mean+1
  
  result = cbind(Irri_Dry,pvalueIrri_Dry)
  result2 = cbind(result,plant)
  result3 = cbind(result2,specie)
  log2foldchange <- rbind(log2foldchange,result3)
}}

colnames(log2foldchange)=c('diff','pvalue','Treatment','ASV')
log2foldchange2 = data.frame(log2foldchange)

####Adjust the p-value
log2foldchange2$Significance <- "No Significant"
pval_thres <- 0.1
log2foldchange2$Significance[which(log2foldchange2$pvalue < pval_thres)] <- "q < 0.05"
log2foldchange2$Significance <- log2foldchange2$Significance %>% factor

log2foldchange2_selected_2 =log2foldchange2[log2foldchange2$Significance != "No Significant",] 

write.table(log2foldchange2, "log2foldchange2 fungi water.txt")
write.table(log2foldchange2_selected_2, "log2foldchange2 fungi water seleccionados.txt")

head(log2foldchange2_selected)

Log_Irrigation = log2foldchange2_selected_2

#Treatment
log2foldchange <- c()

for(water in melted_tableSyncom$Irrigation  %>% unique){
  melted_sub <- melted_tableSyncom %>% subset(Irrigation ==  water) %>% droplevels
  for(specie in melted_sub$variable  %>% unique){
    melted_sub2 <- melted_sub %>% subset(variable ==  specie) %>% droplevels
    for(plant in melted_sub2[melted_sub2$Treatment != 'NCC',]$Treatment  %>% unique){
    melted_sub3 <- melted_sub2[melted_sub2$Treatment != 'NCC',] %>% subset(Treatment ==  plant) %>% droplevels
    
    NCC = data.frame(t(melted_sub2[melted_sub2$Treatment=='NCC',][,4]))
    CC = data.frame(t(melted_sub3[,4]))
    
    NCC_mean = apply(NCC, 1, mean) 
    CC_mean = apply(CC, 1, mean) 
    
    NCC_CC <- log2(NCC_mean+1) - log2(CC_mean+1) 
    
    NCC_CC_statistic = t.test(NCC, CC)
    pvalueNCC_CC = NCC_CC_statistic$p.value#Water_mean+1
    
    result = cbind(NCC_CC,pvalueNCC_CC)
    result2 = cbind(result,plant)
    result3 = cbind(result2,specie)
    result4 = cbind(result3,water)
    log2foldchange <- rbind(log2foldchange,result4)
  }}}

colnames(log2foldchange)=c('diff','pvalue','Treatment','ASV', 'Irrigation')
log2foldchange2 = data.frame(log2foldchange)

####Adjust the p-value
log2foldchange2$Significance <- "No Significant"
pval_thres <- 0.1
log2foldchange2$Significance[which(log2foldchange2$pvalue < pval_thres)] <- "q < 0.05"
log2foldchange2$Significance <- log2foldchange2$Significance %>% factor

log2foldchange2_selected_TREATMENT =log2foldchange2[log2foldchange2$Significance != "No Significant",] 
write.table(log2foldchange2_selected_TREATMENT, "log2foldchange2 fungi Treatment seleccionados.txt")


guild_fg_table <- read.table("Functional Trait Fungi core.txt", header = TRUE)
melted_fungal_traits <- guild_fg_table[,-c(1,4:6)] %>% melt


fungal_traits_mean <- guild_fg_table %>%group_by(Irrigation, Treatment)%>%
  summarise_all(mean)

melted_fungal_traits_Norpa_mean2 <- fungal_traits_mean[,-(3:6)] %>% melt
melted_fungal_traits_Norpa_mean2$ID <- paste(melted_fungal_traits_Norpa_mean2$Irrigation, melted_fungal_traits_Norpa_mean2$Treatment, sep="_")


safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888",
                             "#619333","#CC4677", "#DDCC79", "#149733", "#332850", "#AA4469",
                             "#94AA99","#87CCFE", "#882215")

melted_fungal_traits_Norpa_mean2$Treatment = factor(melted_fungal_traits_Norpa_mean2$Treatment, c('NCC', 'GBL_NR', 'GBL_R','G'))

ggplot(data=melted_fungal_traits_Norpa_mean2, aes(x=Treatment, y=value, fill=variable)) + 
  geom_bar(stat="identity")+ 
  theme_few()+ 
  facet_grid(.~Irrigation,space = "fixed",scales = "fixed")+
  ylab("Relative abundance (%)")+
  scale_fill_manual(values=safe_colorblind_palette)+
  #theme(legend.position = 'none')+
  ylim(0,50)

melted_fungal_traits$Treatment = factor(melted_fungal_traits$Treatment, c('NCC', 'GBL_NR', 'GBL_R','G'))
melted_fungal_traits2 = cbind(melted_fungal_traits,scale(melted_fungal_traits[4]))


colnames(melted_fungal_traits2) = c('Irrigation','Treatment','variable','value','log_value')

melted_fungal_traits2$Significance <- "No Significant"

ggplot(data = melted_fungal_traits2, aes(Treatment,variable)) +
  geom_raster(aes(fill = log_value))+
  theme_few() +
  #facet_grid(.~Irrigation)+
  #geom_text(aes(label = melterd2), color = "black", size = 1.5) +
  geom_tile(aes(color = Significance),fill = '#00000000', linewidth = 0.2,width = 1,height = 1) + #
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-2,4),name = "Corr") + #,na.value = "#D9D9D9"
  scale_color_manual(values = c("black",'grey'),na.value =  "transparent",name = "Significative") +
  labs(y='Microbial paramenters',x='Chemical parameters')+
  theme(axis.text.x = element_text(angle = -45, hjust=-0.05), #element_blank()
        axis.title = element_blank())


log2foldchange <- c()

for(specie in melt.function[melt.function$Treatment!='G',]$variable  %>% unique){
    melted_sub2 <- melt.function[melt.function$Treatment!='G',] %>% subset(variable ==  specie) %>% droplevels
    
    Irrigated = data.frame(t(melted_sub2[melted_sub2$Irrigation=='Irrigated',][,7]))
    Dryland = data.frame(t(melted_sub2[melted_sub2$Irrigation=='Dryland',][,7]))
    
    Irrigated_mean = apply(Irrigated, 1, mean) 
    Drylnad_mean = apply(Dryland, 1, mean) 
    
    Irri_Dry <- log2(Irrigated_mean+1) - log2(Drylnad_mean+1) 
    
    Irri_Dry_statistic = t.test(Irrigated,Dryland)
    pvalueIrri_Dry = Irri_Dry_statistic$p.value#Water_mean+1
    
    result = cbind(Irri_Dry,pvalueIrri_Dry)
    result2 = cbind(result,specie)
    log2foldchange <- rbind(log2foldchange,result2)
  }

colnames(log2foldchange)=c('diff','pvalue','Functional')
log2foldchange2 = data.frame(log2foldchange)

####Adjust the p-value
log2foldchange2$Significance <- "No Significant"
pval_thres <- 0.05
log2foldchange2$Significance[which(log2foldchange2$pvalue < pval_thres)] <- "q < 0.05"
log2foldchange2$Significance <- log2foldchange2$Significance %>% factor

log2foldchange2_selected =log2foldchange2[log2foldchange2$Significance == "q < 0.05",] 


##############RDA###########
par(mfrow=c(1,1))
Fun <- rda(t(ps.perc.F.2@otu_table)~., chemistry[c(33:61,63,64),-(1:5)], scale=T)
summary(Fun)
Fun$CA

f<-data.frame(rownames(Fun$CCA$u),Fun$CCA$u[, 1:2])
colnames(f) =c('ID', 'RDA1','RDA2')
sam = cbind(rownames(ps.perc.F.2@sam_data), ps.perc.F.2@sam_data)
colnames(sam) =c('ID', 'Plot','Irrigation','Treatment','Plants','Residue','Rep')
f <- merge(f,sam, by='ID')
rownames(f)<-r$ID

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499") 

f$Treatment = factor(f$Treatment, c('NCC','GBL_NR','GBL_R','G'))

load.F = data.frame(Fun$CCA$biplot[,1:2])


RDAfUN = ggplot(f, aes(RDA1, RDA2)) + 
  geom_vline(xintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_hline(yintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_point(aes(colour=Treatment, shape=Irrigation),size=5) +
  theme_bw() +
  theme(legend.position="right")+ 
  labs(y=ylab_text, x=xlab_text)+
  scale_colour_manual(values = safe_colorblind_palette)+
  geom_segment(data=load.F,aes(x=0,y=0,xend=RDA1,yend=RDA2),arrow=arrow(length=unit(.3,"lines")),color="black",size=0.8)+
  annotate("text",x=load.F$RDA1, y=load.F$RDA2,label=rownames(load.F),size=4)
  xlim(-0.3,0.9)

anova (Fun, by = 'margin', parallel = 8)
#         Df Variance      F Pr(>F)  
#NT0       1   165.47 1.2403  0.040 *
#PMN       1   164.05 1.2297  0.030 *


###Phylum vs chemistry
ps.phylum.10@otu_table
tableSyncom = cbind(ps.phylum.10@sam_data, ps.phylum.10@otu_table)
tableSyncom2 = cbind(rownames(tableSyncom), tableSyncom)
colnames(tableSyncom2) = c('ID',colnames(tableSyncom))
tableSyncom2$ID = as.numeric(tableSyncom2$ID)
tableSyncom3 = tableSyncom2[order(tableSyncom2$ID),]

sam = cbind(rownames(chemistry), chemistry)
colnames(sam) =c('ID', colnames(chemistry))
tableSyncom3$ID = as.factor(tableSyncom3$ID)

Fun_che <- data.frame(cbind(tableSyncom3,sam[,-c(1)]))

Corr.F <-cor(Fun_che[,-c(1:7)], method = 'spearman')
p.mat.F <-cor_pmat(Fun_che[,-c(1:7)], method = 'spearman')


Corr.F2 = data.frame(Corr.F[c(2:9,11,12),-c(1:12)])
p.mat.F2 = data.frame(p.mat.F[c(2:9,11,12),-c(1:12)])

Corr.F3 = data.frame(Corr.F[c(1:12),-c(1:12)])
p.mat.F3 = data.frame(p.mat.F[c(1:12),-c(1:12)])


ggcorrplot(Corr.F3,
           p.mat = p.mat.F3, insig = "blank",
           method = "circle",
           hc.order = TRUE, #type = "upper",
           outline.color = "white",
           ggtheme = ggplot2::theme_bw,
           colors = c("#6D9EC1", "white", "#E46726"))



Alf = cbind(rownames(Alfadiversity), Alfadiversity)
colnames(Alf) =c('ID', colnames(Alfadiversity))
Alf$ID = as.factor(Alf$ID)

Bac_Alf_quim <- data.frame(cbind(sam[33:64,], Alf[,c(1,8)]))
ncol(Bac_Alf_quim)

Corr.F <-cor(Bac_Alf_quim[,-c(1:5,19)], method = 'spearman')
p.mat.F <-cor_pmat(Bac_Alf_quim[,-c(1:5,19)], method = 'spearman')

ggcorrplot(Corr.F,
           p.mat = p.mat.F, insig = "blank",
           method = "circle",
           hc.order = TRUE, #type = "upper",
           outline.color = "white",
           ggtheme = ggplot2::theme_bw,
           colors = c("#6D9EC1", "white", "#E46726"))

