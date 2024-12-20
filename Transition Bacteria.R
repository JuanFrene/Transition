packages <- c('ggthemes','dplyr', "ape", "ShortRead", "Biostrings", "phyloseq", 
              "corncob","DESeq2", "microbiome", "DECIPHER", "phangorn", "tibble", "lme4", 
              "lmerTest", "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire")
sapply(packages, require, character.only = TRUE)              

# Set working directory
setwd("G:/My Drive/labs/NMSU/Transtion/Bacteria/")#CC Network/

load("G:/My Drive/labs/NMSU/Transtion/Bacteria/transition bacteria.RData")

# Set plotting theme
theme_set(theme_bw())
# Assign variables for imported data
#sharedfile_bac = "final.asv.shared"
#taxfile_bac = "final.asv.ASV.cons.taxonomy"

#mothur_data_bac <- import_mothur(mothur_shared_file = sharedfile_bac,
#                             mothur_constaxonomy_file = taxfile_bac)
#mothur_data_bac@otu_table[5:5]

#Agregar la data
#Metadata <- read.table("Metadata.txt", header = TRUE)
#sample_dat = sample_data(Metadata)
#sample_names(sample_dat) = rownames(Metadata)

## Building the Phyloseq Object
physeq_bac = merge_phyloseq(mothur_data_bac, sample_dat)
ps.1 = prune_taxa(taxa_sums(physeq_bac) >= 5, physeq_bac)
sample_sums(ps.1)

# Now, let's make sure that the correct information is included in our phyloseq object
ps.2 = subset_taxa(ps.1, Rank1 ==  "Bacteria" |Rank1 ==  "Archaea") #Kingdom == "Bacteria" | Kingdom == "Archaea"
ps.perc <- transform_sample_counts(ps.2, function(x) x / sum(x) * 100) 
ps.perc.2 = prune_taxa(taxa_sums(ps.perc) >= 0.1, ps.perc)
ps.2.H <- transform(ps.2, "hellinger")


tax = data.frame(ps.perc.2@tax_table)
head(tax)
tax_bacteroidota = tax[tax$Rank2 =='Bacteroidota',]

ps.2@tax_table[1:25,]


###ANALYSIS###
# 1. PCoA plot####
paleta_alive <- c('#FF0000','#00008B',"#FFB919","#00CC1C")

plot_ordination(ps.perc.2, ordinate(ps.perc.2, "NMDS", "bray"), color = "Treatment", shape = 'Irrigation') + theme_few() +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  geom_point(size = 4) +
  scale_color_manual(values = safe_colorblind_palette)

Kdist <- dist(t(ps.perc.2@otu_table)) 
pc<-ape::pcoa(Kdist)

PCoA <- ordinate(ps.2.H, "PCoA", "bray")
title ="PCoA - PCoA1 vs PCoA2"
ev1 <- (pc$values$Relative_eig[1])*100
ev2 <- (pc$values$Relative_eig[2])*100
xlab_text <- paste("PCoA1 (", round(ev1,2), "%)")
ylab_text <- paste("PCoA2 (", round(ev2,2), "%)")

x<-data.frame(pc$vectors[, 1:2])
x <- merge(x,ps.2@sam_data, by=0)
rownames(x)<-x$Row.names
x$Row.names <- NULL

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499") 

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
  scale_colour_manual(values = safe_colorblind_palette)+
  ylim(-1,1)


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
ps.T3_bray <- phyloseq::distance(ps.2.H, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(ps.2))

# Adonis test
PERMANOVA = adonis2(ps.T3_bray ~ Treatment*Irrigation, data = sampledf)
#Df SumOfSqs      R2      F Pr(>F)   
#Treatment             3   1.4557 0.10380 1.0462  0.002 **
#Irrigation            1   0.5055 0.03604 1.0899  0.002 **
#Treatment:Irrigation  3   1.3952 0.09949 1.0028  0.404   
#Residual             23  10.6670 0.76066                 
#Total                30  14.0233 1.00000  

PERMANOVA2 = data.frame(PERMANOVA)
variable  = c('variable','variable','variable','variable','variable')
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
data_melt <- LC[-(5),]

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

ggsave("Bacteria_variance.pdf", height=5, width=3.5, units='in')



# 3. Simper function
otudat <- t(otu_table(ps.perc.2))
otudat <- otudat[rowSums(otudat) > 0,]
ncol(otudat)
## in case you also have empty otu
otudat <- otudat[, colSums(otudat) > 0]
Simper_IB = simper(otudat, sample_data(ps.perc.2)$Irrigation)
Simper_TB = simper(otudat, sample_data(ps.perc.2)$Treatment)


table_I_genus = summary(Simper_IB)$
table_T_genus_NCC_GBL_R = summary(Simper_TB)$NCC_GBL_R
table_T_genus_NCC_G = summary(Simper_TB)$NCC_G
table_T_genus_NCC_GBL_NR = summary(Simper_TB)$NCC_GBL_NR
table_T_genus_GBL_R_G = summary(Simper_TB)$GBL_R_G
table_T_genus_GBL_R_GBL_NR = summary(Simper_TB)$GBL_R_GBL_NR
table_T_genus_G_GBL_NR = summary(Simper_TB)$G_GBL_NR

only_I_sig = table_I[table_I$p == '0.001',]
only_I_sig_r = only_I_sig[only_I_sig$ratio > '1.5',]

names = rownames(summary(Simper_TB)$GBL_R_GBL_NR)
table = cbind(names, summary(Simper_TB)$GBL_R_GBL_NR)
table[table$names == 'ASV0000022',]

# 4. Plotting Relative Abundance Bar Charts####
# phylum-level
ps.compositional <- microbiome::transform(ps.perc, "compositional")
ps.phyla.perc <-taxa_level(ps.perc.2, "Rank2")
ps.clases.perc <-taxa_level(ps.compositional, "Rank3")
ps.genus.perc <-taxa_level(ps.compositional, "Rank6")

# identify the 10 most abundant phylum
phylum.10 <- names(sort(taxa_sums(ps.phyla.perc), TRUE)[1:12])

clases.10 <- names(sort(taxa_sums(ps.clases.perc), TRUE)[1:15])

# now we can use the list of the top phylum and subset those from the phyloseq object
phylum.10 <- names(sort(taxa_sums(ps.phyla.perc), TRUE)[1:12])
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
ps.genus.perc.bac <-taxa_level(ps.perc, "Rank6")

# identify the 50 most abundant genus
genus.10.bac <- names(sort(taxa_sums(ps.genus.perc.bac), TRUE)[4:53])
# now we can use the list of the top phylum and subset those from the phyloseq object
ps.genus.10.bac <- prune_taxa(genus.10.bac, ps.genus.perc.bac)

Table.bac = cbind(ps.genus.10.bac@sam_data, log10(ps.genus.10.bac@otu_table+1)) 
Table_mean.bac <- Table.bac%>%group_by(Treatment, Irrigation)%>%
  summarise_all(mean)

Mean1 <- data.frame(Table_mean.bac)
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
alpha.div<-estimate_richness(ps.2, measures=c("Shannon", "Observed",'Chao', 'Simpson','Richness'))
even <- evenness(physeq_bac, 'pielou')
write.table(Data, "Alfadiversity en txt.txt")
Alfadiversity = read.table("Alfadiversity en txt.txt")
# it is easiest to append the alpha-diversity estimates to the metadata file for downstream analyses

Data = data.frame(physeq_bac@sam_data)

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

Alfadiversity2 = Alfadiversity[-c(13,17),]
##Shannon
# the below commands will perform anova's based on flavonoid production and plant genotype
# Note that I have included Block as a random factor here, thus I have to make a linear mixed effects model before running the anova
# make sure lmerTest is a loaded package or you will not see p-values

Observed.model.bac <-lm(Observed ~ Irrigation   + (1|Rep), data = Alfadiversity[Alfadiversity$Treatment!='G',])
anova(Observed.model.bac) 
#Df     Sum Sq   Mean Sq F value Pr(>F)
#Irrigation  1  174285651 174285651  1.8406 0.1886
#Residuals  22 2083183708  94690169

Observed.model.bac_I <-lm(Observed ~ Treatment   + (1|Rep), data = Alfadiversity[Alfadiversity$Irrigation!='Dryland',])
anova(Observed.model.bac_I) 
#Df     Sum Sq   Mean Sq F value Pr(>F)
#Irrigation  1   57411104  19137035  0.1579 0.9226
#Residuals 12 1454694278 121224523

Observed.model.bac_D <-lm(Observed ~ Treatment   + (1|Rep), data = Alfadiversity2[Alfadiversity2$Irrigation=='Dryland',])
anova(Observed.model.bac_D) 
#Df     Sum Sq   Mean Sq F value Pr(>F)
#Irrigation  1   645995074 215331691  3.4621 0.05106 .
#Residuals 12 746372675  62197723 

Shannon.model.bac <-lm(Shannon ~ Irrigation   + (1|Rep), data = Alfadiversity[Alfadiversity$Treatment!='G',])
anova(Shannon.model.bac) 
#Df     Sum Sq   Mean Sq F value Pr(>F)
#Irrigation  1  0.07909 0.079088  1.3485 0.2592
#Residuals  20 1.17294 0.058647

Shannon.model.bac_I <-lm(Shannon ~ Treatment   + (1|Rep), data = Alfadiversity2[Alfadiversity2$Irrigation!='Dryland',])
anova(Shannon.model.bac_I) 
#Df     Sum Sq   Mean Sq F value Pr(>F)
#Irrigation  1    0.18955 0.063185   1.066 0.4065
#Residuals 10 0.59275 0.059275

Shannon.model.bac_D <-lm(Shannon ~ Treatment   + (1|Rep), data = Alfadiversity[Alfadiversity$Irrigation=='Dryland',])
anova(Shannon.model.bac_D) 
#Df     Sum Sq   Mean Sq F value Pr(>F)
#Irrigation  1   645995074 215331691  3.4621 0.05106 .
#Residuals 12 746372675  62197723 

Alfadiversity2$Treatment = factor(Alfadiversity2$Treatment, c('NCC','GBL_NR','GBL_R','G'))
ggplot(Alfadiversity, aes(Irrigation, Observed,  fill = Treatment)) + 
  geom_boxplot() + 
  theme_few()
  
ggplot(Alfadiversity2, aes(Irrigation, Shannon,  fill = Treatment)) + 
  geom_boxplot() + 
  theme_bw()+
  theme(legend.position = 'none')

ggplot(Data2, aes(Irrigation, Chao1,  fill = Treatment)) + 
  geom_boxplot() + 
  theme_few()


##Venn Diagram
library(gplots)
library(euler)

ps.perc.N <- subset_samples(ps.2, Treatment == "NCC")
ps.perc.R <- subset_samples(ps.2, Treatment == "GBL_R")
ps.perc.NR <- subset_samples(ps.2, Treatment == "GBL_NR")
ps.perc.G <- subset_samples(ps.2, Treatment == "G")

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
ps.phyla.perc
ncol(ps.phyla.perc@otu_table)
Phylum = data.frame(cbind(ps.phyla.perc@sam_data,ps.phyla.perc@otu_table))
colnames(Phylum)
Phylum$Treatment = factor(Phylum$Treatment, c('NCC','GBL_NR','GBL_R','G'))


#Proteobacteria
summary(aov(Proteobacteria ~ Irrigation   + (1|Rep), data = Phylum[Phylum$Treatment!='G',]))
#NS
summary(aov(Proteobacteria ~ Treatment   + (1|Rep), data = Phylum[Phylum$Irrigation!='Dryland',]))
#NS
summary(aov(Proteobacteria ~ Treatment   + (1|Rep), data = Phylum[Phylum$Irrigation=='Dryland',]))
#NS

ggplot(Phylum, aes(Irrigation, Proteobacteria,  fill = Treatment)) + 
  geom_boxplot() + 
  theme_few()

#Actinobacteriota
summary(aov(Actinobacteriota ~ Irrigation   + (1|Rep), data = Phylum[Phylum$Treatment!='G',]))
#NS
summary(actino_I = aov(Actinobacteriota ~ Treatment   + (1|Rep), data = Phylum[Phylum$Irrigation!='Dryland',]))
#NS

LSDoav <- HSD.test(actino_I, c("Treatment"), alpha = 0.05, group=TRUE) 
LSDoav

#Dryland
summary(actino_D = aov(Actinobacteriota ~ Treatment   + (1|Rep), data = Phylum[Phylum$Irrigation=='Dryland',]))
#                     Df  Sum Sq  Mean Sq F value Pr(>F)  
#Treatment             3 0.019999 0.006666   8.381 0.00284 **
#Residuals   12 0.009545 0.000795           


LSDoav <- HSD.test(actino_D, c("Treatment"), alpha = 0.05, group=TRUE) 
LSDoav

#Actinobacteriota groups
#GBL_R         0.3724595      a
#G             0.3566673      a
#GBL_NR        0.3538103      a
#NCC           0.2809956      b


ggplot(Phylum, aes(Irrigation, Actinobacteriota,  fill = Treatment)) + 
  geom_boxplot() + 
  theme_few()

#Acidobacteriota 
summary(aov(Acidobacteriota ~ Irrigation   + (1|Rep), data = Phylum[Phylum$Treatment!='G',]))
#NS

ggplot(Phylum, aes(Irrigation, Acidobacteriota,  fill = Treatment)) + 
  geom_boxplot() + 
  theme_few()

summary(acido_I = aov(Acidobacteriota ~ Treatment   + (1|Rep), data = Phylum[Phylum$Irrigation!='Dryland',]))
# Df   Sum Sq Mean Sq F value  Pr(>F)   
#Treatment    3 0.013711 0.00457   6.095 0.00922 **
#Residuals   12 0.008999 0.00075

LSDoav <- HSD.test(acido_I, c("Treatment"), alpha = 0.05, group=TRUE) 
LSDoav
#Acidobacteriota groups
#NCC         0.16410674      a
#GBL_NR      0.10716231     ab
#GBL_R       0.09654901      b
#G           0.09019068      b

summary(acido_D = aov(Acidobacteriota ~ Treatment   + (1|Rep), data = Phylum[Phylum$Irrigation=='Dryland',]))
#             Df   Sum Sq  Mean Sq F value  Pr(>F)   
#Treatment    3 0.011349 0.003783   7.379 0.00462 **
#Residuals   12 0.006152 0.000513  

LSDoav <- HSD.test(acido_D, c("Treatment"), alpha = 0.05, group=TRUE) 
LSDoav
#Acidobacteriota groups
#NCC         0.13270126      a
#G           0.07949514      b
#GBL_NR      0.07145318      b
#GBL_R       0.06576203      b

#Bacteroidota  
summary(aov(Bacteroidota ~ Irrigation   + (1|Rep), data = Phylum[Phylum$Treatment!='G',]))
#Df   Sum Sq  Mean Sq F value Pr(>F)  
#Irrigation   1 0.005518 0.005518   6.286 0.0201 *
#Residuals   22 0.019312 0.000878

ggplot(Phylum, aes(Irrigation, Bacteroidota,  fill = Treatment)) + 
  geom_boxplot() + 
  theme_few()

summary(bacte_I = aov(Bacteroidota ~ Treatment   + (1|Rep), data = Phylum[Phylum$Irrigation!='Dryland',]))
# Df   Sum Sq Mean Sq F value  Pr(>F)   
#Treatment    3 0.002728 0.0009094   4.796 0.0202 *
#Residuals   12 0.002275 0.0001896 

LSDoav <- HSD.test(bacte_I, c("Treatment"), alpha = 0.05, group=TRUE) 
LSDoav

#Bacteroidota groups
#G        0.07031161      a
#GBL_R    0.05014057     ab
#GBL_NR   0.03898380      b
#NCC      0.03774725      b

summary(aov(Bacteroidota ~ Treatment   + (1|Rep), data = Phylum[Phylum$Irrigation=='Dryland',]))
#NS

#Chloroflexi   
summary(aov(Chloroflexi ~ Irrigation   + (1|Rep), data = Phylum[Phylum$Treatment!='G',]))
#NS

summary(aov(Chloroflexi ~ Treatment   + (1|Rep), data = Phylum[Phylum$Irrigation!='Dryland',]))
# NS

summary(aov(Chloroflexi ~ Treatment   + (1|Rep), data = Phylum[Phylum$Irrigation=='Dryland',]))
#NS

ggplot(Phylum, aes(Irrigation, Chloroflexi,  fill = Treatment)) + 
  geom_boxplot() + 
  theme_few()

#Verrucomicrobiota    
summary(aov(Verrucomicrobiota ~ Irrigation   + (1|Rep), data = Phylum[Phylum$Treatment!='G',]))
#NS

summary(aov(Verrucomicrobiota ~ Treatment   + (1|Rep), data = Phylum[Phylum$Irrigation!='Dryland',]))
# NS

summary(aov(Verrucomicrobiota ~ Treatment   + (1|Rep), data = Phylum[Phylum$Irrigation=='Dryland',]))
#NS

ggplot(Phylum, aes(Irrigation, Verrucomicrobiota,  fill = Treatment)) + 
  geom_boxplot() + 
  theme_few()

#Patescibacteria    
summary(aov(Patescibacteria ~ Irrigation   + (1|Rep), data = Phylum[Phylum$Treatment!='G',]))
#NS

summary(aov(Patescibacteria ~ Treatment   + (1|Rep), data = Phylum[Phylum$Irrigation!='Dryland',]))
# NS

summary(aov(Patescibacteria ~ Treatment   + (1|Rep), data = Phylum[Phylum$Irrigation=='Dryland',]))
#      Df   Sum Sq   Mean Sq F value  Pr(>F)   
#Treatment    3 0.002868 0.0009559   7.274 0.00488 **
#Residuals   12 0.001577 0.0001314 


ggplot(Phylum, aes(Irrigation, Patescibacteria,  fill = Treatment)) + 
  geom_boxplot() + 
  theme_few()

#Gemmatimonadota     
summary(aov(Gemmatimonadota ~ Irrigation   + (1|Rep), data = Phylum[Phylum$Treatment!='G',]))
# Df    Sum Sq   Mean Sq F value Pr(>F)  
#Irrigation   1 0.0003593 0.0003593   4.623 0.0428 *
#Residuals   22 0.0017100 0.0000777

ggplot(Phylum, aes(Irrigation, Gemmatimonadota,  fill = Treatment)) + 
  geom_boxplot() + 
  theme_few()

summary(aov(Gemmatimonadota ~ Treatment   + (1|Rep), data = Phylum[Phylum$Irrigation!='Dryland',]))
# NS

summary(aov(Gemmatimonadota ~ Treatment   + (1|Rep), data = Phylum[Phylum$Irrigation=='Dryland',]))
#           Df    Sum Sq   Mean Sq F value Pr(>F)  
#Treatment    3 0.0005697 1.899e-04    5.92 0.0102 *
#Residuals   12 0.0003849 3.208e-05 


#Myxococcota      
summary(aov(Myxococcota ~ Irrigation   + (1|Rep), data = Phylum[Phylum$Treatment!='G',]))
# Df    Sum Sq   Mean Sq F value  Pr(>F)   
#Irrigation   1 0.0003041 3.041e-04   12.82 0.00167 **
#Residuals   22 0.0005219 2.372e-05  
ggplot(Phylum, aes(Irrigation, Myxococcota,  fill = Treatment)) + 
  geom_boxplot() + 
  theme_few()

summary(aov(Myxococcota ~ Treatment   + (1|Rep), data = Phylum[Phylum$Irrigation!='Dryland',]))
# NS

summary(aov(Myxococcota ~ Treatment   + (1|Rep), data = Phylum[Phylum$Irrigation=='Dryland',]))
# Df    Sum Sq   Mean Sq F value Pr(>F)  
#Treatment    3 0.0002016 6.721e-05   4.245 0.0292 *
#Residuals   12 0.0001900 1.584e-05 

ggplot(Phylum, aes(Irrigation, Myxococcota,  fill = Treatment)) + 
  geom_boxplot() + 
  theme_few()


Phylum_mean <- Phylum%>%group_by(Plot)%>%
  summarise_all(mean)
Phylum_mean2= data.frame(Phylum_mean)
#Proteobacteria 0.1941714
#Acidobacteriota 0.1009275
#Bacteroidota 0.06053032       
#Verrucomicrobiota 0.05307449  
#Planctomycetota 0.06366277  
#Actinobacteriota 0.3214369
#Chloroflexi 0.09592767 
#Gemmatimonadota 0.029 
#Patescibacteria 0.02906979   
#Gemmatimonadota 0.02236484      

#### log2foldchange
#Root and Front vs water
library(reshape2)
ps.perc
tableSyncom = cbind(ps.perc.2@sam_data, t(ps.perc.2@otu_table))
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


log2foldchange <- c()

for(plant in melted_tableSyncom$Treatment  %>% unique){
  melted_sub <- melted_tableSyncom %>% subset(Treatment ==  plant) %>% droplevels
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
log2foldchange2$Significance[which(log2foldchange2$pvalue < pval_thres)] <- "q < 0.1"
log2foldchange2$Significance <- log2foldchange2$Significance %>% factor

log2foldchange2_selected =log2foldchange2[log2foldchange2$Significance == "q < 0.1",] 


##############RDA###########
par(mfrow=c(1,1))
Bact <- rda(t(ps.perc.2@otu_table)~., chemistry[33:64,-(1:5)], scale=T)
summary(Bact)
Bact$CCA$u

r<-data.frame(rownames(Bact$CCA$u),Bact$CCA$u[, 1:2])
colnames(r) =c('ID', 'RDA1','RDA2')
sam = cbind(rownames(ps.2@sam_data), ps.2@sam_data)
colnames(sam) =c('ID', 'Plot','Irrigation','Treatment','Plants','Residue','Rep')
r <- merge(r,sam, by='ID')
rownames(r)<-r$ID

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499") 

r$Treatment = factor(r$Treatment, c('NCC','GBL_NR','GBL_R','G'))

load = data.frame(Bact$CCA$biplot[,1:2])


RDABact = ggplot(r, aes(RDA1, RDA2)) + 
  geom_vline(xintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_hline(yintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_point(aes(colour=Treatment, shape=Irrigation),size=5) +
  theme_bw() +
  theme(legend.position="right")+ 
  labs(x='RDA1(13.82%)', y='RDA (11.40%)')+
  scale_colour_manual(values = safe_colorblind_palette)+
  geom_segment(data=load,aes(x=0,y=0,xend=RDA1,yend=RDA2),arrow=arrow(length=unit(.3,"lines")),color="black",linewidth=0.8)+
  annotate("text",x=load$RDA1, y=load$RDA2,label=rownames(load),size=4)+
  xlim(-0.3,0.9)
  
anova (Bact, by = 'margin', parallel = 8)
#PMC       1    203.5 1.1827  0.036 *


###Phylum vs chemistry
ps.phyla.perc.B <-taxa_level(ps.perc.2, "Rank2")
phylum.10.B <- names(sort(taxa_sums(ps.phyla.perc.B), TRUE)[1:12])
ps.phylum.B <- prune_taxa(phylum.10.B, ps.phyla.perc)


ps.phylum.10@otu_table
tableSyncom.B = cbind(ps.phylum.B@sam_data, ps.phylum.B@otu_table)
tableSyncom.B2 = cbind(rownames(tableSyncom.B), tableSyncom.B)
colnames(tableSyncom.B2) = c('ID',colnames(tableSyncom.B))
tableSyncom.B2$ID = as.numeric(tableSyncom.B2$ID)
tableSyncom.B3 = tableSyncom.B2[order(tableSyncom.B2$ID),]

sam = cbind(rownames(chemistry[c(33:64),-(1:5)]), chemistry[c(33:64),-(1:5)])
colnames(sam) =c('ID', colnames(chemistry[,-(1:6)]))
tableSyncom.B3$ID = as.factor(tableSyncom.B3$ID)

Bac_che <- data.frame(cbind(tableSyncom.B3,sam[,-c(1)]))

Corr.B <-cor(Bac_che[,-c(1:7)], method = 'spearman')
p.mat.B <-cor_pmat(Bac_che[,-c(1:7)], method = 'spearman')


Corr.B2 = Corr.B[1:12,-c(1:12)]
p.mat.B2 = p.mat.B[1:12,-c(1:12)]

order = colnames(Corr.B2)
ggcorrplot(Corr.B2,
           p.mat = p.mat.B2, insig = "blank",
           method = "circle",
           hc.order = TRUE, #type = "upper",
           outline.color = "white",
           ggtheme = ggplot2::theme_bw,
           colors = c("#6D9EC1", "white", "#E46726"))

ggcorrplot(Corr.F2,
           p.mat = p.mat.B2, insig = "blank",
           method = "circle",
           hc.order = TRUE, #type = "upper",
           outline.color = "white",
           ggtheme = ggplot2::theme_bw,
           colors = c("#6D9EC1", "white", "#E46726"))
