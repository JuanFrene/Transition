library("FactoMineR")
library("factoextra")
library(ade4)
library(MASS)
library(ellipse)
library(ggplot2)
library(agricolae)
library(vegan)
library(ggthemes)
library(RVAideMemoire)
library(reshape2)
library(emmeans)

setwd("G:/My Drive/labs/NMSU/Transtion/Bacteria/")

Funtional_data <- read.table("Functional table.txt", header = TRUE)
Metadata <- read.table("Metadata.txt", header = TRUE)
F_data = cbind(Metadata, t(Funtional_data))
F_data2 = F_data[,c(2,3,7:11,15:20,23:25,33,37:39,48:53)]
colnames(F_data)

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499") 

head(Funtional_data)

melt.function <- melt(F_data[,-(6)])
head(melt.function)
melt.function.mean <- melt.function%>%group_by(Irrigation, Treatment,variable)%>%
  summarise_all(mean)
melt.function.mean$ID <- paste(melt.function.mean$Irrigation, melt.function.mean$Treatment, sep="_")

melt.function.mean$ID= factor(melt.function.mean$ID, levels=c("Dryland_NCC","Dryland_GBL_NR","Dryland_GBL_R","Dryland_G","Irrigated_NCC","Irrigated_GBL_R","Irrigated_GBL_NR","Irrigated_G"))
ggplot(melt.function.mean, aes(x = ID, y = value, fill = variable)) + 
  theme_few() +
  geom_bar(stat = "identity") + 
  #scale_fill_manual(values = safe_colorblind_palette) + 
  labs(fill = "Phylum")+
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5))


#Do a model per ion
Func_em <- NULL
Func_Dryland <- NULL
#

for(specie in melted_fungal_traits[melted_fungal_traits$Irrigation!='Dryland',]$variable  %>% unique){
  melted_sub <- melted_fungal_traits[melted_fungal_traits$Irrigatio!='Dryland',] %>% subset(variable == specie) %>% droplevels
  m1 <- lm(data = melted_sub,
           formula = value ~ Treatment)
  m1_Ans <- emmeans(m1, pairwise~ Treatment, ref =1, adjust = "none")
  m1_em <- m1_Ans$emmeans %>% as.data.frame
  m1_p <- m1_Ans$contrasts %>% as.data.frame
  m1_em$Linea <- specie
  m1_p$Linea <- specie
  Func_em <- rbind(Func_em,m1_em)
  Func_Dryland <- rbind(Func_Dryland,m1_p)
}

####Adjust the p-value
Func_Dryland$Significance <- "NoSignificant"
pval_thres <- 0.05
Func_Dryland$Significance[which(Func_Dryland$p.value < pval_thres)] <- "q < 0.05"
Func_Dryland$Significance <-Func_Dryland$Significance %>% factor

Func_Dryland_selected = Func_Dryland[Func_Dryland$Significance == 'q < 0.05',]

#Func_water -> nitrogen_fixation 12 Dryland - Irrigated  0.1825000000 0.083903730 22  2.1751119 0.04064868                             nitrogen_fixation

#Irrigation
#       contrast estimate         SE df   t.ratio    p.value                              Linea Significance
#11   GBL_NR - G  -0.1350 0.04826144 12 -2.797264 0.01612586                 methanol_oxidation     q < 0.05
#12    GBL_R - G  -0.1250 0.04826144 12 -2.590059 0.02365532                 methanol_oxidation     q < 0.05
#17   GBL_NR - G  -0.1350 0.04826144 12 -2.797264 0.01612586                      methylotrophy     q < 0.05
#18    GBL_R - G  -0.1250 0.04826144 12 -2.590059 0.02365532                      methylotrophy     q < 0.05
#80  NCC - GBL_R  -0.0600 0.02571883 12 -2.332921 0.03786917                       cellulolysis     q < 0.05
#91 NCC - GBL_NR   0.0275 0.01198958 12  2.293659 0.04066154 dark_oxidation_of_sulfur_compounds     q < 0.05

#Dryland
#          contrast estimate         SE df   t.ratio     p.value                                  Linea Significance
#31    NCC - GBL_NR  -0.2325 0.09862598 12 -2.357391 0.036222749                nitrate_denitrification     q < 0.05
#37    NCC - GBL_NR  -0.2325 0.09862598 12 -2.357391 0.036222749                nitrite_denitrification     q < 0.05
#43    NCC - GBL_NR  -0.2325 0.09862598 12 -2.357391 0.036222749          nitrous_oxide_denitrification     q < 0.05
#49    NCC - GBL_NR  -0.2325 0.09862598 12 -2.357391 0.036222749                        denitrification     q < 0.05
#73    NCC - GBL_NR  -0.2325 0.09862598 12 -2.357391 0.036222749                    nitrite_respiration     q < 0.05
#82  GBL_NR - GBL_R   0.0625 0.02033572 12  3.073409 0.009655439                           cellulolysis     q < 0.05
#130 GBL_NR - GBL_R  -0.0975 0.04124369 12 -2.363998 0.035790121                              human_gut     q < 0.05
#132      GBL_R - G   0.0925 0.04124369 12  2.242768 0.044574335                              human_gut     q < 0.05
#142 GBL_NR - GBL_R  -0.0975 0.04124369 12 -2.363998 0.035790121                             mammal_gut     q < 0.05
#144      GBL_R - G   0.0925 0.04124369 12  2.242768 0.044574335                             mammal_gut     q < 0.05
#222      GBL_R - G   0.0725 0.03144108 12  2.305900 0.039770435        nonphotosynthetic_cyanobacteria     q < 0.05
#229   NCC - GBL_NR  -0.2475 0.07908988 12 -3.129351 0.008702338 anoxygenic_photoautotrophy_S_oxidizing     q < 0.05
#232 GBL_NR - GBL_R   0.1750 0.07908988 12  2.212673 0.047053999 anoxygenic_photoautotrophy_S_oxidizing     q < 0.05
#233     GBL_NR - G   0.1800 0.07908988 12  2.275892 0.041988797 anoxygenic_photoautotrophy_S_oxidizing     q < 0.05
#235   NCC - GBL_NR  -0.2475 0.07908988 12 -3.129351 0.008702338             anoxygenic_photoautotrophy     q < 0.05
#238 GBL_NR - GBL_R   0.1750 0.07908988 12  2.212673 0.047053999             anoxygenic_photoautotrophy     q < 0.05
#239     GBL_NR - G   0.1800 0.07908988 12  2.275892 0.041988797             anoxygenic_photoautotrophy     q < 0.05
#247   NCC - GBL_NR  -0.2525 0.07608247 12 -3.318767 0.006124145                        photoautotrophy     q < 0.05
#250 GBL_NR - GBL_R   0.1700 0.07608247 12  2.234418 0.045249599                        photoautotrophy     q < 0.05
#251     GBL_NR - G   0.1700 0.07608247 12  2.234418 0.045249599                        photoautotrophy     q < 0.05
#278    NCC - GBL_R  -8.1300 3.51509913 12 -2.312879 0.039270798            anaerobic_chemoheterotrophy     q < 0.05
#280 GBL_NR - GBL_R  -8.8775 3.51509913 12 -2.525533 0.026637099            anaerobic_chemoheterotrophy     q < 0.05

F_data$Treatment = factor(F_data$Treatment, c('NCC','GBL_NR','GBL_R','G'))
ggplot(F_data, aes(Irrigation, methanotrophy,  fill = Treatment)) + 
  geom_boxplot() + 
  theme_few()+
  scale_fill_manual(values = safe_colorblind_palette)


melted_F_data <- F_data2 %>% melt
melted_F_data$Treatment = factor(melted_F_data$Treatment, c('NCC', 'GBL_NR', 'GBL_R','G'))
melted_F_data = cbind(melted_F_data,log2(melted_F_data[4]+1))

colnames(melted_F_data) = c('Irrigation','Treatment','variable','value','log_value')
melted_F_data$Significance <- "No Significant"

ggplot(data = melted_F_data, aes(Treatment,variable)) +
  geom_raster(aes(fill = log_value))+
  theme_few() +
  facet_grid(.~Irrigation)+
  #geom_text(aes(label = melterd2), color = "black", size = 1.5) +
  geom_tile(aes(color = Significance),fill = '#00000000', linewidth = 0.2,width = 1,height = 1) + #
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(0,4),name = "Corr") + #,na.value = "#D9D9D9"
  scale_color_manual(values = c("black",'grey'),na.value =  "transparent",name = "Significative") +
  labs(y='Microbial paramenters',x='Chemical parameters')+
  theme(axis.text.x = element_text(angle = -45, hjust=-0.05), #element_blank()
        axis.title = element_blank(),legend.position = 'top')


