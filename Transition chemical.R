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

setwd("G:/My Drive/labs/NMSU/Transtion/")

chemistry <- read.table("Chemical properties.txt", header = TRUE)

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499") 

head(chemistry)
aovH2O <- aov(H2O~Treatment*Irrigation*Year, data=chemistry)
summary(aovH2O)

aovH2O_I <- aov(H2O~Irrigation, data=chemistry[chemistry$Treatment!='G',])
summary(aovH2O_I)
#NS

aovH2O_Y <- aov(H2O~Year, data=chemistry)
summary(aovH2O_Y)
# Df  Sum Sq Mean Sq F value Pr(>F)    
#Year         1 0.06238 0.06238   130.8 <2e-16 ***
#Residuals   62 0.02958 0.00048   

aovH2O_Ir <- aov(H2O~Treatment*Year, data=chemistry[chemistry$Irrigation != 'Dryland',])
summary(aovH2O_Ir)
#NS

aovH2O_Dr <- aov(H2O~Treatment*Year, data=chemistry[chemistry$Irrigation == 'Dryland',])
summary(aovH2O_Dr)
#Df  Sum Sq Mean Sq F value   Pr(>F)    
#Treatment       3 0.00652 0.00217   7.084  0.00143 ** 
#Year            1 0.03976 0.03976 129.650 3.68e-11 ***
#Treatment:Year  3 0.00460 0.00153   5.005  0.00777 ** 
#Residuals      24 0.00736 0.00031

LSDoav <- HSD.test(aovH2O_Dr, c("Treatment":"Year"), alpha = 0.05, group=TRUE) 
LSDoav

chemistry$Year = factor(chemistry$Year)
ggplot(data=chemistry, aes(x=Year, y=H2O, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  facet_grid(.~Irrigation,space = "fixed",scales = "free_y") +
  ggtitle("H20")+ ylab("nmol / g soil")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Year")+
  theme_bw()

####SOC
head(chemistry)
aovSOC <- aov(SOC~Treatment*Irrigation*Year, data=chemistry)
summary(aovSOC)

aovSOC_I <- aov(SOC~Irrigation, data=chemistry[chemistry$Treatment!='G',])
summary(aovSOC_I)
# Df Sum Sq Mean Sq F value Pr(>F)  
#Irrigation   1   3.62   3.616   5.919 0.0179 *
#Residuals   62  37.87   0.611    

LSDoav <- HSD.test(aovSOC_I, c("Irrigation"), alpha = 0.05, group=TRUE) 
LSDoav
#SOC groups
#Dryland   7.958531      a
#Irrigated 7.483156      b

aovSOC_Y <- aov(SOC~Year, data=chemistry)
summary(aovSOC_Y)
# NS

aovSOC_Ir <- aov(SOC~Treatment*Year, data=chemistry[chemistry$Irrigation != 'Dryland',])
summary(aovSOC_Ir)
# Df Sum Sq Mean Sq F value  Pr(>F)    
#Treatment       3  6.956  2.3185  11.757 6.2e-05 ***
#Year            1  2.576  2.5759  13.062 0.00139 ** 
#Treatment:Year  3  1.285  0.4284   2.172 0.11754    
#Residuals      24  4.733  0.1972

LSDoav <- HSD.test(aovSOC_Ir, c("Treatment","Year"), alpha = 0.05, group=TRUE) 
LSDoav
#SOC groups
#GO:2023     8.26750      a
#GBL:2023  8.02250     ab
#GBL-Rem:2023 7.98750     ab
#GBL:2024  7.83900     ab 2.3
#GO:2024     7.23050    abc
#GBL-Rem:2024 7.09475     bc 12$
#Fallow:2023 6.79000      c
#Fallow:2024 6.63350      c 2.5

aovSOC_Dr <- aov(SOC~Treatment*Year, data=chemistry[chemistry$Irrigation == 'Dryland',])
summary(aovSOC_Dr)
#Df Sum Sq Mean Sq F value   Pr(>F)    
#Treatment       3 11.391   3.797  12.408 4.24e-05 ***
#Year            1  0.004   0.004   0.013   0.9111    
#Treatment:Year  3  3.582   1.194   3.902   0.0211 *  
#Residuals      24  7.344   0.306   

LSDoav <- HSD.test(aovSOC_Dr, c("Treatment",'Year'), alpha = 0.05, group=TRUE) 
LSDoav

#SOC groups
#GBL-Rem:2023 8.73000      a
#GBL:2024  8.65500      a
#PG:2024     8.49500      a
#GBL:2023  8.45750      a
#GBL-Rem:2024 7.83825     ab 11%
#PG:2023     7.53250     ab
#Fallow:2023 7.07000      b
#Fallow:2024 6.89000      b 2.3

ggplot(data=chemistry, aes(x=Year, y=SOC, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  facet_grid(.~Irrigation,space = "fixed",scales = "free_y") +
  ggtitle("SOC")+ ylab("mg / kg soil")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Year")+
  theme_bw()

####MAOC
head(chemistry)
aovMAOC <- aov(MAOC~Treatment*Irrigation*Year, data=chemistry)
summary(aovMAOC)

aovMAOC_I <- aov(MAOC~Irrigation, data=chemistry[chemistry$Treatment!='G',])
summary(aovMAOC_I)
#NS
aovMAOC_Y <- aov(MAOC~Year, data=chemistry)
summary(aovMAOC_Y)
#Df Sum Sq Mean Sq F value   Pr(>F)    
#Year         1  6.197   6.197   23.04 1.04e-05 ***
#Residuals   62 16.676   0.269

#MAOC Irrigation
aovMAOC_Ir <- aov(MAOC~Treatment*Year, data=chemistry[chemistry$Irrigation != 'Dryland',])
summary(aovMAOC_Ir)
#  Df Sum Sq Mean Sq F value   Pr(>F)    
#Treatment       3  5.175   1.725   9.165 0.000319 ***
#Year            1  4.332   4.332  23.018 6.95e-05 ***
#Treatment:Year  3  0.460   0.153   0.814 0.498702    
#Residuals      24  4.517   0.188  

LSDoav <- HSD.test(aovMAOC_Ir, c("Treatment","Year"), alpha = 0.05, group=TRUE) 
LSDoav
#MAOC groups
#GBL-R:2023  7.693907      a
#GO:2023     7.434131     ab
#GBL-NR:2023 7.280493     ab
#GBL-R:2024  7.075867    abc
#GO:2024     6.649747     bc
#Fallow:2023 6.515293     bc
#GBL-NR:2024 6.184972      c
#Fallow:2024 6.069645      c

#MAOC Dryland
aovMAOC_Dr <- aov(MAOC~Treatment*Year, data=chemistry[chemistry$Irrigation == 'Dryland',])
summary(aovMAOC_Dr)
#Df Sum Sq Mean Sq F value   Pr(>F)    
#Treatment       3  2.8573  0.9524  10.345 0.000147 ***
#Year            1 2.0706  2.0706  22.490 7.99e-05 ***
#  Treatment:Year  3 0.8599  0.2866   3.113 0.05080   
# Residuals      24 2.2096  0.0921*  
  

LSDoav <- HSD.test(aovMAOC_Dr, c("Treatment",'Year'), alpha = 0.05, group=TRUE) 
LSDoav

#MAOC groups
#GBL-NR:2023 8.73000      a
#GBL-R:2024  8.65500      a
#PG:2024     8.49500      a
#GBL-R:2023  8.45750      a
#GBL-NR:2024 7.83825     ab
#PG:2023     7.53250     ab
#Fallow:2023 7.07000      b
#Fallow:2024 6.89000      b

ggplot(data=chemistry, aes(x=Year, y=MAOC, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  facet_grid(.~Irrigation,space = "fixed",scales = "free_y") +
  ggtitle("MAOC")+ ylab("g / kg soil")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Year")+
  theme_bw()

####POC
head(chemistry)
aovPOC <- aov(POC~Treatment*Irrigation*Year, data=chemistry)
summary(aovPOC)

aovPOC_I <- aov(POC~Irrigation, data=chemistry[chemistry$Treatment!='G',])
summary(aovPOC_I)
#Df Sum Sq Mean Sq F value Pr(>F)  
#Irrigation   1  1.131  1.1315   4.789 0.0324 *
#Residuals   62 14.648  0.2363 

LSDoav <- HSD.test(aovPOC_I, c("Irrigation"), alpha = 0.05, group=TRUE) 
LSDoav
#POC groups
#Dryland   0.5691033      a
#Irrigated 0.3031784      b

aovPOC_Y <- aov(POC~Year, data=chemistry)
summary(aovPOC_Y)
#Df Sum Sq Mean Sq F value   Pr(>F)    
#Year         1   6.762   6.762    46.5 4.42e-09 ***
#Residuals   62  9.017   0.145

aovPOC_Ir <- aov(POC~Treatment*Year, data=chemistry[chemistry$Irrigation != 'Dryland',])
summary(aovPOC_Ir)
#  Df Sum Sq Mean Sq F value   Pr(>F)    
#Treatment       3  0.4687  0.1562   2.623   0.0738 .  
#Year            1 1.7334  1.7334  29.097 1.54e-05 ***
#Treatment:Year  3 0.4504  0.1501   2.520   0.0820 .  
#Residuals      24 1.4298  0.0596  


aovPOC_Dr <- aov(POC~Treatment*Year, data=chemistry[chemistry$Irrigation == 'Dryland',])
summary(aovPOC_Dr)
#Df Sum Sq Mean Sq F value   Pr(>F)    
#Treatment       3 0.897   0.299   2.150    0.120    
#Year            1  5.574   5.574  40.101 1.51e-06 ***
#Treatment:Year  3  0.759   0.253   1.819    0.171    
#Residuals      24  3.336   0.139  

LSDoav <- HSD.test(aovPOC_Dr, c("Treatment",'Year'), alpha = 0.05, group=TRUE) 
LSDoav

ggplot(data=chemistry, aes(x=Year, y=POC, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  facet_grid(.~Irrigation,space = "fixed",scales = "free_y") +
  ggtitle("POC")+ ylab("g / kg soil")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Year")+
  theme_bw()

####MBC
head(chemistry)
aovMBC <- aov(MBC~Treatment*Irrigation*Year, data=chemistry)
summary(aovMBC)

aovMBC_I <- aov(MBC~Irrigation, data=chemistry[chemistry$Treatment!='G',])
summary(aovMBC_I)
#NS
aovMBC_Y <- aov(MBC~Year, data=chemistry)
summary(aovMBC_Y)
#NS

#MBC Irrigation
aovMBC_Ir <- aov(MBC~Treatment*Year, data=chemistry[chemistry$Irrigation != 'Dryland',])
summary(aovMBC_Ir)
#  Df Sum Sq Mean Sq F value   Pr(>F)    
#Treatment       26437    8812   7.863 0.000797 ***
#Year            1     34      34   0.030 0.863881    
#Treatment:Year  3   9143    3048   2.719 0.066896 .  
#Residuals      24  26899    1121   

LSDoav <- HSD.test(aovMBC_Ir, c("Treatment","Year"), alpha = 0.05, group=TRUE) 
LSDoav
# MBC groups
#GBL-R:2024  284.0319      a
#GBL-NR:2023 272.4437      a
#GBL-R:2023  266.9401      a
#GO:2024     254.1345     ab
#GO:2023     240.9110     ab
#GBL-NR:2024 212.0784     ab
#Fallow:2024 206.3828     ab
#Fallow:2023 184.5369      b

aovMBC_Dr <- aov(MBC~Treatment*Year, data=chemistry[chemistry$Irrigation == 'Dryland',])
summary(aovMBC_Dr)
#Df Sum Sq Mean Sq F value   Pr(>F)    
#Treatment       3 30335   10112   6.777 0.00181 **
#Year            1   6528    6528   4.375 0.04723 * 
#Treatment:Year  3  20447    6816   4.568 0.01144 * 
#Residuals      24  35811    1492      

LSDoav <- HSD.test(aovMBC_Dr, c("Treatment",'Year'), alpha = 0.05, group=TRUE) 
LSDoav

#MBC groups
#PG:2024     317.5320      a
#GBL-R:2024  250.7599     ab
#GBL-R:2023  240.9251     ab
#GBL-NR:2024 232.6071     ab
#GBL-NR:2023 229.3709     ab
#PG:2023     202.6688      b
#Fallow:2023 185.4665      b
#Fallow:2024 171.7949      b

ggplot(data=chemistry, aes(x=Year, y=MBC, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  facet_grid(.~Irrigation,space = "fixed",scales = "free_y") +
  ggtitle("MBC")+ ylab("g / kg soil")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Year")+
  theme_bw()

####PMC
head(chemistry)
aovPMC <- aov(PMC~Treatment*Irrigation*Year, data=chemistry)
summary(aovPMC)

aovPMC_I <- aov(PMC~Irrigation, data=chemistry[chemistry$Treatment!='G',])
summary(aovPMC_I)
#NS
aovPMC_Y <- aov(PMC~Year, data=chemistry)
summary(aovPMC_Y)
#NS

aovPMC_Ir <- aov(PMC~Treatment*Year, data=chemistry[chemistry$Irrigation != 'Dryland',])
summary(aovPMC_Ir)
#  Df Sum Sq Mean Sq F value   Pr(>F)    
#Treatment       3  749.1  249.68   3.138  0.044 *
#Year            1    3.0    2.98   0.037  0.848  
#Treatment:Year  3  172.7   57.58   0.724  0.548  
#Residuals      24 1909.9   79.58                 

LSDoav <- HSD.test(aovPMC_Ir, c("Treatment","Year"), alpha = 0.05, group=TRUE) 
LSDoav
# PMC groups
#GBL-R:2024  36.27369      a
#GO:2024     34.40181      a
#GBL-NR:2023 32.14879      a
#GBL-R:2023  30.66601      a
#GO:2023     29.50966      a
#GBL-NR:2024 28.09498      a
#Fallow:2023 23.01961      a
#Fallow:2024 19.01333      a

#PMC Dryland
aovPMC_Dr <- aov(PMC~Treatment*Year, data=chemistry[chemistry$Irrigation == 'Dryland',])
summary(aovPMC_Dr)
#Df Sum Sq Mean Sq F value   Pr(>F)    
#Treatment       3 1881.0   627.0   8.061 0.00069 ***
# Year            1  371.2   371.2   4.773 0.03892 *  
#Treatment:Year  3 1742.6   580.9   7.468 0.00107 ** 
#Residuals      24 1866.8    77.8     

LSDoav <- HSD.test(aovPMC_Dr, c("Treatment",'Year'), alpha = 0.05, group=TRUE) 
LSDoav

#PMC groups
#PG:2024     50.27511      a
#GBL-R:2024  48.27461      a
#GBL-NR:2023 34.11306     ab
#GBL-R:2023  32.73804     ab
#GBL-NR:2024 32.24611     ab
#Fallow:2023 26.46152      b
#PG:2023     24.51435      b
#Fallow:2024 14.27879      b

ggplot(data=chemistry, aes(x=Year, y=PMC, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  facet_grid(.~Irrigation,space = "fixed",scales = "free_y") +
  ggtitle("PMC")+ ylab("mg / kg soil")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Year")+
  theme_bw()


####TN
head(chemistry)
aovTN <- aov(TN~Treatment*Irrigation*Year, data=chemistry)
summary(aovTN)

aovTN_I <- aov(TN~Irrigation, data=chemistry[chemistry$Treatment!='G',])
summary(aovTN_I)
#NS
aovTN_Y <- aov(TN~Year, data=chemistry)
summary(aovTN_Y)
#NS

aovTN_Ir <- aov(TN~Treatment*Year, data=chemistry[chemistry$Irrigation != 'Dryland',])
summary(aovTN_Ir)
#   Df Sum Sq Mean Sq F value Pr(>F)  
#Treatment       3 0.0139 0.00465   0.109 0.9537  
#Year            1 0.0043 0.00425   0.100 0.7545  
#Treatment:Year  3 0.3914 0.13048   3.073 0.0469 *
#Residuals      24 1.0192 0.04247                  


aovTN_Dr <- aov(TN~Treatment*Year, data=chemistry[chemistry$Irrigation == 'Dryland',])
summary(aovTN_Dr)
#Df Sum Sq Mean Sq F value   Pr(>F)    
#Treatment       3 0.3087 0.10289   3.497  0.031 *
#Year            1 0.0014 0.00139   0.047  0.830  
#Treatment:Year  3 0.0745 0.02484   0.844  0.483  
#Residuals      24 0.7062 0.02942     

LSDoav <- HSD.test(aovTN_Dr, c("Treatment",'Year'), alpha = 0.05, group=TRUE) 
LSDoav

#TN groups
#PG:2024     50.27511      a
#GBL-R:2024  48.27461      a
#GBL-NR:2023 34.11306     ab
#GBL-R:2023  32.73804     ab
#GBL-NR:2024 32.24611     ab
#Fallow:2023 26.46152      b
#PG:2023     24.51435      b
#Fallow:2024 14.27879      b

ggplot(data=chemistry, aes(x=Year, y=TN, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  facet_grid(.~Irrigation,space = "fixed",scales = "free_y") +
  ggtitle("TN")+ ylab("g / kg soil")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Year")+
  theme_bw()

####PMN
head(chemistry)
aovPMN <- aov(PMN~Treatment*Irrigation*Year, data=chemistry)
summary(aovPMN)

aovPMN_I <- aov(PMN~Irrigation, data=chemistry[chemistry$Treatment!='G',])
summary(aovPMN_I)
# Df Sum Sq Mean Sq F value   Pr(>F)    
#Irrigation   1  62.62   62.62    16.3 0.000151 ***
#Residuals   62 238.26    3.84


aovPMN_Y <- aov(PMN~Year, data=chemistry)
summary(aovPMN_Y)
#NS

aovPMN_Ir <- aov(PMN~Treatment*Year, data=chemistry[chemistry$Irrigation != 'Dryland',])
summary(aovPMN_Ir)
#   Df Sum Sq Mean Sq F value Pr(>F)  
#Treatment       3 12.50    4.17   3.116   0.0449 *  
#Year            1  37.84   37.84  28.300 1.85e-05 ***
#Treatment:Year  3   4.22    1.41   1.053   0.3874    
#Residuals      24  32.09    1.34                   

LSDoav <- HSD.test(aovPMN_Ir, c("Treatment",'Year'), alpha = 0.05, group=TRUE) 
LSDoav

#PMN groups
#GBL-R:2023  7.628055      a
#GBL-NR:2023 6.626617     ab
#GO:2023     6.273610     ab
#Fallow:2023 5.396698    abc
#GBL-R:2024  5.019237    abc
#GBL-NR:2024 4.581805     bc
#Fallow:2024 4.338365     bc
#GO:2024     3.286535      c

aovPMN_Dr <- aov(PMN~Treatment*Year, data=chemistry[chemistry$Irrigation == 'Dryland',])
summary(aovPMN_Dr)
#Df Sum Sq Mean Sq F value   Pr(>F)    
#Treatment       3 34.17   11.39   5.075 0.007309 ** 
#Year            1  47.77   47.77  21.286 0.000111 ***
#Treatment:Year  3  15.82    5.27   2.349 0.097749 .  
#Residuals      24  53.86    2.24    

LSDoav <- HSD.test(aovPMN_Dr, c("Treatment",'Year'), alpha = 0.05, group=TRUE) 
LSDoav

#PMN groups
#GBL-R:2024  10.139408      a
#GBL-NR:2024  9.706489     ab
#Fallow:2024  8.548912    abc
#GBL-NR:2023  6.959133    abc
#GBL-R:2023   6.626461     bc
#PG:2024      5.981004      c
#PG:2023      5.916593      c
#Fallow:2023  5.099350      c

ggplot(data=chemistry, aes(x=Year, y=PMN, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  facet_grid(.~Irrigation,space = "fixed",scales = "free_y") +
  ggtitle("PMN")+ ylab("mg / kg soil")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Year")+
  theme_bw()

####NT0
head(chemistry)
aovNT0 <- aov(NT0~Treatment*Irrigation, data=chemistry)
summary(aovNT0)

aovNT0_I <- aov(NT0~Irrigation, data=chemistry[chemistry$Treatment!='G',])
summary(aovNT0_I)
# Df Sum Sq Mean Sq F value   Pr(>F)    
#Irrigation    1  818.8   818.8   19.09 0.000137 ***
#Residuals   30 1287.0    42.9 

aovNT0_Ir <- aov(NT0~Treatment, data=chemistry[chemistry$Irrigation != 'Dryland',])
summary(aovNT0_Ir)
#   Df Sum Sq Mean Sq F value Pr(>F)  
#Treatment       3 84.20  28.067    3.63 0.0451 *
#Residuals   12  92.77   7.731                  

LSDoav <- HSD.test(aovNT0_Ir, c("Treatment"), alpha = 0.05, group=TRUE) 
LSDoav

#NT0 groups
#Fallow 10.590594      a
#GBL-R   8.336341     ab
#GBL-NR  6.286856     ab
#GO      4.441102      b

aovNT0_Dr <- aov(NT0~Treatment, data=chemistry[chemistry$Irrigation == 'Dryland',])
summary(aovNT0_Dr)
#Df Sum Sq Mean Sq F value   Pr(>F)    
#Treatment        3  799.7  266.57   10.31 0.00122 **
#Residuals   12  310.3   25.86   

LSDoav <- HSD.test(aovNT0_Dr, c("Treatment"), alpha = 0.05, group=TRUE) 
LSDoav

#NT0 groups
#Fallow 25.440235      a
#GBL-NR 20.417338      a
#GBL-R  18.079293      a
#PG      6.184726      b

ggplot(data=chemistry, aes(x=Irrigation, y=NT0, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  #facet_grid(.~Irrigation,space = "fixed",scales = "free_y") +
  ggtitle("NT0")+ ylab("mg / kg soil")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Year")+
  theme_bw()

####MAON
head(chemistry)
aovMAON <- aov(MAON~Treatment*Irrigation*Year, data=chemistry)
summary(aovMAON)

aovMAON_I <- aov(MAON~Irrigation, data=chemistry[chemistry$Treatment!='G',])
summary(aovMAON_I)
# Df Sum Sq Mean Sq F value   Pr(>F)    
#Irrigation     1 0.1179 0.11794   5.325 0.0244 *
#Residuals   62 1.3732 0.02215   

aovMAON_I <- aov(MAON~Year, data=chemistry)
summary(aovMAON_I)
#NS

aovMAON_Ir <- aov(MAON~Treatment*Year, data=chemistry[chemistry$Irrigation != 'Dryland',])
summary(aovMAON_Ir)
#NS

aovMAON_Dr <- aov(MAON~Treatment, data=chemistry[chemistry$Irrigation == 'Dryland',])
summary(aovMAON_Dr)
#NS

ggplot(data=chemistry, aes(x=Year, y=MAON, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  facet_grid(.~Irrigation,space = "fixed",scales = "free_y") +
  ggtitle("MAON")+ ylab("mg / kg soil")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Year")+
  theme_bw()

####PON
head(chemistry)
aovPON <- aov(PON~Treatment*Irrigation*Year, data=chemistry)
summary(aovPON)

aovPON_I <- aov(PON~Irrigation, data=chemistry[chemistry$Treatment!='G',])
summary(aovPON_I)
#NS

aovPON_Y <- aov(PON~Year, data=chemistry)
summary(aovPON_Y)
#NS

aovPON_Ir <- aov(PON~Treatment*Year, data=chemistry[chemistry$Irrigation != 'Dryland',])
summary(aovPON_Ir)
#NS

aovPON_Dr <- aov(PON~Treatment, data=chemistry[chemistry$Irrigation == 'Dryland',])
summary(aovPON_Dr)
#NS

ggplot(data=chemistry, aes(x=Year, y=PON, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  facet_grid(.~Irrigation,space = "fixed",scales = "free_y") +
  ggtitle("PON")+ ylab("mg / kg soil")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Year")+
  theme_bw()

####BG
head(chemistry)
aovBG <- aov(BG~Treatment*Irrigation, data=chemistry)
summary(aovBG)

aovBG_I <- aov(BG~Irrigation, data=chemistry[chemistry$Treatment!='G',])
summary(aovBG_I)
# Df Sum Sq Mean Sq F value   Pr(>F)    
#Irrigation    1  1   4114    4114   5.872 0.0216 *
#Residuals   30  21017     701  

LSDoav <- HSD.test(aovBG_I, c("Irrigation"), alpha = 0.05, group=TRUE) 
LSDoav

aovBG_Ir <- aov(BG~Treatment, data=chemistry[chemistry$Irrigation != 'Dryland',])
summary(aovBG_Ir)
#   Df Sum Sq Mean Sq F value Pr(>F)  
#Treatment       3 3   8742  2914.0   5.118 0.0165 *
#Residuals   12   6832   569.4                   

LSDoav <- HSD.test(aovBG_Ir, c("Treatment"), alpha = 0.05, group=TRUE) 
LSDoav

#BG groups
#BL-R  108.50207      a
#GBL-NR  84.13770     ab
#GO      65.93287     ab
#Fallow  44.99046      b

aovBG_Dr <- aov(BG~Treatment, data=chemistry[chemistry$Irrigation == 'Dryland',])
summary(aovBG_Dr)
#NS


ggplot(data=chemistry, aes(x=Irrigation, y=BG, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  #facet_grid(.~Irrigation,space = "fixed",scales = "free_y") +
  ggtitle("BG")+ ylab("mg / kg soil")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Year")+
  theme_bw()

####NAG
head(chemistry)
aovNAG <- aov(NAG~Treatment*Irrigation, data=chemistry)
summary(aovNAG)

aovNAG_I <- aov(NAG~Irrigation, data=chemistry[chemistry$Treatment!='G',])
summary(aovNAG_I)
# NS

aovNAG_Ir <- aov(NAG~Treatment, data=chemistry[chemistry$Irrigation != 'Dryland',])
summary(aovNAG_Ir)
#NS


aovNAG_Dr <- aov(NAG~Treatment, data=chemistry[chemistry$Irrigation == 'Dryland',])
summary(aovNAG_Dr)
#Df Sum Sq Mean Sq F value   Pr(>F)    
#Treatment    3 2244.0   748.0   14.22 0.000296 ***
#Residuals   12  631.2    52.6 

LSDoav <- HSD.test(aovNAG_Dr, c("Treatment"), alpha = 0.05, group=TRUE) 
LSDoav

#NAG groups
#GBL   58.42219      a
#PG     48.15289     ab
#GBL-Rem 39.95768     bc
#Fallow 26.04607      c



ggplot(data=chemistry, aes(x=Irrigation, y=NAG, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  #facet_grid(.~Irrigation,space = "fixed",scales = "free_y") +
  ggtitle("NAG")+ ylab("mg / kg soil")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Year")+
  theme_bw()


C_N = chemistry$SOC / chemistry$TN

chemistry2 = cbind(chemistry, C_N)

###C_N
head(chemistry)
aovMAON <- aov(C_N~Treatment*Irrigation*Year, data=chemistry2)
summary(aovMAON)

aovMAON_I <- aov(C_N~Irrigation, data=chemistry2[chemistry2$Treatment!='G',])
summary(aovMAON_I)
# Df Sum Sq Mean Sq F value   Pr(>F)    
#Irrigation     1 0.1179 0.11794   5.325 0.0244 *
#Residuals   62 1.3732 0.02215   

aovMAON_I <- aov(C_N~Year, data=chemistry2)
summary(aovMAON_I)
#NS

aovMAON_Ir <- aov(C_N~Treatment*Year, data=chemistry2[chemistry2$Irrigation != 'Dryland',])
summary(aovMAON_Ir)
#NS

aovMAON_Dr <- aov(C_N~Treatment*Year, data=chemistry2[chemistry2$Irrigation == 'Dryland',])
summary(aovMAON_Dr)
#NS

ggplot(data=chemistry2, aes(x=Year, y=C_N, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  facet_grid(.~Irrigation,space = "fixed",scales = "free_y") +
  ggtitle("MAON")+ ylab("mg / kg soil")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Year")+
  theme_bw()

library(ggcorrplot)
##Corr y p-value
Corr <-cor(chemistry[,-c(1:7)], method = 'pearson', use='pairwise.complete.obs')
p.mat <-cor_pmat(chemistry[,-c(1:7)], method = 'pearson')

ggcorrplot(Corr,
           p.mat = p.mat, insig = "blank",
           method = "circle",
           hc.order = TRUE, type = "upper",
           #outline.color = "white",
           ggtheme = ggplot2::theme_bw,
           colors = c("#6D9EC1", "white", "#E46726")
)


##############################################

res.pca <- prcomp(chemistry[33:64,-c(1:5)], scale = TRUE, center = TRUE)
fviz_eig(res.pca)
summary(res.pca)


load = data.frame(res.pca$rotation*2)

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
  annotate("text",x=load$PC1*1.5,y=load$PC2,label=rownames(load),size=2.5)+
  geom_vline(xintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_hline(yintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_segment(data=load,aes(x=0,y=0,xend=PC1,yend=PC2),arrow=arrow(length=unit(.3,"lines")),color="black",size=.4)+
  geom_point(size=3, aes(color = Treatment, shape=Irrigation),stroke = 1) + #size = 4, #,  shape=Year
  labs(x=xlab_text.F, y=ylab_text.F) +
  scale_fill_manual(values = safe_colorblind_palette) +
  scale_color_manual(values = safe_colorblind_palette)+
  theme_bw()+ theme(legend.position = 'right')

