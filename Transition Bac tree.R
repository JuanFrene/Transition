###################################################################################################################################################
#https://yulab-smu.top/treedata-book/
#https://epirhandbook.com/en/phylogenetic-trees-1.html
#https://yulab-smu.top/treedata-book/chapter4.html
#https://yulab-smu.top/treedata-book/chapter10.html
#https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/advanced-models-for-differential-abundance.html#deseq2-differential-abundance-testing-for-sequencing-data
#https://rstudio.com/products/rstudio/download/#download
#https://cran.r-project.org/bin/windows/Rtools/
#####################################################################################################################################################
library(seqinr)
library(adegenet)
library(ape)
library(ggtree)
library(DECIPHER)
library(viridis)
library(ggplot2)
library("BiocParallel")
library("readxl")
library(ggtree)
library(ggtreeExtra)
library(dplyr)
library(reshape2)
library(treeio)
library(readr)
library(ggstar)
library(Biostrings)
library(ggnewscale)
#options(getClass.msg=FALSE)

###################################################################################################################################################
# Set working directories- Windows-HP
setwd("G:/My Drive/labs/NMSU/Transtion/Bacteria/Tree/")#CC Network/
getwd()
register(SnowParam(6))
##################################################################################################################################################
# load the sequences from the file
# change "DNA" to "RNA" or "AA" if necessary
seqs <- readDNAStringSet("final.treatment CC D.fasta", format = "fasta")

# look at some of the sequences (optional)
seqs

# nucleotide sequences need to be in the same orientation
# if they are not, then they can be reoriented (optional)
seqs <- OrientNucleotides(seqs)

# perform the alignment
aligned <- AlignSeqs(seqs)

# view the alignment in a browser (optional)
#BrowseSeqs(aligned, highlight=2)

# write the alignment to a new FASTA file
writeXStringSet(aligned, file="Bacteria_fasta_aligned.fa")

# read in the aligned data
dna <- read.alignment("Bacteria_fasta_aligned.fa", format = "fasta")

# create a distance matrix for the alignment 
D <- dist.alignment(dna, matrix = "similarity")

#Write distance table to file and save
df <- melt(as.matrix(D), varnames = c("row", "col"))
df[df$row > df$col,]
write.csv(df, file ="distance_table.csv")
#grep 'NA' distance_table.csv | cut -d, -f1 --complement | cut -d, -f3 --complement | sed 's/","/ /g' | sed 's/"//g' | tr -s ' ' '\n' | sort | uniq -c -d > problem_fasta_sequences

#View distance matrix. Darker shades of gray mean a larger distance # you can also make cool color plots but they're much more complicated because they use the image() function
temp <- as.data.frame(as.matrix(D))
table.paint(temp, cleg=0, clabel.row=.5, clabel.col=.5) + scale_color_viridis()

#create tree object
treeD <- nj(D)
class(treeD) #all trees created using {ape} package will be of class phylo

treeD <- ladderize(treeD)

#check the formatting of the tip.labels in the tree
treeD$tip.label

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ metadata ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import metadata from excel file- Linux
metadata2 <- read.table("Metadata_tree_CC.txt", header = TRUE)
Samples_ID = row.names(metadata2)
metadata2 = cbind(Samples_ID,metadata2)
# inspect the first 6 Sample_IDs
head(metadata2$Samples_ID) 
nrow(metadata2)
#determine if all samples are present in the tree file and vice versa
metadata2$Samples_ID %in% treeD$tip.label
treeD$tip.label %in% metadata2$Samples_ID

#show any sample IDs that are not on the tree
metadata2$Samples_ID[!tree$tip.label %in% metadata2$Samples_ID]#character(0) = none


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ draw tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mypal=c('#d4c29e','#de5275','#6aceee','#8e67a0','#e9ae8c',"#55003b",'#a0ae00','#e0e004',
        '#195850','#535121','#e0aa02','#f35828',"#c1204c",'#5e3a49',"#88CCEE", "#CC6677",
        "#DDCC77", "#117733","#AA4419", "#44AA59")

#linear with sequence lengths
dat <- data.frame(id=c(1,2,3,5,14,15,17,18,19,23,24,25,26,28,29,30,50,51,52,57,58),type=c("Thermoleophilia_unclassified",'RB41','Candidatus_Yanofskybacteria_ge', 'Actinobacteriota_unclassified',
                                          'Aeromicrobium','Hamadaea','S0134_terrestrial_group_ge','Xanthomonadaceae_unclassified','Segetibacter',
                                          'Propionibacteriaceae_unclassified','Subgroup_17_ge','Rubrobacter','Ilumatobacteraceae_ge','MB-A2-108_ge','TK10_ge',
                                          'Gemmatimonadaceae_ge','MB-A2-108_ge','Vicinamibacteraceae_ge','Gaiellales_ge','TRA3-20_ge','JG30-KF-CM45_ge')) #

ggtree(treeD,  size=0.25) %<+% metadata +  #color='black', gg_mcmc_tree2 <-
  geom_tiplab(aes(label=Genus), offset=0.02, hjust = 0.5, vjust=0.5, size=2.5, align = T, linetype = "dotted", linesize = 0.05) +
  geom_treescale(fontsize=3, linesize=0.4)
    geom_hilight(data=dat, mapping=aes(node=id, fill=type),
               align="right",extend=0.3,show.legend=FALSE,alpha=0.3) 
metadata


gg_mcmc_tree <- ggtree(treeD, layout="fan", color='black', size=0.25) %<+% metadata +  #
  geom_tiplab(aes(label=Genus), offset=0.02, hjust = 0.5, vjust=0.5, size=2.5, align = T, linetype = "dotted", linesize = 0.05) +
  scale_fill_manual(values = mypal, guide=guide_legend(keywidth=0.5, keyheight=0.5, order=6))+ 
  geom_treescale(fontsize=3, linesize=0.4)
  geom_hilight(data=dat, mapping=aes(node=id, fill=type),
             align="right",extend=0.2,show.legend=TRUE,alpha=0.3) 


Species_data <- metadata2[, "Phylum", drop = F]
Species_data <- data.frame(apply(Species_data, 2, as.factor))
colnames(Species_data) <- "Species"
rownames(Species_data) <- metadata2$Samples_ID
Species_colours <-c('#d4c29e','#de5275',"#88CCEE", "#DDCC77", "#117733", "#332288", "#AA4499", 
                    "#44AA99", "#999933")  ##rainbow(length(unique(metadata2$Species)), alpha = 1)
names(Species_colours) <- c(sort(unique(metadata2$Phylum)))

Species_tree = gheatmap(gg_mcmc_tree, Species_data,
                   width = 0.1,
                   offset = 0.21, 
                   colnames_position = "top",
                   colnames_angle = 0.45, 
                   colnames_offset_y = 0.1,
                   hjust = 0,
                   font.size = 2)+ 
                    scale_fill_manual(breaks=names(Species_colours), 
                                      values=Species_colours, name="Species")


Species_tree                    

Soil_data <- metadata2[, "Treatment", drop = F]
colnames(Soil_data) <- "Treatment"
unique(Soil_data)
rownames(Soil_data) <- metadata$Samples_ID
Soil_colours <- c('#8e67a0','#e9ae8c','#528551')
names(Soil_colours) <-  c(sort(unique(metadata$Treatment)))


Soil_hm <- gheatmap(Species_tree, Soil_data,
                         width = 0.1,
                         # Increase offset
                         offset = 0.25, 
                         colnames_position = "bottom",
                         colnames_angle = 0.45, 
                         colnames_offset_y = 1,
                         hjust = 0,
                         font.size = 3) +
                         scale_fill_manual(breaks=names(c(Species_colours,Soil_colours)), 
                                           values=c(Species_colours,Soil_colours), 
                                           name=c('Species', "Treatment"))

Soil_hm

####Adjust the p-value
metadata2$Enhanced <- "Irrigated"
pval_thres <- 0
metadata2$Enhanced[which(metadata2$diff < pval_thres)] <- "Non-Irrigated"
metadata2$Enhanced <- metadata2$Enhanced %>% factor

Irrigation_data <- metadata2[, "Enhanced", drop = F]
colnames(Irrigation_data) <- "Enhanced"
unique(Irrigation_data)
rownames(Irrigation_data) <- metadata2$Samples_ID
Irrigation_colours <- c('#FF7F50',"#00CC7A")
names(Irrigation_colours) <-  c(sort(unique(metadata2$Enhanced)))

Irrigation_hm <- gheatmap(Soil_hm, Irrigation_data,
                           width = 0.1,
                           # Increase offset
                           offset = 0.29, 
                           colnames_position = "bottom",
                           colnames_angle = 0.45, 
                           colnames_offset_y = 1,
                           hjust = 0,
                           font.size = 3) +
  scale_fill_manual(breaks=names(c(Species_colours,Soil_colours,Irrigation_colours)), values=c(Species_colours,Soil_colours,Irrigation_colours), name=c('Treatment'))
  geom_hilight(data=dat, mapping=aes(node=id, fill=type),
               align="right",extend=0.2,show.legend=TRUE,alpha=0.3) 

Irrigation_hm

