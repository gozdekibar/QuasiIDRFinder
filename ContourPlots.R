library(MASS)
library(viridis)
library(ggplot2)
library(gridExtra)  ## for grid.arrange()
library(tidyverse)
library(RColorBrewer)

library("readxl")

#setwd("/Users/kibar/Desktop/QuasiIDRFinder/")

my_data_metapredict <- read_excel("./Dataset/MasterTable_canonicalHuman_qIDR.xlsx")

thre=0.001


## 1) top qIDRs
my_data_metapredict_top=my_data_metapredict[my_data_metapredict$`min(p-val)`< 0.05, ]
my_data_metapredict_top$`min(p-val)`[which(my_data_metapredict_top$`min(p-val)`< thre)]=thre

highlight_df_3 <- my_data_metapredict_top %>% 
  filter(my_data_metapredict_top$overlap_metapredict  =="yes" )

blues <- brewer.pal(9, "Blues")
blue_range <- colorRampPalette(blues)

ratio <- with(my_data_metapredict_top, diff(range(-log10(`min(p-val)`)))/diff(c(90,300)))
ggplot(my_data_metapredict_top, aes(x=-log10(`min(p-val)`), y= `length_IDR(qIDR)`) ) + 
  geom_density2d_filled(bins=9) + geom_point(data=highlight_df_3,  aes(x=-log10(`min(p-val)`),y=`length_IDR(qIDR)`), 
                                             color='black', size=0.4)+ theme_classic() +ylim(90,300)+scale_fill_manual(values = blue_range(11))+ggtitle("top qIDRs with metapredict overlaps")+ coord_fixed(ratio=ratio)+theme(legend.title = element_text( size=9), legend.text=element_text(size=7))+xlab(expression(-log[10]~(p-value)))+ylab("Length (amino acids)") # for the x axis label




##2)all qIDRs  
my_data_all_qIDRs=my_data_metapredict  
my_data_all_qIDRs$`min(p-val)`[which(my_data_all_qIDRs$`min(p-val)`< thre)]=thre


keys <- read.csv(file = './Dataset/HG30_pep_annotation.csv')
TFnames <- read.table("./Dataset/annotations/TF_IDs.txt",header = F)
PLDs <- read.table("./Dataset/annotations/PLD_IDs.txt",header = F)
Pld_99 <- read_excel("./Dataset/annotations/PLD_99p.xlsx")

#convert PLD keys to match
PLD_keys=keys$gene_symbol[match(PLDs$V1, keys$ID )]


highlight_df_TFs <- my_data_all_qIDRs %>% 
  filter(my_data_all_qIDRs$gene_symbol  %in% TFnames$V1 )


highlight_df_PLDs <- my_data_all_qIDRs %>% 
  filter(my_data_all_qIDRs$gene_symbol  %in% PLD_keys )



## 99 peridic PLDs


highlight_df_PLDs99 <- my_data_all_qIDRs %>% 
  filter(my_data_all_qIDRs$gene_symbol  %in% Pld_99$Name )



ratio <- with(my_data_all_qIDRs, diff(range(-log10(`min(p-val)`)))/diff(c(90,300)))

##all qIDRs density
ggplot(my_data_all_qIDRs, aes(x=-log10(`min(p-val)`), y= `length_IDR(qIDR)`) ) + 
  geom_density2d_filled(bins=9) +ylim(90,300)+scale_fill_manual(values = blue_range(11))+ theme_classic()+ggtitle("all IDRs with significant periodicity")+ coord_fixed(ratio=ratio)+theme(legend.title = element_text( size=9), legend.text=element_text(size=7))+xlab(expression(-log[10]~(p-value)))+ylab("Length (amino acids)") # for the x axis label



greens <- brewer.pal(9, "Greens")
greens_range <- colorRampPalette(greens)

ratio <- with(highlight_df_TFs, diff(range(-log10(`min(p-val)`)))/diff(c(90,300)))

ggplot(highlight_df_TFs, aes(x=-log10(`min(p-val)`), y= `length_IDR(qIDR)`) ) + 
  geom_density2d_filled(bins=9) + theme_classic() +ylim(90,300)+scale_fill_manual(values = greens_range(11))+ggtitle("Transcription Factors")+ coord_fixed(ratio=ratio)+theme(legend.title = element_text( size=9), legend.text=element_text(size=7))+xlab(expression(-log[10]~(p-value)))+ylab("Length (amino acids)") # for the x axis label




purples <- brewer.pal(9, "Purples")
purple_range <- colorRampPalette(purples)

ratio <- with(highlight_df_PLDs, diff(range(-log10(`min(p-val)`)))/diff(c(90,300)))


ggplot(highlight_df_PLDs, aes(x=-log10(`min(p-val)`), y= `length_IDR(qIDR)`) ) + 
  geom_density2d_filled(bins=9) + theme_classic() +ylim(90,300)+scale_fill_manual(values = purple_range(11))+ggtitle("Prion-like domains")+ coord_fixed(ratio=ratio)+theme(legend.title = element_text( size=9), legend.text=element_text(size=7))+xlab(expression(-log[10]~(p-value)))+ylab("Length (amino acids)") # for the x axis label




reds <- brewer.pal(9, "Reds")
red_range <- colorRampPalette(reds)

ratio <- with(highlight_df_PLDs99, diff(range(-log10(`min(p-val)`)))/diff(c(90,300)))


ggplot(highlight_df_PLDs99, aes(x=-log10(`min(p-val)`), y= `length_IDR(qIDR)`) ) + 
  geom_density2d_filled(bins=9) + theme_classic() +ylim(90,300)+scale_fill_manual(values = red_range(11))+ggtitle("Aromatic-rich prion-like domains")+ coord_fixed(ratio=ratio)+theme(legend.title = element_text( size=9), legend.text=element_text(size=7))+xlab(expression(-log[10]~(p-value)))+ylab("Length (amino acids)") # for the x axis label




