################################################################################
#                                                                              #
# Assigning AOA Clades Based on amoA Database                                  #
#                                                                              #
# Author: Anna Burns                                                           #
# Last edited: 19.04.2024                                                      #
#                                                                              #
################################################################################

################################################################################
# Required packages                                                            #
################################################################################

library(dplyr)
library(rstatix)
library(ggplot2)
library(phyloseq)

################################################################################
# Data upload                                                                  #
################################################################################

rel.ab <- read.csv('Data/prok.asv.csv')

sample <- read.csv('../wetscapes_amoa/Data/prok.sample.csv')

tax <- read.csv('../wetscapes_amoa/Data/prok.tax.csv', row.names = 1)

aoa.clades <- read.csv('../wetscapes_amoa/Data/16SrRNA.aoa.clades.csv')

aoa.abs <- read_csv("PROK Analysis/qPCR_mRNA/absolute.abundances.csv")

################################################################################
# Code                                                                         #
################################################################################

# this code is for reproducibility purposes; the data is already saved under:
# ./wetscapes_amoa/Data/16srRNA.aoa.clades.csv

## create tidy dataframe of absolute abundance data (16S rRNA)

# abs.ab <- absolute_abundances %>% remove_rownames() %>% column_to_rownames(var = '...1')
# 
# sample <- sample %>% select(-X) %>% remove_rownames() %>% column_to_rownames(var = 'Sample')
# 
# asv.m <- t(data.matrix(abs.ab))
# 
# tax.m <- as.matrix(tax)
# 
# OTU <- otu_table(asv.m, taxa_are_rows = TRUE)
# 
# TAX <- tax_table(tax.m)
# 
# INFO <- sample_data(sample)
# 
# ps <- phyloseq(OTU, TAX, INFO)
# 
# ps.melt <- psmelt(ps)
# 
# ## filter tidy dataframe by AOA ASVs
# 
# aoa.clades <- aoa.clades %>% select(X, Assigned_Clade_Tree2)
# aoa.asvs <- as.vector(aoa.clades$X)
# 
# rel.ab.aoa <- ps.melt[ps.melt$OTU %in% aoa.asvs,]
# 
# aoa.full <- dplyr::left_join(rel.ab.aoa, aoa.clades, by = join_by('OTU' == 'X'))

# write.csv(aoa.full, '../wetscapes_amoa/Data/16srRNA.aoa.absolute.csv')

# clean absolute data results for 2018 drought cycle

aoa.abs.2018 <- aoa.abs %>% filter(depth == '05-10 cm',
                                   year == '18' | season2 == '19-Feb',
                                   season2 != '18-Feb',
                                   loc == 'PW' | loc == 'CW') 

aoa.abs.2018$season2 <- factor(aoa.abs.2018$season2, levels = c('18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'))

################################################################################
# Figures and Statistics                                                       #
################################################################################

aoa.summary <- aoa.abs.2018 %>% group_by(loc, season2, OTU) %>% mutate(mean = mean(Abundance), sd = sd(Abundance), n = n(), se = sd/sqrt(n))

aoa.present <- aoa.present %>% filter(mean != 0)

aoa.calc <- aoa.present %>% group_by(loc, season2, Assigned_Clade_Tree2) %>% summarize(mean = mean(Abundance), sd = sd(Abundance), n = n(), se = sd/sqrt(n))

ggplot(aoa.calc, aes(x = season2, y = Assigned_Clade_Tree2, shape = loc)) +
  geom_point(aes(size = mean)) + 
  scale_shape_manual(values = c(1, 16)) + 
  scale_size_continuous(range = c(5,15), breaks = c(1e6, 5e6, 1e7)) + 
  xlab('') +
  ylab('') + 
  theme_classic() + 
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=15,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(), 
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

aoa.calc2 <- aoa.present %>% select(Assigned_Clade_Tree2, loc, season2, Abundance)

rstatix::levene_test(aoa.calc2[aoa.calc2$OTU == 'Dada1987',], Abundance~season2)

aoa.rel.drought.shallow[aoa.rel.drought.shallow$loc == 'PW' & aoa.rel.drought.shallow$Assigned_Clade_Tree2 == 'UD',] %>% shapiro_test(Abundance)
aoa.rel.drought.shallow[aoa.rel.drought.shallow$loc == 'PW' & aoa.rel.drought.shallow$Assigned_Clade_Tree2 == 'UD',] %>% levene_test(Abundance ~ season2)

aoa.rel.drought.shallow[aoa.rel.drought.shallow$loc == 'PW' & aoa.rel.drought.shallow$Assigned_Clade_Tree2 == 'NP-Eta',] %>% kruskal_test(Abundance ~ season2)

