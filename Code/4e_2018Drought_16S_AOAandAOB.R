################################################################################
#                                                                              #
# 16S rRNA relative abundance of AOA and AOB ASVs                              #
#                                                                              #
# Author: Anna Burns                                                           #
# Last edited: 04.06.2024                                                      #
# Last tested: 04.06.2024                                                      #
#                                                                              #
################################################################################

################################################################################
# Required packages                                                            #
################################################################################

library(dplyr)
library(tidyverse)
library(phyloseq)

# to install phyloseq (do not run if already installed)
# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install(version = "3.18")
# BiocManager::install('phyloseq')

################################################################################
# Data upload                                                                  #
################################################################################

aoa.16s <- read.csv('Data/16srRNA.aoa.csv')
aob.16s <- read.csv('Data/16srRNA.aob.csv')

################################################################################
# Code                                                                         #
################################################################################

asv <- read.csv('Data/prok.absolute.abundances.asv.csv')
sample <- read.csv('Data/prok.sample.csv')
tax <- read.csv('Data/prok.tax.csv')

sample.2018drought <- sample %>% filter(season2 == '18-Apr' | season2 == '18-Jun' | season2 == '18-Aug' | season2 == '18-Oct' | 
                                          season2 == '18-Dec' |season2 == '19-Feb') %>% filter(depth == '05-10 cm')

sample.formatted <- sample.2018drought %>% remove_rownames() %>% column_to_rownames(var = 'Sample')

drought.ids <- as.vector(sample.2018drought$Sample)
drought.asv <- asv[, asv$Sample %in% drought.ids]

drought.asv.formatted <- drought.asv %>% remove_rownames() %>% column_to_rownames(var = 'Sample')
drought.asv.formatted.t <- t(drought.asv.formatted)

# format aoa data (saved as Data/16srRNA.aoa.csv)

aoa.tax <- tax %>% filter(Class == 'Nitrososphaeria')
aoa.ids <- as.vector(aoa.tax$ID)
aoa.tax <- aoa.tax %>% remove_rownames() %>% column_to_rownames(var = 'ID')
aoa.asv <- drought.asv.formatted.t[rownames(drought.asv.formatted.t) %in% aoa.ids, ]

aoa.tax.ps <- tax_table(as.matrix(aoa.tax))
aoa.asv.ps <- otu_table(as.matrix(aoa.asv), taxa_are_rows = TRUE)
sample.ps <- sample_data(sample.formatted)

aoa.ps <- phyloseq(aoa.tax.ps, aoa.asv.ps, sample.ps)
aoa.ps.melt <- psmelt(aoa.ps)

aoa.ps.present <- aoa.ps.melt %>% filter(Abundance > 0)
write.csv(aoa.ps.present, 'Data/16srRNA.aoa.csv')

# format aob data (saved as Data/16srRNA.aob.csv)

aob.tax <- tax %>% filter(Class == 'Nitrospira' | Order == 'Nitrosomonadales')
aob.ids <- as.vector(aob.tax$ID)
aob.tax <- aob.tax %>% remove_rownames() %>% column_to_rownames(var = 'ID')
aob.asv <- drought.asv.formatted.t[rownames(drought.asv.formatted.t) %in% aob.ids, ]

aob.tax.ps <- tax_table(as.matrix(aob.tax))
aob.asv.ps <- otu_table(as.matrix(aob.asv), taxa_are_rows = TRUE)
sample.ps <- sample_data(sample.formatted)

aob.ps <- phyloseq(aob.tax.ps, aob.asv.ps, sample.ps)
aob.ps.melt <- psmelt(aob.ps)

aob.ps.present <- aob.ps.melt %>% filter(Abundance > 0)
write.csv(aob.ps.present, 'Data/16srRNA.aob.csv')

################################################################################
# Figures and statistics                                                       #
################################################################################

aoa.summary <- aoa.16s %>% filter(loc == 'PW' | loc == 'CW') %>% group_by(loc, season2, group) %>% summarise(AOA = sum(Abundance))
aoa.avg <- aoa.summary %>% group_by(loc, season2) %>% summarise(aoa.mean = mean(AOA), aoa.sd = sd(AOA), n.aoa = n(), se.aoa = aoa.sd/sqrt(n.aoa)) 

aob.summary <- aob.16s %>% filter(loc == 'PW' | loc == 'CW') %>% group_by(loc, season2, group) %>% summarise(AOB = sum(Abundance))

summary.16s <- left_join(aob.summary, aoa.summary, by = c('loc', 'season2'))
summary.16s$AOA <- summary.16s$AOA %>% replace_na(0)

summary.long <- pivot_longer(summary.16s, cols = !c('loc', 'season2'), names_to = 'Domain', values_to = 'Abundance')

summary.long$season2 <- factor(summary.long$season2, levels = c('18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'), 
                               labels = c('2018-04', '2018-06', '2018-08', '2018-10', '2018-12', '2019-02'))


ggplot(summary.long[summary.long$loc == 'PW',], aes(x = season2, y = Abundance, group = Domain, col = Domain)) +
  geom_point(aes(shape = Domain)) +
  geom_line()
