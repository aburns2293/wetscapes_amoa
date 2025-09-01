################################################################################
#                                                                              #
# Functional and phylogenetic assignment of ammonia                            #
# oxidizing bacteria and potential comammox bacteria                           #
#                                                                              #
# Author: Anna Burns                                                           #
# Last edited: 28.07.2023                                                      #
# Last tested: 25.04.2024                                                      #
#                                                                              #
################################################################################

################################################################################
# Required packages                                                            #
################################################################################

library(tidyverse)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(tibble)
library(vegan)

################################################################################
# Data upload                                                                  #
################################################################################

absolute.abundance <- read.csv('Data/prok.absolute.abundances.asv.csv') %>% remove_rownames() %>% column_to_rownames(var = 'Sample')

asv <- read.csv('Data/prok.asv.csv')

tax <- read.csv('Data/prok.tax.csv')

sample <- read.csv('Data/prok.sample.csv') %>% remove_rownames() %>% column_to_rownames(var = 'Sample')

potential.comammox <- read.csv('Data/nob.comammox.hits.csv')
  # includes nitrate oxidizing bacteria that had a BLAST ID in the first 100 hits
  # with known ammonia oxidizing bacteria

################################################################################
# Functions                                                                    #
################################################################################

################################################################################
# Code                                                                         #
################################################################################

# ammonia oxidizing bacteria

# potential comammox bacteria

ggplot(potential.comammox) + geom_histogram(aes(x = PercID), binwidth = 0.5) +
  theme_classic() + 
  xlab('% ID Match') + 
  ggtitle('Comammox Blast Results', subtitle ='Suggested ID cutoff at 97.5%') +
  geom_vline(xintercept = 99, col = 'blue', linewidth = 2)

## separating all hits with a Blast ID Match > 97.5% and adding their taxonomic/sample data

comammox_99 <- potential.comammox %>% filter(PercID >= 99) 

comammox_99_list <- as.vector(comammox_99$ASV) %>% unique()

comammox.absolute <- absolute.abundance[,colnames(absolute.abundance) %in% comammox_99_list]

comammox.absolute <- t(comammox.absolute)

comammox.vec <- as.vector(row.names(comammox.absolute))

tax.comammox <- tax[tax$ID %in% comammox.vec,]

tax.comammox <- tax.comammox %>% remove_rownames() %>% column_to_rownames(var = 'ID')

sample.names.comammox <- as.vector(colnames(comammox.absolute))

qpcr.sample.comammox <- sample[rownames(sample) %in% sample.names.comammox,]

## make phyloseq for potential comammox

tax.comammox.m <- as.matrix(tax.comammox)

comammox.absolute.m <- data.matrix(comammox.absolute)

comammox.OTU <- otu_table(comammox.absolute, taxa_are_rows = TRUE)

comammox.TAX <- tax_table(tax.comammox.m)

comammox.INFO <- sample_data(qpcr.sample.comammox)

comammox.ps <- phyloseq(comammox.OTU, comammox.TAX, comammox.INFO)

comammox.ps.melt <- psmelt(comammox.ps)

## analysis of seasonal trends in potential comammox nobs

comammox.ps.present <- comammox.ps.melt %>% filter(Abundance != 0)

comammox.ps.present$season2 <- factor(comammox.ps.present$season2, levels = c('17-Apr', '17-Aug', '17-Nov',
                                              '18-Feb', '18-Apr', '18-Jun', 
                                              '18-Aug', '18-Oct', '18-Dec',
                                              '19-Feb', '19-Jun', '19-Oct'))

comammox.summary <- comammox.ps.melt %>% group_by(loc, season2, OTU) %>% summarise(mean = mean(Abundance), sd = sd(Abundance), n = n(), se = sd/sqrt(n))

# PW

ggplot(comammox.summary[comammox.summary$loc == 'PW',], aes(x = season2, y = mean, col = OTU)) +
  geom_point() + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean, color = OTU),width = 0.1, linewidth = 1) +
  geom_line(aes(group = OTU)) +
  theme_classic() +
  ylab('copies/g DW soil') +
  xlab('') +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=20,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        axis.title.y=element_text(margin = margin(l = 60)),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.7),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

kruskal.test(Abundance~season2, data = comammox.ps.melt[comammox.ps.melt$loc == 'PW',])

# PD
ggplot(comammox.summary[comammox.summary$loc == 'PD',], aes(x = season2, y = mean, col = OTU)) +
  geom_point() + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean, color = OTU),width = 0.1, linewidth = 1) +
  geom_line(aes(group = OTU)) +
  theme_classic() +
  ylab('copies/g DW soil') +
  xlab('') +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=20,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        axis.title.y=element_text(margin = margin(l = 60)),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

kruskal.test(Abundance~season2, data = comammox.ps.melt[comammox.ps.melt$loc == 'PD',])

# CW
ggplot(comammox.summary[comammox.summary$loc == 'CW',], aes(x = season2, y = mean, col = OTU)) +
  geom_point() + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean, color = OTU),width = 0.1, linewidth = 1) +
  geom_line(aes(group = OTU)) +
  theme_classic() +
  ylab('copies/g DW soil') +
  xlab('') +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=20,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        axis.title.y=element_text(margin = margin(l = 60)),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

# CD
ggplot(comammox.summary[comammox.summary$loc == 'CD',], aes(x = season2, y = mean, col = OTU)) +
  geom_point() + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean, color = OTU),width = 0.1, linewidth = 1) +
  geom_line(aes(group = OTU)) +
  theme_classic() +
  ylab('copies/g DW soil') +
  xlab('') +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=20,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        axis.title.y=element_text(margin = margin(l = 60)),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

# relative abundances

asv <- asv %>% remove_rownames() %>% column_to_rownames(var = 'ID')

asv.t <- t(asv)

rel.abundance <- decostand(asv.t, method = 'total')

comammox_99 <- potential.comammox %>% filter(PercID >= 99) 

comammox_99_list <- as.vector(comammox_99$ASV) %>% unique()

comammox.rel <- rel.abundance[,colnames(rel.abundance) %in% comammox_99_list]

max(comammox.rel)

comammox.rel.present <- comammox.rel[,colSums(comammox.rel) > 0]

totals <- rowSums(comammox.rel.present)
max(totals)
mean(totals)
sd(totals)/417
n <- n(totals)
se <- (sd(totals)/n(totals))
# relative abundances for drought period

relevant.samples <- c('X16S.18AprPW105', 'X16S.18AprPW205', 'X16S.18AprPW305', 
                      'X16S.18AprCW105', 'X16S.18AprCW205', 'X16S.18AprCW305',
                      'X16S.18JunPW105', 'X16S.18JunPW205', 'X16S.18JunPW305',
                      'X16S.18JunCW105', 'X16S.18JunCW205', 'X16S.18JunCW305',
                      'X16S.18AugPW105', 'X16S.18AugPW205', 'X16S.18AugPW305',
                      'X16S.18AugCW105', 'X16S.18AugCW205', 'X16S.18AugCW305',
                      'X16S.18OctPW105', 'X16S.18OctPW205', 'X16S.18OctPW305',
                      'X16S.18OctCW105', 'X16S.18OctCW205', 'X16S.18OctCW305',
                      'X16S.18DecPW105', 'X16S.18DecPW205', 'X16S.18DecPW305',
                      'X16S.18DecCW105', 'X16S.18DecCW205', 'X16S.18DecCW305',
                      'X16S.19FebPW105', 'X16S.19FebPW205', 'X16S.19FebPW305',
                      'X16S.19FebCW105', 'X16S.19FebCW205', 'X16S.19FebCW305'
                      )

drought.comammox.rel <- comammox.rel[rownames(comammox.rel) %in% relevant.samples,]

pw.comammox.rel <- drought.comammox.rel[str_detect(rownames(drought.comammox.rel), 'PW') == TRUE,]

max(pw.comammox.rel)
mean(pw.comammox.rel)

cw.comammox.rel <- drought.comammox.rel[str_detect(rownames(drought.comammox.rel), 'CW') == TRUE,]

max(cw.comammox.rel)
mean(cw.comammox.rel)

# comparison between aob, aoa and comammox

asv <- asv %>% column_to_rownames(var = 'ID')
tax <- tax %>% column_to_rownames(var = 'ID')

wetscapes.ps <- phyloseq(tax_table(as.matrix(tax)),
                         otu_table(data.matrix(asv), taxa_are_rows = TRUE),
                         sample_data(sample))

wetscapes.ps.rel <- transform_sample_counts(wetscapes.ps, function(x) x/sum(x))

ao.tax <- tax %>% filter(Order == 'Nitrosomonadales' |
                           Class == 'Nitrososphaeria' |
                           Order == 'Nitrospirales')

ao.ids <- as.vector(rownames(ao.tax))

ao.ps.rel <- subset_taxa(wetscapes.ps.rel, taxa_names(wetscapes.ps.rel) %in% ao.ids)

ao.ps.melt <- psmelt(ao.ps.rel)

ao.ps.melt <- ao.ps.melt %>% filter(Abundance > 0)

ao.ps.melt$AO_Group <- NA

ao.ps.melt[ao.ps.melt$Domain == 'Bacteria',]$AO_Group <- 'AOB'
ao.ps.melt[ao.ps.melt$Domain == 'Archaea',]$AO_Group <- 'AOA'
ao.ps.melt[ao.ps.melt$OTU %in% comammox_975_list,]$AO_Group <- 'Comammox'

ao.ps.melt <- ao.ps.melt %>% mutate(Abundance_Perc = Abundance * 100) %>% filter(depth == '05-10 cm', 
                                                                                 loc %in% c('PW', 'CW'),
                                                                                 season2 %in% c('18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'))

ao.ps.melt.sum <- ao.ps.melt %>% group_by(loc, season2, group, AO_Group) %>% summarise(Total_Abundance = sum(Abundance_Perc))
ao.ps.melt.mean <- ao.ps.melt.sum %>% group_by(loc, season2, AO_Group) %>% summarise(Mean_Abundance = mean(Total_Abundance),
                                                                                            sd = sd(Total_Abundance),
                                                                                            n = n(),
                                                                                            se = sd/sqrt(n))

ao.ps.melt.mean$season2 <- factor(ao.ps.melt.mean$season2, levels = c('18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'))
ao.ps.melt.mean$AO_Group <- factor(ao.ps.melt.mean$AO_Group, levels = c('AOB', 'AOA', 'Comammox'))

ggplot(ao.ps.melt.mean[ao.ps.melt.mean$AO_Group == 'Comammox',], aes(x = season2, y = Mean_Abundance)) + 
  geom_point(size = 6) + 
  geom_line(linewidth = 1.3, aes(group = 1)) +
  geom_errorbar(aes(ymax = Mean_Abundance + se, ymin = Mean_Abundance), width = 0.1, linewidth = 1) + 
  xlab('') + 
  ylab('Relative Abundance [%]') + 
  scale_y_continuous(limits = c(0, 0.65)) +
  ggtitle('Relative abundance of potential comammox ASVs, PW') + 
  theme_classic() + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 20))

comammox.melt <- ao.ps.melt %>% filter(AO_Group == 'Comammox')
comammox.melt$OTU %>% unique()
mean(comammox.melt$Abundance_Perc)
sd(comammox.melt$Abundance_Perc)
