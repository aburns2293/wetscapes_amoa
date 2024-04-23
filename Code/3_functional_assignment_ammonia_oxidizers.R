################################################################################
#                                                                              #
# Functional and phylogenetic assignment of ammonia                            #
# oxidizing bacteria and potential comammox bacteria                           #
#                                                                              #
# Author: Anna Burns                                                           #
# Last edited: 28.07.2023                                                      #
#                                                                              #
################################################################################

################################################################################
# Required packages                                                            #
################################################################################

library(tidyverse)
library(dplyr)
library(ggplot2)
library(phyloseq)

################################################################################
# Data upload                                                                  #
################################################################################

absolute.abundance <- read.csv('Data/prok.sample.csv') %>% remove_rownames() %>% column_to_rownames(var = 'X')

tax <- read.csv('Data/prok.tax.csv')

sample <- read.csv('Data/prok.sample.csv') %>% select(-X) %>% remove_rownames() %>% column_to_rownames(var = 'Sample')

potential.comammox <- read.csv('Data/nob.comammox.hits.csv') %>% select(-X)
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
  geom_vline(xintercept = 97.5, col = 'blue', linewidth = 2)

## separating all hits with a Blast ID Match > 97.5% and adding their taxonomic/sample data

comammox_975 <- potential.comammox %>% filter(PercID >= 97.5) 

comammox_975_list <- as.vector(comammox_975$ASV) %>% unique()

comammox.absolute <- absolute.abundance[,colnames(absolute.abundance) %in% comammox_975_list]

comammox.absolute <- t(comammox.absolute)

comammox.vec <- as.vector(row.names(comammox.absolute))

tax.comammox <- tax[tax$X %in% comammox.vec,]

tax.comammox <- tax.comammox %>% remove_rownames() %>% column_to_rownames(var = 'X')

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

comammox.ps.present <- comammox_absolute %>% filter(Abundance != 0)

comammox.ps.present$season2 <- factor(comammox.ps.present$season2, levels = c('17-Apr', '17-Aug', '17-Nov',
                                              '18-Feb', '18-Apr', '18-Jun', 
                                              '18-Aug', '18-Oct', '18-Dec',
                                              '19-Feb', '19-Jun', '19-Oct'))

comammox.summary <- comammox_absolute %>% group_by(loc, season2, OTU) %>% summarise(mean = mean(Abundance), sd = sd(Abundance), n = n(), se = sd/sqrt(n))

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
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

kruskal.test(Abundance~season2, data = comammox_absolute[comammox_absolute$loc == 'PW',])

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

kruskal.test(Abundance~season2, data = comammox_absolute[comammox_absolute$loc == 'PD',])

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
