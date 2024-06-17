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
library(rstatix)

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
# Functions                                                                    #
################################################################################

gene.test <- function(site, Data) {
  
  distrib <- levene_test(formula = Abundance~season2, data = Data[Data$loc == site,])
  print(paste0('Distribution (Levene) p = ', round(distrib$p, digits = 3)))
  
  norm <- shapiro.test(Data[Data$loc == site,]$Abundance)
  print(paste0('Normality (Shapiro-Wilk) p = ', round(norm$p, digits = 3)))
  
  aov <- summary(aov(Abundance~season2, Data[Data$loc == site,]))
  print(paste0('Anova p = ', round(aov[[1]][[5]][[1]], digits = 3)))
  
  tukey <- TukeyHSD(aov(Abundance~season2, Data[Data$loc == site,]))
  print('Tukey HSD output:')
  print(tukey[['season2']])
}
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
aoa.avg <- aoa.summary %>% group_by(loc, season2) %>% summarise(mean = mean(AOA), sd = sd(AOA), n = n(), se = sd/sqrt(n)) 
aoa.avg$Domain <- 'AOA'

aob.summary <- aob.16s %>% filter(loc == 'PW' | loc == 'CW') %>% group_by(loc, season2, group) %>% summarise(AOB = sum(Abundance))
aob.avg <- aob.summary %>% group_by(loc, season2) %>% summarise(mean = mean(AOB), sd = sd(AOB), n = n(), se = sd/sqrt(n)) 
aob.avg$Domain <- 'AOB'

summary.16s <- rbind(aoa.avg, aob.avg)

summary.16s$season2 <- factor(summary.16s$season2, levels = c('18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'), 
                               labels = c('2018-04', '2018-06', '2018-08', '2018-10', '2018-12', '2019-02'))


ggplot(summary.16s[summary.16s$loc == 'PW',], aes(x = season2, y = mean, group = Domain, col = Domain)) +
  geom_point(aes(shape = Domain), size = 10) +
  geom_line(linewidth = 1.3) +
  theme_classic() + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean, color = Domain), width = 0.1, linewidth = 1) +
  scale_shape_manual(values = c(15,16,17)) + 
  scale_color_manual(values = c("#A04000","#117A65")) + 
  scale_y_continuous(limits = c(0, 3.5e9), breaks = c(0, 1e9, 2e9, 3e9), labels = c('0.0','1.0', '2.0','3.0')) +
  scale_x_discrete(drop = FALSE) + 
  labs(y = expression(paste(x10^{9}, ' copies/g DW soil')), x = '') + 
  theme(axis.title=element_text(size=20,face = "bold"),
        axis.title.y=element_text(margin = margin(l = 60)),
        axis.text=element_text(size=20,face = "bold"),
        axis.text.x = element_text(angle=45, vjust = 0.6),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

gene.test(Data = aoa.16s, site = 'PW')
kruskal.test(Abundance~season2, data = aoa.16s[aoa.16s$loc == 'PW',])
dunn_test(aoa.16s[aoa.16s$loc == 'PW',], Abundance~season2, p.adjust.method = 'bonferroni')

gene.test(Data = aob.16s, site = 'PW')
kruskal.test(Abundance~season2, data = aob.16s[aob.16s$loc == 'PW',])
dunn_test(aob.16s[aob.16s$loc == 'PW',], Abundance~season2, p.adjust.method = 'bonferroni')

ggplot(summary.16s[summary.16s$loc == 'CW',], aes(x = season2, y = mean, group = Domain, col = Domain)) +
  geom_point(aes(shape = Domain), size = 10) +
  geom_line(linewidth = 1.3) +
  theme_classic() + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean, color = Domain), width = 0.1, linewidth = 1) +
  scale_shape_manual(values = c(15,16,17)) + 
  scale_color_manual(values = c("#A04000","#117A65")) + 
  scale_y_continuous(limits = c(0, 1.2e9), breaks = c(0, 3e8, 6e8, 9e8, 12e8), labels = c('0.0','3.0', '6.0','9.0', '12.0')) +
  scale_x_discrete(drop = FALSE) + 
  labs(y = expression(paste(x10^{8}, ' copies/g DW soil')), x = '') + 
  theme(axis.title=element_text(size=20,face = "bold"),
        axis.title.y=element_text(margin = margin(l = 60)),
        axis.text=element_text(size=20,face = "bold"),
        axis.text.x = element_text(angle=45, vjust = 0.6),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

gene.test(Data = aob.16s, site = 'CW')
kruskal.test(Abundance~season2, data = aob.16s[aob.16s$loc == 'CW',])
dunn_test(aob.16s[aob.16s$loc == 'CW',], Abundance~season2, p.adjust.method = 'bonferroni')

