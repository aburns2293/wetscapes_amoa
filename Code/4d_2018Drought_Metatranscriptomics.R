################################################################################
#                                                                              #
# Analysis of metatranscriptomic dataset for PW + CW, 04.2018 - 02.2019        #
#                                                                              #
# Author: Anna Burns                                                           #
# Last edited: 15.01.2024                                                      #
# Last tested: 25.04.2024                                                      #
#                                                                              #
################################################################################

################################################################################
# Required packages                                                            #
################################################################################

library(dplyr)
library(stringr)
library(tidyverse)
library(RColorBrewer)
library(lubridate)

################################################################################
# Data upload                                                                  #
################################################################################

aoa <- read.csv('Data/mRNA.aoa.absolute.csv')

aob <- read.csv('Data/mRNA.aob.absolute.csv')

################################################################################
# Functions                                                                    #
################################################################################

gene.test <- function(site, Data, subunit) {
  
  distrib <- levene_test(formula = Amount~month, data = Data[Data$Site == site & Data$Subunit == subunit,])
  print(paste0('Distribution (Levene) p = ', round(distrib$p, digits = 3)))
  
  norm <- shapiro.test(Data[Data$Site == site & Data$Subunit == subunit,]$Amount)
  print(paste0('Normality (Shapiro-Wilk) p = ', round(norm$p, digits = 3)))
  
  aov <- summary(aov(Amount ~ month, Data[Data$Site == site & Data$Subunit == subunit,]))
  print(paste0('Anova p = ', round(aov[[1]][[5]][[1]], digits = 3)))
  
  tukey <- TukeyHSD(aov(Amount ~ month, Data[Data$Site == site & Data$Subunit == subunit,]))
  print('Tukey HSD output:')
  print(tukey[['month']])
}

################################################################################
# Code                                                                         #
################################################################################

# preparing aoa data

aoa$month <- factor(aoa$month, levels = c('18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'),
                    labels = c('2018-04', '2018-06', '2018-08', '2018-10', '2018-12', '2019-02'))

aoa_summary <- aoa %>% group_by(Site, month, Subunit) %>% summarise(mean = mean(Amount), sd = sd(Amount), n = n(), se = sd/sqrt(n))

aoa_summary$month <- factor(aoa_summary$month, levels = c('18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'), 
                              labels = c('2018-04', '2018-06', '2018-08', '2018-10', '2018-12', '2019-02'))

# preparing aob data

aob$month <- factor(aob$month, levels = c('18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'))

aob_summary <- aob %>% group_by(Site, month, Subunit) %>% summarise(mean = mean(Amount), sd = sd(Amount), n = n(), se = sd/sqrt(n))

aob_summary$month <- factor(aob_summary$month, levels = c('18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'), 
                            labels = c('2018-04', '2018-06', '2018-08', '2018-10', '2018-12', '2019-02'))

################################################################################
# Figures and Statistics                                                       #
################################################################################

# aoa

## coastal wet site

ggplot(aoa_summary[aoa_summary$Site == 'CW',], aes(x = month, y = mean, col = Subunit, group = Subunit)) +
  geom_point(aes(shape = Subunit), size = 10) + 
  geom_line(linewidth = 1.3) + 
  theme_classic() + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean, color = Subunit), width = 0.1, linewidth = 1) +
  scale_shape_manual(values = c(15,16,17)) + 
  scale_color_manual(values = c("#A04000","#117A65","#7D3C98")) + 
  scale_y_continuous(limits = c(0, 1e7), breaks = c(0, 3e6, 6e6, 9e6), labels = c('0.0','3.0', '6.0','9.0')) +
  scale_x_discrete(drop = FALSE) + 
  labs(y = expression(paste(x10^{6}, ' copies/g DW soil')), x = '') + 
  theme(axis.title=element_text(size=20,face = "bold"),
        axis.title.y=element_text(margin = margin(l = 60)),
        axis.text.x = element_text(angle=45, vjust = 0.6),
        axis.text=element_text(size=20,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.85,0.85),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

gene.test('CW', aoa, 'amo-A')
kruskal.test(Amount~month, data = aoa[aoa$Site == 'CW' & aoa$Subunit == 'amo-A',])
dunn_test(aoa[aoa$Site == 'CW' & aoa$Subunit == 'amo-A',], Amount~month, p.adjust.method = 'bonferroni')

## percolation wet site

ggplot(aoa_summary[aoa_summary$Site == 'PW',], aes(x = month, y = mean, col = Subunit, group = Subunit)) +
  geom_point(aes(shape = Subunit), size = 10) + 
  geom_line(size = 1.3) + 
  theme_classic() + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean, color = Subunit), width = 0.1, linewidth = 1) +
  scale_shape_manual(values = c(15,16,17)) + 
  scale_color_manual(values = c("#A04000","#117A65","#7D3C98")) + 
  scale_y_continuous(limits = c(0, 1.7e8), breaks = c(0, 5e7, 1e8, 1.5e8), labels = c('0.0','0.5','1.0','1.5')) +
  scale_x_discrete(drop = FALSE) + 
  labs(x = '', y = expression(paste(x10^{8}, ' copies/g DW soil'))) +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=20,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        axis.text.x = element_text(angle=45, vjust = 0.6),
        axis.title.y=element_text(margin = margin(l = 60)),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.85,0.85),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

gene.test('PW', aoa, 'amo-C')
kruskal.test(Amount~month, data = aoa[aoa$Site == 'PW' & aoa$Subunit == 'amo-C',])
dunn_test(aoa[aoa$Site == 'PW' & aoa$Subunit == 'amo-C',], Amount~month, p.adjust.method = 'bonferroni')

# aob

# coastal wet site

ggplot(aob_summary[aob_summary$Site == 'CW',], aes(x = month, y = mean, col = Subunit, group = Subunit)) +
  geom_point(aes(shape = Subunit), size = 10) + 
  geom_line(size = 1.3) + 
  theme_classic() + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean, color = Subunit), width = 0.1, linewidth = 1) +
  scale_shape_manual(values = c(15,16,17)) + 
  scale_color_manual(values = c("#117A65","#7D3C98")) + 
  scale_y_continuous(limits = c(0, 2.2e6), breaks = c(0, 0.5e6, 1e6, 1.5e6,2e6), labels = c('0.0','0.5','1.0','1.5', '2.0')) +
  scale_x_discrete(drop = FALSE) + 
  labs(x = '', y = expression(paste(x10^{6}, ' copies/g DW soil'))) +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=20,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        axis.text.x = element_text(angle=45, vjust = 0.6),
        axis.title.y=element_text(margin = margin(l = 60)),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.85,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

gene.test('CW', aob, 'amo-C')
kruskal.test(Amount~month, data = aob[aob$Site == 'CW' & aob$Subunit == 'amo-B',])
dunn_test(aob[aob$Site == 'CW' & aob$Subunit == 'amo-B',], Amount~month, p.adjust.method = 'bonferroni')

# percolation wet site

ggplot(aob_summary[aob_summary$Site == 'PW',], aes(x = month, y = mean, col = Subunit, group = Subunit)) +
  geom_point(aes(shape = Subunit), size = 10) + 
  geom_line(size = 1.3) + 
  theme_classic() + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean, color = Subunit), width = 0.1, linewidth = 1) +
  scale_shape_manual(values = c(15,16,17)) + 
  scale_color_manual(values = c("#117A65","#7D3C98")) + 
  scale_y_continuous(limits = c(0,7e6), breaks = c(0, 3e6, 6e6), labels = c('0.0','3.0','6.0')) +
  scale_x_discrete(drop = FALSE) + 
  labs(x = '', y = expression(paste(x10^{6}, ' copies/g DW soil'))) +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=20,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        axis.text.x = element_text(angle=45, vjust = 0.6),
        axis.title.y=element_text(margin = margin(l = 60)),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.85,0.91),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

gene.test('PW', aob, 'amo-B')
kruskal.test(Amount~month, data = aob[aob$Site == 'PW' & aob$Subunit == 'amo-B',])
dunn_test(aob[aob$Site == 'PW' & aob$Subunit == 'amo-B',], Amount~month, p.adjust.method = 'bonferroni')

# summary statistics

aoa$Domain <- 'AOA'
aob$Domain <- 'AOB'

## aoa vs. aob in PW

pw.aom <- rbind(aoa[aoa$Site == 'PW',], aob[aob$Site == 'PW',]) 
levene_test(pw.aom, Amount ~ Domain)
kruskal.test(pw.aom, Amount ~ Domain)

## aoa vs. aob in CW

cw.aom <- rbind(aoa[aoa$Site == 'CW',], aob[aob$Site == 'CW',]) 
levene_test(cw.aom, Amount ~ Domain)
shapiro_test(cw.aom$Amount)
kruskal.test(cw.aom, Amount ~ Domain)

## CW aoa vs. PW aoa

levene_test(aoa, Amount ~ Site)
shapiro_test(aoa$Amount, aoa)
kruskal.test(aoa[aoa$Subunit == 'amo-C',], Amount~Site)
dunn_test(aoa[aoa$Subunit == 'amo-C',], Amount~Site, p.adjust.method = 'bonferroni')

