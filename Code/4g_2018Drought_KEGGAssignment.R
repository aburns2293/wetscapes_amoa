################################################################################
#                                                                              #
# KEGG assignment and gene analysis of 2018 metatranscriptomic data            #
#                                                                              #
# Author: Anna Burns                                                           #
# Last edited: 19.04.2024                                                      #
#                                                                              #
################################################################################

################################################################################
# Required packages                                                            #
################################################################################

library(tidyverse)
library(dplyr)
library(ggplot2)
library(stringr)
library(rstatix)

################################################################################
# Data upload                                                                  #
################################################################################

mrna.tax <- read.csv('./Data/kegg.nitrogencycle.gene.tax.csv')

kegg <- read.csv('./Data/kegg.nitrogencycle.gene.abundance.csv')

mrna.sample <- read.csv('./Data/kegg.nitrogencycle.sample.csv')
mrna.sample$ID2 <- str_c(mrna.sample$sample, mrna.sample$month, sep = '_')

################################################################################
# Functions                                                                    #
################################################################################

gene.test <- function(site, gene, Data) {
  
  distrib <- levene_test(formula = Abundance~Date, data = Data[Data$Site == site & Data$Gene == gene,])
  print(paste0('Distribution (Levene) p = ', round(distrib$p, digits = 3)))
  
  norm <- shapiro.test(Data[Data$Site == site & Data$Gene == gene,]$Abundance)
  print(paste0('Normality (Shapiro-Wilk) p = ', round(norm$p, digits = 3)))
  
  aov <- summary(aov(Abundance~Date, Data[Data$Site == site & Data$Gene == gene,]))
  print(paste0('Anova p = ', round(aov[[1]][[5]][[1]], digits = 3)))
  
  tukey <- TukeyHSD(aov(Abundance~Date, Data[Data$Site == site & Data$Gene == gene,]))
  print('Tukey HSD output:')
  print(tukey[['Date']])
}

################################################################################
# Code                                                                         #
################################################################################

################################################################################
# Figures and Statistics                                                       #
################################################################################

# nitrogen fixation

kegg <- kegg %>% remove_rownames() %>% column_to_rownames(var = 'X')

nit.fix <- kegg %>% select(KO1381, KO1382, KO1383, KO1368, KO1402, KO1403, KO1401)

nit.fix <- nit.fix %>% rownames_to_column(var = 'ID')

nit.fix$Sample <- str_sub(nit.fix$ID, 1, 3)
nit.fix$Site <- str_sub(nit.fix$ID, 1, 2)
nit.fix$Date <- str_sub(nit.fix$ID, 5)
nit.fix$Date <- factor(nit.fix$Date, levels = c('18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'))

ggplot(nit.fix, aes(x = Date, y = KO1382, col = Site, group = Site)) + geom_point()

nifh.sum <- nit.fix %>% group_by(Site, Date) %>% summarise(mean = mean(KO1382), sd = sd(KO1382), n = n(), se = sd/sqrt(n))

ggplot(nifh.sum, aes(x = Date, y = mean, col = Site, group = Site)) + 
  geom_point(aes(shape = Site), size = 10) + 
  geom_line(size = 1.3) + 
  scale_color_manual(values = c("#A04000","#117A65","#7D3C98")) + 
  labs(y = expression(paste(x10^{8}, ' nifH transcripts/g DW soil')), x = '') + 
  scale_y_continuous(limits = c(0, 2.3e8), breaks = c(0, 0.5e8, 1e8, 1.5e8, 2e8), labels = c('0.0', '0.5', '1.0', '1.5', '2.0')) +
  scale_x_discrete(drop = FALSE) + 
  scale_shape_manual(values = c(15,16,17)) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean), width = 0.1, size = 1) +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=15,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

nit.fix.long <- pivot_longer(nit.fix, cols = !c('ID', 'Sample', 'Site', 'Date'), names_to = 'Gene', values_to = 'Abundance')

levene_test(Abundance~Site, data = nit.fix.long)
kruskal.test(Abundance~Site, data = nit.fix.long)

gene.test('CW', 'KO1381', nit.fix.long)

kruskal_test(KO1381~Date, data = nit.fix[nit.fix$Site == 'CW',])
dunn_test(KO1381~Date, data = nit.fix[nit.fix$Site == 'CW',], p.adjust.method = 'bonferroni')

ggplot(nit.fix.long[nit.fix.long$Site == 'CW',], aes(x = Date, y = Abundance)) + 
  geom_point() + 
  facet_wrap(vars(Gene), scales = 'free')

ggplot(nit.fix.long, aes(x = Date, y = Gene, fill = Abundance)) + geom_tile() + facet_wrap(vars(Site)) + theme_classic() +
  scale_y_discrete(labels = c('AntG', 'nifH', 'nifK', 'vnfD', 'vnfK', 'vnfH'))

ggplot(nit.fix.long[nit.fix.long$Gene == 'KO1383' | nit.fix.long$Gene == 'KO1382',], aes(x = Date, y = Abundance, col = Gene)) + 
  geom_point() + 
  facet_wrap(vars(Site), scales = 'free')

nit.sum <- nit.fix.long %>% filter(Gene == 'KO1383' | Gene == 'KO1382') %>% 
  group_by(Gene, Site, Date) %>% 
  summarise(mean = mean(Abundance), sd = sd(Abundance), n=n(), se = sd/sqrt(n))

ggplot(nit.sum[nit.sum$Site == 'PW',], aes(x = Date, y = mean, col = Gene, group = Gene)) + 
  geom_point(aes(shape = Gene), size = 10) + 
  geom_line(size = 1.3) + 
  scale_color_manual(values = c("#A04000","#117A65"), labels = c('nifH', 'nifK')) + 
  labs(y = expression(paste(x10^{8}, ' transcripts/g DW soil')), x = '') + 
  scale_y_continuous(limits = c(0, 2.3e8), breaks = c(0, 0.5e8, 1e8, 1.5e8, 2e8), labels = c('0.0', '0.5', '1.0', '1.5', '2.0')) +
  scale_x_discrete(drop = FALSE) + 
  scale_shape_manual(values = c(15,16), labels = c('nifH', 'nifK')) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean), width = 0.1, size = 1) +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=15,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

kruskal_test(KO1383~Date, data = nit.fix[nit.fix$Site == 'PW',])
dunn_test(KO1383~Date, data = nit.fix[nit.fix$Site == 'PW',], p.adjust.method = 'bonferroni')

kruskal_test(KO1383~Date, data = nit.fix[nit.fix$Site == 'CW',])
dunn_test(KO1383~Date, data = nit.fix[nit.fix$Site == 'CW',], p.adjust.method = 'bonferroni')

ggplot(nit.sum[nit.sum$Site == 'CW',], aes(x = Date, y = mean, col = Gene, group = Gene)) + 
  geom_point(aes(shape = Gene), size = 10) + 
  geom_line(size = 1.3) + 
  scale_color_manual(values = c("#A04000","#117A65"), labels = c('nifH', 'nifK')) + 
  labs(y = expression(paste(x10^{7}, ' transcripts/g DW soil')), x = '') + 
  scale_y_continuous(limits = c(0, 6.5e7), breaks = c(0, 2e7, 4e7, 6e7), labels = c('0.0', '2.0', '4.0', '6.0')) +
  scale_x_discrete(drop = FALSE) + 
  scale_shape_manual(values = c(15,16), labels = c('nifH', 'nifK')) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean), width = 0.1, size = 1) +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=15,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

kruskal.test(KO1383~Site, data = nit.fix)
kruskal.test(KO1382~Site, data = nit.fix)

# denitrification

denit <- kegg %>% select(KO1362, KO1363, KO1365, KO1378, KO1379, KO1361, KO1394, KO1385, KO1377, KO1366)

denit <- denit %>% rownames_to_column(var = 'ID')

denit$Sample <- str_sub(denit$ID, 1, 3)
denit$Site <- str_sub(denit$ID, 1, 2)
denit$Date <- str_sub(denit$ID, 5)
denit$Date <- factor(denit$Date, levels = c('18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'))

denit.long <- pivot_longer(denit, cols = !c('ID', 'Sample', 'Site', 'Date'), names_to = 'KEGG_Code', values_to = 'Abundance')

denit.fix.long$Gene <- 'NA'

denit.long[denit.long$KEGG_Code == 'KO1362',]$Gene <-'NarG'
denit.long[denit.long$KEGG_Code == 'KO1363',]$Gene <-'NarH'
denit.long[denit.long$KEGG_Code == 'KO1365',]$Gene <-'NarI'
denit.long[denit.long$KEGG_Code == 'KO1378',]$Gene <-'NapA'
denit.long[denit.long$KEGG_Code == 'KO1379',]$Gene <-'NapB'
denit.long[denit.long$KEGG_Code == 'KO1361',]$Gene <-'NirK'
denit.long[denit.long$KEGG_Code == 'KO1394',]$Gene <-'NirS'
denit.long[denit.long$KEGG_Code == 'KO1385',]$Gene <-'NorB'
denit.long[denit.long$KEGG_Code == 'KO1377',]$Gene <-'NorC'
denit.long[denit.long$KEGG_Code == 'KO1366',]$Gene <-'NosZ'

ggplot(denit.long, aes(x = Date, y = Gene, fill = Abundance)) + geom_tile() + facet_wrap(vars(Site)) + theme_classic()

ggplot(denit.long[denit.long$Site == 'CW',], aes(x = Date, y = Abundance)) + 
  geom_point() + 
  facet_wrap(vars(Gene), scales = 'free')

ggplot(denit.long[denit.long$Site == 'PW',], aes(x = Date, y = Abundance)) + 
  geom_point() + 
  facet_wrap(vars(Gene), scales = 'free')

levene_test(Abundance~Site, data= denit.long)
kruskal.test(Abundance~Site, denit.long)

denit.sum <- denit.long %>% 
  group_by(Gene, Site, Date) %>% 
  summarise(mean = mean(Abundance), sd = sd(Abundance), n=n(), se = sd/sqrt(n))

ggplot(denit.sum[denit.sum$Site == 'PW',], aes(x = Date, y = mean, col = Gene)) + 
  geom_point(size = 4) + 
  geom_line(aes(group = Gene), size = 1.3) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean), width = 0.1, size = 1) +
  labs(y = expression(paste(x10^{8}, ' transcripts/g DW soil')), x = '') + 
  scale_y_continuous(limits = c(0, 11.5e8), breaks = c(0, 2e8, 4e8, 6e8, 8e8, 10e8), labels = c('0.0', '2.0', '4.0', '6.0', '8.0', '10.0')) +
  scale_x_discrete(drop = FALSE) + 
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=15,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(), legend.position = c(0.9,0.8),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

ggplot(denit.sum[denit.sum$Site == 'CW',], aes(x = Date, y = mean, col = Gene)) + 
  geom_point(size = 4) + 
  geom_line(aes(group = Gene), size = 1.3) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean), width = 0.1, size = 1) +
  labs(y = expression(paste(x10^{8}, ' transcripts/g DW soil')), x = '') + 
  scale_y_continuous(limits = c(0, 1.5e8), breaks = c(0, 0.5e8, 1e8, 1.5e8), labels = c('0.0', '0.5', '1.0', '1.5')) +
  scale_x_discrete(drop = FALSE) + 
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=15,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(), legend.position = c(0.9,0.8),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

denit.dynamic.sum <- denit.sum %>% filter(Gene == 'NarG' | Gene == 'NarH' | Gene == 'NirK' | Gene == 'NorB')
denit.dynamic <- denit.fix.long %>% filter(Gene == 'NarG' | Gene == 'NarH' | Gene == 'NirK' | Gene == 'NorB')

ggplot(denit.dynamic.sum[denit.dynamic.sum$Site == 'PW',], aes(x = Date, y = mean, col = Gene)) + 
  geom_point(size = 4) + 
  geom_line(aes(group = Gene), size = 1.3) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean), width = 0.1, size = 1) +
  labs(y = expression(paste(x10^{8}, ' transcripts/g DW soil')), x = '') + 
  scale_y_continuous(limits = c(0, 11.5e8), breaks = c(0, 2e8, 4e8, 6e8, 8e8, 10e8), labels = c('0.0', '2.0', '4.0', '6.0', '8.0', '10.0')) +
  scale_x_discrete(drop = FALSE) + 
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=15,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(), legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

gene.test('CW', 'NirK', denit.fix.long)
kruskal.test(Abundance~Date, denit.fix.long[denit.fix.long$Site == 'CW' & denit.fix.long$Gene == 'NirK',])
dunn_test(Abundance~Date, data = denit.fix.long[denit.fix.long$Site == 'CW' & denit.fix.long$Gene == 'NirK',])

ggplot(denit.dynamic.sum[denit.dynamic.sum$Site == 'CW',], aes(x = Date, y = mean, col = Gene)) + 
  geom_point(size = 4) + 
  geom_line(aes(group = Gene), size = 1.3) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean), width = 0.1, size = 1) +
  labs(y = expression(paste(x10^{8}, ' transcripts/g DW soil')), x = '') + 
  scale_y_continuous(limits = c(0, 1.5e8), breaks = c(0, 0.5e8, 1e8, 1.5e8), labels = c('0.0', '0.5', '1.0', '1.5')) +
  scale_x_discrete(drop = FALSE) + 
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=15,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(), legend.position = c(0.9,0.8),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

# dnra

dnra <- kegg %>% select(KO1357, KO1358, KO1384, KO1395)

dnra <- dnra %>% rownames_to_column(var = 'ID')

dnra$Sample <- str_sub(dnra$ID, 1, 3)
dnra$Site <- str_sub(dnra$ID, 1, 2)
dnra$Date <- str_sub(dnra$ID, 5)
dnra$Date <- factor(dnra$Date, levels = c('18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'))

dnra.long <- pivot_longer(dnra, cols = !c('ID', 'Sample', 'Site', 'Date'), names_to = 'KEGG_Code', values_to = 'Abundance')

dnra.long$Gene <- 'NA'

dnra.long[dnra.long$KEGG_Code == 'KO1357',]$Gene <-'NirB'
dnra.long[dnra.long$KEGG_Code == 'KO1358',]$Gene <-'NirD'
dnra.long[dnra.long$KEGG_Code == 'KO1384',]$Gene <-'NrfA'
dnra.long[dnra.long$KEGG_Code == 'KO1395',]$Gene <-'NrfH'

ggplot(dnra.long, aes(x = Date, y = Gene, fill = Abundance)) + geom_tile() + facet_wrap(vars(Site)) + theme_classic()

ggplot(dnra.long[dnra.long$Site == 'CW',], aes(x = Date, y = Abundance)) + 
  geom_point() + 
  facet_wrap(vars(Gene), scales = 'free')

ggplot(dnra.long[dnra.long$Site == 'PW',], aes(x = Date, y = Abundance)) + 
  geom_point() + 
  facet_wrap(vars(Gene), scales = 'free')

levene_test(formula = Abundance~Site, data = dnra.long)
shapiro.test(dnra.long$Abundance)

kruskal.test(Abundance~Site, dnra.long)

dnra.sum <- dnra.long %>% 
  group_by(Gene, Site, Date) %>% 
  summarise(mean = mean(Abundance), sd = sd(Abundance), n=n(), se = sd/sqrt(n))

ggplot(dnra.sum[dnra.sum$Site == 'PW',], aes(x = Date, y = mean, col = Gene)) + 
  geom_point(size = 4) + 
  geom_line(aes(group = Gene), size = 1.3) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean), width = 0.1, size = 1) +
  labs(y = expression(paste(x10^{8}, ' transcripts/g DW soil')), x = '') + 
  scale_y_continuous(limits = c(0, 4.5e8), breaks = c(0, 1e8, 2e8, 3e8, 4e8), labels = c('0.0', '1.0', '2.0', '3.0', '4.0')) +
  scale_x_discrete(drop = FALSE) + 
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=15,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(), legend.position = c(0.9,0.8),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

ggplot(dnra.sum[dnra.sum$Site == 'CW',], aes(x = Date, y = mean, col = Gene)) + 
  geom_point(size = 4) + 
  geom_line(aes(group = Gene), size = 1.3) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean), width = 0.1, size = 1) +
  labs(y = expression(paste(x10^{7}, ' transcripts/g DW soil')), x = '') + 
  scale_y_continuous(limits = c(0, 6.7e7), breaks = c(0, 2e7, 4e7, 6e7), labels = c('0.0', '2.0', '4.0', '6.0')) +
  scale_x_discrete(drop = FALSE) + 
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=15,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(), legend.position = c(0.9,0.8),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

gene.test('CW', 'NirB', Data = dnra.long)
kruskal.test(Abundance~Date, dnra.long[dnra.long$Site == 'CW' & dnra.long$Gene == 'NrfA',])
dunn_test(Abundance~Date, data = dnra.long[dnra.long$Site == 'CW' & dnra.long$Gene == 'NrfA',])

# anra

anra <- kegg %>% select(KO1360, KO1387, KO1364, KO1356, KO1397, KO1359)

anra <- anra %>% rownames_to_column(var = 'ID')

anra$Sample <- str_sub(anra$ID, 1, 3)
anra$Site <- str_sub(anra$ID, 1, 2)
anra$Date <- str_sub(anra$ID, 5)
anra$Date <- factor(anra$Date, levels = c('18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'))

anra.long <- pivot_longer(anra, cols = !c('ID', 'Sample', 'Site', 'Date'), names_to = 'KEGG_Code', values_to = 'Abundance')

anra.long$Gene <- 'NA'

anra.long[anra.long$KEGG_Code == 'KO1360',]$Gene <-'NarB'
anra.long[anra.long$KEGG_Code == 'KO1387',]$Gene <-'NR'
anra.long[anra.long$KEGG_Code == 'KO1364',]$Gene <-'NasA'
anra.long[anra.long$KEGG_Code == 'KO1356',]$Gene <-'NasB'
anra.long[anra.long$KEGG_Code == 'KO1397',]$Gene <-'NIT-6'
anra.long[anra.long$KEGG_Code == 'KO1359',]$Gene <-'NirA'

ggplot(anra.long, aes(x = Date, y = Gene, fill = Abundance)) + geom_tile() + facet_wrap(vars(Site)) + theme_classic()

ggplot(anra.long[anra.long$Site == 'CW',], aes(x = Date, y = Abundance)) + 
  geom_point() + 
  facet_wrap(vars(Gene), scales = 'free')

ggplot(anra.long[anra.long$Site == 'PW',], aes(x = Date, y = Abundance)) + 
  geom_point() + 
  facet_wrap(vars(Gene), scales = 'free')

levene_test(formula = Abundance~Site, data = anra.long)
shapiro.test(anra.long$Abundance)

kruskal.test(Abundance~Site, anra.long)

anra.sum <- anra.long %>% 
  group_by(Gene, Site, Date) %>% 
  summarise(mean = mean(Abundance), sd = sd(Abundance), n=n(), se = sd/sqrt(n))

ggplot(anra.sum[anra.sum$Site == 'PW',], aes(x = Date, y = mean, col = Gene)) + 
  geom_point(size = 4) + 
  geom_line(aes(group = Gene), size = 1.3) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean), width = 0.1, size = 1) +
  labs(y = expression(paste(x10^{8}, ' transcripts/g DW soil')), x = '') + 
  scale_y_continuous(limits = c(0, 2.5e8), breaks = c(0, 0.5e8, 1.0e8, 1.5e8, 2.0e8), labels = c('0.0', '0.5', '1.0', '1.5', '2.0')) +
  scale_x_discrete(drop = FALSE) + 
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=15,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(), legend.position = c(0.9,0.8),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

ggplot(anra.sum[anra.sum$Site == 'CW',], aes(x = Date, y = mean, col = Gene)) + 
  geom_point(size = 4) + 
  geom_line(aes(group = Gene), size = 1.3) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean), width = 0.1, size = 1) +
  labs(y = expression(paste(x10^{7}, ' transcripts/g DW soil')), x = '') + 
  scale_y_continuous(limits = c(0, 9e7), breaks = c(0, 2e7, 4e7, 6e7, 8e7), labels = c('0.0', '2.0', '4.0', '6.0', '8.0')) +
  scale_x_discrete(drop = FALSE) + 
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=15,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(), legend.position = c(0.9,0.8),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

gene.test('CW', 'NR', Data = anra.long)

kruskal.test(Abundance~Date, anra.long[anra.long$Site == 'CW' & anra.long$Gene == 'NR',])
dunn_test(Abundance~Date, data = anra.long[anra.long$Site == 'CW' & anra.long$Gene == 'NR',])

# heat map of all interesting nitrogen-cycling genes

dynamic <- kegg %>% select(X, KO1381, KO1383, KO1382, KO1357, KO1364, KO1287, KO1288, KO1289)

dynamic <- rename(dynamic, ID = X)

dynamic$Sample <- str_sub(dynamic$ID, 1, 3)
dynamic$Site <- str_sub(dynamic$ID, 1, 2)
dynamic$Date <- str_sub(dynamic$ID, 5)
dynamic$Date <- factor(dynamic$Date, levels = c('18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'))

dynamic.long <- pivot_longer(dynamic, cols = !c('ID', 'Sample', 'Site', 'Date'), names_to = 'KEGG_Code', values_to = 'Abundance')

dynamic.long$Gene <- 'NA'

dynamic.long[dynamic.long$KEGG_Code == 'KO1381',]$Gene <-'nifD'
dynamic.long[dynamic.long$KEGG_Code == 'KO1382',]$Gene <-'nifH'
dynamic.long[dynamic.long$KEGG_Code == 'KO1383',]$Gene <-'nifK'
dynamic.long[dynamic.long$KEGG_Code == 'KO1357',]$Gene <-'nirB'
dynamic.long[dynamic.long$KEGG_Code == 'KO1364',]$Gene <-'nasA'
dynamic.long[dynamic.long$KEGG_Code == 'KO1287',]$Gene <-'amoA'
dynamic.long[dynamic.long$KEGG_Code == 'KO1288',]$Gene <-'amoB'
dynamic.long[dynamic.long$KEGG_Code == 'KO1289',]$Gene <-'amoC'

dynamic.long$Gene <- factor(dynamic.long$Gene, levels = c('nasA', 'nirB', 'nifK', 'nifH', 'nifD'))

ggplot(dynamic.long, aes(x = Date, y = Gene, fill = Abundance)) + geom_tile() + facet_wrap(vars(Site)) + theme_classic()

dynamic.m <- data.matrix(dynamic)

dynamic.t <- t(dynamic.m)

heatmap(dynamic.t[,str_detect(colnames(dynamic.t), 'PW') == TRUE])

dynamic.mean <- dynamic.long %>% group_by(Site, Date, Gene) %>% summarise(mean = mean(Abundance))

dynamic.mean.wide <- pivot_wider(dynamic.mean, id_cols = c(Site, Date), names_from = Gene, values_from = mean)

pw.mean.wide <- dynamic.mean.wide %>% filter(Site == 'PW') %>% column_to_rownames(var = 'Date') %>% select(!Site)

pw.mean.m <- t(data.matrix(pw.mean.wide))

pw.mean.ordered <- pw.mean.m[c(1,2,3,5,6,7,4,8),]

heatmap(pw.mean.ordered, Colv = NA, Rowv = NA, scale = 'row')

cw.mean.wide <- dynamic.mean.wide %>% filter(Site == 'CW') %>% column_to_rownames(var = 'Date') %>% select(!Site)

cw.mean.m <- t(data.matrix(cw.mean.wide))

cw.mean.ordered <- cw.mean.m[c(1,2,3,5,6,7,4,8),]

heatmap(cw.mean.ordered, Colv = NA, Rowv = NA, scale = 'row')

ggplot(dynamic.mean[dynamic.mean$Site == 'PW',], aes(x = Date, y = Gene, size = mean)) + geom_point() + theme_classic()

# testing amo genes

levene_test(dynamic.long[dynamic.long$Gene == 'nifH',], Abundance~Site)
dynamic.long[dynamic.long$Gene == 'amoA',] %>% shapiro_test(Abundance)

kruskal_test(data = dynamic.long[dynamic.long$Gene == 'nasA',], formula = Abundance ~ Site)

levene_test(dynamic.long[dynamic.long$Site == 'CW' & dynamic.long$Gene == 'amoA',], Abundance~Date)
dynamic.long[dynamic.long$Site == 'CW' & dynamic.long$Gene == 'amoA',] %>% shapiro_test(Abundance)

summary(aov(Abundance~Date, dynamic.long[dynamic.long$Site == 'CW' & dynamic.long$Gene == 'amoA',]))
TukeyHSD(aov(Abundance~Date, dynamic.long[dynamic.long$Site == 'CW' & dynamic.long$Gene == 'amoA',]))

kruskal_test(formula = Abundance~Date, data = dynamic.long[dynamic.long$Site == 'CW' & dynamic.long$Gene == 'amoA',])
dunn_test(formula = Abundance~Date, data = dynamic.long[dynamic.long$Site == 'CW' & dynamic.long$Gene == 'amoA',])

