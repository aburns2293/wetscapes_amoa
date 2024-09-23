################################################################################
#                                                                              #
# 2018 drought cycle: Nutrient analysis                                        #
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
library(tibble)

################################################################################
# Data upload                                                                  #
################################################################################

sample <- read.csv('Data/prok.sample.csv') %>% remove_rownames() %>% column_to_rownames(var = 'Sample')

sample.2018drought <- sample %>% dplyr::filter(season2 == '18-Apr'| season2 == '18-Jun' | 
                                          season2 == '18-Aug' | season2 == '18-Oct' | season2 == '18-Dec' | season2 == '19-Feb')

sample.2018drought.topsoil <- sample.2018drought %>% filter(depth == '05-10 cm')

sample.2018drought.topsoil$season2 <- factor(sample.2018drought.topsoil$season2, levels = c('18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'))

nitrate.ammonium <- sample.2018drought.topsoil %>% select(NO3,NH4,loc,season2) %>% pivot_longer(!c(loc,season2), names_to = 'Nutrient', values_to='Volume')

doc <- sample %>% select(loc, season2, TOM, Bio, hCOM, lCOM) %>% mutate(Other = TOM - Bio - hCOM - lCOM) %>%
  pivot_longer(!c(loc,season2), names_to = 'Fraction', values_to='Volume')

################################################################################
# Code                                                                         #
################################################################################

# NO3

ggplot(sample.2018drought.topsoil, aes(x = season2, y = NO3)) +
  geom_boxplot() +
  facet_wrap(vars(loc)) +
  xlab('') +
  ylab('NO3 (mgC/gDM)') +
  ggtitle('Topsoil nitrate fluxes during 2018 drought cycle') +
  theme_classic() +
  geom_rect(aes(xmin=2.5,
                xmax = 5.5,
                ymin = -Inf,
                ymax = Inf), fill = 'red', alpha = 0.003) +
  guides(col=guide_legend(title="Site"))

## relationship between NO3 and groundwater depth

ggplot(sample.2018drought.topsoil, aes(x = Mean_GW, y = NO3, group = loc, col = loc, shape = Drought_Status)) +
  geom_point() +
  xlab('Groundwater depth (cm)') +
  ylab('NO3 (mgC/gDM)') +
  ggtitle('Relationship between groundwater table depth and NO3 content in topsoil') +
  theme_classic() +
  guides(col=guide_legend(title="Site")) 

## relationship between NO3 and topsoil water content

ggplot(sample.2018drought.topsoil, aes(x = water, y = NO3, group = loc, col = loc, shape = Drought_Status)) +
  geom_point() +
  xlab('Groundwater depth (cm)') +
  ylab('NO3 (mgC/gDM)') +
  ggtitle('Relationship between water and NO3 content in topsoil') +
  theme_classic() +
  guides(col=guide_legend(title="Site")) 

# NH4

ggplot(sample.2018drought.topsoil, aes(x = season2, y = NH4)) +
  geom_boxplot() +
  facet_wrap(vars(loc)) +
  xlab('') +
  ylab('NH4 (mgC/gDM)') +
  ggtitle('Topsoil ammonium fluxes during 2018 drought cycle') +
  theme_classic() +
  geom_rect(aes(xmin=2.5,
                xmax = 5.5,
                ymin = -Inf,
                ymax = Inf), fill = 'red', alpha = 0.003) +
  guides(col=guide_legend(title="Site"))

## relationship between NH4 and groundwater depth

ggplot(sample.2018drought.topsoil, aes(x = Mean_GW, y = NH4, group = loc, col = loc, shape = Drought_Status)) +
  geom_point() +
  xlab('Groundwater depth (cm)') +
  ylab('NH4 (mgC/gDM)') +
  ggtitle('Relationship between groundwater table depth and NH4 content in topsoil') +
  theme_classic() +
  guides(col=guide_legend(title="Site")) 

## relationship between NH4 and topsoil water content

ggplot(sample.2018drought.topsoil, aes(x = water, y = NH4, group = loc, col = loc, shape = Drought_Status)) +
  geom_point() +
  xlab('Groundwater depth (cm)') +
  ylab('NH4 (mgC/gDM)') +
  ggtitle('Relationship between water and NH4 content in topsoil') +
  theme_classic() +
  guides(col=guide_legend(title="Site")) 

# Relationship between NO3/NH4 accumulation

ggplot(nitrate.ammonium, aes(x = season2, y = Volume, group_by = Nutrient)) +
  geom_boxplot(aes(col = Nutrient)) +
  facet_wrap(vars(loc)) +
  theme_classic() +
  xlab('') +
  ylab('Volume (mgC/gDM)') +
  ggtitle('Site-specific nitrate and ammonium accumulation during 2018 drought cycle') +
  geom_rect(aes(xmin=2.5,
                xmax = 5.5,
                ymin = -Inf,
                ymax = Inf), fill = 'red', alpha = 0.003)

# DOC figures and statistics

doc.summary <- doc %>% filter(Fraction == 'TOM') %>% group_by(loc, season2) %>% na.omit() %>% summarise(mean = mean(Volume), sd = sd(Volume), n = n(), se=sd/sqrt(n))

doc.summary$season2 <- factor(doc.summary$season2, levels = c('17-Apr', '17-Aug', '17-Nov', '18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb', '19-Jun', '19-Oct'), 
                              labels = c('2017-04', '2017-08', '2017-11', '2018-04', '2018-06', '2018-08', '2018-10', '2018-12', '2019-02', '2019-10'))

ggplot(doc.summary, aes(x = season2, y = mean, col = loc, group = loc)) +
  geom_point(aes(shape = loc), size = 10) +
  scale_x_discrete(drop = FALSE) + 
  geom_line(linewidth = 1.3) +
  scale_shape_manual(values = c(15,16, 17, 18))


levene_test(doc[doc$Fraction == 'TOM',], Volume~loc)
shapiro.test(doc[doc$Fraction == 'TOM',]$Volume)
kruskal.test(doc[doc$Fraction == 'lCOM',], Volume~loc)

ggplot(doc.summary[doc.summary$loc == 'PW',], aes(x = season2, y = mean, col = Fraction, group = Fraction)) + 
  geom_point(aes(shape = Fraction), size = 10) + 
  geom_line(linewidth = 1.3) + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean, color = Fraction), width = 0.1, linewidth = 1) +
  scale_shape_manual(values = c(15,16, 17, 18, 8)) + 
  scale_color_manual(values = c("#8DDBE0","#117A65","#7D3C98", '#E6C229', '#F46677')) + 
  scale_x_discrete(drop = FALSE) + 
  labs(y = 'mg C/g DW soil', x = '') + 
  theme(axis.title=element_text(size=20,face = "bold", colour = '#36454F'),
        axis.text=element_text(size=20,face = "bold"),
        axis.text.x = element_text(angle=45, vjust = 0.6),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.85,0.8),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

levene_test(doc[doc$loc == 'PW' & doc$Fraction == 'lCOM',], Volume~season2)
shapiro_test(doc[doc$loc == 'PW' & doc$Fraction == 'lCOM',]$Volume)

summary(aov(Volume ~ season2, doc[doc$loc == 'PW' & doc$Fraction == 'hCOM',]))
TukeyHSD((aov(Volume ~ season2, doc[doc$loc == 'PW' & doc$Fraction == 'hCOM',])))

kruskal_test(doc[doc$loc == 'PW' & doc$Fraction == 'lCOM',], Volume~season2)
dunn_test(doc[doc$loc == 'PW' & doc$Fraction == 'Bio',], Volume~season2, p.adjust.method= 'BH')

ggplot(doc.summary[doc.summary$loc == 'CW',], aes(x = season2, y = mean, col = Fraction, group = Fraction)) + 
  geom_point(aes(shape = Fraction), size = 10) + 
  geom_line(linewidth = 1.3) + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean, color = Fraction), width = 0.1, linewidth = 1) +
  scale_shape_manual(values = c(15,16, 17, 18, 8)) + 
  scale_color_manual(values = c("#8DDBE0","#117A65","#7D3C98", '#E6C229', '#F46677')) + 
  scale_x_discrete(drop = FALSE) + 
  labs(y = 'mg C/g DW soil', x = '') + 
  theme(axis.title=element_text(size=20,face = "bold", colour = '#36454F'),
        axis.text=element_text(size=20,face = "bold"),
        axis.text.x = element_text(angle=45, vjust = 0.6),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.87,0.8),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

levene_test(doc[doc$loc == 'CW' & doc$Fraction == 'TOM',], Volume~season2)
shapiro_test(doc[doc$loc == 'CW' & doc$Fraction == 'TOM',]$Volume)

summary(aov(Volume ~ season2, doc[doc$loc == 'CW' & doc$Fraction == 'TOM',]))
TukeyHSD((aov(Volume ~ season2, doc[doc$loc == 'CW' & doc$Fraction == 'TOM',])))

kruskal_test(doc[doc$loc == 'CW' & doc$Fraction == 'lCOM',], Volume~season2)
dunn_test(doc[doc$loc == 'CW' & doc$Fraction == 'lCOM',], Volume~season2, p.adjust.method= 'BH')

## nitrate and ammonium figures

nit.summary <- nitrate.ammonium %>% group_by(loc, season2, Nutrient) %>% summarize(mean = mean(Volume), n = n(), sd = sd(Volume), se = sd/sqrt(n))

nit.summary$season2 <- factor(nit.summary$season2, levels = c('18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'), 
                              labels = c('2018-04', '2018-06', '2018-08', '2018-10', '2018-12', '2019-02'))

ggplot(nit.summary[nit.summary$loc == 'CW',], aes(x = season2, y = mean, col = Nutrient, group = Nutrient)) + 
  geom_point(aes(shape = Nutrient), size = 10) + 
  geom_line(linewidth = 1.3) + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean, color = Nutrient), width = 0.1, linewidth = 1) +
  scale_shape_manual(values = c(15,16)) + 
  scale_color_manual(values = c("#7D3C98", '#E6C229')) + 
  scale_x_discrete(drop = FALSE) + 
  labs(y = 'mg C/g DW soil', x = '') + 
  theme(axis.title=element_text(size=20,face = "bold", colour = '#36454F'),
        axis.text=element_text(size=20,face = "bold"),
        axis.text.x = element_text(angle=45, vjust = 0.6),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

ggplot(nit.summary[nit.summary$loc == 'PW',], aes(x = season2, y = mean, col = Nutrient, group = Nutrient)) + 
  geom_point(aes(shape = Nutrient), size = 10) + 
  geom_line(linewidth = 1.3) + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean, color = Nutrient), width = 0.1, linewidth = 1) +
  scale_shape_manual(values = c(15,16)) + 
  scale_color_manual(values = c("#7D3C98", '#E6C229')) + 
  scale_x_discrete(drop = FALSE) + 
  labs(y = 'mg C/g DW soil', x = '') + 
  theme(axis.title=element_text(size=20,face = "bold", colour = '#36454F'),
        axis.text=element_text(size=20,face = "bold"),
        axis.text.x = element_text(angle=45, vjust = 0.6),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

