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

sample.2018drought <- sample %>% dplyr::filter(season2 == '18-Feb' | season2 == '18-Apr'| season2 == '18-Jun' | 
                                          season2 == '18-Aug' | season2 == '18-Oct' | season2 == '18-Dec' | season2 == '19-Feb')

sample.2018drought.topsoil <- sample.2018drought %>% filter(depth == '05-10 cm')

sample.2018drought.topsoil$season2 <- factor(sample.2018drought.topsoil$season2, levels = c('18-Feb', '18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'))

nitrate.ammonium <- sample.2018drought.topsoil %>% select(NO3,NH4,loc,season2) %>% pivot_longer(!c(loc,season2), names_to = 'Nutrient', values_to='Volume')

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
