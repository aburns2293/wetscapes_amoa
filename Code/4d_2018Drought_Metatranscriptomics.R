################################################################################
#                                                                              #
# Analysis of metatranscriptomic dataset for PW + CW, 04.2018 - 02.2019        #
#                                                                              #
# Author: Anna Burns                                                           #
# Last edited: 15.01.2024                                                      #
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

################################################################################
# Code                                                                         #
################################################################################

# preparing aoa data

aoa$month <- factor(aoa$month, levels = c('18Apr', '18Jun', '18Aug', '18Oct', '18Dec', '19Feb'))

aoa_summary <- aoa %>% group_by(Site, month, Subunit) %>% summarise(mean = mean(Amount), sd = sd(Amount), n = n(), se = sd/sqrt(n))

aoa_summary$month2 <- str_c(str_sub(aoa_summary$month, 1, 2), str_sub(aoa_summary$month, 3,5), sep = '-')

aoa_summary$month2 <- factor(aoa_summary$month2, levels = c('18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'))

# preparing aob data

aob$month <- factor(aob$month, levels = c('18Apr', '18Jun', '18Aug', '18Oct', '18Dec', '19Feb'))

aob_summary <- aob %>% group_by(Site, month, Subunit) %>% summarise(mean = mean(Amount), sd = sd(Amount), n = n(), se = sd/sqrt(n))

aob_summary$month2 <- str_c(str_sub(aob_summary$month, 1, 2), str_sub(aob_summary$month, 3,5), sep = '-')

aob_summary$month2 <- factor(aob_summary$month2, levels = c('18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'))

################################################################################
# Figures and Statistics                                                       #
################################################################################

# aoa

## coastal wet site

ggplot(aoa_summary[aoa_summary$Site == 'CW',], aes(x = month2, y = mean, col = Subunit, group = Subunit)) +
  geom_point(aes(shape = Subunit), size = 10) + 
  geom_line(linewidth = 1.3) + 
  theme_classic() + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean, color = Subunit), width = 0.1, linewidth = 1) +
  scale_shape_manual(values = c(15,16,17)) + 
  scale_color_manual(values = c("#A04000","#117A65","#7D3C98")) + 
  scale_y_continuous(limits = c(0, 1e7), breaks = c(0, 3e6, 6e6, 9e6), labels = c(0.0,3.0, 6.0,9.0)) +
  scale_x_discrete(drop = FALSE) + 
  labs(y = expression(paste(x10^{6}, ' copies/g DW soil')), x = '') + 
  theme(axis.title=element_text(size=20,face = "bold"),
        axis.title.y=element_text(margin = margin(l = 60)),
        axis.text=element_text(size=20,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

## percolation wet site

ggplot(aoa_summary[aoa_summary$Site == 'PW',], aes(x = month2, y = mean, col = Subunit, group = Subunit)) +
  geom_point(aes(shape = Subunit), size = 10) + 
  geom_line(size = 1.3) + 
  theme_classic() + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean, color = Subunit), width = 0.1, linewidth = 1) +
  scale_shape_manual(values = c(15,16,17)) + 
  scale_color_manual(values = c("#A04000","#117A65","#7D3C98")) + 
  scale_y_continuous(limits = c(0, 1.7e8), breaks = c(0, 5e7, 1e8, 1.5e8), labels = c(0,0.5,1.0,1.5)) +
  scale_x_discrete(drop = FALSE) + 
  labs(x = '', y = expression(paste(x10^{8}, ' copies/g DW soil'))) +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=20,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        axis.title.y=element_text(margin = margin(l = 60)),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

# aob

# coastal wet site

ggplot(aob_summary[aob_summary$Site == 'CW',], aes(x = month2, y = mean, col = Subunit, group = Subunit)) +
  geom_point(aes(shape = Subunit), size = 10) + 
  geom_line(size = 1.3) + 
  theme_classic() + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean, color = Subunit), width = 0.1, linewidth = 1) +
  scale_shape_manual(values = c(15,16,17)) + 
  scale_color_manual(values = c("#117A65","#7D3C98")) + 
  scale_y_continuous(limits = c(0, 3e6), breaks = c(0, 1e6, 2e6, 3e6), labels = c(0,1,2,3)) +
  scale_x_discrete(drop = FALSE) + 
  labs(x = '', y = expression(paste(x10^{6}, ' copies/g DW soil'))) +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=20,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        axis.title.y=element_text(margin = margin(l = 60)),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

# percolation wet site

ggplot(aob_summary[aob_summary$Site == 'PW',], aes(x = month2, y = mean, col = Subunit, group = Subunit)) +
  geom_point(aes(shape = Subunit), size = 10) + 
  geom_line(size = 1.3) + 
  theme_classic() + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean, color = Subunit), width = 0.1, linewidth = 1) +
  scale_shape_manual(values = c(15,16,17)) + 
  scale_color_manual(values = c("#117A65","#7D3C98")) + 
  scale_y_continuous(limits = c(0,7e6), breaks = c(0, 3e6, 6e6), labels = c(0,3,6)) +
  scale_x_discrete(drop = FALSE) + 
  labs(x = '', y = expression(paste(x10^{6}, ' copies/g DW soil'))) +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=20,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        axis.title.y=element_text(margin = margin(l = 60)),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())
