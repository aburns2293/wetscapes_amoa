
################################################################################
#                                                                              #
# qPCR Results for 2018 Drought Samples                                        #
#                                                                              #
# Author: Anna Burns                                                           #
# Last edited: 23.04.2023                                                      #
# Last tested: 25.04.2024                                                      #
#                                                                              #
################################################################################

################################################################################
# Required packages                                                            #
################################################################################

library(lubridate)
library(ggplot2)
library(dplyr)
library(stringr)
library(readxl)
library(rstatix)

################################################################################
# Data upload                                                                  #
################################################################################

qpcr <- read.csv('./Data/amoa.qpcr.csv')

# adding drought variable

qpcr$Drought_Status <- 'NA'
qpcr[qpcr$Date %in% c('2018-04', '2018-12', '2019-02'),]$Drought_Status <- 'Non-drought'
qpcr[qpcr$Date %in% c('2018-06', '2018-08', '2018-10'),]$Drought_Status <- 'Drought'

qpcr.summary <- qpcr %>% group_by(Site, Date, Domain) %>% summarise(mean = mean(Abundance),
                                                                    sd = sd(Abundance),
                                                                    n = n(),
                                                                    se = sd/sqrt(n))

################################################################################
# Figures and Statistics                                                       #
################################################################################

## PW

ggplot(qpcr.summary[qpcr.summary$Site == 'PW',], aes(x = Date, y = mean, col = Domain, group = Domain)) +
  geom_point(aes(shape = Domain), size = 10) + 
  geom_line(linewidth = 1.3) + 
  theme_classic() + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean, color = Domain), width = 0.1, linewidth = 1) +
  scale_shape_manual(values = c(15,17)) + 
  scale_color_manual(values = c("#A04000","#117A65","#7D3C98")) + 
  scale_y_continuous(limits = c(0, 3e8), breaks = c(0, 1e8, 2e8, 3e8), labels = c('0.0', '1.0', '2.0', '3.0')) +
  scale_x_discrete(drop = FALSE) + 
  labs(y = expression(paste(x10^{8}, ' copies/g DW soil')), x = '') + 
  theme(axis.title=element_text(size=20,face = "bold"),
        axis.title.y=element_text(margin = margin(l = 60)),
        axis.text=element_text(size=20,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

kruskal.test(Abundance~Date, data = qpcr[qpcr$Domain == 'AOB' & qpcr$Site == 'PW',])
kruskal_effsize(Abundance~Date, data = qpcr[qpcr$Domain == 'AOB' & qpcr$Site == 'PW',])
dunn_test(Abundance~Date, data = qpcr[qpcr$Domain == 'AOB' & qpcr$Site == 'PW',], p.adjust.method = 'bonferroni')

## PD

ggplot(qpcr.summary[qpcr.summary$Site == 'PD',], aes(x = Date, y = mean, col = Domain, group = Domain)) +
  geom_point(aes(shape = Domain), size = 10) + 
  geom_line(linewidth = 1.3) + 
  theme_classic() + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean, color = Domain), width = 0.1, linewidth = 1) +
  scale_shape_manual(values = c(15,17)) + 
  scale_color_manual(values = c("#A04000","#117A65","#7D3C98")) + 
  scale_y_continuous(limits = c(0, 6.6e8), breaks = c(0, 2e8, 4e8, 6e8), labels = c('0.0', '2.0', '4.0', '6.0')) +
  scale_x_discrete(drop = FALSE) + 
  labs(y = expression(paste(x10^{8}, ' copies/g DW soil')), x = '') + 
  theme(axis.title=element_text(size=20,face = "bold"),
        axis.title.y=element_text(margin = margin(l = 60)),
        axis.text=element_text(size=20,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

kruskal.test(Abundance~Date, data = qpcr[qpcr$Domain == 'AOB' & qpcr$Site == 'PD',])
kruskal_effsize(Abundance~Date, data = qpcr[qpcr$Domain == 'AOB' & qpcr$Site == 'PD',])
dunn_test(Abundance~Date, data = qpcr[qpcr$Domain == 'AOB' & qpcr$Site == 'PD',], p.adjust.method = 'bonferroni')

## CW

ggplot(qpcr.summary[qpcr.summary$Site == 'CW',], aes(x = Date, y = mean, col = Domain, group = Domain)) +
  geom_point(aes(shape = Domain), size = 10) + 
  geom_line(linewidth = 1.3) + 
  theme_classic() + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean, color = Domain), width = 0.1, linewidth = 1) +
  scale_shape_manual(values = c(15,17)) + 
  scale_color_manual(values = c("#A04000","#117A65","#7D3C98")) + 
  scale_y_continuous(limits = c(0, 1.2e8), breaks = c(0, 0.5e8, 1e8), labels = c('0.0', '0.5', '1.0')) +
  scale_x_discrete(drop = FALSE) + 
  labs(y = expression(paste(x10^{8}, ' copies/g DW soil')), x = '') + 
  theme(axis.title=element_text(size=20,face = "bold"),
        axis.title.y=element_text(margin = margin(l = 60)),
        axis.text=element_text(size=20,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

kruskal.test(Abundance~Date2, data = qpcr[qpcr$Domain == 'AOA' & qpcr$Site == 'CW',])
kruskal_effsize(Abundance~Date, data = qpcr[qpcr$Domain == 'AOB' & qpcr$Site == 'CW',])
dunn.test(qpcr[qpcr$Domain == 'AOA' & qpcr$Site == 'CW',]$Abundance, qpcr[qpcr$Domain == 'AOA' & qpcr$Site == 'CW',]$Drought_Status, method = 'bonferroni')

## CD

ggplot(qpcr.summary[qpcr.summary$Site == 'CD',], aes(x = Date, y = mean, col = Domain, group = Domain)) +
  geom_point(aes(shape = Domain), size = 10) + 
  geom_line(linewidth = 1.3) + 
  theme_classic() + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean, color = Domain), width = 0.1, linewidth = 1) +
  scale_shape_manual(values = c(15,17)) + 
  scale_color_manual(values = c("#A04000","#117A65","#7D3C98")) + 
  scale_y_continuous(limits = c(0, 1.7e8), breaks = c(0, 0.5e8, 1e8, 1.5e8), labels = c('0.0', '0.5', '1.0', '1.5')) +
  scale_x_discrete(drop = FALSE) + 
  labs(y = expression(paste(x10^{8}, ' copies/g DW soil')), x = '') + 
  theme(axis.title=element_text(size=20,face = "bold"),
        axis.title.y=element_text(margin = margin(l = 60)),
        axis.text=element_text(size=20,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

kruskal.test(Abundance~Date, data = qpcr[qpcr$Domain == 'AOB' & qpcr$Site == 'CD',])
kruskal_effsize(Abundance~Date, data = qpcr[qpcr$Domain == 'AOB' & qpcr$Site == 'CD',])
dunn_test(Abundance~Date, data = qpcr[qpcr$Domain == 'AOB' & qpcr$Site == 'CD',], p.adjust.method = 'bonferroni')
