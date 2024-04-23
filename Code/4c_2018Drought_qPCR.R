
################################################################################
#                                                                              #
# qPCR Results for 2018 Drought Samples                                        #
#                                                                              #
# Author: Anna Burns                                                           #
# Last edited: 19.04.2023                                                      #
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

################################################################################
# Data upload                                                                  #
################################################################################

full_data <- read.csv('../wetscapes/AmoA_qPCR/long_qpcr_wt.csv')
full_data$Date2 <- ymd(full_data$Date2)
long_qpcr_clean <- read_csv("~/wetscapes/AmoA_qPCR/long_qpcr_clean.csv")

full_data2 <- full_data %>% filter(Domain != 'Water')

################################################################################
# Figures and Statistics                                                       #
################################################################################

## PW

ggplot(full_data2[full_data2$Site == 'PW',], aes(x = Date, y = mean, col = Domain, group = Domain)) +
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

kruskal.test(Abundance~Date2, data = long_qpcr_clean[long_qpcr_clean$Domain == 'AOB' & long_qpcr_clean$Site == 'PW',])
kruskal_effsize(Abundance~Date2, data = long_qpcr_clean[long_qpcr_clean$Domain == 'AOB' & long_qpcr_clean$Site == 'PW',])
dunn_test(Abundance~Date2, data = long_qpcr_clean[long_qpcr_clean$Domain == 'AOB' & long_qpcr_clean$Site == 'PW',], p.adjust.method = 'bonferroni')

## PD

ggplot(full_data2[full_data2$Site == 'PD',], aes(x = Date, y = mean, col = Domain, group = Domain)) +
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

kruskal.test(Abundance~Date2, data = long_qpcr_clean[long_qpcr_clean$Domain == 'AOB' & long_qpcr_clean$Site == 'PD',])
kruskal_effsize(Abundance~Date2, data = long_qpcr_clean[long_qpcr_clean$Domain == 'AOB' & long_qpcr_clean$Site == 'PD',])
dunn_test(Abundance~Date2, data = long_qpcr_clean[long_qpcr_clean$Domain == 'AOB' & long_qpcr_clean$Site == 'PD',], p.adjust.method = 'bonferroni')

## CW

ggplot(full_data2[full_data2$Site == 'CW',], aes(x = Date, y = mean, col = Domain, group = Domain)) +
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
kruskal.test(Abundance~Date2, data = long_qpcr_clean[long_qpcr_clean$Domain == 'AOB' & long_qpcr_clean$Site == 'CW',])
kruskal_effsize(Abundance~Date2, data = long_qpcr_clean[long_qpcr_clean$Domain == 'AOB' & long_qpcr_clean$Site == 'CW',])
dunn_test(Abundance~Date2, data = long_qpcr_clean[long_qpcr_clean$Domain == 'AOB' & long_qpcr_clean$Site == 'CW',], p.adjust.method = 'bonferroni')

## CD

ggplot(full_data2[full_data2$Site == 'CD',], aes(x = Date, y = mean, col = Domain, group = Domain)) +
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

kruskal.test(Abundance~Date2, data = long_qpcr_clean[long_qpcr_clean$Domain == 'AOB' & long_qpcr_clean$Site == 'CD',])
kruskal_effsize(Abundance~Date2, data = long_qpcr_clean[long_qpcr_clean$Domain == 'AOB' & long_qpcr_clean$Site == 'CD',])
dunn_test(Abundance~Date2, data = long_qpcr_clean[long_qpcr_clean$Domain == 'AOB' & long_qpcr_clean$Site == 'CD',], p.adjust.method = 'bonferroni')
