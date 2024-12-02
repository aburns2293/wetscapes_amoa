
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

qpcr.summary <- qpcr %>% group_by(Site, Date, Domain, Drought_Status) %>% summarise(mean = mean(Abundance),
                                                                    sd = sd(Abundance),
                                                                    n = n(),
                                                                    se = sd/sqrt(n))

################################################################################
# Figures and Statistics                                                       #
################################################################################

## PW AOA

ggplot(qpcr.summary[qpcr.summary$Site == 'PW' & qpcr.summary$Domain == "AOA",], aes(x = Date, y = mean, group = Domain)) +
  geom_point( size = 10, aes(col = Drought_Status)) + 
  geom_line(linewidth = 1.3) + 
  theme_classic() + 
  scale_color_manual(values = c("#A04000","#117A65")) + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean+5000000, col = Drought_Status), width = 0.1, linewidth = 1) +
  scale_y_continuous(limits = c(0, 1.3e8), breaks = c(0, 3e7, 6e7, 9e7, 12e7), labels = c('0.0', '3.0', '6.0', '9.0', "12.0")) +
  scale_x_discrete(drop = FALSE) + 
  guides(color = "none") +
  labs(y = expression(paste(x10^{7}, ' copies/g DW soil')), x = '') + 
  theme(axis.title=element_text(size=20,face = "bold"),
        axis.text=element_text(size=20,face = "bold"),
        axis.text.x=element_text(angle=45, hjust=1),
        title = element_text(size = 20, face = 'bold'),
        axis.line = element_line(linewidth = 1.3),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

ggplot(qpcr[qpcr$Site == "PW" & qpcr$Domain == "AOA",], aes(x = Drought_Status, y = Abundance)) +
  geom_boxplot(aes(fill = Drought_Status), outlier.size = 4) +
  theme_classic() + 
  xlab("") +
  ylab("") +
  guides(fill = "none") +
  scale_fill_manual(values = c("#A04000","#117A65")) + 
  scale_y_continuous(limits = c(0, 1.3e8), breaks = c(0, 3e7, 6e7, 9e7, 12e7), labels = c('0.0', '3.0', '6.0', '9.0', "12.0")) +
  theme(axis.title=element_text(size=20,face = "bold"),
        axis.text.x=element_blank(),
        axis.line.y=element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        title = element_text(size = 20, face = 'bold'),
        axis.line.x = element_line(linewidth = 1.3, color = "black", linetype = "dashed"),
        aspect.ratio = 1)

kruskal.test(Abundance~Date, data = qpcr[qpcr$Domain == 'AOB' & qpcr$Site == 'PW',])
kruskal_effsize(Abundance~Date, data = qpcr[qpcr$Domain == 'AOB' & qpcr$Site == 'PW',])
dunn_test(Abundance~Date, data = qpcr[qpcr$Domain == 'AOB' & qpcr$Site == 'PW',], p.adjust.method = 'BH')

## PW AOB

ggplot(qpcr.summary[qpcr.summary$Site == 'PW' & qpcr.summary$Domain == "AOB",], aes(x = Date, y = mean, group = Domain)) +
  geom_point( size = 10, aes(col = Drought_Status)) + 
  geom_line(linewidth = 1.3) + 
  theme_classic() + 
  scale_color_manual(values = c("#A04000","#117A65")) + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean+.1e8, col = Drought_Status), width = 0.1, linewidth = 1) +
  scale_y_continuous(limits = c(0, 3.3e8), breaks = c(0, 1e8, 2e8, 3e8), labels = c('0.0', '1.0', '2.0', '3.0')) +
  scale_x_discrete(drop = FALSE) + 
  guides(color = "none") +
  labs(y = expression(paste(x10^{8}, ' copies/g DW soil')), x = '') + 
  theme(axis.title=element_text(size=20,face = "bold"),
        axis.text=element_text(size=20,face = "bold"),
        axis.text.x=element_text(angle=45, hjust=1),
        title = element_text(size = 20, face = 'bold'),
        axis.line = element_line(linewidth = 1.3),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

ggplot(qpcr[qpcr$Site == "PW" & qpcr$Domain == "AOB",], aes(x = Drought_Status, y = Abundance)) +
  geom_boxplot(aes(fill = Drought_Status), outlier.size = 4) +
  theme_classic() + 
  xlab("") +
  ylab("") +
  guides(fill = "none") +
  scale_fill_manual(values = c("#A04000","#117A65")) + 
  scale_y_continuous(limits = c(0, 3.3e8), breaks = c(0, 1e8, 2e8, 3e8), labels = c('0.0', '1.0', '2.0', '3.0')) +
  theme(axis.title=element_text(size=20,face = "bold"),
        axis.text.x=element_blank(),
        axis.line.y=element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        title = element_text(size = 20, face = 'bold'),
        axis.line.x = element_line(linewidth = 1.3, color = "black", linetype = "dashed"),
        aspect.ratio = 1)

## CW AOA

ggplot(qpcr.summary[qpcr.summary$Site == 'CW' & qpcr.summary$Domain == "AOA",], aes(x = Date, y = mean, group = Domain)) +
  geom_point( size = 10, aes(col = Drought_Status)) + 
  geom_line(linewidth = 1.3) + 
  theme_classic() + 
  scale_color_manual(values = c("#A04000","#117A65")) + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean+50000, col = Drought_Status), width = 0.1, linewidth = 1) +
  scale_y_continuous(limits = c(0, 6e6), breaks = c(0, 2e6, 4e6, 6e6), labels = c('0.0', '2.0', '4.0', '6.0')) +
  scale_x_discrete(drop = FALSE) + 
  guides(color = "none") +
  labs(y = expression(paste(x10^{6}, ' copies/g DW soil')), x = '') + 
  theme(axis.title=element_text(size=20,face = "bold"),
        axis.text=element_text(size=20,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        axis.line = element_line(linewidth = 1.3),
        axis.text.x = element_text(angle=45, hjust = 1),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

ggplot(qpcr[qpcr$Site == "CW" & qpcr$Domain == "AOA",], aes(x = Drought_Status, y = Abundance)) +
  geom_boxplot(aes(fill = Drought_Status), outlier.size = 4) +
  theme_classic() + 
  xlab("") +
  ylab("") +
  guides(fill = "none") +
  scale_fill_manual(values = c("#A04000","#117A65")) + 
  scale_y_continuous(limits = c(0, 6e6), breaks = c(0, 2e6, 4e6, 6e6), labels = c('0.0', '2.0', '4.0', '6.0')) +
  theme(axis.title=element_text(size=20,face = "bold"),
        axis.text.x=element_blank(),
        axis.line.y=element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        title = element_text(size = 20, face = 'bold'),
        axis.line.x = element_line(linewidth = 1.3, color = "black", linetype = "dashed"),
        aspect.ratio = 1)



kruskal.test(Abundance~Date2, data = qpcr[qpcr$Domain == 'AOA' & qpcr$Site == 'CW',])
kruskal_effsize(Abundance~Date, data = qpcr[qpcr$Domain == 'AOB' & qpcr$Site == 'CW',])
dunn_test(qpcr[qpcr$Domain == 'AOA' & qpcr$Site == 'CW',], Abundance~Date2, p.adjust.method = 'BH')

## CW AOB

ggplot(qpcr.summary[qpcr.summary$Site == 'CW' & qpcr.summary$Domain == "AOB",], aes(x = Date, y = mean, group = Domain)) +
  geom_point( size = 10, aes(col = Drought_Status)) + 
  geom_line(linewidth = 1.3) + 
  theme_classic() + 
  scale_color_manual(values = c("#A04000","#117A65")) + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean+5000000, col = Drought_Status), width = 0.1, linewidth = 1) +
  scale_y_continuous(limits = c(0, 13e7), breaks = c(0, 3e7, 6e7, 9e7, 12e7), labels = c('0.0', '3.0', '6.0', '9.0', "12.0")) +
  scale_x_discrete(drop = FALSE) + 
  guides(color = "none") +
  labs(y = expression(paste(x10^{6}, ' copies/g DW soil')), x = '') + 
  theme(axis.title=element_text(size=20,face = "bold"),
        axis.text=element_text(size=20,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        axis.line = element_line(linewidth = 1.3),
        axis.text.x = element_text(angle=45, hjust = 1),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

ggplot(qpcr[qpcr$Site == "CW" & qpcr$Domain == "AOB",], aes(x = Drought_Status, y = Abundance)) +
  geom_boxplot(aes(fill = Drought_Status), outlier.size = 4) +
  theme_classic() + 
  xlab("") +
  ylab("") +
  guides(fill = "none") +
  scale_fill_manual(values = c("#A04000","#117A65")) + 
  scale_y_continuous(limits = c(0, 13e7), breaks = c(0, 3e7, 6e7, 9e7, 12e7), labels = c('0.0', '3.0', '6.0', '9.0', "12.0")) +
  theme(axis.title=element_text(size=20,face = "bold"),
        axis.text.x=element_blank(),
        axis.line.y=element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        title = element_text(size = 20, face = 'bold'),
        axis.line.x = element_line(linewidth = 1.3, color = "black", linetype = "dashed"),
        aspect.ratio = 1)

