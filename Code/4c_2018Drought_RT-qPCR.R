################################################################################
#                                                                              #
# RT-qPCR results for 2018 drought cycle:                                      #
# Bacterial and archaeal amoA                                                  #
#                                                                              #  
# Author: Anna Burns                                                           #
# Last edited: 19.04.2024                                                      #
# Last tested: 25.04.2024                                                      #
#                                                                              #
################################################################################

################################################################################
# Required packages                                                            #
################################################################################

library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpattern)
library(rstatix)
library(stringr)

################################################################################
# Data upload                                                                  #
################################################################################

rt.qpcr <- read.csv('Data/amoa.rtqpcr.csv')

rt.summary <- rt.qpcr %>% group_by(Site, Date, Domain) %>% summarise(Mean = mean(Abundance), SD = sd(Abundance), n = n(), SE = SD/sqrt(n))

rt.summary$Date <- factor(rt.summary$Date, levels = c('2018-04', '2018-06', '2018-08', '2018-10', '2018-12', '2019-02'))

# combined qpcr and rt qpcr data

qpcr <- read.csv('./Data/amoa.qpcr.csv')


qpcr <- qpcr %>% select(Sample, Site, Date, Domain, Abundance)

rt.qpcr2 <- rt.qpcr %>% select(Sample, Site, Date, Domain, Abundance)

rt.qpcr2$Sample <- str_sub(rt.qpcr2$Sample, start = 5)
qpcr$Sample <- str_sub(qpcr$Sample, start = 6)

rt.qpcr2$Type <- 'RT-qPCR'
qpcr$Type <- 'qPCR'

full <- rbind(rt.qpcr2, qpcr)

full.summary <- full %>% na.omit() %>% group_by(Site, Type, Domain, Date) %>% summarise(Mean = mean(Abundance), sd = sd(Abundance), n = n(), se = sd/sqrt(n))

################################################################################
# Figures and Statistics                                                       #
################################################################################

# PD
ggplot(rt.summary[rt.summary$Site == 'PD' ,], aes(x = Date, y = Mean, col = Domain, group = Domain)) +
  geom_point(aes(shape = Domain), size = 10) + 
  geom_line(linewidth = 1.3) + 
  theme_classic() + 
  geom_errorbar(aes(ymax = Mean + SE, ymin = Mean, color = Domain), width = 0.1, linewidth = 1) +
  scale_shape_manual(values = c(15,17)) + 
  scale_color_manual(values = c("#A04000","#117A65","#7D3C98")) + 
  scale_y_continuous(limits = c(0, 1.5e8), breaks = c(0, 0.5e8, 1e8, 1.5e8), labels = c(0.0,0.5, 1.0,1.5)) +
  scale_x_discrete(drop = FALSE) + 
  labs(y = expression(paste('amoA  ', x10^{8}, ' copies/g DW soil')), x = '') + 
  theme(axis.title=element_text(size=20,face = "bold"),
        axis.title.y=element_text(margin = margin(l = 60)),
        axis.text=element_text(size=20,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

# PW
ggplot(rt.summary[rt.summary$Site == 'PW',], aes(x = Date, y = Mean, col = Domain, group = Domain)) +
  geom_point(aes(shape = Domain), size = 10) + 
  geom_line(linewidth = 1.3) + 
  theme_classic() + 
  geom_errorbar(aes(ymax = Mean + SE, ymin = Mean, color = Domain), width = 0.1, linewidth = 1) +
  scale_shape_manual(values = c(15,17)) + 
  scale_color_manual(values = c("#A04000","#117A65","#7D3C98")) + 
  scale_y_continuous(limits = c(0, 4.3e8), breaks = c(0, 1e8, 2e8, 3e8, 4e8), labels = c(0.0,1.0,2.0,3.0,4.0)) +
  scale_x_discrete(drop = FALSE) + 
  labs(y = expression(paste('amoA  ', x10^{8}, ' copies/g DW soil')), x = '') + 
  theme(axis.title=element_text(size=20,face = "bold"),
        axis.title.y=element_text(margin = margin(l = 60)),
        axis.text=element_text(size=20,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

# CW
ggplot(rt.summary[rt.summary$Site == 'CW',], aes(x = Date, y = Mean, col = Domain, group = Domain)) +
  geom_point(aes(shape = Domain), size = 10) + 
  geom_line(linewidth = 1.3) + 
  theme_classic() + 
  geom_errorbar(aes(ymax = Mean + SE, ymin = Mean, color = Domain), width = 0.1, linewidth = 1) +
  scale_shape_manual(values = c(15,17)) + 
  scale_color_manual(values = c("#A04000","#117A65","#7D3C98")) + 
  scale_y_continuous(limits = c(0, 2.5e8), breaks = c(0, 0.5e8, 1e8, 1.5e8, 2e8, 2.5e8), labels = c(0.0,0.5,1.0,1.5,2.0, 2.5)) +
  scale_x_discrete(drop = FALSE) + 
  labs(y = expression(paste('amoA  ', x10^{8}, ' copies/g DW soil')), x = '') + 
  theme(axis.title=element_text(size=20,face = "bold"),
        axis.title.y=element_text(margin = margin(l = 60)),
        axis.text=element_text(size=20,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

# summary figure with qPCR

ggplot(full, aes(x = Site, y = Abundance, fill = Domain, pattern = Type)) +
  geom_boxplot_pattern(color = "black", 
                       pattern_fill = "black",
                       pattern_angle = 45,
                       pattern_density = 0.1,
                       pattern_spacing = 0.025,
                       pattern_key_scale_factor = 0.6) + 
  scale_pattern_manual(values = c(RNA = "stripe", DNA = "none")) + 
  labs(y = expression(paste(x10^{8}, ' copies/g DW soil')), x = '') + 
  scale_y_continuous(limits = c(0, 8e8), breaks = c(0, 2e8, 4e8, 6e8, 8e8), labels = c('0.0', '2.0', ' 4.0', '6.0', '8.0')) +
  theme(axis.title=element_text(size=20,face = "bold"),
        axis.title.y=element_text(margin = margin(l = 60)),
        axis.text=element_text(size=20,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

ggplot(full.summary, aes(x = Site, y = Mean)) + 
  geom_col_pattern(aes(fill = Domain, pattern = Type), position = 'dodge', 
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) + 
  scale_pattern_manual(values = c(RNA = "stripe", DNA = "none")) + 
  geom_errorbar(aes(ymax = Mean + se, ymin = Mean), width = 0.1, linewidth = 1, position = position_dodge(width = 0.9)) + 
  labs(y = expression(paste(x10^{8}, ' copies/g DW soil')), x = '') + 
  scale_y_continuous(limits = c(0, 2e8), breaks = c(0, 0.5e8, 1e8, 1.5e8, 2e8), labels = c('0.0', '0.5', ' 1.0', '1.5', '2.0')) +
  theme(axis.title=element_text(size=20,face = "bold"),
        axis.title.y=element_text(margin = margin(l = 60)),
        axis.text=element_text(size=20,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

kruskal.test(Abundance~Type, data = full[full$Domain == 'AOB' & full$Site == 'CD',])
dunn_test(Abundance~Type, data = full[full$Domain == 'AOA' & full$Site == 'PD',], p.adjust.method = 'bonferroni')  
