################################################################################
#                                                                              #
# SSU- and 16S- based taxonomic composition of AOA and AOB                     #
#                                                                              #
# Author: Anna Burns                                                           #
# Last edited: 27.09.2024                                                      #
# Last tested: 27.09.2024                                                      #
#                                                                              #
################################################################################

################################################################################
# Required packages                                                            #
################################################################################

library(dplyr)
library(ggplot2)

################################################################################
# Data upload                                                                  #
################################################################################

aoa.16s <- read.csv('Data/16srRNA.aoa.absolute.csv')
ssu <- read.csv('Data/16srRNA.ssu.csv')
  
aob.16s <- read.csv('Data/16srRNA.aob.csv')
################################################################################
# Functions                                                                    #
################################################################################



################################################################################
# Code                                                                         #
################################################################################

# 16S AOA 

aoa.16s.sum <- aoa.16s %>% group_by(loc, season2, Order, group) %>% summarise(tot_abundance = sum(Abundance))

aoa.16s.sum2 <- aoa.16s.sum %>% group_by(loc, season2, Order) %>% summarise(mean = mean(tot_abundance), sd = sd(tot_abundance), n = n(), se = sd/sqrt(n))

aoa.16s.sum2$season2 <- factor(aoa.16s.sum2$season2, levels = c('18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'),
                              labels = c('2018-04', '2018-06', '2018-08', '2018-10', '2018-12', '2019-02'))

ggplot(aoa.16s.sum2[aoa.16s.sum2$loc == 'PW',], aes(x= season2, y = mean, group = Order)) + 
  geom_col(aes(fill = Order), width = .95) + 
  scale_fill_manual(values = c("#A04000","#117A65","#7D3C98")) + 
  scale_color_manual(values = c("#A04000","#117A65","#7D3C98")) +
  labs(y = expression(paste(x10^{7}, ' 16S rRNA copies/g DW soil')), x = '') + 
  scale_y_continuous(limits = c(0, 3.2e7), breaks = c(0, 1e7, 2e7, 3e7), labels = c('0.0', '1.0', '2.0', '3.0')) +
  scale_x_discrete(drop = FALSE) + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean, color = Order), width = 0.1, linewidth = 1) +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=15,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

ggplot(aoa.16s.sum2[aoa.16s.sum2$loc == 'CW',], aes(x= season2, y = mean, group = Order)) + 
  geom_col(aes(fill = Order), width = .95) + 
  scale_fill_manual(values = c("#A04000","#117A65","#7D3C98")) + 
  scale_color_manual(values = c("#A04000","#117A65","#7D3C98")) +
  labs(y = expression(paste(x10^{6}, ' 16S rRNA copies/g DW soil')), x = '') + 
  scale_y_continuous(limits = c(0, 5e6), breaks = c(0, 1e6, 2e6, 3e6, 4e6, 5e6), labels = c('0.0', '1.0', '2.0', '3.0', '4.0', '5.0')) +
  scale_x_discrete(drop = FALSE) + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean, color = Order), width = 0.1, linewidth = 1) +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=15,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

## 16S AOB

## SSU AOA PW

## SSU AOB PW

