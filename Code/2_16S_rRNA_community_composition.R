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
library(cowplot)

################################################################################
# Data upload                                                                  #
################################################################################

aoa.16s <- read.csv('Data/16srRNA.aoa.absolute.csv')
aob.16s <- read.csv('Data/16srRNA.aob.csv')

aoa.ssu <- read.csv('Data/aoa.ssu.csv') %>% select(-X)
aob.ssu <- read.csv('Data/aob.ssu.csv') %>% select(-X)
################################################################################
# Functions                                                                    #
################################################################################



################################################################################
# Code                                                                         #
################################################################################

# 16S AOA 

aoa.16s.sum <- aoa.16s %>% group_by(loc, season2, Order, group) %>% summarise(tot_abundance = sum(Abundance))

aoa.16s.sum2 <- aoa.16s.sum %>% group_by(loc, season2, Order) %>% summarise(mean = mean(tot_abundance), sd = sd(tot_abundance), n = n(), se = sd/sqrt(n)) %>%
  arrange(desc((Order))) %>% mutate(cumulative = cumsum(mean))

aoa.16s.sum2$season2 <- factor(aoa.16s.sum2$season2, levels = c('18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'),
                              labels = c('2018-04', '2018-06', '2018-08', '2018-10', '2018-12', '2019-02'))

aoa.16s.pw <- ggplot(aoa.16s.sum2[aoa.16s.sum2$loc == 'PW',], aes(x= season2, y = mean, group = Order)) + 
  geom_col(aes(fill = Order), width = .98) + 
  scale_fill_manual(values = c("#A04000","#117A65","#7D3C98")) + 
  scale_color_manual(values = c("#A04000","#117A65","#7D3C98")) +
  labs(y = expression(paste(x10^{7}, ' 16S rRNA copies/g DW soil')), x = '') + 
  scale_y_continuous(limits = c(0, 3.5e7), breaks = c(0, 1e7, 2e7, 3e7), labels = c('0.0', '1.0', '2.0', '3.0'),
                     expand = c(0,0)) +
  scale_x_discrete(drop = FALSE) + 
  guides(col = 'none', fill = 'none') + 
  geom_errorbar(data = aoa.16s.sum2[aoa.16s.sum2$loc == 'PW' & aoa.16s.sum2$Order == 'Nitrososphaerales',],
                aes(ymax = cumulative + se, ymin = cumulative-100000), width = 0.1, linewidth = 1, col = "#7D3C98") +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text.x=element_text(size=13,face = "bold"),
        axis.text.y=element_text(size=20, face='bold'),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

aoa.16s.cw <- ggplot(aoa.16s.sum2[aoa.16s.sum2$loc == 'CW',], aes(x= season2, y = mean, group = Order)) + 
  geom_col(aes(fill = Order), width = .98) + 
  scale_fill_manual(values = c("#A04000","#117A65","#7D3C98")) + 
  scale_color_manual(values = c("#A04000","#117A65","#7D3C98")) +
  labs(y = expression(paste(x10^{6}, ' 16S rRNA copies/g DW soil')), x = '') + 
  scale_y_continuous(limits = c(0, 5e6), breaks = c(0, 1e6, 2e6, 3e6, 4e6, 5e6), labels = c('0.0', '1.0', '2.0', '3.0', '4.0', '5.0'),
                     expand = c(0,0)) +
  scale_x_discrete(drop = FALSE) + 
  guides(col = 'none', fill = 'none') + 
  geom_errorbar(data = aoa.16s.sum2[aoa.16s.sum2$loc == 'CW' & aoa.16s.sum2$Order == 'Ca. Nitrosotaleales',],
                aes(ymax = cumulative + se, ymin = cumulative-10000), width = 0.1, linewidth = 1, col = "#117A65") +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text.x=element_text(size=13,face = "bold"),
        axis.text.y=element_text(size=20, face='bold'),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

## 16S AOB

aob.16s.sum <- aob.16s %>% group_by(loc, season2, Order, group) %>% summarise(tot_abundance = sum(Abundance))

aob.16s.sum2 <- aob.16s.sum %>% group_by(loc, season2, Order) %>% summarise(mean = mean(tot_abundance), sd = sd(tot_abundance), n = n(), se = sd/sqrt(n)) %>%
  arrange(desc((Order))) %>% mutate(cumulative = cumsum(mean))

aob.16s.sum2$season2 <- factor(aob.16s.sum2$season2, levels = c('18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'),
                               labels = c('2018-04', '2018-06', '2018-08', '2018-10', '2018-12', '2019-02'))

aob.16s.pw <- ggplot(aob.16s.sum2[aob.16s.sum2$loc == 'PW',], aes(x= season2, y = mean, group = Order)) + 
  geom_col(aes(fill = Order), width = .98) + 
  scale_fill_manual(values = c("#F0C808",'#499DD4')) + 
  scale_color_manual(values = c("#F0C808",'#499DD4')) +
  labs(y = expression(paste(x10^{9}, ' 16S rRNA copies/g DW soil')), x = '') + 
  scale_y_continuous(limits = c(0, 3.3e9), breaks = c(0, 1e9, 2e9, 3e9), labels = c('0.0', '1.0', '2.0', '3.0'),
                     expand = c(0,0)) +
  guides(fill = 'none', col = 'none') + 
  scale_x_discrete(drop = FALSE) + 
  geom_errorbar(aes(ymax = cumulative + se, ymin = cumulative-10000000, col = Order), width = 0.1, linewidth = 1) +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text.x=element_text(size=13,face = "bold"),
        axis.text.y=element_text(size=20, face='bold'),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

aob.16s.cw <- ggplot(aob.16s.sum2[aob.16s.sum2$loc == 'CW',], aes(x= season2, y = mean, group = Order)) + 
  geom_col(aes(fill = Order), width = .98) + 
  scale_fill_manual(values = c("#F0C808",'#499DD4')) + 
  scale_color_manual(values = c("#F0C808",'#499DD4')) +
  guides(fill = 'none', col = 'none') + 
  labs(y = expression(paste(x10^{8}, ' 16S rRNA copies/g DW soil')), x = '') + 
  scale_y_continuous(limits = c(0, 1.2e9), breaks = c(0, 2e8, 4e8, 6e8, 8e8, 1e9), labels = c('0.0', '2.0', '4.0', '6.0', '8.0', '10.0'),
                     expand = c(0,0)) +
  scale_x_discrete(drop = FALSE) + 
  geom_errorbar(aes(ymax = cumulative + se, ymin = cumulative-5000000, col = Order), width = 0.1, linewidth = 1) +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text.x=element_text(size=13,face = "bold"),
        axis.text.y=element_text(size=20, face='bold'),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

## SSU AOA PW

aoa.ssu.sum <- aoa.ssu %>% group_by(Site, Date, Order, Subsample) %>% summarise(tot_abundance = sum(Abundance))

aoa.ssu.sum2 <- aoa.ssu.sum %>% group_by(Site, Date, Order) %>% summarise(mean = mean(tot_abundance), sd = sd(tot_abundance), n = n(), se = sd/sqrt(n)) %>%
  arrange(Order) %>% mutate(cumulative = cumsum(mean))

aoa.ssu.sum2$Date <- factor(aoa.ssu.sum2$Date, levels = c('18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'),
                               labels = c('2018-04', '2018-06', '2018-08', '2018-10', '2018-12', '2019-02'))

aoa.ssu.pw <- ggplot(aoa.ssu.sum2[aoa.ssu.sum2$Site == 'PW',], aes(x= Date, y = mean, group = Order)) + 
  geom_col(aes(fill = Order), width = .98) + 
  scale_fill_manual(values = c("#7D3C98", "#A04000", "#117A65")) + 
  scale_color_manual(values = c("#7D3C98", "#A04000", "#117A65")) +
  labs(y = expression(paste(x10^{10}, ' SSU transcripts/g DW soil')), x = '') + 
  scale_y_continuous(limits = c(0, 3.5e10), breaks = c(0, 1e10, 2e10, 3e10), labels = c('0.0', '1.0', '2.0', '3.0'),
                     expand = c(0,0)) +
  scale_x_discrete(drop = FALSE) + 
  guides(col = 'none', fill = 'none') + 
  geom_errorbar(data = aoa.ssu.sum2[aoa.ssu.sum2$Site == 'PW' & aoa.ssu.sum2$Order == 'Nitrososphaerales' | aoa.ssu.sum2$Site == 'PW' & aoa.ssu.sum2$Order == 'Ca. Nitrosopumilales' & aoa.ssu.sum2$Date == '2018-10',],
                aes(ymax = cumulative + se, ymin = cumulative-100000000, col = Order), width = 0.1, linewidth = 1) +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text.x=element_text(size=13,face = "bold"),
        axis.text.y=element_text(size=20, face='bold'),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

aoa.ssu.sum3 <- aoa.ssu.sum2
aoa.ssu.sum3$Order <- factor(aoa.ssu.sum3$Order, levels = c('Nitrososphaerales', 'Ca. Nitrosotaleales', 'Ca. Nitrosopumilales'))

aoa.ssu.cw <- ggplot(aoa.ssu.sum3[aoa.ssu.sum3$Site == 'CW',], aes(x= Date, y = mean, group = Order)) + 
  geom_col(aes(fill = Order), width = .98) + 
  scale_fill_manual(values = c("#7D3C98", "#117A65", "#A04000")) + 
  scale_color_manual(values = c("#7D3C98", "#117A65", "#A04000")) +
  labs(y = expression(paste(x10^{8}, ' SSU transcripts/g DW soil')), x = '') + 
  scale_y_continuous(limits = c(0, 6.2e8), breaks = c(0, 2e8, 4e8, 6e8), labels = c('0.0', '2.0', '4.0', '6.0'),
                     expand = c(0,0)) +
  scale_x_discrete(drop = FALSE) + 
  guides(col = 'none', fill = 'none') + 
  geom_errorbar(data = aoa.ssu.sum3[aoa.ssu.sum3$Site == 'CW' & aoa.ssu.sum3$Order != 'Ca. Nitrosopumilales' | 
                                      aoa.ssu.sum3$Site == 'CW' & aoa.ssu.sum3$Order == 'Ca. Nitrosopumilales' & aoa.ssu.sum3$Date == '2018-08' | 
                                      aoa.ssu.sum3$Site == 'CW' & aoa.ssu.sum3$Order == 'Ca. Nitrosopumilales' & aoa.ssu.sum3$Date == '2018-10' ,],
                  aes(ymax = cumulative + se, ymin = cumulative-2000000, col = Order), width = 0.1, linewidth = 1) +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text.x=element_text(size=13,face = "bold"),
        axis.text.y=element_text(size=20, face='bold'),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)


## SSU AOB PW

aob.ssu.sum <- aob.ssu %>% group_by(Site, Date, Order, Subsample) %>% summarise(tot_abundance = sum(Abundance))

aob.ssu.sum2 <- aob.ssu.sum %>% group_by(Site, Date, Order) %>% summarise(mean = mean(tot_abundance), sd = sd(tot_abundance), n = n(), se = sd/sqrt(n))

aob.ssu.sum2$Date <- factor(aob.ssu.sum2$Date, levels = c('18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'),
                            labels = c('2018-04', '2018-06', '2018-08', '2018-10', '2018-12', '2019-02'))

aob.ssu.sum2$Order <- factor(aob.ssu.sum2$Order, levels = c('Nitrosomonadales', 'Nitrospirales'))

aob.ssu.sum2 <- aob.ssu.sum2 %>% arrange(desc(Order)) %>% mutate(cumulative = cumsum(mean))

aob.ssu.pw <- ggplot(aob.ssu.sum2[aob.ssu.sum2$Site == 'PW',], aes(x= Date, y = mean, group = Order)) + 
  geom_col(aes(fill = Order), width = .98) + 
  scale_fill_manual(values = c("#F0C808",'#499DD4')) + 
  scale_color_manual(values = c("#F0C808",'#499DD4')) +
  labs(y = expression(paste(x10^{11}, ' SSU transcripts/g DW soil')), x = '') + 
  scale_y_continuous(limits = c(0, 8.1e11), breaks = c(0, 2e11, 4e11, 6e11, 8e11), labels = c('0.0', '2.0', '4.0', '6.0', '8.0'),
                     expand = c(0,0)) +
  guides(fill = 'none', col = 'none') + 
  scale_x_discrete(drop = FALSE) + 
  geom_errorbar(aes(ymax = cumulative + se, ymin = cumulative-5000000000, col = Order), width = 0.1, linewidth = 1) +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text.x=element_text(size=13,face = "bold"),
        axis.text.y=element_text(size=20, face='bold'),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

aob.ssu.cw <- ggplot(aob.ssu.sum2[aob.ssu.sum2$Site == 'CW',], aes(x= Date, y = mean, group = Order)) + 
  geom_col(aes(fill = Order), width = .98) + 
  scale_fill_manual(values = c("#F0C808",'#499DD4')) + 
  scale_color_manual(values = c("#F0C808",'#499DD4')) +
  guides(fill = 'none', col = 'none') + 
  labs(y = expression(paste(x10^{11}, ' SSU transcripts/g DW soil')), x = '') + 
  scale_y_continuous(limits = c(0, 4.5e11), breaks = c(0, 1e11, 2e11, 3e11, 4e11), labels = c('0.0', '1.0', '2.0', '3.0', '4.0'),
                     expand = c(0,0)) +
  scale_x_discrete(drop = FALSE) + 
  geom_errorbar(aes(ymax = cumulative + se, ymin = cumulative-1000000000, col = Order), width = 0.1, linewidth = 1) +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text.x=element_text(size=13,face = "bold"),
        axis.text.y=element_text(size=20, face='bold'),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

## arrange figures

pw.16s <- plot_grid(aoa.16s.pw, aob.16s.pw)
cw.16s <- plot_grid(aoa.16s.cw, aob.16s.cw)

pw.ssu <- plot_grid(aoa.ssu.pw, aob.ssu.pw)
cw.ssu <- plot_grid(aoa.ssu.cw, aob.ssu.cw)
