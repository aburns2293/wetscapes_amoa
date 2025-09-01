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
library(rstatix)

################################################################################
# Data upload                                                                  #
################################################################################

aoa.16s <- read.csv('Data/16srRNA.aoa.absolute.csv')
aob.16s <- read.csv('Data/16srRNA.aob.csv') %>% filter(Family == 'Nitrosomonadaceae')

aoa.ssu <- read.csv('Data/aoa.ssu.csv') %>% select(-X)
aob.ssu <- read.csv('Data/aob.ssu.csv') %>% select(-X) %>% filter(Family == 'Nitrosomonadaceae')

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
  theme(axis.title=element_text(size=20,face = "bold"),axis.text.x=element_text(size=20,face = "bold", angle=45, hjust = 1),
        axis.text.y=element_text(size=20, face='bold'),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

aoa.16s.pw

## statistics

aoa.16s.sum3 <- aoa.16s.sum %>% group_by(loc, season2, group) %>% dplyr::summarise(sum = sum(tot_abundance))

leveneTest(sum~season2, aoa.16s.sum3[aoa.16s.sum3$loc == 'PW',])
shapiro.test(aoa.16s.sum3[aoa.16s.sum3$loc == 'PW',]$sum)
figs <- lm(sum~season2, aoa.16s.sum3[aoa.16s.sum3$loc == 'PW',])
plot(figs)

summary(aov(sum~season2, aoa.16s.sum3[aoa.16s.sum3$loc == 'PW',]))
TukeyHSD(aov(sum~season2, aoa.16s.sum3[aoa.16s.sum3$loc == 'PW',]))

## cw
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
  theme(axis.title=element_text(size=20,face = "bold"),axis.text.x=element_text(size=20,face = "bold", angle=45, hjust = 1),
        axis.text.y=element_text(size=20, face='bold'),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

## no statistics necessary for CW AOA, only one time point

## 16S AOB

aob.16s.sum <- aob.16s %>% group_by(loc, season2, group) %>% summarise(tot_abundance = sum(Abundance))

aob.16s.sum2 <- aob.16s.sum %>% group_by(loc, season2) %>% summarise(mean = mean(tot_abundance), sd = sd(tot_abundance), n = n(), se = sd/sqrt(n))

aob.16s.sum2$season2 <- factor(aob.16s.sum2$season2, levels = c('18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'),
                               labels = c('2018-04', '2018-06', '2018-08', '2018-10', '2018-12', '2019-02'))

aob.16s.pw <- ggplot(aob.16s.sum2[aob.16s.sum2$loc == 'PW',], aes(x= season2, y = mean)) + 
  geom_col(width = .95, col = '#499DD4', fill = '#499DD4') + 
  scale_fill_manual(values = c('#499DD4')) + 
  scale_color_manual(values = c('#499DD4')) +
  labs(y = expression(paste(x10^{9}, ' 16S rRNA copies/g DW soil')), x = '') + 
  scale_y_continuous(limits = c(0, 2.8e9), breaks = c(0, 1e9, 2e9), labels = c('0.0', '1.0', '2.0'),
                     expand = c(0,0)) +
  guides(fill = 'none', col = 'none') + 
  scale_x_discrete(drop = FALSE) + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean, col = '#499DD4'), width = 0.1, linewidth = 1) +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text.x=element_text(size=20,face = "bold", angle=45, hjust = 1),
        axis.text.y=element_text(size=20, face='bold'),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

aob.16s.pw

## statistics
aob.16s.sum <- aob.16s.sum %>% ungroup()
levene_test(aob.16s.sum[aob.16s.sum$loc == 'PW',], tot_abundance ~ season2)
shapiro.test(aob.16s.sum[aob.16s.sum$loc == 'PW' & aob.16s.sum$season2 == '18-Dec',]$tot_abundance)

summary(aov(tot_abundance~season2, aob.16s.sum[aob.16s.sum$loc == 'PW',]))
TukeyHSD(aov(tot_abundance~season2, aob.16s.sum[aob.16s.sum$loc == 'PW',]))

aob.16s.cw <- ggplot(aob.16s.sum2[aob.16s.sum2$loc == 'CW',], aes(x= season2, y = mean)) + 
  geom_col(fill = '#499DD4', col = '#499DD4', width = .95) + 
  guides(fill = 'none', col = 'none') + 
  labs(y = expression(paste(x10^{8}, ' 16S rRNA copies/g DW soil')), x = '') + 
  scale_y_continuous(limits = c(0, 8.2e8), breaks = c(0, 2e8, 4e8, 6e8, 8e8), labels = c('0.0', '2.0', '4.0', '6.0', '8.0'),
                     expand = c(0,0)) +
  scale_x_discrete(drop = FALSE) + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean), col = '#499DD4', width = 0.1, linewidth = 1) +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text.x=element_text(size=20,face = "bold", angle=45, hjust = 1),
        axis.text.y=element_text(size=20, face='bold'),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

aob.16s.cw

levene_test(aob.16s.sum[aob.16s.sum$loc == 'CW',], tot_abundance ~ season2)
shapiro.test(aob.16s.sum[aob.16s.sum$loc == 'CW' & aob.16s.sum$season2 == '19-Feb',]$tot_abundance)

summary(aov(tot_abundance~season2, aob.16s.sum[aob.16s.sum$loc == 'CW',]))
TukeyHSD(aov(tot_abundance~season2, aob.16s.sum[aob.16s.sum$loc == 'CW',]))

## overall summary statistics

### aoa pw vs. aoa cw

leveneTest(sum~loc, aoa.16s.sum3[aoa.16s.sum3$loc == 'CW' | aoa.16s.sum3$loc == 'PW',])
shapiro.test(aoa.16s.sum3[aoa.16s.sum3$loc == 'CW' | aoa.16s.sum3$loc == 'PW',]$sum)

kruskal.test(aoa.16s.sum3[aoa.16s.sum3$loc == 'CW' | aoa.16s.sum3$loc == 'PW',], sum~loc)

### aob pw vs. aob cw

levene_test(aob.16s.sum[aob.16s.sum$loc == 'CW' | aob.16s.sum$loc == 'PW',], tot_abundance~loc)
shapiro.test(aob.16s.sum[aob.16s.sum$loc == 'CW',]$tot_abundance)
shapiro.test(aob.16s.sum[aob.16s.sum$loc == 'PW',]$tot_abundance)

kruskal.test(aob.16s.sum[aob.16s.sum$loc == 'CW' | aob.16s.sum$loc == 'PW',], tot_abundance~loc)

### aoa pw vs. aob pw

aob.16s.sum$Organism <- 'AOB'
aob.16s.sum3 <- aob.16s.sum %>% rename(sum=tot_abundance)
aoa.16s.sum3$Organism <- 'AOA'

pw.16s <- rbind(aoa.16s.sum3[aoa.16s.sum3$loc == 'PW',], aob.16s.sum3[aob.16s.sum3$loc == 'PW',])

pw.16s <- pw.16s %>% ungroup()
levene_test(pw.16s, sum~Organism)
shapiro.test(pw.16s[pw.16s$Organism == 'AOA',]$sum)
shapiro.test(pw.16s[pw.16s$Organism == 'AOB',]$sum)

kruskal.test(pw.16s, sum~Organism)

### aoa cw vs. aob cw

cw.16s <- rbind(aoa.16s.sum3[aoa.16s.sum3$loc == 'CW',], aob.16s.sum2[aob.16s.sum2$loc == 'CW',])
cw.16s <- cw.16s %>% ungroup()
levene_test(cw.16s, sum~Organism)
shapiro.test(cw.16s[cw.16s$Organism == 'AOA',]$sum)
shapiro.test(cw.16s[cw.16s$Organism == 'AOB',]$sum)

kruskal.test(cw.16s, sum~Organism)

## SSU AOA PW

aoa.ssu.sum <- aoa.ssu %>% group_by(Site, Date, Order, Subsample) %>% summarise(tot_abundance = sum(Abundance))

aoa.ssu.sum2 <- aoa.ssu.sum %>% group_by(Site, Date, Order) %>% summarise(mean = mean(tot_abundance), sd = sd(tot_abundance), n = n(), se = sd/sqrt(n)) %>%
  arrange(Order) %>% mutate(cumulative = cumsum(mean))

aoa.ssu.sum2$Date <- factor(aoa.ssu.sum2$Date, levels = c('18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'),
                               labels = c('2018-04', '2018-06', '2018-08', '2018-10', '2018-12', '2019-02'))

aoa.ssu.sum2$Order <- factor(aoa.ssu.sum2$Order, levels = c('Nitrososphaerales', 'Ca. Nitrosopumilales', 'Ca. Nitrosotaleales'))

aoa.ssu.pw <- ggplot(aoa.ssu.sum2[aoa.ssu.sum2$Site == 'PW',], aes(x= Date, y = mean, group = Order)) + 
  geom_col(aes(fill = Order), width = .98) + 
  scale_fill_manual(values = c("#7D3C98", "#A04000", "#117A65")) + 
  scale_color_manual(values = c("#7D3C98", "#A04000", "#117A65")) +
  labs(y = expression(paste(x10^{10}, ' SSU transcripts/g DW soil')), x = '') + 
  scale_y_continuous(limits = c(0, 3.5e10), breaks = c(0, 1e10, 2e10, 3e10), labels = c('0.0', '1.0', '2.0', '3.0'),
                     expand = c(0,0)) +
  scale_x_discrete(drop = FALSE) + 
  guides(col = 'none', fill = 'none') + 
  geom_errorbar(data = aoa.ssu.sum2[aoa.ssu.sum2$Site == 'PW' & aoa.ssu.sum2$Order == 'Nitrososphaerales' | 
                                      aoa.ssu.sum2$Site == 'PW' & aoa.ssu.sum2$Order == 'Ca. Nitrosopumilales' & aoa.ssu.sum2$Date == '2018-10',],
                aes(ymax = cumulative + se, ymin = cumulative-100000000, col = Order), width = 0.1, linewidth = 1) +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text.x=element_text(size=20,face = "bold", angle=45, hjust = 1),
        axis.text.y=element_text(size=20, face='bold'),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

aoa.ssu.pw
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
  theme(axis.title=element_text(size=20,face = "bold"),axis.text.x=element_text(size=20,face = "bold", angle=45, hjust = 1),
        axis.text.y=element_text(size=20, face='bold'),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)


## SSU AOB PW

aob.ssu.sum <- aob.ssu %>% group_by(Site, Date, Sample) %>% summarise(tot_abundance = sum(Abundance))

aob.ssu.sum2 <- aob.ssu.sum %>% group_by(Site, Date) %>% summarise(mean = mean(tot_abundance), sd = sd(tot_abundance), n = n(), se = sd/sqrt(n))

aob.ssu.sum2$Date <- factor(aob.ssu.sum2$Date, levels = c('18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'),
                            labels = c('2018-04', '2018-06', '2018-08', '2018-10', '2018-12', '2019-02'))

aob.ssu.pw <- ggplot(aob.ssu.sum2[aob.ssu.sum2$Site == 'PW',], aes(x= Date, y = mean)) + 
  geom_col(col = '#499DD4', fill = '#499DD4', width = .95) + 
  labs(y = expression(paste(x10^{11}, ' SSU transcripts/g DW soil')), x = '') + 
  scale_y_continuous(limits = c(0, 6.1e11), breaks = c(0, 2e11, 4e11, 6e11), labels = c('0.0', '2.0', '4.0', '6.0'),
                     expand = c(0,0)) +
  guides(fill = 'none', col = 'none') + 
  scale_x_discrete(drop = FALSE) + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean), col = '#499DD4', width = 0.1, linewidth = 1) +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text.x=element_text(size=20,face = "bold", angle=45, hjust = 1),
        axis.text.y=element_text(size=20, face='bold'),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

aob.ssu.pw

aob.ssu.cw <- ggplot(aob.ssu.sum2[aob.ssu.sum2$Site == 'CW',], aes(x= Date, y = mean)) + 
  geom_col(col = '#499DD4', fill = '#499DD4', width = .95) + 
  labs(y = expression(paste(x10^{10}, ' SSU transcripts/g DW soil')), x = '') + 
  scale_y_continuous(limits = c(0, 8.1e10), breaks = c(0, 2e10, 4e10, 6e10, 8e10), labels = c('0.0', '2.0', '4.0', '6.0', '8.0'),
                     expand = c(0,0)) +
  scale_x_discrete(drop = FALSE) + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean), col = '#499DD4', width = 0.1, linewidth = 1) +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text.x=element_text(size=20,face = "bold", angle=45, hjust = 1),
        axis.text.y=element_text(size=20, face='bold'),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

aob.ssu.cw

## statistics

aob.ssu.sum4 <- aob.ssu.sum %>% ungroup()

levene_test(aob.ssu.sum4[aob.ssu.sum4$Site == 'PW',], tot_abundance ~ Date)
shapiro.test(aob.ssu.sum4[aob.ssu.sum4$Site == 'PW' & aob.ssu.sum4$Date == '19-Feb',]$tot_abundance)

summary(aov(tot_abundance ~ Date, data = aob.ssu.sum4[aob.ssu.sum4$Site == 'PW',]))
TukeyHSD(aov(tot_abundance ~ Date, data = aob.ssu.sum4[aob.ssu.sum4$Site == 'PW',]))

levene_test(aob.ssu.sum4[aob.ssu.sum4$Site == 'CW',], tot_abundance ~ Date)
shapiro.test(aob.ssu.sum4[aob.ssu.sum4$Site == 'CW' & aob.ssu.sum4$Date == '18-Aug',]$tot_abundance)

kruskal.test(aob.ssu.sum4[aob.ssu.sum4$Site == 'CW',], tot_abundance ~ Date)
dunn_test(aob.ssu.sum4[aob.ssu.sum4$Site == 'CW',], tot_abundance ~ Date)

levene_test(aob.ssu.sum4, tot_abundance ~ Site)
shapiro.test(aob.ssu.sum4[aob.ssu.sum4$Site == 'PW',]$tot_abundance)
shapiro.test(aob.ssu.sum4[aob.ssu.sum4$Site == 'CW',]$tot_abundance)
kruskal.test(aob.ssu.sum4, tot_abundance ~ Site)

aoa.ssu.sum4 <- aoa.ssu.sum %>% group_by(Site, Date, Subsample) %>% summarise(total = sum(tot_abundance)) %>% mutate(Organism = 'AOA')
aob.ssu.sum4 <- aob.ssu.sum4 %>% mutate(Organism = 'AOB') %>% rename(total = tot_abundance)

pw.ssu.data <- rbind(aoa.ssu.sum4[aoa.ssu.sum4$Site == 'PW',], aob.ssu.sum4[aob.ssu.sum4$Site == 'PW',]) %>% ungroup()
cw.ssu.data <- rbind(aoa.ssu.sum4[aoa.ssu.sum4$Site == 'CW',], aob.ssu.sum4[aob.ssu.sum4$Site == 'CW',]) %>% ungroup()

levene_test(pw.ssu.data, total ~ Organism)
shapiro.test(pw.ssu.data[pw.ssu.data$Organism == 'AOB',]$total)
shapiro.test(pw.ssu.data[pw.ssu.data$Organism == 'AOA',]$total)
kruskal.test(pw.ssu.data, total ~ Organism)

levene_test(cw.ssu.data, total ~ Organism)
shapiro.test(cw.ssu.data[cw.ssu.data$Organism == 'AOB',]$total)
shapiro.test(cw.ssu.data[cw.ssu.data$Organism == 'AOA',]$total)
kruskal.test(cw.ssu.data, total ~ Organism)
## arrange figures

pw.16s <- plot_grid(aoa.16s.pw, aob.16s.pw)
cw.16s <- plot_grid(aoa.16s.cw, aob.16s.cw)

pw.ssu <- plot_grid(aoa.ssu.pw, aob.ssu.pw)
cw.ssu <- plot_grid(aoa.ssu.cw, aob.ssu.cw)
