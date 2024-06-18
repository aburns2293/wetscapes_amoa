
################################################################################
#                                                                              #
# Defining drought thresholds via k-means clustering                           #
#                                                                              #
# Author: Anna Burns                                                           #
# Last edited: 26.07.2023                                                      #
# Last tested: 25.04.2024                                                      #
#                                                                              #
################################################################################

################################################################################
# Required packages                                                            #
################################################################################

library(ggplot2)
library(dplyr)
library(tidyverse)
library(cluster)
library(stringr)
library(lubridate)
library(rstatix)
library(readr)

################################################################################
# Data upload                                                                  #
################################################################################

# daily water table measurements from 06.11.2017 to 18.02.2019

water_table_summary <- read.csv("Data/water_table_summary.csv")

water_table_summary$Year <- as.numeric(water_table_summary$Year)

# sample data

sample_data <- read.csv('Data/prok.sample.csv')

# precipitation

greifswald_precip <- read_delim("Data/precipitation.data.txt", 
                                delim = ";", escape_double = FALSE, trim_ws = TRUE)
# temperature

greifswald_temp <- read_delim("../wetscapes_amoa/Data/temperature.data.txt", delim = ";", escape_double = FALSE, trim_ws = TRUE)

################################################################################
# Functions                                                                    #
################################################################################

# evaluates the optimal number of clusters for k-means algorithm

wt.cluster.optimal_number <- function(wt_data, site){
  
  site_wt <- wt_data %>% filter(Site == site) %>% select(GW_level) %>% na.omit()
  
  site_wt_matrix <- as.matrix(site_wt)
  
  site_clus <- clusGap(site_wt_matrix, kmeans, K.max = 10, B = 500)
  
  print(plot(site_clus)) 
  
  return(site_clus)
  }

# takes input of optimal cluster numbers and runs clustering algorithm

wt.cluster.full <- function(wt_data, optimal_number, site) {
  
  site_opt_num <- as.integer(optimal_number)
  
  site_wt <- wt_data %>% filter(Site == site) %>% select(GW_level) %>% na.omit()
  
  site_wt_matrix <- as.matrix(site_wt)
  
  set.seed(11)
  
  site_wt_km <- kmeans(site_wt_matrix, site_opt_num)
  
  site_wt$clusters <- as.factor(site_wt_km$cluster)
  
  return(site_wt) }

################################################################################
# Code                                                                         #
################################################################################

# weather - precipitation and temperature

## weather data has already been added to the 'prok.sample.csv' file, so code here
## is not necessary beyond reproducibility purposes

## cleaning water table data

water_table_summary$Date <- mdy(water_table_summary$Date)
  
## cleaning precipitation data

greifswald_precip$MESS_DATUM <- as.character(greifswald_precip$MESS_DATUM)

greifswald_precip$MESS_DATUM <- as.Date(greifswald_precip$MESS_DATUM, format = '%Y%m%d')

greifswald_precip$JAHR <- as.numeric(format(greifswald_precip$MESS_DATUM, "%Y"))

greifswald_precip$MONAT <- as.numeric(format(greifswald_precip$MESS_DATUM, "%m"))

greifswald_precip$TAG <- as.numeric(format(greifswald_precip$MESS_DATUM, "%d"))

## cleaning temperature data

greifswald_temp$MESS_DATUM <- as.character(greifswald_temp$MESS_DATUM)

greifswald_temp$MESS_DATUM <- substr(greifswald_temp$MESS_DATUM, 1, nchar(greifswald_temp$MESS_DATUM)-2)

greifswald_temp$MESS_DATUM <- as.Date(greifswald_temp$MESS_DATUM, format = '%Y%m%d')

greifswald_temp <- greifswald_temp %>% group_by(MESS_DATUM) %>% summarise(MEAN_TEMP = mean(TT_TER), MEAN_HUM = mean(RF_TER))

## combining precipitation and temperature data (21st century)

greifswald_weather <- left_join(greifswald_precip, greifswald_temp, by = 'MESS_DATUM')

greifswald_weather <- greifswald_weather %>% filter(JAHR >= 2000)

## monthly weather averages

greifswald_weather_sum <- greifswald_weather %>% group_by(JAHR, MONAT) %>% 
  summarise(mean_temp = mean(MEAN_TEMP), sd_temp = sd(MEAN_TEMP),
            mean_hum = mean(MEAN_HUM), sd_hum = sd(MEAN_HUM), mean_precip = mean(RS), sd_precip = sd(RS)) %>%
  arrange(MONAT)

# Clustering water table values

## PW drought thresholds

pw_wt_opt <- wt.cluster.optimal_number(water_table_summary, 'PW')

pw.k <- as.data.frame(pw_wt_opt[["Tab"]])

pw.k$k <- row.names(pw.k)

pw.k$k <- factor(pw.k$k, levels = seq(1, 10, 1))

ggplot(pw.k, aes(x = k, y = gap, group = 1)) + geom_line() + geom_point() + 
  ylab('Gap(k)') +
  xlab('Clusters k') +
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=20,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 20),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1) + 
  geom_point(aes(x = 2, y = pw.k[pw.k$k == 2,]$gap), shape = 18, size = 10, col = '#7D3C98')

ggplot(water_table_summary[water_table_summary$Site == 'PW',]) + geom_histogram(aes(x = GW_level)) + geom_vline(aes(xintercept = -5.45), linetype = 2, size = 2) + 
  ggtitle('a) PW') + 
  xlab('Groundwater Depth (cm)') + 
  ylab('Count') + 
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=15,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

### identify the 'best' cluster from this figure and enter it in the wt.cluster.full function
### 'best' - balance between minimal cluster silhouette and minimal number of clusters for interpretation

pw_wt <- wt.cluster.full(water_table_summary, 2, 'PW')

ggplot(pw_wt) + geom_boxplot(aes(x = clusters, y = GW_level))

pw_wt$Drought <- 'NA'

### before entering these variables, make sure that the cluster with the higher average water table in the plot corresponds to 'non-drought'

pw_wt$Drought_Status[pw_wt$clusters == 1] <- 'Non-Drought'
pw_wt$Drought_Status[pw_wt$clusters == 2] <- 'Drought'

### the drought threshold is the maximum water table level during drought conditions

pw_wt %>% group_by(Drought_Status) %>% summarise(max = max(GW_level), min = min(GW_level), mean = mean(GW_level))

## PD drought thresholds

pd_wt_opt <- wt.cluster.optimal_number(water_table_summary, 'PD')

pd_wt <- wt.cluster.full(water_table_summary, 2, 'PD')

ggplot(pd_wt) + geom_boxplot(aes(x = clusters, y = GW_level))

ggplot(water_table_summary[water_table_summary$Site == 'PD',]) + geom_histogram(aes(x = GW_level)) + geom_vline(aes(xintercept = -31.05), linetype = 2, size = 2) + 
  ggtitle('b) PD') + 
  xlab('Groundwater Depth (cm)') + 
  ylab('Count') + 
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=15,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

pd_wt <- pd_wt %>% add_column(Drought_Status = NA)

pd.k <- as.data.frame(pd_wt_opt[["Tab"]])

pd.k$k <- row.names(pd.k)

pd.k$k <- factor(pd.k$k, levels = seq(1, 10, 1))

ggplot(pd.k, aes(x = k, y = gap, group = 1)) + geom_line() + geom_point() + 
  ylab('Gap') +
  ggtitle('b) PD') + 
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=15,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank()) + 
  geom_point(aes(x = 2, y = pd.k[pd.k$k == 2,]$gap), shape = 18, size = 5, size = 2)

### before entering these variables, make sure that the cluster with the higher average water table in the plot corresponds to 'non-drought'

pd_wt$Drought_Status[pd_wt$clusters == 2] <- 'Non-Drought'
pd_wt$Drought_Status[pd_wt$clusters == 1] <- 'Drought'

pd_wt %>% group_by(Drought_Status) %>% summarise(max = max(GW_level), min = min(GW_level), mean = mean(GW_level))

## CW drought thresholds

cw_wt_opt <- wt.cluster.optimal_number(water_table_summary, 'CW')

cw_wt <- wt.cluster.full(water_table_summary, 2, 'CW')

ggplot(cw_wt) + geom_boxplot(aes(x = clusters, y = GW_level))

cw_wt <- cw_wt %>% add_column(Drought_Status = NA)

cw.k <- as.data.frame(cw_wt_opt[["Tab"]])

cw.k$k <- row.names(cw.k)

cw.k$k <- factor(cw.k$k, levels = seq(1, 10, 1))

ggplot(cw.k, aes(x = k, y = gap, group = 1)) + geom_line() + geom_point() + 
  ylab('Gap') +
  ggtitle('c) CW') + 
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=15,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank()) + 
  geom_point(aes(x = 2, y = cw.k[cw.k$k == 2,]$gap), shape = 18, size = 5, size = 2)

ggplot(water_table_summary[water_table_summary$Site == 'CW',]) + geom_histogram(aes(x = GW_level)) + geom_vline(aes(xintercept = -30.24), linetype = 2, size = 2) + 
  ggtitle('c) CW') + 
  xlab('Groundwater Depth (cm)') + 
  ylab('Count') + 
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=15,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

### before entering these variables, make sure that the cluster with the higher average water table in the plot corresponds to 'non-drought'

cw_wt$Drought_Status[cw_wt$clusters == 1] <- 'Non-Drought'
cw_wt$Drought_Status[cw_wt$clusters == 2] <- 'Drought'

cw_wt %>% group_by(Drought_Status) %>% summarise(max = max(GW_level), min = min(GW_level), mean = mean(GW_level))

## CD drought thresholds

cd_wt_opt <- wt.cluster.optimal_number(water_table_summary, 'CD')

cd_wt <- wt.cluster.full(water_table_summary, 2, 'CD')

ggplot(cd_wt) + geom_boxplot(aes(x = clusters, y = GW_level))

cd_wt <- cd_wt %>% add_column(Drought_Status = NA)

cd.k <- as.data.frame(cd_wt_opt[["Tab"]])

cd.k$k <- row.names(cd.k)

cd.k$k <- factor(cd.k$k, levels = seq(1, 10, 1))

ggplot(cd.k, aes(x = k, y = gap, group = 1)) + geom_line() + geom_point() + 
  ylab('Gap') +
  ggtitle('d) CD') + 
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=15,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank()) + 
  geom_point(aes(x = 2, y = cd.k[cd.k$k == 2,]$gap), shape = 18, size = 5, size = 2)

ggplot(water_table_summary[water_table_summary$Site == 'CD',]) + geom_histogram(aes(x = GW_level)) + geom_vline(aes(xintercept = -45.88), linetype = 2, size = 2) + 
  ggtitle('d) CD') + 
  xlab('Groundwater Depth (cm)') + 
  ylab('Count') + 
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=15,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())

### before entering these variables, make sure that the cluster with the higher average water table in the plot corresponds to 'non-drought'

cd_wt$Drought_Status[cd_wt$clusters == 2] <- 'Non-Drought'
cd_wt$Drought_Status[cd_wt$clusters == 1] <- 'Drought'

cd_wt %>% group_by(Drought_Status) %>% summarise(max = max(GW_level), min = min(GW_level), mean = mean(GW_level))

## summarize and combine drought data

full_drought_wt <- rbind(pw_wt, pd_wt, cd_wt, cw_wt) %>% select(-clusters)

## add site

full_drought_wt <- full_drought_wt %>% mutate(Site = row.names(full_drought_wt))

full_drought_wt$Site <- str_sub(full_drought_wt$Site, 0, 2)

## add dates and format

full_drought_wt <- full_drought_wt %>% mutate(Date = row.names(full_drought_wt))

full_drought_wt$Date <- str_sub(full_drought_wt$Date,6)

full_drought_wt$Date <- str_c('20',full_drought_wt$Date)

full_drought_wt$Date <- as.Date(full_drought_wt$Date)

################################################################################
# Figures and Statistics                                                       #
################################################################################

# Groundwater table depth with drought thresholds

water_table_summary$Date <- mdy(water_table_summary$Date)
drought2018 <- interval(start = "2018-03-15", end = "2019-02-28")
water_table_drought2018 <- water_table_summary[which(water_table_summary$Date %within% drought2018),]

ggplot(water_table_drought2018[water_table_drought2018$Site == 'CW' | water_table_drought2018$Site == 'PW', ]) + 
  geom_line(aes(x = Date, y = GW_level, col = Site), size = 1.3) +
  theme_classic() +
  ylab('Depth below surface (cm)') +
  xlab('') +
  geom_hline(yintercept = -5.45, col = "#117A65", linetype = 'dashed', size = 1) +
  geom_hline(yintercept = -30.24, col = "#A04000", linetype = 'dashed', size = 1) +
  geom_hline(yintercept = 0, col = "black", linetype = 'solid', size = 1, alpha = 0.3) +
  #scale_y_continuous(limits = c(-30, 15)) + 
  scale_color_manual(values = c("#A04000", "#117A65")) + 
  scale_x_date(breaks = seq(as.Date('2017-12-01'), as.Date('2020-02-21'), by = "2 months"), date_labels = '%Y-%m') +
  theme(axis.title=element_text(size=24,face = "bold"),
        axis.text=element_text(size=20,face = "bold"),
        axis.text.x = element_text(angle=45, vjust = 0.6),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, linewidth=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1) + 
  guides(color = guide_legend(list(linetype = c(1, 1, 1, 1) ) ) )

# Relationship between groundwater depth and drought status

wilcox.test(GW_level ~ Drought_Status, data = pw_wt)
wilcox.test(GW_level ~ Drought_Status, data = pd_wt)
wilcox.test(GW_level ~ Drought_Status, data = cw_wt)
wilcox.test(GW_level ~ Drought_Status, data = cd_wt)

ggplot(full_drought_wt) +
  geom_boxplot(aes(x = Site, y = GW_level, col = Drought_Status)) +
  theme_classic() +
  ggtitle('Relationship between drought status and water table depth') +
  ylab('Water table depth from surface (cm)') +
  xlab('') +
  annotate('text', x=c(1, 2, 2.9, 3.9), y=c(5,48,10,20), label='***') +
  labs(color = NULL)

# Topsoil water content

### separate topsoil data

sample_topsoil <- sample_data %>% filter(depth == '05-10 cm',
                                         loc == 'PW' | loc == 'CW',
                                         year == '18' | season2 == '19-Feb',
                                         season2 != '18-Feb') %>% select(water, loc, season2)

sample_topsoil$season2 <- factor(sample_topsoil$season2, levels = c('18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'), labels = c('2018-04', '2018-06', '2018-08', '2018-10', '2018-12', '2019-02'))

sample_topsoil_sum <- sample_topsoil %>% group_by(loc, season2) %>% summarise(mean = mean(water), sd = sd(water), n = n(), se = sd/sqrt(n))

ggplot(sample_topsoil_sum, aes(x=season2, y = mean, col = loc, group=loc)) +
  geom_point(aes(shape = loc), size = 10) + 
  scale_y_continuous(limits = c(30,100)) +
  geom_line(size = 1.3) + 
  theme_classic() +
  ylab('Relative water weight (%)') +
  xlab('') + 
  scale_shape_manual(values = c(15,16)) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean), width = 0.1, size = 1) +
  scale_color_manual(values = c("#A04000", "#117A65")) + 
  theme(axis.title=element_text(size=24,face = "bold"),
        axis.text.x = element_text(angle=45, vjust = 0.6),
        axis.text=element_text(size=20,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, linewidth=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        aspect.ratio = 1)

levene_test(formula = water~season2, data = sample_topsoil[sample_topsoil$loc == 'PW',])
shapiro_test(sample_topsoil[sample_topsoil$loc == 'PW',]$water)
summary(aov(water~season2, sample_topsoil[sample_topsoil$loc == 'PW',]))
mean(sample_topsoil[sample_topsoil$loc == 'PW',]$water)

levene_test(formula = water~season2, data = sample_topsoil[sample_topsoil$loc == 'CW',])
shapiro_test(sample_topsoil[sample_topsoil$loc == 'CW',]$water)
kruskal_test(formula = water~season2, data = sample_topsoil[sample_topsoil$loc == 'CW',])
dunn_test(formula = water~season2, data = sample_topsoil[sample_topsoil$loc == 'CW',], p.adjust.method = 'BH')

# Topsoil water content and drought

wilcox.test(water ~ Drought_Status, data = sample_topsoil[sample_topsoil$loc == 'PW',])
wilcox.test(water ~ Drought_Status, data = sample_topsoil[sample_topsoil$loc == 'PD',])
wilcox.test(water ~ Drought_Status, data = sample_topsoil[sample_topsoil$loc == 'CD',])
wilcox.test(water ~ Drought_Status, data = sample_topsoil[sample_topsoil$loc == 'CW',])

ggplot(sample_topsoil, aes(x=loc, y=water, col = Drought_Status)) +
  geom_boxplot() +
  theme_classic() +
  labs(color = NULL) +
  xlab('') +
  ylab('Relative water weight (%)') +
  annotate('text', x=c(1, 2, 2.9, 3.9), y=c(78,78,95,100), label=c('p=0.03','p=0.12','p=0.01','p=0.28')) +
  ggtitle('Topsoil water content during drought phases')
