
################################################################################
#                                                                              #
# Defining drought thresholds via k-means clustering                           #
#                                                                              #
# Author: Anna Burns                                                           #
# Last edited: 26.07.2023                                                      #
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

################################################################################
# Data upload                                                                  #
################################################################################

# daily water table measurements from 06.11.2017 to 18.02.2019

water_table_summary <- read_csv("Data/water_table_summary.csv") %>% select(-'...1') %>% 
  remove_rownames() %>% column_to_rownames(var='Group.1') %>%
  filter(Site != 'AW', Site != 'AD') %>% separate(Date, c('Year', 'Month', 'Date'), '-' )

water_table_summary$Year <- as.numeric(water_table_summary$Year)

# sample data

sample_data <- read.csv('Data/prok.sample.csv') %>% select(-X)

# precipitation

greifswald_precip <- read_delim("Data/precipitation.data.txt", 
                                delim = ";", escape_double = FALSE, trim_ws = TRUE)
# temperature

greifswald_temp <- read_delim("Data/temperature.data.txt", delim = ";", escape_double = FALSE, trim_ws = TRUE)

################################################################################
# Functions                                                                    #
################################################################################

# evaluates the optimal number of clusters for k-means algorithm

wt.cluster.optimal_number <- function(wt_data, site){
  
  site_wt <- wt_data %>% filter(Site == site) %>% select(GW_level) %>% na.omit()
  
  site_wt_matrix <- as.matrix(site_wt)
  
  site_clus <- clusGap(site_wt_matrix, kmeans, K.max = 10, B = 500)
  
  print(plot(site_clus)) }

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

### identify the 'best' cluster from this figure and enter it in the wt.cluster.full function
### 'best' - balance between minimal cluster silhouette and minimal number of clusters for interpretation

pw_wt <- wt.cluster.full(water_table_summary, 2, 'PW')

ggplot(pw_wt) + geom_boxplot(aes(x = clusters, y = GW_level))

pw_wt <- pw_wt %>% add_column(Drought_Status = NA)

### before entering these variables, make sure that the cluster with the higher average water table in the plot corresponds to 'non-drought'

pw_wt$Drought_Status[pw_wt$clusters == 1] <- 'Non-Drought'
pw_wt$Drought_Status[pw_wt$clusters == 2] <- 'Drought'

### the drought threshold is the maximum water table level during drought conditions

pw_wt %>% group_by(Drought_Status) %>% summarise(max = max(GW_level), min = min(GW_level), mean = mean(GW_level))

## PD drought thresholds

pd_wt_opt <- wt.cluster.optimal_number(water_table_summary, 'PD')

pd_wt <- wt.cluster.full(water_table_summary, 2, 'PD')

ggplot(pd_wt) + geom_boxplot(aes(x = clusters, y = GW_level))

pd_wt <- pd_wt %>% add_column(Drought_Status = NA)

### before entering these variables, make sure that the cluster with the higher average water table in the plot corresponds to 'non-drought'

pd_wt$Drought_Status[pd_wt$clusters == 2] <- 'Non-Drought'
pd_wt$Drought_Status[pd_wt$clusters == 1] <- 'Drought'

pd_wt %>% group_by(Drought_Status) %>% summarise(max = max(GW_level), min = min(GW_level), mean = mean(GW_level))

## CW drought thresholds

cw_wt_opt <- wt.cluster.optimal_number(water_table_summary, 'CW')

cw_wt <- wt.cluster.full(water_table_summary, 2, 'CW')

ggplot(cw_wt) + geom_boxplot(aes(x = clusters, y = GW_level))

cw_wt <- cw_wt %>% add_column(Drought_Status = NA)

### before entering these variables, make sure that the cluster with the higher average water table in the plot corresponds to 'non-drought'

cw_wt$Drought_Status[cw_wt$clusters == 1] <- 'Non-Drought'
cw_wt$Drought_Status[cw_wt$clusters == 2] <- 'Drought'

cw_wt %>% group_by(Drought_Status) %>% summarise(max = max(GW_level), min = min(GW_level), mean = mean(GW_level))

## CD drought thresholds

cd_wt_opt <- wt.cluster.optimal_number(water_table_summary, 'CD')

cd_wt <- wt.cluster.full(water_table_summary, 2, 'CD')

ggplot(cd_wt) + geom_boxplot(aes(x = clusters, y = GW_level))

cd_wt <- cd_wt %>% add_column(Drought_Status = NA)

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

ggplot(full_drought_wt) + 
  geom_line(aes(x = Date, y = GW_level, col = Site)) +
  theme_classic() +
  ylab('Water table depth from surface (cm)') +
  xlab('') +
  geom_hline(yintercept = -5.45, col = 'purple', linetype = 'dashed') +
  geom_hline(yintercept = -31.05, col = '#00BFC4', linetype = 'dashed') +
  geom_hline(yintercept = -30.24, col = '#7CAE00', linetype = 'dashed') +
  geom_hline(yintercept = -45.88, col = 'red', linetype = 'dashed') +
  ggtitle('Daily water table depth values with drought thresholds, 2017-09 to 2020-02') +
  scale_y_continuous(breaks = seq(-100,50,25))

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

# Topsoil water content and groundwater depth

### separate topsoil data

sample_topsoil <- sample_data %>% filter(depth == '05-10 cm')

ggplot(sample_topsoil, aes(x=Mean_GW, y = water, col = loc, group=loc)) +
  geom_point() +
  geom_smooth(se=FALSE) +
  theme_classic() +
  ylab('Relative water weight (%)') +
  xlab('Water table depth from surface (cm)') +
  labs(color = 'Site') +
  ggtitle('Relationship between water table depth and topsoil water content')

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
