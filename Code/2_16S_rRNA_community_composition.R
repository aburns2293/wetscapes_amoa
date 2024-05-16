################################################################################
#                                                                              #
# 16S rRNA relative abundance ASVs and community composition                   #
#                                                                              #
# Author: Anna Burns                                                           #
# Last edited: 26.07.2023                                                      #
# Last tested: 25.04.2024                                                      #
#                                                                              #
################################################################################

################################################################################
# Required packages                                                            #
################################################################################

library(tidyverse)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(vegan)
library(tibble)

################################################################################
# Data upload                                                                  #
################################################################################

asv <- read.csv('Data/prok.asv.csv') %>% remove_rownames() %>% column_to_rownames(var = 'ID') 
sample <- read.csv('Data/prok.sample.csv')

################################################################################
# Functions                                                                    #
################################################################################

# calculates nmds data for desired site

nmds.data <- function(asv, sample, site) {
  
  if(site != 'all') {
    
    site.asv <- asv[, str_detect(colnames(asv), site)]
    
    site.asv$row.sums <- rowSums(site.asv)
    
    site.asv.clean <- site.asv %>% filter(row.sums != 0) %>% select(-row.sums)
    
    site.asv.t <- as.data.frame(t(site.asv.clean))
    
    site.asv.rel <- decostand(site.asv.t, method = 'total')
    
    site.mds <- metaMDS(site.asv.t, distance = 'bray', autotransform = FALSE)
    
    return(site.mds)
    
  }else{ 
    
    asv$row.sums <- rowSums(asv)
    
    asv.clean <- asv %>% filter(row.sums != 0) %>% select(-row.sums)
    
    asv.t <- as.data.frame(t(asv.clean))
    
    asv.rel <- decostand(asv.t, method = 'total')
    
    mds <- metaMDS(asv.rel, distance = 'bray', autotransform = FALSE)
    
    return(mds)
    
  }
}

# adds nmds site ordination to sample data

sample.nmds.data <- function(mds, sample, site) {
  
  if(site != 'all') {
    
    site.mds <- mds
    
    site.scores <- scores(site.mds)
    
    site.scores <- as.data.frame(site.scores$sites)
    
    site.scores$Sample <- rownames(site.scores)
    
    site.sample <- sample[sample$loc == site,]
    
    site.sample.new <- left_join(site.sample, site.scores, by = 'Sample')
    
    return(site.sample.new)
    
  }else{ 
    
    scores <- scores(mds)
    
    site.scores <- as.data.frame(scores$sites)
    
    site.scores$Sample <- rownames(site.scores)
    
    sample.new <- left_join(sample, site.scores, by = 'Sample')
    
    return(sample.new)
    
  }
}

################################################################################
# Code                                                                         #
################################################################################

# full data-set community composition

sample.short <- sample %>% filter(season2 == '19-Feb' | year == '18')
sample.short <- sample.short %>% filter(season2 != '18-Feb')
sample.short.shallow <- sample.short %>% filter(depth == '05-10 cm')

sample.short.shallow.ids <- as.vector(sample.short.shallow$Sample) %>% unique()

asv.short.shallow <- asv[,colnames(asv) %in% sample.short.shallow.ids]

nmds <- nmds.data(asv.short.shallow, sample.short.shallow, 'all')


scores <- scores(nmds)

site.scores <- as.data.frame(scores$sites)

site.scores$Sample <- rownames(site.scores)

sample.new <- dplyr::left_join(sample.short.shallow, site.scores, by = 'Sample')

#full.sample <- sample.nmds.data(nmds, sample.short, 'all')

env.var <- sample.short.shallow %>% select(Sample, water, NO3, NH4, P) %>% remove_rownames() %>% column_to_rownames(var = 'Sample')

env <- envfit(nmds, env.var, na.rm = TRUE)

env

scores <- as.data.frame(scores(nmds)$sites)

cont <- as.data.frame(scores(env, "vectors")) * ordiArrowMul(env)
cat <- as.data.frame(scores(env, "factors")) * ordiArrowMul(env)

ggplot(sample.new, aes(x = NMDS1, y = NMDS2)) + 
  stat_ellipse(aes(group = loc, col = loc), linewidth = 2, linetype = 'dashed') +
  geom_point(aes(col = as.factor(loc)), size = 5) +
  theme_classic() +
  ggtitle('') +  
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E6AB02")) + 
  scale_y_continuous(limits = c(-2,2)) + 
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=15,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, linewidth=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.95,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank())  

## statistics

asv.short.shallow$row.sums <- rowSums(asv.short.shallow)

asv.clean <- asv.short.shallow %>% filter(row.sums != 0) %>% select(-row.sums)

asv.t <- as.data.frame(t(asv.clean))

asv.rel <- decostand(asv.t, method = 'total')

dist <- vegdist(asv.rel, method = "bray")

beta <- betadisper(dist, group = sample.short.shallow$loc)

anova(beta)
boxplot(beta)

TukeyHSD(beta)

perm <- adonis2(dist ~ sample.short.shallow$loc, permutations = 999)

perm

# PW ordination with environmental variables

pw.asv <- asv.short.shallow[,str_detect(colnames(asv.short.shallow), 'P') == TRUE]

pw.asv <- pw.asv[!apply(pw.asv == "0", 1, all),] 

pw.short.shallow <- sample.short.shallow %>% filter(loc == 'PW' | loc == 'PD')

pw.mds <- nmds.data(pw.asv, pw.short.shallow, 'all')

pw.sample <- sample.nmds.data(pw.mds, pw.short.shallow, 'all')

pw.env.var <- pw.sample %>% select(Sample, NO3, NH4, P, water, TOM, Drought_Status) %>% remove_rownames() %>% column_to_rownames(var = 'Sample')

pw.env <- envfit(pw.mds, pw.env.var, na.rm = TRUE)

pw.env

pw.scores <- as.data.frame(scores(pw.mds)$sites)

pw.cont <- as.data.frame(scores(pw.env, "vectors")) * ordiArrowMul(pw.env)
pw.cat <- as.data.frame(scores(pw.env, "factors")) * ordiArrowMul(pw.env)

pw.dist <- as.matrix(vegdist(pw.asv, method = "bray"))

set.seed(111)

pw.clus <- clusGap(pw.dist, kmeans, K.max = 10, B = 500)

print(plot(pw.clus))

set.seed(112)

pw.km <- kmeans(pw.dist, 4, nstart = 25)

pw.scores$cluster <- as.factor(pw.km$cluster)

sample$clusters <- as.factor(km$cluster)

ggplot(pw.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point() #+
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = pw.cont) +
  geom_point(data = pw.cat, aes(x = NMDS1, y = NMDS2), 
             shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  geom_text_repel(data = pw.cat, aes(x = NMDS1, y = NMDS2+0.04), 
                  label = row.names(pw.cat), colour = "navy", fontface = "bold") + 
  geom_text_repel(data = pw.cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
                  fontface = "bold", label = row.names(pw.cont)) + 
  theme_classic() +
  labs(col = 'Depth') +
  ggtitle('PW: 16S ASV NMDS with Envfit')

# PD ordination with environmental variables

pd.asv <- absolute.asv[str_detect(row.names(absolute.asv), 'PD') == TRUE,]

pd.mds <- metaMDS(pd.asv, distance = 'bray', autotransform = FALSE)

pd.sample <- sample.nmds.data(pd.mds, sample.short[sample.short$loc == 'PD',], 'PD')

pd.env.var <- pd.sample %>% select(Sample, NO3, NH4, Drought_Status) %>% remove_rownames() %>% column_to_rownames(var = 'Sample')

pd.env <- envfit(pd.mds, pd.env.var, na.rm = TRUE)

pd.env

pd.scores <- as.data.frame(scores(pd.mds)$sites)

pd.cont <- as.data.frame(scores(pd.env, "vectors")) * ordiArrowMul(pd.env)
pd.cat <- as.data.frame(scores(pd.env, "factors")) * ordiArrowMul(pd.env)

pd.dist <- as.matrix(vegdist(pd.asv, method = "bray"))

set.seed(211)

pd.clus <- clusGap(pd.dist, kmeans, K.max = 10, B = 500)

print(plot(pd.clus))

set.seed(212)

pd.km <- kmeans(pd.dist, 3 , nstart = 25)

pd.scores$cluster <- as.factor(pd.km$cluster)

ggplot(pd.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(col = as.factor(cluster))) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = pd.cont) +
  geom_point(data = pd.cat, aes(x = NMDS1, y = NMDS2), 
             shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  geom_text_repel(data = pd.cat, aes(x = NMDS1, y = NMDS2+0.04), 
                  label = row.names(pd.cat), colour = "navy", fontface = "bold") + 
  geom_text_repel(data = pd.cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
                  fontface = "bold", label = row.names(pd.cont)) + 
  theme_classic() +
  labs(col = 'Depth') +
  ggtitle('PD: 16S ASV NMDS with Envfit')

# CW ordination with environmental variables

cw.nmds <- nmds.data(asv, sample, 'CW')

cw.sample <- sample.nmds.data(cw.nmds, sample, 'CW')

cw.env.var <- cw.sample %>% select(Sample, depth, NO3, NH4, Drought_Status) %>% remove_rownames() %>% column_to_rownames(var = 'Sample')

cw.env <- envfit(cw.nmds, cw.env.var, na.rm = TRUE)

cw.env

cw.scale <- ordiArrowMul(cw.env)

cw.scores <- as.data.frame(scores(cw.nmds)$sites)

cw.scores$depth <- cw.sample$depth

cw.cont <- as.data.frame(scores(cw.env, "vectors"))*cw.scale
cw.cat <- as.data.frame(scores(cw.env, "factors"))*cw.scale

ggplot(cw.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(col = as.factor(depth))) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = cw.cont) +
  geom_point(data = cw.cat, aes(x = NMDS1, y = NMDS2), 
             shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  geom_text_repel(data = cw.cat, aes(x = NMDS1, y = NMDS2+0.04), 
                  label = row.names(cw.cat), colour = "navy", fontface = "bold") + 
  geom_text_repel(data = cw.cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
                  fontface = "bold", label = row.names(cw.cont)) + 
  theme_classic() +
  labs(col = 'Depth') +
  ggtitle('CW: 16S ASV NMDS with Envfit')

# CD ordination with environmental variables

cd.nmds <- nmds.data(asv, sample, 'CD')

cd.sample <- sample.nmds.data(cd.nmds, sample, 'CD')

cd.env.var <- cd.sample %>% select(Sample, depth, NO3, NH4, Drought_Status) %>% remove_rownames() %>% column_to_rownames(var = 'Sample')

cd.env <- envfit(cd.nmds, cd.env.var, na.rm = TRUE)

cd.env

cd.scores <- as.data.frame(scores(cd.nmds)$sites)

cd.depth <- cd.sample %>% remove_rownames() %>% column_to_rownames(var = 'Sample') %>% select(depth)

cd.scores.2 <- merge(cd.scores, cd.depth, by='row.names')

cd.cont <- as.data.frame(scores(cd.env, "vectors")) * ordiArrowMul(cd.env)
cd.cat <- as.data.frame(scores(cd.env, "factors")) * ordiArrowMul(cd.env)

ggplot(cd.scores.2, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(col = as.factor(depth))) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = cd.cont) +
  geom_point(data = cd.cat, aes(x = NMDS1, y = NMDS2), 
             shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  geom_text_repel(data = cd.cat, aes(x = NMDS1, y = NMDS2+0.04), 
                  label = row.names(cd.cat), colour = "navy", fontface = "bold") + 
  geom_text_repel(data = cd.cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
                  fontface = "bold", label = row.names(cd.cont)) + 
  theme_classic() +
  labs(col = 'Depth') +
  ggtitle('CD: 16S ASV NMDS with Envfit')


# sample data

sample.short <- sample %>% filter(year == '18' | season2 == '19-Feb') %>% filter(depth == '05-10 cm') 

sample.short$season2 <- factor(sample.short$season2, levels = c('18-Feb', '18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'))

sample.short[sample.short$group == '18OctPD1',]$NO3 <- 296.74
sample.short[sample.short$group == '18OctPD2',]$NO3 <- 373.76
sample.short[sample.short$group == '18OctPD3',]$NO3 <- 787.33

sample.nitrogen <- sample.short %>% select(loc, season2, NH4, NO3)

new_row1 <- c(loc = 'PD', season2 = '18-Dec', NH4 = 6.1, NO3 = 170.9)
new_row2 <- c(loc = 'PD', season2 = '18-Dec', NH4 = 5.9, NO3 = 102.5)
new_row3 <- c(loc = 'PD', season2 = '18-Dec', NH4 = 4.4, NO3 = 193.1)

sample.nitrogen <- rbind(sample.nitrogen, new_row1, new_row2, new_row3)

sample.nitrogen <- sample.nitrogen %>% filter(season2 != '18-Feb')

sample_nitrogen_long <- sample.nitrogen %>% pivot_longer(!c(loc, season2), names_to = 'Nutrient', values_to = 'Volume')

sample_nitrogen_long$Volume <- as.numeric(sample_nitrogen_long$Volume)

sample_nitrogen_summary <- sample_nitrogen_long %>% group_by(loc, season2, Nutrient) %>% summarise(mean = mean(Volume), n = n(), sd = sd(Volume), se = sd/n)

sample_nitrogen_summary$date <- sample_nitrogen_summary$season2

sample_nitrogen_summary[sample_nitrogen_summary$season2 == '18-Feb',]$date <- as.factor('2018-02')

ggplot(sample_nitrogen_summary[sample_nitrogen_summary$loc == 'CW',], aes(x = season2, y = mean, group = Nutrient)) +
  geom_line(aes(linetype = Nutrient), linewidth = 1.2) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean), width = 0.3, linewidth = 1) +
  scale_linetype_manual(values = c('solid', 'dotdash')) +
  xlab('Time') +
  ylab('mg N/g DW') + 
  theme(axis.title=element_text(size=20,face = "bold"),axis.text=element_text(size=15,face = "bold"),
        title = element_text(size = 20, face = 'bold'),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, linewidth=2),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),legend.position = c(0.9,0.9),
        legend.text = element_text(size = 15),legend.key = element_rect(fill = NA),legend.title = element_blank(),
        legend.background = element_blank(),
        strip.background = element_rect(fill = NA, color = 'black', linewidth = 2),
        strip.text = element_text(size = 20, face = 'bold'))

