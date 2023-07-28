################################################################################
#                                                                              #
# 16S rRNA relative abundance ASVs and community composition                   #
#                                                                              #
# Author: Anna Burns                                                           #
# Last edited: 26.07.2023                                                      #
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

################################################################################
# Data upload                                                                  #
################################################################################

asv <- read.csv('Data/prok.asv.csv') %>% remove_rownames() %>% column_to_rownames('X') 

sample <- read.csv('Data/prok.sample.csv') %>% select(-X)

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
    
    site.mds <- metaMDS(site.asv.rel, distance = 'bray', autotransform = FALSE)
    
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

nmds <- nmds.data(asv, sample, 'all')

full.sample <- sample.nmds.data(nmds, sample, 'all')

ggplot(full.sample, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(col = as.factor(loc))) +
  theme_classic() +
  labs(col = 'Site') +
  ggtitle('All sites: 16S ASV Community Composition')

## statistics

asv$row.sums <- rowSums(asv)

asv.clean <- asv %>% filter(row.sums != 0) %>% select(-row.sums)

asv.t <- as.data.frame(t(asv.clean))

asv.rel <- decostand(asv.t, method = 'total')

dist <- vegdist(asv.rel, method = "bray")

beta <- betadisper(dist, group = sample$loc)

anova(beta) 

perm <- adonis2(dist ~ sample$loc, permutations = 999)

perm

# PW ordination with environmental variables

pw.nmds <- nmds.data(asv, sample, 'PW')

pw.sample <- sample.nmds.data(pw.nmds, sample, 'PW')

pw.env.var <- pw.sample %>% select(Sample, depth, NO3, NH4, Drought_Status) %>% remove_rownames() %>% column_to_rownames(var = 'Sample')

pw.env <- envfit(pw.nmds, pw.env.var, na.rm = TRUE)

pw.env

pw.scores <- as.data.frame(scores(pw.nmds)$sites)

pw.scores$depth <- pw.sample$depth

pw.cont <- as.data.frame(scores(pw.env, "vectors")) * ordiArrowMul(pw.env)
pw.cat <- as.data.frame(scores(pw.env, "factors")) * ordiArrowMul(pw.env)

ggplot(pw.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(col = as.factor(depth))) +
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

pd.nmds <- nmds.data(asv, sample, 'PD')

pd.sample <- sample.nmds.data(pd.nmds, sample, 'PD')

pd.env.var <- pd.sample %>% select(Sample, depth, NO3, NH4, Drought_Status) %>% remove_rownames() %>% column_to_rownames(var = 'Sample')

pd.env <- envfit(pd.nmds, pd.env.var, na.rm = TRUE)

pd.env

pd.scores <- as.data.frame(scores(pd.nmds)$sites)

pd.scores$depth <- pd.sample$depth

pd.cont <- as.data.frame(scores(pd.env, "vectors")) * ordiArrowMul(pd.env)
pd.cat <- as.data.frame(scores(pd.env, "factors")) * ordiArrowMul(pd.env)

ggplot(pd.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(col = as.factor(depth))) +
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

