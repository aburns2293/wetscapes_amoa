


rel.ab <- read.csv('Data/prok.asv.csv')

aoa.clades <- read.csv('Data/16SrRNA.aoa.clades.csv')

sample <- read.csv('Data/prok.sample.csv')

tax <- read.csv('Data/prok.tax.csv', row.names = 1)

rel.ab <- rel.ab %>% remove_rownames() %>% column_to_rownames(var = 'X')

sample <- sample %>% select(-X) %>% remove_rownames() %>% column_to_rownames(var = 'Sample')

asv.m <- data.matrix(rel.ab)

tax.m <- as.matrix(tax)

OTU <- otu_table(asv.m, taxa_are_rows = TRUE)

TAX <- tax_table(tax.m)

INFO <- sample_data(sample)

ps <- phyloseq(OTU, TAX, INFO)

ps.relab <- ps %>% transform_sample_counts(function(x) {x/sum(x)}) 

ps.melt <- psmelt(ps.relab)





aoa.clades <- aoa.clades %>% select(X, Assigned_Clade_Tree2)
aoa.asvs <- as.vector(aoa.clades$X)

rel.ab.aoa <- ps.melt[ps.melt$OTU %in% aoa.asvs,]


aoa.full <- dplyr::left_join(rel.ab.aoa, aoa.clades, by = join_by('OTU' == 'X'))








aoa.rel.ab <- read.csv('Data/16SrRNA.aoa.relative.csv')

aoa.rel.shallow <- aoa.rel.ab %>% filter(depth == '05-10 cm') 

aoa.rel.drought.shallow <- aoa.rel.shallow %>% filter(year == '18' | season2 == '19-Feb')

aoa.rel.drought.shallow <- aoa.rel.drought.shallow %>% filter(season2 != '18-Feb')

aoa.rel.drought.shallow$season2 <- factor(aoa.rel.drought.shallow$season2, levels = c('18-Apr', '18-Jun', '18-Aug', '18-Oct', '18-Dec', '19-Feb'))

aoa.present <- aoa.rel.drought.shallow %>% group_by(loc, season2, Assigned_Clade_Tree2) %>% mutate(groupsum = sum(Abundance))

aoa.present <- aoa.present %>% filter(groupsum != 0)

ggplot(aoa.summary, aes(x = season2, y = Assigned_Clade_Tree2)) +
  geom_point(aes(size = mean)) + 
  #geom_point(aes(size = (mean + se)), shape = 1) + 
  facet_wrap(~loc) + 
  xlab('') +
  ylab('') + 
  theme_classic() + 
  theme(axis.text = element_text(size = 13, face = 'bold'),
        title = element_text(size = 20, face = 'bold'),
        legend.text = element_text(size = 13),
        legend.title = element_blank(),
        panel.border = element_rect(linetype = 'solid', color = 'black', fill = NA, size = 2),
        strip.text = element_text(size = 13, face = 'bold'),
        strip.background = element_rect(linewidth = 2))


aoa.summary <- aoa.present %>% group_by(loc, season2, Assigned_Clade_Tree2) %>% summarise(mean = mean(Abundance), n = n(), sd = sd(Abundance), se = sd/sqrt(n))

ggplot(aoa.summary, aes(x = season2,y = mean, col = Assigned_Clade_Tree2, group = Assigned_Clade_Tree2)) +
  geom_point() + 
  geom_line() + 
  facet_wrap(~loc) + 
  xlab('') +
  ylab('') + 
  theme_classic() + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean), width = 0.1) + 
  theme(axis.text = element_text(size = 13, face = 'bold'),
        title = element_text(size = 20, face = 'bold'),
        legend.text = element_text(size = 13),
        legend.title = element_blank(),
        panel.border = element_rect(linetype = 'solid', color = 'black', fill = NA, size = 2),
        strip.text = element_text(size = 13, face = 'bold'),
        strip.background = element_rect(linewidth = 2))
