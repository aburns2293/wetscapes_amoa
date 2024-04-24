library(EML)
#remotes::install_github("EDIorg/EMLassemblyline")
library(EMLassemblyline)

tables <- c("16srRNA.aoa.absolute.csv",
            "16SrRNA.aoa.clades.csv",
            "16SrRNA.aoa.relative.csv",
            "amoa.qpcr.csv",
            'aoa.asv.csv',
            'kegg.nitrogencycle.gene.abundance.csv',
            'kegg.nitrogencycle.gene.tax.csv',
            'kegg.nitrogencycle.sample.csv',
            'mRNA.aoa.absolute.csv',
            'mRNA.aob.absolute.csv',
            'nob.comammox.hits.csv',
            'precipitation.data.txt',
            'prok.asv.csv',
            'prok.sample.csv',
            'prok.tax.csv',
            'temperature.data.txt',
            'water_table_summary.csv')

tables.names <- c("AOA 16s rRNA absolute abundances",
                  "AOA 16S rRNA clade assignments",
                  "AOA 16S rRNA relative abundances",
                  "qPCR results for B and A amoA",
                  'Abundance and sample data table for AOA (clean)',
                  'Abundance of nitrogen cycling genes (KEGG)',
                  'KEGG taxonomies for nitrogen cycling genes',
                  'Sample information for KEGG outputs',
                  'Absolute mRNA copies for AOA',
                  'Absoltue mRNA copies for AOB',
                  'Potential comammox NOB',
                  'Regional precipitation data',
                  'Full ASV of prokaryotic 16S rRNA',
                  'Full prokaryote sample data',
                  'Full prokaryote taxonomy',
                  'Regional temperature data',
                  'Daily water table data')

tables.descs <- c('Absolute abundance of 16s rRNA gene copies for AOA',
                  'Clade identification for AOA ASVs',
                  'Relative abundance of 16S rRNA gene copies for AOA',
                  'qPCR results for bacterial and archaeal amoA',
                  'AOA relative abundance ASV table',
                  'KEGG-identified nitrification gene abundances from metatranscriptome',
                  'KEGG-identified gene taxonomy',
                  'KEGG-assosciated sample information',
                  'AOA absolute metatranscriptomic data',
                  'AOB absolute metatranscriptomic data',
                  'NOB BLASTn results for comammox ID',
                  'Site precipitation data',
                  'Full prokaryote ASV table for 16S rRNA metagenome',
                  'Full sample data',
                  'Full prokaryote taxonomy',
                  'Site temperature data',
                  'Water table depth below surface')

template_core_metadata(path = './Data/Metadata',
                       license = 'CCBY',
                       file.type = '.txt')

template_table_attributes(
  path = './Data',
  data.table = tables)

template_categorical_variables(path = './Data')

template_geographic_coverage(path = './Data',
                             data.path = './Data',
                             data.table = 'geo.data.csv',
                             site.col = 'site',
                             lat.col = 'lat',
                             lon.col = 'long')

make_eml(path = './Data/EML_Data',
         data.path = './Data',
         eml.path = './Data',
         dataset.title = 'Divergent responses of soil AOA and AOB to drought in rewetted fen peatlands',
         temporal.coverage = c('2017-04-01', '2020-02-01'),
         maintenance.description = 'Complete',
         data.table = tables,
         data.table.name = tables.names, 
         data.table.description = table.descs,
         return.obj = TRUE)

# styling

eml.aoa <- read_xml('./Data/metadata_formatted.xml')
style.sheet <- read_xml('./Data/eml_appinfo2documentation.xsl')
style.sheet2 <- read_xml('./Data/ascii-treeview.xsl')

out.doc <- xml_xslt(eml.aoa, style.sheet2)

library(htmltools)
save_html(xml_xslt(eml.aoa, style.sheet2), './Data/formatted.html')

