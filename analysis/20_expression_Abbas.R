# plot UMAPs for Abbas dataset

library(dplyr)
library(ggplot2)

genes = c('CEACAM1', 'HAVCR2', 'HMGB1', 'LGALS9', 'PTDSS1')

AML.Abbas = readRDS('./data/objects/AML_Abbas.rds')
AML.Abbas = ScaleData(AML.Abbas, features = genes)

for (g in genes) {
  message(g)
  p=FeaturePlot(subset(AML.Abbas, orig.ident %in% c('PT1A', 'PT2A', 'PT3A', 'PT4A', 'PT5A', 
                                                    'PT6A', 'PT7A', 'PT8A')), 
                reduction = 'ref.umap', features = g, slot = 'scale.data', order=T) +
    scale_color_gradientn(colours = BuenColors::jdb_palette(name = 'brewer_celsius')) + 
    NoLegend() + NoAxes() + 
    theme(plot.title = element_blank())
  ggsave(paste0('./figures/UMAP/20230508_Abbas_',g,'.png'), width = 4, height = 4, dpi = 600, plot = p)
}

for (g in genes) {
  message(g)
  p=FeaturePlot(subset(AML.Abbas, orig.ident %in% c('NL1', 'NL2')), 
                reduction = 'ref.umap', features = g, slot = 'scale.data', order=T) +
    scale_color_gradientn(colours = BuenColors::jdb_palette(name = 'brewer_celsius')) + 
    NoLegend() + NoAxes() + 
    theme(plot.title = element_blank())
  ggsave(paste0('./figures/UMAP/20230508_Abbas_normal_',g,'.png'), width = 4, height = 4, dpi = 600, plot = p)  
}
