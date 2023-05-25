# plot UMAPs for van Galen dataset

library(dplyr)
library(ggplot2)

genes = c('CEACAM1', 'HAVCR2', 'HMGB1', 'LGALS9', 'PTDSS1')

AML.vanGalen = readRDS('./data/objects/AML_van_Galen.rds')
AML.vanGalen = ScaleData(AML.vanGalen, features = genes)

for (g in genes) {
  message(g)
  p=FeaturePlot(subset(AML.vanGalen, orig.ident %in% unique(AML.vanGalen$orig.ident)[grepl('D0', unique(AML.vanGalen$orig.ident))]), 
                reduction = 'ref.umap', features = g, slot = 'scale.data', order=T) +
    scale_color_gradientn(colours = BuenColors::jdb_palette(name = 'brewer_celsius')) + 
    NoLegend() + NoAxes() + 
    theme(plot.title = element_blank())
  ggsave(paste0('./figures/UMAP/20230508_vanGalen_',g,'.png'), width = 4, height = 4, dpi = 600, plot = p)
}

for (g in genes) {
  message(g)
  p=FeaturePlot(subset(AML.vanGalen, orig.ident %in% c('BM1', 'BM2', 'BM3', 'BM4')), 
                reduction = 'ref.umap', features = g, slot = 'scale.data', order=T) +
    scale_color_gradientn(colours = BuenColors::jdb_palette(name = 'brewer_celsius')) + 
    NoLegend() + NoAxes() + 
    theme(plot.title = element_blank())
  ggsave(paste0('./figures/UMAP/20230508_vanGalen_normal_',g,'.png'), width = 4, height = 4, dpi = 600, plot = p)  
}
