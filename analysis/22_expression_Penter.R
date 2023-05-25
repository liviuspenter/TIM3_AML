# plot UMAPs for Penter dataset

library(dplyr)
library(ggplot2)

genes = c('CEACAM1', 'HAVCR2', 'HMGB1', 'LGALS9', 'PTDSS1')

AML.Penter = readRDS('./data/Penter/AML.10026.rds')
AML.Penter = ScaleData(AML.Penter, features = genes)

for (g in genes) {
  message(g)
  p=FeaturePlot(subset(AML.Penter, orig.ident %in% c('AML1007.1', 'AML1010.1', 'AML1012.1', 'AML1016.1', 
                                                     'AML1019.1', 'AML1022.1', 'AML1026.2')), 
                reduction = 'ref.umap', features = g, slot = 'scale.data', order=T) +
    scale_color_gradientn(colours = BuenColors::jdb_palette(name = 'brewer_celsius')) + 
    NoLegend() + NoAxes() + 
    theme(plot.title = element_blank())
  ggsave(paste0('./figures/UMAP/20230508_Penter_',g,'_transplant.png'), width = 4, height = 4, dpi = 600, plot = p)
}

for (g in genes) {
  message(g)
  p=FeaturePlot(subset(AML.Penter, orig.ident %in% c('AML1002.1', 'AML1005.1', 'AML1008.1', 'AML3001.1', 
                                                     'AML3003.1', 'AML3005.1', 'AML5003.1', 'AML8002.1', 
                                                     'AML8004.1', 'AML8007.1')), 
                reduction = 'ref.umap', features = g, slot = 'scale.data', order=T) +
    scale_color_gradientn(colours = BuenColors::jdb_palette(name = 'brewer_celsius')) + 
    NoLegend() + NoAxes() + 
    theme(plot.title = element_blank())
  ggsave(paste0('./figures/UMAP/20230508_Penter_',g,'_naive.png'), width = 4, height = 4, dpi = 600, plot = p)  
}
