# plot UMAPs for Beneyto-Calabuig dataset

library(dplyr)
library(ggplot2)

genes = c('CEACAM1', 'HAVCR2', 'HMGB1', 'LGALS9', 'PTDSS1')

AML.BC = readRDS('./data/Beneyto_Calabuig/seurat_main_cohort.rds')
AML.BC = ScaleData(AML.BC, features = genes)

for (g in genes) {
  message(g)
  p=FeaturePlot(subset(AML.BC, sample %in% c('A.1.d0', 'A.2.d0', 'A.3.d0', 'A.4.d0', 'A.5.d0', 'A.6.d0',
                                             'A.7.d0', 'A.8.d0', 'A.9.d0', 'A.10.d0', 'A.11.d0', 'A.12.d0',
                                             'A.12.d0', 'A.13.d0', 'A.14.d0', 'A.15.d0', 'B.1.d0', 'B.2.d0', 'B.3.d0', 'B.4.d0')), 
                reduction = 'ref.umap', features = g, slot = 'scale.data', order=T) +
    scale_color_gradientn(colours = BuenColors::jdb_palette(name = 'brewer_celsius')) + 
    NoLegend() + NoAxes() + 
    theme(plot.title = element_blank())
  ggsave(paste0('./figures/UMAP/20230508_BC_',g,'.png'), width = 4, height = 4, dpi = 600, plot = p)
}

for (g in genes) {
  message(g)
  p=FeaturePlot(subset(AML.BC, sample %in% c('A.0.d0', 'Reference.d0')), 
                reduction = 'ref.umap', features = g, slot = 'scale.data', order=T) +
    scale_color_gradientn(colours = BuenColors::jdb_palette(name = 'brewer_celsius')) + 
    NoLegend() + NoAxes() + 
    theme(plot.title = element_blank())
  ggsave(paste0('./figures/UMAP/20230508_BC_normal_',g,'.png'), width = 4, height = 4, dpi = 600, plot = p)  
}
