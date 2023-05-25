# plot UMAPs for Huo dataset

library(dplyr)
library(ggplot2)

genes = c('CEACAM1', 'HAVCR2', 'HMGB1', 'LGALS9', 'PTDSS1')

AML.Huo = readRDS('./data/objects/Normal_Huo.rds')
AML.Huo = ScaleData(AML.Huo, features = genes)

for (g in genes) {
  message(g)
  p=FeaturePlot(AML.Huo, reduction = 'ref.umap', features = g, slot = 'scale.data', order=T) +
    scale_color_gradientn(colours = BuenColors::jdb_palette(name = 'brewer_celsius')) + 
    NoLegend() + NoAxes() + 
    theme(plot.title = element_blank())
  ggsave(paste0('./figures/UMAP/20230508_Huo_',g,'.png'), width = 4, height = 4, dpi = 600, plot = p)
}
