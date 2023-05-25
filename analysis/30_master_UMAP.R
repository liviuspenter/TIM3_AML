# plot master UMAP with all datasets integrated

library(dplyr)
library(Seurat)

source('./analysis/definitions.R')

AML.Abbas = readRDS('./data/objects/AML_Abbas.rds')
Abbas.embedding = Embeddings(AML.Abbas, reduction = 'ref.umap') %>% as.data.frame()
Abbas.embedding$dataset = 'Abbas'
Abbas.embedding$predicted.celltype = AML.Abbas$predicted.celltype
rm(AML.Abbas)

AML.BC = readRDS('./data/Beneyto_Calabuig/seurat_main_cohort.rds')
BC.embedding = Embeddings(AML.BC, reduction = 'ref.umap') %>% as.data.frame()
BC.embedding$dataset = 'BC'
BC.embedding$predicted.celltype = AML.BC$predicted.celltype
rm(AML.BC)

AML.Huo = readRDS('./data/objects/Normal_Huo.rds')
Huo.embedding = Embeddings(AML.Huo, reduction = 'ref.umap') %>% as.data.frame()
Huo.embedding$dataset = 'Huo'
Huo.embedding$predicted.celltype = AML.Huo$predicted.celltype
rm(AML.Huo)

AML.vanGalen = readRDS('./data/objects/AML_van_Galen.rds')
vanGalen.embedding = Embeddings(AML.vanGalen, reduction = 'ref.umap') %>% as.data.frame()
vanGalen.embedding$dataset = 'vanGalen'
vanGalen.embedding$predicted.celltype = AML.vanGalen$predicted.celltype
rm(AML.vanGalen)

AML.Penter = readRDS('./data/Penter/AML.10026.rds')
Penter.embedding = Embeddings(AML.Penter, reduction = 'ref.umap') %>% as.data.frame()
Penter.embedding$dataset = 'Penter'
Penter.embedding$predicted.celltype = AML.Penter$predicted.celltype
rm(AML.Penter)

embedding = bind_rows(Abbas.embedding, BC.embedding, Huo.embedding, vanGalen.embedding, Penter.embedding)

for (dataset in unique(embedding$dataset)) {
  message(dataset)
  subset = embedding %>% filter(dataset == dataset)
  p=ggplot() +
    geom_point(data=embedding, 
               aes(x=refUMAP_1, y=refUMAP_2), color = 'grey95', size=0.5) + 
    geom_point(data=subset,
               aes(x=refUMAP_1, y=refUMAP_2), color = 'firebrick', size=0.5) + 
    theme_classic() +
    NoLegend() + 
    NoAxes() 
  ggsave(paste0('./figures/UMAP/20230525_UMAP_all_', dataset, '.png'), width = 6, height = 6, dpi = 600, plot = p)
}

p=ggplot(embedding, aes(x=refUMAP_1, y=refUMAP_2, color=predicted.celltype)) +
  geom_point(size=0.5) + 
  scale_color_manual(values = AML.combined.colors) + 
  theme_classic() +
  NoLegend() + 
  NoAxes() 
ggsave('./figures/UMAP/20230525_UMAP_all.png', width = 6, height = 6, dpi = 600, plot = p)

subset = embedding %>% filter(predicted.celltype %in% c('MAIT', 'gdT')) #%>% sample_n(10000)
p=ggplot() +
  geom_point(data=embedding, 
             aes(x=refUMAP_1, y=refUMAP_2), color = 'grey95', size=0.5) + 
  geom_point(data=subset,
             aes(x=refUMAP_1, y=refUMAP_2, color=predicted.celltype), size=0.5) + 
  scale_color_manual(values = AML.combined.colors) + 
  theme_classic() +
  NoLegend() + 
  NoAxes() 
ggsave('./figures/UMAP/20230525_UMAP_all_MAIT_gdT.png', width = 6, height = 6, dpi = 600, plot = p)

subset = embedding %>% filter(predicted.celltype %in% c('HSC', 'LMPP', 'GMP', 'CD14 Mono', 'CD16 Mono', 
                                                        'pDC', 'cDC2', 'Prog_DC')) %>% sample_n(10000)
p=ggplot() +
  geom_point(data=embedding, 
             aes(x=refUMAP_1, y=refUMAP_2), color = 'grey95', size=0.5) + 
  geom_point(data=subset, 
             aes(x=refUMAP_1, y=refUMAP_2, color=predicted.celltype), size=0.5) + 
  scale_color_manual(values = AML.combined.colors) + 
  theme_classic() +
  NoLegend() + 
  NoAxes() 
ggsave('./figures/UMAP/20230525_UMAP_all_myelopoiesis.png', width = 6, height = 6, dpi = 600, plot = p)

subset = embedding %>% filter(predicted.celltype %in% c('Prog_RBC', 'Prog_Mk')) %>% sample_n(10000)
p=ggplot() +
  geom_point(data=embedding, 
             aes(x=refUMAP_1, y=refUMAP_2), color = 'grey95', size=0.5) + 
  geom_point(data=subset,
             aes(x=refUMAP_1, y=refUMAP_2, color=predicted.celltype), size=0.5) + 
  scale_color_manual(values = AML.combined.colors) + 
  theme_classic() +
  NoLegend() + 
  NoAxes() 
ggsave('./figures/UMAP/20230525_UMAP_all_RBC_Mk.png', width = 6, height = 6, dpi = 600, plot = p)

subset = embedding %>% filter(predicted.celltype %in% c('CD56 bright NK', 'NK')) %>% sample_n(10000)
p=ggplot() +
  geom_point(data=embedding, 
             aes(x=refUMAP_1, y=refUMAP_2), color = 'grey95', size=0.5) + 
  geom_point(data=subset,
             aes(x=refUMAP_1, y=refUMAP_2, color=predicted.celltype), size=0.5) + 
  scale_color_manual(values = AML.combined.colors) + 
  theme_classic() +
  NoLegend() + 
  NoAxes() 
ggsave('./figures/UMAP/20230525_UMAP_all_NK.png', width = 6, height = 6, dpi = 600, plot = p)

subset = embedding %>% filter(predicted.celltype %in% c('CD4 Naive', 'CD4 Memory', 'Treg')) %>% sample_n(10000)
p=ggplot() +
  geom_point(data=embedding, 
             aes(x=refUMAP_1, y=refUMAP_2), color = 'grey95', size=0.5) + 
  geom_point(data=subset,
             aes(x=refUMAP_1, y=refUMAP_2, color=predicted.celltype), size=0.5) + 
  scale_color_manual(values = AML.combined.colors) + 
  theme_classic() +
  NoLegend() + 
  NoAxes() 
ggsave('./figures/UMAP/20230525_UMAP_all_CD4.png', width = 6, height = 6, dpi = 600, plot = p)

subset = embedding %>% filter(predicted.celltype %in% c('CD8 Effector_1', 'CD8 Effector_2', 
                                                        'CD8 Memory_1', 'CD8 Memory_2', 'CD8 Naive')) %>% 
  sample_n(10000)
p=ggplot() +
  geom_point(data=embedding, 
             aes(x=refUMAP_1, y=refUMAP_2), color = 'grey95', size=0.5) + 
  geom_point(data=subset,
             aes(x=refUMAP_1, y=refUMAP_2, color=predicted.celltype), size=0.5) + 
  scale_color_manual(values = AML.combined.colors) + 
  theme_classic() +
  NoLegend() + 
  NoAxes() 
ggsave('./figures/UMAP/20230525_UMAP_all_CD8.png', width = 6, height = 6, dpi = 600, plot = p)


subset = embedding %>% filter(predicted.celltype %in% c('Prog_B 1', 'Prog_B 2', 'Plasmablast')) %>% 
  sample_n(10000)
p=ggplot() +
  geom_point(data=embedding, 
             aes(x=refUMAP_1, y=refUMAP_2), color = 'grey95', size=0.5) + 
  geom_point(data=subset,
             aes(x=refUMAP_1, y=refUMAP_2, color=predicted.celltype), size=0.5) + 
  scale_color_manual(values = AML.combined.colors) + 
  theme_classic() +
  NoLegend() + 
  NoAxes() 
ggsave('./figures/UMAP/20230525_UMAP_all_B_progenitor.png', width = 6, height = 6, dpi = 600, plot = p)

subset = embedding %>% filter(predicted.celltype %in% c('Naive B', 'Memory B')) %>% 
  sample_n(10000)
p=ggplot() +
  geom_point(data=embedding, 
             aes(x=refUMAP_1, y=refUMAP_2), color = 'grey95', size=0.5) + 
  geom_point(data=subset,
             aes(x=refUMAP_1, y=refUMAP_2, color=predicted.celltype), size=0.5) + 
  scale_color_manual(values = AML.combined.colors) + 
  theme_classic() +
  NoLegend() + 
  NoAxes() 
ggsave('./figures/UMAP/20230525_UMAP_all_B_cell.png', width = 6, height = 6, dpi = 600, plot = p)

