library(ComplexHeatmap)
library(ggplot2)
library(ggrepel)
library(ggsignif)

source('./analysis/definitions.R')

df.Abbas = read.csv2(file = './data/expression/20230504_expression_Abbas.csv', sep = '\t')
df.BC = read.csv2(file = './data/expression/20230504_expression_BC.csv', sep = '\t')
df.Huo = read.csv2(file = './data/expression/20230504_expression_Huo.csv', sep = '\t')
df.Penter = read.csv2(file = './data/expression/20230504_expression_Penter.csv', sep = '\t')
df.vanGalen = read.csv2(file = './data/expression/20230504_expression_vanGalen.csv', sep = '\t')

df.combined = dplyr::bind_rows(df.Abbas, df.BC, df.Huo, df.Penter, df.vanGalen)
df.combined$predicted.celltype = factor(df.combined$predicted.celltype, levels = names(AML.combined.colors))

df.combined$meta = factor(df.combined$meta, levels = c('normal', 'AML', 'naive', 'transplant'))
# percentage of NK cells amongst T/NK cells
df.freq = df.combined %>% 
  group_by(sample, meta, dataset, predicted.celltype) %>% 
  filter(predicted.celltype %in% TNK.clusters) %>%
  tally()
df.freq$freq = apply(df.freq, MARGIN = 1, FUN = function(x) { 
  as.numeric(x['n']) /  
    as.numeric(colSums(df.freq[which(df.freq$sample == as.character(x['sample'])), 'n']))
})

df.freq$category = paste0(df.freq$dataset, '.', df.freq$meta)

write.table(df.freq %>% filter(predicted.celltype == 'NK') %>%
              arrange(category), file = './data/20230525_NK_cells_per_sample.csv', 
            sep = '\t', quote = F)

ggplot(df.freq[which(df.freq$predicted.celltype == 'NK'),], 
       aes(x=reorder(category, freq), y=100*freq, color=category)) + 
  geom_boxplot(outlier.colour = NA) + 
  geom_jitter(size=0.5) +
  scale_y_continuous('% NK cells') + 
  scale_color_manual(values = dataset.colors) +
  geom_signif(comparisons = list(c('Penter.transplant', 'Abbas.AML'),
                                 c('Penter.transplant', 'BC.AML'),
                                 #c('Penter.transplant', 'Penter.naive'),
                                 c('Penter.transplant', 'vanGalen.AML')), 
              step_increase = 0.1, color = 'black', test = 't.test', textsize = 3) + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        axis.title.x = element_blank())
ggsave('./figures/plots/20230507_NK_cells_freq.svg', width = 2, height = 3)

# expression of TIM3 
df.HAVCR2 = df.combined %>% 
  group_by(sample, meta, dataset, predicted.celltype) %>% 
  summarize(HAVCR2.expr = mean(as.numeric(HAVCR2)))

df.HAVCR2$category = paste0(df.HAVCR2$dataset, '.', df.HAVCR2$meta)

# export for prism
write.table(df.HAVCR2 %>% 
              filter(predicted.celltype == 'HSC') %>%
              filter(dataset == 'Abbas'), file = './data/20230525_HAVCR2_HSC.csv', 
            sep = '\t', quote = F)


ggplot(df.HAVCR2[which(df.HAVCR2$predicted.celltype == 'NK' & df.HAVCR2$dataset %in% c('Penter', 'Abbas')),], 
       aes(x=reorder(category, HAVCR2.expr), y=HAVCR2.expr, color=category)) + 
  geom_boxplot(outlier.colour = NA) + 
  geom_jitter(size=0.5) +
  scale_y_continuous('HAVCR2') + 
  scale_color_manual(values = dataset.colors) +
  geom_signif(comparisons = list(c('Penter.transplant', 'Abbas.AML'),
                                 c('Penter.transplant', 'Penter.naive'),
                                 c('Penter.transplant', 'Abbas.normal')), 
              step_increase = 0.1, color = 'black', test = 't.test', textsize = 3) + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        axis.title.x = element_blank())
ggsave('./figures/plots/20230507_HAVCR2_expr_Penter_Abbas_NK.svg', width = 1.5, height = 3)

ggplot(df.HAVCR2[which(df.HAVCR2$predicted.celltype == 'HSC' & df.HAVCR2$dataset %in% c('Penter', 'Abbas')),], 
       aes(x=reorder(category, HAVCR2.expr), y=HAVCR2.expr, color=category)) + 
  geom_boxplot(outlier.color =  NA) + 
  geom_jitter(size=0.5) +
  scale_y_continuous('HAVCR2') + 
  scale_color_manual(values = dataset.colors) +
  geom_signif(comparisons = list(c('Abbas.normal', 'Penter.transplant'),
                                 c('Abbas.normal', 'Penter.naive'),
                                 c('Abbas.normal', 'Abbas.AML')), 
              step_increase = 0.1, color = 'black', test = 't.test', textsize = 3) + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        axis.title.x = element_blank())
ggsave('./figures/plots/20230507_HAVCR2_expr_Penter_Abbas_HSC.svg', width = 1.5, height = 3)

