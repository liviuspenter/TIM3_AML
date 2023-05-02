setwd('/Users/liviuspenter/dfci/TIM3/')

AML.10026 = readRDS('./data/objects/AML.10026.rds')

AML.10026 = ScaleData(AML.10026, features = c('CEACAM1', 'LGALS9' ,'HAVCR2', 'PTDSS1'))

df = as.data.frame(t(GetAssayData(AML.10026, slot = 'scale.data', 
                                  assay = 'RNA')[c('CEACAM1', 'LGALS9' ,'HAVCR2', 'PTDSS1'),]))

df$predicted.celltype = factor(AML.10026$predicted.celltype[rownames(df)], levels = names(nanoranger.R::AML.combined.colors))
df$sample = AML.10026$orig.ident[rownames(df)]

ggplot(df %>% group_by(sample, predicted.celltype) %>% 
  summarize(HAVCR2.expr = mean(HAVCR2)), aes(x=predicted.celltype, y=HAVCR2.expr)) + 
  geom_jitter() + 
  stat_summary(geom='crossbar', linewidth=0.5, size=0.5, fun=median, fun.min = median, fun.max = median, color='firebrick') + 
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        axis.title.x = element_blank())

ggplot(df %>% group_by(sample, predicted.celltype) %>% 
         summarize(expr = mean(CEACAM1)), aes(x=predicted.celltype, y=expr)) + 
  geom_jitter() + 
  stat_summary(geom='crossbar', linewidth=0.5, size=0.5, fun=median, fun.min = median, fun.max = median, color='firebrick') + 
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        axis.title.x = element_blank())

ggplot(df %>% group_by(sample, predicted.celltype) %>% 
         summarize(expr = mean(LGALS9)), aes(x=predicted.celltype, y=expr)) + 
  geom_jitter() + 
  stat_summary(geom='crossbar', linewidth=0.5, size=0.5, fun=median, fun.min = median, fun.max = median, color='firebrick') + 
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        axis.title.x = element_blank())
