# extract expression for CEACAM1, HMGB1, LGALS9, HAVCR2 and PTDSS1 from each dataset

genes = c('CEACAM1', 'HAVCR2', 'HMGB1', 'LGALS9', 'PTDSS1')

### pull data from Seurat objects

# Abbas
AML.Abbas = readRDS('./data/objects/AML_Abbas.rds')
AML.Abbas = ScaleData(AML.Abbas, features = genes)

df.Abbas = as.data.frame(t(GetAssayData(AML.Abbas, slot = 'scale.data', 
                                         assay = 'RNA')[genes,]))
df.Abbas$predicted.celltype = AML.Abbas$predicted.celltype[rownames(df.Abbas)]
df.Abbas$sample = AML.Abbas$orig.ident[rownames(df.Abbas)]
rm(AML.Abbas)

# Beneyto-Calabuig
AML.BC = readRDS('./data/Beneyto_Calabuig/seurat_main_cohort.rds')
AML.BC = ScaleData(AML.BC, features = genes)

df.BC = as.data.frame(t(GetAssayData(AML.BC, slot = 'scale.data', 
                                        assay = 'RNA')[genes,]))
df.BC$predicted.celltype = AML.BC$predicted.celltype[rownames(df.BC)]
df.BC$sample = AML.BC$sample[rownames(df.BC)]
rm(AML.BC)

# Huo

AML.Huo = readRDS('./data/objects/Normal_Huo.rds')
AML.Huo = ScaleData(AML.Huo, features = genes)

df.Huo = as.data.frame(t(GetAssayData(AML.Huo, slot = 'scale.data', 
                                     assay = 'RNA')[genes,]))
df.Huo$predicted.celltype = AML.Huo$predicted.celltype[rownames(df.Huo)]
df.Huo$sample = AML.Huo$orig.ident[rownames(df.Huo)]
rm(AML.Huo)

# van Galen
AML.vanGalen = readRDS('./data/objects/AML_van_Galen.rds')
AML.vanGalen = ScaleData(AML.vanGalen, features = genes)

df.vanGalen = as.data.frame(t(GetAssayData(AML.vanGalen, slot = 'scale.data', 
                                           assay = 'RNA')[genes,]))
df.vanGalen$predicted.celltype = AML.vanGalen$predicted.celltype[rownames(df.vanGalen)]
df.vanGalen$sample = AML.vanGalen$orig.ident[rownames(df.vanGalen)]
rm(AML.vanGalen)

# Penter
AML.Penter = readRDS('./data/Penter/AML.10026.rds')
AML.Penter = ScaleData(AML.Penter, features = genes)

df.Penter = as.data.frame(t(GetAssayData(AML.Penter, slot = 'scale.data', 
                                           assay = 'RNA')[genes,]))
df.Penter$predicted.celltype = AML.Penter$predicted.celltype[rownames(df.Penter)]
df.Penter$sample = AML.Penter$orig.ident[rownames(df.Penter)]
rm(AML.Penter)

### write output to csv files along with basic metadata

df.Abbas = df.Abbas[which(grepl('NL|PT.A', df.Abbas$sample)),]
df.Abbas$meta = ifelse(grepl('NL', df.Abbas$sample), 'normal', 'AML')
df.Abbas$dataset = 'Abbas'
write.table(df.Abbas, file = './data/expression/20230504_expression_Abbas.csv', sep = '\t', quote = F)

df.BC = df.BC[which(grepl('d0', df.BC$sample)),]
df.BC$meta = ifelse(df.BC$sample %in% c('A.0.d0', 'Reference.d0'), 'normal', 'AML')
df.BC$dataset = 'BC'
write.table(df.BC, file = './data/expression/20230504_expression_BC.csv', sep = '\t', quote = F)

df.Huo$meta = 'normal'
df.Huo$dataset = 'Huo'
write.table(df.Huo, file = './data/expression/20230504_expression_Huo.csv', sep = '\t', quote = F)

df.Penter = df.Penter[which(grepl('\\.1', df.Penter$sample)),]
df.Penter$meta = ifelse(df.Penter$sample %in% c('AML1007.1', 'AML1010.1', 'AML1012.1', 'AML1016.1', 'AML1019.1',
                                                'AML1022.1', 'AML1026.1'), 'transplant', 'naive')
df.Penter$dataset = 'Penter'
write.table(df.Penter, file = './data/expression/20230504_expression_Penter.csv', sep = '\t', quote = F)

df.vanGalen$meta = ifelse(grepl('BM', df.vanGalen$sample), 'normal', 'AML')
df.vanGalen$dataset = 'vanGalen'
write.table(df.vanGalen, file = './data/expression/20230504_expression_vanGalen.csv', sep = '\t', quote = F)

