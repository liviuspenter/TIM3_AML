TNK.clusters = c('CD8 Naive','CD8 Effector_1','CD8 Effector_2','CD8 Memory_1','CD8 Memory_2',
                 'CD4 Naive','CD4 Memory','Treg','CD56 bright NK','NK','gdT','MAIT')

myeloid.clusters = c('HSC', 'LMPP', 'GMP', 'CD14 Mono', 'CD16 Mono', 'pDC', 'cDC2', 'Prog_DC')

dataset.colors = c('Huo.normal' = 'blue', 'vanGalen.normal' = 'blue', 'BC.normal' = 'blue','Abbas.normal' = 'blue',
                   'vanGalen.AML' = 'orange', 'BC.AML' = 'orange', 'Abbas.AML' = 'red', 'Penter.naive' = 'red',
                   'Penter.transplant' = 'firebrick')

AML.combined.colors = c('HSC' = RColorBrewer::brewer.pal(8, name = 'YlOrRd')[8],
                        'LMPP' = RColorBrewer::brewer.pal(8, name = 'YlOrRd')[7],
                        'GMP' = RColorBrewer::brewer.pal(8, name = 'YlOrRd')[6],
                        'CD14 Mono' = RColorBrewer::brewer.pal(8, name = 'YlOrRd')[5],
                        'CD16 Mono' = RColorBrewer::brewer.pal(8, name = 'YlOrRd')[4],
                        'pDC' = RColorBrewer::brewer.pal(8, name = 'YlOrRd')[3],
                        'cDC2' = RColorBrewer::brewer.pal(8, name = 'YlOrRd')[2],
                        'Prog_DC' = RColorBrewer::brewer.pal(8, name = 'YlOrRd')[1],
                        # CD8 T cells
                        'CD8 Naive' = RColorBrewer::brewer.pal(8, name = 'Purples')[2],
                        'CD8 Effector_1' = RColorBrewer::brewer.pal(8, name = 'Purples')[3],
                        'CD8 Effector_2' = RColorBrewer::brewer.pal(8, name = 'Purples')[4],
                        'CD8 Memory_1' = RColorBrewer::brewer.pal(8, name = 'Purples')[5],
                        'CD8 Memory_2' = RColorBrewer::brewer.pal(8, name = 'Purples')[6],
                        'CD56 bright NK' = RColorBrewer::brewer.pal(8, name = 'Purples')[7],
                        'NK' = RColorBrewer::brewer.pal(8, name = 'Purples')[8],
                        # CD4 T cells
                        'gdT' = RColorBrewer::brewer.pal(6, name = 'Blues')[2],
                        'MAIT' = RColorBrewer::brewer.pal(6, name = 'Blues')[3],
                        'CD4 Naive' = RColorBrewer::brewer.pal(6, name = 'Blues')[4],
                        'CD4 Memory' = RColorBrewer::brewer.pal(6, name = 'Blues')[5],
                        'Treg' = RColorBrewer::brewer.pal(6, name = 'Blues')[6],
                        # B cells
                        'Prog_B 1' = RColorBrewer::brewer.pal(5, name = 'Greens')[1],
                        'Prog_B 2' = RColorBrewer::brewer.pal(5, name = 'Greens')[2],
                        'Naive B' = RColorBrewer::brewer.pal(5, name = 'Greens')[3],
                        'Memory B' = RColorBrewer::brewer.pal(5, name = 'Greens')[4],
                        'Plasmablast' = RColorBrewer::brewer.pal(5, name = 'Greens')[5],
                        # other
                        'Prog_Mk' = RColorBrewer::brewer.pal(3, name = 'Greys')[2],
                        'Prog_RBC' = RColorBrewer::brewer.pal(3, name = 'Greys')[3])