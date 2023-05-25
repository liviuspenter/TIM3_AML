# map dataset from Beneyto-Calabuig et al., Cell Stem Cell 2023 to Seurat reference map
# needs to be downloaded from https://doi.org/10.6084/m9.figshare.20291628.v2

library(Seurat)

AML.BC = readRDS('./data/Beneyto_Calabuig/seurat_main_cohort.rds')

# map to bone marrow reference - can be downloaded from Seurat 
# https://satijalab.org/seurat/articles/multimodal_reference_mapping.html
bm = readRDS(file = './data/objects/bm.reference.rds')
bm[["spca.annoy.neighbors"]] <- LoadAnnoyIndex(object = bm[["spca.annoy.neighbors"]], file = "./data/objects/reftmp.idx")

anchors <- FindTransferAnchors(reference = bm, query = AML.BC, k.filter = NA, reference.reduction = "spca", 
                               reference.neighbors = "spca.annoy.neighbors", dims = 1:50)
AML.BC <- MapQuery(anchorset = anchors, AML.BC, reference = bm, 
                      refdata = list(celltype = "celltype.l2", predicted_ADT = "ADT"),
                      reference.reduction = "spca", reduction.model = "wnn.umap")

AML.BC$sample = paste0(AML.BC$patient, '.', AML.BC$day)
AML.BC$orig.ident = paste0(AML.BC$patient, '.', AML.BC$day)
saveRDS(file = './data/Beneyto_Calabuig/seurat_main_cohort.rds', AML.BC)