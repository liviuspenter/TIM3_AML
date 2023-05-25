# process data from Huo et al., Science Immunology 2023 downloaded from GSE224714
# data needs to be downloaded into ./data/Huo/ 

library(Seurat)

seurat.data = Read10X('./data/Huo/')
seurat.data = seurat.data[,which(grepl('N._', colnames(seurat.data)))]

Huo = CreateSeuratObject(seurat.data, project = 'Huo', min.cells = 3, min.features = 200)
Huo[['percent.mt']] <- PercentageFeatureSet(Huo, pattern = '^MT-')

Huo = subset(Huo, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)

Huo = NormalizeData(Huo)
Huo = FindVariableFeatures(Huo)
Huo = ScaleData(Huo)
Huo = RunPCA(Huo)
Huo = RunUMAP(Huo, dims = 1:30)
Huo = FindNeighbors(Huo)
Huo = FindClusters(Huo, resolution = 0.3)

# map to bone marrow reference - can be downloaded from Seurat 
# https://satijalab.org/seurat/articles/multimodal_reference_mapping.html
bm = readRDS(file = './data/objects/bm.reference.rds')
bm[["spca.annoy.neighbors"]] <- LoadAnnoyIndex(object = bm[["spca.annoy.neighbors"]], file = "./data/objects/reftmp.idx")

anchors <- FindTransferAnchors(reference = bm, query = Huo, k.filter = NA, reference.reduction = "spca", 
                               reference.neighbors = "spca.annoy.neighbors", dims = 1:50)
Huo.mapped <- MapQuery(anchorset = anchors, Huo, reference = bm, 
                refdata = list(celltype = "celltype.l2", predicted_ADT = "ADT"),
                reference.reduction = "spca", reduction.model = "wnn.umap")

saveRDS('./data/objects/Normal_Huo.rds', object = Huo.mapped)
