# process data from Abbas et al., Nature Communications 2021 downloaded from NCBI Geo GSM5936941
# data needs to be downloaded into ./data/Abbas/ 

library(Seurat)

load('./data/Abbas/GSM5936941_readCounts_AML.rda')
boo = as.data.frame(readCounts)
colnames(boo) = gsub(colnames(boo), pattern = './GRCh38.RNA/', replacement = '')
load('./data/Abbas/GSM5936942_readCounts_TME_NL.rda')
boo2 = as.data.frame(readCounts)
colnames(boo2) = gsub(colnames(boo2), pattern = './GRCh38.RNA/', replacement = '')

AML.Abbas = CreateSeuratObject(cbind(boo, boo2))
rm(boo)
rm(boo2)

AML.Abbas[['percent.mt']] <- PercentageFeatureSet(AML.Abbas, pattern = '^MT-')

AML.Abbas = NormalizeData(AML.Abbas)
AML.Abbas = FindVariableFeatures(AML.Abbas)
AML.Abbas = ScaleData(AML.Abbas)
AML.Abbas = RunPCA(AML.Abbas)
AML.Abbas = RunUMAP(AML.Abbas, dims = 1:30)
AML.Abbas = FindNeighbors(AML.Abbas)
AML.Abbas = FindClusters(AML.Abbas, resolution = 0.3)

# map to bone marrow reference - can be downloaded from Seurat 
# https://satijalab.org/seurat/articles/multimodal_reference_mapping.html
bm = readRDS(file = './data/objects/bm.reference.rds')
bm[["spca.annoy.neighbors"]] <- LoadAnnoyIndex(object = bm[["spca.annoy.neighbors"]], file = "./data/objects/reftmp.idx")

anchors <- FindTransferAnchors(reference = bm, query = AML.Abbas, k.filter = NA, reference.reduction = "spca", 
                               reference.neighbors = "spca.annoy.neighbors", dims = 1:50)
AML.Abbas <- MapQuery(anchorset = anchors, AML.Abbas, reference = bm, 
                      refdata = list(celltype = "celltype.l2", predicted_ADT = "ADT"),
                      reference.reduction = "spca", reduction.model = "wnn.umap")

saveRDS(AML.Abbas, file = './data/objects/AML_Abbas.rds')