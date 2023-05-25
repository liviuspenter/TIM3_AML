# process data from van Galen et al., Cell 2019 downloaded from NCBI Geo GSE116256
# data needs to be downloaded into ./data/van_Galen/ 

library(Seurat)

files = list.files('./data/van_Galen/', pattern = 'D0*.dem.txt.gz')
files = c(files, list.files('./data/van_Galen/', pattern = '*BM[1-4]*'))

# gather data
so.merged = NULL
for (f in files) {
  message(f)
  identifier = stringr::str_split_fixed(f, pattern = '_', n=2)[,2]
  identifier = gsub(identifier, pattern = '.dem.txt.gz', replacement = '')
  patient = stringr::str_split_fixed(identifier, pattern = '-', n=2)[,1]
  timepoint = stringr::str_split_fixed(identifier, pattern = '-', n=2)[,2]
  
  mat = as.data.frame(data.table::fread(paste0('./data/van_Galen/',f)))
  rownames(mat) = mat$Gene
  mat= mat[,-1]
  so = CreateSeuratObject(mat, project = identifier, min.cells = 3, min.features = 200)
  so$patient = patient
  so$timepoint = timepoint
  if (is.null(so.merged)) {
    so.merged = so
  } else {
    so.merged = merge(so.merged, so)
  }
}

# map to bone marrow reference - can be downloaded from Seurat 
# https://satijalab.org/seurat/articles/multimodal_reference_mapping.html
bm = readRDS(file = './data/objects/bm.reference.rds')
bm[["spca.annoy.neighbors"]] <- LoadAnnoyIndex(object = bm[["spca.annoy.neighbors"]], file = "./data/objects/reftmp.idx")

so.merged = NormalizeData(so.merged)
so.merged = FindVariableFeatures(so.merged)
so.merged = ScaleData(so.merged)
so.merged = RunPCA(so.merged)
so.merged = RunUMAP(so.merged, dims = 1:30)
so.merged = FindNeighbors(so.merged)
so.merged = FindClusters(so.merged, reso.mergedlution = 0.3)

anchors <- FindTransferAnchors(reference = bm, query = so, k.filter = NA, reference.reduction = "spca", 
                               reference.neighbors = "spca.annoy.neighbors", dims = 1:50)
AML <- MapQuery(anchorset = anchors, so, reference = bm, 
                refdata = list(celltype = "celltype.l2", predicted_ADT = "ADT"),
                reference.reduction = "spca", reduction.model = "wnn.umap")

saveRDS('./data/objects/AML_van_Galen.rds', object = AML)
