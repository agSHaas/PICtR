library(Seurat)
library(tidyverse)
#singlet.list <- readRDS("/omics/groups/OE0285/internal/joppsail/Doublets/Data/singlet.list.DV49.DV50")
# singlet.list <- list(
#   large=sob_concat_1.list$singlet,
#   pb=sob_concat_pbs.list$singlet,
#   inf=sob_concat_inf.list$singlet)

#singlet.list$large <- subset(singlet.list$large, feature=c("CD8a", "CD44","CD45","CD19","CD4","TCRb","CD274","CD25","CD127","CD11c"))
#singlet.list$large <- subset(singlet.list$large, downsample=40000)

#singlet.list <- lapply(singlet.list, function(list){ 
#  list.obj<-RunUMAP(list, features = rownames(list), verbose = T)
#  Idents(list.obj) <- rownames(list.obj@meta.data)
#  list.obj <- FindNeighbors(list.obj, features =rownames(list),dims = NULL,verbose = T)
#  list.obj <- FindClusters(list.obj, resolution = c(1.0), verbose = T)
  
#  return(list.obj)
#})

#saveRDS(singlet.list, "/omics/groups/OE0285/internal/joppsail/Doublets/Data/singlet.list.DV49.DV50")


doublet.list <- readRDS("/omics/groups/OE0285/internal/joppsail/Doublets/Data/doublet.list.DV49.DV50")
# doublet.list <- list(
#   large=sob_concat_1.list$doublet,
#   pb=sob_concat_pbs.list$doublet,
#   inf=sob_concat_inf.list$doublet)
doublet.list$large <- subset(doublet.list$large, feature=c("CD8a", "CD44","CD45","CD19","CD4","TCRb","CD274","CD25","CD127","CD11c"))

doublet.list <- lapply(doublet.list, function(list){ 
  list.obj<-RunUMAP(list, features = rownames(list), verbose = T)
  Idents(list.obj) <- rownames(list.obj@meta.data)
  list.obj <- FindNeighbors(list.obj, features =rownames(list),dims = NULL,verbose = T)
  list.obj <- FindClusters(list.obj, resolution = c(1.0), verbose = T)
  
  return(list.obj)
})

saveRDS(doublet.list, "/omics/groups/OE0285/internal/joppsail/Doublets/Data/doublet.list.DV49.DV50")
