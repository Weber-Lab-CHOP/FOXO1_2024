# Code written by Alice Wang: Alice.Wang@pennmedicine.upenn.edu Github: alice-y-wang # 
# load libraries # ----

library(data.table)
library(magrittr)
library(ggplot2)
library(Seurat)
library(RColorBrewer)
library(viridis)
library(compiler)
library(readr)
library(matrixStats)
library(Matrix)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(R.utils)
library(limma)
library(tidyverse)
library(DoubletFinder)
library(plyr)
library(rlist)
library(hdf5r)
library(readr)

# load source file for relevant single-cell functions and set working directory # -----
source('./sc-utilities.R')
setwd('./')

# 0) load and create seurat objects for all samples # ------

# set data path and sample names 
raw.data.path <- "./"
sample.names <- c("FOXO1_tumor", "NGFR_tumor")


# iterate through each sample and create and save seurat object. 
for (i in 1:length(sample.names)){
  
  #load h5 file into seurat 
  h5.path <- paste0(raw.data.path, sample.names[i],"/outs/filtered_feature_bc_matrix.h5")
  sample.h5 <- Read10X_h5(h5.path)
  
  # create seurat object
  seurat.obj <- CreateSeuratObject(counts = sample.h5, min.cells=3,  min.features=200, project = sample.names[i])
  
  #add Mito
  seurat.obj[["percent.mt"]] <- PercentageFeatureSet(object=seurat.obj, pattern="^MT-")
  
  #add ribo 
  seurat.obj <- PercentageFeatureSet(seurat.obj, "^RP[SL]", col.name = "percent.ribo")
  
  
  #rename cells
  seurat.obj <- RenameCells(seurat.obj, new.names = sub("-1", paste0("-", sample.names[i]), Cells(seurat.obj)))
  
  # save in a seurat_objects folder
  saveRDS(seurat.obj, paste0("seurat_objects/", sample.names[i], "_unfiltered.RDS"))
}



# 1) Computationally identify Doublets ####
seurat.obj.path <- "seurat_objects/"
sample.names <- c("FOXO1_tumor","NGFR_tumor")

for (i in 1: length(sample.objects)){
  seurat.obj <- readRDS(paste0(seurat.obj.path, sample.names[i], "_unfiltered.RDS"))
  
  #default parameters for now. 
  seurat.obj <- NormalizeData(seurat.obj, assay="RNA", normalization.method = "LogNormalize", scale.factor=10000) %>%
    FindVariableFeatures (assay="RNA", selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(assay = "RNA", verbose=TRUE) %>%
    RunPCA(npcs = 30, verbose = FALSE)%>%
    RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
    FindClusters(resolution = c(0.4, 0.8), verbose = FALSE)
  
  #find doublet
  seurat.obj <- FindDoublets(seurat.obj, exp_rate = 0.05)
  
  # save object 
  saveRDS(seurat.obj, paste0(seurat.obj.path, sample.names[i], "_unfiltered_wDoublets-Reductions.RDS"))
}

# 2) Filter and do Dimensionality Reductions ### ------

# filtering samples
seurat.obj.path <- "seurat_objects/"
sample.names <- c("FOXO1_tumor", "NGFR_tumor")


for (i in 1:length(sample.objects)){
  if (sample.names[i] == "NGFR_tumor"){
    print(sample.names[i])
    seurat.obj <- readRDS(sample.objects[i])
    seurat.obj.filtered <- subset(seurat.obj, subset = nFeature_RNA > 300 & 
                                  nFeature_RNA < 9000 & percent.mt < 10)
    saveRDS(seurat.obj.filtered, paste0(seurat.obj.path, sample.names[i], "_filtered_noReductions.RDS"))
}
  else {
    print(sample.names[i])
    seurat.obj <- readRDS(sample.objects[i])
    seurat.obj.filtered <- subset(seurat.obj, subset = nFeature_RNA > 300 & 
                                  nFeature_RNA < 7500 & percent.mt < 10, Doublet_Singlet == "Singlet")
  saveRDS(seurat.obj.filtered, paste0(seurat.obj.path, sample.names[i], "_filtered_noReductions.RDS"))
  }
  
}


# perform dimensionality reductions, regressing out cell cycle
sample.names <- c("FOXO1_tumor","NGFR_tumor")
cell_cycle=FALSE

for(i in  1: length(sample.objects)){
  seurat.obj <- readRDS(sample.objects[i])
  
  # normalizing, finding variable features, Use lognorm2 #
  DefaultAssay(seurat.obj) <- "RNA"
  seurat.obj <- NormalizeData(seurat.obj) %>%
    FindVariableFeatures(nfeatures = 2000)
  
  # Remove the genes that were express in less than 1% cells
  cell.express.perc = rowSums( seurat.obj@assays$RNA@counts >= 1) /ncol( seurat.obj@assays$RNA@counts)
  table(cell.express.perc > 0.01)
  
  sel.features <- setdiff(VariableFeatures(seurat.obj), names(cell.express.perc)[cell.express.perc < 0.01])
  
  VariableFeatures(seurat.obj) <-  sel.features
  
  #add cell cycle scoring genes 
  cycle3 = fread('regev_lab_cell_cycle_genes.txt', header = F)$V1 #replace with your own list 
  s.genes <- cycle3[1:43]
  g2m.genes <- cycle3[44:97]
  seurat.obj <- CellCycleScoring(object = seurat.obj, s.features = s.genes,
                                 g2m.features = g2m.genes, set.ident = TRUE)
  
  # heat shock regression
  heat_shock_gene <- fread('heat_shock_geneList.txt') #replace with own list 
  heat_shock_gene <- heat_shock_gene$`Approved symbol`
  seurat.obj<- AddModuleScore(seurat.obj, features = list(heat_shock_gene),
                              name = 'HeatShock.Score')
  
  if (cell_cycle == FALSE){
    print("false cell cycle")
    # REGRESS Cell cycle
    pc = 30
    DefaultAssay(seurat.obj) <- "RNA"
    seurat.obj<- ScaleData(seurat.obj, features = c(VariableFeatures(seurat.obj), c("CD3E", "CD4", "muCD19BBz-WPRE-sv40", "CD8A")), 
                           vars.to.regress = c("percent.mt","HeatShock.Score1", "nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score")) %>%
      RunPCA(npc = 30, verbose = F, features = VariableFeatures(seurat.obj))  %>%
      RunUMAP(dims = 1:pc, reduction = 'pca', reduction.name = paste0("regressCC_PC", pc, 'umap.rna'), reduction.key = paste0("regressCC_PC", pc, 'rnaUMAP_')) %>%
      FindNeighbors(dims = 1:pc, reduction = 'pca', graph.name = paste0("regressCC_snn", "PC", pc)) %>%
      FindClusters(res = c(0.2, 0.4, 0.6, 0.8), graph.name =  paste0("regressCC_snn", "PC", pc))
    
    saveRDS(seurat.obj, paste0(seurat.obj.path, sample.names[i], "_filtered_wReductions_regressCC.RDS"))
  }
}


# 3) SCT normalization # ----
sample.names <- c("FOXO1_tumor","NGFR_tumor")
seurat.obj.path <- "seurat_objects/"

for (i in 1:length(sample.objects)){
  
  seurat.obj <- readRDS(paste0(seurat.obj.path, sample.names[i], "_filtered_wReductions_regressCC.RDS"))
  
  #sct transform
  seurat.obj <- SCTransform(seurat.obj, assay="RNA", new.assay.name = "SCT", 
                            vars.to.regress = c("percent.mt", "nFeature_RNA", "nCount_RNA")) 
  
  saveRDS(seurat.obj, paste0(seurat.obj.path, sample.names[i], "_filtered_wSCT.RDS"))
  
}

# 4) Integrate by tissue type for tumor # -------------
ngfr.tumor.obj <- readRDS(paste0(seurat.obj.path, "NGFR_tumor", "_filtered_wSCT.RDS"))
foxo1.tumor.obj <-  readRDS(paste0(seurat.obj.path, "FOXO1_tumor", "_filtered_wSCT.RDS"))

tumor.list <- c(ngfr.tumor.obj,foxo1.tumor.obj)

tumor.features <- SelectIntegrationFeatures(object.list = tumor.list, nfeatures = 3000)
tumor.list <- PrepSCTIntegration(object.list =tumor.list, anchor.features = tumor.features)

tumor.anchors <- FindIntegrationAnchors(object.list = tumor.list, normalization.method = "SCT",
                                         anchor.features = tumor.features)

tumor.integrated.sct <- IntegrateData(anchorset = tumor.anchors, normalization.method = "SCT")

saveRDS(tumor.integrated.sct, "seurat_objects/Tumor_integrated_SCT.RDS")


