# Function modified from Gregory M Chen, Written By Alice Wang

FindDoublets <- function(seurat.rna, PCs = 1:50, exp_rate = 0.02, sct = FALSE, multiome=FALSE, v5=FALSE){
  # sct--do SCTransform or not
  DefaultAssay(seurat.rna) <- "RNA"
  
  ## pK identification
  sweep.res.list <- paramSweep_v3(seurat.rna, PCs = PCs, sct = sct)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pk_infer <- 0.09
  
  if (multiome == TRUE)  {
    pK_infer <- bcmvn[bcmvn$BCmetric == max(bcmvn$BCmetric), ]$pK
    pK_infer <- as.numeric(levels(pK_infer))[pK_infer]
  }
  
  
  ## Homotypic Doublet proportion Estimate
  annotations <- seurat.rna@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)          
  exp_rate <- (dim(seurat.rna)[2] / 2000) * 0.016
  nExp_poi <- round(exp_rate * length(seurat.rna$seurat_clusters))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies
  seurat.rna <- doubletFinder_v3(seurat.rna, PCs = PCs, pN = 0.25,
                                 pK = pk_infer, nExp = nExp_poi, reuse.pANN = FALSE, 
                                 sct = sct)
  
  seurat.rna <- doubletFinder_v3(seurat.rna, PCs = PCs, pN = 0.25, 
                                 pK = pk_infer, nExp = nExp_poi.adj,
                                 reuse.pANN = paste0("pANN_0.25_0.09_", nExp_poi), 
                                 sct = sct)
  doublet_var = paste0('DF.classifications_0.25_0.09_', nExp_poi.adj)
  seurat.rna[['Doublet_Singlet']] = seurat.rna[[doublet_var]]
  
  mnames = names(seurat.rna@meta.data)
  seurat.rna@meta.data[, grep(mnames, pattern = '0.25_0.09')] <- NULL
  return(seurat.rna)
}

