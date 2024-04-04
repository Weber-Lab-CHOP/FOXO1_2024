# Code written by Kat Mueller: muellerkp@chop.edu Github: muellerkp

# scRNA-seq analysis
# This script follows the scRNA-seq preprocessing script and contains analyses for data in Figure 4 and Extended Data Figure 9

# Script Contents:
# 1. Clean up object
#       Subset to remove cluster 11 (not T cells)
#       Reorder clusters numerically and recolor UMAP
# 2. Density plots
# 3. Add metadata for sample/cell type
#       Sample + Cluster ID
#       CD4/CD8 gates
# 4. Subset object by CD4/CD8 cell type
# 5. Add Cell Cycle scores to clusters at resolution 0.6
# 6. Add T Effector, Tissue Resident Memory, T regulatory and FOXO1-associated gene signatures
#       FOXO1 gene signature is derived from bulk RNA-seq shown in Fig. 2 using CD8+ HA CAR T cells with FOXO1 vs tNGFR. 
#       The top 100 differentially expressed genes from this comparison were used to define the signature.
# 7. FindAllMarkers and FindMarkers to identify DE genes
# 8. Gene Ontology analyses
# 9. KEGG pathway analyses
# 10. Wilcoxon Rank Sum test comparing gene signatures in FOXO1 vs. tNGFR cells

#Load Libraries-----
library(celldex)
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(rmarkdown)
library(ProjecTILs)
library(STACAS)
library(patchwork)
library(ggplot2)
library(reshape2)
library(scRNAseq)
library(SingleR)
library(EnsDb.Hsapiens.v86)
library(cetcolor)
library(Nebulosa)
library(RColorBrewer)
library(magrittr)
library(cowplot)
library(viridis)
library(devtools)
library(wesanderson)
library(tidyseurat)
library(rstatix)
library(ggpubr)
library(coin)
library(scCustomize)

#Load Object-----
setwd("~/Desktop/scRNA_Alice")
TUMOR_integrated <- readRDS("~/Desktop/scRNA_Alice/seurat_objects/Tumor_integrated_SCT.RDS")

#Clean Up Object-----
#Subset object to remove cluster 11, non-T cells
Idents(TUMOR_integrated) <- "int.snn_res.0.6"
TUMOR_integrated<- subset(TUMOR_integrated, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10))

#Reorder clusters in the UMAP
TUMOR_integrated@meta.data$int.snn_res.0.6 <- factor(TUMOR_integrated@meta.data$int.snn_res.0.6, 
                                                     levels = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
#Plot UMAP
DimPlot(TUMOR_integrated, reduction = "umap.int",
        group.by = "int.snn_res.0.6", split.by = "orig.ident") + theme(aspect.ratio = 1)

colors <- c('#332288', '#661100','#117733', '#6699CC','#44AA99', '#882255','#999933', '#CC6677','#AA4499', '#88CCEE','#DDCC77')
#Color palette chosen because it is colorblind friendly - https://stackoverflow.com/questions/9563711/r-color-palettes-for-many-data-classes
#Note that in figure 4, cluster numbering starts at 1 for ease of visualization (i.e. C0 in script is C1 in figure, etc.)

#Plot re-colored UMAP
DimPlot(TUMOR_integrated, reduction = "umap.int",
        group.by = "int.snn_res.0.6", cols = colors) + theme(aspect.ratio = 1)

#Visualize sample distribution across clusters
DimPlot(TUMOR_integrated, reduction = "umap.int",
        group.by = "int.snn_res.0.6", split.by = "orig.ident", cols = colors) + theme(aspect.ratio = 1)

#Count cells per cluster
row_count <- as.data.frame(x = table(TUMOR_integrated$int.snn_res.0.6))
row_count

#Density plots to visualize sample distribution -----
#Set embeddings for UMAP coordinates
data <- Embeddings(object = TUMOR_integrated[["umap.int"]])[(colnames(TUMOR_integrated)), c(1, 2)]
data <- as.data.frame(data)
data$cluster <- TUMOR_integrated@active.ident

##Note UMAP colors are set to WHITE so indiviual cells are not visible in density plots created in ggplot2
colors <- c('#FFFFFF', '#FFFFFF','#FFFFFF', '#FFFFFF','#FFFFFF', '#FFFFFF','#FFFFFF', '#FFFFFF','#FFFFFF', '#FFFFFF','#FFFFFF')

pl1 <- DimPlot(TUMOR_integrated, split.by = "orig.ident", cols = colors,
               combine = FALSE)

# Plots with only contours, no shading
pl1[[1]] + geom_density_2d(aes_string(x = "IntUMAP_1", y = "IntUMAP_2"), 
                           colour = "black") + xlim(-9, 9) + ylim(-7, 7) + theme(aspect.ratio = 1)

##Filled viridis version, shown in Fig. 4
pl1[[1]] + geom_density_2d_filled(aes_string(x = "IntUMAP_1", y = "IntUMAP_2"), 
                                  colour = "black") + xlim(-9, 9) + ylim(-7, 7) + theme(aspect.ratio = 1)

#Add metadata for cell types-----

#Add orig.ident to res.0.6 clusters
head(TUMOR_integrated)
Idents(TUMOR_integrated) <- "orig.ident"
TUMOR_integrated$orig.ident_0.6 <- paste(Idents(TUMOR_integrated), TUMOR_integrated$int.snn_res.0.6, sep = "_")

# Count cells from each sample in all clusters
row_count <- as.data.frame(x = table(TUMOR_integrated$orig.ident_0.6))
row_count
write.csv(row_count, "Res0.6_bysample.csv")

#Add CD4/CD8 Double negative (DN), Single Positive (SP) and Double Positive (DP) gates to orig.ident
CD8s <- WhichCells(TUMOR_integrated, expression = CD8A>0)
TUMOR_integrated$CD8s<- ifelse(colnames(TUMOR_integrated) %in% CD8s, "CD8+", "CD8-")

CD4s <- WhichCells(TUMOR_integrated, expression = CD4>0)
TUMOR_integrated$CD4s<- ifelse(colnames(TUMOR_integrated) %in% CD4s, "CD4+", "CD4-")

Idents(TUMOR_integrated) <- "CD8s"
TUMOR_integrated$CD4s.CD8s <- paste(Idents(TUMOR_integrated), TUMOR_integrated$CD4s, sep = "_")

row_count <- as.data.frame(x = table(TUMOR_integrated$CD4s.CD8s))
row_count

Idents(TUMOR_integrated) <- "orig.ident"
TUMOR_integrated$orig.ident_CD4s.CD8s <- paste(Idents(TUMOR_integrated), TUMOR_integrated$CD4s.CD8s, sep = "_")

row_count <- as.data.frame(x = table(TUMOR_integrated$orig.ident_CD4s.CD8s))
write.csv(row_count, "CD4s.CD8s_bysample.csv")

Idents(TUMOR_integrated) <- "CD4s.CD8s"
TUMOR_integrated$CD4s.CD8s_res0.6 <- paste(Idents(TUMOR_integrated), TUMOR_integrated$int.snn_res.0.6, sep = "_")

row_count <- as.data.frame(x = table(TUMOR_integrated$CD4s.CD8s_res0.6))
write.csv(row_count, "CD4s.CD8s_res0.6.csv")

#Subset object for CD4 or CD8 cells-----
Idents(TUMOR_integrated) <- "CD4s.CD8s"
CD4.SP<- subset(TUMOR_integrated, idents = "CD8-_CD4+")
DimPlot(CD4.SP, reduction = "umap.int",
        group.by = "int.snn_res.0.6", split.by = "orig.ident") + theme(aspect.ratio = 1)
row_count <- as.data.frame(x = table(CD4.SP$int.snn_res.0.6))
row_count

CD8.SP<- subset(TUMOR_integrated, idents = "CD8+_CD4-")
DimPlot(CD8.SP, reduction = "umap.int",
        group.by = "int.snn_res.0.6", split.by = "orig.ident") + theme(aspect.ratio = 1)
row_count <- as.data.frame(x = table(CD8.SP$int.snn_res.0.6))
row_count

#Add Cell Cycle score to each cluster-----
Idents(TUMOR_integrated) <- "int.snn_res.0.6"
TUMOR_integrated$Phase_res0.6 <- paste(Idents(TUMOR_integrated), TUMOR_integrated$Phase, sep = "_")

row_count <- as.data.frame(x = table(TUMOR_integrated$Phase_res0.6))
write.csv(row_count, "Phase_res0.6.csv")

#Add T effector (Teff), and Tissue Resident Memory (TRM) signatures-----
DefaultAssay(TUMOR_integrated) <- "RNA"

#Teff signature
Teff.sig <- list(c("IL2RA", "TNFRSF8", "TNFRSF4", "ICOS", "KLRG1", "HAVCR2", "TBX21", "IL2", "GZMA", 
                   "TNF", "CCL3", "CCL4", "CCL5", "PRF1", "GZMB", "IFNG", "CD69"))
TUMOR_integrated <- AddModuleScore(object = TUMOR_integrated, features = Teff.sig, name = "Teff.sig")

FeaturePlot(object = TUMOR_integrated, features = "Teff.sig1")
#module score range: -0.5 to 1.5
#Show -0.5 to 1.3
FeaturePlot(object = TUMOR_integrated, features = "Teff.sig1", split.by = "orig.ident") + theme(aspect.ratio = 1) & scale_color_gradientn(colours = rev(brewer.pal(n = 7, "RdBu")), limits = c(-0.5, 1.3)) & theme(legend.position = c(0.9, 1.1))

colors <- c('#021EA6', '#000000')
VlnPlot(TUMOR_integrated, features = c("Teff.sig1"), group.by = "orig.ident", ncol=2, pt.size=0, cols = colors) + theme(aspect.ratio = 1)

VlnPlot(TUMOR_integrated, features = c("IL2RA", "TNFRSF8", "TNFRSF4", "ICOS", "KLRG1", "HAVCR2", "TBX21", "IL2", "GZMA", 
                                       "TNF", "CCL3", "CCL4", "CCL5", "PRF1", "GZMB", "IFNG", "CD69"), group.by = "orig.ident", stack = TRUE) + theme(aspect.ratio = 1) + NoLegend() + coord_flip()

VlnPlot(TUMOR_integrated, features = c("IL2RA", "TNFRSF8", "TNFRSF4", "ICOS", "KLRG1", "HAVCR2", "TBX21", "IL2", "GZMA", 
                                       "TNF", "CCL3", "CCL4", "CCL5", "PRF1", "GZMB", "IFNG", "CD69"), group.by = "orig.ident", stack = TRUE) + theme(aspect.ratio = 1) + NoLegend()
#TRM signature
TRM.sig <- list(c("CD69", "CD101", "PDCD1", "CXCR6", "CCR8", "ITGA1", "ITGAE", "SELPLG", "ZNF683", "PRDM1", "RUNX3", "NOTCH1", "RBPJ"))
#Note on gene names:
#CD103 is ITGAE
#CD49A is ITGA1
#CLA is SELPLG
#HOBIT is ZNF683
#BLIMP1 is PRDM1
#NOTCH is NOTCH1
TUMOR_integrated <- AddModuleScore(object = TUMOR_integrated, features = TRM.sig, name = "TRM.sig")

FeaturePlot(object = TUMOR_integrated, features = "TRM.sig1")
#module score range: -0.5 to 0.5
#Showing -0.5 to 0.4

FeaturePlot(object = TUMOR_integrated, features = "TRM.sig1", split.by = "orig.ident") + theme(aspect.ratio = 1) & scale_color_gradientn(colours = rev(brewer.pal(n = 11, "PRGn")), limits = c(-0.5, 0.4)) & theme(legend.position = c(0.9, 0.9))

colors <- c('#021EA6', '#000000')
VlnPlot(TUMOR_integrated, features = c("TRM.sig1"), group.by = "orig.ident", ncol=2, pt.size=0, cols = colors) + theme(aspect.ratio = 1)

VlnPlot(TUMOR_integrated, features = c("CD69", "CD101", "PDCD1", "CXCR6", "CCR8", "ITGA1", "ITGAE", "SELPLG", "ZNF683", "PRDM1", "RUNX3", "NOTCH1", "RBPJ"), group.by = "orig.ident", stack = TRUE) + theme(aspect.ratio = 1) + NoLegend() + coord_flip()

VlnPlot(TUMOR_integrated, features = c("CD69", "CD101", "PDCD1", "CXCR6", "CCR8", "ITGA1", "ITGAE", "SELPLG", "ZNF683", "PRDM1", "RUNX3", "NOTCH1", "RBPJ"), group.by = "orig.ident", stack = TRUE) + theme(aspect.ratio = 1) + NoLegend()


#Add Treg signature-----
Treg.sig <- list(c("CTLA4", "FOXP3", "IL2RA"))
TUMOR_integrated <- AddModuleScore(object = TUMOR_integrated, features = Treg.sig, name = "Treg.sig")

VlnPlot(TUMOR_integrated, features = c("Treg.sig1"), group.by = "orig.ident", ncol=2, pt.size=0) + theme(aspect.ratio = 1)

VlnPlot(TUMOR_integrated, features = c("CTLA4", "FOXP3", "IL2RA"), group.by = "orig.ident", stack = TRUE) + theme(aspect.ratio = 1) + NoLegend() + coord_flip()

VlnPlot(TUMOR_integrated, features = c("CTLA4", "FOXP3", "IL2RA"), group.by = "orig.ident", stack = TRUE) + theme(aspect.ratio = 1) + NoLegend()


#Add FOXO1 bulk CD8 HA signatures-----
#Add FOXO1 CD8 signature from HA CAR T cells - all genes
FOXO1_HA.CD8sig<- read.table("~/Desktop/scRNA_Alice/FOXO1.genelist.csv")
head(FOXO1_HA.CD8sig)
FOXO1.list <- list(FOXO1_HA.CD8sig$V1)
FOXO1.list
TUMOR_integrated <- AddModuleScore(object = TUMOR_integrated, features = FOXO1.list, name = "FOXO1_HA.CD8sig")
head(TUMOR_integrated)
FeaturePlot(object = TUMOR_integrated, features = "FOXO1_HA.CD8sig1")
#3200 genes in this list 
FeaturePlot(object = TUMOR_integrated, features = "FOXO1_HA.CD8sig1", split.by = "orig.ident") + theme(aspect.ratio = 1) & scale_color_gradientn(colours = rev(brewer.pal(n = 7, "RdBu")), limits = c(0, 0.2)) & theme(legend.position = c(0.9, 1.1))
colors <- c('#021EA6', '#000000')
VlnPlot(TUMOR_integrated, features = c("FOXO1_HA.CD8sig1"), group.by = "orig.ident", ncol=2, pt.size=0, cols = colors) + theme(aspect.ratio = 1)
VlnPlot(TUMOR_integrated, features = c("FOXO1_HA.CD8sig1"), group.by = "int.snn_res.0.6", ncol=2, pt.size=0) + theme(aspect.ratio = 1)

#FOXO1 CD8 signature from HA CAR T cells - top 100 genes
FOXO1_HA.CD8_top100<- read.table("~/Desktop/scRNA_Alice/FOXO1_top100genes.csv")
head(FOXO1_HA.CD8_top100)
FOXO1_HA.CD8_top100 <- list(FOXO1_HA.CD8_top100$V1)
FOXO1_HA.CD8_top100
head(TUMOR_integrated)
TUMOR_integrated <- AddModuleScore(object = TUMOR_integrated, features = FOXO1_HA.CD8_top100, name = "FOXO1_HA.CD8_top100genes")
head(TUMOR_integrated)
FeaturePlot(object = TUMOR_integrated, features = "FOXO1_HA.CD8_top100genes1")
FeaturePlot(object = TUMOR_integrated, features = "FOXO1_HA.CD8_top200genes1", split.by = "orig.ident") + theme(aspect.ratio = 1) & scale_color_gradientn(colours = rev(brewer.pal(n = 7, "RdBu")), limits = c(0, 0.3)) & theme(legend.position = c(0.8, 0.9))
colors <- c('#021EA6', '#000000')
VlnPlot(TUMOR_integrated, features = c("FOXO1_HA.CD8_top200genes1"), group.by = "orig.ident", ncol=2, pt.size=0, cols = colors) + theme(aspect.ratio = 1)
VlnPlot(TUMOR_integrated, features = c("FOXO1_HA.CD8_top200genes1"), group.by = "int.snn_res.0.6", ncol=2, pt.size=0) + theme(aspect.ratio = 1)

FeaturePlot(object = TUMOR_integrated, features = "FOXO1_HA.CD8_top100genes1", split.by = "orig.ident") + theme(aspect.ratio = 1) & scale_color_gradientn(colours = rev(brewer.pal(n = 11, "PRGn"))) & theme(legend.position = c(0.8, 1.2))
FeaturePlot(object = CD4.SP, features = "FOXO1_HA.CD8_top100genes1", split.by = "orig.ident") + theme(aspect.ratio = 1) & scale_color_gradientn(colours = rev(brewer.pal(n = 11, "PRGn")), limits = c(0, .7)) & theme(legend.position = c(0.8, 1))
FeaturePlot(object = CD8.SP, features = "FOXO1_HA.CD8_top100genes1", split.by = "orig.ident") + theme(aspect.ratio = 1) & scale_color_gradientn(colours = rev(brewer.pal(n = 11, "PRGn")), limits = c(0, .7)) & theme(legend.position = c(0.8, 1))

#Save the object with added metadata-----
saveRDS(TUMOR_integrated, "seurat_objects/TUMOR_integrated_annotated.RDS")

#Extract metadata
write.csv(TUMOR_integrated@meta.data,"TUMOR.metadata.csv")
write.csv(CD4.SP@meta.data,"TUMOR.CD4.metadata.csv")
write.csv(CD8.SP@meta.data,"TUMOR.CD8.metadata.csv")


#FindAllMarkers - Upregulated genes by cluster using integrated resolution 0.6-----
Idents(TUMOR_integrated) <- "int.snn_res.0.6"
Res0.6.markers <- FindAllMarkers(TUMOR_integrated, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25  )
write.csv(Res0.6.markers, "TUMOR_integrated.Res0.6genes.csv")
markers <-Res0.6.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
print(markers, n=110)
write.csv(markers,"Res0.6_top10genes.csv")

#FindMarkers (This compares an individual cluster to all other clusters, wherease FindAllMarkers compares all clusters simultaneously)-----
#Analysis for Cluster 0 is showed in Extended Data Fig. 9. All DE analyses were used to inform cluster annotations in Fig. 4.
head(TUMOR_integrated)
Idents(TUMOR_integrated) <- "int.snn_res.0.6"
cluster0.markers <- FindMarkers(TUMOR_integrated, ident.1 = 0)
head(cluster0.markers, n=10)
#write.csv(cluster0.markers,"cluster0.markers.csv")

head(TUMOR_integrated)
Idents(TUMOR_integrated) <- "int.snn_res.0.6"
cluster1.markers <- FindMarkers(TUMOR_integrated, ident.1 = 1)
head(cluster1.markers, n=10)
#write.csv(cluster1.markers,"cluster1.markers.csv")

head(TUMOR_integrated)
Idents(TUMOR_integrated) <- "int.snn_res.0.6"
cluster2.markers <- FindMarkers(TUMOR_integrated, ident.1 = 2)
head(cluster2.markers, n=10)
#write.csv(cluster2.markers,"cluster2.markers.csv")

head(TUMOR_integrated)
Idents(TUMOR_integrated) <- "int.snn_res.0.6"
cluster3.markers <- FindMarkers(TUMOR_integrated, ident.1 = 3)
head(cluster3.markers, n=10)
#write.csv(cluster3.markers,"cluster3.markers.csv")

head(TUMOR_integrated)
Idents(TUMOR_integrated) <- "int.snn_res.0.6"
cluster4.markers <- FindMarkers(TUMOR_integrated, ident.1 = 4)
head(cluster4.markers, n=10)
#write.csv(cluster4.markers,"cluster4.markers.csv")

head(TUMOR_integrated)
Idents(TUMOR_integrated) <- "int.snn_res.0.6"
cluster5.markers <- FindMarkers(TUMOR_integrated, ident.1 = 5)
head(cluster5.markers, n=10)
#write.csv(cluster5.markers,"cluster5.markers.csv")

head(TUMOR_integrated)
Idents(TUMOR_integrated) <- "int.snn_res.0.6"
cluster6.markers <- FindMarkers(TUMOR_integrated, ident.1 = 6)
head(cluster6.markers, n=10)
#write.csv(cluster6.markers,"cluster6.markers.csv")

head(TUMOR_integrated)
Idents(TUMOR_integrated) <- "int.snn_res.0.6"
cluster7.markers <- FindMarkers(TUMOR_integrated, ident.1 = 7)
head(cluster7.markers, n=10)
#write.csv(cluster7.markers,"cluster7.markers.csv")

head(TUMOR_integrated)
Idents(TUMOR_integrated) <- "int.snn_res.0.6"
cluster8.markers <- FindMarkers(TUMOR_integrated, ident.1 = 8)
head(cluster8.markers, n=10)
#write.csv(cluster8.markers,"cluster8.markers.csv")

head(TUMOR_integrated)
Idents(TUMOR_integrated) <- "int.snn_res.0.6"
cluster9.markers <- FindMarkers(TUMOR_integrated, ident.1 = 9)
head(cluster9.markers, n=10)
#write.csv(cluster9.markers,"cluster9.markers.csv")

head(TUMOR_integrated)
Idents(TUMOR_integrated) <- "int.snn_res.0.6"
cluster10.markers <- FindMarkers(TUMOR_integrated, ident.1 = 10)
head(cluster10.markers, n=10)
#write.csv(cluster10.markers,"cluster10.markers.csv")

#Extract metadata
write.csv(TUMOR_integrated@meta.data,"TUMOR.metadata.csv")
write.csv(TUMOR_integrated@meta.data$orig.ident,"TUMOR.metadata.short.csv")
write.csv(CD4.SP@meta.data, "CD4.metadata.csv")
write.csv(CD8.SP@meta.data, "CD8.metadata.csv")


#Gene Ontology analyses-----
#Change DE genes to each set of markers to look at each cluster; Cluster 0 is shown here.
DEgenes <- cluster0.markers

#GO without filtering, by p_adj
DEgenes <- tibble::rownames_to_column(DEgenes)
geneList <- DEgenes$rowname
convertedList <- AnnotationDbi::select(org.Hs.eg.db,
                                       keys = geneList,
                                       columns = c("ENTREZID", "SYMBOL"),
                                       keytype = "SYMBOL")
geneList <- convertedList$ENTREZID
geneList <- geneList[!is.na(geneList)]

ego <- enrichGO(gene = geneList,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05,
                readable = TRUE)
head(ego)
result_table <- ego@result
result_table[1:10, 1:3]
dotplot(ego)


#Kegg pathway analyses-----
DEgenes <- cluster0.markers
head(DEgenes)
DEgenes <- tibble::rownames_to_column(DEgenes)
geneList <- DEgenes$rowname
convertedList <- AnnotationDbi::select(org.Hs.eg.db,
                                       keys = geneList,
                                       columns = c("ENTREZID", "SYMBOL"),
                                       keytype = "SYMBOL")
geneList <- convertedList$ENTREZID
geneList <- geneList[!is.na(geneList)]

#kegg pathway overrepresentation analysis
kk <- enrichKEGG(gene = geneList,
                 organism = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)
result_table <- kk@result
result_table[1:10, 1:3]
dotplot(kk)

#Wilcoxon rank sum test and effect size-----
Tumor.metadata<- read.csv("Tumor.sigs.csv")
Tumor.metadata
head(Tumor.metadata)

bxp <- ggboxplot(
  Tumor.metadata, x = "orig.ident", y = "TRM.sig1", 
  ylab = "TRM signature", xlab = "Sample", add = "jitter"
)
bxp

#wilcoxon rank sum test (This is the same as a Mann Whitney U test)
#Note both coin and rstatix packages have a wilcox_test function; it is necessary to specify which to use (rstatix).
#Bonferroni, Holm, etc. are only necessary for multiple comparisons
#Teff statistical significance
stat.test.Teff <- Tumor.metadata %>% 
  rstatix::wilcox_test(Teff.sig1 ~ orig.ident) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test.Teff

Tumor.metadata %>% wilcox_effsize(Teff.sig1 ~ orig.ident)

#TRM statistical significance
stat.test.TRM <- Tumor.metadata %>% 
  rstatix::wilcox_test(TRM.sig1 ~ orig.ident) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test.TRM

Tumor.metadata %>% wilcox_effsize(TRM.sig1 ~ orig.ident)

#Treg statistical significance
stat.test.Treg <- Tumor.metadata %>% 
  rstatix::wilcox_test(Treg.sig1 ~ orig.ident) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test.Treg

Tumor.metadata %>% wilcox_effsize(Treg.sig1 ~ orig.ident)

#Bulk RNA-seq FOXO1 gene signature (all genes): statistical significance
stat.test.FOXO1sig <- Tumor.metadata %>% 
  rstatix::wilcox_test(FOXO1_HA.CD8sig1 ~ orig.ident) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test.FOXO1sig

Tumor.metadata %>% wilcox_effsize(FOXO1_HA.CD8sig1 ~ orig.ident)

#Bulk RNA-seq FOXO1 gene signature (top 100 genes): statistical significance
FOXO1_HA.CD8_top100genes1.Score <- Tumor.metadata %>% 
  rstatix::wilcox_test(FOXO1_HA.CD8_top100genes1 ~ orig.ident) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
FOXO1_HA.CD8_top100genes1.Score

Tumor.metadata %>% wilcox_effsize(FOXO1_HA.CD8_top100genes1 ~ orig.ident)
