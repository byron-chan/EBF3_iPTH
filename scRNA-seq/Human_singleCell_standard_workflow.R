# script to perform standard workflow steps to analyze single cell RNA-Seq data

# setwd("~/Documents/Ebf3_PTHScl_OVX_scRNA_analysis")


# load libraries
library(Seurat)
library(tidyverse)
library(hdf5r)
library(DoubletFinder)
library(harmony)
library(clustree)
library(DESeq2)
library(scCustomize)
library(SingleR)
library(celldex)
library(MAST)
library(monocle3)
library(dittoSeq)
library(SeuratDisk)
library(reticulate)
library(sceasy)
library(RColorBrewer)
library(ggplot2)
library(reshape)
library(forcats)
library(dplyr)
library(zellkonverter)
library(EnhancedVolcano)
library(effsize)
library(rliger)
library(hdf5r)


# Load the dataset

veh1 = Read10X_h5(filename = '~/Documents/Human_P50_Scl_scRNA_analysis/Matrix/PTH/25/count/sample_filtered_feature_bc_matrix.h5')
veh2 = Read10X_h5(filename = '~/Documents/Human_P50_Scl_scRNA_analysis/Matrix/PTH/26/count/sample_filtered_feature_bc_matrix.h5')
pth1 = Read10X_h5(filename = '~/Documents/Human_P50_Scl_scRNA_analysis/Matrix/PTH/10/count/sample_filtered_feature_bc_matrix.h5')
pth2 = Read10X_h5(filename = '~/Documents/Human_P50_Scl_scRNA_analysis/Matrix/PTH/34/count/sample_filtered_feature_bc_matrix.h5')

# importing human atlas
Raw_data <- Read10X(data.dir = '~/Documents/Human_P50_Scl_scRNA_analysis/matrix_files/')
metadata <- read.csv('metadata.csv')
rownames(metadata) <- metadata$X
atlas <- CreateSeuratObject(counts = Raw_data, meta.data = metadata,project = "atlas", min.cells = 3, min.features = 200)
DefaultAssay(atlas) <- "RNA"

#PRESERVE ATLAS STRUCTURE

library(Seurat)
library(Matrix)
library(data.table)

# -------------------------
# Read files
# -------------------------

mat <- readMM("atlas_counts.mtx")

genes <- fread("atlas_genes.csv")$gene
cells <- fread("atlas_cells.csv")$cell

rownames(mat) <- genes
colnames(mat) <- cells

mat <- as(mat, "dgCMatrix")

meta <- fread("atlas_metadata.csv")

meta <- as.data.frame(meta)

rownames(meta) <- meta[,1]
meta <- meta[,-1]

head(rownames(meta))
head(colnames(mat))

atlas <- CreateSeuratObject(
  counts = mat,
  meta.data = meta,
  project = "Atlas"
)

atlas[["umap"]] <- CreateDimReducObject(
  embeddings = as.matrix(
    atlas@meta.data[,c("UMAP1","UMAP2")]
  ),
  key = "UMAP_"
)

table(atlas$author_celltype)

DimPlot(
  atlas,
  reduction = "umap",
  group.by = "author_celltype",
  label = TRUE
)

#### Make sure atlas has PCA + UMAP
atlas <- NormalizeData(atlas)
atlas <- FindVariableFeatures(atlas)
atlas <- ScaleData(atlas, verbose = FALSE)
atlas <- RunPCA(atlas, npcs = 30, verbose = FALSE)
atlas <- RunUMAP(atlas, dims = 1:30, reduction.name = "umap")


#filter out the atlas metadata
cols_to_keep <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "author_celltype")

atlas_filtered@meta.data <- atlas_filtered@meta.data[, cols_to_keep, drop = FALSE]




# Create one merged Seurat object, initialize the Seurat object with the raw (non-normalized data).

veh1 = CreateSeuratObject(counts = veh1, project = "vehicle1", min.cells = 3, min.features = 200)
veh2 = CreateSeuratObject(counts = veh2, project = "vehicle2", min.cells = 3, min.features = 200)
pth1 = CreateSeuratObject(counts = pth1, project = "pth1", min.cells = 3, min.features = 200)
pth2 = CreateSeuratObject(counts = pth2, project = "pth2", min.cells = 3, min.features = 200)



# 1. QC -------
# % MT reads

veh1$percent.mt = PercentageFeatureSet(veh1, pattern = "^MT-")
veh2$percent.mt = PercentageFeatureSet(veh2, pattern = "^MT-")
pth1$percent.mt = PercentageFeatureSet(pth1, pattern = "^MT-")
pth2$percent.mt = PercentageFeatureSet(pth2, pattern = "^MT-")

atlas$percent.mt = PercentageFeatureSet(atlas, pattern = "^MT-")




VlnPlot(veh1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(veh1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

VlnPlot(veh2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(veh2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

VlnPlot(pth1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(pth1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

VlnPlot(pth2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(pth2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

VlnPlot(atlas, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(pth2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')



# 2. Filtering ----------------- you want cells that express more than 200 genes and less than 30 percent mito genes
veh1_filtered <- subset(veh1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 &
                            percent.mt < 20)

veh2_filtered <- subset(veh2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 &
                          percent.mt < 20)

pth1_filtered <- subset(pth1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 &
                          percent.mt < 20)

pth2_filtered <- subset(pth2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 &
                          percent.mt < 20)

atlas_filtered <- subset(atlas, subset = nFeature_RNA > 200 & nFeature_RNA < 4000)


# merge only after QC and filtering
# merging all the files

merged_seurat <- merge(veh1_filtered, y = c(veh2_filtered, pth1_filtered, pth2_filtered),
                       add.cell.ids = c("veh1","veh2","pth1", "pth2"))

merged_seurat <- merge(pth1_filtered, y = c(pth2_filtered),
                       add.cell.ids = c("pth1", "pth2"))

merged_seurat$dataset <- "Teriparatide"

#merging again

merged_seurat <- merge(atlas, y = merged_seurat, add.cell.ids = c("Atlas", "Teriparatide"))

merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat, nfeatures = 3000)
merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- RunPCA(merged_seurat, npc = 15)
merged_seurat = FindNeighbors(merged_seurat, dims = 1:15)
merged_seurat = FindClusters(merged_seurat,resolution = c(0.1, 0.2))
merged_seurat <- RunUMAP(merged_seurat, dims = 1:20)


# Create sample column

merged_seurat$sample = rownames(merged_seurat@meta.data)
merged_seurat@meta.data = separate(merged_seurat@meta.data, col = 'sample', into = c('Condition', 'Barcode'), sep = '_')


# 3. Normalize data ----------
merged_seurat = NormalizeData(merged_seurat)


# 4. Identify highly variable features -------------- genes that are highly variable between cells

merged_seurat = FindVariableFeatures(merged_seurat)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(merged_seurat), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(merged_seurat)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# 5. Scaling -------------
merged_seurat = ScaleData(merged_seurat)

# 6. Perform Linear dimensionality reduction -------------- you have to do this then find neighbors and cluster and do non-linear (umap) after

merged_seurat = RunPCA(merged_seurat, npcs = 20)

# determine dimensionality of the data
ElbowPlot(merged_seurat, ndims = 50)

# 7. Clustering ------------ choose the dimension based on your elbow plot
merged_seurat = FindNeighbors(merged_seurat, dims = 1:20)


# understanding resolution (lower the number the lower the clusters) - gives you resolution for distinct clusters you can play around with this

merged_seurat = FindClusters(merged_seurat,resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))


#8 non-linear dimensionality reduction -------------- (umap) ----- visualizing the clusters

merged_seurat <- RunUMAP(merged_seurat, dims = 1:20)

Idents(merged_seurat) = "RNA_snn_res.0.2"

DimPlot(merged_seurat, reduction = "umap", label = T, label.size = 4, split.by = "orig.ident")
DimPlot(merged_seurat, reduction = "umap", label = T)

# look at whether we need batch correction or not
DimPlot(merged_seurat, group.by = "orig.ident", label = T)
DimPlot(merged_seurat, reduction = "umap", label = T, label.size = 8)|DimPlot(merged_seurat, group.by = "orig.ident", label = T)

# setting your sample names as idents so you can use the idents argument to call your sample names
Idents(merged_seurat) <- "orig.ident"
Idents(merged_seurat_harmony) <- "orig.ident"



# clustree to determine optimal clustering

clustree(merged_seurat_harmony, prefix = "RNA_snn_res.", node_colour = "sc3_stability", label.size = 4)

# running doublet finder ---- this should be ran before harmony! but you can decide if its worth doing doublet finder first


# running harmony ---------

merged_seurat_harmony = RunHarmony(merged_seurat, group.by.vars = 'orig.ident', max.iter = 20, plot_convergence = T)

# you can see now harmony as a new method of reduction, similar to PCA and UMAP
merged_seurat_harmony@reductions

# adding harmony embeddings for UMAP and clustering

merged_seurat_harmony_embed = Embeddings(merged_seurat_harmony, "harmony")
merged_seurat_harmony_embed[1:10,1:10]

# clustering using harmony batch correction
merged_seurat_harmony = merged_seurat_harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
  FindClusters(resolution = c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))

Idents(merged_seurat_harmony) = "cell_type"

# visualize the integration
DimPlot(merged_seurat_harmony, reduction = 'umap', label = T, split.by = 'orig.ident')

DimPlot(merged_seurat_harmony, reduction = 'umap', label = T, group.by = "RNA_snn_res.0.4")

# before and after batch correction
p2|p5
p1|p4
p5|p6

# visualize PCA results
print(merged_seurat_filtered[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(merged_seurat_filtered, dims = 1, cells = 500, balanced = TRUE)


# setting identity of clusters ----- default was 0.8 resolution

Idents(merged_seurat_harmony) = "RNA_snn_res.0.2"
Idents(merged_seurat_harmony)

# Finding DEGs (cluster biomarkers) ----------------- Find all markers will try to define differentially expressed markers for each clusters compared to rest
merged_seurat_harmony[["RNA"]] <- JoinLayers(merged_seurat_harmony[["RNA"]]) # join the layers this is new in seurat v5
recluster[["RNA"]] <- JoinLayers(recluster[["RNA"]])
merged_seurat[["RNA"]] <- JoinLayers(merged_seurat[["RNA"]])
atlas[["RNA"]] <- JoinLayers(atlas[["RNA"]])
CAR[["RNA"]] <- JoinLayers(CAR[["RNA"]])

stromal_combined[["RNA"]] <- JoinLayers(stromal_combined[["RNA"]])

markers = FindAllMarkers(merged_seurat_harmony, only.pos = T, min.pct = 0.5, logfc.threshold = 0.25)
markers.clusters = markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)

write.csv(markers.clusters, file = "recluster.markers.clusters.csv")



# Finding markers within specific clusters - in this case it is 4 and 8 (4 is Ebf3 cells)
markers_cluster4 = FindConservedMarkers(merged_seurat_filtered,
                                        ident.1 = 4,
                                        grouping.var = 'Condition')

markers_cluster8 = FindConservedMarkers(merged_seurat_filtered, ident.1 = 8, grouping.var = 'Condition')

head(markers_cluster4)
head(markers_cluster8)

# rename cluster 4 to CAR identify
Idents(merged_seurat_filtered)
merged_seurat_filtered = RenameIdents(merged_seurat_filtered, '4' = 'CAR cells')

DimPlot(merged_seurat_filtered, reduction = 'umap', label = T)

# Findmarkers between conditions - male vehicle vs male pth etc
merged_seurat_filtered$cluster.cond = paste0(merged_seurat_filtered$RNA_snn_res.0.1, '_', merged_seurat_filtered$Condition)
View(merged_seurat_filtered@meta.data)
Idents(merged_seurat_filtered) = merged_seurat_filtered$cluster.cond

DimPlot(merged_seurat_filtered, reduction = 'umap', label = T)

# find markers

m_pth_response = FindMarkers(merged_seurat_filtered, ident.1 = '4_mpth2m', ident.2 = '4_mvch2m')
f_pth_response = FindMarkers(merged_seurat_filtered, ident.1 = '4_fpth2m', ident.2 = '4_fvch2m')
head(m_pth_response)
head(f_pth_response)

# plots
DimPlot(merged_seurat_harmony, reduction = 'umap', group.by = 'Condition', label = T)
DimPlot(merged_seurat_harmony, reduction = 'umap', split.by = 'Condition', label = T)
DimPlot(merged_seurat_harmony, reduction = 'umap', label = T)

FeaturePlot_scCustom(merged_seurat_harmony, features = c("CXCL12", "LEPR" ,"ADIPOQ","KITLG", "EBF3", "FOXC1"))
FeaturePlot_scCustom(atlas, features = c("CXCL12", "EBF3", "LEPR", "KITLG"))

VlnPlot(merged_seurat_harmony, features = c("Ebf3","Cxcl12","Lepr", "Sp7", "tdTomato","Adipoq", "Bglap", "Col1a1", "Sox9"))

test <- subset(merged_seurat_harmony, idents  = "2m")
VlnPlot(test, features = "Xist", split.by = "orig.ident", group.by = "Condition")

# trying to split the data to plot by individual conditions
pdf(file = "conditions_merged.pdf",width = 6,height = 4)
conditions.merged_seurat_harmony = SplitObject(merged_seurat_harmony, split.by = "orig.ident")
plot.conditions.merged = lapply(X = conditions.merged_seurat_harmony, FUN = function(x) {
  DimPlot(x, reduction = "umap", label = T, label.size = 4)
})
plot.conditions.merged
dev.off()



# trying to split the data to plot by individual conditions
pdf(file = "conditions.pdf",width = 6,height = 4)
conditions.tdtomato = SplitObject(tdtomato, split.by = "orig.ident")
plot.conditions.tdtomato = lapply(X = conditions.tdtomato, FUN = function(x) {
  DimPlot(x, reduction = "umap", label = T, label.size = 4)
})
plot.test
dev.off()


# finding markers, only.pos = T makes it higher expression compared to lower

markers.clusters = markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)
write.csv(markers.clusters, file = "new.markers.clusters.csv")


view(pth.markers)

# finding subcluster in specific cluster

merged_seurat_harmony <- FindSubCluster(merged_seurat_harmony, "3", "RNA_snn", subcluster.name = "test_sub", resolution = 0.02, algorithm = 1)
DimPlot(merged_seurat_harmony, reduction = "umap", group.by = "test_sub", label = T)
DimPlot(merged_seurat_harmony, reduction = "umap", label = T)



# how many cells are in each cluster
Idents(merged_seurat_harmony) = "RNA_snn_res.0.2"
table(Idents(merged_seurat_harmony))

Idents(merged_seurat_harmony) = "time"
table(Idents(merged_seurat_harmony))

table(Idents(merged_seurat_harmony), merged_seurat_harmony$stim)


# how many cells are in each cluster per time point

Idents(merged_seurat_harmony) = "orig.ident"
cluster2d = subset(merged_seurat_harmony, idents = c("male_vehicle_2d", "female_vehicle_2d", "female_pth_2d"))
cluster2m = subset(merged_seurat_harmony, idents = c("male_vehicle_2m", "female_vehicle_2m", "male_pth_2m", "female_pth_2m"))

Idents(cluster2d) = "RNA_snn_res.0.2"
table(Idents(cluster2d))

Idents(cluster2m) = "RNA_snn_res.0.2"
table(Idents(cluster2m))

# proper plotting with feature plots

FeaturePlot_scCustom(merged_seurat_harmony, features = c("CXCL12","LEPR","EBF3","ADIPOQ"), colors_use = viridis_plasma_dark_high, split.by = "stim")

VlnPlot_scCustom(recluster, features = "tdTomato", ggplot_default_colors = T, pt.size = 0, add.noise = F)

FeaturePlot_scCustom(merged_seurat_harmony, features = c("Pax7"), colors_use = viridis_plasma_dark_high, split.by = "orig.ident")

VlnPlot_scCustom(merged_seurat_harmony, features = c("Ebf3","Cxcl12", "Lepr", "Adipoq", "Lpl", "Runx2", "Bglap", "Col1a1", "Sp7", "Spp1", "Dmp1", "Alpl"), idents = c("Osteo-CAR", "Adipo-CAR"), pt.size = 0)

Idents(merged_seurat_harmony) = "RNA_snn_res.0.2"
DimPlot_scCustom(merged_seurat_harmony, reduction = 'umap', colors_use = pal, label = F) %>%
  LabelClusters(id = "ident", fontface = "bold", repel = F)


# custom palette list for umap
color_pal_list <- c("#d741a7", "#afa2ff", "#5398be", "#dea54b", "#ef2917", "#a30b37", "#87b38d","#aef3e7", "#8ee3ef", "#72788d", "#bbb6df", "#f2b7c6", "#f2dc5d", "#f18701", "#b56b45", "#37ff8b")

# rename clusters and add it to metadata under column cluster_ids

merged_seurat_harmony <- RenameIdents(merged_seurat_harmony, '0' = "Erythroid progenitor", '1' = "EC1", '2' = "CAR1", '3' = "Myeloid", '4' = "Myeloid progenitor", '5' = "Neutrophil", '6' = "Erythroid precursor", '7' = "B", '8' = "Skeletal muscle", '9' = "Erythroblast", '10' = "Osteoblast", '11' = "CAR2", '12' = "Mast", '13' = "Skeletal progenitor", '14' = "EC2", '15' = "EC3")
merged_seurat_harmony$cluster_ids <- merged_seurat_harmony@active.ident

merged_seurat_harmony <- RenameIdents(merged_seurat_harmony, '4' = "CAR1", '5' = "CAR2", '3_0' = "CAR3", '3_1' = "Pre-osteoblast")
merged_seurat_harmony$recluster_ids <- merged_seurat_harmony@active.ident


#rename clusters by cell annotation

cell.annotation <- c("Erythroid1", "EC1", "CAR1", "Myeloid", "Myeloid progenitor", "Neutrophil", "Erythroid precursor", "B", "Skeletal muscle", "Erythroblast", "Osteoblast", "CAR2", "Mast", "Skeletal progenitor", "EC2", "EC3")
names(cell.annotation) <- levels(merged_seurat_harmony)
merged_seurat_harmony <- RenameIdents(merged_seurat_harmony, cell.annotation)
DimPlot(merged_seurat_harmony, reduction = "umap", label = T)

export_df <- merged_seurat_harmony@meta.data

write.csv(merged_seurat_harmony@meta.data, "cell_annotation_metadata.csv")

# Reclustering data to remove HSCs and low tdtomato count --- we want clusters EC1, CAR1, myeloid, myeloid progenitors, neutrophil, osteoblasts, car2, mast, skeletal mscs, ec2 and ec3

recluster <- subset(merged_seurat_harmony, idents = c("0", "1", "2", "3", "4", "5","6", "7", "8", "10"))

recluster <- RenameIdents(recluster, '0' = "Erythroid progenitor", '1' = "EC1", '2' = "CAR1", '3' = "Myeloid", '4' = "Myeloid progenitor", '5' = "Neutrophil", '6' = "Erythroid precursor", '7' = "B", '8' = "Skeletal muscle", '9' = "Erythroblast", '10' = "Osteoblast", '11' = "CAR2", '12' = "Mast", '13' = "Skeletal progenitor", '14' = "EC2", '15' = "EC3")
recluster$new_cluster_ids <- recluster@active.ident

# finding subclusters =====================================
#1 look at what clusters you want to subclusters at different resolutions
#2 proceed to use findsubcluster() to create new meta data column, then set identity to new meta data column
#3 combine certain clusters if need be by using renameidents()
#4 once you are happy with the clusters, add this ident into new meta data column

recluster2 <- FindSubCluster(recluster2, c("4"), "RNA_snn", subcluster.name = "new_clusters", resolution = 0.1, algorithm = 1)
Idents(recluster2) = "new_clusters"
DimPlot(recluster2, reduction = "umap", group.by = "new_clusters", label = T)

# merging some subclusters =======================

recluster1 <- RenameIdents(recluster1, '12' = "4")


recluster <- RenameIdents(recluster, '2_2' = "2_1")
table(Idents(recluster))

recluster$new_clusters <- recluster@active.ident
Idents(recluster) = "new_clusters"


# saving old recluster incase anything goes wrong
recluster.old <- recluster

# no need for reclustering ====================== only necessary if you remove cells - but in this case we removed whole clusters
recluster1 <- FindVariableFeatures(recluster1)
recluster1 <- ScaleData(recluster1)
recluster1 <- RunPCA(recluster1, npcs = 30)
recluster1 <- FindNeighbors(recluster1, dims = 1:30)
recluster1 <- FindClusters(recluster1, resolution = c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))
recluster1 <- RunUMAP(recluster1, dims = 1:30)

# need to rerun harmony again....
recluster2 <- RunHarmony(recluster1, group.by.vars = 'orig.ident', max.iter = 20, plot_convergence = T)

recluster2 <- recluster2 %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
  FindClusters(resolution = c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))

Idents(recluster2) = "RNA_snn_res.0.2"

DimPlot_scCustom(recluster1, reduction = 'umap', pt.size = 0.1, label = F, ggplot_default_colors = T) %>%
  LabelClusters(id = "ident", fontface = "bold", repel = F, size = 5)

recluster1 <- recluster2

Idents(recluster) = "new_clusters"
DimPlot_scCustom(recluster1, reduction = 'umap', ggplot_default_colors = T, pt.size = 0.5, label = F) %>%
  LabelClusters(id = "ident", fontface = "bold", repel = F, size = 5)

Idents(recluster) = "orig.ident"

DimPlot(recluster1, label = T, group.by = "RNA_snn_res.0.4")


FeaturePlot_scCustom(CAR, features = "Il31ra", colors_use = viridis_plasma_dark_high, split.by = "stim")

FeaturePlot_scCustom(recluster, features = c("Scx"), colors_use = viridis_plasma_dark_high)

VlnPlot_scCustom(recluster, features = "Scx", ggplot_default_colors = T, add.noise = F, pt.size = 0)

#"Cxcl12", "Lepr", "Adipoq","Ebf3", "Pth1r","Runx2","Sp7","Alpl","Bglap", "Pparg", "Plin2"
Stacked_VlnPlot(merged_seurat_harmony, features = c("CXCL12", "LEPR","EBF3","FOXC1","KITLG", "LIMCH1", "SP7", "RUNX2","SPP1","IBSP", "COL1A1","BGLAP", "PTH1R"), ggplot_default_colors = F, colors_use = pal, pt.size = 0, plot_legend = T, x_lab_rotate = TRUE, add.noise = F, raster = F)
Stacked_VlnPlot(atlas, features = c("CXCL12", "LEPR", "EBF3","FOXC1","KITLG", "COL1A1","BGLAP", "PTH1R", "ALPL"), ggplot_default_colors = T, pt.size = 0, plot_legend = T, x_lab_rotate = TRUE, add.noise = F, group.by = "author_celltype")


# number of cells per cluster at each condition

table(Idents(merged_seurat), merged_seurat$orig.ident)
write.csv(table(Idents(recluster), recluster$orig.ident), file = "cellnumbers.csv")


# number of cells per cluster at each condition

test <- table(Idents(merged_seurat_harmony), merged_seurat_harmony$orig.ident)
write.csv(test, file = "testchart.csv", row.names = T)

# merging all the downsampled clusters back to one object for stacked violin plots
Stacked_VlnPlot(merged_seurat_harmony, features = c("EBF3","CXCL12"), ggplot_default_colors = T, pt.size = 0, plot_legend = F)


#PTH receptor expressed in cluster 1, 3, 4, and 7 - look at response of PTH in these first, others will be secondary

# getting markers for reclusters
Idents(recluster1) = "RNA_snn_res.0.4"
new.markers = FindAllMarkers(recluster, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)

# getting CAR markers

CAR <- subset(recluster, idents = c("CAR1","CAR2","CAR3","CAR4","AdipoCAR", "Pre-osteoblast", "Mature osteoblast"))
CARDEG <- subset(recluster, idents = c("CAR1","CAR2","CAR3"))

Idents(CAR) = "stim"
Idents(CAR) = "cluster_ids"


markers = FindAllMarkers(CARDEG, only.pos = T, min.pct = 0.4, min.diff.pct = 0.2, logfc.threshold = 0.5, test.use = "MAST")
markers.clusters = markers %>%
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC)

write.csv(markers.clusters, file = "reclustered.markers.csv")



#sorting the markers
new.markers.clusters = new.markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC) #### can also use slice_max(n = 5, order_by = avg_log2FC)
write.csv(new.markers.clusters, file = "recluster.markers.csv")

heatmap.clusters <- read_csv("recluster.markers.heatmap.csv")

maxcells <- min(table(Idents(recluster))) ## figure out what is the minimum cells for each cluster
DoHeatmap(subset(recluster_heatmap, downsample = maxcells), features = heatmap.clusters$gene, disp.max = 2, size = 4) + guides(colour=F) + theme(axis.text.y = element_text(size = 10), text = element_text(size = 15))

recluster_heatmap <- subset(recluster, idents = c("osteoCAR1", "osteoCAR2", "adipoCAR", "osteoblast1", "osteoblast2", "immature osteocyte"))


# compare cluster markers = will do this between car and osteoblast clusters (2 and 10)

Idents(recluster) = "new_clusters"
# start with only.pos = F to see what is different between genes and then you can adjust to only see which is only positive
car.markers <- FindMarkers(recluster, ident.1 = "4", ident.2 = c("7"), min.pct = 0.5, logfc.threshold = 1, only.pos = F)
head(car.markers, n = 10)
write.csv(car.markers, file = "CARDEG10v1.csv")

cluster11.markers <- FindMarkers(recluster, ident.1 = "11", ident.2 = c("2_0", "2_1", "1", "10_0", "10_1", "10_2", "13", "14_0", "14_1", "15"), min.pct = 0.5, min.diff.pct = 0.4, only.pos = T)
head(cluster11.markers, n = 10)
write.csv(cluster11.markers, file = "DEGS2vs11.csv")

cluster10.2v10.markers <- FindMarkers(recluster, ident.1 = c("10_0", "10_1"), ident.2 = "10_2", min.pct = 0.5, min.diff.pct = 0.5, only.pos = T)
head(cluster10.2v10.markers, n = 20)
write.csv(cluster10.2v10.markers, file = "DEGS10_2vs10.csv")

cluster10.0v1.markers <- FindMarkers(recluster, ident.1 = "10_0", ident.2 = "10_1", min.pct = 0.5, min.diff.pct = 0.5, only.pos = T)
head(cluster10.0v1.markers, n = 20)
write.csv(cluster10.0v1.markers, file = "DEGS10_0vs10_1.csv")

cluster10.1v0.markers <- FindMarkers(recluster, ident.1 = "10_1", ident.2 = "10_0", min.pct = 0.5, min.diff.pct = 0.5, only.pos = T)
head(cluster10.1v0.markers, n = 20)
write.csv(cluster10.1v0.markers, file = "DEGS10_0vs10_1.csv")


# cluster re-order in ascending numeric order
levels(recluster2) <- c("0", "1", "3", "7", "4_0", "4_1", "13","2", "5", "11", "14", "6", "9", "10", "8", "12")
levels(merged_seurat_harmony2) <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
levels(recluster) <- c("OsteoCAR1", "OsteoCAR2", "OsteoCAR3", "OsteoCAR4","AdipoCAR", "Pre-osteoblast", "Mature osteoblast", "Periosteal1", "Periosteal2", "Periosteal3", "Endothelial", "Pericyte1", "Pericyte2", "Myogenic")
levels(osteo.10d) <- c("OsteoCAR1", "OsteoCAR2", "AdipoCAR", "Osteoblast1", "Osteoblast2")

# reodering clusters to be compatible with dittoplot - cluster_order
cluster_order <- match(levels(merged_seurat_harmony2@meta.data[["RNA_snn_res.0.2"]]), metaLevels("RNA_snn_res.0.2", merged_seurat_harmony2))

Idents(recluster2) = "new_clusters"

# how many cells are in each cluster and per time point
Idents(recluster) = "new_clusters"
table(Idents(recluster))

Idents(recluster) = "time"
Idents(recluster) = "stim"



# looking at DEGs within each clusters across conditions =============================================================
# current idents 1 - EPC, 2_0 - CAR1, 2_1 - CAR2, 10_0 - OB1, 10_1 - OB2, 10_2 - OB3, 11 - CAR3, 13 - MSC progenitor, 14_0 - EC1, 14_1 - EC2, 15 - EC3
Idents(recluster) = "cluster_ids"

test <- subset(recluster2, idents = "4_0")

car1 <- subset(CAR, idents = "CAR1") #pth receptor
car2 <- subset(CAR, idents = "CAR2") #pth receptor
car3 <- subset(CAR, idents = "CAR3") #some cells have pthr receptor
car4 <- subset(recluster, idents = "CAR4") #some cells have pthr receptor
car5 <- subset(recluster, idents = "AdipoCAR") #some cells have pthr receptor
ob1 <- subset(recluster, idents = "Pre-osteoblast") #pth receptor
ob2 <- subset(recluster, idents = "Mature osteoblast")
per1 <- subset(recluster, idents = "Periosteal1") #
per2 <- subset(recluster, idents = "Periosteal2") #
pc <- subset(recluster, idents = "Pericyte")
ec1 <- subset(recluster, idents = "Endothelial1")
ec2 <- subset(recluster, idents = "Endothelial2")
myo <- subset(recluster, idents = "Myogenic")

# DEGS for 2d vs 2m ---- 2m done already ---- starting 2d
Idents(car1) = "stim" #
Idents(car2) = "stim" #
Idents(car3) = "stim" #
Idents(car4) = "time" #
Idents(car5) = "time" #
Idents(ob1) = "time" #
Idents(ob2) = "time" #
Idents(per1) = "time" #
Idents(per2) = "time"

Idents(test) = "time"




FeaturePlot_scCustom(subset_ob2.2m, features = c("Bglap"), max.cutoff = 4, split.by = "stim", colors_use = viridis_plasma_dark_high)
VlnPlot_scCustom(subset_car2.2d, features = c("Lars2"), pt.size = 0)

# subsetting clusters with only 2d - osteoblasts first
test.2w <- subset(test, idents = "2w")
DimPlot(test.2w, split.by = "stim", group.by = "Condition")


table(test.2w@meta.data$Condition)
table(test.2w@meta.data$time)

#============================================================================================
stim_cells = Cells(car3)[which(car3$stim == "Teriparatide")]
ctrl_cells = Cells(car3)[which(car3$stim == "Control")]

# for 2 day analysis only --- mvch2d is actually fvch2d
female_ctrl_cells = Cells(epc.2d)[which(epc.2d$Condition == "mvch2d")]

# use later ctrl_cells = Cells(ob.2d)[which(ob.2d$stim == "Control")]
slct_stim_cells = sample(stim_cells, size = 68, replace = F)
slct_ctrl_cells = sample(ctrl_cells, size = 68, replace = F)

# use for 2 day analysis
slct_ctrl_cells = sample(female_ctrl_cells, size = 156, replace = F)

# continue for analysis
subset_car3 = subset(car3, cells = c(slct_stim_cells, slct_ctrl_cells))

table(subset_epc.2m@meta.data$Condition)
table(subset_car3@meta.data$stim)

DimPlot(subset_car2, group.by = "Condition", split.by = "stim")
DimPlot(subset_car3, group.by = "stim")

#============================================================================================

Idents(subset_car3) = "stim"

### no thresholding or mini.pct for GSEA analysis - needs background genes
car3.de <- FindMarkers(subset_car3, ident.1 = "Teriparatide", ident.2 = "Control", logfc.threshold = 0, min.pct = 0, test.use = "MAST")
head(car3.de, n=1000)
write.csv(car3.de, file = "car3PTHvsatlasDEG.csv")

### for over-representation analysis - you should have thresholding
car2.de.2m <- FindMarkers(subset_car1.2m, ident.1 = "PTH", ident.2 = "Control", logfc.threshold = 0.5 , min.pct = 0.5, only.pos = F)
head(car1.de.2m, n=10)
write.csv(car1.de.2m, file = "2mcar1degsORA.csv")


FeaturePlot_scCustom(subset_car1.2m, features = c("Gli3"), max.cutoff = "4", split.by = "stim", colors_use = viridis_plasma_dark_high)
VlnPlot(subset_car1.2m, features = c("Lrp6"), split.by = "stim", pt.size = 0)


# plotting volcanoplots with enhanced volcano
library(EnhancedVolcano)

# processing the DEG dataframe computed from FindMarkers()

files <- list.files(path = "volcanodeg", pattern = "\\.csv$", full.names = TRUE)
df <- do.call(rbind, lapply(files, read.csv))


### better to use this
df <- df %>%
  group_by(X) %>%
  #filter(abs(avg_log2FC) > 0.5, p_val_adj < 0.05) %>%   # only keep rows with |avg_log2FC| > 0.5
  slice_min(order_by = if_else(p_val_adj == 0, Inf, p_val_adj),
            n = 1, with_ties = FALSE) %>%
  ungroup()

write.csv(df, file = "combinedCAR123DEGS.csv", row.names = T)


# set the base colour as 'black'
keyvals <- rep('#000000', nrow(df))

# set the base name/label as 'NS'
names(keyvals) <- rep('NS', nrow(df))

# modify keyvals for variables with fold change > 0.5
keyvals[which(df$avg_log2FC > 0.5 &
                df$p_val_adj < 0.05)] <- 'red3'
names(keyvals)[which(df$avg_log2FC > 0.5 &
                       df$p_val_adj < 0.05)] <- 'Upregulated'

# modify keyvals for variables with fold change < -0.5
keyvals[which(df$avg_log2FC < -0.5 &
                df$p_val_adj < 0.05)] <- '#3C84FF'
names(keyvals)[which(df$avg_log2FC < -0.5 &
                       df$p_val_adj < 0.05)] <- 'Downregulated'

EnhancedVolcano(df,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                lab =df$X,
                pCutoff = 0.05,
                FCcutoff = 0.5,
                drawConnectors = T,
                lengthConnectors = unit(0.02, 'npc'),
                boxedLabels = T,
                ylab = bquote(-~Log[10]~ "(adjusted P value)"),
                xlab = bquote(~Log[2]~"fold change"),
                col = c('black', 'black', 'black', 'red3'),
                colAlpha = 0.6,
                colCustom = keyvals,
                #ylim = c(0,10),
                #xlim = c(-13,20.5),
                cutoffLineWidth = 1,
                labSize = 6,
                gridlines.minor = F, gridlines.major = F,
                pointSize = c(ifelse(abs(df$avg_log2FC) > 0.5 & df$p_val_adj < 0.05, 2, 2)),
                #selectLab = c("Col1a1", "Spp1", "Alpl", "Grem1"),
                selectLab = c("EBF3", "FOXC1", "PPARG", "APOE","CXCL12", "KITLG", "COL1A1","PDE4D"),
                labFace = "bold")



#------- SingleR automated single cell annotation
# will probably want to use ImmGenData (mouse bulk RNAseq) and MouseRNAseqData

ref.data1 <- celldex::ImmGenData()
ref.data2 <- celldex::MouseRNAseqData()
seurat.counts <- GetAssayData(merged_seurat_harmony, assay = 'RNA', layer = 'counts')

singleR.annotate <- SingleR(seurat.counts, ref = ref.data2, labels = ref.data2$label.main)
merged_seurat_harmony$singleR.labels2 <- singleR.annotate$labels[match(rownames(merged_seurat_harmony@meta.data), rownames(singleR.annotate))]


## ----- clusterprofiler setup-------------------
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(AnnotationHub)
library(ReactomePA)


#### read in data

df <- read.csv('ob1DEG.csv')
df <- read.csv('DEGS10_0vs10_1.csv')

#ranking by log2fc --- i dont like this one since it doesnt take significance into consideration

genelist <- df$avg_log2FC
names(genelist) <- df$X
head(genelist)
genelist <- sort(genelist, decreasing = T)
plot(genelist)

# for ORA
gene <- df$X
head(gene)

#ranking instead with foldchange and pvalue

genelist <- (-sign(df$avg_log2FC)*(log10(df$p_val)))
names(genelist) <- df$X
head(genelist)
genelist <- sort(genelist, decreasing = T)
plot(genelist)


## gene ontology analysis - gene set enrichment analysis
fgseaRes <- gseGO(geneList = genelist,
                  OrgDb = org.Mm.eg.db,
                  keyType = "SYMBOL",
                  ont = "BP",
                  nPermSimple = 10000,
                  minGSSize = 10,
                  maxGSSize = 1000,
                  pvalueCutoff = 0.05,
                  eps = 0,
                  verbose = F
)

write.csv(fgseaRes, file = "GSEABP2mcar1.csv")

dotplot(fgseaRes, showCategory = 20, split = ".sign", font.size = 8) + facet_grid(.~.sign) ####### find a. way to plot NES for GSEA

View(as.data.frame(fgseaRes))

### gene ontology but with over-representation

ora <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                keyType = 'SYMBOL',
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05)

write.csv(ora, file = "ORABP2mcar1.csv")

dotplot(ora, showCategory = 20)

## KEGG pathway enrichment analysis ### WIP need to convert symbols into entrezid

kegg <- gseKEGG(gene = genelist,
                organism = "mmu",
                nPermSimple = 10000,
                minGSSize = 10,
                maxGSSize = 1000,
                pvalueCutoff = 0.05,
                verbose = F
)

## reactome pathway analysis

reactome <- gsePathway(geneList = genelist,
                       organism = "mmu",
                       OrgDb = org.Mm.eg.db,
                       keyType = "SYMBOL",
                       nPermSimple = 10000,
                       minGSSize = 10,
                       maxGSSize = 1000,
                       pvalueCutoff = 0.05,
                       verbose = F
)

### number of cells within cluster that express double positive cells
Idents(recluster) = "new_clusters"
test <- subset(recluster, idents = c("2", "3"))
Idents(test) = "new_clusters"
FeatureScatter_scCustom(test, feature1 = "Pecam1", feature2 = "Bglap", pt.size = 1.5, colors_use = c("dodgerblue", "red"))

(sum(FetchData(test, vars = "Pth1r") > 0.5 & FetchData(test, vars = "Ebf3") > 0.5)) / ncol(test)

Percent_Expressing(test, features = c("Pth1r", "Ebf3"))


# umap can label in better details - first assign the variable for DimPlot then apply it into LabelClusters()
p <- DimPlot(recluster)
LabelClusters(p, id = "ident",  fontface = "bold", color = "black", repel = F, size = 4)


DimPlot(merged_seurat_harmony, label = T, split.by = "orig.ident")

FeaturePlot_scCustom(merged_seurat_harmony, features = "tdTomato", colors_use = viridis_plasma_dark_high, max.cutoff = 4)
VlnPlot(merged_seurat_harmony, features = "Tnfsf11", pt.size = 0)


#convert seurat object to anndata
SaveH5Seurat(recluster, filename = "recluster.h5Seurat")
Convert("recluster.h5Seurat", dest = "h5ad")

SaveH5Seurat(merged_seurat_harmony, filename = "merged.h5Seurat")
Convert("merged.h5Seurat", dest = "h5ad")



test <- subset(recluster, idents = c("0_0", "0_1", "4", "7", "9", "11", "12"))

test <- FindVariableFeatures(test)
test <- ScaleData(test, features = rownames(test))
test <- RunPCA(test, npcs = 30)
test <- FindNeighbors(test, dims = 1:30, reduction = 'umap')
test <- FindClusters(test, resolution = c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4), reduction = 'cca.aligned')
test <- RunUMAP(test, dims = 1:30, reduction = 'umap')
DimPlot(test, label = T)




#----- analysis
Idents(atlas) = "cell_type"

atlas_stromal <- subset(
  atlas,
  subset = cell_type == "stromal cell of bone marrow"
)

DefaultAssay(atlas_stromal) <- "RNA"

query_stromal <- subset(merged_seurat_harmony, idents = c("9"))
Idents(query_stromal) = "orig.ident"
query_pth_stromal <- subset(query_stromal, idents = c("pth1", "pth2"))
Idents(query_pth_stromal) = "RNA_snn_res.0.4"




# Make sure we're plotting real expression
DefaultAssay(atlas) <- "RNA"

# Confirm all genes exist in the atlas
ens_ids %in% rownames(atlas[["RNA"]])

# convert Ensembl to gene symbol in the ATLAS
library(org.Hs.eg.db)

ens_ids <- rownames(atlas[["RNA"]])

symbols <- mapIds(
  org.Hs.eg.db,
  keys = ens_ids,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

keep <- !is.na(symbols) & !duplicated(symbols)

atlas <- subset(atlas, features = ens_ids[keep])
rownames(atlas) <- symbols[keep]


atlas_stromal$dataset <- "atlas"
query_pth_stromal$dataset <- "query"

stromal_combined <- merge(atlas_stromal, query_pth_stromal)

Stacked_VlnPlot(stromal_combined, features = c("KITLG","EBF3", "EBF1", "FOXC1", "GNAS", "ADCY3", "ATF3"), ggplot_default_colors = T, pt.size = 0, plot_legend = T, x_lab_rotate = TRUE, add.noise = F, group.by = "dataset")


#--- reference mapping

query <- merged_seurat_harmony
reference <- atlas

DefaultAssay(reference) <- "RNA"
reference

# Find anchors to map query onto reference
anchors <- FindTransferAnchors(
  reference = atlas,
  query = merged_seurat_harmony,
  dims = 1:30  # use the PCs from atlas
)

# Map labels from atlas to query
predictions <- TransferData(
  anchorset = anchors,
  refdata = atlas$cell_type,  # or l1 / cell_type depending on resolution
  dims = 1:30
)

# Add predicted labels to your dataset
merged_seurat_harmony <- AddMetaData(merged_seurat_harmony, metadata = predictions)


# UMAP colored by your cluster
DimPlot(merged_seurat_harmony, group.by = "RNA_snn_res.0.4")

# UMAP colored by predicted atlas cell type
DimPlot(merged_seurat_harmony, group.by = "predicted.id")



merged_seurat <- merge(reference, y = c(merged_seurat_harmony),
                       add.cell.ids = c("atlas", "query"))

test <- merged_seurat_harmony

# merge the idents together from atlas and query
merged_seurat_harmony$new_id <- ifelse(
  !is.na(merged_seurat_harmony$cell_type) & merged_seurat_harmony$cell_type != "",
  merged_seurat_harmony$cell_type,
  merged_seurat_harmony$predicted.id
)

Idents(merged_seurat_harmony) = "new_id"


merged_seurat_harmony$stim <- ifelse(
  merged_seurat_harmony$orig.ident == "atlas",
  "Control",
  ifelse(
    merged_seurat_harmony$orig.ident %in% c("pth1", "pth2"),
    "Teriparatide",
    NA
  )
)


CAR <- subset(merged_seurat_harmony, idents = c("CAR1", "CAR2", "CAR3", "Pre-osteoblast"))
CARDEG <- subset(merged_seurat_harmony, idents = c("CAR1", "CAR2", "CAR3"))

Stacked_VlnPlot(CARDEG, features = c("SERPINE1", "TGFBI", "CCN2", "CCN1", "INHBA","POSTN", "FN1", "ACTA2", "TAGLN", "PLOD2", "ZEB2"), ggplot_default_colors = F, pt.size = 0, plot_legend = T, x_lab_rotate = TRUE, add.noise = F, split.by = "stim", colors_use = deg_colors)

Stacked_VlnPlot(CAR, features = c("PDE1C", "PDE7B", "PDE5A", "PDE7A", "PDE4D",), ggplot_default_colors = F, pt.size = 0, plot_legend = T, x_lab_rotate = TRUE, add.noise = F, split.by = "stim", colors_use = deg_colors)



merged_seurat_harmony <- RenameIdents(
  merged_seurat_harmony,
  "bone cell" = "Osteo-lineage",
  "CD4-positive, alpha-beta T cell" = "CD4 T",
  "CD8-positive, alpha-beta T cell" = "CD8 T",
  "dendritic cell, human" = "DC",
  "endothelial cell" = "Endothelial",
  "hematopoietic precursor cell" = "HSPC",
  "lymphocyte of B lineage" = "B",
  "mature NK T cell" = "NK",
  "monocyte" = "Monocyte",
  "smooth muscle cell" = "Smooth-muscle",
  "stromal cell of bone marrow" = "Stromal",
  "T cell" = "Other T"
)



merged_seurat_harmony$cluster_id <- merged_seurat_harmony@active.ident


Idents(merged_seurat_harmony) = "new_id"




# figures

DimPlot_scCustom(atlas, reduction = 'umap', pt.size = 0.1, label = F, ggplot_default_colors = T, raster = F, group.by = "author_celltype") %>%
  LabelClusters(id = "author_celltype", fontface = "bold", repel = F, size = 4)

DimPlot_scCustom(merged_seurat_harmony, reduction = 'umap', pt.size = 0.1, label = F, ggplot_default_colors = T, raster = F, group.by = "cluster_id") %>%
  LabelClusters(id = "cluster_id", fontface = "bold", repel = F, size = 4)

DimPlot(merged_seurat_harmony, reduction = 'umap', pt.size = 0.1, label = F, raster = F, split.by = "stim") %>%
  LabelClusters(id = "ident", fontface = "bold", repel = F, size = 4)


pal <- DiscretePalette_scCustomize(num_colors = 30, palette = "varibow")
pal2 <- DiscretePalette_scCustomize(num_colors = 40, palette = "varibow")
my_colors <- pal[c(1,2,14,15,12,6,7,8,9,17,11,5,13,18,19,16,10)]

deg_colors <- pal2[c(25,40)]


car3 <- subset(CAR, idents = "3")
car4 <- subset(CAR, idents = "4")
car5 <- subset(CAR, idents = "5")


# ----------------------------
# PTH response module
# ----------------------------
pth_genes <- c(
  "TNFSF11",
  "FOS",
  "FOSL1",
  "JUNB",
  "EGR1",
  "CREM",
  "NR4A2",
  "NR4A3",
  "NFKBIA",
  "IL6",
  "PDE4D",
  "IGF1",
  "PDE1C"
)


# ----------------------------
# TGF-beta signaling module
# ----------------------------
tgfb_signal <- c("Tgfb1","Tgfbr1","Ltbp1","Itgb8","Fstl1","Ccn5","Mmp13","Loxl2","Sfrp4","Smad7","Smad2","Smad3","Col1a1","Col3a1","Postn","Timp1","Fn1","Acta2","Ccn2","Serpine1")

tgfb_emt <- c(
  "SERPINE1",
  "FN1",
  "COL1A1",
  "COL1A2",
  "POSTN",
  "TAGLN",
  "TIMP1",
  "VIM",
  "CCN2",
  "ITGA5",
  "MMP2",
  "MMP14",
  "PLOD2",
  "TGM2",
  "ROCK2",
  "LTBP1",
  "LOXL2"
)


CAR <- AddModuleScore(
  CAR,
  features = list(pth_genes),
  name = "PTH_Score",
  ctrl = 50
)

CARDEG <- AddModuleScore(
  CARDEG,
  features = list(tgfb_emt),
  name = "TGFB_Score",
  ctrl = 50
)


Stacked_VlnPlot(
  CARDEG,
  features = c("TGFB_Score"),
  group.by = "recluster_ids",
  split.by = "stim",
  add.noise = F, colors_use = deg_colors,
  pt.size = 0
)


FeaturePlot_scCustom(
  CAR,
  features = c("PTH_Score1"),
  split.by = "stim",
  reduction = "umap",
  raster = F,
  label = T,
  repel = T
)

# ---- renaming clusters

current_ids <- levels(Idents(merged_seurat_harmony))
current_ids

protected_ids <- c("CAR1", "CAR2", "CAR3", "Pre-osteoblast")


remaining_ids <- setdiff(current_ids, protected_ids)

# sort for reproducibility (important)
remaining_ids <- sort(remaining_ids)

new_numeric_labels <- as.character(seq_along(remaining_ids) - 1)

rename_map <- setNames(new_numeric_labels, remaining_ids)

merged_seurat_harmony <- RenameIdents(merged_seurat_harmony, rename_map)

final_levels <- c(
  protected_ids,
  as.character(seq_along(remaining_ids) - 1)
)

Idents(merged_seurat_harmony) <- factor(
  Idents(merged_seurat_harmony),
  levels = final_levels
)


DimPlot(
  merged_seurat_harmony,
  cells.highlight = list(
    pth1 = WhichCells(merged_seurat_harmony, expression = orig.ident == "pth1"),
    pth2 = WhichCells(merged_seurat_harmony, expression = orig.ident == "pth2")
  ),
  cols.highlight = c("#E41A1C", "#1F78B4"),
  cols = "grey85",
  pt.size = .5,
  raster = F
)

DimPlot(
  merged_seurat_harmony,
  cells.highlight = list(
    Teriparatide = WhichCells(merged_seurat_harmony, expression = stim == "Teriparatide")
    #pth2 = WhichCells(merged_seurat_harmony, expression = orig.ident == "pth2")
  ),
  cols.highlight = c("#E41A1C"),
  cols = "grey85",
  pt.size = .5,
  raster = F
)




#=======

merged_seurat_harmony$pth_highlight <- "Atlas"
merged_seurat_harmony$pth_highlight[merged_seurat_harmony$orig.ident == "pth1"] <- "pth1"
merged_seurat_harmony$pth_highlight[merged_seurat_harmony$orig.ident == "pth2"] <- "pth2"


DimPlot_scCustom(
  merged_seurat_harmony,
  group.by = "pth_highlight",
  colors_use = c(
    "pth1"  = "#D55E00",
    "pth2"  = "#0072B2",
    "Other" = "grey100"
  )
)



# module score effect size testing just between two comparisons======
library(dplyr)

df <- CARDEG@meta.data %>%
  select(
    cluster = recluster_ids,
    stim,
    TGFB_Score
  ) %>%
  filter(stim %in% c("Control", "Teriparatide"))


# --- running it by downsampling

# Step 1: get min cells per cluster
min_n <- df %>%
  dplyr::count(cluster, stim) %>%
  group_by(cluster) %>%
  summarise(n_min = min(n), .groups = "drop")

# Step 2: downsample per cluster
df_balanced <- df %>%
  left_join(min_n, by = "cluster") %>%
  group_by(cluster) %>%
  group_modify(~ {
    n_use <- unique(.x$n_min)
    bind_rows(
      slice_sample(filter(.x, stim == "Control"), n = n_use),
      slice_sample(filter(.x, stim == "Teriparatide"),     n = n_use)
    )
  }) %>%
  ungroup() %>%
  select(-n_min)

#sanity check
df_balanced %>% dplyr::count(cluster, stim)


# =============================
# 1️⃣ Wilcoxon test per cluster
# =============================
wilcox_results <- df_balanced %>%
  group_by(cluster) %>%
  summarise(
    # Wilcoxon rank-sum test p-value
    p_value = tryCatch(
      wilcox.test(TGFB_Score ~ stim, exact = FALSE)$p.value,
      error = function(e) NA
    ),
    # median scores per condition
    median_control = median(TGFB_Score[stim == "Control"], na.rm = TRUE),
    median_PTH     = median(TGFB_Score[stim == "Teriparatide"], na.rm = TRUE),
    # delta median
    delta_median   = median_PTH - median_control,
    # number of cells per condition (should be equal after downsampling)
    n_control = sum(stim == "Control"),
    n_PTH     = sum(stim == "Teriparatide"),
    .groups = "drop"
  ) %>%
  # optional: adjust p-values for multiple clusters
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# =============================
# 2️⃣ Cliff's delta per cluster
# =============================
cliff_results <- df_balanced %>%
  group_by(cluster) %>%
  summarise(
    cliffs_delta = tryCatch(
      cliff.delta(
        TGFB_Score[stim == "Teriparatide"],
        TGFB_Score[stim == "Control"]
      )$estimate,
      error = function(e) NA
    ),
    .groups = "drop"
  )

# =============================
# 3️⃣ Merge Wilcoxon and Cliff's delta results
# =============================
final_results <- wilcox_results %>%
  left_join(cliff_results, by = "cluster") %>%
  mutate(
    # categorize effect size
    cliffs_magnitude = case_when(
      abs(cliffs_delta) < 0.147 ~ "negligible",
      abs(cliffs_delta) < 0.33  ~ "small",
      abs(cliffs_delta) < 0.474 ~ "medium",
      TRUE                      ~ "large"
    )
  )

# =============================
# 4️⃣ Inspect final results
# =============================
final_results

write.csv(final_results, file = "final_TGFB_module_score.csv")







## LIGER integration

library(rliger)

expr_list <- list(
  atlas = GetAssayData(atlas_filtered, assay = "RNA", slot = "counts"),
  pth1  = GetAssayData(pth1_filtered,  assay = "RNA", slot = "counts"),
  pth2  = GetAssayData(pth2_filtered,  assay = "RNA", slot = "counts"),
  veh1  = GetAssayData(veh1_filtered,  assay = "RNA", slot = "counts"),
  veh2  = GetAssayData(veh2_filtered,  assay = "RNA", slot = "counts")
)

liger_obj <- rliger::createLiger(expr_list)

liger_obj <- rliger::normalize(liger_obj)

liger_obj <- rliger::selectGenes(
  liger_obj,
  var.thresh = 0,
  combine = "union"
)




liger_obj <- rliger::scaleNotCenter(liger_obj)


liger_obj <- runIntegration(
  liger_obj,
  k = 15,        # similar to ~15 PCs
  lambda = 4    # moderate integration strength
)

liger_obj <- quantile_norm(liger_obj)

liger_obj <- runUMAP(
  liger_obj,
  distance = "cosine",
  n_neighbors = 30,
  min_dist = 0.3
)

liger_obj <- louvainCluster(
  liger_obj,
  resolution = 1
)


plotByDatasetAndCluster(liger_obj, axis.labels = c("UMAP1", "UMAP2"))

plotByCluster(liger_obj)

# View gene loadings for a factor of interest
getFactorGeneLoading(liger_obj, factor = 5)


liger_seurat <- ligerToSeurat(liger_obj)


DimPlot(liger_seurat, group.by = "dataset")
DimPlot(liger_seurat, group.by = "liger_clusters")




getMatrix(liger_obj, slot = "norm")



# Seurat > GEO submission workflow

merged_seurat_harmony <- JoinLayers(merged_seurat_harmony, assay = "RNA")
counts <- GetAssayData(merged_seurat_harmony, assay = "RNA", layer = "counts")

library(Matrix)

counts <- as(counts, "dgCMatrix")

# Matrix
writeMM(counts, file = "matrix.mtx")

# Genes (features)
write.table(
  rownames(counts),
  file = "genes.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

# Barcodes (cells)
write.table(
  colnames(counts),
  file = "barcodes.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

meta <- merged_seurat_harmony@meta.data

write.csv(meta, file = "cell_metadata.csv")










