# script to perform standard workflow steps to analyze single cell RNA-Seq data
# data: 2 days and 2 months Ebf3-creERT;tdt chase with PTH vs vehicle bone marrow cells

# setwd("~/Documents/scRNA_seq_analysis")
# setwd("~/Documents/scRNA_seq_analysis/data/new_files/")


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
library(cowplot)
library(effsize)
library(UCell)




# files legends
# HKd1 - female vehicle 2 days
# HKd2 - male vehicle 2 days
# HKd4 - female PTH 2 days - i suspect vehicle and PTH may have been mixed
# HKd5 - male PTH 2 months
# HKd6 - female PTH 2 months
# HKd7 - male vehicle 2 months
# HKd8 - female vehicle 2 months

#we lost male PTH and i know this is HKd3
# potentially it is HKd1 male vehicle!!!!!
# HKd2 female PTH!!!
# HKd3 male PTH
# HKd4 female vehicle!!!

#repeat 2 days and 10 days vehicle and PTH chase

# Load the dataset

mveh2d = Read10X_h5(filename = '~/Documents/Sox9_PTH_scRNA_analysis/Matrix/HKd17/filtered_feature_bc_matrix.h5')
fveh2d = Read10X_h5(filename = '~/Documents/Sox9_PTH_scRNA_analysis/Matrix/HKd18/filtered_feature_bc_matrix.h5')
mpth2d = Read10X_h5(filename = '~/Documents/Sox9_PTH_scRNA_analysis/Matrix/HKd19/filtered_feature_bc_matrix.h5')
fpth2d = Read10X_h5(filename = '~/Documents/Sox9_PTH_scRNA_analysis/Matrix/HKd20/filtered_feature_bc_matrix.h5')

mveh3w = Read10X_h5(filename = '~/Documents/Sox9_PTH_scRNA_analysis/Matrix/HKd13/filtered_feature_bc_matrix.h5')
fveh3w = Read10X_h5(filename = '~/Documents/Sox9_PTH_scRNA_analysis/Matrix/HKd14/filtered_feature_bc_matrix.h5')
mpth3w = Read10X_h5(filename = '~/Documents/Sox9_PTH_scRNA_analysis/Matrix/HKd15/filtered_feature_bc_matrix.h5')
fpth3w = Read10X_h5(filename = '~/Documents/Sox9_PTH_scRNA_analysis/Matrix/HKd16/filtered_feature_bc_matrix.h5')


mveh2m = Read10X_h5(filename = '~/Documents/Sox9_PTH_scRNA_analysis/Matrix/HKd9/filtered_feature_bc_matrix.h5')
fveh2m = Read10X_h5(filename = '~/Documents/Sox9_PTH_scRNA_analysis/Matrix/HKd10/filtered_feature_bc_matrix.h5')
mpth2m = Read10X_h5(filename = '~/Documents/Sox9_PTH_scRNA_analysis/Matrix/HKd11/filtered_feature_bc_matrix.h5')
fpth2m = Read10X_h5(filename = '~/Documents/Sox9_PTH_scRNA_analysis/Matrix/HKd12/filtered_feature_bc_matrix.h5')

veh



#after sorting out genes - looks like HKd1 and HKd2 are swapped - HKd2 is actually female - based on DEG showing Xist is upregulated


# Create one merged Seurat object, initialize the Seurat object with the raw (non-normalized data).

mveh2d = CreateSeuratObject(counts = mveh2d, project = "male_vehicle_2d", min.cells = 3, min.features = 200)
fveh2d = CreateSeuratObject(counts = fveh2d, project = "female_vehicle_2d", min.cells = 3, min.features = 200)
mpth2d = CreateSeuratObject(counts = mpth2d, project = "male_pth_2d", min.cells = 3, min.features = 200)
fpth2d = CreateSeuratObject(counts = fpth2d, project = "female_pth_2d", min.cells = 3, min.features = 200)

mveh3w = CreateSeuratObject(counts = mveh3w, project = "male_vehicle_3w", min.cells = 3, min.features = 200)
fveh3w = CreateSeuratObject(counts = fveh3w, project = "female_vehicle_3w", min.cells = 3, min.features = 200)
mpth3w = CreateSeuratObject(counts = mpth3w, project = "male_pth_3w", min.cells = 3, min.features = 200)
fpth3w = CreateSeuratObject(counts = fpth3w, project = "female_pth_3w", min.cells = 3, min.features = 200)

mveh2m = CreateSeuratObject(counts = mveh2m, project = "male_vehicle_2m", min.cells = 3, min.features = 200)
fveh2m = CreateSeuratObject(counts = fveh2m, project = "female_vehicle_2m", min.cells = 3, min.features = 200)
mpth2m = CreateSeuratObject(counts = mpth2m, project = "male_pth_2m", min.cells = 3, min.features = 200)
fpth2m = CreateSeuratObject(counts = fpth2m, project = "female_pth_2m", min.cells = 3, min.features = 200)







# 1. QC -------
# % MT reads
mveh2d$percent.mt = PercentageFeatureSet(mveh2d, pattern = "^mt-")
fveh2d$percent.mt = PercentageFeatureSet(fveh2d, pattern = "^mt-")
mpth2d$percent.mt = PercentageFeatureSet(mpth2d, pattern = "^mt-")
fpth2d$percent.mt = PercentageFeatureSet(fpth2d, pattern = "^mt-")

mveh3w$percent.mt = PercentageFeatureSet(mveh3w, pattern = "^mt-")
fveh3w$percent.mt = PercentageFeatureSet(fveh3w, pattern = "^mt-")
mpth3w$percent.mt = PercentageFeatureSet(mpth3w, pattern = "^mt-")
fpth3w$percent.mt = PercentageFeatureSet(fpth3w, pattern = "^mt-")

mveh2m$percent.mt = PercentageFeatureSet(mveh2m, pattern = "^mt-")
fveh2m$percent.mt = PercentageFeatureSet(fveh2m, pattern = "^mt-")
mpth2m$percent.mt = PercentageFeatureSet(mpth2m, pattern = "^mt-")
fpth2m$percent.mt = PercentageFeatureSet(fpth2m, pattern = "^mt-")



VlnPlot(mveh2d, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(mveh2d, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

VlnPlot(fveh2d, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(fveh2d, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

VlnPlot(mpth2d, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(mpth2d, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

VlnPlot(fpth2d, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(fpth2d, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')


VlnPlot(mveh3w, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(mveh3w, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

VlnPlot(fveh3w, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(fveh3w, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

VlnPlot(mpth3w, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(mpth3w, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

VlnPlot(fpth3w, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(fpth3w, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')



VlnPlot(mveh2m, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(mveh2m, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

VlnPlot(fveh2m, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(fveh2m, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

VlnPlot(mpth2m, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(mpth2m, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

VlnPlot(fpth2m, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(fpth2m, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')







# 2. Filtering ----------------- you want cells that express more than 200 genes and less than 20 percent mito genes


mveh2d_filtered <- subset(mveh2d, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
fveh2d_filtered <- subset(fveh2d, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
mpth2d_filtered <- subset(mpth2d, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
fpth2d_filtered <- subset(fpth2d, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)

mveh3w_filtered <- subset(mveh3w, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
fveh3w_filtered <- subset(fveh3w, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
mpth3w_filtered <- subset(mpth3w, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
fpth3w_filtered <- subset(fpth3w, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)

mveh2m_filtered <- subset(mveh2m, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
fveh2m_filtered <- subset(fveh2m, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
mpth2m_filtered <- subset(mpth2m, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
fpth2m_filtered <- subset(fpth2m, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)

young.recluster <- recluster


# merge only after QC and filtering
# merging all the files
merged_seurat = merge(mveh2d_filtered, y = c(fveh2d_filtered, mpth2d_filtered, mveh3w_filtered, fveh3w_filtered, mpth3w_filtered, fpth3w_filtered, mveh2m_filtered, fveh2m_filtered, mpth2m_filtered, fpth2m_filtered),
                      add.cell.ids = c("mveh2d", "fveh2d","mpth2d","mveh3w","fveh3w","mpth3w","fpth3w", "mveh2m", "fveh2m","mpth2m", "fpth2m"))

merged_seurat = merge(veh10d_filtered, y = c(pth10d_filtered),
                      add.cell.ids = c("veh", "pth"))




# Create sample column

merged_seurat_harmony$sample = rownames(merged_seurat_harmony@meta.data)
merged_seurat_harmony@meta.data = separate(merged_seurat_harmony@meta.data, col = 'sample', into = c('Condition', 'Barcode'), sep = '_')


# Create a stimu column that groups your vehicle versus PTH --- from now you can split by "stim"
control_samples <- c("mveh2d","fveh2d","mveh3w", "fveh3w","mveh2m", "fveh2m")
pth_samples <- c("mpth2d", "mpth3w", "fpth3w", "mpth2m", "fpth2m")

control_samples <- c("veh")
pth_samples <- c("pth")

merged_seurat_harmony@meta.data$stim <- sapply(merged_seurat_harmony@meta.data$Condition, function(ita) ifelse(ita %in% control_samples, "Control", "PTH"))

merged_seurat_harmony2@meta.data$stim <- sapply(merged_seurat_harmony2@meta.data$Condition, function(ita) ifelse(ita %in% control_samples, "Control", "PTH"))

# create a time column that groups 2 days and 2 months of vehicle --- from now you can split by "time"

samples_2d <- c("mvch2d", "fvch2d", "fpth2d")
samples_2m <- c("mvch2m", "fvch2m", "mpth2m", "fpth2m")

samples_2d <- c("veh2d", "pth2d")
samples_10d <- c("veh10d", "pth10d")

merged_seurat_harmony@meta.data$time <- sapply(merged_seurat_harmony@meta.data$Condition, function(ita) ifelse(ita %in% samples_2d, "2d", "2m"))

merged_seurat_harmony2@meta.data$time <- sapply(merged_seurat_harmony2@meta.data$Condition, function(ita) ifelse(ita %in% samples_2d, "2d", "10d"))

recluster@meta.data$time <- sapply(recluster@meta.data$Condition, function(ita) ifelse(ita %in% samples_2d, "2d", "10d"))

# Create a time column that groups your vehicle versus PTH --- from now you can split by "stim"
samples_2d <- c("mveh2d", "fveh2d", "mpth2d")
samples_3w <- c("mveh3w", "fveh3w", "mpth3w", "fpth3w")
samples_2m <- c("mveh2m", "fveh2m", "mpth2m", "fpth2m")


recluster@meta.data$time <- sapply(
  recluster@meta.data$Condition,
  function(ita) {
    if (ita %in% samples_2d) {
      "2d"
    } else if (ita %in% samples_3w) {
      "3w"
    } else if (ita %in% samples_2m) {
      "2m"
    } else {
      NA
    }
  }
)



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

merged_seurat = RunPCA(merged_seurat, npcs = 50)

# determine dimensionality of the data
ElbowPlot(merged_seurat, ndims = 50)

# 7. Clustering ------------ choose the dimension based on your elbow plot
merged_seurat = FindNeighbors(merged_seurat, dims = 1:30)

# understanding resolution (lower the number the lower the clusters) - gives you resolution for distinct clusters you can play around with this

merged_seurat = FindClusters(merged_seurat, resolution = c(0.1,0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))



#8 non-linear dimensionality reduction -------------- (umap) ----- visualizing the clusters

merged_seurat <- RunUMAP(merged_seurat, dims = 1:30)

Idents(merged_seurat) = "RNA_snn_res.0.2"

p1 = DimPlot(merged_seurat, reduction = "umap", label = T, label.size = 8)

# look at whether we need batch correction or not
p2 = DimPlot(merged_seurat, split.by = "orig.ident", label = T)
p1|p2

# setting your sample names as idents so you can use the idents argument to call your sample names
Idents(merged_seurat) <- "orig.ident"
Idents(merged_seurat_harmony) <- "orig.ident"

plot.male.vehicle = subset(merged_seurat_harmony, idents = c("male_vehicle_2d", "male_vehicle_2m"))
plot.female.vehicle = subset(merged_seurat_harmony, idents = c("female_vehicle_2d", "female_vehicle_2m"))
plot.male.female.2m = subset(merged_seurat_harmony, idents = c("male_vehicle_2m", "female_vehicle_2m"))
plot.vehicle.2days = subset(merged_seurat_harmony, idents = c("female_vehicle_2d", "female_pth_2d"))
plot.male.pth.2m = subset(merged_seurat_harmony, idents = c("male_vehicle_2m", "male_pth_2m"))
plot.female.pth.2m = subset(merged_seurat_harmony, idents = c("female_vehicle_2m", "female_pth_2m"))
plot.vehicle.2m = subset(merged_seurat_harmony, i)

plot.2d = subset(merged_seurat_harmony, idents = c("male_vehicle_2d", "female_vehicle_2d", "female_pth_2d"))
umap.2d = DimPlot(plot.2d, reduction = "umap", group.by = "RNA_snn_res.0.2", label = T)

plot.2m = subset(merged_seurat_harmony, idents = c("male_vehicle_2m", "female_vehicle_2m", "male_pth_2m", "female_pth_2m"))
umap.2m = DimPlot(plot.2m, reduction = "umap", group.by = "RNA_snn_res.0.2", label = T)

umap.2d|umap.2m

DimPlot(plot.male.vehicle, reduction = "umap", label = T)
DimPlot(plot.female.vehicle, reduction = "umap", label = T)
DimPlot(plot.male.female.2m, reduction = "umap", label = T)
DimPlot(plot.vehicle.2days, reduction = "umap", label = T)
DimPlot(plot.male.pth.2m, reduction = "umap", label = T)
DimPlot(plot.female.pth.2m, reduction = "umap", label = T)

split_merged_seurat = SplitObject(merged_seurat, split.by = "orig.ident")
p3 = DimPlot(split_merged_seurat$female_vehicle_2d, reduction = "umap", label = T)
p4 = DimPlot(split_merged_seurat$female_vehicle_2m, reduction = "umap", label = T)
p3|p4

FeaturePlot(merged_seurat, features = "tdTomato")
VlnPlot(merged_seurat, features = "tdTomato")

# clustree to determine optimal clustering

clustree(merged_seurat_harmony, prefix = "RNA_snn_res.", node_colour = "sc3_stability", label.size = 4)

# running doublet finder ---- this should be ran before harmony! but you can decide if its worth doing doublet finder first


# number of cells per cluster at each condition

test1 <- View(table(Idents(sub_ob), sub_ob$orig.ident))
write.csv(test, file = "testchart.csv", row.names = T)


write.csv(table(Idents(recluster), recluster$orig.ident), file = "cellnumber.csv", row.names = T)

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

Idents(merged_seurat_harmony) = "RNA_snn_res.0.2"

p4 = DimPlot(merged_seurat_harmony, reduction = 'umap', group.by = 'RNA_snn_res.1.4', label = T)
p5 = DimPlot(merged_seurat_harmony, reduction = 'umap', split.by = 'orig.ident', label = T)
p6 = DimPlot(merged_seurat, reduction = 'umap', group.by = 'orig.ident', label = T)

# before and after batch correction
p2|p5
p1|p4
p5|p6

# visualize PCA results
print(merged_seurat_filtered[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(merged_seurat_filtered, dims = 1, cells = 500, balanced = TRUE)

#filtering object
merged_seurat_harmony2_filtered <- subset(merged_seurat_harmony2, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 20)


# Finding DEGs (cluster biomarkers) ----------------- Find all markers will try to define differentially expressed markers for each clusters compared to rest
merged_seurat_harmony[["RNA"]] <- JoinLayers(merged_seurat_harmony[["RNA"]]) # join the layers this is new in seurat v5

markers = FindAllMarkers(recluster, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
markers.clusters = markers %>%
  group_by(cluster) %>%
  slice_max(n = 500, order_by = avg_log2FC)

write.csv(markers.clusters, file = "all.recluster.markers.csv")

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

FeaturePlot(merged_seurat_harmony, features = c("Ebf3","Cxcl12","Lepr","Sp7", "tdTomato","Adipoq", "Bglap", "Col1a1", "Sox9"))

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


# Find cells that have the gene in this case tdTomato
tdtomato_1_more = WhichCells(merged_seurat_harmony, slot = "data", expression = tdTomato>1)
tdtomato_1_equal = WhichCells(merged_seurat_harmony, slot = "data", expression = tdTomato==1)
tdtomato_1_less = WhichCells(merged_seurat_harmony, slot = "data", expression = tdTomato<1)

# Look at the gene "tdTomato" distribution of your data
summary(FetchData(merged_seurat_harmony, vars = "tdTomato"))
RidgePlot(merged_seurat_harmony, features = 'tdTomato')

# look at number of cells expressing a gene 'tdTomato'

sum(GetAssayData(object = merged_seurat_harmony, slot = "data")['tdTomato',])

sum(GetAssayData(object = merged_seurat_harmony, slot = "data")['tdTomato',]>0)/nrow(merged_seurat_harmony@meta.data)

sum(nrow(merged_seurat_harmony))

sum(GetAssayData(object = merged_seurat, slot = "data")['tdTomato',]>0)/nrow(merged_seurat_harmony@meta.data)


# trying to split the data to plot by individual conditions
pdf(file = "conditions.pdf",width = 6,height = 4)
conditions.tdtomato = SplitObject(tdtomato, split.by = "orig.ident")
plot.conditions.tdtomato = lapply(X = conditions.tdtomato, FUN = function(x) {
  DimPlot(x, reduction = "umap", label = T, label.size = 4)
})
plot.test
dev.off()


# finding markers, only.pos = T makes it higher expression compared to lower

markers = FindAllMarkers(merged_seurat_harmony, only.pos = T, min.pct = 0.5, min.diff.pct = 0.5, logfc.threshold = 0.25)
markers.clusters = markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)
write.csv(markers.clusters, file = "old.markers.clusters.csv")



# finding subcluster in specific cluster

merged_seurat_harmony <- FindSubCluster(merged_seurat_harmony, "6", "RNA_snn", subcluster.name = "test_sub", resolution = 0.05, algorithm = 1)
DimPlot(merged_seurat_harmony, reduction = "umap", group.by = "test_sub", label = T)


test <- recluster
test <- test_backup
Idents(recluster) = "test_sub"

# how many cells are in each cluster
Idents(merged_seurat_harmony) = "RNA_snn_res.0.2"
table(Idents(merged_seurat_harmony))

Idents(merged_seurat_harmony) = "time"
table(Idents(merged_seurat_harmony))

table(Idents(merged_seurat_harmony), merged_seurat_harmony$orig.ident)


# how many cells are in each cluster per time point

Idents(merged_seurat_harmony) = "orig.ident"
cluster2d = subset(merged_seurat_harmony, idents = c("male_vehicle_2d", "female_vehicle_2d", "female_pth_2d"))
cluster2m = subset(merged_seurat_harmony, idents = c("male_vehicle_2m", "female_vehicle_2m", "male_pth_2m", "female_pth_2m"))

Idents(cluster2d) = "RNA_snn_res.0.2"
table(Idents(cluster2d))

Idents(cluster2m) = "RNA_snn_res.0.2"
table(Idents(cluster2m))

# proper plotting with feature plots
FeaturePlot_scCustom(merged_seurat, features = "Vdr", colors_use = viridis_plasma_dark_high)

VlnPlot(merged_seurat, features = "tdTomato")


FeaturePlot_scCustom(merged_seurat_harmony, features = "Vdr", colors_use = viridis_plasma_dark_high)

VlnPlot_scCustom(merged_seurat_harmony, features = "tdTomato", ggplot_default_colors = T, pt.size = 0)

FeaturePlot_scCustom(merged_seurat_harmony, features = c("Col1a1", "Runx2", "Bglap", "Col1a1", "Sp7", "Spp1", "Alpl","Entpd3", "Car3", "Dmp1", "Mepe", "Phex"), colors_use = viridis_plasma_dark_high, max.cutoff = 4)

VlnPlot_scCustom(merged_seurat_harmony, features = c("Ebf3","Cxcl12", "Lepr", "Adipoq", "Lpl", "Runx2", "Bglap", "Col1a1", "Sp7", "Spp1", "Dmp1", "Alpl"), idents = c("Osteo-CAR", "Adipo-CAR"), pt.size = 0)

Idents(merged_seurat_harmony1) = "RNA_snn_res.0.2"
DimPlot_scCustom(merged_seurat_harmony, reduction = 'umap', colors_use = color_pal_list, pt.size = 0.5, label = F) %>%
  LabelClusters(id = "ident", fontface = "bold", repel = F)

plot.2d = subset(merged_seurat_harmony, idents = c("male_vehicle_2d", "female_vehicle_2d", "female_pth_2d"))
Idents(plot.2d) = "RNA_snn_res.0.2"
umap2d <- DimPlot_scCustom(plot.2d, reduction = "umap", colors_use = color_pal_list, pt.size = 0.5, label = F) %>%
  LabelClusters(id = "ident", fontface = "bold", repel = F)

plot.2m = subset(merged_seurat_harmony, idents = c("male_vehicle_2m", "female_vehicle_2m", "male_pth_2m", "female_pth_2m"))
Idents(plot.2m) = "RNA_snn_res.0.2"
umap2m <- DimPlot_scCustom(plot.2m, reduction = "umap", colors_use = color_pal_list, pt.size = 0.5, label = F) %>%
  LabelClusters(id = "ident", fontface = "bold", repel = F)


# custom palette list for umap
color_pal_list <- c("#d741a7", "#afa2ff", "#5398be", "#dea54b", "#ef2917", "#a30b37", "#87b38d","#aef3e7", "#8ee3ef", "#72788d", "#bbb6df", "#f2b7c6", "#f2dc5d", "#f18701", "#b56b45", "#37ff8b")

# rename clusters and add it to metadata under column cluster_ids

merged_seurat_harmony <- RenameIdents(merged_seurat_harmony, '0' = "Erythroid progenitor", '1' = "EC1", '2' = "CAR1", '3' = "Myeloid", '4' = "Myeloid progenitor", '5' = "Neutrophil", '6' = "Erythroid precursor", '7' = "B", '8' = "Skeletal muscle", '9' = "Erythroblast", '10' = "Osteoblast", '11' = "CAR2", '12' = "Mast", '13' = "Skeletal progenitor", '14' = "EC2", '15' = "EC3")
merged_seurat_harmony$cluster_ids <- merged_seurat_harmony@active.ident

merged_seurat_harmony2 <- RenameIdents(merged_seurat_harmony2, '0' = "OsteoCAR1", '1' = "OsteoCAR2", '2' = "Pericytes", '3' = "Periosteal1", '4' = "Periosteal2", '5' = "Osteoblast1", '6' = "Endothelial", '7' = "AdipoCAR", '8' = "MSC progenitor", '9' = "Osteoblast2", '10' = "Periosteal3")
merged_seurat_harmony2$cluster_ids <- merged_seurat_harmony2@active.ident

recluster1 <- RenameIdents(recluster1, '0' = "OsteoCAR1", '1' = "OsteoCAR2", '2' = "OsteoCAR3", '3' = "OsteoCAR4", '4' = "Periosteal1", '5' = "Periosteal2", '6' = "Pre-osteoblast", '7' = "Endothelial", '8' = "AdipoCAR", '9' = "Pericyte1", '10' = "Pericyte2", '11' = "Periosteal3", '12' = "Myogenic", '13' = "Mature osteoblast")
recluster1$cluster_ids <- recluster1@active.ident
Idents(recluster1) <- "cluster_ids"

levels(recluster1) <- c("OsteoCAR1", "OsteoCAR2", "OsteoCAR3", "OsteoCAR4", "AdipoCAR", "Pre-osteoblast", "Mature osteoblast", "Periosteal1", "Periosteal2", "Periosteal3", "Pericyte1", "Pericyte2", "Endothelial", "Myogenic")

recluster <- RenameIdents(recluster, 'AdipoCAR' = "CAR1", 'OsteoCAR4' = "CAR2", 'OsteoCAR1' = "CAR3", 'OsteoCAR3' = "CAR4", 'OsteoCAR2' = "CAR5", 'Pre-osteoblast_0' = "CAR6", 'Pre-osteoblast_1' = "Pre-osteoblast")


#rename clusters by cell annotation

cell.annotation <- c("Erythroid1", "EC1", "CAR1", "Myeloid", "Myeloid progenitor", "Neutrophil", "Erythroid precursor", "B", "Skeletal muscle", "Erythroblast", "Osteoblast", "CAR2", "Mast", "Skeletal progenitor", "EC2", "EC3")
names(cell.annotation) <- levels(merged_seurat_harmony)
merged_seurat_harmony <- RenameIdents(merged_seurat_harmony, cell.annotation)
DimPlot(merged_seurat_harmony, reduction = "umap", label = T)

export_df <- merged_seurat_harmony@meta.data

write.csv(merged_seurat_harmony@meta.data, "cell_annotation_metadata.csv")

# Reclustering data to remove HSCs and low tdtomato count --- we want clusters EC1, CAR1, myeloid, myeloid progenitors, neutrophil, osteoblasts, car2, mast, skeletal mscs, ec2 and ec3

recluster <- subset(merged_seurat_harmony, idents = c("3", "5", "6", "7","8","9","10","11", "12", "13"))

recluster <- RenameIdents(recluster, '0' = "Erythroid progenitor", '1' = "EC1", '2' = "CAR1", '3' = "Myeloid", '4' = "Myeloid progenitor", '5' = "Neutrophil", '6' = "Erythroid precursor", '7' = "B", '8' = "Skeletal muscle", '9' = "Erythroblast", '10' = "Osteoblast", '11' = "CAR2", '12' = "Mast", '13' = "Skeletal progenitor", '14' = "EC2", '15' = "EC3")
recluster$new_cluster_ids <- recluster@active.ident
recluster <- RenameIdents(recluster, '7_3' = "7_2")


# removing the metadata column

recluster$new_clusters_ids <- NULL

# finding subclusters =====================================
#1 look at what clusters you want to subclusters at different resolutions
#2 proceed to use findsubcluster() to create new meta data column, then set identity to new meta data column
#3 combine certain clusters if need be by using renameidents()
#4 once you are happy with the clusters, add this ident into new meta data column

recluster <- FindSubCluster(recluster, "14", "RNA_snn", subcluster.name = "new_clusters", resolution = 0.1, algorithm = 1)
Idents(recluster) = "new_clusters"
DimPlot(recluster, reduction = "umap", group.by = "new_clusters", label = T)

# merging some subclusters =======================

recluster <- RenameIdents(recluster, 'Pre-osteoblast_3' = "Pre-osteoblast")
recluster <- RenameIdents(recluster, '2_2' = "2_1")
table(Idents(recluster))

recluster$new_clusters <- recluster@active.ident
Idents(recluster) = "new_clusters"


# saving old recluster incase anything goes wrong
recluster.old <- recluster

# no need for reclustering ====================== only necessary if you remove cells - but in this case we removed whole clusters
recluster1 <- FindVariableFeatures(recluster)
recluster1 <- ScaleData(recluster1)
recluster1 <- RunPCA(recluster1, npcs = 20)
recluster1 <- FindNeighbors(recluster1, dims = 1:20)
recluster1 <- FindClusters(recluster1, resolution = c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))
recluster1 <- RunUMAP(recluster1, dims = 1:20)

# need to rerun harmony again....
recluster1 <- RunHarmony(recluster1, group.by.vars = 'orig.ident', max.iter = 20, plot_convergence = T)

recluster1 <- recluster1 %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
  FindClusters(resolution = c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))

Idents(recluster1) = "RNA_snn_res.0.1"

DimPlot_scCustom(recluster, reduction = 'umap', pt.size = 0.1, label = F, ggplot_default_colors = F, colors_use = my_colors) %>%
  LabelClusters(id = "ident", fontface = "bold", repel = T, size = 0)

Stacked_VlnPlot(CAR, features = c("Tnfrsf11a"), ggplot_default_colors = F , pt.size = 0, x_lab_rotate = T, add.noise = F, colors_use = my_colors)

DotPlot_scCustom(CAR, features = "Pth1r", x_lab_rotate=T, colors_use = viridis_plasma_dark_high) + scale_size(range = c(4,8)) + coord_flip()

Stacked_VlnPlot(merged_seurat_harmony, features = c("Cxcl12", "Lepr", "Ebf3", "Adipoq", "Foxc1", "Wif1", "Alpl", "Spp1","Col1a1", "Bglap", "Dmp1","Pth1r", "Prg4", "Acan", "Col2a1","Scx"), ggplot_default_colors = T , pt.size = 0, x_lab_rotate = T)

pal <- DiscretePalette_scCustomize(num_colors = 30, palette = "varibow")
pal2 <- DiscretePalette_scCustomize(num_colors = 40, palette = "varibow")
my_colors <- pal[c(1,2,14,15,12,6,7,8,9,17,11,5,13,18,19,16,10)]

deg_colors <- pal2[c(25,40)]

my_colors <- c(
  "CAR1" = pal[11],
  "CAR2" = pal[13],
  "CAR3" = pal[15],
  "CAR4" = pal[17],
  "CAR5" = pal[18],
  "CAR6" = pal[19],
  "Pre-osteoblast" = pal[1],
  "Mature osteoblast" = pal[4],
  "Periosteal1" = pal[21],
  "Periosteal2" = pal[22],
  "Periosteal3" = pal[24],
  "Pericyte1" = pal[29],
  "Pericyte2" = pal[30],
  "Endothelial" = pal[3],
  "Myogenic" = pal[6]
)



Idents(recluster) = "orig.ident"
recluster.2d = subset(recluster, idents = c("male_vehicle_2d", "female_vehicle_2d", "female_pth_2d"))
Idents(recluster.2d) = "RNA_snn_res.0.1"
umap2d <- DimPlot_scCustom(recluster.2d, reduction = "umap", colors_use = color_pal_list, pt.size = 0.5, label = F) %>%
  LabelClusters(id = "ident", fontface = "bold", repel = F)

recluster.2m = subset(recluster, idents = c("male_vehicle_2m", "female_vehicle_2m", "male_pth_2m", "female_pth_2m"))
Idents(recluster.2m) = "RNA_snn_res.0.1"
umap2m <- DimPlot_scCustom(recluster.2m, reduction = "umap", colors_use = color_pal_list, pt.size = 0.5, label = F) %>%
  LabelClusters(id = "ident", fontface = "bold", repel = F)


FeaturePlot_scCustom(CAR, features = c("Tnfrsf11a"), colors_use = viridis_plasma_dark_high)

FeaturePlot_scCustom(recluster, features = "Pth1r", colors_use = viridis_plasma_dark_high)

VlnPlot_scCustom(recluster, features = c("Enpep", "Thy1"), ggplot_default_colors = T, pt.size = 0, add.noise = F)

Stacked_VlnPlot(recluster, features = c("tdTomato", "Prg4","Col2a1", "Col10a1", "Acan", "Clec3b", "Mfap5"), ggplot_default_colors = T, pt.size = 0, add.noise = F, x_lab_rotate = T)

#pth target -- "Tnfsf11", "Mmp13", "Vdr", "Cebpb", "Igf1", "Il6", "Nfil3", "Atf4", "Plaur", "Nr4a2", "Nr4a3", "Fosl1", "Fosl2", "Crem", "Rgs2", "Pde4d", "Pde7b", "Pde10a"
# -- "Pparg", "Cebpa", "Zfp423", "Adipoq", "Apoe", "Nfia","Stat5a", "Stat5b","Insig2", "Fto"
# -- "Tgfb1", "Tgfbr1", "Smad3", "Acvr1b", "Ltbp1","Cebpb", "Inhba","Bmp4", "Bmp5", "Bmp6"
# -- "Wnt4", "Wnt5a", "Wnt5b", "Sfrp4", "Sfrp2", "Fzd1", "Fzd3", "Fzd5", "Lrp4","Lrp5","Lrp6", "Notum", "Nkd2", "Ctnnb1"
# -- "Notch1","Notch2","Notch3","Jag1", "Hes1", "Heyl"
# -- "Cxcl14", "Kitl", "Ngf", "Angpt1", "Vegfa","Vegfc", "Thbd", "Pdgfra", "Pdgfrb", "Vcam1", "Vav3"
# -- "Kitl", "Angpt1", "Vcam1", "Vegfa", "Vegfc","Pdgfra","Pdgfrb","Cxcl14","Vav3", "Thbd"
# -- "Ebf3", "Ebf1","Foxc1","Runx1", "Runx2"
# -- "Cebpb", "Tead1", "Tead4", "Hdac4", "Cthrc1", "Col3a1", "Col4a1", "Col4a2","Col6a1", "Col6a2", "Col6a3", "Col18a1","Col24a1"
# -- "Tnfsf11", "Mmp13", "Vdr", "Cebpb", "Igf1", "Il6", "Atf4", "Nr4a2", "Nr4a3", "Pde10a"
# -- "Dlx5","Runx2","Sp7","Bglap","Pth1r","Cebpb","Alpl","Wif1","Spp1","Col1a1","Col3a1", "Col4a1", "Col6a1", "Col18a1", "Col24a1", "Tead1", "Hdac4"
# -- "Dlx5", "Sp7", "Bglap", "Pth1r", "Col1a1", "Col3a1", "Col4a1", "Col6a1", "Col18a1", "Col24a1"
# -- "Casp9", "Dffb", "Casp8", "Bcl2l11", "Cradd", "Casp6", "Bad", "Irf1", "Cflar", "Bid", "Bok", "Bcl2l1"

selected_genes <- c("Epo")

Stacked_VlnPlot(CAR, features = selected_genes, ggplot_default_colors = F, pt.size = 0, plot_legend = T, x_lab_rotate = TRUE, add.noise = F, colors_use = deg_colors, split.by = "stim")

Stacked_VlnPlot(CAR.2d, features = selected_genes, ggplot_default_colors = F, pt.size = 0, plot_legend = T, x_lab_rotate = TRUE, add.noise = F, colors_use = deg_colors, split.by = "stim") | Stacked_VlnPlot(CAR.10d, features = selected_genes, ggplot_default_colors = F, pt.size = 0, plot_legend = T, x_lab_rotate = TRUE, add.noise = F, colors_use = deg_colors, split.by = "stim")


# CAR DEGS lists
# CAR1 -- "Chchd2", "Tmem59", "Arpc1b", "Naca", "Lrpap1", "Ywhab", "Cd81", "Oaz1", "Srpr", "C1ra", "Septin7", "Btf3", "Rarres2", "Dynll1", "Nfe2l1", "P4hb", "Set", "Os9", "Pcbp1", "Sar1a"
# CAR2 -- "Mllt10", "Kirrel", "Dock1", "Phactr2", "Rbm6", "Nfib", "Ccny", "Akap13", "Trio", "Arhgap35", "Snd1", "Kat6b", "Atp9b", "Med13l", "Trim2", "Tnrc6b", "Zfp638", "Ctif", "Prkce", "Mdfic"
# CAR3 -- "Clec2d", "Snai2", "Col6a5", "Sned1", "Ednra", "Cbln1", "Ccn5", "Abcg2", "Prelp", "Ttc14", "Fam13a", "Trim2", "Fmo2", "Esm1", "Foxc1", "Arhgap29", "Pnisr", "Plac9", "Mxi1", "Soat1"
# CAR4 -- "Dnah12", "Pecam1", "Sema4a", "Cx3cr1", "Ttc23l", "Nup107", "Ppm1n", "Nup210l", "Dnaja3", "B230369F24Rik", "Agbl3", "Gm47507", "Gm56599", "Bphl", "Nfs1", "Ccdc138", "Lsmem1", "Fbf1", "Vasp", "Plekhg2"
# CAR5 -- "Ccl2", "Cxcl1", "Vasn", "Niban2", "Sgk1", "Ier5", "Myc", "Mt2", "Tnfaip3", "Gm16685", "Mt1", "Ngf", "Klf4", "Cebpb", "Nfkbia", "Gfpt2", "Trib1", "Rcl1", "Dnajb1", "Rnf149"
# CAR6 -- "Kcnk2", "Alpl", "Wif1", "Slc20a2", "Ndnf", "Tenm4", "Limch1", "Lrp4", "Gm20619", "Spp1", "Nav2", "Adamts5", "Ncam1", "Nhsl1", "Smad6", "Kcnh1", "Col8a1", "Olfml3", "Ak5", "Tmtc2"

# Cluster DEGS -- "Cxcl12", "Lepr", "Ebf3", "Adipoq", "Foxc1", "Limch1", "Kcnk2", "Wif1", "Runx2", "Sp7", "Dlx5", "Lrp4","Spp1", "Alpl" ,"Col1a1", "Bglap", "Pth1r", "Prg4", "Mfap5", "Clec3b","Tnmd", "Fmod", "Tnfrsf11b", "Rgs5", "Myh11", "Kcnj8", "Pln", "Cdh5", "Pecam1", "Pax7", "Myf5"

list_genes <- c("Cxcl12", "Lepr", "Ebf3", "Adipoq", "Foxc1", "Limch1", "Kcnk2", "Wif1", "Runx2", "Sp7", "Dlx5", "Lrp4","Spp1", "Alpl" ,"Col1a1", "Bglap", "Dmp1","Pth1r", "Prg4", "Dcn", "Mgp","Tnmd", "Fmod", "Tnfrsf11b", "Rgs5", "Myh11", "Cdh5", "Pecam1", "Emcn")
plot6 <- DotPlot_scCustom(recluster, features = list_genes, x_lab_rotate=T, colors_use = viridis_plasma_dark_high) + scale_y_discrete(limits = rev) + scale_size(range = c(4,8))
plot1+plot2+plot3+plot4+plot5+plot6 + scale_y_discrete(limits = rev)

selected_genes <- c("Casp9", "Casp8", "Casp6", "Dffb", "Bid", "Bcl2l11", "Bad", "Bok", "Cradd", "Irf1", "Cflar")
DotPlot(CAR, features = selected_genes, dot.scale = 8, cols= "RdBu", split.by = "stim") + scale_size(range = c(4,8)) + coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_x_discrete(limits = rev)



#PTH receptor expressed in cluster 1, 3, 4, and 7 - look at response of PTH in these first, others will be secondary

# getting markers for reclusters
Idents(recluster) = "cluster_ids"
new.markers = FindAllMarkers(recluster, only.pos = T, min.pct = 0.2, min.diff.pct = 0.2, logfc.threshold = 0.25)

new.markers = FindAllMarkers(recluster, only.pos = T, min.pct = 0.4, min.diff.pct = 0.2, logfc.threshold = 0.5, test.use = "MAST")

# getting CAR markers

CAR <- subset(recluster, idents = c("CAR1","CAR2", "CAR3", "CAR4", "CAR5", "CAR6", "Pre-osteoblast", "Mature osteoblast"))
CARDEG <- subset(recluster, idents = c("CAR1","CAR2", "CAR3", "CAR4", "CAR5", "CAR6"))
ecpc <- subset(recluster, idents = c("Pericyte1", "Pericyte2", "Endothelial"))
sub_ob <- subset(merged_seurat_harmony, idents = c("4"))


recluster <- subset(recluster, idents = c("0", "1", "2", "3", "4", "5", "6", "7_0", "7_1", "7_2", "8", "9"))

Idents(peri) = "cluster_ids"
peri.2d <- subset(peri, idents = "2d")
peri.10d <- subset(peri, idents = "10d")

## "Pdzrn4", "Hp", "Cxcl12", "Lamp1", "Spacdr", "Vdr", "Cebpb", "Cd34", "Phex", "Angpt1", "Pparg", "Kitl", "Adipoq", "Pdgfd", "Ebf3", "Tnc", "Lpl", "Bglap", "Pth1r", "Cdh11", "Pde10a", "Pde3a", "Pde7b", "Ibsp", "Grem1", "Shroom3", "Bmp6", "Plxna2", "Cdh2", "Lepr"
selected_genes <- c("Cxcl12", "Kitl", "Plxna2", "Phex", "Pde10a", "Pde3a", "Pde7b", "Cd34")
Idents(peri.2d) = "cluster_ids"
Idents(peri.10d) = "cluster_ids"

levels(peri.10d) <- c("Periosteal1", "Periosteal2", "Periosteal3")
levels(recluster) <- c("0","1","2","3","4","5","6","7_0","7_1","7_2","8","9", "10")
levels(merged_seurat_harmony) <- c("0","1","2","3","4","5","6","7","8","9", "10")

VlnPlot(peri, features = c("Pdzrn4"), split.by = "stim", pt.size = 0)
Stacked_VlnPlot(peri.2d, features = selected_genes, ggplot_default_colors = T, pt.size = 0, plot_legend = T, x_lab_rotate = TRUE, split.by = "stim", add.noise = F) | Stacked_VlnPlot(peri.10d, features = selected_genes, ggplot_default_colors = T, pt.size = 0, plot_legend = T, x_lab_rotate = TRUE, split.by = "stim", add.noise = F)
FeaturePlot_scCustom(peri.2d, features = c("Cd34"), split.by = "stim", colors_use = viridis_plasma_dark_high)
FeaturePlot_scCustom(peri.10d, features = c("Cd34"), split.by = "stim", colors_use = viridis_plasma_dark_high)

Idents(CAR) = "time"
CAR.2d <- subset(CAR, idents = "2d")
CAR.10d <- subset(CAR, idents = "10d")

markers = FindAllMarkers(CARDEG, only.pos = T, min.pct = 0.4, min.diff.pct = 0.2, logfc.threshold = 0.5, test.use = "MAST")
markers.clusters = markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)

write.csv(markers.clusters, file = "clustermarkerswithCARonly.csv")


#sorting the markers
heatmap.markers = new.markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC) #### can also use slice_max(n = 5, order_by = avg_log2FC)
write.csv(heatmap.markers, file = "clusterspecificmarkers.csv")

heatmap.clusters <- read_csv("heatmapclustermarkerswithCARonly.csv")

test <- ScaleData(CARDEG, features = rownames(CARDEG)) ### need to scale the data for all features to display heatmap properly

maxcells <- min(table(Idents(test))) ## figure out what is the minimum cells for each cluster
DoHeatmap(subset(test, downsample = maxcells), features = heatmap.clusters$gene, disp.max = 2, size = 4) + guides(colour=F) + theme(axis.text.y = element_text(size = 10), text = element_text(size = 15,  angle = 0))



# compare cluster markers = will do this between car and osteoblast clusters - whats different between CAR

Idents(recluster) = "cluster_ids"
# start with only.pos = F to see what is different between genes and then you can adjust to only see which is only positive
test.markers <- FindMarkers(sub_ob, ident.1 = "Pre-osteoblast_0", ident.2 = c("Pre-osteoblast_1"), min.pct = 0.2, logfc.threshold = 0.5, only.pos = T)
head(test.markers, n = 20)
write.csv(test.markers, file = "AdipoCARvsOsteoCAR4volcano.csv")

myo.de.10d <- FindMarkers(myo.10d, ident.1 = "PTH", ident.2 = "Control", logfc.threshold = 0, min.pct = 0, test.use = "MAST")


# cluster re-order in ascending numeric order
levels(recluster) <- c("1", "2_0", "2_1", "10_0", "10_1", "10_2", "11", "13", "14_0", "14_1", "15")
levels(merged_seurat_harmony2) <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
levels(recluster) <- c("OsteoCAR1", "OsteoCAR2", "OsteoCAR3", "OsteoCAR4","AdipoCAR", "Pre-osteoblast_0", "Pre-osteoblast_1","Mature osteoblast", "Periosteal1", "Periosteal2", "Periosteal3", "Endothelial", "Pericyte1", "Pericyte2", "Myogenic")

#reorder split.by plots
recluster$orig.ident <- factor(x=recluster$orig.ident, levels = c("vehicle_2d", "pth_2d", "vehicle_10d", "pth_10d"))

levels(osteo.10d) <- c("OsteoCAR1", "OsteoCAR2", "AdipoCAR", "Osteoblast1", "Osteoblast2")

# reodering clusters to be compatible with dittoplot - cluster_order
cluster_order <- match(levels(merged_seurat_harmony2@meta.data[["RNA_snn_res.0.2"]]), metaLevels("RNA_snn_res.0.2", merged_seurat_harmony2))



# how many cells are in each cluster and per time point
Idents(recluster) = "new_clusters"
table(Idents(recluster))

Idents(recluster) = "time"
Idents(recluster) = "stim"


# looking at DEGs within each clusters across conditions =============================================================
# current idents 1 - EPC, 2_0 - CAR1, 2_1 - CAR2, 10_0 - OB1, 10_1 - OB2, 10_2 - OB3, 11 - CAR3, 13 - MSC progenitor, 14_0 - EC1, 14_1 - EC2, 15 - EC3
Idents(recluster) = "cluster_ids"


car1 <- subset(recluster, idents = "CAR1") #pth receptor
car2 <- subset(recluster, idents = "CAR2") #pth receptor
car3 <- subset(recluster, idents = "CAR3") #some cells have pthr receptor
car4 <- subset(recluster, idents = "CAR4") #some cells have pthr receptor
car5 <- subset(recluster, idents = "CAR5") #some cells have pthr receptor
car6 <- subset(recluster, idents = "CAR6")

ob1 <- subset(recluster, idents = "Pre-osteoblast") #pth receptor
ob2 <- subset(recluster, idents = "Mature osteoblast")

per1 <- subset(recluster, idents = "Periosteal1") #
per2 <- subset(recluster, idents = "Periosteal2") #
per3 <- subset(recluster, idents = "Periosteal3")
ec <- subset(recluster, idents = "Endothelial")
pc1 <- subset(recluster, idents = "Pericyte1")
pc2 <- subset(recluster, idents = "Pericyte2")
myo <- subset(recluster, idents = "Myogenic")

# DEGS for 2d vs 2m ---- 2m done already ---- starting 2d
Idents(car1) = "time" #
Idents(car2) = "time" #
Idents(car3) = "time" #
Idents(car4) = "time" #
Idents(car5) = "time" #
Idents(car6) = "time"
Idents(ob1) = "time" #
Idents(ob2) = "time" #
Idents(per1) = "time" #
Idents(per2) = "time"
Idents(per3) = "time"
Idents(ec) = "time"
Idents(pc1) = "time"
Idents(pc2) = "time"
Idents(myo) = "time"




# subsetting clusters with only 2d - osteoblasts first
car2.2d <- subset(car2, idents = "2d")
DimPlot(car2.2d, split.by = "stim", group.by = "Condition")

table(ob1.2d@meta.data$Condition)
table(ob1.2d@meta.data$stim)

stim_cells = Cells(epc.2m)[which(epc.2m$stim == "PTH")]
ctrl_cells = Cells(epc.2m)[which(epc.2m$stim == "Control")]

# for 2 day analysis only --- mvch2d is actually fvch2d
female_ctrl_cells = Cells(epc.2d)[which(epc.2d$Condition == "mvch2d")]

# use later ctrl_cells = Cells(ob.2d)[which(ob.2d$stim == "Control")]
slct_stim_cells = sample(stim_cells, size = 604, replace = F)
slct_ctrl_cells = sample(ctrl_cells, size = 604, replace = F)

# use for 2 day analysis
slct_ctrl_cells = sample(female_ctrl_cells, size = 156, replace = F)

# continue for analysis
subset_epc.2m = subset(epc.2m, cells = c(slct_stim_cells, slct_ctrl_cells))

table(subset_epc.2m@meta.data$Condition)
table(subset_epc.2m@meta.data$stim)

DimPlot(subset_epc.2m, group.by = "Condition", split.by = "stim")
DimPlot(subset_epc.2m, group.by = "stim")


Idents(car2.2d) = "stim"

### no thresholding or mini.pct for GSEA analysis - needs background genes
pc2.de <- FindMarkers(car2.2d, ident.1 = "PTH", ident.2 = "Control", logfc.threshold = 0, min.pct = 0, test.use = "MAST")
head(pc2.de, n=1000)
write.csv(pc2.de, file = "2dcar2PTHDEG.csv")

### for over-representation analysis - you should have thresholding
car2.de.2m <- FindMarkers(subset_car1.2m, ident.1 = "PTH", ident.2 = "Control", logfc.threshold = 0.5 , min.pct = 0.5, only.pos = F)
head(car1.de.2m, n=10)
write.csv(car1.de.2m, file = "2mcar1degsORA.csv")


FeaturePlot_scCustom(subset_car1.2m, features = c("Gli3"), max.cutoff = "4", split.by = "stim", colors_use = viridis_plasma_dark_high)
VlnPlot(subset_car1.2m, features = c("Lrp6"), split.by = "stim", pt.size = 0)


#combining dataframes
df1 <- read.csv("2d DEG analysis/2dcar1DEG.csv")
df1$cluster <- "OsteoCAR1"
df2 <- read.csv("2d DEG analysis/2dcar2DEG.csv")
df2$cluster <- "OsteoCAR2"
df3 <- read.csv("2d DEG analysis/2dcar3DEG.csv")
df3$cluster <- "OsteoCAR3"
df4 <- read.csv("2d DEG analysis/2dcar4DEG.csv")
df4$cluster <- "OsteoCAR4"
df5 <- read.csv("2d DEG analysis/2dcar5DEG.csv")
df5$cluster <- "AdipoCAR"
df6 <- read.csv("2d DEG analysis/2dob1DEG.csv")
df6$cluster <- "Pre-osteoblast"
df7 <- read.csv("2d DEG analysis/2dob2DEG.csv")
df7$cluster <- "Mature osteoblast"
df8 <- read.csv("2d DEG analysis/2dper1DEG.csv")
df8$cluster <- "Periosteal1"
df9 <- read.csv("2d DEG analysis/2dper2DEG.csv")
df9$cluster <- "Periosteal2"
df10 <- read.csv("2d DEG analysis/2dper3DEG.csv")
df10$cluster <- "Periosteal3"
df11 <- read.csv("2d DEG analysis/2decDEG.csv")
df11$cluster <- "Endothelial"
df12 <- read.csv("2d DEG analysis/2dpc1DEG.csv")
df12$cluster <- "Pericyte1"
df13 <- read.csv("2d DEG analysis/2dpc2DEG.csv")
df13$cluster <- "Pericyte2"
df14 <- read.csv("2d DEG analysis/2dmyoDEG.csv")
df14$cluster <- "Myogenic"

df15 <- read.csv("10d DEG analysis/10dcar1DEG.csv")
df15$cluster <- "OsteoCAR1"
df16 <- read.csv("10d DEG analysis/10dcar2DEG.csv")
df16$cluster <- "OsteoCAR2"
df17 <- read.csv("10d DEG analysis/10dcar3DEG.csv")
df17$cluster <- "OsteoCAR3"
df18 <- read.csv("10d DEG analysis/10dcar4DEG.csv")
df18$cluster <- "OsteoCAR4"
df19 <- read.csv("10d DEG analysis/10dcar5DEG.csv")
df19$cluster <- "AdipoCAR"
df20 <- read.csv("10d DEG analysis/10dob1DEG.csv")
df20$cluster <- "Pre-osteoblast"
df21 <- read.csv("10d DEG analysis/10dob2DEG.csv")
df21$cluster <- "Mature osteoblast"
df22 <- read.csv("10d DEG analysis/10dper1DEG.csv")
df22$cluster <- "Periosteal1"
df23 <- read.csv("10d DEG analysis/10dper2DEG.csv")
df23$cluster <- "Periosteal2"
df24 <- read.csv("10d DEG analysis/10dper3DEG.csv")
df24$cluster <- "Periosteal3"
df25 <- read.csv("10d DEG analysis/10decDEG.csv")
df25$cluster <- "Endothelial"
df26 <- read.csv("10d DEG analysis/10dpc1DEG.csv")
df26$cluster <- "Pericyte1"
df27 <- read.csv("10d DEG analysis/10dpc2DEG.csv")
df27$cluster <- "Pericyte2"
df28 <- read.csv("10d DEG analysis/10dmyoDEG.csv")
df28$cluster <- "Myogenic"

# apply 2 days to all the df
df_list <- list(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14)

df_list <- lapply(df_list, function(df) {
  df$time <- "2d"
  return(df)
})

list2env(setNames(df_list, paste0("df", 1:14)), envir = .GlobalEnv)

# apply 10 days to all the df
df_list <- list(df15, df16, df17, df18, df19, df20, df21, df22, df23, df24, df25, df26, df27, df28)

df_list <- lapply(df_list, function(df) {
  df$time <- "10d"
  return(df)
})

list2env(setNames(df_list, paste0("df", 15:28)), envir = .GlobalEnv)

# Create a list of your data frames
df_list <- list(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14, df15, df16, df17, df18, df19, df20, df21, df22, df23, df24, df25, df26, df27, df28)


# Apply filtering for significant genes
df_list_filtered <- lapply(df_list, function(df) {
  df %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
})

# Apply z-score scaling to avg_log2FC for each data frame
# Function to calculate the z-score for avg_log2FC
calc_zscore <- function(df) {
  df$z_score <- (df$avg_log2FC - mean(df$avg_log2FC, na.rm = TRUE)) / sd(df$avg_log2FC, na.rm = TRUE)
  return(df)
}

df_list_filtered_scaled <- lapply(df_list_filtered, calc_zscore)

#combining the df

combined_df <- bind_rows(df_list_filtered_scaled)

colnames(combined_df)[colnames(combined_df) == 'X'] <- 'gene'

write.csv(combined_df, file = "clusterDEGsummarized.csv")

combined_df$cluster <- factor(combined_df$cluster, levels = c("OsteoCAR1", "OsteoCAR2", "OsteoCAR3", "OsteoCAR4","AdipoCAR", "Pre-osteoblast", "Mature osteoblast", "Periosteal1", "Periosteal2", "Periosteal3", "Endothelial", "Pericyte1", "Pericyte2", "Myogenic"))
combined_df$cluster <- fct_rev(combined_df$cluster)
combined_df$time <- factor(combined_df$time, levels = c("2d", "10d"))
ggplot(combined_df, aes(x=gene, y=cluster, fill=z_score)) + geom_tile() + 
  facet_grid(~time, labeller = labeller(time = c("2d" = "2 days DEGs", "10d" = "10 days DEGs"))) + scale_fill_gradient2(high = "#FF0000", mid = "#FFFFFF", low = "#0000FF", midpoint = 0, limits = c(-3, 3), oob = scales::squish) +
  guides(fill = guide_colourbar(title = "Z-score")) + theme(panel.border = element_rect(color = 'black', fill = NA, size = 1)) + theme(axis.text.y = element_text(face="bold", size = 10)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank())


# finding out number of genes that are upregulated and downregulated

nrow(df1[df1$avg_log2FC > 0.5 & df1$p_val_adj < 0.05, ]) #upregulated
nrow(df1[df1$avg_log2FC < -0.5 & df1$p_val_adj < 0.05, ]) #downregulated



# filtering genes

filtered <- df[abs(df$avg_log2FC) > 0.5 & df$p_val_adj < 0.05,]

filtered <- combined_df %>%
  filter(cluster == "Pre-osteoblast", avg_log2FC < -0.5, p_val_adj < 0.05, time == "2d")

write.csv(filtered, file = "test.csv", row.names = F)

# figure out the number of DEGS that meet your conditions

nrow(df[df$avg_log2FC > 0.5 & df$p_val_adj < 0.05, ])
nrow(df[df$avg_log2FC < -0.5 & df$p_val_adj < 0.05, ])
nrow(df[abs(df$avg_log2FC) > 0.5 & df$p_val_adj < 0.05,])


# plotting volcanoplots with enhanced volcano
library(EnhancedVolcano)


# processing the DEG dataframe computed from FindMarkers()
df <- read.csv("volcanodeg/2dob1DEG.csv") #clusterdegds.csv is your list of differentially expressed genes from FindMarkers()

files <- list.files(path = "volcanodeg", pattern = "\\.csv$", full.names = TRUE)
df <- do.call(rbind, lapply(files, read.csv))

### better to use this
df <- df %>%
  group_by(X) %>%
  #filter(abs(avg_log2FC) > 0.5, p_val_adj < 0.05) %>%   # only keep rows with |avg_log2FC| > 0.5
  slice_min(order_by = if_else(p_val_adj == 0, Inf, p_val_adj),
            n = 1, with_ties = FALSE) %>%
  ungroup()

#write.csv(df, file = "combinedCAROBDEG.csv", row.names = T)



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
                xlim = c(-7,11),
                cutoffLineWidth = 1,
                labSize = 6,
                gridlines.minor = F, gridlines.major = F,
                pointSize = c(ifelse(abs(df$avg_log2FC) > 0.5 & df$p_val_adj < 0.05, 2, 2)),
                selectLab = c("Vdr", "Tnfsf11", "Mmp13","Pparg", "Cebpa", "Zfp467"),
                labFace = "bold")
                #selectLab = c("Hbb-bt", "Hba-a1", "B2m", "Cst3", "Rppl18", "Rps3a1", "Rps2", "Rpl27a", "Clec12a", "Pabpc4", "Gpx1", "Col3a1", "Col1a1", "Gli3", "Arhgap21", "Sulf1", "Itga8", "Lama3", "Col6a3", "Dbp", "Nrg3", "Ptn", "Stk39", "Ppp3ca", "Rnf150", "Ppard", "Pdgfd", "Tnfsf11", "Il1rn", "Sned1", "Zfp148", "Zfp148", "Ctsb", "Per3", "Hk2", "Il4ra", "Ctnnb1", "Il17ra"))
                
                #selectLab = c("Stat3", "Jak1", "Cebpb", "Per3", "Rplp0", "Dbp", "mt-Nd3", "Rps3a3", "Mme", "Fosb", "Irs2", "Hspa5", "Cebpb", "Socs3", "Casp4", "Adamts1", "Pim1", "Nr4a1"))

#                labFace = "bold",


#------- SingleR automated single cell annotation
# will probably want to use ImmGenData (mouse bulk RNAseq) and MouseRNAseqData

ref.data1 <- celldex::ImmGenData()
ref.data2 <- celldex::MouseRNAseqData()
seurat.counts <- GetAssayData(merged_seurat_harmony, assay = 'RNA', layer = 'counts')

singleR.annotate <- SingleR(seurat.counts, ref = ref.data1, labels = ref.data$label.main)
merged_seurat_harmony$singleR.labels2 <- singleR.annotate$labels[match(rownames(merged_seurat_harmony@meta.data), rownames(singleR.annotate))]


## ----- clusterprofiler setup-------------------
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(AnnotationHub)
library(ReactomePA)


#### read in data

df <- read.csv('2mcar1DEG.csv')
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
Idents(recluster) = "recluster_ids"
test <- subset(recluster, idents = "EPC")
Idents(test) = "recluster_ids"
FeatureScatter_scCustom(test, feature1 = "Ebf3", feature2 = "Pecam1", pt.size = 1.5, colors_use = "dodgerblue")

(sum(FetchData(test, vars = "Pecam1") > 0.5 & FetchData(test, vars = "Ebf3") > 0.5)) / ncol(test)

Percent_Expressing(CARDEG, features = c("Ebf3"), split_by = c("stim"))


# cell cycle analysis ==========================

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

cellcycle <- CellCycleScoring(recluster.2d.veh, s.features = s.genes, g2m.features = g2m.genes)
Idents(cellcycle) = "new_clusters"
Idents(cellcycle) = "Phase"
write.csv(table(Idents(cellcycle), cellcycle$Phase), file = "cellcycle.csv")

Idents(recluster) = "time"
recluster.2d <- subset(recluster, idents = "2d")
recluster.10d <- subset(recluster, idents = "10d")

Idents(recluster.2d) = "stim"
Idents(recluster.10d) = "stim"

recluster.10d.veh <- subset(recluster.10d, idents = "Control")
recluster.10d.pth <- subset(recluster.10d, idents = "PTH")

Idents(recluster.10d.veh) = "cluster_ids"
Idents(recluster.10d.pth) = "cluster_ids"

cellcycle <- CellCycleScoring(recluster.10d.pth, s.features = s.genes, g2m.features = g2m.genes)
write.csv(table(Idents(cellcycle), cellcycle$Phase), file = "10dpthcellcycle.csv")

#if you want to look at the number of cells in each cluster per time
table(Idents(recluster), recluster$time)
Idents(recluster) = "cluster_ids"
test <- subset(recluster, idents = "2m")
Idents(test) = "recluster_ids"
table(Idents(test), test$stim)


# renaming cluster IDs after reclustering
Idents(recluster) = "new_clusters"
levels(recluster) <- c("2_0", "2_1", "11", "10_0", "10_1", "10_2", "13", "1" ,"14_0", "14_1", "15")

recluster <- RenameIdents(recluster,'2_0' = "osteoCAR1", '2_1' = "osteoCAR2", '11' = "adipoCAR", '10_0' = "osteoblast1", '10_1' = "osteoblast2", '10_2' = "immature osteocyte", '13' = "MSC progenitor", '1' = "EPC", '14_0' = "EC1", '14_1' = "EC2", '15' = "EC3")
recluster$recluster_ids <- recluster@active.ident
write.csv(recluster@meta.data, "recluster_cell_annotation_metadata.csv")


# umap can label in better details - first assign the variable for DimPlot then apply it into LabelClusters()
p <- DimPlot_scCustom(test, pt.size = 0.1, label = F)
LabelClusters(p, id = "ident",  fontface = "bold", color = "black", repel = T)




# trying monocle3 cd

cds <- SeuratWrappers::as.cell_data_set(recluster) #change to cds

cds <- cluster_cells(cds, reduction_method = "UMAP")

plot_cells(cds, show_trajectory_graph = F, color_cells_by = "partition")

cds <- learn_graph(cds, use_partition = T)

cds <- order_cells(cds)

plot_cells(cds, color_cells_by = "pseudotime", label_branch_points = F, label_leaves = F, graph_label_size = 0, show_trajectory_graph = T)


# dittoplotting

Idents(merged_seurat_harmony2) = "time"
test <- subset(merged_seurat_harmony2, idents = "2d")
dittoBarPlot(test, "RNA_snn_res.0.2", group.by = "stim", var.labels.reorder = cluster_order)


# to convert to anndata, first change RNA5 to RNA3 (seurat v5 feature!)
recluster1[["RNA3"]] <- as(recluster1[["RNA"]], "Assay")
DefaultAssay(recluster1) <- "RNA3"
recluster1[["RNA"]] <- NULL
recluster1 <- RenameAssays(recluster1, RNA3='RNA')


#convert seurat object to anndata
SaveH5Seurat(recluster1, filename = "ebf3PTH.h5Seurat")
Convert("ebf3PTH.h5Seurat", dest = "h5ad")


Convert("ebf3PTH.h5ad", dest = "h5Seurat", overwrite = T)
test <- LoadH5Seurat("ebf3PTH.h5Seurat")
View(test)

list_genes <- c("tdTomato", "Cxcl12", "Lepr", "Ebf3")
DotPlot_scCustom(recluster, features = list_genes, cols = viridis_plasma_dark_high, dot.scale = 8) + 
  guides(color = guide_colorbar(title = 'Scaled Average Expression')) + 
  RotatedAxis()

"Tnfsf11", "Vdr", "Mmp13", "If1", "Il6", "Pde4d", "Adcy2", "Crem"
"Ebf3", "Ebf1", "Foxc1", "Runx1", "Runx2"
"Kitl", "Cxcl14", "Ngf", "Angpt1", "Vegfa", "Vegfc", "Thbd"

selected_genes <- c("Snai1", "Snai2", "Twist1", "Zeb1", "Zeb2")
Idents(CAR.2d) = "cluster_ids"
Idents(CAR.10d) = "cluster_ids"

levels(CAR) <- c("CAR1","CAR2", "CAR3","CAR4", "CAR5","CAR6","Pre-osteoblast", "Mature osteoblast")
levels(recluster) <- c("CAR1","CAR2", "CAR3","CAR4", "CAR5","CAR6","Pre-osteoblast", "Mature osteoblast","Periosteal1", "Periosteal2", "Periosteal3", "Pericyte1", "Pericyte2", "Endothelial", "Myogenic")

recluster$recluster_ids <- factor(
  recluster$recluster_ids,
  levels = c(
    "CAR1","CAR2","CAR3","CAR4","CAR5","CAR6",
    "Pre-osteoblast","Mature osteoblast",
    "Periosteal1","Periosteal2","Periosteal3",
    "Pericyte1","Pericyte2",
    "Endothelial","Myogenic"
  )
)


Idents(recluster) <- factor(Idents(recluster), levels = rev(levels(Idents(recluster))))
DotPlot_scCustom(recluster, features = c("Cxcl12", "Lepr", "Ebf3", "Adipoq", "Foxc1", "Limch1", "Kcnk2", "Wif1", "Runx2", "Sp7", "Spp1", "Alpl", "Col1a1", "Bglap", "Pth1r", "Prg4", "Mfap5", "Clec3b","Tnmd", "Fmod", "Tnfrsf11b", "Rgs5", "Myh11", "Kcnj8", "Pln", "Cdh5", "Pecam1"), cols = viridis_plasma_dark_high, dot.scale = 8)

Stacked_VlnPlot(CAR.2d, features = selected_genes, ggplot_default_colors = T, pt.size = 0, plot_legend = T, x_lab_rotate = TRUE, split.by = "stim", add.noise = F) | Stacked_VlnPlot(CAR.10d, features = selected_genes, ggplot_default_colors = T, pt.size = 0, plot_legend = T, x_lab_rotate = TRUE, split.by = "stim", add.noise = F)

Stacked_VlnPlot(merged_seurat_harmony, features = c("Nr4a2", "Pde3a", "Pde10a"),  ggplot_default_colors = F, pt.size = 0, plot_legend = T, x_lab_rotate = TRUE, add.noise = F, colors_use = deg_colors, split.by = "stim")


#panel of genes
# adipogenic - Pparg Cebpa Adipoq Lpin1 Esr1 Ebf1
# osteogenic - Runx2 Sp7 Alpl Pth1r Col1a1 Bglap Spp1 
# osteoclastic - Tnfsf11 Csf1 Cthrc1 S1p Sema4d Ctf1 Pdgfb Atp6v0d2
# PTH signal - Mef2c Efnb2 Igf1 Fgf2 Tnfsf10 Vdr Dgat1 Mmp13 PKA Hdac4 Hdac5 Tnfrsf11b Zfp467 Cav1
# Wnt signal - Lrp4 Lrp5 Lrp6 Lrp8 Fzd1 Dkk1 Wif1 Axin2 Wnt5b Wnt4 Wnt2b Wnt11 Wnt5a
# CAR genes - Cxcl12 Kitl Lepr Foxc1, Pdgfrb Ebf1 Ebf3 Il7
# T lymphocyte and macrophage - Cd40lg Wnt10b
# Inflammation Tnfa Il17a Nlrp3 Flt3 Il7 Il1b Il6 Thpo
# calcium channels - Cacna1b Cacna1c Cacna1d Cacna1h Cacna1g Cacna1l
# Bmps - Bmp1 Bmp2 Bmp3 Bmp4 Bmp5 Bmp6
# Plexin B1/Sema4d pathway - Plxnb1 Sema4d Vegfa Akt1 Rhod
# Igf pathway - Igf1 Igf1r Irs1 Fos Jun
# recruitment of CAR cells - Il1 Pdgf Vegf Bmp Cxcl12/Cxcr4/Ccl7 Tgfb
# Fgfs - Fgf1 Fgf2 Fgf7 Fgf10 Fgf11 Fgf13 Fgf18 Fgfr1 Fgfr2 Fgfr3
# Cxcls/Ccl/Cxcr - Ccl2 Ccl9 Ccl19 Ccl11 Ccl7 Cxcl1 Cxcl5 Cxcl9 Cxcl10 Cxcl14
# Interleukins - Il1rn Il15 Il34 Il20ra Il4ra Il17rd Il16 Il31ra Il15ra Il1r1 Il17ra Il
# fibrotic panel of CAR - col1a1 col3a1 col6a3 gli1 pdpn pdgfra


# testing how to do apoptosis scoring

library(UCell)

test <- recluster
RE_AP <- read.gmt('~/Documents/Ebf3_PTH_scRNA_analysis/gmt/REACTOME_APOPTOSIS.v2024.1.Mm.gmt')
AP_genes <- c(RE_AP$gene)
AP_genes <- stringr::str_to_title(AP_genes)
AP_genes <- list(AP_genes[AP_genes%in%rownames(test)])
test <- AddModuleScore_UCell(test, features = AP_genes, name = 'AP.Score')


levels(test) <- c("OsteoCAR1", "OsteoCAR2", "OsteoCAR3", "OsteoCAR4","AdipoCAR", "Pre-osteoblast", "Mature osteoblast", "Periosteal1", "Periosteal2", "Periosteal3", "Endothelial", "Pericyte1", "Pericyte2", "Myogenic")


Idents(test) = "cluster_ids"

FeaturePlot(recluster, features = c("Pth1r"), label = T) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

RidgePlot(CARDEG, features = 'signature_1AP.Score', group.by = "cluster_ids")

VlnPlot(test, features = "signature_1AP.Score", pt.size = 0)

test1 <- subset(test, idents = "Mature osteoblast")

VlnPlot(test1, feature ="signature_1AP.Score", split.by = "stim", pt.size =)

Idents(test1) = "stim"

ctrl_scores <- test1$signature_1AP.Score[test1$stim == "Control"]
pth_scores <- test1$signature_1AP.Score[test1$stim == "PTH"]
t_test_results <- t.test(pth_scores, ctrl_scores)
print(t_test_results)



# trying scVelo ===============================================

# to convert to anndata, first change RNA5 to RNA3 (seurat v5 feature!)
library(sceasy)

use_condaenv('/opt/miniconda3/')
loompy <- reticulate::import('loompy')

test[["RNA3"]] <- as(test[["RNA"]], "Assay")
DefaultAssay(test) <- "RNA3"
test[["RNA"]] <- NULL
test <- RenameAssays(test, RNA3='RNA')


#convert seurat object to h5ad

library(reticulate)

# save your seurat object as test so you dont affect your main analysis
test <- recluster
recluster <- recluster1

recluster.old <- recluster

Idents(test) = "stim"
test1 <- subset(test, ident = "Control")
test2 <- subset(test, ident = "PTH")
Idents(test2) = "cluster_ids"

# this is to convert from seurat to h5ad
as.anndata(test2, file_path = "~/Documents/Ebf3_PTH_scRNA_analysis/", file_name = "pth.h5ad")


# Store all dataframes in a list
df_list <- list(df1 = df15, df2 = df16, df3 = df17, df4 = df18, df5 = df19)

# filtered
df_list <- lapply(df_list, function(df) {
  df %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
})

# Create a unique gene list from all dataframes
all_genes <- unique(unlist(map(df_list, ~ .x$X)))

# Create a presence/absence matrix
presence_matrix <- data.frame(Gene = all_genes)

for (i in seq_along(df_list)) {
  df_name <- names(df_list)[i]
  presence_matrix[[df_name]] <- presence_matrix$Gene %in% df_list[[i]]$X
}

# Convert TRUE/FALSE to 1/0 for better readability
presence_matrix[,-1] <- lapply(presence_matrix[,-1], as.integer)

# Print the presence/absence matrix
print(presence_matrix)

# -------------------------------
# FIND UNIQUE GENES PER DATAFRAME
# -------------------------------

unique_genes_list <- list()

for (i in seq_along(df_list)) {
  df_name <- names(df_list)[i]
  others <- df_list[-i]  # All dataframes except the current one
  unique_genes_list[[df_name]] <- setdiff(df_list[[i]]$X, unlist(map(others, ~ .x$X)))
}

# Print unique genes in each dataframe
for (name in names(unique_genes_list)) {
  cat("\nUnique genes in", name, ":\n")
  print(unique_genes_list[[name]])
}

# -------------------------------
# EXPORT THE UNIQUE GENES IN LONG FORMAT WITH ADDITIONAL STATISTICS
# -------------------------------
unique_genes_long <- do.call(rbind, lapply(names(unique_genes_list), function(df_name) {
  # Get the current dataframe
  df_current <- df_list[[df_name]]
  # Unique genes for this dataframe
  unique_genes <- unique_genes_list[[df_name]]
  # Subset the dataframe to the unique genes and select desired columns:
  # Rename the X column to Gene
  df_subset <- df_current %>%
    filter(X %in% unique_genes) %>%
    select(Gene = X, avg_log2FC, p_val_adj)
  # Add a column for the source
  df_subset$Source <- df_name
  df_subset
}))

# Write the long-format table as a CSV file
write.csv(unique_genes_long, file = "2daysuniqueDEGs.csv", row.names = FALSE)


# -------------------------------
# FIND COMMON GENES (IN ALL 5 DF)
# -------------------------------

common_genes <- Reduce(intersect, map(df_list, ~ .x$X))
cat("\nGenes present in ALL 5 dataframes:\n")
print(common_genes)

combined_df <- do.call(rbind, lapply(names(df_list), function(df_name) {
  # Get the current dataframe
  df_current <- df_list[[df_name]]
  # Subset the dataframe to select desired columns:
  # Rename the X column to Gene
  df_subset <- df_current %>%
    filter(X %in% common_genes) %>%
    select(Gene = X, avg_log2FC, p_val_adj)
  # Add a column for the source
  df_subset$Source <- df_name
  df_subset
}))



# Write the long-format table as a CSV file
write.csv(combined_df, file = "10dayscommonDEGs.csv", row.names = FALSE)

# -------------------------------
# FIND GENES PRESENT IN AT LEAST 2 DATAFRAMES
# -------------------------------
gene_counts <- table(unlist(map(df_list, ~ .x$X)))
genes_in_2_or_more <- names(gene_counts[gene_counts >= 2])

cat("\nGenes present in at least 2 dataframes:\n")
print(genes_in_2_or_more)




# correlation matrix

normalized_data <- GetAssayData(CARDEG, assay = "RNA", layer = "data")
# Get cluster identities
clusters <- Idents(CARDEG)

# Helper function to compute average expression for each cluster
average_cluster_expression <- function(data, clusters) {
  cluster_ids <- unique(clusters)
  cluster_avg_list <- list()
  
  for (cluster in cluster_ids) {
    cells_in_cluster <- names(clusters[clusters == cluster])
    avg_expr <- rowMeans(data[, cells_in_cluster])
    cluster_avg_list[[as.character(cluster)]] <- avg_expr
  }
  return(cluster_avg_list)
}

# Calculate average expression for each cluster
cluster_avg_expression <- average_cluster_expression(normalized_data, clusters)

# Combine all cluster averages into one matrix
expr_matrix <- do.call(cbind, cluster_avg_expression)

# Set appropriate row and column names for the matrix
colnames(expr_matrix) <- names(cluster_avg_expression)
rownames(expr_matrix) <- rownames(normalized_data)

# Compute the Pearson correlation matrix
cor_matrix <- cor(expr_matrix, method = "pearson")

# Visualize the Pearson correlation matrix with values on the heatmap
library(pheatmap)

cluster_annotation <- data.frame(Cluster = factor(c("OsteoCAR1", "OsteoCAR2", "OsteoCAR3", "OsteoCAR4", "AdipoCAR","Pre-osteoblast"),
                                                  levels = c("OsteoCAR1", "OsteoCAR2", "OsteoCAR3", "OsteoCAR4", "AdipoCAR","Pre-osteoblast")))
pheatmap(cor_matrix,
         clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation",
         display_numbers = TRUE,   # Display correlation values
         number_format = "%.2f",   # Format of the displayed numbers (2 decimal places)
         fontsize_number = 10,     # Font size of the numbers
         color = colorRampPalette(c("blue", "white", "red"))(100),
         treeheight_row = 0,
         treeheight_col = 0)  # Color scale






library(cowplot)
library(ggplot2)

# -------------------------------
# CAR DEGs gene lists
# -------------------------------
CAR1_list <- c("Chchd2", "Tmem59", "Arpc1b", "Naca", "Lrpap1", "Ywhab", "Cd81", "Oaz1", "Srpr", "C1ra",
               "Septin7", "Btf3", "Rarres2", "Dynll1", "Nfe2l1", "P4hb", "Set", "Os9", "Pcbp1", "Sar1a")

CAR2_list <- c("Mllt10", "Kirrel", "Dock1", "Phactr2", "Rbm6", "Nfib", "Ccny", "Akap13", "Trio", "Arhgap35",
               "Snd1", "Kat6b", "Atp9b", "Med13l", "Trim2", "Tnrc6b", "Zfp638", "Ctif", "Prkce", "Mdfic")

CAR3_list <- c("Clec2d", "Snai2", "Col6a5", "Sned1", "Ednra", "Cbln1", "Ccn5", "Abcg2", "Prelp", "Ttc14",
               "Fam13a", "Trim2", "Fmo2", "Esm1", "Foxc1", "Arhgap29", "Pnisr", "Plac9", "Mxi1", "Soat1")

CAR4_list <- c("Dnah12", "Pecam1", "Sema4a", "Cx3cr1", "Ttc23l", "Nup107", "Ppm1n", "Nup210l", "Dnaja3",
               "B230369F24Rik", "Agbl3", "Gm47507", "Gm56599", "Bphl", "Nfs1", "Ccdc138", "Lsmem1", "Fbf1",
               "Vasp", "Plekhg2")

CAR5_list <- c("Ccl2", "Cxcl1", "Vasn", "Niban2", "Sgk1", "Ier5", "Myc", "Mt2", "Tnfaip3", "Gm16685",
               "Mt1", "Ngf", "Klf4", "Cebpb", "Nfkbia", "Gfpt2", "Trib1", "Rcl1", "Dnajb1", "Rnf149")

CAR6_list <- c("Kcnk2", "Alpl", "Wif1", "Slc20a2", "Ndnf", "Tenm4", "Limch1", "Lrp4", "Gm20619", "Spp1",
               "Nav2", "Adamts5", "Ncam1", "Nhsl1", "Smad6", "Kcnh1", "Col8a1", "Olfml3", "Ak5", "Tmtc2")

# -------------------------------
# Common color and dot size scale
# -------------------------------
dot_scale <- 6
size_range <- c(1, 10)

# -------------------------------
# Create dot plots with reversed y-axis
# -------------------------------
p1 <- DotPlot_scCustom(CARDEG, features = CAR1_list, dot.scale = dot_scale, x_lab_rotate=T, colors_use = viridis_plasma_dark_high) + 
  scale_size(range = size_range) + scale_y_discrete(limits = rev)

p2 <- DotPlot_scCustom(CARDEG, features = CAR2_list, dot.scale = dot_scale, x_lab_rotate=T, colors_use = viridis_plasma_dark_high) + 
  scale_size(range = size_range) + scale_y_discrete(limits = rev)

p3 <- DotPlot_scCustom(CARDEG, features = CAR3_list, dot.scale = dot_scale, x_lab_rotate=T, colors_use = viridis_plasma_dark_high) + 
  scale_size(range = size_range) + scale_y_discrete(limits = rev)

p4 <- DotPlot_scCustom(CARDEG, features = CAR4_list, dot.scale = dot_scale, x_lab_rotate=T, colors_use = viridis_plasma_dark_high) + 
  scale_size(range = size_range) + scale_y_discrete(limits = rev)

p5 <- DotPlot_scCustom(CARDEG, features = CAR5_list, dot.scale = dot_scale, x_lab_rotate=T, colors_use = viridis_plasma_dark_high) + 
  scale_size(range = size_range) + scale_y_discrete(limits = rev)

p6 <- DotPlot_scCustom(CARDEG, features = CAR6_list, dot.scale = dot_scale, x_lab_rotate=T, colors_use = viridis_plasma_dark_high) + 
  scale_size(range = size_range) + scale_y_discrete(limits = rev)

# -------------------------------
# Combine all plots in a 2x3 grid
# -------------------------------
plot_grid(p1, p2, p3, p4, p5, p6, ncol = 2, align = "hv")



test <- merged_seurat_harmony

merged_seurat_harmony <- subset(merged_seurat_harmony, subset = orig.ident != "female_pth_2d")

#==========================
# Mouse ortholog signatures
#==========================

CAR1.genes <- c(
  "Frzb",
  "Pdk4",
  "Mgst1",
  "Tagln",
  "Gpx3",
  "Abca8a",      # use Abca8 if that's the symbol in your dataset
  "Plaat4",
  "Snai2",
  "Fmo2",
  "Itgb2",
  "Tsc22d1",
  "Cyb5a",
  "Chmp5",
  "Brk1",
  "Cstb",
  "Aldoa",
  "U2af1",
  "Slc25a6",
  "Gadd45gip1",
  "Septin7"
)

CAR2.genes <- c(
  "Cd302",        # LY75-CD302
  "Mybpc1",
  "Hfm1",
  "Plcb4",
  "Larp1b",
  "Galnt17",
  "Rgs6",
  "Gcnt1",
  "Tex14",
  "Kazn",
  "Wwox",
  "Sik3",
  "Depp1",
  "Lrmda",
  "Auts2",
  "Nrxn3",
  "Asxl1",
  "Zfand3",
  "Cemip2",
  "Sypl1"
)

CAR3.genes <- c(
  "Ccn2",
  "Slc5a3",
  "Jun",
  "Id4",
  "Scara5",
  "Mrc2",
  "Igfbp4",
  "Mesd",
  "Thbs1",
  "Col1a1",
  "Sfrp1",
  "Col3a1",
  "Edil3",
  "Pltp",
  "Apoe",
  "Itgbl1",
  "Gas1",
  "Col8a1",
  "Col18a1",
  "Thy1"
)

CAR1.genes <- intersect(CAR1.genes, rownames(test))
CAR2.genes <- intersect(CAR2.genes, rownames(test))
CAR3.genes <- intersect(CAR3.genes, rownames(test))

test <- AddModuleScore(
  test,
  features = list(
    CAR1.genes,
    CAR2.genes,
    CAR3.genes
  ),
  name = c("HumanCAR1","HumanCAR2","HumanCAR3")
)

FeaturePlot(
  test,
  features = c("HumanCAR11","HumanCAR22","HumanCAR33"),
  min.cutoff = "q05",
  max.cutoff = "q95",
  ncol = 3
)
VlnPlot(
  test,
  features = c("HumanCAR33"),
  group.by = "recluster_ids",
  pt.size = 0
)


# ----------------------------
# PTH response module
# ----------------------------
pth_genes <- c(
  "Tnfsf11",
  "Fos",
  "Fosl1",
  "Junb",
  "Egr1",
  "Crem",
  "Nr4a2",
  "Nr4a3",
  "Nfkbia",
  "Il6",
  "Pde4d",
  "Igf1",
  "Pde10a"
)

tgfb_emt <- c(
  "Serpine1",   # PAI-1 – classic TGFβ EMT/migration target
  "Fn1",        # Fibronectin
  "Col1a1",
  "Postn",
  "Tagln",      # Actin remodeling
  "Timp1",      # EMT / myofibroblast-like shift
  "Vim",        # Mesenchymal intermediate filament
  "Ccn2",       # Matrix + migration
  "Itga5",      # Fibronectin receptor
  "Mmp2",
  "Mmp14",
  "Plod2",
  "Tgm2",
  "Rock2", ##
  "Ltbp1", ##
  "Loxl2"##
)

apoptosis_genes <- c(
  "Bcl2","Bcl2l1","Mcl1","Bcl2a1a","Bcl2a1b",
  "Cflar","Xiap","Birc2","Birc3",
  "Traf1","Traf2",
  "Akt1","Akt2",
  "Rel","Rela","Nfkb1","Nfkbia"
  )




recluster <- AddModuleScore(
  recluster,
  features = list(pth_genes),
  name = "PTH_Score",
  ctrl = 50
)

recluster <- AddModuleScore(
  recluster,
  features = list(tgfb_emt),
  name = "TGFB_EMT",
  ctrl = 50
)

recluster <- AddModuleScore(
  recluster,
  features = list(apoptosis_genes),
  name = "Apoptosis_Score",
  ctrl = 50
)


Stacked_VlnPlot(
  CAR,
  features = c("Survival_Score"),
  group.by = "recluster_ids",
  split.by = "stim",
  add.noise = F, colors_use = deg_colors,
  pt.size = 0
)

FeaturePlot_scCustom(
  CAR,
  features = c("Apoptosis_Score"),
  split.by = "stim",
  reduction = "umap",
  raster = F,
  label = T,
  repel = F,
  pt.size = 1
)

# module score effect size testing just between two comparisons======
library(dplyr)

df <- recluster@meta.data %>%
  select(
    cluster = recluster_ids,
    stim,
    Apoptosis_Score1
  ) %>%
  filter(stim %in% c("Control", "PTH"))


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
      slice_sample(filter(.x, stim == "PTH"),     n = n_use)
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
      wilcox.test(Apoptosis_Score1 ~ stim, exact = FALSE)$p.value,
      error = function(e) NA
    ),
    # median scores per condition
    median_control = median(Apoptosis_Score1[stim == "Control"], na.rm = TRUE),
    median_PTH     = median(Apoptosis_Score1[stim == "PTH"], na.rm = TRUE),
    # delta median
    delta_median   = median_PTH - median_control,
    # number of cells per condition (should be equal after downsampling)
    n_control = sum(stim == "Control"),
    n_PTH     = sum(stim == "PTH"),
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
        Apoptosis_Score1[stim == "PTH"],
        Apoptosis_Score1[stim == "Control"]
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

write.csv(final_results, file = "final_apoptosis_module_score.csv")


#### ====== testing for multiple comparisons ====

comparisons <- list(
  c("Control", "PTH"),
  c("Control", "OPG-Fc"),
  c("Control", "OPG-Fc + PTH"),
  c("PTH","OPG-Fc + PTH")
)


df <- recluster@meta.data %>%
  select(
    cluster = recluster_ids,
    stim,
    TGFB_EMT1
  ) %>%
  filter(stim %in% c("Control", "PTH", "OPG-Fc", "OPG-Fc + PTH"))


analyze_pair <- function(df, stim_pair) {
  
  # subset only the two conditions
  df_pair <- df %>% filter(stim %in% stim_pair)
  
  # Step 1: downsample per cluster
  min_n <- df_pair %>%
    dplyr::count(cluster, stim) %>%
    group_by(cluster) %>%
    summarise(n_min = min(n), .groups = "drop")
  
  df_balanced <- df_pair %>%
    left_join(min_n, by = "cluster") %>%
    group_by(cluster) %>%
    group_modify(~ {
      n_use <- unique(.x$n_min)
      bind_rows(
        slice_sample(filter(.x, stim == stim_pair[1]), n = n_use),
        slice_sample(filter(.x, stim == stim_pair[2]), n = n_use)
      )
    }) %>%
    ungroup() %>%
    select(-n_min)
  
  # Step 2: Wilcoxon + effect size
  wilcox_results <- df_balanced %>%
    group_by(cluster) %>%
    summarise(
      p_value = tryCatch(
        wilcox.test(TGFB_EMT1 ~ stim, exact = FALSE)$p.value,
        error = function(e) NA
      ),
      median_1 = median(TGFB_EMT1[stim == stim_pair[1]], na.rm = TRUE),
      median_2 = median(TGFB_EMT1[stim == stim_pair[2]], na.rm = TRUE),
      delta_median = median_2 - median_1,
      n_1 = sum(stim == stim_pair[1]),
      n_2 = sum(stim == stim_pair[2]),
      .groups = "drop"
    ) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH"))
  
  # Step 3: Cliff's delta
  cliff_results <- df_balanced %>%
    group_by(cluster) %>%
    summarise(
      cliffs_delta = tryCatch(
        cliff.delta(
          TGFB_EMT1[stim == stim_pair[2]],
          TGFB_EMT1[stim == stim_pair[1]]
        )$estimate,
        error = function(e) NA
      ),
      .groups = "drop"
    )
  
  # Step 4: merge
  final_results <- wilcox_results %>%
    left_join(cliff_results, by = "cluster") %>%
    mutate(
      comparison = paste(stim_pair[2], "vs", stim_pair[1]),
      cliffs_magnitude = case_when(
        abs(cliffs_delta) < 0.147 ~ "negligible",
        abs(cliffs_delta) < 0.33  ~ "small",
        abs(cliffs_delta) < 0.474 ~ "medium",
        TRUE                       ~ "large"
      )
    )
  
  return(final_results)
}

results_all <- lapply(comparisons, function(pair) analyze_pair(df, pair)) %>%
  bind_rows()

# inspect
results_all

write.csv(results_all, file = "final_TGFB_module_score.csv")



#extracting gene universe

gene.universe <- rownames(recluster[["RNA"]]@features)

write.csv(gene.universe, "gene_universe.csv", row.names = F)


# Seurat > GEO submission workflow

counts <- GetAssayData(recluster, assay = "RNA", layer = "counts")

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

meta <- recluster@meta.data

write.csv(meta, file = "cell_metadata.csv")
