silibrary(Seurat)
library(DoubletFinder)
library(dplyr)
library(ggplot2)
library(SingleR)
library(celldex)
library(SeuratObject)

setwd("/media/ayat/Study/Master_1st_sem/single_cell/project1")

# Load data
BMMC1_D1T1 <- readRDS("GSM4138872_scRNA_BMMC_D1T1.rds")
BMMC_D1T2 <- readRDS("GSM4138873_scRNA_BMMC_D1T2.rds")
CD34_D2T1 <- readRDS("GSM4138874_scRNA_CD34_D2T1.rds")
CD34_D3T1 <- readRDS("GSM4138875_scRNA_CD34_D3T1.rds")

# Create Seurat objects
samples <- list(
  seurat_object_bmmc_d1t1 = CreateSeuratObject(counts = BMMC1_D1T1, min.cells = 3, min.features = 200),
  seurat_object_bmmc_d2t1 = CreateSeuratObject(counts = BMMC_D1T2, min.cells = 3, min.features = 200),
  seurat_object_cd34_d2t1 = CreateSeuratObject(counts = CD34_D2T1, min.cells = 3, min.features = 200),
  seurat_object_cd34_d3t1 = CreateSeuratObject(counts = CD34_D3T1, min.cells = 3, min.features = 200)
)

# Add metadata
metadata <- data.frame(
  Donor = c("D1", "D1", "D2", "D3"),
  Replicate = c("T1", "T2", "T1", "T1"),
  Sex = c("F", "F", "M", "F"),
  row.names = names(samples)
)

for (sample_name in names(samples)) {
  samples[[sample_name]]$Donor <- metadata[sample_name, "Donor"]
  samples[[sample_name]]$Replicate <- metadata[sample_name, "Replicate"]
  samples[[sample_name]]$Sex <- metadata[sample_name, "Sex"]
}

# Assign each Seurat object to its own variable
for (name in names(samples)) {
  assign(name, samples[[name]])
}
rm(samples)
## How many cells are in each sample?
n_BMMC1_D1T1 <- ncol(BMMC1_D1T1)
n_BMMC_D1T2 <- ncol(BMMC_D1T2)
n_CD34_D2T1 <- ncol(CD34_D2T1)
n_CD34_D3T1 <- ncol(CD34_D3T1)

## How many genes are in the expression matrices?
g_BMMC1_D1T1 <- nrow(BMMC1_D1T1)
g_BMMC_D1T2 <- nrow(BMMC_D1T2)
g_CD34_D2T1 <- nrow(CD34_D2T1)
g_CD34_D3T1 <- nrow(CD34_D3T1)

## What information is now part of the meta-data of the objects?
head(seurat_object_bmmc_d1t1@meta.data)
head(seurat_object_bmmc_d2t1@meta.data)
head(seurat_object_cd34_d2t1@meta.data)
head(seurat_object_cd34_d3t1@meta.data)

#------------------week 2 ------------------------

##Quality Control 

# BMMC1_D1T1 object analysis

# Add mitochondrial gene percentage to metadata
# This calculates the percentage of mitochondrial gene expression for each cell
seurat_object_bmmc_d1t1[["percent.mt"]] <- PercentageFeatureSet(seurat_object_bmmc_d1t1, pattern = "^MT")

# Plot distribution of QC metrics: number of features, counts, and mitochondrial gene percentages
png("d1t1_distr.png")
VlnPlot(seurat_object_bmmc_d1t1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Generate scatter plots for QC metrics relationships
# Scatter plot of RNA counts vs mitochondrial gene percentage and RNA counts vs number of features
png("d1t1_scatter.png")
plot_1 <- FeatureScatter(seurat_object_bmmc_d1t1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot_2 <- FeatureScatter(seurat_object_bmmc_d1t1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot_1, plot_2))
dev.off()

# Filter cells based on feature count and mitochondrial percentage
# Retain cells with 200-2500 genes and less than 5% mitochondrial content
seurat_object_bmmc_d1t1 <- subset(seurat_object_bmmc_d1t1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt <= 5)

# Normalize the data
# Log normalize the data to make it comparable across cells
seurat_object_bmmc_d1t1 <- NormalizeData(seurat_object_bmmc_d1t1)

# Find highly variable features
# Identify genes with high variance for downstream analyses
seurat_object_bmmc_d1t1 <- FindVariableFeatures(seurat_object_bmmc_d1t1, selection.method = "vst", nfeatures = 2000)

# Plot the top 10 most variable genes
top10_bmmc <- head(VariableFeatures(seurat_object_bmmc_d1t1), 10)
plot1_bmmc <- VariableFeaturePlot(seurat_object_bmmc_d1t1)
png("d1t1_top_variable_genes.png")
LabelPoints(plot = plot1_bmmc, points = top10_bmmc, repel = TRUE)
dev.off()

# Scale the data to remove technical noise
seurat_object_bmmc_d1t1 <- ScaleData(seurat_object_bmmc_d1t1)

# Perform PCA for dimensionality reduction
seurat_object_bmmc_d1t1 <- RunPCA(seurat_object_bmmc_d1t1)

# Plot an elbow plot to determine the optimal number of principal components to use
png("d1t1_elbow_plot.png")
ElbowPlot(seurat_object_bmmc_d1t1)
dev.off()

# Plot a heatmap of the top PCs for visualization
png("d1t1_pc_heatmap.png")
PCHeatmap(seurat_object_bmmc_d1t1, dims = 1:20, cells = 500, balanced = TRUE, ncol = 4)
dev.off()

# Find neighbors and clusters based on PCA results
seurat_object_bmmc_d1t1 <- FindNeighbors(seurat_object_bmmc_d1t1, dims = 1:20)
seurat_object_bmmc_d1t1 <- FindClusters(seurat_object_bmmc_d1t1)

# Perform UMAP for visualization in a 2D space
seurat_object_bmmc_d1t1 <- RunUMAP(seurat_object_bmmc_d1t1, dims = 1:20)
png("d1t1_umap.png")
DimPlot(seurat_object_bmmc_d1t1, reduction = "umap", label = TRUE)
dev.off()

# DoubletFinder: Perform parameter sweep to find the optimal pK
sweep.res_bmmc1 <- paramSweep(seurat_object_bmmc_d1t1, PCs = 1:20, sct = FALSE)
sweep.stats_bmmc1 <- summarizeSweep(sweep.res_bmmc1, GT = FALSE)
bcmvn_BMMCT1 <- find.pK(sweep.stats_bmmc1)

# Plot BC metric vs pK to visualize the optimal parameter for doublet detection
png("d1t1_bcmetric_vs_pk.png")
ggplot(bcmvn_BMMCT1, aes(x = pK, y = BCmetric)) + 
  geom_point() + 
  labs(title = "Optimal pK Selection", x = "pK", y = "BC metric")
dev.off()

# Extract the optimal pK value for DoubletFinder
pK_bmmc1 <- bcmvn_BMMCT1 %>%
  filter(BCmetric == max(BCmetric)) %>%
  dplyr::select(pK)
pK_bmmc1 <- as.numeric(as.character(pK_bmmc1[[1]]))

# Estimate the expected number of doublets based on the dataset
expected_doublet_rate <- 0.05
nExp <- round(ncol(seurat_object_bmmc_d1t1) * expected_doublet_rate)

# Model homotypic doublets based on cluster annotations
annotations_bmmc1 <- seurat_object_bmmc_d1t1@meta.data$seurat_clusters
homotypic.prop_bmmc1 <- modelHomotypic(annotations_bmmc1)

# Adjust the number of expected doublets based on homotypic proportion
nExp_poi_bmmc1 <- round(0.75 * nrow(seurat_object_bmmc_d1t1@meta.data))
nExp_poi.adj_bmmc1 <- round(nExp_poi_bmmc1 * (1 - homotypic.prop_bmmc1))

# Run DoubletFinder to identify and annotate doublets in the dataset
seurat_object_bmmc_d1t1 <- doubletFinder(
  seurat_object_bmmc_d1t1,
  PCs = 1:20,
  pN = 0.25,
  pK = pK_bmmc1,
  nExp = nExp_poi_bmmc1,
  reuse.pANN = FALSE,
  sct = FALSE
)

# Plot UMAP before doublet removal
png("bmmc_d1t1_umap_before_doublet_removal.png")
# Detect the column name matching the pattern
df_class_col <- grep("^DF\\.classifications_\\d+.*", colnames(seurat_object_bmmc_d1t1@meta.data), value = TRUE)[1]

# Use the detected column in DimPlot
DimPlot(
  seurat_object_bmmc_d1t1,
  reduction = "umap",
  group.by = df_class_col,
  label = TRUE
) + ggtitle("UMAP Before Doublet Removal (BMMC_D1T1)")
dev.off()

table(seurat_object_bmmc_d1t1$DF.classifications_0.25_0.005_4276)
table(seurat_object_bmmc_d1t1$seurat_clusters, seurat_object_bmmc_d1t1$DF.classifications_0.25_0.005_4276)

# Remove doublets
seurat_object_bmmc_d1t1 <- subset(seurat_object_bmmc_d1t1, subset = DF.classifications_0.25_0.005_4276 == "Singlet")
# Plot UMAP after doublet removal
png("bmmc_d1t1_umap_after_doublet_removal.png")
DimPlot(seurat_object_bmmc_d1t1, reduction = "umap", label = TRUE) +
  ggtitle("UMAP After Doublet Removal (BMMC_D1T1)")
dev.off()

# ----------BMMC_D1T2---------
# Add mitochondrial gene percentage to metadata
seurat_object_bmmc_d2t1[["percent.mt"]] <- PercentageFeatureSet(seurat_object_bmmc_d2t1, pattern = "^MT")

# Plot distribution of QC metrics
png("d1t2_distr.png")
VlnPlot(seurat_object_bmmc_d2t1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Generate scatter plots for QC metrics relationships
png("d1t2_scatter.png")
plot1_bmmc_d2t1 <- FeatureScatter(seurat_object_bmmc_d2t1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_bmmc_d2t1 <- FeatureScatter(seurat_object_bmmc_d2t1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1_bmmc_d2t1, plot2_bmmc_d2t1))
dev.off()

# Filter cells based on feature count
seurat_object_bmmc_d2t1 <- subset(seurat_object_bmmc_d2t1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize the data
seurat_object_bmmc_d2t1 <- NormalizeData(seurat_object_bmmc_d2t1)

# Find highly variable features
seurat_object_bmmc_d2t1 <- FindVariableFeatures(seurat_object_bmmc_d2t1, selection.method = "vst", nfeatures = 2000 )

# Plot the top 10 variable genes
top10_bmmc_d2t1 <- head(VariableFeatures(seurat_object_bmmc_d2t1), 10)
plot1_bmmc_d2t1 <- VariableFeaturePlot(seurat_object_bmmc_d2t1)
png("d1t2_top_variable_genes.png")
LabelPoints(plot = plot1_bmmc_d2t1, points = top10_bmmc_d2t1, repel = TRUE)
dev.off()

# Scale the data and run PCA
seurat_object_bmmc_d2t1 <- ScaleData(seurat_object_bmmc_d2t1)
seurat_object_bmmc_d2t1 <- RunPCA(seurat_object_bmmc_d2t1)

# Plot elbow plot
png("d1t2_elbow_plot.png")
ElbowPlot(seurat_object_bmmc_d2t1)
dev.off()

# Plot PC heatmap
png("d1t2_pc_heatmap.png")
PCHeatmap(seurat_object_bmmc_d2t1, dims = 1:20, cells = 500, balanced = TRUE, ncol = 4)
dev.off()

# Find neighbors and clusters
seurat_object_bmmc_d2t1 <- FindNeighbors(seurat_object_bmmc_d2t1, dims = 1:20)
seurat_object_bmmc_d2t1 <- FindClusters(seurat_object_bmmc_d2t1)

# Run UMAP
seurat_object_bmmc_d2t1 <- RunUMAP(seurat_object_bmmc_d2t1, dims = 1:20)
png("d1t2_umap.png")
DimPlot(seurat_object_bmmc_d2t1, reduction = "umap", label = TRUE)
dev.off()

# Perform DoubletFinder parameter sweep
sweep.res_bmmc2 <- paramSweep(seurat_object_bmmc_d2t1, PCs = 1:20, sct = FALSE)
sweep.stats_bmmc2 <- summarizeSweep(sweep.res_bmmc2, GT = FALSE)
bcmvn_BMMCT2 <- find.pK(sweep.stats_bmmc2)

# Plot BCmetric vs pK
png("d1t2_bcmetric_vs_pk.png")
ggplot(bcmvn_BMMCT2, aes(x = as.numeric(as.character(pK)), y = BCmetric)) + 
  geom_point() + 
  geom_line() + 
  labs(title = "Optimal pK Selection (BMMC_D1T2)", x = "pK", y = "BC metric") + 
  theme_minimal()
dev.off()

# Extract optimal pK
pK_bmmc2 <- bcmvn_BMMCT2 %>%
  filter(BCmetric == max(BCmetric)) %>%
  dplyr::select(pK)
pK_bmmc2 <- as.numeric(as.character(pK_bmmc2[[1]]))

# Model homotypic doublets
annotations_bmmc2 <- seurat_object_bmmc_d2t1@meta.data$seurat_clusters
homotypic.prop_bmmc2 <- modelHomotypic(annotations_bmmc2)

# Calculate expected doublets
nExp_poi_bmmc2 <- round(0.75 * nrow(seurat_object_bmmc_d2t1@meta.data))
nExp_poi.adj_bmmc2 <- round(nExp_poi_bmmc2 * (1 - homotypic.prop_bmmc2))

# Run DoubletFinder
seurat_object_bmmc_d2t1 <- doubletFinder(seurat_object_bmmc_d2t1, PCs = 1:20, pN = 0.25, pK = pK_bmmc2, nExp = nExp_poi_bmmc2, reuse.pANN = FALSE, sct = FALSE)

png("UMAP Before Doublet Removal (BMMC_D2T1).png")
# Plot UMAP before doublet removal
df_class_col2 <- grep("^DF\\.classifications_\\d+.*", colnames(seurat_object_bmmc_d2t1@meta.data), value = TRUE)[1]

# Use the detected column in DimPlot
DimPlot(
  seurat_object_bmmc_d2t1,
  reduction = "umap",
  group.by = df_class_col2,
  label = TRUE
) + ggtitle("UMAP Before Doublet Removal (BMMC_D2T1)")
dev.off()

table(seurat_object_bmmc_d2t1$DF.classifications_0.25_0.005_4330)
table(seurat_object_bmmc_d2t1$seurat_clusters, seurat_object_bmmc_d2t1$DF.classifications_0.25_0.005_4330)


# Remove doublets
seurat_object_bmmc_d2t1 <- subset(seurat_object_bmmc_d2t1, subset = DF.classifications_0.25_0.005_4330 == "Singlet")

# Plot UMAP after doublet removal
png("bmmc_d2t1_umap_after_doublet_removal.png")
DimPlot(seurat_object_bmmc_d2t1, reduction = "umap", label = TRUE) +
  ggtitle("UMAP After Doublet Removal (BMMC_D1T2)")
dev.off()

# ------------CD34_D2T1-----------------------

# Add mitochondrial gene percentage to metadata
seurat_object_cd34_d2t1[["percent.mt"]] <- PercentageFeatureSet(seurat_object_cd34_d2t1, pattern = "^MT")

# Plot distribution of QC metrics
png("cd34_d2t1_distr.png")
VlnPlot(seurat_object_cd34_d2t1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Scatter plots for QC metrics
png("cd34_d2t1_scatter.png")
plot1_cd34_d2t1 <- FeatureScatter(seurat_object_cd34_d2t1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_cd34_d2t1 <- FeatureScatter(seurat_object_cd34_d2t1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1_cd34_d2t1, plot2_cd34_d2t1))
dev.off()

# Filter cells
seurat_object_cd34_d2t1 <- subset(seurat_object_cd34_d2t1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500& percent.mt < 5)

# Normalize the data
seurat_object_cd34_d2t1 <- NormalizeData(seurat_object_cd34_d2t1)

# Find highly variable features
seurat_object_cd34_d2t1 <- FindVariableFeatures(seurat_object_cd34_d2t1, selection.method = "vst", nfeatures = 2000)

# Plot the top 10 variable genes
top10_cd34_d2t1 <- head(VariableFeatures(seurat_object_cd34_d2t1), 10)
plot1_cd34_d2t1 <- VariableFeaturePlot(seurat_object_cd34_d2t1)
png("cd34_d2t1_top_variable_genes.png")
LabelPoints(plot = plot1_cd34_d2t1, points = top10_cd34_d2t1, repel = TRUE)
dev.off()

# Scale data and run PCA
seurat_object_cd34_d2t1 <- ScaleData(seurat_object_cd34_d2t1)
seurat_object_cd34_d2t1 <- RunPCA(seurat_object_cd34_d2t1)

# Plot elbow plot
png("cd34_d2t1_elbow_plot.png")
ElbowPlot(seurat_object_cd34_d2t1)
dev.off()

# Plot PC heatmap
png("cd34_d2t1_pc_heatmap.png")
PCHeatmap(seurat_object_cd34_d2t1, dims = 1:20, cells = 500, balanced = TRUE, ncol = 4)
dev.off()

# Find neighbors and clusters
seurat_object_cd34_d2t1 <- FindNeighbors(seurat_object_cd34_d2t1, dims = 1:20)
seurat_object_cd34_d2t1 <- FindClusters(seurat_object_cd34_d2t1)

# Run UMAP
seurat_object_cd34_d2t1 <- RunUMAP(seurat_object_cd34_d2t1, dims = 1:20)
png("cd34_d2t1_umap.png")
DimPlot(seurat_object_cd34_d2t1, reduction = "umap", label = TRUE)
dev.off()

# Perform DoubletFinder parameter sweep
sweep.res_cd34_d2t1 <- paramSweep(seurat_object_cd34_d2t1, PCs = 1:20, sct = FALSE)
sweep.stats_cd34_d2t1 <- summarizeSweep(sweep.res_cd34_d2t1, GT = FALSE)
bcmvn_cd34_d2t1 <- find.pK(sweep.stats_cd34_d2t1)

# Plot BCmetric vs pK
png("cd34_d2t1_bcmetric_vs_pk.png")
ggplot(bcmvn_cd34_d2t1, aes(x = as.numeric(as.character(pK)), y = BCmetric)) + 
  geom_point() + 
  geom_line() + 
  labs(title = "Optimal pK Selection (CD34_D2T1)", x = "pK", y = "BC metric") + 
  theme_minimal()
dev.off()

# Extract optimal pK
pK_cd34_d2t1 <- bcmvn_cd34_d2t1 %>%
  filter(BCmetric == max(BCmetric)) %>%
  dplyr::select(pK)
pK_cd34_d2t1 <- as.numeric(as.character(pK_cd34_d2t1[[1]]))

# Model homotypic doublets
annotations_cd34_d2t1 <- seurat_object_cd34_d2t1@meta.data$seurat_clusters
homotypic.prop_cd34_d2t1 <- modelHomotypic(annotations_cd34_d2t1)

# Calculate expected doublets
nExp_poi_cd34_d2t1 <- round(0.75 * nrow(seurat_object_cd34_d2t1@meta.data))
nExp_poi.adj_cd34_d2t1 <- round(nExp_poi_cd34_d2t1 * (1 - homotypic.prop_cd34_d2t1))

# Run DoubletFinder
seurat_object_cd34_d2t1 <- doubletFinder(seurat_object_cd34_d2t1, PCs = 1:20, pN = 0.25, pK = pK_cd34_d2t1, nExp = nExp_poi_cd34_d2t1, reuse.pANN = FALSE, sct = FALSE)

# Plot UMAP before doublet removal
png("cd34_d2t1_umap_before_doublet_removal.png")
DimPlot(seurat_object_cd34_d2t1, reduction = "umap", group.by = "DF.classifications_0.25_0.3_1159", label = TRUE) +
  ggtitle("UMAP Before Doublet Removal (CD34_D2T1)")
dev.off()

# Remove doublets
seurat_object_cd34_d2t1 <- subset(seurat_object_cd34_d2t1, subset = DF.classifications_0.25_0.3_1159 == "Singlet")

# Plot UMAP after doublet removal
png("cd34_d2t1_umap_after_doublet_removal.png")
DimPlot(seurat_object_cd34_d2t1, reduction = "umap", label = TRUE) +
  ggtitle("UMAP After Doublet Removal (CD34_D2T1)")
dev.off()


#-------------- CD34_D3T1 -----------

# Add mitochondrial gene percentage to metadata
seurat_object_cd34_d3t1[["percent.mt"]] <- PercentageFeatureSet(seurat_object_cd34_d3t1, pattern = "^MT")

# Plot distribution of QC metrics: nFeature_RNA, nCount_RNA, and percent.mt
png("cd34_d3t1_distr.png")
VlnPlot(seurat_object_cd34_d3t1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Generate scatter plots for QC metrics relationships
png("cd34_d3t1_scatter.png")
plot1_cd34_d3t1 <- FeatureScatter(seurat_object_cd34_d3t1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_cd34_d3t1 <- FeatureScatter(seurat_object_cd34_d3t1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1_cd34_d3t1, plot2_cd34_d3t1))
dev.off()

# Filter cells based on feature count
# Retain cells with 200-2500 genes and less than 5% mitochondrial content
seurat_object_cd34_d3t1 <- subset(seurat_object_cd34_d3t1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize the data
seurat_object_cd34_d3t1 <- NormalizeData(seurat_object_cd34_d3t1)

# Find highly variable features
seurat_object_cd34_d3t1 <- FindVariableFeatures(seurat_object_cd34_d3t1, selection.method = "vst", nfeatures = 2000 )

# Plot the top 10 variable genes
top10_cd34_d3t1 <- head(VariableFeatures(seurat_object_cd34_d3t1), 10)
plot1_cd34_d3t1 <- VariableFeaturePlot(seurat_object_cd34_d3t1)
png("cd34_d3t1_top_variable_genes.png")
LabelPoints(plot = plot1_cd34_d3t1, points = top10_cd34_d3t1, repel = TRUE)
dev.off()

# Scale the data and run PCA
seurat_object_cd34_d3t1 <- ScaleData(seurat_object_cd34_d3t1)
seurat_object_cd34_d3t1 <- RunPCA(seurat_object_cd34_d3t1)

# Elbow plot to determine the number of PCs to use
png("cd34_d3t1_elbow_plot.png")
ElbowPlot(seurat_object_cd34_d3t1)
dev.off()

# Plot PC heatmap
png("cd34_d3t1_pc_heatmap.png")
PCHeatmap(seurat_object_cd34_d3t1, dims = 1:20, cells = 500, balanced = TRUE, ncol = 4)
dev.off()

# Find neighbors and clusters
seurat_object_cd34_d3t1 <- FindNeighbors(seurat_object_cd34_d3t1, dims = 1:20)
seurat_object_cd34_d3t1 <- FindClusters(seurat_object_cd34_d3t1)

# Run UMAP for dimensionality reduction
seurat_object_cd34_d3t1 <- RunUMAP(seurat_object_cd34_d3t1, dims = 1:20)
png("cd34_d3t1_umap.png")
DimPlot(seurat_object_cd34_d3t1, reduction = "umap", label = TRUE)
dev.off()

# Perform DoubletFinder parameter sweep
sweep.res_cd34_d3t1 <- paramSweep(seurat_object_cd34_d3t1, PCs = 1:20, sct = FALSE)
sweep.stats_cd34_d3t1 <- summarizeSweep(sweep.res_cd34_d3t1, GT = FALSE)
bcmvn_cd34_d3t1 <- find.pK(sweep.stats_cd34_d3t1)

# Plot BCmetric vs pK
png("cd34_d3t1_bcmetric_vs_pk.png")
ggplot(bcmvn_cd34_d3t1, aes(x = as.numeric(as.character(pK)), y = BCmetric)) + 
  geom_point() + 
  geom_line() + 
  labs(title = "Optimal pK Selection (CD34_D3T1)", x = "pK", y = "BC Metric") + 
  theme_minimal()
dev.off()

# Extract optimal pK
pK_cd34_d3t1 <- bcmvn_cd34_d3t1 %>%
  filter(BCmetric == max(BCmetric)) %>%
  dplyr::select(pK)
pK_cd34_d3t1 <- as.numeric(as.character(pK_cd34_d3t1[[1]]))

# Model homotypic doublets
annotations_cd34_d3t1 <- seurat_object_cd34_d3t1@meta.data$seurat_clusters
homotypic.prop_cd34_d3t1 <- modelHomotypic(annotations_cd34_d3t1)

# Calculate expected doublets
nExp_poi_cd34_d3t1 <- round(0.75 * nrow(seurat_object_cd34_d3t1@meta.data))
nExp_poi.adj_cd34_d3t1 <- round(nExp_poi_cd34_d3t1 * (1 - homotypic.prop_cd34_d3t1))

# Run DoubletFinder
seurat_object_cd34_d3t1 <- doubletFinder(
  seurat_object_cd34_d3t1,
  PCs = 1:20,
  pN = 0.25,
  pK = pK_cd34_d3t1,
  nExp = nExp_poi_cd34_d3t1,
  reuse.pANN = FALSE,
  sct = FALSE
)

# Plot UMAP before doublet removal
png("cd34_d3t1_umap_before_doublet_removal.png")
DimPlot(seurat_object_cd34_d3t1, reduction = "umap", group.by = "DF.classifications_0.25_0.03_4026", label = TRUE) +
  ggtitle("UMAP Before Doublet Removal (CD34_D3T1)")
dev.off()

# Remove doublets
seurat_object_cd34_d3t1 <- subset(seurat_object_cd34_d3t1, subset = DF.classifications_0.25_0.03_4026 == "Singlet")

# Plot UMAP after doublet removal
png("cd34_d3t1_umap_after_doublet_removal.png")
DimPlot(seurat_object_cd34_d3t1, reduction = "umap", label = TRUE) +
  ggtitle("UMAP After Doublet Removal (CD34_D3T1)")
dev.off()

#-------------------------Batch correction-------------------------------------
# Merge datasets (before batch correction)
merged_data <- merge(
  seurat_object_bmmc_d1t1,
  y = c(seurat_object_bmmc_d2t1, seurat_object_cd34_d2t1, seurat_object_cd34_d3t1),
  add.cell.ids = c("BMMC_D1T1", "BMMC_D1T2", "CD34_D2T1", "CD34_D3T1")
)
#Step 1: Perform Clustering and UMAP Without Batch Correction# Scale the data
merged_data <- ScaleData(merged_data)

# Run PCA
merged_data <- RunPCA(merged_data)

# Run UMAP
merged_data <- RunUMAP(merged_data, dims = 1:20)

# Find neighbors and clusters
merged_data <- FindNeighbors(merged_data, dims = 1:20)
merged_data <- FindClusters(merged_data, resolution = 0.5)

saveRDS(merged_data,"merged_data_postprocessed.RDS")
# Plot UMAP
png("merged_no_batch_umap.png")
DimPlot(merged_data, reduction = "umap", group.by = "orig.ident", label = TRUE)
dev.off()

png("merged_no_batch_umap_sex.png")
DimPlot(merged_data, reduction = "umap", group.by = "Sex", label = TRUE) +
  ggtitle("UMAP Without Batch Correction (Grouped by Sex)")
dev.off()

png("merged_no_batch_umap_donor.png")
DimPlot(merged_data, reduction = "umap", group.by = "Donor", label = TRUE) +
  ggtitle("UMAP Without Batch Correction (Grouped by Donor)")
dev.off()


png("merged_no_batch_umap_replicate.png")
DimPlot(merged_data, reduction = "umap", group.by = "Replicate", label = TRUE) +
  ggtitle("UMAP Without Batch Correction (Grouped by Replicate)")
dev.off()

#integration steps
#identify integration features
integration_features <- SelectIntegrationFeatures(object.list = list(
  seurat_object_bmmc_d1t1,
  seurat_object_bmmc_d2t1,
  seurat_object_cd34_d2t1,
  seurat_object_cd34_d3t1
))

#identify integration anchors

integration_anchors <- FindIntegrationAnchors(
  object.list = list(
    seurat_object_bmmc_d1t1,
    seurat_object_bmmc_d2t1,
    seurat_object_cd34_d2t1,
    seurat_object_cd34_d3t1
  ),
  anchor.features = integration_features
)

#integrate data
integrated_data <- IntegrateData(anchorset = integration_anchors)

#Steps After Integration
#scaling,PCA,UMAP
integrated_data <- ScaleData(integrated_data)
integrated_data <- RunPCA(integrated_data)
integrated_data <- RunUMAP(integrated_data, dims = 1:20)
#clustering
integrated_data <- FindNeighbors(integrated_data, dims = 1:20)
integrated_data <- FindClusters(integrated_data, resolution = 0.5)

#visualization
png("integrated_batch_umap.png")
DimPlot(integrated_data, reduction = "umap", group.by = "orig.ident", label = TRUE)
dev.off()

png("integrated_batch_umap_sex.png")
DimPlot(integrated_data, reduction = "umap", group.by = "Sex", label = TRUE) +
  ggtitle("UMAP With Batch Correction (Grouped by Sex)")
dev.off()

png("integrated_batch_umap_donor.png")
DimPlot(integrated_data, reduction = "umap", group.by = "Donor", label = TRUE) +
  ggtitle("UMAP With Batch Correction (Grouped by Donor)")
dev.off()

png("integrated_batch_umap_replicate.png")
DimPlot(integrated_data, reduction = "umap", group.by = "Replicate", label = TRUE) +
  ggtitle("UMAP With Batch Correction (Grouped by Replicate)")
dev.off()

length(unique(integrated_data$seurat_clusters))

saveRDS(integrated_data,"integrated_data_1.RDS")
#---------------------------------------

#-----------week 3----------
#Automatic annotation
# Load the Human Primary Cell Atlas reference data
reference <- celldex::HumanPrimaryCellAtlasData()

# Extract the normalized expression matrix from the Seurat object
data_matrix <- GetAssayData(integrated_data, layer = "data")

# Perform automatic annotation
singler_annotations <- SingleR(
  test = data_matrix,
  ref = reference,
  labels = reference$label.main
)

# Add SingleR annotations to Seurat metadata
integrated_data$SingleR <- singler_annotations$labels


# Plot UMAP with annotations
png("automatic_annotation_umap.png")
DimPlot(integrated_data, reduction = "umap", group.by = "SingleR", label = TRUE) +
  ggtitle("UMAP with Automatic Cell Type Annotation")
dev.off()


png("umap_seurat_clusters.png")
DimPlot(integrated_data, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("UMAP with Seurat Clusters")
dev.off()


#Manual annotation
# Identify markers for all clusters
markers <- FindAllMarkers(
  object = integrated_data,
  only.pos = TRUE,  # Keep only positive markers
  min.pct = 0.25,  # Minimum percentage of cells expressing the gene
  logfc.threshold = 0.25  # Minimum log fold-change
)

# Save markers for reference
saveRDS(markers, "markers_all_clusters.rds")

# Check the top markers for each cluster
top_markers<-markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)


# Step 3: Define Feature Markers Based on Table 2
feature_HSC <- c("CD34", "CD38", "Sca1", "Kit")
feature_LMPP <- c("CD38", "CD52", "CSF3R", "ca1", "Kit", "CD34", "Flk2")
feature_CLP <- c("IL7R")
feature_GMPN <- c("ELANE")
feature_CMP <- c("IL3", "GM-CSF", "M-CSF")
feature_B <- c("CD19", "CD20", "CD38")
feature_PreB <- c("CD19", "CD34")
feature_Plasma <- c("SDC1", "IGHA1", "IGLC1", "MZB1", "JCHAIN")
feature_Cd8 <- c("CD3D", "CD3E", "CD8A", "CD8B")
feature_Cd4 <- c("CD3D", "CD3E", "CD4")
feature_NK <- c("FCGR3A", "NCAM1", "NKG7", "KLRB1")
feature_Erythrocytes <- c("GATA1", "HBB", "HBA1", "HBA2")
feature_pDC <- c("IRF8", "IRF4", "IRF7")
feature_cDC <- c("CD1C", "CD207", "ITGAM", "NOTCH2", "SIRPA")
feature_CD14 <- c("CD14", "CCL3", "CCL4", "IL1B")
feature_CD16 <- c("FCGR3A", "CD68", "S100A12")
feature_Basophils <- c("GATA2")

# UMAP Feature Plot for GMPN
png("umap_feature_GMPN.png", width = 800, height = 600)
FeaturePlot(integrated_data, features = feature_GMPN, label = TRUE, min.cutoff = "q10") 
dev.off()

# Violin Plot for GMPN
png("violin_feature_GMPN.png", width = 800, height = 600)
VlnPlot(integrated_data, features = feature_GMPN, group.by = "seurat_clusters", pt.size = 0.5) 
dev.off()


# UMAP Feature Plot for B Cells
png("umap_feature_B.png", width = 800, height = 600)
FeaturePlot(integrated_data, features = feature_B, label = TRUE, min.cutoff = "q10") 
dev.off()

# Violin Plot for B Cells
png("violin_feature_B.png", width = 800, height = 600)
VlnPlot(integrated_data, features = feature_B, group.by = "seurat_clusters", pt.size = 0.5) 
dev.off()

# UMAP Feature Plot for CD8+ T Cells
png("umap_feature_Cd8.png", width = 800, height = 600)
FeaturePlot(integrated_data, features = feature_Cd8, label = TRUE, min.cutoff = "q10")
dev.off()

# Violin Plot for CD8+ T Cells
png("violin_feature_Cd8.png", width = 800, height = 600)
VlnPlot(integrated_data, features = feature_Cd8, group.by = "seurat_clusters", pt.size = 0.5) 
dev.off()

# UMAP Feature Plot for CD4+ T Cells
png("umap_feature_Cd4.png", width = 800, height = 600)
FeaturePlot(integrated_data, features = feature_Cd4, label = TRUE, min.cutoff = "q10")
dev.off()

# Violin Plot for CD4+ T Cells
png("violin_feature_Cd4.png", width = 800, height = 600)
VlnPlot(integrated_data, features = feature_Cd4, group.by = "seurat_clusters", pt.size = 0.5)
dev.off()


# UMAP Feature Plot for pDC
png("umap_feature_pDC.png", width = 800, height = 600)
FeaturePlot(integrated_data, features = feature_pDC, label = TRUE, min.cutoff = "q10") 
dev.off()

# Violin Plot for pDC
png("violin_feature_pDC.png", width = 800, height = 600)
VlnPlot(integrated_data, features = feature_pDC, group.by = "seurat_clusters", pt.size = 0.5)
dev.off()


#  Manually Assign Cluster Labels
cluster_ids <- c(
  "0" = "CD4+ T Cells (CD4)",
  "1" = "Granulocyte-Monocyte Progenitors (GMP)",
  "2" = "Common Lymphoid Progenitors (CLP)",
  "3" = "Hematopoietic Stem Cells (HSC)",
  "4" = "Natural Killer Cells (NK)",
  "5" = "Erythrocytes (Ery)",
  "6" = "Plasma Cells",
  "7" = "CD8+ T Cells (CD8)",
  "8" = "B Cells (B)",
  "9" = "Plasmacytoid Dendritic Cells (pDC)",
  "10" = "Pre B Cells",
  "11" = "Conventional Dendritic Cells (cDC)",
  "12" = "CD14+ Monocytes",
  "13" = "Basophils"
)

# Apply manual annotations to metadata
integrated_data$manual_annotation <- factor(integrated_data$seurat_clusters, levels = names(cluster_ids))
levels(integrated_data$manual_annotation) <- cluster_ids

# Save the updated Seurat object
saveRDS(integrated_data, "integrated_data_with_manual_annotation_final.rds")

# Step 6: Visualization of Manual Annotations
# UMAP with Manual Annotations
png("manual_annotation_umap.png", width = 800, height = 600)
DimPlot(integrated_data, reduction = "umap", group.by = "manual_annotation", label = TRUE) +
  ggtitle("UMAP with Manual Annotations")
dev.off()

# Step 7: Compare Manual and Automatic Annotations
# UMAP with Automatic Annotations (if available)
png("automatic_annotation_umap.png", width = 800, height = 600)
DimPlot(integrated_data, reduction = "umap", group.by = "SingleR", label = TRUE) +
  ggtitle("UMAP with Automatic Annotations")
dev.off()

# Step 8: Validate Key Marker Genes Across Clusters
# Violin Plot for Selected Markers
selected_markers <- c("CD34", "S100A9", "HBA1") # Adjust based on key markers
png("violin_plot_markers.png", width = 800, height = 600)
VlnPlot(integrated_data, features = selected_markers, group.by = "manual_annotation", pt.size = 0.5) +
  ggtitle("Violin Plot for Selected Marker Genes")
dev.off()

# UMAP Feature Plots for Selected Marker Genes
for (gene in selected_markers) {
  png(paste0(gene, "_feature_umap.png"), width = 800, height = 600)
  FeaturePlot(integrated_data, features = gene, reduction = "umap") +
    ggtitle(paste("UMAP Plot of", gene))
  dev.off()
}

# Step 9: (Optional) Merge Similar Clusters for Final Labels
final_labels <- c(
  "CD4+ T Cells (CD4)" = "T Cells",
  "CD8+ T Cells (CD8)" = "T Cells",
  "Granulocyte-Monocyte Progenitors (GMP)" = "Myeloid Progenitors",
  "Erythrocytes (Ery)" = "Erythrocytes",
  "B Cells (B)" = "B Cells",
  "Natural Killer Cells (NK)" = "NK Cells",
  "Plasma Cells" = "Plasma Cells",
  "Plasmacytoid Dendritic Cells (pDC)" = "Dendritic Cells",
  "Conventional Dendritic Cells (cDC)" = "Dendritic Cells",
  "CD14+ Monocytes" = "Monocytes",
  "Pre B Cells" = "B Cell Progenitors",
  "Basophils" = "Basophils"
)

integrated_data$final_annotation <- recode(integrated_data$manual_annotation, !!!final_labels)

# Save final annotated dataset
saveRDS(integrated_data, "integrated_data_with_final_annotations.rds")

# Final UMAP with Combined Labels
png("final_annotation_umap.png", width = 800, height = 600)
DimPlot(integrated_data, reduction = "umap", group.by = "final_annotation", label = TRUE) +
  ggtitle("UMAP with Final Annotations")
dev.off()

#------------
# Compute proportions of cell types for each sample
cell_type_proportions <- integrated_data@meta.data %>%
  group_by(orig.ident, SingleR) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count))

# Save proportions to a CSV file for reference
write.csv(cell_type_proportions, "cell_type_proportions.csv")

# Plot cell-type proportions

png("cell_type_proportions_barplot.png", width = 1200, height = 800)
ggplot(cell_type_proportions, aes(x = orig.ident, y = proportion, fill = SingleR)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Cell-Type Proportions by Sample",
    x = "Sample",
    y = "Proportion",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#-------------Differential expression analysis----------

# Perform differential expression analysis for B cells vs T cells
B_vs_T_DE <- FindMarkers(
  object = integrated_data,
  ident.1 = "B_cell",
  ident.2 = "T_cells",
  group.by = "SingleR",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

# Save results for B cells vs T cells
write.csv(B_vs_T_DE, "B_vs_T_DE.csv")

# Perform differential expression analysis for T cells vs Monocytes
T_vs_mono_DE <- FindMarkers(
  object = integrated_data,
  ident.1 = "T_cells",
  ident.2 = "Monocyte",
  group.by = "SingleR",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

# Save results for T cells vs Monocytes
write.csv(T_vs_mono_DE, "T_vs_Monocytes_DE.csv")

# Add significance column
B_vs_T_DE$significance <- ifelse(
  B_vs_T_DE$p_val_adj < 0.05 & abs(B_vs_T_DE$avg_log2FC) > 0.25,
  ifelse(B_vs_T_DE$avg_log2FC > 0, "Upregulated in B Cells", "Upregulated in T Cells"),
  "Nonsignificant"
)

# Create volcano plot
png("volcano_B_vs_T.png", width = 1000, height = 800)
ggplot(B_vs_T_DE, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Upregulated in B Cells" = "red", "Upregulated in T Cells" = "blue", "Nonsignificant" = "grey")) +
  labs(
    title = "Differential Expression: B Cells vs T Cells",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Significance"
  ) +
  theme_minimal()
dev.off()

# Add significance column
T_vs_mono_DE$significance <- ifelse(
  T_vs_mono_DE$p_val_adj < 0.05 & abs(T_vs_mono_DE$avg_log2FC) > 0.25,
  ifelse(T_vs_mono_DE$avg_log2FC > 0, "Upregulated in T Cells", "Upregulated in Monocytes"),
  "Nonsignificant"
)

# Create volcano plot
png("volcano_T_vs_Monocytes.png", width = 1000, height = 800)
ggplot(T_vs_mono_DE, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Upregulated in T Cells" = "blue", "Upregulated in Monocytes" = "green", "Nonsignificant" = "grey")) +
  labs(
    title = "Differential Expression: T Cells vs Monocytes",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Significance"
  ) +
  theme_minimal()
dev.off()

# Extract top 5 genes for B Cells vs T Cells
top5_b_vs_t <- B_vs_T_DE %>%
  arrange(desc(abs(avg_log2FC))) %>%
  head(5) %>%
  mutate(Comparison = "B Cells vs T Cells", Gene = rownames(.))

# Extract top 5 genes for T Cells vs Monocytes
top5_t_vs_mono <- T_vs_mono_DE %>%
  arrange(desc(abs(avg_log2FC))) %>%
  head(5) %>%
  mutate(Comparison = "T Cells vs Monocytes", Gene = rownames(.))

# Combine the data
top_genes <- bind_rows(top5_b_vs_t, top5_t_vs_mono)

# Create a combined column for x-axis labeling
top_genes <- top_genes %>%
  mutate(Cell_Type = ifelse(Comparison == "B Cells vs T Cells", 
                            "B Cells vs T Cells", "T Cells vs Monocytes"))
# Generate the dot plot
png("differential_expression_dotplot.png", width = 1200, height = 800)
ggplot(top_genes, aes(x = Cell_Type, y = Gene, size = -log10(p_val_adj), color = avg_log2FC)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(
    title = "Differentially Expressed Genes: Comparison of Top 5 Genes",
    x = "Cell Type Comparison",
    y = "Genes",
    size = "-Log10 Adjusted P-value",
    color = "Log2 Fold Change"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
#-------------
#Differential Expression Between BMMC and CD34 (Independent of Cell Types)
# Add sample group metadata for BMMC and CD34
integrated_data$Sample_Group <- ifelse(grepl("BMMC", integrated_data$orig.ident), "BMMC", "CD34")

# Perform differential expression analysis
bmmc_vs_cd34_DE <- FindMarkers(
  object = integrated_data,
  ident.1 = "BMMC",
  ident.2 = "CD34",
  group.by = "Sample_Group",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

bmmc_vs_cd34_DE<-read.csv("DE_results_BMMC_vs_CD34.csv")

# Extract top 5 DEGs
top5_bmmc_vs_cd34 <- bmmc_vs_cd34_DE %>%
  arrange(desc(abs(avg_log2FC))) %>%
  head(5)


# Save the results
write.csv(bmmc_vs_cd34_DE, "BMMC_vs_CD34_DEGs.csv")
print(top5_bmmc_vs_cd34)
#    X         p_val avg_log2FC pct.1 pct.2     p_val_adj
#1    GFI1B 1.207239e-177  -8.376365 0.204 0.141 2.414479e-174
#2    RBPMS 1.624770e-181  -7.976753 0.154 0.176 3.249540e-178
#3 SERPING1 5.755628e-236  -7.610669 0.205 0.178 1.151126e-232
#4  CCDC175 1.750546e-235  -7.574297 0.181 0.148 3.501092e-232
#5    ABCG2 7.933035e-242  -7.567490 0.008 0.049 1.586607e-238
 
# Add significance column
bmmc_vs_cd34_DE$significance <- ifelse(
  bmmc_vs_cd34_DE$p_val_adj < 0.05 & abs(bmmc_vs_cd34_DE$avg_log2FC) > 0.25,
  ifelse(bmmc_vs_cd34_DE$avg_log2FC > 0, "Upregulated in BMMC", "Upregulated in CD34"),
  "Nonsignificant"
)

# Create volcano plot
png("volcano_BMMC_vs_CD34.png", width = 1000, height = 800)
ggplot(bmmc_vs_cd34_DE, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Upregulated in BMMC" = "red", "Upregulated in CD34" = "blue", "Nonsignificant" = "grey")) +
  labs(
    title = "Differential Expression: BMMC vs CD34",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Significance"
  ) +
  theme_minimal()
dev.off()


#Differential Expression Between Monocytes in BMMC and CD34
# Subset Monocytes
monocytes <- subset(integrated_data, subset = SingleR == "Monocyte")

# Add group metadata for BMMC and CD34 within Monocytes
monocytes$Sample_Group <- ifelse(grepl("BMMC", monocytes$orig.ident), "BMMC", "CD34")

# Perform differential expression analysis for Monocytes
mono_bmmc_vs_cd34_DE <- FindMarkers(
  object = monocytes,
  ident.1 = "BMMC",
  ident.2 = "CD34",
  group.by = "Sample_Group",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

# Extract top 5 DEGs
top5_mono_bmmc_vs_cd34 <- mono_bmmc_vs_cd34_DE %>%
  arrange(desc(abs(avg_log2FC))) %>%
  head(5)

# Save the results
write.csv(mono_bmmc_vs_cd34_DE, "Monocytes_BMMC_vs_CD34_DEGs.csv")
print(top5_mono_bmmc_vs_cd34)
#p_val avg_log2FC pct.1 pct.2    p_val_adj
#CPA3    2.853123e-29 -12.843008 0.170 0.198 5.706246e-26
#NPW     1.658065e-27 -10.668568 0.162 0.207 3.316129e-24
#CALN1   3.302685e-37 -10.652328 0.168 0.204 6.605370e-34
#GAL3ST4 6.914244e-18   9.508985 0.127 0.043 1.382849e-14
#HOXA9   9.343806e-32  -8.743432 0.125 0.057 1.868761e-28
# Add significance column
mono_bmmc_vs_cd34_DE$significance <- ifelse(
  mono_bmmc_vs_cd34_DE$p_val_adj < 0.05 & abs(mono_bmmc_vs_cd34_DE$avg_log2FC) > 0.25,
  ifelse(mono_bmmc_vs_cd34_DE$avg_log2FC > 0, "Upregulated in BMMC", "Upregulated in CD34"),
  "Nonsignificant"
)

# Create volcano plot
png("volcano_Monocytes_BMMC_vs_CD34.png", width = 1000, height = 800)
ggplot(mono_bmmc_vs_cd34_DE, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Upregulated in BMMC" = "red", "Upregulated in CD34" = "blue", "Nonsignificant" = "grey")) +
  labs(
    title = "Differential Expression: Monocytes (BMMC vs CD34)",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Significance"
  ) +
  theme_minimal()
dev.off()
#------------
#Pathway analysis
# Perform differential expression analysis (if not already done)
bmmc_vs_cd34_DE <- FindMarkers(
  object = integrated_data,
  ident.1 = "BMMC",
  ident.2 = "CD34",
  group.by = "Sample_Group",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

# Filter significant DE genes
# Perform differential expression analysis
de_results <- FindMarkers(
  integrated_data,
  ident.1 = "BMMC",
  ident.2 = "CD34",
  group.by = "Sample_Group",      # Use the Sample_Group column
  logfc.threshold = 0.25,         # Minimum log fold change
  test.use = "wilcox"             # Wilcoxon test
)

# View the results
head(de_results)

# Save the results to a file
write.csv(de_results, "DE_results_BMMC_vs_CD34.csv", row.names = TRUE)

# Filter significant genes
significant_genes <- rownames(de_results[de_results$p_val_adj < 0.05 & abs(de_results$avg_log2FC) > 0.25, ])

# Print number of significant genes
length(significant_genes)

# Save significant genes
write.csv(significant_genes, "Significant_genes_BMMC_vs_CD34.csv", row.names = FALSE)

# Set active.ident to Sample_Group
Idents(integrated_data) <- "Sample_Group"

# Confirm the change
table(Idents(integrated_data))

#Pathway analysis
# Load necessary libraries
library(enrichR)

# GO Biological Process
png("GO_Biological_Process_BMMC_vs_CD34.png", width = 1200, height = 800)
DEenrichRPlot(
  object = integrated_data,
  ident.1 = "BMMC",
  ident.2 = "CD34",
  enrich.database = "GO_Biological_Process_2023",
  assays = "RNA",
  logfc.threshold = 0.25,
  max.genes = 100,
  display.rows = 15
) +
  ggtitle("GO Biological Process: BMMC vs CD34")
dev.off()

# GO Cellular Component
png("GO_Cellular_Component_BMMC_vs_CD34.png", width = 1200, height = 800)
DEenrichRPlot(
  object = integrated_data,
  ident.1 = "BMMC",
  ident.2 = "CD34",
  enrich.database = "GO_Cellular_Component_2023",
  assays = "RNA",
  logfc.threshold = 0.25,
  max.genes = 100,
  display.rows = 15
) +
  ggtitle("GO Cellular Component: BMMC vs CD34")
dev.off()

# GO Molecular Function
png("GO_Molecular_Function_BMMC_vs_CD34.png", width = 1200, height = 800)
DEenrichRPlot(
  object = integrated_data,
  ident.1 = "BMMC",
  ident.2 = "CD34",
  enrich.database = "GO_Molecular_Function_2023",
  assays = "RNA",
  logfc.threshold = 0.25,
  max.genes = 100,
  display.rows = 15
) +
  ggtitle("GO Molecular Function: BMMC vs CD34")
dev.off()



# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)

# Prepare the gene list
# First, let's get the significant DEGs
de_results <- de_results[order(de_results$p_val_adj), ]
significant_genes <- rownames(de_results)[de_results$p_val_adj < 0.05]

# Perform GO enrichment analysis
go_enrichment <- enrichGO(gene = significant_genes,
                          OrgDb = org.Hs.eg.db,
                          keyType = "SYMBOL",
                          ont = "BP",  # Biological Process
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)

# Display top enriched pathways
head(go_enrichment@result, n=20)
# Create a visualization of the enrichment results
# Dotplot of top 20 enriched terms
png("Top 20 Enriched GO Terms.png")
dotplot(go_enrichment, showCategory=20) +
  theme_minimal() +
  ggtitle("Top 20 Enriched GO Terms")
dev.off()

# Extract the pathway with the lowest p-value
top_pathway <- go_enrichment@result[which.min(go_enrichment@result$p.adjust), ]

# Display the top pathway
print(top_pathway)

# Description GeneRatio   BgRatio RichFactor FoldEnrichment   zScore       pvalue
#GO:1903131 GO:1903131 mononuclear cell differentiation  113/1407 484/18888  0.2334711       3.134187 13.49412 1.108577e-28
#p.adjust       qvalue
#GO:1903131 6.21025e-25 4.525329e-25
#The pathway with the lowest p-value is "mononuclear cell differentiation," 
#which is statistically significant in the GO enrichment analysis. This biological process involves the development and differentiation of mononuclear cells, 
#which are critical components of the immune system, including lymphocytes and monocytes. 
#These cells play essential roles in immune response, inflammation, and tissue repair.
#---------------------week4-----------------------------
