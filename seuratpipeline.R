#setwd("D:/baylor/TRIGEMINAL/b6tgnuintegrateddim20res01annotated")
library(Seurat)
library(SeuratDisk)
library(hdf5r)
library(dplyr)
library(patchwork)
library(ggplot2)
library(mclust)
import(Scanpy)
# load data


#b6tgnuintegrateddim20res01annotated <- Convert("spinalcord_single_nuclei.h5ad", dest = "h5seurat", overwrite = TRUE)
b6tgnuintegrateddim20res01annotated <- Read10X_h5(filename = "spinalcord_single_nuclei.h5ad")
#b6tgnuintegrateddim20res01annotated <-LoadH5Seurat("D:/baylor/TRIGEMINAL/b6tgnuintegrateddim20res01annotated/Neuron_paper_cell_label/objintegrated.h5seurat")


b6tgnuintegrateddim20res01annotated <- CreateSeuratObject(counts = b6tgnuintegrateddim20res01annotated, project = "b6tgnuintegrateddim20res01annotated", min.cells = 3, min.features = 200)
b6tgnuintegrateddim20res01annotated

dense.size <- object.size(b6tgnuintegrateddim20res01annotated)
dense.size

sparse.size <- object.size(b6tgnuintegrateddim20res01annotated)
sparse.size

b6tgnuintegrateddim20res01annotated[["percent.mt"]] <- PercentageFeatureSet(b6tgnuintegrateddim20res01annotated, pattern = "^MT-|^Mt-")

# Visualize QC metrics as a violin plot
jpeg(file="b6tgnuintegrateddim20res01annotatedPercentageFeatureSet.jpeg")
VlnPlot(b6tgnuintegrateddim20res01annotated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

jpeg(file="b6tgnuintegrateddim20res01annotatedFeatureScatter.jpeg")
plot1 <- FeatureScatter(b6tgnuintegrateddim20res01annotated, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(b6tgnuintegrateddim20res01annotated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()


b6tgnuintegrateddim20res01annotated <- NormalizeData(b6tgnuintegrateddim20res01annotated, normalization.method = "LogNormalize", scale.factor = 10000)

b6tgnuintegrateddim20res01annotated <- FindVariableFeatures(b6tgnuintegrateddim20res01annotated, selection.method = "vst", nfeatures = 32000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(b6tgnuintegrateddim20res01annotated), 10)
write.csv(top10,"b6tgnuintegrateddim20res01annotatedtophighvariablegenes.csv")

top5 <- head(VariableFeatures(b6tgnuintegrateddim20res01annotated), 5)
write.csv(top5,"b6tgnuintegrateddim20res01annotatedtophighvariablegenes5.csv")

# plot variable features with and without labels
jpeg(file="b6tgnuintegrateddim20res01annotatedvariable_plot1.jpeg")
plot1 <- VariableFeaturePlot(b6tgnuintegrateddim20res01annotated)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()

# plot variable features with and without labels
jpeg(file="b6tgnuintegrateddim20res01annotatedvariable5_plot1.jpeg")
plot1 <- VariableFeaturePlot(b6tgnuintegrateddim20res01annotated)
plot2 <- LabelPoints(plot = plot1, points = top5, repel = TRUE)
plot1 + plot2
dev.off()

all.genes <- rownames(b6tgnuintegrateddim20res01annotated)
b6tgnuintegrateddim20res01annotated <- ScaleData(b6tgnuintegrateddim20res01annotated, features = all.genes)

b6tgnuintegrateddim20res01annotated <- RunPCA(b6tgnuintegrateddim20res01annotated, features = VariableFeatures(object = b6tgnuintegrateddim20res01annotated))

jpeg(file="b6tgnuintegrateddim20res01annotatedvizdim3load.jpeg")
VizDimLoadings(b6tgnuintegrateddim20res01annotated, dims = 1:2, reduction = "pca")
dev.off()

b6tgnuintegrateddim20res01annotated <- FindNeighbors(b6tgnuintegrateddim20res01annotated, dims = 1:30)
b6tgnuintegrateddim20res01annotated <- FindClusters(b6tgnuintegrateddim20res01annotated, resolution = 0.1)
b6tgnuintegrateddim20res01annotated <- RunUMAP(b6tgnuintegrateddim20res01annotated, dims = 1:30)

jpeg(file="b6tgnuintegrateddim20res01annotatedaged0.530dims.jpeg")
DimPlot(b6tgnuintegrateddim20res01annotated, group.by = "celltype",reduction = "umap")
dev.off()

jpeg(file="b6tgnuintegrateddim20res01annotatedgender0.130dims.jpeg")
DimPlot(b6tgnuintegrateddim20res01annotated, group.by = "gender",reduction = "umap")
dev.off()

jpeg(file="b6tgnuintegrateddim20res01annotateddonor0.130dims.jpeg")
DimPlot(b6tgnuintegrateddim20res01annotated, group.by = "donor",reduction = "umap")
dev.off()

jpeg(file="b6tgnuintegrateddim20res01annotatednotes0.130dims.jpeg")
DimPlot(b6tgnuintegrateddim20res01annotated, group.by = "notes",reduction = "umap")
dev.off()

jpeg(file="b6tgnuintegrateddim20res01annotateddisease0.130dims.jpeg")
DimPlot(b6tgnuintegrateddim20res01annotated, group.by = "disease",reduction = "umap")
dev.off()

jpeg(file="b6tgnuintegrateddim20res01annotatedid0.130dims.jpeg")
DimPlot(b6tgnuintegrateddim20res01annotated, group.by = "sampleid",reduction = "umap")
dev.off()

jpeg(file="b6tgnuintegrateddim20res01annotatedauthors0.130dims.jpeg")
DimPlot(b6tgnuintegrateddim20res01annotated, group.by = "foldername",reduction = "umap")
dev.off()

jpeg(file="b6tgnuintegrateddim20res01annotated0.110dims.jpeg")
DimPlot(b6tgnuintegrateddim20res01annotated,reduction = "umap")
dev.off()


b6tgnuintegrateddim20res01annotated <- FindClusters(b6tgnuintegrateddim20res01annotated, resolution = 0.1)
b6tgnuintegrateddim20res01annotated <- RunUMAP(b6tgnuintegrateddim20res01annotated, dims = 1:10)
jpeg(file="b6tgnuintegrateddim20res01annotated0.310dims.jpeg")
DimPlot(b6tgnuintegrateddim20res01annotated, reduction = "umap")
dev.off()

b6tgnuintegrateddim20res01annotated <- FindClusters(b6tgnuintegrateddim20res01annotated, resolution = 0.4)
b6tgnuintegrateddim20res01annotated <- RunUMAP(b6tgnuintegrateddim20res01annotated, dims = 1:10)
jpeg(file="b6tgnuintegrateddim20res01annotated0.410dims.jpeg")
DimPlot(b6tgnuintegrateddim20res01annotated, reduction = "umap", label=T)
dev.off()

# find all markers of cluster 2
cluster2.markers <- FindMarkers(b6tgnuintegrateddim20res01annotated, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
write.csv(cluster2.markers, "b6tgnuintegrateddim20res01annotatedmarkers.csv")

# find markers for every cluster compared to all remaining cells, report only the positive - generate the Differentially expressed genes - DEGs list
# ones
b6tgnuintegrateddim20res01annotated.markers <- FindAllMarkers(b6tgnuintegrateddim20res01annotated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
b6tgnuintegrateddim20res01annotated.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(b6tgnuintegrateddim20res01annotated.markers, "tg_combinedclustermarkersres01.csv")



#Aditional steps
# Load the Seurat package and your Seurat object
library(Seurat)
my_seurat <- Read10X(data.dir = "path/to/data")
my_seurat <- CreateSeuratObject(counts = my_seurat$`Gene Expression`)

# Identify cell types using clustering and/or manual annotation
# For example:
my_seurat <- FindNeighbors(my_seurat, dims = 1:10)
my_seurat <- FindClusters(my_seurat, resolution = 0.5)
my_seurat$cell_type <- Idents(my_seurat)

# Create the cell type proportion barplot
VlnPlot(my_seurat, features = "cell_type", pt.size = 0, group.by = "cell_type")
VlnPlot(my_seurat, features = "cell_type", pt.size = 0, group.by = "cell_type", geom = "bar")

#Feature expression heatmap

jpeg(file="objB6TGNuheatmap.jpeg")
DoHeatmap(
  b6tgnuintegrateddim20res01annotated,
  features = NULL,
  cells = NULL,
  group.by = "name",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = NULL,
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)
dev.off()

#Heatmap with a given gene list

counts <- GetAssyData(b6tgnuintegrateddim20res01annotated, assay="RNA", slot="data")
genes <- c("Iba1","Gfap","Iba1","Atf3","cFos","Fos","Il1b","Il6","Nos2","Nox","Ccl2","Cd68","Itgam","Tac1","Calca","Trpm8","Trpv1","S100b","Gfra2",
           "Pou4F2","Gal","Cd55","Scn11a","Fxyd7","Ngfr","Nefh","Hapln4","Cbln2","Kcnab1","Sst","Il31ra","Apoe","Fabp7","Mpz","Gldn","Scn7a","Dcn","Pdgfra","Mgp","Alpl",
           "Cd74","Igfbp7","Tinagl1","Htr1f","Klf6","Klf9","Mt1","Egr1","Egr2","Cyr61","Nr4a1","Ctgf","Jag1","Lrp1")
counts <- as.matrix(counts[rownames(counts) %in% genes, ])

jpeg(file="genelistB6TGNuheatmap.jpeg")
DoHeatmap(
  counts,
  features = NULL,
  cells = NULL,
  group.by = "ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = NULL,
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)


dev.off()



saveRDS(b6tgnuintegrateddim20res01annotated, file = "b6tgnuintegrateddim20res01annotated.rds")


