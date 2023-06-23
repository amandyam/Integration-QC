library(Seurat)
#library(SeuratData)
#importscanpy
#data <- ('/storage/chentemp/u249905/seurat/drgfilter.h5ad')
data <-'/storage/chentemp/u249905/seurat/drgfilter.h5ad'
#seuratobj <- CreateSeuratObject(counts=data, project="/storage/chentemp/u249905/seurat/drgfilter.h5ad")
sceasy::convertFormat(data, from="anndata", to="seurat", outFile='drgfilter.rds')
#write10X_h5(seurat_obj, filename = "h5_file.h5")




