library(Signac)
library(Seurat)
library(data.table)

setwd("E:/GRE_screen")

# read seurat object
seurat_mat=readRDS("Seurat.Rds")

# group cells by celltype
Idents(seurat_mat)='celltype'

# calculate average expression by celltype
seurat_mat@active.assay='RNA'
matrix_gene <- AverageExpression(seurat_mat)
saveRDS(matrix_gene,'S1_averageexpression_gene_celltype.Rds')

seurat_mat@active.assay='ATAC'
matrix_peak <- AverageExpression(seurat_mat)
saveRDS(matrix_peak,'S1_averageexpression_peak_celltype.Rds')

# identify DE genes and DA peaks
seurat_mat@active.assay='RNA'
markers_gene=FindAllMarkers(seurat_mat)
fwrite(markers_gene,'S1_markers_gene_celltype.csv')

seurat_mat@active.assay='ATAC'
markers_peak=FindAllMarkers(seurat_mat)
fwrite(markers_gene,'S1_markers_peak_celltype.csv')
