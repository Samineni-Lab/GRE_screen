## GRE_screen
This repository describes the R scripts used for identifying putative celltype-specific genomic regulatory elements (GREs) using snRNA-seq and snATAC-seq data.

Important packages in R:
R version 4.4.1
Seurat 5.1.0
Signac 1.13.0


Step 1: Get expression/accessibility data of differentially expressed (DE) genes and differentially accessible (DA) peaks for each cell type.

Step 2: Calculate the correlation coefficients between any DE genes to any DA peaks, and annotate peaks and genes.

Step 3: Identify putative celltype-specific enhancers and silencers.

Step 4: Generate plots.
