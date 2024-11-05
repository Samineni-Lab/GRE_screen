## Mapping Projection-TAG reads
This section describes the command line scripts used for mapping Projection-TAG reads using cellranger.

Step 1: Generate a custom mm10 reference genome including Projection-TAG features.

Step 2: Map the sequencing reads to the custom genome.


## Plots
This section includes the R scripts for generating key plots in the manuscript.


## Putative GRE screen
This section describes the R scripts used for identifying putative celltype-specific genomic regulatory elements (GREs) using snRNA-seq and snATAC-seq data.

Step 1: Get expression/accessibility data of differentially expressed (DE) genes and differentially accessible (DA) peaks for each cell type.

Step 2: Calculate the correlation coefficients between any DE genes to any DA peaks, and annotate peaks and genes.

Step 3: Identify putative celltype-specific enhancers and silencers.

Step 4: Generate plots.
