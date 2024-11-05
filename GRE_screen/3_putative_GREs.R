library(dplyr)
library(tibble)
library(data.table)

setwd("E:/GRE_screen")

# load FindMarker results and correlation table
markers_gene=fread("S1_markers_gene_celltype.csv")
markers_peak=fread("S1_markers_peak_celltype.csv")
peak_gene_correlation=fread("S2_peak_gene_correlation_distance_TSS.csv")

# Select and rename columns
markers_DE_gene <- markers_gene %>%
  select(-p_val)%>%
  rename(avg_log2FC_RNA = avg_log2FC,
         pct.1_RNA = pct.1,
         pct.2_RNA = pct.2,
         p_val_adj_RNA = p_val_adj,
         cluster_RNA = cluster)

markers_DA_peak <- markers_peak %>%
  select(-p_val)%>%
  rename(avg_log2FC_ATAC = avg_log2FC,
         pct.1_ATAC = pct.1,
         pct.2_ATAC = pct.2,
         p_val_adj_ATAC = p_val_adj,
         cluster_ATAC = cluster,
         peak = gene)

# join tables
final_df <- peak_gene_correlation%>%
  inner_join(markers_DE_gene, by = c("gene" = "gene"))%>%
  inner_join(markers_DA_peak, by = c("peak" = "peak"))


# Identify putative enhancers and regulated genes
enhancer_df <- final_df%>%
  filter(
    avg_log2FC_ATAC > 0.5, p_val_adj_ATAC < 0.05, # peak criteria
    avg_log2FC_RNA > 0.5, p_val_adj_RNA < 0.05, # gene criteria
    cor > 0.75,cluster_ATAC==cluster_RNA,abs(distance)<5e6
  )

fwrite(enhancer_df,file='S3_celltype_puEnhancer_gene.csv')
    
# Identify putative silencers and regulated genes
enhancer_df <- final_df%>%
  filter(
    avg_log2FC_ATAC > 0.5, p_val_adj_ATAC < 0.05, # peak criteria
    avg_log2FC_RNA < (-0.5), p_val_adj_RNA < 0.05, # gene criteria
    cor < (-0.75),cluster_ATAC==cluster_RNA,abs(distance)<5e6
  )

fwrite(enhancer_df,file='S3_celltype_puSilencer_gene.csv')
