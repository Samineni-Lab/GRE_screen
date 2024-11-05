library(dplyr)
library(tibble)
library(biomaRt)
library(GenomicRanges)
library(data.table)

setwd("E:/GRE_screen")

# Identify DE genes with cutoffs
# Positive DE genes (avg_log2FC>0.5) will be used for putative enhancers
# negative DE genes (avg_log2FC<-0.5) will be used for putative silencers
markers_gene=fread("S1_markers_gene_celltype.csv")
DE_gene=markers_gene%>%filter(abs(avg_log2FC)>0.5,p_val_adj<0.05)%>%pull(gene)%>%unique()

# filter average expression matrix for DE genes
matrix_gene=readRDS('S1_averageexpression_gene_celltype.Rds')
matrix_DE_gene=matrix_gene[DE_gene,]


# Identify DA peaks with cutoffs
markers_peak=fread("S1_markers_peak_celltype.csv")
DA_peak=markers_peak%>%filter(avg_log2FC>0.5,p_val_adj<0.05)%>%pull(gene)%>%unique()

# filter average expression matrix for DA peaks
matrix_peak=readRDS('S1_averageexpression_peak_celltype.Rds')
matrix_DA_peak=matrix_peak[DA_peak,]

# caluclate coorelation between any DE genes and any DA peaks
correlation=cor(matrix_DE_gene,matrix_DA_peak)%>%round(3)
saveRDS(correlation,'S2_correlation_DE_gene_DA_peak.Rds')

# reformat outputs
correlation.melt=correlation%>%
  as.data.frame()%>%
  rownames_to_column('peak')%>%
  melt(id.vars='peak',measure.vars=colnames(correlation))%>%
  dplyr::rename(gene=variable,cor=value)

# Get peak information
peak=stringr::str_split(rownames(correlation),'-',simplify=T)%>%
  as_tibble()%>%
  dplyr::rename(chromosome=V1,peak.start=V2,peak.end=V3)%>%
  mutate(peak=paste(chromosome,peak.start,peak.end,sep='-'),peak.start=as.numeric(peak.start),peak.end=as.numeric(peak.end))%>%
  mutate(peak.center=(peak.start+peak.end)/2)%>%
  dplyr::select(peak,chromosome,peak.start,peak.end,peak.center)

# Get gene information
mouse=useMart("ensembl", dataset = "mmusculus_gene_ensembl")
gene <- getBM(attributes = c("mgi_symbol","chromosome_name",'start_position','end_position','transcription_start_site','strand'), 
              mart =  mouse)
gene=gene%>%
  dplyr::select(gene=mgi_symbol,chromosome=chromosome_name,TSS=transcription_start_site,strand)%>%
  mutate(chromosome=paste0('chr',chromosome))%>%
  group_by(gene,chromosome,strand)%>%
  summarise(TSS=round(median(TSS)))%>%
  filter(gene %in% colnames(correlation))%>%
  dplyr::select(gene,chromosome,TSS,gene.strand=strand)

# join data
correlation.merge=correlation.melt%>%
  inner_join(peak)%>%
  inner_join(gene) # join by gene and chromosome, so that only peak-gene pair on the same chromosome will be kept

# calculate distance
correlation.merge$distance=ifelse(correlation.merge$gene.strand=='+',
                                  correlation.merge$peak.center-correlation.merge$TSS,
                                  correlation.merge$TSS-correlation.merge$peak.center)

fwrite(correlation.merge,file='S2_peak_gene_correlation_distance_TSS.csv')


