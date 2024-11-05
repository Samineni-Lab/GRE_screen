# add cellranger to PATH
export PATH=/ref/rglab/software/cellranger-6.1.2:$PATH
which cellranger

cd /data

#============ count GEX library ============#
# get reference genome
reference=/ref/refdata-cellranger-arc-mm10-2020-A-2.0.0-Projection-TAG
sample='20230418_mouse_SSC-16_SSC_BC_RNA'

cellranger count --id="$sample" --sample="$sample" --transcriptome="$reference" --include-introns=TRUE --fastqs=/data/fastq/"$sample" --expect-cells=5000"
