# add cellranger-arc to PATH
export PATH=/ref/software/cellranger-6.1.2:$PATH


cd /ref


# prepare files for custome genome by using the Projection-TAG.fa and Projection-TAG.gtf
cp -r refdata-cellranger-arc-mm10-2020-A-2.0.0 refdata-cellranger-arc-mm10-2020-A-2.0.0-custom

cd refdata-cellranger-arc-mm10-2020-A-2.0.0-custom

# copy and update fa file
cp fasta/genome.fa fasta/genome.custom.fa 
cat Projection-TAG.fa >> fasta/genome.custom.fa 
grep "^>" fasta/genome.custom.fa 

# copy and update gtf file
cp genes/genes.gtf genes/genes.custom.gtf
cat Projection-TAG.gtf >> genes/genes.custom.gtf
tail genes/genes.custom.gtf

# build reference genome
cd /ref
cellranger mkref --genome=refdata-cellranger-arc-mm10-2020-A-2.0.0-Projection-TAG --fasta=/ref/refdata-cellranger-arc-mm10-2020-A-2.0.0-custom/fasta/genome.custom.fa --genes=/ref/refdata-cellranger-arc-mm10-2020-A-2.0.0-custom/genes/genes.custom.gtf






