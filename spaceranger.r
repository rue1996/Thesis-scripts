#First_step
#code adapted from 10X Genomics (https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/using/count)

module load software/spaceranger-1.3.0
spaceranger count --id=run_spaceranger_count_DH1 \
--transcriptome=/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_01/raw_data/refdata-gex-GRCh38-2020-A \
--fastqs=/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_01/raw_data/DH1 \
--sample=DH1-SCI7T001-SCI5T001_HK33HDSX2 \
--image=/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_01/raw_data/15_yr_old_DH1.tiff \
--slide=V10Y04-037 \
--area=D1 \
--localcores=8 \
--localmem=64
