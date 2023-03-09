#step 1
#Nuclei Segmentation

module load software/matlab-R2019b
fname='/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_01/raw_data/15_yr_old_DH1.tiff';
N=5;
matlab -nodisplay -nosplash -nodesktop -r "VNS('/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_01_/raw_data/15_yr_old_DH1.tiff, 2)"

#step 2
#Refining Segmentation

module load software/matlab-R2019b
fname='/scratch/mshruv003/spatial_tutorial/spatial/VistoSeg/code/Lieber_Institute_OTS-20-7690_rush_anterior_A1.tif';
M=3;
matlab -nodisplay -nosplash -nodesktop -r "refineVNS('/scratch/mshruv003/spatial_tutorial/spatial/VistoSeg/code/Lieber_Institute_OTS-20-7690_rush_anterior_A1.tif', 3)";

#step 3
#Obtaining the nuclei count

module load software/matlab-R2020b
mask='/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_01/raw_data/15_yr_old_DH1._nuclei.mat';
jsonname='/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_01/raw_data/run_spaceranger_count_DH1/outs/spatial/scalefactors_json.json';
posname='/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_01/raw_data/run_spaceranger_count_DH1/outs/spatial/tissue_positions_list.csv';
matlab -nodisplay -nosplash -nodesktop -r "countNuclei('/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_01/raw_data/15_yr_old_DH1._nuclei.mat', '/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_01/raw_data/run_spaceranger_count_DH1/outs/spatial/scalefactors_json.json', '/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_01/raw_data/run_spaceranger_count_DH1/outs/spatial/tissue_positions_list.csv')"


