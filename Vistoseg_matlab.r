#2nd_Step
#Using matlab
#Section: 2.5.2 Histology image and segmentation (Vistoseg)
#code adapted from: http://research.libd.org/VistoSeg/

fname='/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_1/raw_data/Visto/4_year_old_T2.tif';
N=5;
VNS('/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_1/raw_data/Visto/4_year_old_T2.tif', 5)


fname='/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_1/raw_data/Visto/4_year_old_T2.tif';
M=1;
refineVNS('/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_1/raw_data/Visto/4_year_old_T2.tif', 1)

mask='/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_1/raw_data/Visto/4_year_old_T2_nuclei.mat';
jsonname='/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_1/raw_data/run_spaceranger_count_DH1/outs/spatial/scalefactors_json.json';
posname='/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_1/raw_data/run_spaceranger_count_DH1/outs/spatial/tissue_positions_list.csv';
countNuclei(mask,jsonname,posname)


fname='/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_2/raw_data/Visto/4_year_old_T1.tif';
N=5;
VNS(fname,N)


fname='/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_2/raw_data/Visto/4_year_old_T1.tif';
M=5;
refineVNS(fname,M)

mask='/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_2/raw_data/Visto/4_year_old_T1_nuclei.mat';
jsonname='/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_2/raw_data/run_spaceranger_count_DH2/outs/spatial/scalefactors_json.json';
posname='/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_2/raw_data/run_spaceranger_count_DH2/outs/spatial/tissue_positions_list.csv';
countNuclei(mask,jsonname,posname)

fname='/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_2/raw_data/Visto/15_year_old_T1.tif';
N=5;
VNS(fname,N)


fname='/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_2/raw_data/Visto/15_year_old_T1.tif';
M=4;
refineVNS(fname,M)

mask='/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_2/raw_data/Visto/15_year_old_T1_nuclei.mat';
jsonname='/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_2/raw_data/run_spaceranger_count_DH4/outs/spatial/scalefactors_json.json';
posname='/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_2/raw_data/run_spaceranger_count_DH4/outs/spatial/tissue_positions_list.csv';
countNuclei(mask,jsonname,posname)

fname='/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_3/raw_data/Visto/15_year_old_T2.tif';
N=5;
VNS(fname,N)


fname='/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_3/raw_data/Visto/15_year_old_T2.tif';
M=2;
refineVNS(fname,M)

mask='/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_3/raw_data/Visto/15_year_old_T2_nuclei.mat';
jsonname='/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_3/raw_data/run_spaceranger_count_DH3/outs/spatial/scalefactors_json.json';
posname='/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_3/raw_data/run_spaceranger_count_DH3/outs/spatial/tissue_positions_list.csv';
countNuclei(mask,jsonname,posname)


fname='/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_01/raw_data/Visto/15_yr_old_DH1.tif';
N=5;
VNS(fname,N)


fname='/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_01/raw_data/Visto/15_yr_old_DH1.tif';
M=5;
refineVNS(fname,M)

mask='/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_01/raw_data/Visto/15_yr_old_DH1_nuclei.mat';
jsonname='/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_01/raw_data/run_spaceranger_count_DH1a/outs/spatial/scalefactors_json.json';
posname='/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_01/raw_data/run_spaceranger_count_DH1a/outs/spatial/tissue_positions_list.csv';
countNuclei(mask,jsonname,posname)
#Deleted in run spaceranger DH1a but it is saved in vistoseg in gui

fname='/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_01/raw_data/Visto/15_yr_old_DH2.tif';
N=5;
VNS(fname,N)


fname='/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_01/raw_data/Visto/15_yr_old_DH2.tif';
M=1;
refineVNS(fname,M)

mask='/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_01/raw_data/Visto/15_yr_old_DH2_nuclei.mat';
jsonname='/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_01/raw_data/run_spaceranger_count_DH2a/outs/spatial/scalefactors_json.json';
posname='/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_01/raw_data/run_spaceranger_count_DH2a/outs/spatial/tissue_positions_list.csv';
countNuclei(mask,jsonname,posname)

fname='/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_02/raw_data/Visto/31_yr_old_DH3.tif';
N=5;
VNS(fname,N)


fname='/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_02/raw_data/Visto/31_yr_old_DH3.tif';
M=5;
refineVNS(fname,M)

mask='/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_02/raw_data/Visto/31_yr_old_DH3_nuclei.mat';
jsonname='/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_02/raw_data/run_spaceranger_count_DH3a/outs/spatial/scalefactors_json.json';
posname='/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_02/raw_data/run_spaceranger_count_DH3a/outs/spatial/tissue_positions_list.csv';
countNuclei(mask,jsonname,posname)

fname='/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_02/raw_data/Visto/31_yr_old_DH4.tif';
N=5;
VNS(fname,N)


fname='/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_02/raw_data/Visto/31_yr_old_DH4.tif';
M=4;
refineVNS(fname,M)

mask='/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_02/raw_data/Visto/31_yr_old_DH4_nuclei.mat';
jsonname='/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_02/raw_data/run_spaceranger_count_DH4a/outs/spatial/scalefactors_json.json';
posname='/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_02/raw_data/run_spaceranger_count_DH4a/outs/spatial/tissue_positions_list.csv';
countNuclei(mask,jsonname,posname)
