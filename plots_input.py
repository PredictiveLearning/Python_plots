import numpy as np
# coding: utf-8

# In[ ]:
#RedshiftsToRead = [True,True, True,True,True,True,True, True, True, True,True, True]
#RedshiftsToRead = [False, True, True,True,True,True,False,False,False,False,False,False]
#RedshiftsToRead = [True, True, True, True,True,True,True,True,True,True,False,False]
RedshiftsToRead = [True, True, True, True,True,True,False,False,False,False,False,False]
#RedshiftsToRead = [True, True, True,True,True, True, True, True,True,True,True,True,True,True,True,False,False]
#RedshiftsToRead = [True,False, False,False,False,False,False,False,False,False,False,False]
#RedshiftsToRead = [False,False,False,False,True,True,True,True,False,False,False,False]
#RedshiftsToRead = [True,False,False,False,False,True, False,False,False, False,False, False]

    
CatalogType='snap'
#CatalogType='tree'

#COSMOLOGIES & DARK MATTER SIMS
PLANCK=1
CATERPILLAR_PLANCK=0
WMAP1=0
WMAP7=0

#RUN OPTIONS
MR=1
MRII=1
opt_plot_MCMC_sample=0
opt_detailed_enrichment=1
opt_individual_elements=0
opt_rings=1
opt_rings_in_bulges=1

Datadir = '/net/bootes/export/data1/data/'
MCMCdir = '/net/bootes/export/data1/Workspace/LGal_Development_Branch/MCMC/'
#MCMCSampledir = '/net/bootes//export/data1/Workspace/LGal_Development_Branch/output/bestfit/final_'
#MCMCSampledir = '/net/bootes//export/data1/Workspace/LGal_Development_Branch/output/task0_'
#MCMCSampledir = '/net/bootes//export/data1/Workspace/LGal_Development_Branch/output/bestfit/task0_'
#MCMCSampledir = '/net/bootes//export/data1/Workspace/1/output/bestfit/final_'
#MCMCSampledir = '/net/bootes//export/data1/Workspace/cosma/task0_'
MCMCSampledir = '/net/bootes/export/data1/Workspace/cosma/final_'
#MCMCSampledir = '/net/bootes/export/data1/Workspace/1/output/task0_'

DirName_MR = '/net/bootes/scratch-ssd/SAM/test3/MR/'
DirName_MRII = '/net/bootes/scratch-ssd/SAM/test3/MRII/'
#DirName_MR = '/net/bootes/scratch-ssd/SAM/test95_sigma_pow/MR/'
#DirName_MRII = '/net/bootes/scratch-ssd/SAM/test95_sigma_pow/MRII/'
#DirName_MR = '/net/bootes/export/data1/SAM/test57/MR/'
#DirName_MRII = '/net/bootes/export/data1/SAM/test57/MRII/'
#DirName_MR = '/data/astrodata/Bruno/SAM/test57/MR/'
#DirName_MRII = '/data/astrodata/Bruno/SAM/test57/MRII/'
             
#DirName_MR = '/net/bootes/scratch-ssd/SAM/test0/MR/'
#DirName_MR = '/net/bootes/scratch2/SAM/Henriques2015a/snaps/MR/'
#DirName_MRII = '/net/bootes/scratch2/SAM/Henriques2015a/snaps/MRII/'
#DirName_MR = '/net/bootes/scratch2/SAM/Henriques2015a/GalTree/MR/'
#DirName_MRII = '/net/bootes/scratch2/SAM/Henriques2015a/GalTree/MRII/'
#DirName_MR = '/net/bootes/scratch2/SAM/Henriques2015a/snaps/nmerge_10/MR/'

#DirName_MR = '/net/bootes/scratch2/SAM/HWL18/MR/'
#DirName_MRII = '/net/bootes/scratch2/SAM/HWL18/MRII/'
#DirName_MR = '/net/bootes/scratch2/SAM/HWL18/only_SN/MR/'
#DirName_MRII = '/net/bootes/scratch2/SAM/HWL18/only_SN/MRII/'
#DirName_MR = '/net/bootes/scratch2/SAM/HWL18/only_AGN/MR/'
#DirName_MRII = '/net/bootes/scratch2/SAM/HWL18/only_AGN/MRII/'
#DirName_MR = '/net/bootes/scratch2/SAM/HWL18/no_feedback/MR/'
#DirName_MRII = '/net/bootes/scratch2/SAM/HWL18/no_feedback/MRII/'

#DirName_MR = '/net/bootes/scratch2/SAM/HWL18/ejection_cut_low_mass/MR/' 
#DirName_MR = '/net/bootes/scratch2/SAM/HWL18/ejection_cut_very_low_mass/MR/'
#DirName_MR = '/net/bootes/scratch2/SAM/HWL18/ejection_cut_high_mass/MR/' 
#DirName_MR = '/net/bootes/scratch2/SAM/HWL18/ejection_cut_very_high_mass/MR/' 
  

#DirName_MR = '/net/bootes/scratch2/SAM/Hen15_tests/Guo_reinc/MR/'
#DirName_MR = '/net/bootes/scratch2/SAM/Hen15_tests/bh_in_instabilities/MR/'

#DirName_MR = '/net/bootes/scratch2/SAM/Henriques2015a/snaps/new_BH_growth/MR/'
#DirName_MRII = '/net/bootes/scratch2/SAM/Henriques2015a/snaps/new_BH_growth/MRII/'

#DirName_MR = '/net/bootes/scratch2/SAM/Henriques2015a/snaps/MR/'
#DirName_MRII = '/net/bootes/scratch2/SAM/Henriques2015a/snaps/MRII/'

#DirName_MR = '/net/bootes/export/data1/Workspace/LGal_Development_Branch/output/'
#DirName_MR = '/net/bootes/scratch2/SAM/Henriques2015a/GalTree/MR/'
#DirName_MR = '/Users/BrunoHenriques/Desktop/Work/SAM/Henriques2015a/GalTree/MR/'


prefix_this_model='This Work'
file_this_model='ThisWork'

do_previous_model1=1
file_previous_model1=Datadir+'Guo2013a_m05'
prefix_previous_model1='Guo2013a'
linestyle_previous_model1=':'

#do_previous_model2=1
#file_previous_model2=Datadir+'Henriques2013a'
#prefix_previous_model2='Henriques2013a - WMAP7'
#linestyle_previous_model2='--'

do_previous_model2=1
file_previous_model2=Datadir+'/Henriques2014a/Henriques2014a'
prefix_previous_model2='Henriques2015'
linestyle_previous_model2='--'
    
#HENRIQUES15
slope_color_cut=[0.075,0.275, 0.3, 0.32,0.38,0.38,0.38,0.38]
offset_color_cut=[1.85,1.213, 1.18,0.99,0.79,0.79,0.79,0.79]     
 
#HENRIQUES17
#slope_color_cut=[0.05,0.275, 0.3, 0.32,0.38,  0.38,0.38,0.38]
#offset_color_cut=[1.85,1.213, 1.18,0.99,0.79,  0.79,0.79,0.79]

#NO AGN
#slope_color_cut=[0.075, 0.5, 0.48, 0.38, 0.18];
#offset_color_cut=[1.95, 1.085, 1.0, 1.0, 1.15];  

#MCMC    
#slope_color_cut=[0.075, 0.5, 0.48, 0.38, 0.18];
#offset_color_cut=[1.8, 1.085, 1.1, 1.0, 1.15];

#slope_color_cut=[0.5, 0.5, 0.48, 0.38, 0.18];
#offset_color_cut=[1.085, 1.085, 1.0, 1.0, 1.15];  
 
minimum_y_color_cut=[0.0,1.3,1.3,1.3,1.3,1.3,1.3,1.3]

MCMC_slope_color_cut=[0.075, 0.5, 0.48, 0.38, 0.18];
MCMC_offset_color_cut=[1.8, 1.085, 1.1, 1.0, 1.15];    
   
obs_offset_color_cut=[0.244,0.69,0.59,0.59,0.59]
obs_slope_color_cut=[1.9,0.88,0.88,0.88,0.88]

one_one_size_small=[5,4]
one_two_size_small=[7,4]
one_three_size_large=[10,4]
one_four_size_large=[12,4]
one_five_size_large=[14,4]

two_one_size_small=[5,6]
two_two_size_small=[7,6]
two_four_size_large=[11,6]
two_two_size_large=[7,6]

three_one_size_small=[5,8]
three_two_size_small=[7,8]


RNUM=12
#in units of kpc
RingRadius=np.array([0.029, 0.059, 0.12, 0.24, 0.47, 0.94, 1.9, 3.8, 7.5, 15., 30., 60.])

#RingRadius=np.array([0.00030, 0.00090, 0.0027, 0.0081, 0.024, 0.073, 0.22, 0.66, 1.97, 5.91, 17.72, 53.14])
#RingRadius=np.array([0.0025, 0.0063, 0.016, 0.039, 0.098, 0.24, 0.61, 1.53, 3.82, 9.54, 23.84, 59.61])
#RingRadius=np.array([0.020, 0.040, 0.080, 0.16, 0.32, 0.64, 1.28, 2.56, 5.12, 10.24, 20.48, 40.96])
#RingRadius=np.array([0.025, 0.063, 0.16, 0.39, 0.98, 2.4, 6.1, 15., 38., 95., 240., 600.])

#RingRadius[i]= 0.44*pow(1.5,i+1)/1000.; 
#RingRadius=np.array([0.660,0.990,1.485,2.228,3.341,5.012,7.518,11.277,16.915,25.373,38.059,57.088])
#RingRadius[i]= 0.1*pow(1.6,i+1)/1000.;
#RingRadius=np.array([0.160,0.256,0.410,0.655,1.049,1.678,2.684,4.295,6.872,10.995,17.592,28.147])
 
#RingRadius[i]= 0.1*pow(1.7,i+1)/1000.;
#RingRadius=np.array([0.170,0.289,0.491,0.835,1.420,2.414,4.103,6.976,11.859,20.160,34.272,58.262])

#RingRadius[i]= 0.2*pow(1.6,i+1)/1000.;
#RingRadius=np.array([0.320,0.512,0.819,1.311,2.097,3.355,5.369,8.590,13.744,21.990,35.184,56.295])
# RingRadius[i]= 0.34*pow(1.52,i+1)/1000.; 
#RingRadius=np.array([0.517,0.786,1.194,1.815,2.759,4.193,6.374,9.688,14.726,22.383,34.022,51.713])


SDSS_min_z=0.005     
SDSS_max_z=0.2

Hubble_h_WMAP1 = 0.732
Hubble_h_WMAP7 = 0.704

    
if WMAP1: 
    FullRedshiftList=[0.00,0.09,0.41,0.99,2.07,3.06,3.87] 
    FullSnapshotList_MR=[63,59,50,41,32,27,24] 
    FullSnapshotList_MRII=[67,63,54,45,36,31,28]
    BoxSize_MR    = 500. #full MR 
    BoxSize_MRII  = 100. #full MRII      
    Hubble_h      = 0.732
    Omega_M       = 0.25 
    Omega_Lambda  = 0.75
    MaxTreeFiles  = 512
    
if WMAP7: 
    FullRedshiftList=[0.00,0.08,0.39,1.02,1.92,2.92,4.05] 
    FullSnapshotList_MR=[53,51,45,37,30, 25,21]     
    FullSnapshotList_MRII=[57,55,49,41,34, 29,25]
    BoxSize_MR    = 521.555 #full MR 
    BoxSize_MRII  = 104.3 #full MRII      
    Hubble_h      = 0.704
    Omega_M       = 0.272 
    Omega_Lambda  = 0.728
    MaxTreeFiles  = 512
 


if PLANCK:      
   
    FullRedshiftList=[0.00,0.11,0.40,1.04,2.07,3.11,3.95,5.03,5.92,6.97,8.22,8.93] 
    FullSnapshotList_MR=[58,54,47,38,30,25,22,19,17,15,13,12]
    #FullRedshiftList=[0.00,0.11,0.26,0.40,0.51,0.57,0.71,0.86,1.04,2.07,3.11,3.95,5.03,5.92,6.97,8.22,8.93] 
    #FullSnapshotList_MR=[58,54,50,47,45,44,42,40,38,30,25,22,19,17,15,13,12]
    #FullRedshiftList=[1.04,1.48,2.07] 
    #FullSnapshotList_MR=[38,34,30]
    #FullRedshiftList=[0.00,0.51,1.04,2.07,3.11,5.03] 
    #FullSnapshotList_MR=[58,45,38,30,25,19] 
    
    FullSnapshotList_MRII=[62,58,51,42,34,29,26,23,21,19,17,16]
    BoxSize_MR    = 500.* 0.960558 #full MR 
    BoxSize_MRII  = 100.* 0.960558 #full MRII      
    Hubble_h      = 0.673
    Omega_M       = 0.315 
    Omega_Lambda  = 0.683
    MaxTreeFiles  = 512
    
if CATERPILLAR_PLANCK:
    FullRedshiftList=[0.00,0.10,0.40,1.00,2.01,2.98,3.96]
    FullSnapshotList=[319,295,245,189,145,124,111]
    BoxSize_MR    = 100.
    BoxSize_MRII    = 100.
    Hubble_h      = 0.673
    Omega_M       = 0.315
    Omega_Lambda  = 0.683
    MaxTreeFiles  = 8
    
    
TNG100_color =  '#1F77B4'
TNG300_color =  '#FF7F0E'
Illustris_color =  '#D62728'
Magneticum_Box4_color = '#0000FF'
Magneticum_Box2_color = '#007F7F'
LGalaxies_color =  '#FF00FF'
EAGLE_color = '#00FFFF'
    
