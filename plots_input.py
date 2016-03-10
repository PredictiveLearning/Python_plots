
# coding: utf-8

# In[ ]:
RedshiftsToRead = [True,True,True,True,True,True,False]

CatalogType='snap'
#CatalogType='tree'

#COSMOLOGIES & DARK MATTER SIMS
PLANCK=1
CATERPILLAR_PLANCK=0
WMAP1=0

#RUN OPTIONS
MRII=1
opt_plot_MCMC_sample=0
opt_detailed_enrichment=0

#PLOT OPTIONS
opt_stellar_mass_vs_halo_mass=0
opt_stellar_mass_function=1
opt_metals_vs_stellarmass=0
opt_BHBM=0
opt_SFRF=1
opt_gas_fraction=0
opt_HI_MF=1
opt_sfr_vs_stellar_mass=0
opt_ur_vs_r=0
opt_UVJ_colour=0
opt_redfraction_color_cut=0
    
opt_gas_metallicity_gradients=0
opt_SFR_gradients=0

opt_bluck_red_fractions=0
opt_sat_fraction=0
opt_HotGas_fraction=0
opt_BHmass_in_radio=0

opt_misc_plots=0
opt_test_resolution_rings=1
    
#when tree on
opt_simple_tree_map=0
opt_full_tree_map=0

#ANIMATIONS
opt_anime_mass_gr=0




Datadir = '/net/bootes/export/data1/data/'
MCMCdir = '/net/bootes//export/data1/Workspace/LGal_Development_Branch/MCMC/'
MCMCSampledir = '/net/bootes//export/data1/Workspace/LGal_Development_Branch/output/'

DirName_MR = '/net/bootes/export/data1/SAM/test2/MR/'
DirName_MRII = '/net/bootes/export/data1/SAM/test2/MRII/'

#DirName_MR = '/net/bootes/scratch2/SAM/test2/MR/'
#DirName_MRII = '/net/bootes/scratch2/SAM/test2/MRII/'

#DirName_MR = '/net/bootes/scratch2/SAM/Henriques2015a/snaps/MR/'
#DirName_MRII = '/net/bootes/scratch2/SAM/Henriques2015a/snaps/MRII/'

#DirName_MR = '/net/bootes/export/data1/Workspace/LGal_Development_Branch/output/'
#DirName_MR = '/net/bootes/scratch2/SAM/Henriques2015a/GalTree/MR/'
#DirName_MR = '/Users/BrunoHenriques/Desktop/Work/SAM/Henriques2015a/GalTree/MR/'



prefix_this_model='This Work - PLANCK1'
file_this_model='ThisWork'

do_previous_model1=1
file_previous_model1=Datadir+'Guo2013a_m05'
prefix_previous_model1='Guo2013a - WMAP7'
linestyle_previous_model1=':'

#do_previous_model2=1
#file_previous_model2=Datadir+'Henriques2013a'
#prefix_previous_model2='Henriques2013a - WMAP7'
#linestyle_previous_model2='--'

do_previous_model2=1
file_previous_model2=Datadir+'/Henriques2014a/Henriques2014a'
prefix_previous_model2='Henriques2015 - PLANCK1'
linestyle_previous_model2='--'


slope_red_fraction=[0.075,0.275, 0.3, 0.32,0.38]
offset_red_fraction=[1.85,1.213, 1.18,0.99,0.79]
minimum_y_red_fraction=[0.0,1.3,1.3,1.3,1.3]



Hubble_h_WMAP1 = 0.732
Hubble_h_WMAP7 = 0.704

    
if WMAP1: 
    FullRedshiftList=[0.00,0.41,0.99,2.07,3.06,3.87] 
    FullSnapshotList=[63,50,41,32,27,24]  
    BoxSize_MR    = 500. #full MR 
    BoxSize_MRII  = 100. #full MRII      
    Hubble_h      = 0.73
    Omega_M       = 0.25 
    Omega_Lambda  = 0.75
    MaxTreeFiles  = 512
    

if PLANCK: 
    FullRedshiftList=[0.00,0.11,0.40,1.04,2.07,3.11,3.95] 
    FullSnapshotList=[58,53, 47,38,30,25,22]  
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
