
# coding: utf-8

# In[ ]:
RedshiftsToRead = [True,True,True,True,True,True,False]
#RedshiftsToRead = [True,True,True,True,True,True]

CatalogType='snap'
#CatalogType='tree'

#COSMOLOGIES & DARK MATTER SIMS
PLANCK=1
CATERPILLAR_PLANCK=0
WMAP1=0
WMAP7=0

#RUN OPTIONS
MRII=0
opt_plot_MCMC_sample=0
opt_detailed_enrichment=0
opt_rings=0

#PLOT OPTIONS
opt_stellar_mass_vs_halo_mass=0
opt_stellar_mass_function=0
opt_redfraction_color_cut_cuts=0
opt_redfraction_color_cut=0
opt_metals_vs_stellarmass=0
opt_gasmetals_vs_stellarmass=0
opt_BHBM=1
opt_SFRF=0
opt_gas_fraction=0
opt_HI_fraction=0
opt_HI_MF=0
opt_morphology_vs_stellarmass=0

opt_sfr_vs_stellar_mass=0
opt_ur_vs_r=0
opt_UVJ_colour=0
opt_UVJ_grid=0
opt_sizes_vs_stellarmass=0
    
#options for H2_AND_RINGS
opt_milkyway_sfr_and_gas_profiles=0
opt_gas_metallicity_gradients=0
opt_SFR_gradients=0
opt_gasfractions_vs_stellarmass=0
opt_H2fraction_vs_stellarmass=0

opt_evo_milkyway_gas_profile=0
opt_test_H2_prescriptions=0

#MISC
opt_BHBM_by_sfr=1
opt_AGN_quenching=1
opt_bluck_red_fractions=0
opt_sat_fraction=0
opt_HotGas_fraction=0
opt_BHmass_in_radio=0
opt_fabian_fb=0

#TESTS
opt_misc_plots=0
opt_test_resolution=0
opt_test_rings=0   
    
#when tree on
opt_simple_tree_map=0
opt_full_tree_map=0

#ANIMATIONS
opt_anime_mass_gr=0


Datadir = '/net/bootes/export/data1/data/'
MCMCdir = '/net/bootes//export/data1/Workspace/LGal_Development_Branch/MCMC/'
MCMCSampledir = '/net/bootes//export/data1/Workspace/LGal_Development_Branch/output/bestfit/'

#DirName_MR = '/net/bootes/export/data1/SAM/test5/MR/'
#DirName_MRII = '/net/bootes/export/data1/SAM/test5/MRII/'

#DirName_MR = '/net/bootes/scratch2/SAM/Henriques2015a/snaps/new_BH_growth/MR/'
#DirName_MRII = '/net/bootes/scratch2/SAM/Henriques2015a/snaps/new_BH_growth/MRII/'

DirName_MR = '/net/bootes/scratch2/SAM/Henriques2015a/snaps/MR/'
DirName_MRII = '/net/bootes/scratch2/SAM/Henriques2015a/snaps/MRII/'

#DirName_MR = '/net/bootes/scratch2/SAM/Henriques2015a/snaps/new_MR_fb_plus_001/'
#DirName_MR = '/net/bootes/scratch2/SAM/Henriques2015a/snaps/new_MR_fb_minus_001/'

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
slope_color_cut=[0.075,0.275, 0.3, 0.32,0.38]
offset_color_cut=[1.85,1.213, 1.18,0.99,0.79]     
 
#HENRIQUES16
#slope_color_cut=[0.075,0.275, 0.3, 0.32,0.38]
#offset_color_cut=[1.92,1.213, 1.18,0.99,0.79]

#slope_color_cut=[0.075, 0.5, 0.48, 0.38, 0.18];
#offset_color_cut=[2.05, 1.085, 1.0, 1.0, 1.15];    
 
minimum_y_color_cut=[0.0,1.3,1.3,1.3,1.3]
 
MCMC_slope_color_cut=[0.075, 0.5, 0.48, 0.38, 0.18];
MCMC_offset_color_cut=[1.8, 1.085, 1.1, 1.0, 1.15];    
   
obs_offset_color_cut=[0.244,0.69,0.59,0.59,0.59]
obs_slope_color_cut=[1.9,0.88,0.88,0.88,0.88]

Hubble_h_WMAP1 = 0.732
Hubble_h_WMAP7 = 0.704

    
if WMAP1: 
    FullRedshiftList=[0.00,0.41,0.99,2.07,3.06,3.87] 
    FullSnapshotList_MR=[63,50,41,32,27,24] 
    FullSnapshotList_MRII=[67,54,45,36,31,28]
    BoxSize_MR    = 500. #full MR 
    BoxSize_MRII  = 100. #full MRII      
    Hubble_h      = 0.73
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
    FullRedshiftList=[0.00,0.11,0.40,1.04,2.07,3.11,3.95] 
    FullSnapshotList_MR=[58,54,47,38,30,25,22]
    #FullRedshiftList=[1.04,1.48,2.07] 
    #FullSnapshotList_MR=[38,34,30]
    #FullRedshiftList=[0.00,0.51,1.04,2.07,3.11,5.03] 
    #FullSnapshotList_MR=[58,45,38,30,25,19] 
    
    FullSnapshotList_MRII=[62,58,51,42,34,29,26]
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
