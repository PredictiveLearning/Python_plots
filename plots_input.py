
# coding: utf-8

# In[ ]:


opt_stellar_mass_function=0
opt_bluck_red_fractions=0
opt_positions=1
opt_test_plots=0

opt_anime_mass_gr=0

RedshiftsToRead = [True,True,True,True,True,False]

#CatalogType='snap'
CatalogType='tree'

PLANCK=1
if PLANCK: 
    RedshiftList=[0.00,0.40,1.04,2.07,3.11,3.95]   
    BoxSize_MR    = 500.* 0.960558 #full MR 
    BoxSize_MRII  = 100.* 0.960558 #full MRII      
    Hubble_h      = 0.673
    Omega_M       = 0.315 
    Omega_Lambda  = 0.683
    MaxTreeFiles  = 512
    
Datadir = '/net/bootes/export/data1/data/'

#DirName_MR = '/net/bootes/scratch2/SAM/test1/MR/'
#DirName_MRII = '/net/bootes/scratch2/SAM/test1/MRII/'

#DirName_MR = '/net/bootes/scratch2/SAM/Henriques2015a/snaps/MR/'
#DirName_MRII = '/net/bootes/scratch2/SAM/Henriques2015a/snaps/MRII/'

#DirName_MR = '/net/bootes/scratch2/SAM/Henriques2015a/GalTree/MR/'
DirName_MR = '/Users/BrunoHenriques/Desktop/Work/SAM/Henriques2015a/GalTree/MR/'


