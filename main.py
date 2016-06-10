
# coding: utf-8

# # Python Plots for LGalaxies

# ## Import Libraries and Read Catalogs

# <p>Use functions read_snap or read_tree to read catalogs. These are both defined in procedures.py. In case of read_snap, SnapshotList will be returned containing the list of snapshots read (usefull to later select galaxies in a given redshift).<p>

# In[148]:

import numpy as np
get_ipython().magic('matplotlib inline')

import pandas as pd

get_ipython().magic('pylab inline')
#import seaborn as sns
#sns.set_style('darkgrid')

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.table import Table
from astropy.io import fits
from importlib import reload
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
from decimal import *
import sys
from scipy.ndimage import zoom 
from mpl_toolkits.axes_grid.inset_locator import inset_axes

import procedures
reload (procedures)
from procedures import *
import plots_input
reload (plots_input)
from plots_input import *

FirstFile = 40
LastFile =  49

Volume_MR = (BoxSize_MR**3.0) * (LastFile - FirstFile + 1) / MaxTreeFiles 
Volume_MRII = (BoxSize_MRII**3.0) * (LastFile - FirstFile + 1) / MaxTreeFiles 

print('Reading started')

if CatalogType=='snap':       
    from LGalaxies_Henriques2015a_struct import LGalaxiesStruct
    from LGalaxies_Henriques2015a_struct import PropertiesToRead
    #from LGalaxies_Henriques2015a_metals_struct import LGalaxiesStruct
    #from LGalaxies_Henriques2015a_metals_struct import PropertiesToRead
    #from LGalaxies_Henriques2015a_Elements_struct import LGalaxiesStruct
    #from LGalaxies_Henriques2015a_Elements_struct import PropertiesToRead   
    #from LGalaxies_fu13_Rings_struct import LGalaxiesStruct
    #from LGalaxies_fu13_Rings_struct import PropertiesToRead
    #from LGalaxies_Henriques2015a_Rings_struct import LGalaxiesStruct
    #from LGalaxies_Henriques2015a_Rings_struct import PropertiesToRead     
    #from LGalaxies_Henriques2015a_Caterpillar_struct import LGalaxiesStruct
    #from LGalaxies_Henriques2015a_Caterpillar_struct import PropertiesToRead
    #from LGalaxies_Henriques2015a_Elements_Rings_struct import LGalaxiesStruct
    #from LGalaxies_Henriques2015a_Elements_Rings_struct import PropertiesToRead
    
    print('\n\nDoing MR')
    (G_MR, SnapshotList_MR) = read_snap(DirName_MR,FirstFile,LastFile,
                     PropertiesToRead,LGalaxiesStruct,RedshiftsToRead,FullRedshiftList)
    if(MRII==1):
        print('\n\nDoing MRII')
        (G_MRII, SnapshotList_MRII) = read_snap(DirName_MRII,FirstFile,LastFile,
                         PropertiesToRead,LGalaxiesStruct,RedshiftsToRead,FullRedshiftList)
if CatalogType=='tree':    
    from LGalaxies_tree_Henriques2015a_struct import LGalaxiesStruct
    from LGalaxies_tree_Henriques2015a_struct import PropertiesToRead_tree    
    #from LGalaxies_tree_Henriques2015a_Rings_struct import LGalaxiesStruct
    #from LGalaxies_tree_Henriques2015a_Rings_struct import PropertiesToRead_tree
    #from LGalaxies_tree_Henriques2015a_Elements_Rings_struct import LGalaxiesStruct
    #from LGalaxies_tree_Henriques2015a_Elements_Rings_struct import PropertiesToRead_tree
    #from LGalaxies_tree_Henriques2015a_Caterpillar_struct import LGalaxiesStruct
    #from LGalaxies_tree_Henriques2015a_Caterpillar_struct import PropertiesToRead_tree   
    
    (G_MR) = read_tree(DirName_MR,FirstFile,LastFile,
                     PropertiesToRead_tree,LGalaxiesStruct)
    SnapshotList = np.zeros(len(FullRedshiftList),dtype=np.int32)
    
    if(MRII==1):
        (G_MRII) = read_tree(DirName_MRII,FirstFile,LastFile,
                     PropertiesToRead_tree,LGalaxiesStruct)    
        SnapshotList_MRII = np.zeros(len(FullRedshiftList),dtype=np.int32)
        
    for ii in range(0,len(FullRedshiftList)):                  
        G0=G_MR[ np.rint(G_MR['Redshift']*100.) == FullRedshiftList[ii]*100. ]             
        SnapshotList[ii]=G0['SnapNum'][0]
        
        if(MRII==1):
            G0=G_MRII[ np.rint(G_MRII['Redshift']*100.) == FullRedshiftList[ii]*100. ]             
            SnapshotList_MRII[ii]=G0['SnapNum'][0]
#endif      

print('reading done\n')
#print (np.log10(G_MR['DiskMass'][0:99]*1.e10))
#print (np.log10(G_MR['BulgeMass'][0:99]*1.e10))
#print (np.log10(G_MR['StellarMass'][0:99]*1.e10))
#print (np.log10(G_MR['MetalsStellarMass'][0:99]*1.e10))
#print (G_MR[0:5])
#help(G_MR)


#DEFINE SOME VARIABLE
if(opt_rings==1):
    RNUM=12
    RingRadius=np.zeros(RNUM,dtype=np.float32)
    for ii in range(0,RNUM):
        RingRadius[ii]= 0.44*pow(1.5,ii+1)/Hubble_h;       
        #print(RingRadius[ii])

plt.rcParams.update({'font.size': 18, 'xtick.labelsize': 18, 'ytick.labelsize': 18, 'axes.linewidth': 2, 
                     'xtick.major.size': 6, 'xtick.major.width': 1.5, 
                     'ytick.major.size': 6, 'ytick.major.width': 1.5, 
                     'xtick.minor.size': 3, 'xtick.minor.width': 1.,                   
                     'ytick.minor.size': 3, 'ytick.minor.width': 1.})




# ## Plots

# In[ ]:

import plots
reload (plots)

if(opt_rings==1):
    import plots_rings
    reload (plots_rings)   


# In[14]:




# # Plots for snapshot output

# In[152]:


with PdfPages('./fig/plots.pdf') as pdf:  
    import procedures
    reload (procedures)
    from procedures import *
    import plots_input
    reload (plots_input)
    from plots_input import *
    
    import plots
    reload (plots)
    if(opt_rings==1):
        import plots_rings
        reload (plots_rings) 
    
    #G0_MR=G_MR[(G_MR['SnapNum']==319) & (G_MR['StellarMass']>0.) & (G_MR['Type']==0)]               
    #G0_MR=G_MR[(G_MR['SnapNum']==319) & (G_MR['StellarMass']>0.) & (G_MR['Mvir']>0.) & (G_MR['Type']==1)]
    #G0_MR=G_MR[(G_MR['StellarMass']>0.) & (G_MR['Mvir']>0.) & (G_MR['Type']==1)]   
    #print(len(G0_MR),len(G_MR))
       
    if(MRII==0):
        G_MRII=G_MR  
        
    opt_test_plot=0  
        
    if opt_test_plot==1:
        fig = plt.figure(figsize=(10,10))
        subplot=plt.subplot()
        subplot.set_xlim([6.0,12.]), subplot.set_ylim([0.0, 16.]) 
        subplot.set_xlabel('x', fontsize=16), subplot.set_ylabel('y', fontsize=16)
               
        G0_MR=G_MR[(G_MR['SnapNum']==255) & (G_MR['StellarMass']>0.) & (G_MR['Type']==0)]               
        StellarMass=np.log10(G0_MR['StellarMass']*1.e10*Hubble_h)
        HaloMass=np.log10(G0_MR['Mvir']*1.e10*Hubble_h)   
        BHMass=np.log10(G0_MR['BlackHoleMass']*1.e10*Hubble_h)
        subplot.scatter(StellarMass,HaloMass,s=5, color='black')
        subplot.scatter(StellarMass,BHMass,s=5, color='blue')
                          
    if opt_test_resolution==1:
        print('Doing test resolution')
        from plots import test_resolution
        ThisRedshiftList=[0.0]        
        test_resolution(G_MR, Volume_MR, G_MRII, Volume_MRII, ThisRedshiftList, pdf)
           
    if opt_stellar_mass_vs_halo_mass==1:
        print('Doing SMHM')
        from plots import stellar_mass_vs_halo_mass
        ThisRedshiftList=[0.0]        
        stellar_mass_vs_halo_mass(G_MR, ThisRedshiftList, pdf)
            
    if opt_stellar_mass_function==1:
        print('Doing SMF')
        from plots import stellar_mass_function
        ThisRedshiftList=[0.0,1.0,2.0,3.0]                 
        stellar_mass_function(G_MR, Volume_MR, G_MRII, Volume_MRII, ThisRedshiftList, pdf)
    
    if opt_redfraction_color_cut_cuts==1:
        print('Doing redfraction_color_cut_cuts')
        from plots import redfraction_color_cut_cuts
        ThisRedshiftList=[0.0,0.4,1.0,2.0,3.0]        
        redfraction_color_cut_cuts(G_MR, ThisRedshiftList, pdf)
    
    if opt_redfraction_color_cut==1:
        print('Doing redfraction_color_cut')
        from plots import redfraction_color_cut
        ThisRedshiftList=[0.0,0.4,1.0,2.0,3.0]        
        redfraction_color_cut(G_MR, ThisRedshiftList, pdf)
    
    if opt_metals_vs_stellarmass==1:
        print('Doing metals_vs_stellarmass')
        from plots import metals_vs_stellarmass
        ThisRedshiftList=[0.0,3.]
        metals_vs_stellarmass(G_MR, ThisRedshiftList, pdf)
        
    if opt_gasmetals_vs_stellarmass==1:
        print('Doing gasmetals_vs_stellarmass')
        from plots import gasmetals_vs_stellarmass
        ThisRedshiftList=[0.0,3.]
        gasmetals_vs_stellarmass(G_MR, ThisRedshiftList, pdf)    
        
    if opt_BHBM==1:
        print('Doing BHBM')
        from plots import BHBM
        ThisRedshiftList=[0.0]        
        BHBM(G_MR, ThisRedshiftList, pdf)    
        
    if opt_SFRF==1:
        print('Doing SFRF')
        from plots import SFRF
        ThisRedshiftList=[0.0]        
        SFRF(G_MR, Volume_MR, G_MRII, Volume_MRII, ThisRedshiftList, pdf)
        
    if opt_gas_fraction==1:
        print('Doing gas_fraction')
        from plots import gas_fraction
        ThisRedshiftList=[0.0]        
        gas_fraction(G_MR, ThisRedshiftList, pdf)
        
    if opt_HI_fraction==1:
        print('Doing HI_fraction')
        from plots import HI_fraction
        ThisRedshiftList=[0.0]        
        HI_fraction(G_MR, ThisRedshiftList, pdf)
        
    if opt_HI_MF==1:
        print('Doing HI_MF')
        from plots import HI_MF
        ThisRedshiftList=[0.0]        
        HI_MF(G_MR, Volume_MR, G_MRII, Volume_MRII, ThisRedshiftList, pdf)
        
    if opt_sfr_vs_stellar_mass==1:
        print('Doing sfr_vs_stellar_mass')
        from plots import sfr_vs_stellar_mass
        ThisRedshiftList=[0.0,0.4,1.0,2.0,3.0]        
        sfr_vs_stellar_mass(G_MR, ThisRedshiftList, pdf)
        
    if opt_ur_vs_r==1:
        print('Doing ur_vs_r')
        from plots import ur_vs_r
        ThisRedshiftList=[0.0]        
        ur_vs_r(G_MR, ThisRedshiftList, pdf)
        
    if opt_UVJ_colour==1:
        print('Doing UVJ_colour')
        from plots import UVJ_colour
        ThisRedshiftList=[0.4,1.0,2.0,3.0]        
        UVJ_colour(G_MR, ThisRedshiftList, pdf)
        
    if opt_UVJ_grid==1:
        print('Doing UVJ_grid')
        from plots import UVJ_grid
        ThisRedshiftList=[1.0,1.5,2.0]        
        UVJ_grid(G_MR, ThisRedshiftList, pdf)
        
    if opt_morphology_vs_stellarmass==1:
        print('Doing morphology_vs_stellarmass')
        from plots import morphology_vs_stellarmass
        ThisRedshiftList=[0.0]        
        morphology_vs_stellarmass(G_MR, G_MRII, ThisRedshiftList, pdf)       
            
    if opt_sizes_vs_stellarmass==1:
        print('Doing sizes_vs_stellarmass')
        from plots import sizes_vs_stellarmass
        ThisRedshiftList=[0.0]        
        sizes_vs_stellarmass(G_MR, ThisRedshiftList, pdf)        
            
            
           
    #ADDITIONAL PLOTS   
    if opt_BHBM_by_sfr==1:
        print('Doing BHBM_by_sfr')
        from plots import BHBM_by_sfr
        ThisRedshiftList=[0.0,1.0,2.0]        
        BHBM_by_sfr(G_MR, ThisRedshiftList, pdf)
    
    if opt_AGN_quenching==1:
        print('Doing AGN_quenching')
        from plots import AGN_quenching
        ThisRedshiftList=[0.0,1.0,2.0]        
        AGN_quenching(G_MR, ThisRedshiftList, pdf)
    
    if opt_bluck_red_fractions==1:
        print('Doing bluck_red_fractions')
        from plots import bluck_red_fractions
        ThisRedshiftList=[0.0]
        bluck_red_fractions(G_MR, ThisRedshiftList, pdf)             
        
    if opt_sat_fraction==1:
        print('Doing Sat Fraction')
        from plots import sat_fraction
        ThisRedshiftList=[0.0,0.4,1.0,2.0,3.0]        
        sat_fraction(G_MR, ThisRedshiftList, pdf)
                
    if opt_HotGas_fraction==1:
        print('Doing HotGas_fraction')
        from plots import HotGas_fraction
        ThisRedshiftList=[0.0]        
        HotGas_fraction(G_MR, ThisRedshiftList, pdf) 
        
    if opt_BHmass_in_radio==1:
        print('Doing BHmass_in_radio')
        from plots import BHmass_in_radio
        ThisRedshiftList=[0.0]        
        BHmass_in_radio(G_MR, ThisRedshiftList, pdf) 
        
    if opt_fabian_fb==1:
        print('Doing fabian_fb')
        from plots import fabian_fb
        ThisRedshiftList=[0.0,0.5,1.0,2.0,3.0,5.0]        
        fabian_fb(G_MR, Volume_MR, ThisRedshiftList, pdf)  
        
    if opt_misc_plots==1:
        print('Doing test_plots')
        from plots import test_plots     
        test_plots(G_MR, SnapshotList, pdf)
           
            
         
       
            
    #PLOTS FOR H2_AND_RINGS
    if(opt_rings==1): 
        print('')
        print('')
        print('Doing plots for H2_and_Rings')
        print('')
        if opt_gasfractions_vs_stellarmass==1:
            print('Doing gasfractions_vs_stellarmass')
            from plots_rings import gasfractions_vs_stellarmass
            ThisRedshiftList=[0.0]        
            gasfractions_vs_stellarmass(G_MR, ThisRedshiftList, pdf)
            
        if opt_H2fraction_vs_stellarmass==1:
            print('Doing H2fraction_vs_stellarmass')
            from plots_rings import H2fraction_vs_stellarmass
            ThisRedshiftList=[0.0]        
            H2fraction_vs_stellarmass(G_MR, ThisRedshiftList, pdf)
            
        if opt_milkyway_sfr_and_gas_profiles==1:
            print('Doing milkyway_sfr_and_gas_profiles')
            from plots_rings import milkyway_sfr_and_gas_profiles
            ThisRedshiftList=[0.0]        
            milkyway_sfr_and_gas_profiles(G_MR, RingRadius, RNUM, ThisRedshiftList, pdf)
            
        if opt_evo_milkyway_gas_profile==1:
            print('Doing evo_milkyway_gas_profile')
            from plots_rings import evo_milkyway_gas_profile
            ThisRedshiftList=[0.0]        
            evo_milkyway_gas_profile(G_MR, RingRadius, RNUM, ThisRedshiftList, pdf)    
        
        if opt_test_H2_prescriptions==1:
            print('Doing test_H2_prescriptions')
            from plots_rings import test_H2_prescriptions
            ThisRedshiftList=[0.0]        
            test_H2_prescriptions(G_MR, RingRadius, RNUM, ThisRedshiftList, pdf)
        
        if opt_gas_metallicity_gradients==1:
            print('Doing gas_metallicity_gradients')
            from plots_rings import gas_metallicity_gradients
            ThisRedshiftList=[0.1]
            gas_metallicity_gradients(G_MR, RingRadius, RNUM, ThisRedshiftList, pdf) 
        
        if opt_SFR_gradients==1:
            print('Doing SFR_gradients')
            from plots_rings import SFR_gradients
            ThisRedshiftList=[0.1]
            SFR_gradients(G_MR, RingRadius, RNUM, ThisRedshiftList, pdf) 
            
        if opt_test_rings==1:
            print('Doing test_rings')
            from plots_rings import test_rings
            ThisRedshiftList=[0.0]
            test_rings(G_MR, RingRadius, RNUM, ThisRedshiftList, pdf) 
            
    print('')
    print('All plots done')
           
#end with PdfPages('./fig/plots.pdf') as pdf: 


# # Plots for tree output

# In[205]:

with PdfPages('./fig/plots.pdf') as pdf:  
    import procedures
    reload (procedures)
    from procedures import *
    import plots_input
    reload (plots_input)
    from plots_input import *
    import plots
    reload (plots)
        
    if CatalogType=='tree':    
        if opt_simple_tree_map==1:
            print('Doing simple tree map')
            from plots import simple_tree_map                   
            simple_tree_map(G_MR, pdf)          
    
    if CatalogType=='tree':    
        if opt_full_tree_map==1:
            print('Doing ful tree map')
            from plots import full_tree_map                   
            full_tree_map(G_MR, pdf, object_type='haloes')       
    
    
    print('')
    print('All plots done')
        
#end with PdfPages('./fig/plots.pdf') as pdf: 


# ## Animations

# In[ ]:

import animations
reload (animations)
    
if opt_anime_mass_gr==1:   
    from animations import anime_mass_gr
    anime_mass_gr(G_MR)


# In[ ]:



