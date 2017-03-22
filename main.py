
# coding: utf-8

# # Python Plots for LGalaxies

# ## Import Libraries and Read Catalogs

# <p>Use functions read_snap or read_tree to read catalogs. These are both defined in procedures.py. In case of read_snap, SnapshotList will be returned containing the list of snapshots read (usefull to later select galaxies in a given redshift).<p>

# In[ ]:

import numpy as np
get_ipython().magic('matplotlib inline')

import pandas as pd

get_ipython().magic('pylab inline')
#import seaborn as sns

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

FirstFile = 0
LastFile =  511


Volume_MR = (BoxSize_MR**3.0) * (LastFile - FirstFile + 1) / MaxTreeFiles / Hubble_h**3
Volume_MRII = (BoxSize_MRII**3.0) * (LastFile - FirstFile + 1) / MaxTreeFiles  / Hubble_h**3

print('Reading started')

if CatalogType=='snap':       
    #from LGalaxies_Henriques2015a_struct import LGalaxiesStruct
    #from LGalaxies_Henriques2015a_struct import PropertiesToRead   
    from LGalaxies_HWT16_struct import LGalaxiesStruct
    from LGalaxies_HWT16_struct import PropertiesToRead   
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
    
    #from LGalaxies_Henriques2016a_struct import LGalaxiesStruct
    #from LGalaxies_Henriques2016a_struct import PropertiesToRead
    #from LGalaxies_Henriques2016a_Rings_struct import LGalaxiesStruct
    #from LGalaxies_Henriques2016a_Rings_struct import PropertiesToRead
    #from LGalaxies_Henriques2016a_Elements_struct import LGalaxiesStruct
    #from LGalaxies_Henriques2016a_Elements_struct import PropertiesToRead    
    #from LGalaxies_Henriques2016a_Elements_Rings_struct import LGalaxiesStruct
    #from LGalaxies_Henriques2016a_Elements_Rings_struct import PropertiesToRead
    
    print('\n\nDoing MR')
    (G_MR, SnapshotList_MR) = read_snap(DirName_MR,FirstFile,LastFile,
                     PropertiesToRead,LGalaxiesStruct,RedshiftsToRead,FullRedshiftList)
    
    if(MRII==1):
        print('\n\nDoing MRII')
        (G_MRII, SnapshotList_MRII) = read_snap(DirName_MRII,FirstFile,LastFile,
                         PropertiesToRead,LGalaxiesStruct,RedshiftsToRead,FullRedshiftList)
if CatalogType=='tree':    
    #from LGalaxies_tree_Henriques2015a_struct import LGalaxiesStruct
    #from LGalaxies_tree_Henriques2015a_struct import PropertiesToRead
    from LGalaxies_tree_HWT16_struct import LGalaxiesStruct
    from LGalaxies_tree_HWT16_struct import PropertiesToRead
    #from LGalaxies_tree_Henriques2015a_Rings_struct import LGalaxiesStruct
    #from LGalaxies_tree_Henriques2015a_Rings_struct import PropertiesToRead_tree
    #from LGalaxies_tree_Henriques2015a_Elements_Rings_struct import LGalaxiesStruct
    #from LGalaxies_tree_Henriques2015a_Elements_Rings_struct import PropertiesToRead_tree
    #from LGalaxies_tree_Henriques2015a_Caterpillar_struct import LGalaxiesStruct
    #from LGalaxies_tree_Henriques2015a_Caterpillar_struct import PropertiesToRead_tree   
    #from LGalaxies_tree_Henriques2016a_struct import LGalaxiesStruct
    #from LGalaxies_tree_Henriques2016a_struct import PropertiesToRead_tree    
    #from LGalaxies_tree_Henriques2016a_Elements_Rings_struct import LGalaxiesStruct
    #from LGalaxies_tree_Henriques2016a_Elements_Rings_struct import PropertiesToRead
    
    (G_MR) = read_tree(DirName_MR,FirstFile,LastFile,PropertiesToRead,LGalaxiesStruct)
    SnapshotList = np.zeros(len(FullRedshiftList),dtype=np.int32)
    
    if(MRII==1):
        (G_MRII) = read_tree(DirName_MRII,FirstFile,LastFile,PropertiesToRead,LGalaxiesStruct)    
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
#print (np.log10(G_MR['BulgeMass'][0:99]*1.e10_tree    ))
#print (np.log10(G_MR['StellarMass'][0:99]*1.e10))
#print (np.log10(G_MR['MetalsStellarMass'][0:99]*1.e10))
#print (G_MR[0:5])
#help(G_MR)


#DEFINE SOME VARIABLE
#if(opt_rings==1):
RNUM=12
RingRadius=np.zeros(RNUM,dtype=np.float32)
#in Kpc
for ii in range(0,RNUM):    
    RingRadius[ii]= 0.44*pow(1.5,ii+1)/Hubble_h;   
    #RingRadius[ii]= 0.3*pow(1.2,ii+1)/Hubble_h;   
    #print(RingRadius[ii])


## Plots


# In[ ]:

import plots 
reload (plots) 
import plots_ringsPropertiesToRead['StellarDiskRadius'] = True
PropertiesToRead['StellarHalfMassRadius'] = True
reload (plots)   


# In[14]:




# # Plots for snapshot output

# In[6]:

plt.rcParams.update({'font.size': 14, 'xtick.labelsize': 14, 'ytick.labelsize': 14, 'axes.linewidth': 2, 
                     'xtick.major.size': 6, 'xtick.major.width': 1.5, 
                     'ytick.major.size': 6, 'ytick.major.width': 1.5, 
                     'xtick.minor.size': 3, 'xtick.minor.width': 1.,                   
                     'ytick.minor.size': 3, 'ytick.minor.width': 1.})  
#plt.rcParams.update({'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 
#                             'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0})

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
        
    if opt_stellar_mass_vs_halo_mass_fractional==1:
        print('Doing SMHM_fractional')
        from plots import stellar_mass_vs_halo_mass_fractional
        ThisRedshiftList=[0.0]        
        stellar_mass_vs_halo_mass_fractional(G_MR, G_MRII, ThisRedshiftList, pdf)
            
    if opt_stellar_mass_function==1:
        print('Doing SMF')
        from plots import stellar_mass_function
        ThisRedshiftList=[0.0,1.0,2.0,3.0]                 
        stellar_mass_function(G_MR, Volume_MR, G_MRII, Volume_MRII, ThisRedshiftList, pdf)
    
    if opt_stellar_mass_function_z0_overplot==1:
        print('Doing SMF_z0_overplot')
        from plots import stellar_mass_function_z0_overplot 
        ThisRedshiftList=[0.0] 
        stellar_mass_function_z0_overplot(G_MR, Volume_MR, G_MRII, Volume_MRII, ThisRedshiftList, pdf)
        
    if opt_stellar_mass_function_allz_overplot==1:
        print('Doing SMF allz_overplot')
        from plots import stellar_mass_function_allz_overplot 
        ThisRedshiftList=[0.4,1.0,2.,3.,4.0,5.0] 
        stellar_mass_function_allz_overplot(G_MR, Volume_MR, G_MRII, Volume_MRII, ThisRedshiftList, pdf)    
        
    if opt_stellar_mass_function_feedback_overplot==1:
        print('Doing SMF_feedback_overplot')
        from plots import stellar_mass_function_feedback_overplot 
        ThisRedshiftList=[0.0,1.0,2.0,3.0] 
        stellar_mass_function_feedback_overplot(G_MR, Volume_MR, G_MRII, Volume_MRII, ThisRedshiftList, pdf)
    
    if opt_redfraction_color_cut_cuts==1:
        print('Doing redfraction_color_cut_cuts')
        from plots import redfraction_color_cut_cuts
        ThisRedshiftList=[0.0,0.4,1.0,2.0,3.0]        
        redfraction_color_cut_cuts(G_MR,G_MRII,ThisRedshiftList,pdf)
    
    if opt_redfraction_color_cut==1:
        print('Doing redfraction_color_cut')
        from plots import redfraction_color_cut
        ThisRedshiftList=[0.0,0.4,1.0,2.0,3.0]        
        redfraction_color_cut(G_MR,G_MRII,ThisRedshiftList,pdf)
    
    if opt_metals_vs_stellarmass==1:
        print('Doing metals_vs_stellarmass')
        from plots import metals_vs_stellarmass
        ThisRedshiftList=[0.0,3.]
        metals_vs_stellarmass(G_MR, ThisRedshiftList, pdf)
        
    if opt_gasmetals_vs_stellarmass==1:
        print('Doing gasmetals_vs_stellarmass')
        from plots import gasmetals_vs_stellarmass
        ThisRedshiftList=[0.0,3.]
        gasmetals_vs_stellarmass(G_MR, G_MRII, RingRadius, RNUM, ThisRedshiftList, pdf)    
        
    if opt_BHBM==1:
        print('Doing BHBM')
        from plots import BHBM
        ThisRedshiftList=[0.0]        
        BHBM(G_MR, ThisRedshiftList, pdf)    
        
    if opt_SFRF==1:
        print('Doing SFRF')
        from plots import SFRF
        ThisRedshiftList=[0.0, 1.0, 2.0, 3.0]        
        SFRF(G_MR, Volume_MR, G_MRII, Volume_MRII, ThisRedshiftList, pdf)
        
    if opt_gas_fraction==1:
        print('Doing gas_fraction')
        from plots import gas_fraction
        ThisRedshiftList=[0.0]        
        gas_fraction(G_MR, ThisRedshiftList, pdf)
        
    if opt_HI_over_Lr_vs_HI==1:
        print('Doing HI_over_Lr_vs_HI')
        from plots import HI_over_Lr_vs_HI
        ThisRedshiftList=[0.0]        
        HI_over_Lr_vs_HI(G_MR, ThisRedshiftList, pdf)
        
    if opt_HI_over_Lr_vs_HI_bins==1:
        print('Doing HI_over_Lr_vs_HI_bins')
        from plots import HI_over_Lr_vs_HI_bins
        ThisRedshiftList=[0.0]        
        HI_over_Lr_vs_HI_bins(G_MR, G_MRII, ThisRedshiftList, pdf)
        
    if opt_HI_MF==1:
        print('Doing HI_MF')
        from plots import HI_MF
        ThisRedshiftList=[0.0]        
        HI_MF(G_MR, Volume_MR, G_MRII, Volume_MRII, ThisRedshiftList, pdf)
    
    if opt_coldgas_vs_stellarmass==1:
        print('Doing coldgas_vs_stellarmass')
        from plots import coldgas_vs_stellarmass
        ThisRedshiftList=[0.0]        
        coldgas_vs_stellarmass(G_MR, ThisRedshiftList, pdf)
    
    if opt_main_sequence==1:
        print('Doing main_sequence')
        from plots import main_sequence
        ThisRedshiftList=[0.0,0.4,1.0,2.0,3.0]        
        main_sequence(G_MR, ThisRedshiftList, pdf)
        
    if opt_ur_vs_r==1:
        print('Doing ur_vs_r')
        from plots import ur_vs_r
        ThisRedshiftList=[0.0]        
        ur_vs_r(G_MRII, ThisRedshiftList, pdf)
        
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
    if opt_ssfr_hist==1:
        print('Doing ssfr_hist')
        from plots import ssfr_hist
        ThisRedshiftList=[0.0]        
        ssfr_hist(G_MR, G_MRII, ThisRedshiftList, pdf)
        
    if opt_SFH==1:
        print('Doing SFH')
        from plots import SFH
        ThisRedshiftList=[0.4]        
        SFH(G_MR, ThisRedshiftList, pdf)
        
    if opt_cooling_heating==1:
        print('Doing cooling_heating')
        from plots import cooling_heating
        ThisRedshiftList=[0.0,1.0,2.0,3.0]        
        #BHBM_by_sfr(G_MR,G_MRII, ThisRedshiftList, pdf)
        cooling_heating(G_MR, ThisRedshiftList, pdf)
    
    if opt_BHBM_by_sfr==1:
        print('Doing BHBM_by_sfr')
        from plots import BHBM_by_sfr
        ThisRedshiftList=[0.0,1.0,2.0]        
        #BHBM_by_sfr(G_MR,G_MRII, ThisRedshiftList, pdf)
        BHBM_by_sfr(G_MR, ThisRedshiftList, pdf)
    
    if opt_AGN_quenching==1:
        print('Doing AGN_quenching')
        from plots import AGN_quenching
        ThisRedshiftList=[0.0,1.0,2.0]        
        AGN_quenching(G_MR, ThisRedshiftList, pdf)
    
    if opt_growth_channels==1:
        print('Doing growth_channels')
        from plots import growth_channels
        ThisRedshiftList=[0.0]        
        growth_channels(G_MR, G_MRII, ThisRedshiftList, pdf)       
        
    if opt_bluck_red_fractions==1:
        print('Doing bluck_red_fractions')
        from plots import bluck_red_fractions
        ThisRedshiftList=[0.0]
        bluck_red_fractions(G_MR, ThisRedshiftList, pdf)       
        
    if opt_satellite_quench_times==1:
        print('Doing satellite_quench_times')
        from plots import satellite_quench_times
        ThisRedshiftList=[0.0]
        satellite_quench_times(G_MR, ThisRedshiftList, pdf)             
        
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
    
    if opt_sfr_massive_galaxies==1:
        print('Doing sfr_massive_galaxies')
        from plots import sfr_massive_galaxies      
        sfr_massive_galaxies(G_MR, pdf) 
        
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
            
        if opt_evo_milkyway_stellar_profiles==1:
            print('Doing evo_milkyway_stellar_profiles')
            from plots_rings import evo_milkyway_stellar_profiles
            ThisRedshiftList=[0.0]        
            evo_milkyway_stellar_profiles(G_MR, RingRadius, RNUM, ThisRedshiftList, pdf)       
        
        if opt_test_H2_prescriptions==1:
            print('Doing test_H2_prescriptions')
            from plots_rings import test_H2_prescriptions
            ThisRedshiftList=[0.0]        
            test_H2_prescriptions(G_MR, RingRadius, RNUM, ThisRedshiftList, pdf)
      
        if opt_milkyway_gradients==1:
            print('Doing milkyway_gradients')
            from plots_rings import milkyway_gradients
            ThisRedshiftList=[0.1]
            milkyway_gradients(G_MR, RingRadius, RNUM, ThisRedshiftList, pdf)
    
        if opt_gas_metallicity_gradients_mass_bins==1:
            print('Doing gas_metallicity_gradients_mass_bins')
            from plots_rings import gas_metallicity_gradients_mass_bins
            ThisRedshiftList=[0.1]
            gas_metallicity_gradients_mass_bins(G_MR, RingRadius, RNUM, ThisRedshiftList, pdf)
            
        if opt_stellar_metallicity_gradients_mass_bins==1:
            print('Doing stellar_metallicity_gradients_mass_bins')
            from plots_rings import stellar_metallicity_gradients_mass_bins
            ThisRedshiftList=[0.0]
            stellar_metallicity_gradients_mass_bins(G_MR, RingRadius, RNUM, ThisRedshiftList, pdf)
            
        if opt_CALIFA_gradients_morph_types==1:
            print('Doing CALIFA_gradients_morph_types')
            from plots_rings import CALIFA_gradients_morph_types
            ThisRedshiftList=[0.0]
            CALIFA_gradients_morph_types(G_MR, RingRadius, RNUM, ThisRedshiftList, pdf)
            
        if opt_CALIFA_gradients_mass_bins==1:
            print('Doing CALIFA_gradients_mass_bins')
            from plots_rings import CALIFA_gradients_mass_bins
            ThisRedshiftList=[0.0]
            CALIFA_gradients_mass_bins(G_MR, RingRadius, RNUM, ThisRedshiftList, pdf)    
            
        if opt_MANGA_gradients_late_types==1:
            print('Doing MANGA_gradients_late_types')
            from plots_rings import MANGA_gradients_late_types
            ThisRedshiftList=[0.0]
            MANGA_gradients_late_types(G_MR, RingRadius, RNUM, ThisRedshiftList, pdf)
                                    
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
            simple_tree_lmap(G_MR, pdf)          
    
    if CatalogType=='tree':    
        if opt_full_tree_map==1:
            print('Doing ful tree map')
            from plots import full_tree_map                   
            full_tree_map(G_MR, pdf, object_type='haloes')       
    
    
    print('')
    print('All plots done')
        
#end with PdfPages('./fig/plots.pdf') as pdf: 


# # Compare 2 files

# In[103]:

#RedshiftsToRead = [False,False,False,False,False,False,True,False]
#RedshiftsToRead = [False,False,False,False,True,False,False,False]
#RedshiftsToRead = [False,False,False,True,False,False,False,False]
RedshiftsToRead = [True,False,False,False,False,False,False,False]

DirName_1 = '/net/bootes/scratch-ssd/SAM/test0/MR/'
DirName_2 = '/net/bootes/scratch-ssd/SAM/test1/MR/'

FirstFile = 40
LastFile  = 40

(G1, SnapshotList_MR) = read_snap(DirName_1,FirstFile,LastFile,
                                    PropertiesToRead,LGalaxiesStruct,RedshiftsToRead,FullRedshiftList)
(G2, SnapshotList_MR) = read_snap(DirName_2,FirstFile,LastFile,
                                    PropertiesToRead,LGalaxiesStruct,RedshiftsToRead,FullRedshiftList)
   

N=200
sel=np.random.rand(len(G1))<float(N/len(G1))
G1=G1[sel]
G2=G2[sel]

if(len(G1)>len(G2)):
    G1=G1[0:len(G2)]
else:    
    G2=G2[0:len(G1)]

print(len(G1), len(G2))
    
with PdfPages('./fig/plots.pdf') as pdf:  
   
    fig = plt.figure(figsize=(15,15))
    grid = gridspec.GridSpec(6, 2)
    
    xlim=[0.001,100.]
    ylim=[0.001,100.]    
    subplot=plt.subplot(grid[0])
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)   
    subplot.set_xscale("log")
    subplot.set_yscale("log")
    plt.scatter(G1['StellarMass'],G2['StellarMass'],s=2,color='black')
    
    xlim=[0.0001,0.1]
    ylim=[0.0001,0.1]    
    subplot=plt.subplot(grid[1])
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)   
    subplot.set_xscale("log")
    subplot.set_yscale("log")
    plt.scatter(G1['StellarHalfMassRadius'],G2['StellarHalfMassRadius'],s=2,color='black')
   
    '''xlim=[0.1,1e12]
    ylim=[0.1,1e12]      
    subplot=plt.subplot(grid[1])
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
    subplot.set_xscale("log")
    subplot.set_yscale("log")  
    plt.scatter(G1['DiskMass_elements'][:,0],G2['DiskMass_elements'][:,0],s=2,color='black')
    
    subplot=plt.subplot(grid[2])
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
    subplot.set_xscale("log")
    subplot.set_yscale("log")  
    plt.scatter(G1['DiskMass_elements'][:,5],G2['DiskMass_elements'][:,5],s=2,color='black')'''
    
    '''xlim=[0.1,1e15]
    ylim=[0.1,1e15]
    subplot=plt.subplot(grid[0])
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
    subplot.set_xscale("log")
    subplot.set_yscale("log")  
    plt.scatter(G1['HotGas_elements'][:,0],G2['HotGas_elements'][:,0],s=2,color='black')
    
    subplot=plt.subplot(grid[1])
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
    subplot.set_xscale("log")
    subplot.set_yscale("log")  
    plt.scatter(G1['HotGas_elements'][:,5],G2['HotGas_elements'][:,5],s=2,color='black')'''
    
    '''xlim=[0.1,1e12]
    ylim=[0.1,1e12]      
    subplot=plt.subplot(grid[0])
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
    subplot.set_xscale("log")
    subplot.set_yscale("log")  
    plt.scatter(G1['sfh_ElementsDiskMass'][:,0,0],G2['sfh_ElementsDiskMass'][:,0,0],s=2,color='black')
    
    subplot=plt.subplot(grid[1])
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
    subplot.set_xscale("log")
    subplot.set_yscale("log")  
    plt.scatter(G1['sfh_ElementsDiskMass'][:,5,5],G2['sfh_ElementsDiskMass'][:,5,5],s=2,color='black')'''
    
    xlim=[1e5,1e15]
    ylim=[1e5,1e15]
    subplot=plt.subplot(grid[2])
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
    subplot.set_xscale("log")
    subplot.set_yscale("log")  
    x=np.arange(xlim[0],xlim[1],(xlim[1]-xlim[0])/100.)
    y=x
    subplot.plot(x,y)
    subplot.scatter(G1['HotGas_elements'][:,0],G2['HotGas_elements'][:,0],s=2,color='black')
    
    xlim=[1e2,1e15]
    ylim=[1e2,1e15]
    subplot=plt.subplot(grid[3])
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
    subplot.set_xscale("log")
    subplot.set_yscale("log")  
    x=np.arange(xlim[0],xlim[1],(xlim[1]-xlim[0])/100.)
    y=x
    subplot.plot(x,y)
    subplot.scatter(G1['HotGas_elements'][:,2],G2['HotGas_elements'][:,2],s=2,color='black')
        
    print("HOT GAS")
    sel=(G1['HotGas_elements'][:,0]!=G2['HotGas_elements'][:,0])
    print(G1['HaloIndex'][sel])
    print(G2['HaloIndex'][sel])
    print(G1['HotGas_elements'][sel,0]/G2['HotGas_elements'][sel,0])    
    print("")
    
    
    
    xlim=[0.1,1e12]
    ylim=[0.1,1e12]
    subplot=plt.subplot(grid[4])
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
    subplot.set_xscale("log")
    subplot.set_yscale("log")  
    x=np.arange(xlim[0],xlim[1],(xlim[1]-xlim[0])/100.)
    y=x
    subplot.plot(x,y)
    subplot.scatter(G1['ColdGas_elements'][:,0],G2['ColdGas_elements'][:,0],s=2,color='black')

    
    subplot=plt.subplot(grid[5])
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
    subplot.set_xscale("log")
    subplot.set_yscale("log")  
    subplot.plot(x,y)
    plt.scatter(G1['ColdGas_elements'][:,5],G2['ColdGas_elements'][:,5],s=2,color='black')
  
    subplot=plt.subplot(grid[6])
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
    subplot.set_xscale("log")
    subplot.set_yscale("log")  
    subplot.plot(x,y)
    plt.scatter(G1['DiskMass_elements'][:,0],G2['DiskMass_elements'][:,0],s=2,color='black')

    subplot=plt.subplot(grid[7])
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
    subplot.set_xscale("log")
    subplot.set_yscale("log")
    subplot.plot(x,y)
    plt.scatter(G1['DiskMass_elements'][:,5],G2['DiskMass_elements'][:,5],s=2,color='black')
    
    '''xlim=[0.1,1e12]
    ylim=[0.1,1e12]      
    subplot=plt.subplot(grid[5])
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
    subplot.set_xscale("log")
    subplot.set_yscale("log")  
    plt.scatter(G1['sfh_ElementsDiskMass'][:,0,0],G2['sfh_ElementsDiskMass'][:,0,0],s=2,color='black')
    
    subplot=plt.subplot(grid[6])
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
    subplot.set_xscale("log")
    subplot.set_yscale("log")  
    plt.scatter(G1['sfh_ElementsDiskMass'][:,5,5],G2['sfh_ElementsDiskMass'][:,5,5],s=2,color='black')'''
    
    xlim=[0.00001,1e4]
    ylim=[0.00001,1e4]   
    subplot=plt.subplot(grid[8])
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
    subplot.set_xscale("log")
    subplot.set_yscale("log")  
    x=np.arange(xlim[0],xlim[1],(xlim[1]-xlim[0])/100.)
    y=x
    subplot.plot(x,y)
    plt.scatter(G1['MetalsColdGas'][:,0],G2['MetalsColdGas'][:,0],s=2,color='black')
    
    
    
    xlim=[0.0000001,1e4]
    ylim=[0.0000001,1e4]
    subplot=plt.subplot(grid[9])
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
    subplot.set_xscale("log")
    subplot.set_yscale("log")  
    x=np.arange(xlim[0],xlim[1],(xlim[1]-xlim[0])/100.)
    y=x
    subplot.plot(x,y)
    plt.scatter(G1['MetalsColdGas'][:,1],G2['MetalsColdGas'][:,1],s=2,color='black')
    
    sel=(G1['MetalsColdGas'][:,0]!=G2['MetalsColdGas'][:,0])
    print(G1['HaloIndex'][sel])
    print(G2['HaloIndex'][sel])
    print(G1['ColdGas_elements'][sel,0]/G2['ColdGas_elements'][sel,0])
        
    '''xlim=[0.1,1e12]
    ylim=[0.1,1e12]   
    subplot=plt.subplot(grid[4])
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
    subplot.set_xscale("log")
    subplot.set_yscale("log")  
    x=np.arange(xlim[0],xlim[1],(xlim[1]-xlim[0])/100.)
    y=x
    subplot.plot(x,y)
    plt.scatter(G1['ColdGasRings_elements'][:,0,0],G2['ColdGasRings_elements'][:,0,0],s=2,color='black')
    
    
    
    xlim=[0.0000001,1e12]
    ylim=[0.0000001,1e12]
    subplot=plt.subplot(grid[5])
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
    subplot.set_xscale("log")
    subplot.set_yscale("log")  
    x=np.arange(xlim[0],xlim[1],(xlim[1]-xlim[0])/100.)
    y=x
    subplot.plot(x,y)
    plt.scatter(G1['ColdGasRings_elements'][:,2,5],G2['ColdGasRings_elements'][:,2,5],s=2,color='black')'''
    
    xlim=[0.1,1e12]
    ylim=[0.1,1e12]      
    subplot=plt.subplot(grid[10])
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
    subplot.set_xscale("log")
    subplot.set_yscale("log")  
    x=np.arange(xlim[0],xlim[1],(xlim[1]-xlim[0])/100.)
    y=x
    subplot.plot(x,y)
    plt.scatter(G1['DiskMassRings_elements'][:,0,0],G2['DiskMassRings_elements'][:,0,0],s=2,color='black')
    
    xlim=[0.0000001,1e12]
    ylim=[0.0000001,1e12] 
    subplot=plt.subplot(grid[11])
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
    subplot.set_xscale("log")
    subplot.set_yscale("log")
    x=np.arange(xlim[0],xlim[1],(xlim[1]-xlim[0])/100.)
    y=x
    subplot.plot(x,y)
    plt.scatter(G1['DiskMassRings_elements'][:,10,0],G2['DiskMassRings_elements'][:,10,0],s=2,color='black')
        
   
    
    
    pdf.savefig()
    plt.close()
    
    
   
	
    
    print("All Plots Done!")


# # Create Fits files

# In[27]:

ThisRedshiftList=[0.0]
(sel)=select_current_redshift(G_MR, ThisRedshiftList, 0, FullSnapshotList_MR)
G0_MR=G_MR[sel] 

G=G0_MR

# L: Logical (Boolean)
# B: Unsigned Byte
# I: 16-bit Integer
# J: 32-bit Integer
# K: 64-bit Integer
# E: Single-precision Floating Point
# D: Double-precision Floating Point
# C: Single-precision Complex
# M: Double-precision Complex
# A: Character
col1 = fits.Column(name='Type',                   format='I', array=G['Type'])
col2 = fits.Column(name='SnapNum',                format='I', array=G['SnapNum'])
col3 = fits.Column(name='CentralMvir',            format='F', array=G['CentralMvir'])
col4 = fits.Column(name='DistanceToCentralGalX',  format='F', array=G['DistanceToCentralGal'][:,0])
col5 = fits.Column(name='DistanceToCentralGalY',  format='F', array=G['DistanceToCentralGal'][:,1])
col6 = fits.Column(name='DistanceToCentralGalZ',  format='F', array=G['DistanceToCentralGal'][:,2])
col7 = fits.Column(name='PosX',                   format='F', array=G['Pos'][:,0])
col8 = fits.Column(name='PosY',                   format='F', array=G['Pos'][:,1])
col9 = fits.Column(name='PosZ',                   format='F', array=G['Pos'][:,2])
col10 = fits.Column(name='VelX',                  format='F', array=G['Vel'][:,0])
col11 = fits.Column(name='VelY',                  format='F', array=G['Vel'][:,1])
col12 = fits.Column(name='VelZ',                  format='F', array=G['Vel'][:,2])
col13 = fits.Column(name='Mvir',                  format='F', array=G['Mvir'])
col14 = fits.Column(name='Rvir',                  format='F', array=G['Rvir'])
col15 = fits.Column(name='Vvir',                  format='F', array=G['Vvir'])
col16 = fits.Column(name='flagSplashBack',        format='I', array=G['flagSplashBack'])
col17 = fits.Column(name='TimeSinceSplashBack',   format='F', array=G['TimeSinceSplashBack'])
col18 = fits.Column(name='ColdGas',               format='F', array=G['ColdGas'])
col19 = fits.Column(name='StellarMass',           format='F', array=G['StellarMass'])
col20 = fits.Column(name='BulgeMass',             format='F', array=G['BulgeMass'])
col21 = fits.Column(name='DiskMass',              format='F', array=G['DiskMass'])
col22 = fits.Column(name='BlackHoleMass',         format='F', array=G['BlackHoleMass'])
col23 = fits.Column(name='Sfr',                   format='F', array=G['Sfr'])
col24 = fits.Column(name='BulgeSize',             format='F', array=G['BulgeSize'])
col25 = fits.Column(name='StellarHalfMassRadius', format='F', array=G['StellarHalfMassRadius'])
col26 = fits.Column(name='MagDust_sdssu',         format='F', array=G['MagDust'][:,15])
col27 = fits.Column(name='MagDust_sdssg',         format='F', array=G['MagDust'][:,16])
col28 = fits.Column(name='MagDust_sdssr',         format='F', array=G['MagDust'][:,17])
col29 = fits.Column(name='MagDust_sdssi',         format='F', array=G['MagDust'][:,18])
col30 = fits.Column(name='MagDust_sdssz',         format='F', array=G['MagDust'][:,19])
col31 = fits.Column(name='MassWeightAge',         format='F', array=G['MassWeightAge'])
col32 = fits.Column(name='rBandWeightAge',        format='F', array=G['rBandWeightAge'])


cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,
                     col11,col12,col13,col14,col15,col16,col17,col18,col19,col20,
                     col21,col22,col23,col24,col25,col26,col27,col28,col29,col30,
                     col31,col32])

tbhdu = fits.new_table(cols)


file=DirName_MR+'snap_MR_H15_BackSplash_z0.00.fits'

prihdr = fits.Header()
prihdr['COMMENT'] = "File containing the galaxy catalog for the Henriques et al. 2015 model"
prihdu = fits.PrimaryHDU(header=prihdr)
thdulist = fits.HDUList([prihdu, tbhdu])
#thdulist.writeto(file, clobber=True)

#header_0 = fits.Header()
#hdu = fits.PrimaryHDU(data=data_arr) header=header_0)
#hdulist = fits.HDUList([hdu])
#hdulist.writeto(file, clobber=True)


print("")
print("file written")
print("")


# In[28]:

DirName_MR = '/net/bootes/scratch2/SAM/Henriques2015a/snaps/BackSplashFlag/MR/'
file=DirName_MR+'snap_MR_H15_BackSplash_z0.00.fits'
HDUList=fits.open(file)
print("")
print("file read")
print("")

#HDUList.info()

model = HDUList[1].data 
#print(model['StellarMass'])



sel=((np.log10(model['Stellarmass']*1.e10/0.673)>9.25) & (model['Stellarmass']>0.) & 
     (model['PosX']/0.673<350) & (model['PosY']/0.673<350) & (model['PosZ']/0.673<350))


print(len(model[sel]))


HDUList.close()
print("")
print("all done")
print("")


# ## Animations

# In[ ]:

import animations
reload (animations)
    
if opt_anime_mass_gr==1:   
    from animations import anime_mass_gr
    anime_mass_gr(G_MR)


# In[ ]:



