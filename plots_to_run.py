import sys
import gc
import time
import numpy as np
import psutil

import pandas as pd


#import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import concurrent.futures
from astropy.table import Table
from astropy.io import fits
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
from decimal import *

from scipy.ndimage import zoom 
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from importlib import reload



#PLOT OPTIONS
plot_stellar_mass_vs_halo_mass_vandokkum=0
plot_stellar_mass_vs_halo_mass_fractional=0
plot_stellar_mass_vs_halo_mass_fractional_allz=0
plot_all_masses_vs_halo_mass_fractional_allz=0
plot_eject_hot_masses_vs_halo_mass_fractional_allz=0
plot_all_masses_vs_halo_mass_fractional_z0=0
plot_stellar_mass_vs_halo_mass_fractional_models_overplot=0
plot_cooling_radius=0
plot_stellar_mass_function_z0_overplot=0
plot_stellar_mass_function_feedback_overplot=0
plot_redfraction_color_cut_cuts=0
plot_stellar_mass_function_allz_overplot=0
plot_cumulative_ND=0
plot_halo_mass_function=0
plot_halo_mass_density_above_mass=0

plot_stellar_mass_function=0
plot_redfraction_color_cut=0
plot_redfraction_SFR_cut=0
plot_wetzel_passive_fraction_vs_stellar_mass=0

plot_metals_vs_stellarmass=0
plot_morphology_vs_stellarmass=0
plot_HI_over_Lr_vs_HI_bins=0
plot_sizes_vs_stellarmass=0
plot_sizes_vs_stellarmass_allz=0
plot_gasmetals_vs_stellarmass=0
plot_gasfractions_vs_stellarmass=0
plot_HI_MF=0

plot_BHBM=0
plot_BHMvir=0
plot_SFRF=0
plot_SFRD=0
plot_main_sequence=0
plot_SSFR_mass = 0
plot_ssfr_hist=0

plot_schmidt_kenn=0
plot_test_metal_evo=0



plot_effective_yield=0
plot_metals_mass=0
plot_coldgas_vs_stellarmass=0
plot_H2fraction_vs_stellarmass=0
plot_gas_fraction=0

plot_HI_over_Lr_vs_HI=0
plot_ur_vs_r=0
plot_UVJ_colour=0
plot_UVJ_grid=0
    
#options for H2_AND_RINGS
plot_milkyway_sfr_and_gas_profiles=0
plot_milkyway_gradients=0
plot_gas_metallicity_gradients_mass_bins=0
plot_MANGA_gradients_late_types=0
plot_CALIFA_gradients_morph_types=0
plot_CALIFA_gradients_mass_bins=0
plot_stellar_metallicity_gradients_mass_bins=0
plot_SFR_gradients=0

plot_evo_milkyway_gas_profile=0

plot_evo_milkyway_stellar_profiles=0
plot_test_H2_prescriptions=0

      
#MISC
plot_SFH=0

plot_cooling_heating=1
plot_BHBM_by_sfr=0
plot_AGN_quenching=0
plot_growth_channels=0
plot_bluck_red_fractions=0
plot_satellite_quench_times=0
plot_sat_fraction=0
plot_number_density_massive_gals=0
plot_HotGas_fraction=0
plot_BHmass_in_radio=0
plot_fabian_fb=0 
plot_H2fraction_fits=0
plot_surface_density_vs_stellar_mass=0

#TESTS
plot_test_resolution=0
plot_misc_plots=0
plot_test_rings=0
    
#when tree on
plot_simple_tree_map=0
plot_full_tree_map=0
plot_sfr_massive_galaxies=0
plot_all_masses_evo=0
plot_mass_fractions_evo_all_single_gal=0
plot_mass_fractions_evo_all_multiple_gals=0
plot_properties_evo=0
plot_halo_growth=0
plot_halo_growth_normalized=0
plot_baryon_fraction=0
plot_halo_growth_rate=0
plot_accretion_history=0
plot_cooling_heating_evo=0
plot_cooling_heating_evo_halo_mass=0

#Favignana
#plot_stellar_mass_function=1
#plot_cumulative_ND=1
#plot_stellar_mass_vs_halo_mass_fractional_allz=1
#plot_sizes_vs_stellarmass_allz=1
#plot_redfraction_SFR_cut=1






#HWB17
#plot_stellar_mass_function=0
#plot_stellar_mass_function_z0_overplot=0
#plot_stellar_mass_vs_halo_mass_fractional=0
#plot_stellar_mass_function_feedback_overplot=0
#plot_all_masses_vs_halo_mass_fractional_z0=0
#plot_cooling_heating=0
#plot_growth_channels=1
#Tree
#plot_all_masses_evo=0
#plot_mass_fractions_evo_all_single_gal=0

#HYJ18
#plot_test_resolution=1
#plot_stellar_mass_function=1
#plot_redfraction_color_cut=1
#plot_metals_vs_stellarmass=1
#plot_morphology_vs_stellarmass=1
#plot_HI_over_Lr_vs_HI_bins=1
#plot_sizes_vs_stellarmass=1
#plot_gasmetals_vs_stellarmass=1
#plot_gasfractions_vs_stellarmass=1
#plot_HI_MF=1
#plot_BHBM=1
#plot_SFRF=1
#plot_SFRD=0
#plot_main_sequence=1
#plot_ssfr_hist=1
#plot_milkyway_sfr_and_gas_profiles=1
#plot_milkyway_gradients=1
#plot_gas_metallicity_gradients_mass_bins=0
#plot_MANGA_gradients_late_types=1
#plot_CALIFA_gradients_morph_types=1
#plot_CALIFA_gradients_mass_bins=1


#ANIMATIONS
plot_anime_mass_gr=0

import plots
reload (plots)
from plots import *

import plots_rings
reload (plots_rings) 
from plots_rings import *

#end switches
def run_plots(plot_to_run):   
    
    print('doing %s\n' %plot_to_run)  
    
    output=0
    
    if plot_to_run == 'test_resolution':      
        ThisRedshiftList=[0.0]        
        output = test_resolution(ThisRedshiftList)
    
    if plot_to_run == 'stellar_mass_vs_halo_mass_vandokkum':    
        ThisRedshiftList=[0.0]        
        output = stellar_mass_vs_halo_mass_vandokkum(ThisRedshiftList)
        
    if plot_to_run == 'stellar_mass_vs_halo_mass_fractional':   
        ThisRedshiftList=[0.0]        
        output = stellar_mass_vs_halo_mass_fractional(ThisRedshiftList)
        
    if plot_to_run == 'stellar_mass_vs_halo_mass_fractional_allz':   
        ThisRedshiftList=[0.0, 1.0, 2.0, 3.0]        
        output = stellar_mass_vs_halo_mass_fractional_allz(ThisRedshiftList)    
        
    if plot_to_run == 'all_masses_vs_halo_mass_fractional_allz':    
        ThisRedshiftList=[0.0, 1.0, 2.0, 3.0]        
        output = all_masses_vs_halo_mass_fractional_allz(ThisRedshiftList)
        
    if plot_to_run == 'eject_hot_masses_vs_halo_mass_fractional_allz':    
        ThisRedshiftList=[0.0, 1.0, 2.0, 3.0]        
        output = eject_hot_masses_vs_halo_mass_fractional_allz(ThisRedshiftList)
                
    if plot_to_run == 'all_masses_vs_halo_mass_fractional_z0':    
        ThisRedshiftList=[0.0]        
        output = all_masses_vs_halo_mass_fractional_z0(ThisRedshiftList)
        
    if plot_to_run == 'stellar_mass_vs_halo_mass_fractional_models_overplot':    
        ThisRedshiftList=[0.0]        
        output = stellar_mass_vs_halo_mass_fractional_models_overplot(ThisRedshiftList)
    
    if plot_to_run == 'cumulative_ND':    
        ThisRedshiftList=[0.1, 0.4, 1.0, 2.0, 3.0] 
        output = cumulative_ND(ThisRedshiftList)
    
    if plot_to_run == 'cooling_radius':    
        ThisRedshiftList=[0.0]        
        output = cooling_radius(ThisRedshiftList)        
                    
    if plot_to_run == 'halo_mass_function':            
        ThisRedshiftList=[0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0]                 
        output = halo_mass_function(ThisRedshiftList)        
          
    if plot_to_run == 'halo_mass_density_above_mass':            
        ThisRedshiftList=[7.0,6.0,5.0,4.0,3.0,2.0,1.0,0.86,0.71,0.57,0.51,0.4, 0.26, 0.1, 0.0]                 
        output = halo_mass_density_above_mass(ThisRedshiftList)             
                        
    if plot_to_run == 'stellar_mass_function':            
        ThisRedshiftList=[0.0,1.0,2.0,3.0]                 
        output = stellar_mass_function(ThisRedshiftList)
    
    if plot_to_run == 'stellar_mass_function_z0_overplot':    
        ThisRedshiftList=[0.0] 
        output = stellar_mass_function_z0_overplot(ThisRedshiftList)
        
    if plot_to_run == 'stellar_mass_function_allz_overplot':    
        ThisRedshiftList=[0.4,1.0,2.,3.,4.0,5.0] 
        output = stellar_mass_function_allz_overplot(ThisRedshiftList)    
        
    if plot_to_run == 'stellar_mass_function_feedback_overplot':    
        ThisRedshiftList=[0.0,1.0,2.0,3.0] 
        output = stellar_mass_function_feedback_overplot(ThisRedshiftList)
    
    if plot_to_run == 'redfraction_color_cut_cuts':    
        ThisRedshiftList=[0.0,0.4,1.0,2.0,3.0]        
        output = redfraction_color_cut_cuts(ThisRedshiftList)
    
    if plot_to_run == 'redfraction_color_cut':    
        ThisRedshiftList=[0.0,0.4,1.0,2.0,3.0]        
        output = redfraction_color_cut(ThisRedshiftList)
        
    if plot_to_run == 'redfraction_SFR_cut':    
        ThisRedshiftList=[0.1,1.0,2.0]        
        output = redfraction_SFR_cut(ThisRedshiftList)
        
    if plot_to_run == 'wetzel_passive_fraction_vs_stellar_mass':    
        ThisRedshiftList=[0.1]        
        output = wetzel_passive_fraction_vs_stellar_mass(ThisRedshiftList)
    
    if plot_to_run == 'metals_vs_stellarmass':    
        ThisRedshiftList=[0.0,3.]
        output = metals_vs_stellarmass(ThisRedshiftList)
        
    if plot_to_run == 'gasmetals_vs_stellarmass':    
        ThisRedshiftList=[0.1,2.]
        output = gasmetals_vs_stellarmass(ThisRedshiftList)    
        
    if plot_to_run == 'effective_yield':    
        ThisRedshiftList=[0.1,2.]
        output = effective_yield(ThisRedshiftList)    
        
    if plot_to_run == 'metals_mass':    
        ThisRedshiftList=[0.1,3.]
        output = metals_mass(ThisRedshiftList) 
        
    if plot_to_run == 'test_metal_evo ':   
        ThisRedshiftList=[0.0, 1.0, 2.0, 3., 5.0, 7.0]
        output = test_metal_evo( ThisRedshiftList)     
          
    if plot_to_run == 'BHBM':    
        ThisRedshiftList=[0.0]        
        output = BHBM(ThisRedshiftList) 
        
    if plot_to_run == 'BHMvir':    
        ThisRedshiftList=[0.0]        
        output = BHMvir(ThisRedshiftList)     
        
    if plot_to_run == 'SFRF':    
        ThisRedshiftList=[0.0, 1.0, 2.0, 3.0]        
        output = SFRF(ThisRedshiftList)
        
    if plot_to_run == 'SFRD':   
        ThisRedshiftList=[0.1, 0.4, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]        
        output = SFRD(ThisRedshiftList)
        
    if plot_to_run == 'gas_fraction':   
        ThisRedshiftList=[0.0]        
        output = gas_fraction(ThisRedshiftList)
        
    if plot_to_run == 'HI_over_Lr_vs_HI':    
        ThisRedshiftList=[0.0]        
        output = HI_over_Lr_vs_HI(ThisRedshiftList)
        
    if plot_to_run == 'HI_over_Lr_vs_HI_bins':   
        ThisRedshiftList=[0.0]        
        output = HI_over_Lr_vs_HI_bins(ThisRedshiftList)
        
    if plot_to_run == 'HI_MF':    
        ThisRedshiftList=[0.0]        
        output = HI_MF(ThisRedshiftList)
    
    if plot_to_run == 'coldgas_vs_stellarmass':    
        ThisRedshiftList=[0.0]        
        output = coldgas_vs_stellarmass(ThisRedshiftList)
    
    if plot_to_run == 'main_sequence':   
        ThisRedshiftList=[0.0,0.4,1.0,2.0,3.0]        
        output = main_sequence(ThisRedshiftList)
        
    if plot_to_run == 'SSFR_mass':    
        ThisRedshiftList=[0.0]        
        output = SSFR_mass(ThisRedshiftList)
        
    if plot_to_run == 'ur_vs_r':  
        ThisRedshiftList=[0.0]        
        output = ur_vs_r(ThisRedshiftList)
        
    if plot_to_run == 'UVJ_colour':   
        ThisRedshiftList=[0.4,1.0,2.0,3.0]        
        output = UVJ_colour(ThisRedshiftList)
        
    if plot_to_run == 'UVJ_grid':  
        ThisRedshiftList=[1.0,1.5,2.0]        
        output = UVJ_grid(ThisRedshiftList)
      
    if plot_to_run == 'morphology_vs_stellarmass':  
        ThisRedshiftList=[0.0]        
        output = morphology_vs_stellarmass(ThisRedshiftList)       
            
    if plot_to_run == 'sizes_vs_stellarmass':    
        ThisRedshiftList=[0.0]        
        output = sizes_vs_stellarmass(ThisRedshiftList)
        
    if plot_to_run == 'sizes_vs_stellarmass_allz':    
        ThisRedshiftList=[0.0, 1.0, 2.0, 3.0]        
        output = sizes_vs_stellarmass_allz(ThisRedshiftList)
                 
    #ADDITIONAL PLOTS   
    if plot_to_run == 'ssfr_hist':   
        ThisRedshiftList=[0.0]        
        output = ssfr_hist(ThisRedshiftList)
        
    if plot_to_run == 'SFH':  
        ThisRedshiftList=[0.4]        
        output = SFH(ThisRedshiftList)
        
    if plot_to_run == 'cooling_heating': 
        ThisRedshiftList=[0.0,1.0,2.0,3.0]        
        output = cooling_heating(ThisRedshiftList)
    
    if plot_to_run == 'BHBM_by_sfr':  
        ThisRedshiftList=[0.0,1.0,2.0]  
        output = BHBM_by_sfr(ThisRedshiftList)
    
    if plot_to_run == 'AGN_quenching':   
        ThisRedshiftList=[0.0,1.0,2.0]        
        output = AGN_quenching(ThisRedshiftList)
    
    if plot_to_run == 'growth_channels':   
        ThisRedshiftList=[0.0]        
        output = growth_channels(ThisRedshiftList)       
        
    if plot_to_run == 'bluck_red_fractions':   
        ThisRedshiftList=[0.0]
        output = bluck_red_fractions(ThisRedshiftList)       
        
    if plot_to_run == 'satellite_quench_times':    
        ThisRedshiftList=[0.0]
        output = satellite_quench_times(ThisRedshiftList)             
        
    if plot_to_run == 'sat_fraction':   
        ThisRedshiftList=[0.0,0.4,1.0,2.0,3.0]        
        output = sat_fraction(ThisRedshiftList)
       
    if plot_to_run == 'number_density_massive_gals':    
        ThisRedshiftList=[2.0,3.0,4.0,5.0]        
        output = number_density_massive_gals(ThisRedshiftList)
        
    if plot_to_run == 'HotGas_fraction':  
        ThisRedshiftList=[0.0]        
        output = HotGas_fraction(ThisRedshiftList) 
        
    if plot_to_run == 'BHmass_in_radio':   
        ThisRedshiftList=[0.0]        
        output = BHmass_in_radio(ThisRedshiftList) 
    
    if plot_to_run == 'sfr_massive_galaxies':  
        output = sfr_massive_galaxies()
        
    if plot_to_run == 'fabian_fb':  
        ThisRedshiftList=[0.0,0.5,1.0,2.0,3.0,5.0]        
        output = fabian_fb(ThisRedshiftList)  
    
    if plot_to_run == 'schmidt_kenn':    
        ThisRedshiftList=[2.0]        
        output = schmidt_kenn(ThisRedshiftList)
          
    if plot_to_run == 'misc_plots':   
        output = misc_plots(SnapshotList)
           
    
    if plot_to_run == 'H2fraction_fits': 
        output = H2fraction_fits()  
        
    if plot_to_run == 'surface_density_vs_stellar_mass': 
        ThisRedshiftList=[0.0, 1.0, 2.0]
        output = surface_density_vs_stellar_mass(ThisRedshiftList)  
          
    #PLOTS FOR H2_AND_RINGS
    if plot_to_run == 'gasfractions_vs_stellarmass':       
        ThisRedshiftList=[0.0]        
        output = gasfractions_vs_stellarmass(ThisRedshiftList)
           
    if plot_to_run == 'H2fraction_vs_stellarmass':     
        ThisRedshiftList=[0.0]        
        output = H2fraction_vs_stellarmass(ThisRedshiftList)
            
    if plot_to_run == 'milkyway_sfr_and_gas_profiles':     
        ThisRedshiftList=[0.0]        
        output = milkyway_sfr_and_gas_profiles(ThisRedshiftList)
            
    if plot_to_run == 'evo_milkyway_gas_profile':       
        ThisRedshiftList=[0.0]        
        output = evo_milkyway_gas_profile(ThisRedshiftList)   
            
    if plot_to_run == 'evo_milkyway_stellar_profiles':       
        ThisRedshiftList=[0.0]        
        output = evo_milkyway_stellar_profiles(ThisRedshiftList)       
        
    if plot_to_run == 'test_H2_prescriptions':   
        ThisRedshiftList=[0.0]        
        output = test_H2_prescriptions(ThisRedshiftList)
      
    if plot_to_run == 'milkyway_gradients':      
        ThisRedshiftList=[0.1]
        output = milkyway_gradients(ThisRedshiftList)
    
    if plot_to_run == 'gas_metallicity_gradients_mass_bins':      
        ThisRedshiftList=[0.1]
        output = gas_metallicity_gradients_mass_bins(ThisRedshiftList)
            
    if plot_to_run == 'stellar_metallicity_gradients_mass_bins':      
        ThisRedshiftList=[0.0]
        output = stellar_metallicity_gradients_mass_bins(ThisRedshiftList)
            
    if plot_to_run == 'CALIFA_gradients_morph_types':      
        ThisRedshiftList=[0.0]
        output = CALIFA_gradients_morph_types(ThisRedshiftList)
            
    if plot_to_run == 'CALIFA_gradients_mass_bins':      
        ThisRedshiftList=[0.0]
        output = CALIFA_gradients_mass_bins(ThisRedshiftList)    
             
    if plot_to_run == 'MANGA_gradients_late_types':        
        ThisRedshiftList=[0.0]
        output = MANGA_gradients_late_types(ThisRedshiftList)
                                    
    if plot_to_run == 'SFR_gradients':      
        ThisRedshiftList=[0.1]
        output = SFR_gradients(ThisRedshiftList) 
            
    if plot_to_run == 'test_rings':       
        ThisRedshiftList=[0.0]
        output = test_rings(ThisRedshiftList) 
        
    return output