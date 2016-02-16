
# coding: utf-8

# # Python Plots for LGalaxies

# ## Import Libraries and Read Catalogs

# <p>Use functions read_snap or read_tree to read catalogs. These are both defined in procedures.py. In case of read_snap, SnapshotList will be returned containing the list of snapshots read (usefull to later select galaxies in a given redshift).<p>

# In[1]:

import numpy as np
get_ipython().magic('matplotlib inline')

import pandas as pd

get_ipython().magic('pylab inline')
#import seaborn as sns
#sns.set_style('darkgrid')

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.table import Table
from importlib import reload
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
from decimal import *
import sys
from scipy.ndimage import zoom 

import procedures
reload (procedures)
from procedures import *
import plots_input
reload (plots_input)
from plots_input import *

FirstFile = 5
LastFile =  5

Volume_MR = (BoxSize_MR**3.0) * (LastFile - FirstFile + 1) / MaxTreeFiles 
Volume_MRII = (BoxSize_MRII**3.0) * (LastFile - FirstFile + 1) / MaxTreeFiles 

if CatalogType=='snap':
    from LGalaxies_Henriques2015a_struct import LGalaxiesStruct
    from LGalaxies_Henriques2015a_struct import PropertiesToRead
    #from LGalaxies_Henriques2015a_Rings_struct import LGalaxiesStruct
    #from LGalaxies_Henriques2015a_Rings_struct import PropertiesToRead
    (G_MR, SnapshotList) = read_snap(DirName_MR,FirstFile,LastFile,
                     PropertiesToRead,LGalaxiesStruct,RedshiftsToRead,FullRedshiftList)
    
if CatalogType=='tree':    
    from LGalaxies_tree_Henriques2015a_struct import LGalaxiesStruct
    from LGalaxies_tree_Henriques2015a_struct import PropertiesToRead_tree    
    (G_MR) = read_tree(DirName_MR,FirstFile,LastFile,
                     PropertiesToRead_tree,LGalaxiesStruct)    
     
    SnapshotList = np.zeros(len(RedshiftList),dtype=np.int32)
    for ii in range(0,len(RedshiftList)):                  
        G0=G_MR[ np.rint(G_MR['Redshift']*100.) == RedshiftList[ii]*100. ]             
        SnapshotList[ii]=G0['SnapNum'][0]
#endif      

print('reading done\n')
#print (np.log10(G_MR['DiskMass'][0:99]*1.e10))
#print (np.log10(G_MR['BulgeMass'][0:99]*1.e10))
#print (np.log10(G_MR['StellarMass'][0:99]*1.e10))
#print (np.log10(G_MR['MetalsStellarMass'][0:99]*1.e10))
#print (G_MR[0:5])
#help(G_MR)


plt.rcParams.update({'font.size': 18, 'xtick.labelsize': 18, 'ytick.labelsize': 18, 'axes.linewidth': 2, 
                     'xtick.major.size': 6, 'xtick.major.width': 1.5, 
                     'ytick.major.size': 6, 'ytick.major.width': 1.5, 
                     'xtick.minor.size': 3, 'xtick.minor.width': 1.,                   
                     'ytick.minor.size': 3, 'ytick.minor.width': 1.})


# ## Plots

# In[2]:

import plots
reload (plots)


# In[3]:


with PdfPages('./fig/plots.pdf') as pdf:  
    import procedures
    reload (procedures)
    from procedures import *
    import plots_input
    reload (plots_input)
    from plots_input import *
    import plots
    reload (plots)

    print('Doing SatFraction')
    
    ThisRedshiftList=[0.0,1.0,2.0,3.0]  
    
    xlim=[8.5,11.5]
    ylim=[0., 1.] 
           
    plot_color=['red','purple']        
    plt.rcParams.update({'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 
                             'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0})
    
    fig = plt.figure(figsize=(15,4))
    grid = gridspec.GridSpec(1, 5)
    grid.update(wspace=0.0, hspace=0.0)
    
    for ii in range (0,len(ThisRedshiftList)):
           
        char_redshift="%0.2f" % ThisRedshiftList[ii]
            
        subplot=plt.subplot(grid[ii])
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        
        xlab='$\mathrm{log_{10}}(M_*[h^{-2}M_{\odot}])$'      
        if ii==0:
            ylab='Satellite Fraction'
        else:
            ylab=''      
        subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)
        
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
        subplot.yaxis.set_minor_locator(MultipleLocator(0.1))
        
        if ii>0:
            plt.tick_params(axis='y', which='both', left='on', labelleft='off')
        
        bin=0.25        
        Mass_arr=np.arange(xlim[0],xlim[1],bin)
        SatFraction=np.zeros(len(Mass_arr),dtype=np.float32)
        HaloSatFraction=np.zeros(len(Mass_arr),dtype=np.float32)
        
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii)        
        G0_MR=G_MR[sel]   
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
        Type=G0_MR['Type']
        
        for ll in range(0,len(Mass_arr)):
            sel_sat=G0_MR[(Type>0) & (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)]
            sel_halosat=G0_MR[(Type==1) & (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)]
            sel_all=G0_MR[(StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)]
                
            if len(sel_all)>0.:
                SatFraction[ll]=float(len(sel_sat))/float(len(sel_all))
                HaloSatFraction[ll]=float(len(sel_halosat))/float(len(sel_all))
            else:               
                SatFraction[ll]=0.
                HaloSatFraction[ll]=0.
                      
        subplot.plot(Mass_arr, SatFraction, color='red', linestyle='-', linewidth=2) 
        subplot.plot(Mass_arr, HaloSatFraction, color='green', linestyle='-', linewidth=2) 
       
        #labels
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.55, 
                    color='black', xlog=0, ylog=0, label='z='+char_redshift[:-1], 
                    fontsize=14, fontweight='normal') 
    
        if(ii==0):
            plot_label (subplot, 'label', xlim, ylim, 
                    x_percentage=0.15, y_percentage=0.90, color='black', xlog=0, ylog=0, 
                    label='galaxies', fontsize=13, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim,
                    x_percentage=0.04, y_percentage=0.92, color='red', x2_percentage=0.13, 
                    xlog=0, ylog=0, linestyle='-', linewidth=2)
            plot_label (subplot, 'label', xlim, ylim, 
                    x_percentage=0.15, y_percentage=0.8, color='black', xlog=0, ylog=0, 
                    label='haloes', fontsize=13, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim,
                    x_percentage=0.04, y_percentage=0.82, color='green', x2_percentage=0.13, 
                    xlog=0, ylog=0, linestyle='-', linewidth=2)
            
    plt.tight_layout()
    plt.savefig('./fig/plots_sat_fraction.pdf')
    pdf.savefig()
    plt.close()
    
    if opt_stellar_mass_vs_halo_mass==1:
        print('Doing SMHM')
        from plots import stellar_mass_vs_halo_mass
        ThisRedshiftList=[0.0]        
        stellar_mass_vs_halo_mass(G_MR, ThisRedshiftList, pdf)
            
    if opt_stellar_mass_function==1:
        print('Doing SMF')
        from plots import stellar_mass_function
        ThisRedshiftList=[0.0,1.0,2.0,3.0]        
        stellar_mass_function(G_MR, Volume_MR, ThisRedshiftList, pdf)
    
    if opt_redfraction_color_cut==1:
        print('Doing redfraction_color_cut')
        from plots import redfraction_color_cut
        ThisRedshiftList=[0.0,0.4,1.0,2.0,3.0]        
        redfraction_color_cut(G_MR, ThisRedshiftList, pdf)
    
    if opt_metals_vs_stellarmass==1:
        print('Doing metals_vs_stellarmass')
        from plots import metals_vs_stellarmass
        ThisRedshiftList=[0.1,3.]
        metals_vs_stellarmass(G_MR, ThisRedshiftList, pdf)
        
    if opt_BHBM==1:
        print('Doing BHBM')
        from plots import BHBM
        ThisRedshiftList=[0.0]        
        BHBM(G_MR, ThisRedshiftList, pdf)    
        
    if opt_SFRF==1:
        print('Doing SFRF')
        from plots import SFRF
        ThisRedshiftList=[0.0]        
        SFRF(G_MR, Volume_MR, ThisRedshiftList, pdf)
        
    if opt_gas_fraction==1:
        print('Doing gas_fraction')
        from plots import gas_fraction
        ThisRedshiftList=[0.0]        
        gas_fraction(G_MR, ThisRedshiftList, pdf)
        
    if opt_HI_MF==1:
        print('Doing HI_MF')
        from plots import HI_MF
        ThisRedshiftList=[0.0]        
        HI_MF(G_MR, Volume_MR, ThisRedshiftList, pdf)
        
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
                        
    if opt_gas_metallicity_gradients==1:
        print('Doing gas_metallicity_gradients')
        from plots import gas_metallicity_gradients
        ThisRedshiftList=[0.1]
        gas_metallicity_gradients(G_MR, ThisRedshiftList, pdf) 
        
    if opt_SFR_gradients==1:
        print('Doing SFR_gradients')
        from plots import SFR_gradients
        ThisRedshiftList=[0.1]
        SFR_gradients(G_MR, ThisRedshiftList, pdf) 
    
    if opt_bluck_red_fractions==1:
        print('Doing bluck_red_fractions')
        from plots import bluck_red_fractions
        ThisRedshiftList=[0.0]
        bluck_red_fractions(G_MR, ThisRedshiftList, pdf)             
        
    if opt_test_plots==1:
        print('Doing test_plots')
        from plots import test_plots     
        test_plots(G_MR, SnapshotList, pdf)
    
    print('')
    print('All plots done')
        
#end with PdfPages('./fig/plots.pdf') as pdf: 


# ## Animations

# In[5]:

import animations
reload (animations)
    
if opt_anime_mass_gr==1:   
    from animations import anime_mass_gr
    anime_mass_gr(G_MR)


# In[ ]:



