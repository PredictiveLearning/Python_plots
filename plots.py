''' 
stellar_mass_vs_halo_mass
stellar_mass_vs_halo_mass_fractional
stellar_mass_function
opt_stellar_mass_function_z0_overplot
opt_stellar_mass_function_feedback_overplot
stellar_mass_vs_halo_mass
redfraction_color_cut
metals_vs_stellarmass
BHBM
SFRF
gas_fraction
HI_MF
sfr_vs_stellar_masssfh_b
ur_vs_r
UVJ_colour
morphology_vs_stellarmass
sizes_vs_stellarmass
    
cooling_heating    
    
bluck_red_fractions
sat_fraction
BHmass_in_radio
sfr_massive_galaxies
misc_plots()
test_resolution_rings


simple_tree_map
full_tree_map
'''

import numpy as np
import pandas as pd
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
from scipy import interpolate
from scipy.stats import binned_statistic_2d
import sys
from scipy.ndimage import zoom
import os.path
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from matplotlib.patches import Ellipse, Circle
from matplotlib.collections import PatchCollection

import procedures
reload (procedures)
from procedures import *
import plots_input
reload (plots_input)
from plots_input import *


def stellar_mass_vs_halo_mass(G_MR, ThisRedshiftList, pdf):

    ylim=[2.0,12.5]
    xlim=[2.5, 15.]
    
    fig = plt.figure(figsize=(10,10))
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim) 

    ylab='$\mathrm{log_{10}}(M_*[M_{\odot}])$'       
    xlab='$\mathrm{log_{10}}(M_{200c}[M_{\odot}])$'
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)
                  
    for ii in range (0,len(ThisRedshiftList)):
  
        #MODEL
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)
        
        G0_MR=G_MR[sel] 
        G0_MR=np.random.choice(G0_MR, size=10000)        
        G0_MR=G0_MR[G0_MR['StellarMass']>0.]
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])          
        HaloMass=np.log10(G0_MR['Mvir']*1.e10/Hubble_h) 
        sel=G0_MR['Type']==2
        subplot.scatter(HaloMass[sel],StellarMass[sel],s=5, color='red')        
        sel=G0_MR['Type']==1
        subplot.scatter(HaloMass[sel],StellarMass[sel],s=5, color='green')
        sel=G0_MR['Type']==0
        subplot.scatter(HaloMass[sel],StellarMass[sel],s=5, color='blue')
        
    #endfor
    
    plt.tight_layout()
    plt.savefig('./fig/plots_smhm.pdf')
    pdf.savefig()
    plt.close()
    
    
    
def stellar_mass_vs_halo_mass_fractional(G_MR, G_MRII, ThisRedshiftList, pdf):

    xlim=[11., 14.5]
    ylim=[-3.,0.]    
        
        
    model_to_print='Hen15_only_AGN'    
    #TO WRITE OUT FOR A GIVEN MODEL
    #PLOT IS ACTUALLY MADE FROM FILES PREVIOUSLY READ
    for ii in range (0,len(ThisRedshiftList)):
  
        #MODEL
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)
        
        G0_MR=G_MR[sel]             
        G0_MR=G0_MR[(G0_MR['StellarMass']>0.0) & (G0_MR['Type']==0)]
        #G0_MR=G0_MR[(G0_MR['StellarMass']>0.0)]
        G0_MR=np.random.choice(G0_MR, size=2000)   
        log_StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])          
        log_HaloMass=np.log10(G0_MR['Mvir']*1.e10/Hubble_h) 
        log_fraction=log_StellarMass-log_HaloMass
        #log_fraction=log_StellarMass
        SSFR=G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h)
        
        #WRITE OUTPUT
        
        #PASSIVE
        sel=np.log10(SSFR)<-11.
        #subplot.scatter(log_HaloMass[sel],log_fraction[sel],s=5, color='red') 
        x_axis=log_HaloMass[sel]
        y_axis=log_fraction[sel]
        fa = open(Datadir+"SMHM_"+model_to_print+"_passive_z0.0.txt", "w")
        fa.write("%d\n" % len(x_axis))
        for kk in range (0,len(x_axis)):               
            fa.write("%0.2f " % x_axis[kk] + "%0.2f\n" % y_axis[kk])         
        fa.close()   
        
        #SF
        sel=np.log10(SSFR)>-11.
        #subplot.scatter(log_HaloMass[sel],log_fraction[sel],s=5, color='blue') 
        x_axis=log_HaloMass[sel]
        y_axis=log_fraction[sel]
        fa = open(Datadir+"SMHM_"+model_to_print+"_SF_z0.0.txt", "w")
        fa.write("%d\n" % len(x_axis))
        for kk in range (0,len(x_axis)):               
            fa.write("%0.2f " % x_axis[kk] + "%0.2f\n" % y_axis[kk])         
        fa.close() 
        
        #MEDIAN
        bin=0.25
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]             
        G0_MR=G0_MR[(G0_MR['StellarMass']>0.0) & (G0_MR['Type']==0)]       
        #G0_MR=G0_MR[(G0_MR['StellarMass']>0.0)]       
        log_StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])       
        log_HaloMass=np.log10(G0_MR['Mvir']*1.e10/Hubble_h)  
        log_fraction=log_StellarMass-log_HaloMass
        #log_fraction=log_StellarMass
        (x_binned, median, mean, pc16, pc84,rms)= median_and_percentiles(bin, xlim[0], xlim[1],log_HaloMass,log_fraction)    
        x_axis=x_binned
        y_axis=median
        fa = open(Datadir+"SMHM_"+model_to_print+"_median_z0.0.txt", "w")
        fa.write("%d\n" % len(x_axis))
        for kk in range (0,len(x_axis)):               
            fa.write("%0.2f " % x_axis[kk] + "%0.2f\n" % y_axis[kk])     
        fa.close() 
        
        y_axis=pc16
        fa = open(Datadir+"SMHM_"+model_to_print+"_pc16_z0.0.txt", "w")
        fa.write("%d\n" % len(x_axis))
        for kk in range (0,len(x_axis)):               
            fa.write("%0.2f " % x_axis[kk] + "%0.2f\n" % y_axis[kk])     
        fa.close()
        
        y_axis=pc84
        fa = open(Datadir+"SMHM_"+model_to_print+"_pc84_z0.0.txt", "w")
        fa.write("%d\n" % len(x_axis))
        for kk in range (0,len(x_axis)):               
            fa.write("%0.2f " % x_axis[kk] + "%0.2f\n" % y_axis[kk])     
        fa.close()
               
        
    
    #PLOTS   
  
    fig = plt.figure(figsize=(two_two_size_small[0],two_two_size_small[1]))
    grid = gridspec.GridSpec(2, 2)
    grid.update(wspace=0.0, hspace=0.0)

    
    model_names=['Hen15_no_feedback','Hen15_only_AGN','Hen15_only_SN','Hen15']
    model_label=['no feedback','only AGN','only SN','Hen15']
    
    for ii in range(0,4):
        
        subplot=plt.subplot(grid[ii])

        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        if ii==2 or ii== 3: 
            xlab='$\mathrm{log_{10}}(M_{200c}[M_{\odot}])$' 
            subplot.set_xlabel(xlab, fontsize=14)
      
        if ii==0 or ii== 2:
            ylab='$\mathrm{log_{10}}(M_*/M_{200c})$'
            subplot.set_ylabel(ylab, fontsize=14)
                       
    
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
        subplot.yaxis.set_minor_locator(MultipleLocator(0.25))
  
        if ii==1 or ii==3:
            plt.tick_params(axis='y', which='both', left='on', labelleft='off')
        if ii==0 or ii==1:    
            plt.tick_params(axis='x', which='both', bottom='on', labelbottom='off')
    
    
        fa = Datadir+"SMHM_"+model_names[ii]+"_passive_z0.0.txt"
        (x_axis,y_axis)=read_file(fa)    
        subplot.scatter(x_axis,y_axis,s=5, color='red') 
        
        fa = Datadir+"SMHM_"+model_names[ii]+"_SF_z0.0.txt"
        (x_axis,y_axis)=read_file(fa)    
        subplot.scatter(x_axis,y_axis,s=5, color='blue') 
            
        fa = Datadir+"SMHM_"+model_names[ii]+"_median_z0.0.txt"
        (x_axis,y_axis)=read_file(fa)
        subplot.plot(x_axis, y_axis,color='black', linewidth=2)   
         
        fa = Datadir+"SMHM_"+model_names[ii]+"_pc16_z0.0.txt"
        (x_axis,y_axis)=read_file(fa)
        subplot.plot(x_axis, y_axis,color='black', linewidth=2, linestyle='--') 
         
        fa = Datadir+"SMHM_"+model_names[ii]+"_pc84_z0.0.txt"
        (x_axis,y_axis)=read_file(fa)
        subplot.plot(x_axis, y_axis,color='black', linewidth=2, linestyle='--')  
       
    
        #MOSTER2012
        '''z=ThisRedshiftList[ii]      
        M_10=11.59
        M_11=1.195
        N_10=0.0351
        N_11=-0.0247
        beta_10=1.376
        beta_11=-0.826
        gamma_10=0.608
        gamma_11=0.329
        moster_mvir=10**np.arange(xlim[0]-1.0,xlim[1]+1.0,0.05)  
        M_1=10**(M_10+M_11*(z/(z+1)))
        N=N_10+N_11*(z/(z+1)) 
        beta=beta_10+beta_11*(z/(z+1))
        gamma=gamma_10+gamma_11*(z/(z+1))
        moster_mstar=moster_mvir*2.*N * ( (moster_mvir/M_1)**(-1.*beta) + (moster_mvir/M_1)**(gamma) )**(-1)      
        subplot.plot(np.log10(moster_mvir), np.log10(moster_mstar/moster_mvir),color='green', linewidth=2)''' 
        
        #Behroozi2013
        z=ThisRedshiftList[0]   
        a=1/(1+z)
        neu=np.exp(-4*a*a)
        log_epsilon=-1.777+(-0.006*(a-1)-0.*z)*neu-0.119*(a-1)
        log_M1=11.514+(-1.793*(a-1)-0.251*z)*neu
        alpha=-1.412+0.731*(a-1.)*neu
        delta=3.508+(2.608*(a-1)-0.043*z)*neu
        gamma=0.316+(1.319*(a-1)+0.279*z)*neu 
        
        '''z=ThisRedshiftList[ii]   
        a=1/(1+z)
        neu=np.exp(-4*a*a)
        log_epsilon=-1.777
        log_M1=11.514
        alpha=-1.412
        delta=3.508
        gamma=0.316'''
       
        x=0.
        f_0=-np.log10(10**(alpha*x) + 1) + delta * ((np.log10(1 + np.exp(x)))**gamma)/(1 + np.exp(10**(-x)))             
        log_Mh=np.arange(8.0,16.0,0.01) 
        x=log_Mh-log_M1
        f_log_Mh=-np.log10(10**(alpha*x) + 1) + delta * ((np.log10(1 + np.exp(x)))**gamma)/(1 + np.exp(10**(-x)))      
        log_mstar=log_epsilon+log_M1+f_log_Mh-f_0     
        subplot.plot(log_Mh, log_mstar-log_Mh,color='limegreen', linewidth=3)  
        
        
        
    
        x_arr=np.arange(xlim[0],xlim[1]+0.05,0.05)           
        y_arr=x_arr*0.+np.log10(0.155)      
        subplot.plot(x_arr,y_arr,color='black', linestyle='--', linewidth=3)         
        
        if(ii==0):
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.15, y_percentage=0.17, 
                        color='black', xlog=0, ylog=0, label='Model Median', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.06, y_percentage=0.19, x2_percentage=0.12, 
                        color='black', xlog=0, ylog=0, linestyle='-', linewidth=2)  
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.15, y_percentage=0.10, 
                        color='black', xlog=0, ylog=0, label='Behroozi et al. (2013)', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.06, y_percentage=0.12, x2_percentage=0.12, 
                        color='limegreen', xlog=0, ylog=0, linestyle='-', linewidth=2) 
        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.85, y_percentage=0.6, 
                    color='black', xlog=0, ylog=0, label='$f_b$', 
                    fontsize=15, fontweight='normal')
        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.1, y_percentage=0.85, 
                    color='black', xlog=0, ylog=0, label=model_label[ii], 
                    fontsize=12, fontweight='normal')
        
    #endfor
    
    
    '''
    Berhoozi2010
    z=0.
    M_10=12.35
    M_1a=0.28
    mstar_00=10.72
    mstar_0a=0.55
    beta_0=0.44
    beta_a=0.18
    delta_0=0.57
    delta_a=0.17
    gamma_0=1.56
    gamma_a=2.51
    log_M1=M_10+M_1a*(-z/(z+1))
    log_mstar0=mstar_00+mstar_0a*(-z/(z+1))
    beta=beta_0+beta_a*(-z/(z+1))
    delta=delta_0+delta_a*(-z/(z+1))
    gamma=gamma_0+gamma_a*(-z/(z+1))        
    log_mstar=np.arange(8.0,12.0,0.01)        
    log_mstar_mstar0=log_mstar-log_mstar0     
    log_Mh=log_M1+beta*(log_mstar_mstar0) + (10**log_mstar_mstar0)**delta/(1+(10**log_mstar_mstar0)**(-gamma))-0.5    
    subplot.plot(log_Mh, log_mstar-log_Mh,color='green', linewidth=2)  '''
       
    
    plt.tight_layout()
    plt.savefig('./fig/plots_smhm_fractional.pdf')
    pdf.savefig()
    plt.close()
    

def stellar_mass_function(G_MR, Volume_MR, G_MRII, Volume_MRII, ThisRedshiftList, pdf):
           
    xlim=[8.0,12.5]
    ylim=[-6.5, 0.5]
    bin=0.25

    model_to_print='Hen15_no_feedback'
    
    fig = plt.figure(figsize=(two_two_size_small[0],two_two_size_small[1]))
    grid = gridspec.GridSpec(2, 2)
    grid.update(wspace=0.0, hspace=0.0)

    for i_z in range(0,len(ThisRedshiftList)):
        
        char_redshift="%0.2f" % ThisRedshiftList[i_z]
        
        subplot=plt.subplot(grid[i_z])

        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        if i_z==2 or i_z == 3: 
            xlab='$\mathrm{log_{10}}(M_*[M_{\odot}])$'
            subplot.set_xlabel(xlab, fontsize=14)
      
        if i_z==0 or i_z == 2:
            ylab='$\mathrm{log_{10}}(\phi [\mathrm{Mpc^{-3}} \mathrm{log_{10}}(M^{-1})])$'
            subplot.set_ylabel(ylab, fontsize=14)
                       
    
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
        subplot.yaxis.set_minor_locator(MultipleLocator(0.25))
  
        if i_z==1 or i_z==3:
            plt.tick_params(axis='y', which='both', left='on', labelleft='off')
        if i_z==0 or i_z==1:    
            plt.tick_params(axis='x', which='both', bottom='on', labelbottom='off')
           
    
        #OBSERVATIONS             
        file = MCMCdir + '/ObsConstraints/StellarMassFunction_z'+char_redshift+'.txt'        
        f = open(file, 'r')     
        line = int(f.readline())     
        obs = Table.read(file, format='ascii', data_start=1, data_end=line+1)
        
        obs_xbin=obs['col1']+(obs['col2']-obs['col1'])/2.
        asy_yerror = [np.log10(obs['col3']/(obs['col3']-obs['col4'])), 
                      np.log10((obs['col3']+obs['col4'])/obs['col3'])]       
        subplot.errorbar(obs_xbin-2.*np.log10(Hubble_h), np.log10(obs['col3'])+3.*np.log10(Hubble_h),
                         yerr=asy_yerror,fmt='o', markersize=5, ecolor='blue', color='blue')
        #sub = plt.subplot(111)
    
        #PREVIOUS MODELS
        RedshiftList_OldModels=[0.1,1.0,2.0,3.0]
        if do_previous_model1==1:           
            char_old_redshift="%0.2f" % RedshiftList_OldModels[i_z]
            file = file_previous_model1+'_smf_z'+char_old_redshift+'.txt' 
            model = Table.read(file, format='ascii')
            subplot.plot(model['col1']-2.*np.log10(Hubble_h_WMAP7), model['col2']+3.*np.log10(Hubble_h_WMAP7),
                         color='red',linestyle=linestyle_previous_model1, linewidth=2)
      
        if do_previous_model2==1:           
            file = file_previous_model2+'_smf_z'+char_redshift+'.txt' 
            model = Table.read(file, format='ascii')
            subplot.plot(model['col1']-2.*np.log10(Hubble_h),model['col2']+3.*np.log10(Hubble_h),
                         color='red',linestyle=linestyle_previous_model2, linewidth=2)
      
    
        #MODEL
        #MR
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)
        
        G0_MR=G_MR[sel]   
        G0_MR=G0_MR[G0_MR['StellarMass']>0.]
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[i_z])
        #StellarMass=np.log10((G0_MR['MassFromInSitu']+G0_MR['MassFromMergers']+G0_MR['MassFromBursts'])*1.e10*Hubble_h)
        
        bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
        hist_MR=np.histogram(StellarMass, bins=bin_arr, range=(xlim[0],xlim[1]))   
       
        #MRII
        if(MRII==1):
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, i_z, FullSnapshotList_MRII)
        
            G0_MRII=G_MRII[sel]   
            G0_MRII=G0_MRII[G0_MRII['StellarMass']>0.]
            StellarMass=stellar_mass_with_err(G0_MRII, Hubble_h, ThisRedshiftList[i_z])
        
            bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
            hist_MRII=np.histogram(StellarMass, bins=bin_arr, range=(xlim[0],xlim[1]))   
           
        
        #join MR+MRII & plot     
        if(MRII==1):
            cut_MR_MRII=10.
            (x_axis,y_axis)=plot_joint_MR_MRII(hist_MR, hist_MRII, cut_MR_MRII, Volume_MR, Volume_MRII, 
                                               bin, subplot, color='red',linewidth=2, linestyle='-')
        else:
            x_axis=hist_MR[1][0:len(hist_MR[1][:])-1]+bin/2.           
            hist_MR=hist_MR[0]       
            y_axis=np.log10(hist_MR/(Volume_MR*bin))
            subplot.plot(x_axis,y_axis, color='red', linewidth=2, linestyle='-') 
         
        
        #WRITE OUTPUT
        fa = open(Datadir+"SMF_"+model_to_print+"_z"+char_redshift+".txt", "w")
        fa.write("%d\n" % len(x_axis))
        for kk in range (0,len(x_axis)):               
            fa.write("%0.2f " % x_axis[kk] + "%0.2f\n" % y_axis[kk])         
        fa.close()    
        
    
        #MCMC sample
        if opt_plot_MCMC_sample==1:
            file = MCMCSampledir + 'mcmc_plus_obs_StellarMassFunction_z'+char_redshift+'.txt' 
            if os.path.isfile(file):
                obs = Table.read(file, format='ascii')      
                subplot.plot(obs['col1']-2.*np.log10(Hubble_h),np.log10(obs['col4'])+3.*np.log10(Hubble_h), 
                             color='black', linewidth=2)               
                if i_z==len(ThisRedshiftList)-1:
                    plot_label (subplot, 'label', xlim, ylim, 
                        x_percentage=0.55, y_percentage=0.85, color='black', xlog=0, ylog=0, 
                        label='MCMC sample', fontsize=13, fontweight='normal') 
                    plot_label (subplot, 'line', xlim, ylim,
                        x_percentage=0.44, y_percentage=0.87, color='black', x2_percentage=0.53, 
                        xlog=0, ylog=0, linestyle='-', linewidth=2)
                
                
        #LABELS
        if i_z==0:
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.91, 
                        color='black', xlog=0, ylog=0, label='Observations used in MCMC', 
                        fontsize=13, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.935, 
                        color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15) 
            
        if i_z==len(ThisRedshiftList)-1:
            plot_label_three_models (subplot, xlim, ylim, position='bottom_left')
                                  
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.35, 
                    color='black', xlog=0, ylog=0, label='z='+char_redshift[:-1], 
                    fontsize=14, fontweight='normal')      
        if i_z==2:
             plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.075, 
                         color='black', xlog=0, ylog=0, label='SMF evo', 
                         fontsize=16, fontweight='normal') 
           
    #endfor

 
    plt.tight_layout()
    plt.savefig('./fig/plots_smf_evo.pdf')
    plt.savefig('./fig/HYW17_plots_smf_evo.pdf')
    pdf.savefig()
    plt.close()
#endif stellar_mass_function==1:


def stellar_mass_function_z0_overplot(G_MR, Volume_MR, G_MRII, Volume_MRII, ThisRedshiftList, pdf):
           
    xlim=[8.0,12.5]
    ylim=[-6.5, 1.0]
        
    plot_color=['red','purple']        
   
    fig = plt.figure(figsize=(6,5))
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1.0))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
    subplot.yaxis.set_minor_locator(MultipleLocator(0.5)) 
    
    xlab='$\mathrm{log_{10}}(M_*[M_{\odot}])$'      
    ylab='$\mathrm{log_{10}}(\phi [\mathrm{Mpc^{-3}} \mathrm{log_{10}}(M^{-1})])$'         
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)    
    
    
  
    #HALO MASS FUNCTION 
    bin=0.1
    for ii in range(0,len(ThisRedshiftList)):
        #MR
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)
        
        G0_MR=G_MR[sel]   
        G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (G0_MR['Type']==0)]   
        Mvir=np.log10(G0_MR['Mvir']*1.e10/Hubble_h*0.155)
        
        bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
        hist_MR=np.histogram(Mvir, bins=bin_arr, range=(xlim[0],xlim[1]))   
        
        #MRII
        if(MRII==1):
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, ii, FullSnapshotList_MRII)
        
            G0_MRII=G_MRII[sel]   
            G0_MRII=G0_MRII[(G0_MRII['StellarMass']>0.) & (G0_MRII['Type']==0)]
            Mvir=np.log10(G0_MRII['Mvir']*1.e10/Hubble_h*0.155)
            
            bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
            hist_MRII=np.histogram(Mvir, bins=bin_arr, range=(xlim[0],xlim[1]))   
           
        
        #join MR+MRII & plot        
        if(MRII==1):
            cut_MR_MRII=10.1
            (x_axis,y_axis)=plot_joint_MR_MRII(hist_MR, hist_MRII, cut_MR_MRII, Volume_MR, Volume_MRII, 
                                               bin, subplot, color='black',linewidth=2, linestyle=':')
        else:
            x_axis=hist_MR[1][0:len(hist_MR[1][:])-1]+bin/2.           
            hist_MR=hist_MR[0]       
            y_axis=np.log10(hist_MR/(Volume_MR*bin))
            subplot.plot(x_axis,y_axis, color='black', linewidth=2, linestyle=':')         
        
        #SECOND AXIS WITH MVIR
        bin=1.0            
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
        Mvir=np.log10(G0_MR['Mvir']*1.e10/Hubble_h)
        (x_binned,median_MR,mean,pc16,pc84,rms)=median_and_percentiles(bin,xlim[0],xlim[1],StellarMass,Mvir) 
      
        if(MRII==1):
            StellarMass=stellar_mass_with_err(G0_MRII, Hubble_h, ThisRedshiftList[ii])
            Mvir=np.log10(G0_MRII['Mvir']*1.e10/Hubble_h)
            (x_binned,median_MRII,mean,pc16,pc84,rms)=median_and_percentiles(bin,xlim[0],xlim[1],StellarMass,Mvir) 
       
        if(MRII==1):
            y_axis=np.concatenate((median_MRII[x_binned<=cut_MR_MRII],median_MR[x_binned>cut_MR_MRII]), axis=0) 
        else:
            y_axis=median_MR
       
        for ii in range(0,len(y_axis)):
            y_axis[ii] = int((y_axis[ii] * 10) + 0.5) / 10.0            
                   
                
                
        ax2 = subplot.twiny()        
        ax2.set_xbound(subplot.get_xbound())        
        ax2.set_xticks(x_binned)
        ax2.set_xticklabels(y_axis)
     
        xlab='$\mathrm{log_{10}}(<M_{\mathrm{vir}}>[M_{\odot}])$'  
        ax2.set_xlabel(xlab, fontsize=14)
        ax2.xaxis.set_label_position('top')  
    
    
    char_redshift="%0.2f" % 0.00         
          
    fa = Datadir+"SMF_Hen15_only_SN_z"+char_redshift+".txt"
    (x_axis,y_axis)=read_file(fa)    
    subplot.plot(x_axis,y_axis, color='blue', linewidth=2, linestyle='-')  
    
    fa = Datadir+"SMF_Hen15_only_AGN_z"+char_redshift+".txt"
    (x_axis,y_axis)=read_file(fa)    
    subplot.plot(x_axis,y_axis, color='brown', linewidth=2, linestyle='-')  
    
    fa = Datadir+"SMF_Hen15_no_feedback_z"+char_redshift+".txt"
    (x_axis,y_axis)=read_file(fa)    
    subplot.plot(x_axis,y_axis, color='black', linewidth=2, linestyle='-')  
        
    fa = Datadir+"SMF_Hen15_z"+char_redshift+".txt"
    (x_axis,y_axis)=read_file(fa)    
    subplot.plot(x_axis,y_axis, color='red', linewidth=2, linestyle='-')  
        
    #OBSERVATIONS             
    file = MCMCdir + '/ObsConstraints/StellarMassFunction_z'+char_redshift+'.txt'        
    f = open(file, 'r')     
    line = int(f.readline())     
    obs = Table.read(file, format='ascii', data_start=1, data_end=line+1)
        
    obs_xbin=obs['col1']+(obs['col2']-obs['col1'])/2.
    asy_yerror = [np.log10(obs['col3']/(obs['col3']-obs['col4'])), 
                  np.log10((obs['col3']+obs['col4'])/obs['col3'])]
    subplot.errorbar(obs_xbin-2.*np.log10(Hubble_h), np.log10(obs['col3'])+3.*np.log10(Hubble_h),
                     yerr=asy_yerror,fmt='o', markersize=5, ecolor='blue', color='blue')
    #sub = plt.subplot(111)
   

    #LABELS
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.91, 
                color='black', xlog=0, ylog=0, label='Observations used in MCMC', 
                fontsize=13, fontweight='normal') 
    plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.935, 
                color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15)   
        
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.12, y_percentage=0.45, 
                color='black', xlog=0, ylog=0, label='Hen15', fontsize=13, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, y_percentage=0.465, 
                color='red', x2_percentage=0.1, xlog=0, ylog=0, linestyle='-', linewidth=2)  
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.12, y_percentage=0.38, 
                color='black', xlog=0, ylog=0, label='only SN', fontsize=13, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, y_percentage=0.395, 
                color='blue', x2_percentage=0.1, xlog=0, ylog=0, linestyle='-', linewidth=2)
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.12, y_percentage=0.31, 
                color='black', xlog=0, ylog=0, label='only AGN', fontsize=13, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, y_percentage=0.325, 
                color='brown', x2_percentage=0.1, xlog=0, ylog=0, linestyle='-', linewidth=2)
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.12, y_percentage=0.24, 
                color='black', xlog=0, ylog=0, label='no feedback', fontsize=13, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, y_percentage=0.255, 
                color='black', x2_percentage=0.1, xlog=0, ylog=0, linestyle='-', linewidth=2)
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.12, y_percentage=0.17, 
                color='black', xlog=0, ylog=0, label=r'$f_b \times M_{\mathrm{Halo}}$', fontsize=13, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, y_percentage=0.185, 
                color='black', x2_percentage=0.1, xlog=0, ylog=0, linestyle=':', linewidth=2)
    
    plt.tight_layout()
    plt.savefig('./fig/plots_stellar_mass_function_z0_overplot.pdf')
    pdf.savefig()
    plt.close()

#end stellar_mass_function_z0_overplot

  
def stellar_mass_function_allz_overplot(G_MR, Volume_MR, G_MRII, Volume_MRII, ThisRedshiftList, pdf):
           
    xlim=[8.0,12.5]
    ylim=[-6.5, -0.5]
    bin=0.25

    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
    subplot.yaxis.set_minor_locator(MultipleLocator(0.5)) 
    
    xlab='$\mathrm{log_{10}}(M_*[M_{\odot}])$'
    ylab='$\mathrm{log_{10}}(\phi [\mathrm{Mpc^{-3}} \mathrm{log_{10}}(M^{-1})])$'
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14) 
    
    z_model=[0.4,1.0,2.0,3.0,4.0,5.0] 
    z_obs=[    0.35,  0.95,   1.3,  1.75, 2.25,   2.75,  3.25,   4.0,   5.0]
    #z_obs=[3.5]
    M_star = [10.78, 10.56, 10.62, 10.51, 10.60, 10.59, 10.83, 11.10, 11.30]
    alpha_1= [-1.38, -1.31, -1.28, -1.28, -1.57, -1.67, -1.76, -1.98, -2.11]
    phi_1  = [1.187, 1.428, 1.069, 0.969, 0.295, 0.228, 0.090, 0.016, 0.003]
    alpha_2= [-0.43,  0.51,  0.29,  0.82,  0.07, -0.08,    1.,    1.,    1.]
    phi_2  = [1.92,   2.19,  1.21,  0.64,  0.45,  0.21,    0.,    0.,    0.]
    
    obs_color=['grey','violet','darkorange','orange','yellow','green','blue','purple','pink','red']
    model_color=['grey','darkorange','green','purple','pink','red']
    
    for i_z in range(0,10):
    #for i_z in range(6,8):
        
        #xx=np.arange(xlim[0],xlim[1]+bin,0.1)
        #S1=phi_1[i_z]*0.001*(10**(xx-M_star[i_z]))**(alpha_1[i_z]+1.)
        #S2=phi_2[i_z]*0.001*(10**(xx-M_star[i_z]))**(alpha_2[i_z]+1.)
        #LLF1 = (S1+S2) * np.exp(-10**(xx-M_star[i_z])) * np.log10(10)
        #subplot.plot(xx,np.log10(LLF1), color='blue', linewidth=2, linestyle='-')                   
          
        
        file = Datadir+"/Davidzon/mf_mass2b_fl5b_tot_Vmax"+"%0.0d" % i_z + ".dat"
        Davidzon = Table.read(file, format='ascii')  
      
        #subplot.plot(Davidzon['col1'], Davidzon['col2'], color='black', linewidth=2, linestyle='-')
        subplot.errorbar(Davidzon['col1'], Davidzon['col2'],[Davidzon['col4'],Davidzon['col3']],
                         fmt='o', markersize=5, ecolor='blue', color=obs_color[i_z])
        
    for i_z in range(0,len(ThisRedshiftList)):
        
        char_redshift="%0.2f" % ThisRedshiftList[i_z]
      
                            
        #MODEL
        #MR
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)
        
        G0_MR=G_MR[sel]   
        G0_MR=G0_MR[G0_MR['StellarMass']>0.]
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[i_z])
        #StellarMass=np.log10((G0_MR['MassFromInSitu']+G0_MR['MassFromMergers']+G0_MR['MassFromBursts'])*1.e10*Hubble_h)
        
        bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
        hist_MR=np.histogram(StellarMass, bins=bin_arr, range=(xlim[0],xlim[1]))   
       
        #MRII
        if(MRII==1):
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, i_z, FullSnapshotList_MRII)
        
            G0_MRII=G_MRII[sel]   
            G0_MRII=G0_MRII[G0_MRII['StellarMass']>0.]
            StellarMass=stellar_mass_with_err(G0_MRII, Hubble_h, ThisRedshiftList[i_z])
        
            bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
            hist_MRII=np.histogram(StellarMass, bins=bin_arr, range=(xlim[0],xlim[1]))   
           
        
        #join MR+MRII & plot     
        if(MRII==1):
            cut_MR_MRII=10.
            (x_axis,y_axis)=plot_joint_MR_MRII(hist_MR, hist_MRII, cut_MR_MRII, Volume_MR, Volume_MRII, 
                                               bin, subplot, color=model_color[i_z],linewidth=2, linestyle='-')
        else:
            x_axis=hist_MR[1][0:len(hist_MR[1][:])-1]+bin/2.           
            hist_MR=hist_MR[0]       
            y_axis=np.log10(hist_MR/(Volume_MR*bin))
            subplot.plot(x_axis,y_axis, color=model_color[i_z], linewidth=2, linestyle='-') 
           
        
        #fa = open(Datadir+"SMF_Hen15_total_z"+char_redshift+".txt", "w")
        #fa.write("%d\n" % len(x_axis))
        #for kk in range (0,len(x_axis)):               
        #    fa.write("%0.2f " % x_axis[kk] + "%0.2f\n" % y_axis[kk])         
        #fa.close()
        
        #fa = Datadir+"SMF_Hen15_total_z"+char_redshift+".txt"
        #(x_axis,y_axis)=read_file(fa)    
        #subplot.plot(x_axis,y_axis, color='black', linewidth=2, linestyle='-') 
        
        #LABELS
        if i_z==0:
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.91, 
                        color='black', xlog=0, ylog=0, label='Davidzon et al. (2017)', 
                        fontsize=13, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.935, 
                        color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15) 
      
                                  
        #plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.35, 
        #            color='black', xlog=0, ylog=0, label='z='+char_redshift[:-1], 
         #           fontsize=14, fontweight='normal')      
        if i_z==2:
             plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.075, 
                         color='black', xlog=0, ylog=0, label='SMF evo', 
                         fontsize=16, fontweight='normal') 
          
    #endfor

 
    plt.tight_layout()
    plt.savefig('./fig/plots_stellar_mass_function_allz_overplot.pdf')
    pdf.savefig()
    plt.close()
    
    
    
    
    
    
    #Active    
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
    subplot.yaxis.set_minor_locator(MultipleLocator(0.5)) 
    
    xlab='$\mathrm{log_{10}}(M_*[M_{\odot}])$'
    ylab='$\mathrm{log_{10}}(\phi [\mathrm{Mpc^{-3}} \mathrm{log_{10}}(M^{-1})])$'
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14) 
    
    for i_z in range(0,10):
   
        file = Datadir+"/Davidzon/mf_mass2b_fl5b_act_Vmax"+"%0.0d" % i_z + ".dat"
        Davidzon = Table.read(file, format='ascii')        
        subplot.errorbar(Davidzon['col1'], Davidzon['col2'],[Davidzon['col4'],Davidzon['col3']],
                         fmt='o', markersize=5, ecolor='blue', color=obs_color[i_z])
        
    for i_z in range(0,len(ThisRedshiftList)):
        
        char_redshift="%0.2f" % ThisRedshiftList[i_z]
      
                            
        #MODEL
        #MR
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)
        
        G0_MR=G_MR[sel]   
    
        color_UV=G0_MR['MagDust'][:,0]-G0_MR['MagDust'][:,2]  
        color_VJ=G0_MR['MagDust'][:,2]-G0_MR['MagDust'][:,7]     
        G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (color_UV < minimum_y_color_cut[i_z+1]) | 
                    (color_UV < (color_VJ*slope_color_cut[i_z+1])+offset_color_cut[i_z+1])]
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[i_z])    
        #StellarMass=np.log10((G0_MR['MassFromInSitu']+G0_MR['MassFromMergers']+G0_MR['MassFromBursts'])*1.e10*Hubble_h)
        
        bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
        hist_MR=np.histogram(StellarMass, bins=bin_arr, range=(xlim[0],xlim[1]))   
       
        #MRII
        if(MRII==1):
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, i_z, FullSnapshotList_MRII)
        
            G0_MRII=G_MRII[sel]        
            color_UV=G0_MRII['MagDust'][:,0]-G0_MRII['MagDust'][:,2]  
            color_VJ=G0_MRII['MagDust'][:,2]-G0_MRII['MagDust'][:,7]       
            G0_MRII=G0_MRII[(G0_MRII['StellarMass']>0.) & (color_UV < minimum_y_color_cut[i_z+1]) | 
                            (color_UV < (color_VJ*slope_color_cut[i_z+1])+offset_color_cut[i_z+1])]
            
            StellarMass=stellar_mass_with_err(G0_MRII, Hubble_h, ThisRedshiftList[i_z])
        
            bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
            hist_MRII=np.histogram(StellarMass, bins=bin_arr, range=(xlim[0],xlim[1]))   
        
        #join MR+MRII & plot     
        if(MRII==1):
            cut_MR_MRII=10.
            (x_axis,y_axis)=plot_joint_MR_MRII(hist_MR, hist_MRII, cut_MR_MRII, Volume_MR, Volume_MRII, 
                                               bin, subplot, color=model_color[i_z],linewidth=2, linestyle='-')
        else:
            x_axis=hist_MR[1][0:len(hist_MR[1][:])-1]+bin/2.           
            hist_MR=hist_MR[0]       
            y_axis=np.log10(hist_MR/(Volume_MR*bin))
            subplot.plot(x_axis,y_axis, color=model_color[i_z], linewidth=2, linestyle='-') 
           
        #fa = open(Datadir+"SMF_Hen15_active_z"+char_redshift+".txt", "w")
        #fa.write("%d\n" % len(x_axis))
        #for kk in range (0,len(x_axis)):               
        #    fa.write("%0.2f " % x_axis[kk] + "%0.2f\n" % y_axis[kk])         
        #fa.close()
        
        #fa = Datadir+"SMF_Hen15_active_z"+char_redshift+".txt"
        #(x_axis,y_axis)=read_file(fa)    
        #subplot.plot(x_axis,y_axis, color='black', linewidth=2, linestyle='-')
        
        #LABELS
        if i_z==0:
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.91, 
                        color='black', xlog=0, ylog=0, label='Davidzon et al. (2017)', 
                        fontsize=13, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.935, 
                        color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15) 
                               
        #plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.35, 
        #            color='black', xlog=0, ylog=0, label='z='+char_redshift[:-1], 
         #           fontsize=14, fontweight='normal')      
        if i_z==2:
             plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.075, 
                         color='black', xlog=0, ylog=0, label='Active', 
                         fontsize=16, fontweight='normal') 
           
    #endfor

 
    plt.tight_layout()
    plt.savefig('./fig/plots_stellar_mass_function_active_allz_overplot.pdf')
    pdf.savefig()
    plt.close()
    
    
    
    
    #PASSIVE    
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
    subplot.yaxis.set_minor_locator(MultipleLocator(0.5)) 
    
    xlab='$\mathrm{log_{10}}(M_*[M_{\odot}])$'
    ylab='$\mathrm{log_{10}}(\phi [\mathrm{Mpc^{-3}} \mathrm{log_{10}}(M^{-1})])$'
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14) 
    
    for i_z in range(0,10):
   
        file = Datadir+"/Davidzon/mf_mass2b_fl5b_pas_Vmax"+"%0.0d" % i_z + ".dat"
        Davidzon = Table.read(file, format='ascii')        
        subplot.errorbar(Davidzon['col1'], Davidzon['col2'],[Davidzon['col4'],Davidzon['col3']],
                         fmt='o', markersize=5, ecolor='blue', color=obs_color[i_z])
        
    for i_z in range(0,len(ThisRedshiftList)):
        
        char_redshift="%0.2f" % ThisRedshiftList[i_z]
      
                            
        #MODEL
        #MR
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)
        
        G0_MR=G_MR[sel]   
    
        color_UV=G0_MR['MagDust'][:,0]-G0_MR['MagDust'][:,2]  
        color_VJ=G0_MR['MagDust'][:,2]-G0_MR['MagDust'][:,7]     
        G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (color_UV > minimum_y_color_cut[i_z+1]) & 
                    (color_UV > (color_VJ*slope_color_cut[i_z+1])+offset_color_cut[i_z+1])]
        #G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h))<-12.)]
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[i_z])    
        #StellarMass=np.log10((G0_MR['MassFromInSitu']+G0_MR['MassFromMergers']+G0_MR['MassFromBursts'])*1.e10*Hubble_h)
        
        bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
        hist_MR=np.histogram(StellarMass, bins=bin_arr, range=(xlim[0],xlim[1]))   
       
        #MRII
        if(MRII==1):
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, i_z, FullSnapshotList_MRII)
        
            G0_MRII=G_MRII[sel]        
            color_UV=G0_MRII['MagDust'][:,0]-G0_MRII['MagDust'][:,2]  
            color_VJ=G0_MRII['MagDust'][:,2]-G0_MRII['MagDust'][:,7]       
            G0_MRII=G0_MRII[(G0_MRII['StellarMass']>0.) & (color_UV > minimum_y_color_cut[i_z+1]) & 
                    (color_UV > (color_VJ*slope_color_cut[i_z+1])+offset_color_cut[i_z+1])]       
            #G0_MRII=G0_MRII[(G0_MRII['StellarMass']>0.) & (np.log10(G0_MRII['Sfr']/(G0_MRII['StellarMass']*1.e10/Hubble_h))<-12.)]
            StellarMass=stellar_mass_with_err(G0_MRII, Hubble_h, ThisRedshiftList[i_z])
        
            bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
            hist_MRII=np.histogram(StellarMass, bins=bin_arr, range=(xlim[0],xlim[1]))   
        
        #join MR+MRII & plot     
        if(MRII==1):
            cut_MR_MRII=10.
            (x_axis,y_axis)=plot_joint_MR_MRII(hist_MR, hist_MRII, cut_MR_MRII, Volume_MR, Volume_MRII, 
                                               bin, subplot, color=model_color[i_z],linewidth=2, linestyle='-')
        else:
            x_axis=hist_MR[1][0:len(hist_MR[1][:])-1]+bin/2.           
            hist_MR=hist_MR[0]       
            y_axis=np.log10(hist_MR/(Volume_MR*bin))
            subplot.plot(x_axis,y_axis, color=model_color[i_z], linewidth=2, linestyle='-') 
           
        
        #fa = open(Datadir+"SMF_Hen15_passive_z"+char_redshift+".txt", "w")
        #fa.write("%d\n" % len(x_axis))
        #for kk in range (0,len(x_axis)):               
        #    fa.write("%0.2f " % x_axis[kk] + "%0.2f\n" % y_axis[kk])         
        #fa.close()
        
        #fa = Datadir+"SMF_Hen15_passive_z"+char_redshift+".txt"
        #(x_axis,y_axis)=read_file(fa)    
        #subplot.plot(x_axis,y_axis, color='black', linewidth=2, linestyle='-')
        
        #LABELS
        if i_z==0:
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.91, 
                        color='black', xlog=0, ylog=0, label='Davidzon et al. (2017)', 
                        fontsize=13, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.935, 
                        color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15) 
                               
        #plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.35, 
        #            color='black', xlog=0, ylog=0, label='z='+char_redshift[:-1], 
         #           fontsize=14, fontweight='normal')      
        if i_z==2:
             plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.075, 
                         color='black', xlog=0, ylog=0, label='Passive', 
                         fontsize=16, fontweight='normal') 
           
    #endfor
 
    plt.tight_layout()
    plt.savefig('./fig/plots_stellar_mass_function_passive_allz_overplot.pdf')
    pdf.savefig()
    plt.close()
#endif stellar_mass_function_allz_overplot==1:
           
      

def stellar_mass_function_feedback_overplot(G_MR, Volume_MR, G_MRII, Volume_MRII, ThisRedshiftList, pdf):
           
    xlim=[8.5,13.5]
    ylim=[-5.5, 1.5]
        
    plot_color=['red','purple']        
   
    fig = plt.figure(figsize=(15,5))
    grid = gridspec.GridSpec(1, 4)
    grid.update(wspace=0.0, hspace=0.0)
         
    model_name=['SMF_Hen15_no_feedback_z','SMF_Hen15_only_AGN_z','SMF_Hen15_only_SN_z','SMF_Hen15_z']    
    label_name=['no feedback','only AGN','only SN','Hen15']    
    
    cmap = plt.get_cmap('winter')
    colors = [cmap(i) for i in np.linspace(0, 0.8, len(ThisRedshiftList))]
    
    for jj in range (0,4):    
        
        subplot=plt.subplot(grid[jj])
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
        subplot.yaxis.set_minor_locator(MultipleLocator(0.5)) 
        
        if jj==0:
            ylab='$\mathrm{log_{10}}(\phi [\mathrm{Mpc^{-3}} \mathrm{log_{10}}(M^{-1})])$'  
        else:
            ylab=''
        xlab='$\mathrm{log_{10}}(M_*[M_{\odot}])$'      
               
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14) 
        
        if jj>0:
            subplot.tick_params(axis='y', which='both', left='on', labelleft='off')
        
        #HALO MASS FUNCTION 
        bin=0.1
        for ii in range(0,len(ThisRedshiftList)):
            #MR
            (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)
        
            G0_MR=G_MR[sel]   
            G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (G0_MR['Type']==0)]   
            Mvir=np.log10(G0_MR['Mvir']*1.e10/Hubble_h*0.155)
            
            bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
            hist_MR=np.histogram(Mvir, bins=bin_arr, range=(xlim[0],xlim[1]))   
        
            #MRII
            if(MRII==1):
                (sel)=select_current_redshift(G_MRII, ThisRedshiftList, ii, FullSnapshotList_MRII)
        
                G0_MRII=G_MRII[sel]   
                G0_MRII=G0_MRII[(G0_MRII['StellarMass']>0.) & (G0_MRII['Type']==0)]
                Mvir=np.log10(G0_MRII['Mvir']*1.e10/Hubble_h*0.155)
        
                bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
                hist_MRII=np.histogram(Mvir, bins=bin_arr, range=(xlim[0],xlim[1]))   
           
        
            #join MR+MRII & plot        
            if(MRII==1):
                cut_MR_MRII=10.1
                (x_axis,y_axis)=plot_joint_MR_MRII(hist_MR, hist_MRII, cut_MR_MRII, Volume_MR, Volume_MRII, 
                                                   bin, subplot, color=colors[ii],linewidth=2, linestyle='--')
            else:
                x_axis=hist_MR[1][0:len(hist_MR[1][:])-1]+bin/2.           
                hist_MR=hist_MR[0]       
                y_axis=np.log10(hist_MR/(Volume_MR*bin))
                subplot.plot(x_axis,y_axis, color=colors[ii], linewidth=2, linestyle='--')         
                  
                
            #STELLAR MASS FUNCTION           
            '''#MR
            (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
            G0_MR=G_MR[sel]               
            StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
            
            bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
            hist_MR=np.histogram(StellarMass, bins=bin_arr, range=(xlim[0],xlim[1]))   
        
            #MRII
            if(MRII==1):
                (sel)=select_current_redshift(G_MRII, ThisRedshiftList, ii, FullSnapshotList_MRII)        
                G0_MRII=G_MRII[sel]                 
                StellarMass=stellar_mass_with_err(G0_MRII, Hubble_h, ThisRedshiftList[ii])
        
                bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
                hist_MRII=np.histogram(StellarMass, bins=bin_arr, range=(xlim[0],xlim[1]))   
           
        
            #join MR+MRII & plot        
            if(MRII==1):
                cut_MR_MRII=10.0
                (x_axis,y_axis)=plot_joint_MR_MRII(hist_MR, hist_MRII, cut_MR_MRII, Volume_MR, Volume_MRII, 
                                                   bin, subplot, color=colors[ii],linewidth=2, linestyle=':')
            else:
                x_axis=hist_MR[1][0:len(hist_MR[1][:])-1]+bin/2.           
                hist_MR=hist_MR[0]       
                y_axis=np.log10(hist_MR/(Volume_MR*bin))
                subplot.plot(x_axis,y_axis, color=colors[ii], linewidth=2, linestyle=':')'''             
                
                
                
                
                
            #STELLAR MASS FUNCTION  FROM FILES    
            char_redshift="%0.2f" %  ThisRedshiftList[ii]      
        
            fa = Datadir+model_name[jj]+char_redshift+".txt"
            (x_axis,y_axis)=read_file(fa)    
            subplot.plot(x_axis,y_axis, color=colors[ii], linewidth=2, linestyle='-')  
       
            if(jj<3):                   
                fa = Datadir+model_name[3]+char_redshift+".txt"
                (x_axis,y_axis)=read_file(fa)    
                subplot.plot(x_axis,y_axis, color=colors[ii], linewidth=2, linestyle=':') 
           
          
                
                
                
                
            #LABELS       
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.9, 
                        color='black', xlog=0, ylog=0, label=label_name[jj], fontsize=12, fontweight='normal') 
           
            if jj==0:
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.13, y_percentage=0.2, 
                            color='blue', xlog=0, ylog=0, label=r'$f_b \times M_{\mathrm{Halo}}$', 
                            fontsize=12, fontweight='normal') 
                plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, y_percentage=0.215, 
                            color='blue', x2_percentage=0.11, xlog=0, ylog=0, linestyle='--', linewidth=2)
                
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.13, y_percentage=0.135,
                            color='blue', xlog=0, ylog=0, label=r'$M_{*}$', 
                            fontsize=12, fontweight='normal') 
                plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, y_percentage=0.155, 
                            color='blue', x2_percentage=0.11, xlog=0, ylog=0, linestyle='-', linewidth=2)
                
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.13, y_percentage=0.07, 
                            color='blue', xlog=0, ylog=0, label='Hen15', 
                            fontsize=12, fontweight='normal') 
                plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, y_percentage=0.085, 
                            color='blue', x2_percentage=0.11, xlog=0, ylog=0, linestyle=':', linewidth=2)
    
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.43, 
                            color=colors[0], xlog=0, ylog=0, label='z=0',fontsize=12, fontweight='normal') 
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.38, 
                            color=colors[1], xlog=0, ylog=0, label='z=1',fontsize=12, fontweight='normal') 
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.33, 
                            color=colors[2], xlog=0, ylog=0, label='z=2',fontsize=12, fontweight='normal') 
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.28, 
                            color=colors[3], xlog=0, ylog=0, label='z=3',fontsize=12, fontweight='normal') 
        
        
    plt.tight_layout()
    plt.savefig('./fig/plots_stellar_mass_function_feedback_overplot.pdf')
    pdf.savefig()
    plt.close()

#end stellar_mass_function_feedback_overplot



def redfraction_color_cut_cuts(G_MR,G_MRII,ThisRedshiftList,pdf):
  
    fig = plt.figure(figsize=(20,4))
    

    for ii in range(0,len(ThisRedshiftList)):    
            
        char_redshift="%0.2f" % ThisRedshiftList[ii]
        
        if ii==0:
            gs1 = gridspec.GridSpec(1, 1)    
            gs1.update(left=0.065, right=0.25, top=0.9, bottom=0.2, wspace=0.0)
            
        if ii==1:          
            gs2 = gridspec.GridSpec(1, 4)  
            gs2.update(left=0.32, right=0.97, top=0.9, bottom=0.2, wspace=0.15, hspace=0.0)
                       
        if ii==0:            
            subplot=plt.subplot(gs1[ii])     
            xlim=[-23.,-13.]
            ylim=[0.7, 2.5]                
            xlab='$\mathrm{M}_r-5\,log_{10}(h)$' 
            ylab='u-r'          
            subplot.xaxis.set_major_locator(MultipleLocator(2))    
            subplot.xaxis.set_minor_locator(MultipleLocator(1))      
            subplot.yaxis.set_minor_locator(MultipleLocator(0.1))
            
        else:
            subplot=plt.subplot(gs2[ii-1])     
            xlim=[0.0, 2.5]
            ylim=[-0.5, 2.5]              
            xlab='V-J' 
            if (ii==1):
                ylab='U-V' 
                plt.tick_params(axis='y', which='both', left='on', labelleft='on')
            else:
                ylab='' 
                plt.tick_params(axis='y', which='both', left='on', labelleft='off')
            subplot.xaxis.set_major_locator(MultipleLocator(0.5))    
            subplot.xaxis.set_minor_locator(MultipleLocator(0.1))      
            subplot.yaxis.set_minor_locator(MultipleLocator(0.1))
    
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)
        
        #MODEL
        bin=0.25          
        Mass_arr=np.arange(xlim[0],xlim[1],bin)
        RedFraction=np.zeros(len(Mass_arr),dtype=np.float32)
        
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]   
        G0_MR=G0_MR[G0_MR['Type']==0]
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
        
        #z=0.
        if ThisRedshiftList[ii]==0:
            color_ur=G0_MR['MagDust'][:,15]-G0_MR['MagDust'][:,17]  
            Magr=G0_MR['MagDust'][:,17]-5.*np.log10(Hubble_h)
            
            Ngals=len(G0_MR)            
            bin_2d=[0.5,0.1]
            Nbins_2d=[int((xlim[1]-xlim[0])/bin_2d[0]),int((ylim[1]-ylim[0])/bin_2d[1])]
            H, xedges, yedges = np.histogram2d(Magr, color_ur, bins=Nbins_2d)            
            extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
            plt.subplots_adjust(bottom=0.15, left=0.15)                   
            mylevels = np.log10(np.logspace(-2.0, 1.0, num=15))  
            H = zoom(H, 20)        
            H=np.log10(H/np.amax(H))
            cont=plt.contourf(H.transpose()[::], origin='lower', cmap='Greys_r', levels=mylevels, extent=extent)        
          
                
            #Current CUT 
            x_arr=np.arange(xlim[0]-1., xlim[1]+1., 0.1)
            y_arr=offset_color_cut[ii]-slope_color_cut[ii]*np.tanh((x_arr+20.07)/1.09)
            subplot.plot(x_arr,y_arr, color='red', linestyle='-', linewidth=2) 
                
            #CUT in MCMC
            x_arr=np.arange(xlim[0]-1., xlim[1]+1., 0.1)
            y_arr=MCMC_offset_color_cut[ii]-MCMC_slope_color_cut[ii]*np.tanh((x_arr+20.07)/1.09)
            subplot.plot(x_arr,y_arr, color='red', linestyle='--', linewidth=2) 
            
            #ORIGINAL BALDRY CUT
            x_arr=np.arange(xlim[0]-1., xlim[1]+1., 0.1)
            y_arr=obs_slope_color_cut[ii]-obs_offset_color_cut[ii]*np.tanh((x_arr+20.07)/1.09)
            subplot.plot(x_arr,y_arr, color='royalblue', linestyle='-', linewidth=2) 
          
            #labels
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.18, y_percentage=0.15, 
                        color='black', xlog=0, ylog=0, label='Current cut', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim,
                x_percentage=0.04, y_percentage=0.17, color='red', x2_percentage=0.15, 
                xlog=0, ylog=0, linestyle='-', linewidth=2)  
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.18, y_percentage=0.1, 
                        color='black', xlog=0, ylog=0, label='OBS cut', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim,
                x_percentage=0.04, y_percentage=0.12, color='royalblue', x2_percentage=0.15, 
                xlog=0, ylog=0, linestyle='-', linewidth=2)  
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.18, y_percentage=0.05, 
                        color='black', xlog=0, ylog=0, label='MCMC', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim,
                x_percentage=0.04, y_percentage=0.07, color='red', x2_percentage=0.15, 
                xlog=0, ylog=0, linestyle='--', linewidth=2)                             
        #z>0.
        else:
            color_VJ=G0_MR['MagDust'][:,2]-G0_MR['MagDust'][:,7] 
            color_UV=G0_MR['MagDust'][:,0]-G0_MR['MagDust'][:,2]                         
           
            #sel=(color_UV < minimum_y_color_cut[ii]) | (color_UV < (color_VJ*slope_color_cut[ii])+offset_color_cut[ii])
           
            Ngals=len(G0_MR)            
            bin_2d=[0.1,0.1]
            Nbins_2d=[int((xlim[1]-xlim[0])/bin_2d[0]),int((ylim[1]-ylim[0])/bin_2d[1])]
            H, xedges, yedges = np.histogram2d(color_VJ[sel], color_UV[sel], bins=Nbins_2d)            
            extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
            plt.subplots_adjust(bottom=0.15, left=0.15) 
            mylevels = np.log10(np.logspace(-2.0, 1.0, num=15))  
            H = zoom(H, 20)        
            H=np.log10(H/np.amax(H))          
            cont=plt.contourf(H.transpose()[::], origin='lower', cmap='Greys_r', levels=mylevels, extent=extent)        
                      
            offset=offset_color_cut[ii]
            slope=slope_color_cut[ii]
            
            x_arr=np.arange(xlim[0]-1., xlim[1]+1., 0.01)             
            cut1= np.empty(len(x_arr))
            cut1.fill(1.3)  
            sel=(x_arr < (1.3-offset)/slope)
            subplot.plot(x_arr[sel],cut1[sel], color='red', linestyle='-', linewidth=2)             
            sel=(x_arr > (1.3-offset)/slope)
            cut2=slope*x_arr+offset
            subplot.plot(x_arr[sel],cut2[sel], color='red', linestyle='-', linewidth=2) 
            
            #MCMC CUT               
            offset=MCMC_offset_color_cut[ii]
            slope=MCMC_slope_color_cut[ii]
            
            x_arr=np.arange(xlim[0]-1., xlim[1]+1., 0.01)             
            cut1= np.empty(len(x_arr))
            cut1.fill(1.3)  
            sel=(x_arr < (1.3-offset)/slope)
            subplot.plot(x_arr[sel],cut1[sel], color='red', linestyle='--', linewidth=2)             
            sel=(x_arr > (1.3-offset)/slope)
            cut2=slope*x_arr+offset
            subplot.plot(x_arr[sel],cut2[sel], color='red', linestyle='--', linewidth=2)
            
            #ORIGINAL OBS CUT               
            offset=obs_offset_color_cut[ii]
            slope=obs_slope_color_cut[ii]
            
            x_arr=np.arange(xlim[0]-1., xlim[1]+1., 0.01)             
            cut1= np.empty(len(x_arr))
            cut1.fill(1.3)  
            sel=(x_arr < (1.3-offset)/slope)
            subplot.plot(x_arr[sel],cut1[sel], color='blue', linestyle='-', linewidth=2)             
            sel=(x_arr > (1.3-offset)/slope)
            cut2=slope*x_arr+offset
            subplot.plot(x_arr[sel],cut2[sel], color='blue', linestyle='-', linewidth=2)
            
        #LABELS          
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.9, 
                    color='black', xlog=0, ylog=0, label='z='+char_redshift[:-1], 
                    fontsize=14, fontweight='normal') 
           
    #endfor
        
    #plt.tight_layout()
    plt.savefig('./fig/plots_redfraction_color_cut_cuts.pdf')   
    pdf.savefig()
    plt.close()
#endif redfraction_color_cut_cuts


def redfraction_color_cut(G_MR,G_MRII,ThisRedshiftList,pdf):
           
    xlim=[7.5,11.5]
    ylim=[0., 1.2]   
    bin=[0.1,0.05]
    Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]

    fig = plt.figure(figsize=(one_five_size_large[0],one_five_size_large[1]))
    grid = gridspec.GridSpec(1, 5)
    grid.update(wspace=0.0, hspace=0.0)

    for ii in range(0,len(ThisRedshiftList)):
        
        char_redshift="%0.2f" % ThisRedshiftList[ii]
        
        subplot=plt.subplot(grid[ii])
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        
        xlab='$\mathrm{log_{10}}(M_*[M_{\odot}])$'      
        if ii==0:
            ylab='Red Fraction'
        else:
            ylab=''      
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)
        
        #subplot.text(xlim[0]+0.1,ylim[0]+.875,'z='+char_redshift, fontsize=16, fontweight='normal')
                        
    
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
        subplot.yaxis.set_minor_locator(MultipleLocator(0.1))
    
        if ii>0:
            plt.tick_params(axis='y', which='both', left='on', labelleft='off')
                   
        #OBSERVATIONS             
        file = MCMCdir + '/ObsConstraints/RedFraction_z'+char_redshift+'.txt'        
        f = open(file, 'r')     
        line = int(f.readline())     
        obs = Table.read(file, format='ascii', data_start=1, data_end=line+1)
        
        obs_xbin=obs['col1']+(obs['col2']-obs['col1'])/2.       
        subplot.errorbar(obs_xbin-2*np.log10(Hubble_h), obs['col3'],obs['col4'],
                 fmt='o', markersize=5, ecolor='blue', color='blue')
        #sub = plt.subplot(111)
        
        #PREVIOUS MODELS 
        '''RedshiftList_OldModels=[0.1,0.4,1.,2.,3.0]
        old_char_redshift="%0.2f" % RedshiftList_OldModels[ii]
        if do_previous_model1==1: 
            file = file_previous_model1+'_redfrac_colorcut_z'+old_char_redshift+'.txt' 
            model = Table.read(file, format='ascii')
            subplot.plot(model['col1']-2*np.log10(Hubble_h_WMAP7),model['col2'],color='red',
                         linestyle=linestyle_previous_model1, linewidth=2)           
        if do_previous_model2==1: 
            file = file_previous_model2+'_redfrac_colorcut_z'+char_redshift+'.txt'  
            model = Table.read(file, format='ascii')
            subplot.plot(model['col1']-2*np.log10(Hubble_h),model['col2'],color='red',
                         linestyle=linestyle_previous_model2, linewidth=2)'''
        
        #MODEL
        bin=0.25        
        Mass_arr=np.arange(xlim[0],xlim[1],bin)
        RedFraction=np.zeros(len(Mass_arr),dtype=np.float32)
        
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]   
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
        
        if(MRII==1):
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, ii, FullSnapshotList_MRII)        
            G0_MRII=G_MRII[sel]   
            StellarMass_MRII=stellar_mass_with_err(G0_MRII, Hubble_h, ThisRedshiftList[ii])
        
        #z=0.
        if ThisRedshiftList[ii]==0:
            color_ur=G0_MR['MagDust'][:,15]-G0_MR['MagDust'][:,17]  
            Magr=G0_MR['MagDust'][:,17]-5.*np.log10(Hubble_h)
            
            if(MRII==1):
                color_ur_MRII=G0_MRII['MagDust'][:,15]-G0_MRII['MagDust'][:,17]  
                Magr_MRII=G0_MRII['MagDust'][:,17]-5.*np.log10(Hubble_h)
            
            for ll in range(0,len(Mass_arr)):
                #sel_red=G0_MR[(color_ur>(offset_color_cut[ii]-slope_color_cut[ii]*np.tanh((Magr+18.07)/1.09))) &
                sel_red=G0_MR[(color_ur>(offset_color_cut[ii]-slope_color_cut[ii]*np.tanh((Magr+20.07)/1.09))) &
                    (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)]
                sel_all=G0_MR[(StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)]     
                if len(sel_all)>0.:
                    RedFraction[ll]=float(len(sel_red))/float(len(sel_all))  
                else:
                    RedFraction[ll]=0.
                    
                if((MRII==1) & (Mass_arr[ll]<9.5)):  
                    sel_red=G0_MRII[(color_ur_MRII>(offset_color_cut[ii]-slope_color_cut[ii] *
                                                    np.tanh((Magr_MRII+20.07)/1.09))) &
                                    (StellarMass_MRII>Mass_arr[ll]-bin/2.) & (StellarMass_MRII<Mass_arr[ll]+bin/2.)]
                    sel_all=G0_MRII[(StellarMass_MRII>Mass_arr[ll]-bin/2.) & (StellarMass_MRII<Mass_arr[ll]+bin/2.)]     
                    if len(sel_all)>0.:
                        RedFraction[ll]=float(len(sel_red))/float(len(sel_all))  
                    else:
                        RedFraction[ll]=0.
                        
        #z>0.
        else:
            color_UV=G0_MR['MagDust'][:,0]-G0_MR['MagDust'][:,2]  
            color_VJ=G0_MR['MagDust'][:,2]-G0_MR['MagDust'][:,7]       
                  
            for ll in range(0,len(Mass_arr)):               
                sel_red=G0_MR[(color_UV > minimum_y_color_cut[ii]) & 
                              (color_UV > (color_VJ*slope_color_cut[ii])+offset_color_cut[ii]) &
                              (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.) & (G0_MR['Type']==0)]           
                sel_all=G0_MR[(StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.) & (G0_MR['Type']==0)]
                if len(sel_all)>0.:
                    RedFraction[ll]=float(len(sel_red))/float(len(sel_all))
                else:
                    RedFraction[ll]=0.
                    
        subplot.plot(Mass_arr, RedFraction, color='red', linestyle='-', linewidth=2) 
      
    
        #MCMC sample
        if opt_plot_MCMC_sample==1:
            file = MCMCSampledir + 'mcmc_plus_obs_RedFraction_z'+char_redshift+'.txt' 
            if os.path.isfile(file):
                obs = Table.read(file, format='ascii')      
                subplot.plot(obs['col1']-2*np.log10(Hubble_h),obs['col4'], color='black', linewidth=2)
                #subplot.plot(obs['col1'],obs['col2'], color='black', linewidth=2)
            
                if ii==len(ThisRedshiftList)-1:
                    plot_label (subplot, 'label', xlim, ylim, 
                        x_percentage=0.55, y_percentage=0.85, color='black', xlog=0, ylog=0, 
                        label='MCMC sample', fontsize=13, fontweight='normal') 
                    plot_label (subplot, 'line', xlim, ylim,
                        x_percentage=0.44, y_percentage=0.87, color='black', x2_percentage=0.53, 
                        xlog=0, ylog=0, linestyle='-', linewidth=2)
    
    
        #LABELS    
        if ii==0:
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.91, 
                        color='black', xlog=0, ylog=0, label='Observations used in MCMC', 
                        fontsize=12, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.935, 
                        color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.03) 
        
            plot_label (subplot,'label',xlim,ylim,x_percentage=0.12,y_percentage=0.83, 
                        color='black',xlog=0,ylog=0,label=prefix_this_model, fontsize=12, fontweight='normal')             
            plot_label (subplot, 'line', xlim, ylim,x_percentage=0.02,y_percentage=0.855,
                        color='red', x2_percentage=0.1,xlog=0, ylog=0, linestyle='-', linewidth=2)
        
        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.65, 
                    color='black', xlog=0, ylog=0, label='z='+char_redshift[:-1], 
                    fontsize=14, fontweight='normal') 
            
            
            
        #CHANGE THIS*****************    
        #if ii==len(ThisRedshiftList)-1:
        #    plot_label_three_models (subplot, xlim, ylim, position='top_left')
                        
    #endfor
        
    plt.tight_layout()
    plt.savefig('./fig/plots_redfraction_color_cut.pdf')
    plt.savefig('./fig/HYW17_plots_redfraction_color_cut.pdf')  
    pdf.savefig()
    plt.close()
#endif redfraction_color_cut

  
def metals_vs_stellarmass(G_MR, ThisRedshiftList, pdf):
           
    xlim=[7.0,12.0]
    ylim=[-1.5, 1.0]   
    bin=0.2
        
    plot_color=['red','purple']        
   
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.25)) 
    
    xlab='$\mathrm{log_{10}}(M_*[M_{\odot}])$'       
    ylab='$\mathrm{log_{10}}(Z/Z_{\odot})$'         
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)    
        
    for ii in range (0,len(ThisRedshiftList)):
        
        char_redshift="%0.2f" % ThisRedshiftList[ii]
        
        #PREVIOUS MODELS    
        '''if ThisRedshiftList[ii]==0.1:
            char_redshift="%0.2f" % ThisRedshiftList[ii]
            if do_previous_model1==1: 
                file = file_previous_model1+'_metals_median_z'+char_redshift+'.txt' 
                model = Table.read(file, format='ascii')
                subplot.plot(model['col1'],model['col2'],color='red',linestyle=linestyle_previous_model1, linewidth=2)

            if do_previous_model2==1: 
                file = file_previous_model2+'_metals_z'+char_redshift+'.txt' 
                model = Table.read(file, format='ascii')
                subplot.plot(model['col1'],model['col2'],color='red',linestyle=linestyle_previous_model2, linewidth=2)'''

        if ii==0:        
            #observations from GALLAZI   
            Nbins=16
            obs_bin=0.2                 
            xmass=np.arange(8.91,12.11,obs_bin)#+2*np.log10(Hubble_h_WMAP7)      
            obsp50=np.array([-0.60,-0.61,-0.65,-0.61,-0.52,-0.41,-0.23,-0.11,-0.01,0.04,0.07,0.10,0.12,0.13,0.14,0.15])
            obsp16=np.array([-1.11,-1.07,-1.10,-1.03,-0.97,-0.90,-0.80,-0.65,-0.41,-0.24,-0.14,-0.09,-0.06,
                            -0.04,-0.03,-0.03])
            obsp84=np.array([-0.00,-0.00,-0.05,-0.01,0.05,0.09,0.14,0.17,0.20,0.22,0.24,0.25,0.26,0.28,0.29,0.30])         
            #subplot.errorbar(xmass, obsp50, yerr=[(obsp50-obsp16),(obsp84-obsp50)], color='blue', fmt='o') 
            subplot.plot(xmass,obsp50,color='blue',linestyle='-',linewidth=2)
            subplot.plot(xmass,obsp16,color='blue',linestyle='--',linewidth=2)
            subplot.plot(xmass,obsp84,color='blue',linestyle='--',linewidth=2)   
            
            #create OBSconstraint file for MCMC 
            '''xmass=np.arange(8.91,12.11,obs_bin)#+2*np.log10(Hubble_h_WMAP7)      
            obs_y=np.zeros(Nbins,dtype=np.float64)
            obs_y_err=np.zeros(Nbins,dtype=np.float64)
            for kk in range (0,len(xmass)):
                print('%0.2f %0.2f %0.2f %0.2f' % ((xmass[kk]-obs_bin/2.),(xmass[kk]+obs_bin/2),
                      (10.**obsp84[kk]+10.**obsp16[kk])/2.,(10.**obsp84[kk]-10.**obsp16[kk])/2.))                           
                obs_y[kk]=(10.**obsp84[kk]+10.**obsp16[kk])/2.
                obs_y_err[kk]=(10.**obsp84[kk]-10.**obsp16[kk])/2.              
            obs_y_err = [np.log10(obs_y/(obs_y-obs_y_err)),np.log10((obs_y+obs_y_err)/obs_y)]
            subplot.errorbar(xmass, np.log10(obs_y), yerr=obs_y_err, color='black', fmt='o')'''  
            
            #z=0
            label=prefix_this_model+', z=0'
            plot_label (subplot, 'label', xlim, ylim, 
                        x_percentage=0.12, y_percentage=0.9, color='black', xlog=0, ylog=0, 
                        label=label, fontsize=13, fontweight='normal')             
            plot_label (subplot,'line',xlim,ylim,x_percentage=0.02,y_percentage=0.92,
                        color=plot_color[0],x2_percentage=0.10,xlog=0,ylog=0,linestyle='-',linewidth=2)
            #z=3
            label=prefix_this_model+', z=3'
            plot_label (subplot, 'label', xlim, ylim, 
                        x_percentage=0.12, y_percentage=0.82, color='black', xlog=0, ylog=0, 
                        label=label, fontsize=13, fontweight='normal')             
            plot_label (subplot,'line',xlim,ylim,x_percentage=0.02,y_percentage=0.84,
                        color=plot_color[1],x2_percentage=0.10,xlog=0,ylog=0,linestyle='-',linewidth=2)
            
            #LABELS
            plot_label (subplot, 'label', xlim, ylim, 
                        x_percentage=0.12, y_percentage=0.74, color='black', xlog=0, ylog=0, 
                        label='Gallazzi 2005', fontsize=13, fontweight='normal')             
            plot_label (subplot,'line',xlim,ylim,x_percentage=0.02,y_percentage=0.76,
                        color='blue',x2_percentage=0.10,xlog=0,ylog=0,linestyle='-',linewidth=2)
            
                               
            #plot_label_three_models (subplot, xlim, ylim, position='bottom_left')
                     
    
        #Model
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]         
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])    
        if(opt_detailed_enrichment==0):                                 
            Metallicity=(G0_MR['MetalsStellarMass'])/((G0_MR['StellarMass'])*0.02)    
        else:
            MassInMetals=G0_MR['MetalsStellarMass'][:,0]+G0_MR['MetalsStellarMass'][:,1]+G0_MR['MetalsStellarMass'][:,2]      
            Metallicity=MassInMetals/(G0_MR['StellarMass']*0.02) 
            
        StellarMass=StellarMass[Metallicity>0.]    
        Metallicity=Metallicity[Metallicity>0.]       
        
        (x_binned, median, mean, pc16, pc84,rms)=median_and_percentiles (bin, xlim[0], xlim[1], StellarMass, Metallicity)    
        subplot.plot(x_binned, np.log10(median),color=plot_color[ii], linewidth=2)
        if (ii==0):
            subplot.plot(x_binned, np.log10(pc16),color=plot_color[ii], linewidth=2, linestyle='--')
            subplot.plot(x_binned, np.log10(pc84),color=plot_color[ii], linewidth=2, linestyle='--')
                   
       
        #MCMC sample
        if opt_plot_MCMC_sample==1:
            file = MCMCSampledir + 'mcmc_plus_obs_StellarMetallicityvsStellarMass_z'+char_redshift+'.txt' 
           
            if os.path.isfile(file):
                obs = Table.read(file, format='ascii')      
                subplot.plot(obs['col1']-2*np.log10(Hubble_h),np.log10(obs['col4']), color='black', linewidth=2)  
                if ii==len(ThisRedshiftList)-1:
                    plot_label (subplot, 'label', xlim, ylim, 
                        x_percentage=0.55, y_percentage=0.85, color='black', xlog=0, ylog=0, 
                        label='MCMC sample', fontsize=13, fontweight='normal') 
                    plot_label (subplot, 'line', xlim, ylim,
                        x_percentage=0.44, y_percentage=0.87, color='black', x2_percentage=0.53, 
                        xlog=0, ylog=0, linestyle='-', linewidth=2)    
    #endfor
        
    plt.tight_layout()
    plt.savefig('./fig/plots_metals_vs_stellarmass.pdf')
    plt.savefig('./fig/HYW17_plots_metals_vs_stellarmass.pdf')
    pdf.savefig()
    plt.close()

#end metals_vs_stellarmass


def gasmetals_vs_stellarmass(G_MR, G_MRII, RingRadius, RNUM, ThisRedshiftList,  pdf):
           
    #xlim=[7.5,12.0]
    #ylim=[8.0, 9.5] 
    xlim=[7.,12.0]
    ylim=[7., 10.] 
    bin=0.2
        
    plot_color=['red','purple']        
   
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.25)) 
    #subplot.yaxis.set_major_locator(MultipleLocator(0.5))    
    subplot.yaxis.set_major_locator(MultipleLocator(0.5))    
    subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
    
    xlab='$\mathrm{log_{10}}(M_*[M_{\odot}])$'       
    ylab='$12 + log_{10}$(O/H)$_{gas}$'         
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)    
        
    for ii in range (0,len(ThisRedshiftList)):
           
        char_redshift="%0.2f" % ThisRedshiftList[ii]    
            
        #Model
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0=G_MR[sel]        
        #G0=G0[(np.log10(G0['Sfr'])>-2.0) & (np.log10(G0['Sfr'])<-1.6)]        
        #G0=G0[(np.log10(G0['Sfr'])>-2.0)]
        '''The selection on type<2 is needed because these galaxies have their cooling artifically 
        shut down and therefore erroneously enriched by the SNIA and AGB channels without 
        cooling to balance it out'''
        G0=G0[(np.log10(G0['Sfr']/(G0['StellarMass']*1.e10/Hubble_h)) >-11.0) & (G0['ColdGas']>0.) &
                   (G0['Type']<2)]
      
        #G0_MR=np.random.choice(G0_MR, size=2000)     
        StellarMass=stellar_mass_with_err(G0, Hubble_h, ThisRedshiftList[ii])  
        #StellarMass=np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)
        if(opt_detailed_enrichment==0):                                 
            Metallicity=G0['MetalsColdGas']/G0['ColdGas']/0.0134
            StellarMass=StellarMass[Metallicity>0.]    
            Metallicity=np.log10(Metallicity[Metallicity>0.])+8.69
        else:
            #if(opt_individual_elements==0):
            if(opt_individual_elements<10):
                MassInGasMetals=G0['MetalsColdGas'][:,0]+G0['MetalsColdGas'][:,1]+G0['MetalsColdGas'][:,2]        
                Metallicity=MassInGasMetals/G0['ColdGas']/0.0134           
                StellarMass=StellarMass[Metallicity>0.]    
                Metallicity=np.log10(Metallicity[Metallicity>0.])+8.69           
                #Metallicity=np.log10(Metallicity[Metallicity>0.])+9.
                #(x_binned, median_MR, mean_MR, pc16, pc84,rms)=median_and_percentiles (bin, xlim[0], xlim[1], StellarMass1, Metallicity)         
                #subplot.plot(x_binned, median,color="brown", linewidth=2)
                #if(ii==0):
                #    plt.scatter(StellarMass1, Metallicity, s=8, color='red')
            
            #Mass * atomic weight
            else:
                N_H=G0['ColdGas_elements'][:,0]/1.
                N_O=G0['ColdGas_elements'][:,4]/16.            
                '''r_eff=GR['GasDiskRadius']/3.*1000./Hubble_h
                N_H=G0['ColdGasRings_elements'][:,0,0]/1.
                N_O=G0['ColdGasRings_elements'][:,0,4]/16.           
                for jj in range (1, RNUM):
                    sel=r_eff>RingRadius[jj]
                    N_H[sel]+=G0['ColdGasRings_elements'][sel,jj,0]/1.
                    N_O[sel]+=G0['ColdGasRings_elements'][sel,jj,4]/16.'''
                        
                sel=(N_O>0.) & (N_H>0.) 
                Metallicity=12.+np.log10(N_O[sel]/N_H[sel])            
                StellarMass=StellarMass[sel] 
            
       
        #Metallicity=np.log10(Metallicity)
        #plt.scatter(StellarMass, Metallicity, s=2, color='black')
          
        (x_binned_MR, median_MR, mean_MR, pc16_MR, pc84_MR,rms_MR)=median_and_percentiles (bin, xlim[0], xlim[1], 
                                                                               StellarMass, Metallicity)    
       
        if(MRII==1):
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, ii, FullSnapshotList_MRII)        
            G0=G_MRII[sel]    
            G0=G0[(np.log10(G0['Sfr']/(G0['StellarMass']*1.e10/Hubble_h)) >-11.0) & 
                  (G0['ColdGas']>0.) & (G0['Type']<2)]
       
            StellarMass=stellar_mass_with_err(G0, Hubble_h, ThisRedshiftList[ii])  
            if(opt_detailed_enrichment==0):                                 
                Metallicity=G0['MetalsColdGas']/G0['ColdGas']/0.0134
                StellarMass=StellarMass[Metallicity>0.]    
                Metallicity=np.log10(Metallicity[Metallicity>0.])+8.69
            else:          
                #Mass * atomic weight
                N_H=G0['ColdGas_elements'][:,0]/1.
                N_O=G0['ColdGas_elements'][:,4]/16.        
                sel=(N_O>0.) & (N_H>0.) 
                Metallicity=12.+np.log10(N_O[sel]/N_H[sel])            
                StellarMass=StellarMass[sel]    
                
            (x_binned_MRII,median_MRII,mean_MRII,pc16_MRII,pc84_MRII,rms_MRII)=median_and_percentiles(bin,xlim[0],xlim[1],
                                                                                       StellarMass,Metallicity) 
              
            cut_MR_MRII=9.0    
            sel=x_binned_MRII<cut_MR_MRII
            x_binned_MRII=x_binned_MRII[sel]
            median_MRII=median_MRII[sel] 
            pc16_MRII=pc16_MRII[sel]
            pc84_MRII=pc84_MRII[sel]
            sel=x_binned_MR>cut_MR_MRII
            x_binned_MR=x_binned_MR[sel]
            median_MR=median_MR[sel] 
            pc16_MR=pc16_MR[sel]
            pc84_MR=pc84_MR[sel]
            x_binned=np.concatenate((x_binned_MRII,x_binned_MR), axis=0)
            median=np.concatenate((median_MRII,median_MR), axis=0)
            pc16=np.concatenate((pc16_MRII,pc16_MR), axis=0)
            pc84=np.concatenate((pc84_MRII,pc84_MR), axis=0)
        else:
            x_binned=x_binned_MR
            median=median_MR
            pc16=pc16_MR
            pc84=pc84_MR
        
        sel=(pc16>ylim[0])
        subplot.plot(x_binned[sel], median[sel],color=plot_color[ii], linewidth=2)
        if (ii==0):
            subplot.plot(x_binned[sel], pc16[sel],color=plot_color[ii], linewidth=2, linestyle='--')
            subplot.plot(x_binned[sel], pc84[sel],color=plot_color[ii], linewidth=2, linestyle='--')
        
    
        '''#guo10    
        if (ii==0):
            x=[9.25,9.75,10.25,10.75,11.25,11.75,12.25]
            y=[8.71,8.89,9.07,9.15,9.12,9.24,9.40] 
            subplot.plot(x,y,color='blue', linewidth=2)'''
            
        #Observations
        obs_bin=0.5
        #Maiolino
        if(ii<2):    
        #if(ii==0):    
            log_mstar=np.arange(xlim[0], xlim[1],obs_bin)
            logM0=np.array([11.18,12.87])
            K0=np.array([9.04,8.90])
            gas_metallicity=-0.0864* ((log_mstar-logM0[ii])**2)+K0[ii]       
            subplot.plot(log_mstar, gas_metallicity,color=plot_color[ii], linewidth=2, linestyle=':')
        
        #create OBSconstraint file for MCMC 
        '''if (ii==0):
            log_mstar=np.arange(xlim[0], xlim[1],obs_bin) 
            obs_y=np.zeros(len(log_mstar),dtype=np.float64)
            obs_y_err=np.zeros(len(log_mstar),dtype=np.float64)
            for kk in range (0,len(log_mstar)):
                print('%0.2f %0.2f %0.2f %0.2f' % ((log_mstar[kk]+2*np.log10(Hubble_h_WMAP7)-obs_bin/2.)
                                                   ,(log_mstar[kk]+2*np.log10(Hubble_h_WMAP7)+obs_bin/2),
                                                   10**(gas_metallicity[kk]),10**(gas_metallicity[kk])*0.1 ))  
                      
                obs_y[kk]=10**(gas_metallicity[kk])
                obs_y_err[kk]=10**((gas_metallicity[kk])*0.1)            
            obs_y_err = [np.log10(obs_y/(obs_y-obs_y_err)),np.log10((obs_y+obs_y_err)/obs_y)]         
            subplot.errorbar(log_mstar, np.log10(obs_y), yerr=obs_y_err, color='black', fmt='o')'''
          
        
        #Tremonti2004
        if (ii==0):
            log_mstar=np.arange(xlim[0], xlim[1],obs_bin)
            log_gas_metallicity=-1.492+1.847*log_mstar-0.08026*(log_mstar**2)
            #subplot.plot(log_mstar, log_gas_metallicity,color='black', linewidth=2, linestyle=':')
            subplot.errorbar(log_mstar, log_gas_metallicity, yerr=0.1, color='blue', markeredgecolor='blue', fmt='o')
        
        #Andrews2012
        file = Datadir + '/brett_gas_metallicity2012.txt' 
        brett12 = Table.read(file, format='ascii')
        err=[-1.*brett12['err_down'],brett12['err_up']]
        subplot.errorbar(brett12['log_mass'], brett12['log_gas_metallicity'], yerr=err, color='black', fmt='o')
         
       
        
        
        
        
         #MCMC sample
        if opt_plot_MCMC_sample==1:
            file = MCMCSampledir + 'mcmc_plus_obs_ColdGasMetallicityvsStellarMass_z'+char_redshift+'.txt'            
            if os.path.isfile(file):
                obs = Table.read(file, format='ascii')                 
                subplot.plot(obs['col1']-2*np.log10(Hubble_h),np.log10(obs['col4']), color='black', linewidth=2)     
                if ii==len(ThisRedshiftList)-1:
                    plot_label (subplot, 'label', xlim, ylim, 
                        x_percentage=0.55, y_percentage=0.85, color='black', xlog=0, ylog=0, 
                        label='MCMC sample', fontsize=13, fontweight='normal') 
                    plot_label (subplot, 'line', xlim, ylim,
                        x_percentage=0.44, y_percentage=0.87, color='black', x2_percentage=0.53, 
                        xlog=0, ylog=0, linestyle='-', linewidth=2)    
        
        #LABELS        
        #z=0
        label=prefix_this_model+', z=0'
        plot_label (subplot, 'label', xlim, ylim, 
                    x_percentage=0.09, y_percentage=0.9, color='black', xlog=0, ylog=0, 
                    label=label, fontsize=13, fontweight='normal')             
        plot_label (subplot,'line',xlim,ylim,x_percentage=0.02,y_percentage=0.92,
                    color=plot_color[0],x2_percentage=0.07,xlog=0,ylog=0,linestyle='-',linewidth=2)
            
        '''#z=0      
        plot_label (subplot, 'label', xlim, ylim, 
                    x_percentage=0.12, y_percentage=0.9, color='black', xlog=0, ylog=0, 
                    label='Maiolino 2008 - z=0.07', fontsize=13, fontweight='normal')             
        plot_label (subplot,'line',xlim,ylim,x_percentage=0.02,y_percentage=0.92,
                    color=plot_color[0],x2_percentage=0.10,xlog=0,ylog=0,linestyle=':',linewidth=2)'''
        
          
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.82, 
                    color='black', xlog=0, ylog=0, label='Tremonti 2004', 
                    fontsize=13, fontweight='normal') 
        plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.84, 
                    color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.06) 
        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.74, 
                    color='black', xlog=0, ylog=0, label='Brett 2012', 
                    fontsize=13, fontweight='normal') 
        plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.76, 
                    color='black', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.06) 
        
        
        #z=3
        label=prefix_this_model+', z=3'
        plot_label (subplot, 'label', xlim, ylim, 
                    x_percentage=0.44, y_percentage=0.12, color='black', xlog=0, ylog=0, 
                    label=label, fontsize=13, fontweight='normal')             
        plot_label (subplot,'line',xlim,ylim,x_percentage=0.37,y_percentage=0.14,
                    color=plot_color[1],x2_percentage=0.42,xlog=0,ylog=0,linestyle='-',linewidth=2)
                
        plot_label (subplot, 'label', xlim, ylim, 
                    x_percentage=0.44, y_percentage=0.04, color='black', xlog=0, ylog=0, 
                    label='Maiolino 2008 - z=3.5', fontsize=13, fontweight='normal')             
        plot_label (subplot,'line',xlim,ylim,x_percentage=0.37,y_percentage=0.06,
                    color=plot_color[1],x2_percentage=0.42,xlog=0,ylog=0,linestyle=':',linewidth=2)
    #endfor
     
        
    '''#problematic objects      
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, 0, FullSnapshotList_MRII)        
    G0_MR=G_MR[sel]     
    G0_MR=G0_MR[(np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h)) >-11.0) & (G0_MR['ColdGas']>0.)]   
    StellarMass=np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)
    MassInGasMetals=G0_MR['MetalsColdGas'][:,0]+G0_MR['MetalsColdGas'][:,1]+G0_MR['MetalsColdGas'][:,2]     
    Metallicity=np.log10(MassInGasMetals/G0_MR['ColdGas']/0.0134)+8.69     
  
    sel=((StellarMass>7.5) & (StellarMass<8.5) & (Metallicity>8.5))
    G_sel=G0_MR[sel]
    print(len(G_sel))    
    print(len(G0_MR))
    
    print(StellarMass[sel])
    print(Metallicity[sel])
    
    MassInGasMetals=G_sel['MetalsColdGas'][:,1]
    Metallicity=np.log10(MassInGasMetals/G_sel['ColdGas']/0.0134)+8.69  
    print(G_sel['Type'])'''
    
    
    plt.tight_layout()
    plt.savefig('./fig/plots_gas_metals_vs_stellarmass.pdf')
    plt.savefig('./fig/HYW17_plots_gas_metals_vs_stellarmass.pdf')
    pdf.savefig()
    plt.close()

#end gasmetals_vs_stellarmass 

def BHBM(G_MR, ThisRedshiftList, pdf):
   
    for ii in range(0,len(ThisRedshiftList)):        
        
        xlim=[8.5,12.5]
        ylim=[5.0, 10.5]
        bin=[0.25,0.25]
        Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]
    
        plot_color=['red','purple']        
     
        fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
        subplot=plt.subplot()
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
            
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25)) 
            
        xlab='$\mathrm{log_{10}}(M_{\mathrm{Bulge}}[M_{\odot}])$'       
        ylab='$\mathrm{log_{10}}(M_{\mathrm{BH}}[M_{\odot}])$' 
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   
            
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR_unsel=G_MR[sel]   
        G0_MR=G0_MR_unsel[(G0_MR_unsel['BulgeMass'] > 0.) & (G0_MR_unsel['BlackHoleMass'] > 0.)]
        Ngals=len(G0_MR) 
       
        BulgeMass=(np.log10(G0_MR['BulgeMass']*1.e10/Hubble_h)) 
        BHMass=(np.log10(G0_MR['BlackHoleMass']*1.e10/Hubble_h)) 
               
        #plt.scatter(BulgeMass, BHMass, s=5, color='black')          
        H, xedges, yedges = np.histogram2d(BulgeMass, BHMass, bins=Nbins,
                                           range=[[xlim[0]-bin[0]/2., xlim[1]+bin[0]/2.], 
                                                  [ylim[0]-bin[1]/2., ylim[1]+bin[1]/2.]])            
        extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
        plt.subplots_adjust(bottom=0.15, left=0.15)     
        mylevels = np.log10(np.logspace(-2.0, 0.0, num=10))  
        H = zoom(H, 20)        
        H=np.log10(H/np.amax(H))       
        cont=plt.contourf(H.transpose()[::], origin='lower', cmap='Greys_r', levels=mylevels, extent=extent)  
        #cont=plt.contourf(H.transpose()[::], origin='lower', cmap='rainbow', levels=mylevels, extent=extent)     
        plt.colorbar(format='%0.1f') 
        #(ax, cmap=None, norm=None, alpha=None, values=None, boundaries=None, orientation='vertical', 
        #ticklocation='auto', extend='neither', spacing='uniform', ticks=None, format=None, 
        #drawedges=False, filled=True, extendfrac=None, extendrect=False, label='')
       
        file = Datadir + 'mcconnel2012.dat'
        obs = Table.read(file, format='ascii', data_start=20)     
        obs_x = np.log10(obs['col14'])
        obs_y = np.log10(obs['col3'])        
        obs_x_err=np.zeros(len(obs_x),dtype=np.float64)+0.24 
        obs_y_err = [np.log10(obs['col3']/obs['col4']),np.log10(obs['col5']/obs['col3'])]
       
        subplot.errorbar(obs_x, obs_y,xerr= obs_x_err, yerr=obs_y_err,
                 fmt='o', markersize=5, ecolor='blue', color='blue')
       
        '''x_arr=np.arange(xlim[0],xlim[1],0.01)
        y_arr=-0.952381*x_arr+17.2        
        subplot.plot(x_arr,y_arr)
        y_arr=-0.952381*x_arr+18.88
        subplot.plot(x_arr,y_arr)
        y_arr=1.05*x_arr-2.91961
        subplot.plot(x_arr,y_arr)'''
        
        #LABELS        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.9, 
                    color='black', xlog=0, ylog=0, label='McConnel 2012', 
                    fontsize=13, fontweight='normal') 
        plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.925, 
                    color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15) 
        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.7, y_percentage=0.1, 
                    color='black', xlog=0, ylog=0, label='BHBM', 
                    fontsize=13, fontweight='normal') 
        
    plt.tight_layout()
    plt.savefig('./fig/plots_bhbm.pdf')
    plt.savefig('./fig/HYW17_plots_bhbm.pdf')
    pdf.savefig()
    plt.close()

#end BHBM




def SFRF(G_MR, Volume_MR, G_MRII, Volume_MRII, ThisRedshiftList, pdf):
       
    fig = plt.figure(figsize=(two_two_size_small[0],two_two_size_small[1]))
    grid = gridspec.GridSpec(2, 2)
    grid.update(wspace=0.0, hspace=0.0)

    xlim=[-1.5,4.0]
    ylim=[-7.0,0.5]
    bin=0.2
        
    for i_z in range(0,len(ThisRedshiftList)):
        
        char_redshift="%0.2f" % ThisRedshiftList[i_z]
        
        subplot=plt.subplot(grid[i_z])

        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        if i_z==2 or i_z == 3: 
            xlab='$\mathrm{log_{10}}(\mathrm{SFR}[M_{\odot}yr^{-1}])$' 
            subplot.set_xlabel(xlab, fontsize=14)
       
        if i_z==0 or i_z == 2:
            ylab='$\mathrm{log_{10}}(\phi [\mathrm{Mpc^{-3}} \mathrm{log_{10}}(SFR^{-1})])$'
    
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
        subplot.yaxis.set_minor_locator(MultipleLocator(0.25))
    
        if i_z==1 or i_z==3:
            plt.tick_params(axis='y', which='both', left='on', labelleft='off')
        if i_z==0 or i_z==1:    
            plt.tick_params(axis='x', which='both', bottom='on', labelbottom='off')

        #MR
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)         
        G0_MR=G_MR[sel]          
        G0_MR=G0_MR[G0_MR['Sfr']>0.]
        #SFR=(np.log10(G0_MR['Sfr']))
        SFR=sfr_with_err(G0_MR, Hubble_h, ThisRedshiftList[i_z])
                
        bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
        hist_MR=np.histogram(SFR, bins=bin_arr, range=(xlim[0],xlim[1]))   
         
        #MRII    
        if(MRII==1):
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, i_z, FullSnapshotList_MRII)         
            G0_MRII=G_MRII[sel]   
            G0_MRII=G0_MRII[G0_MRII['Sfr']>0.]
            #SFR=(np.log10(G0_MRII['Sfr']))
            SFR=sfr_with_err(G0_MRII, Hubble_h, ThisRedshiftList[i_z])
                    
            bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
            hist_MRII=np.histogram(SFR, bins=bin_arr, range=(xlim[0],xlim[1]))   
                        
        #join MR+MRII & plot     
        if(MRII==1):
            cut_MR_MRII=1.0
            (x_axis,y_axis)=plot_joint_MR_MRII(hist_MR, hist_MRII, cut_MR_MRII, Volume_MR, Volume_MRII, 
                                               bin, subplot, color='red',linewidth=2, linestyle='-')
        else:
            x_axis=hist_MR[1][0:len(hist_MR[1][:])-1]+bin/2.           
            hist_MR=hist_MR[0]       
            y_axis=np.log10(hist_MR/(Volume_MR*bin))
            subplot.plot(x_axis,y_axis, color='red', linewidth=2, linestyle='-')           
                
                
        #OBSERVATIONS
        plot_sfrf_obs(i_z,subplot,xlim,ylim)    
    
            
        #MCMC sample
        if i_z==0:
            if opt_plot_MCMC_sample==1:
                file = MCMCSampledir + 'mcmc_plus_obs_SFRF_z'+char_redshift+'.txt' 
                if os.path.isfile(file):
                    obs = Table.read(file, format='ascii')      
                    subplot.plot(obs['col1']-2.*np.log10(Hubble_h_WMAP7),
                                 np.log10(obs['col4'])+3.*np.log10(Hubble_h_WMAP7), color='black', linewidth=2)
            
                    if i_z==len(ThisRedshiftList)-1:
                        plot_label(subplot, 'label', xlim, ylim, x_percentage=0.55, y_percentage=0.85, 
                                    color='black', xlog=0, ylog=0, label='MCMC sample', fontsize=13, fontweight='normal') 
                        plot_label(subplot, 'line', xlim, ylim, x_percentage=0.44, y_percentage=0.87, 
                                   color='black', x2_percentage=0.53, xlog=0, ylog=0, linestyle='-', linewidth=2)    
            
            
            
        #LABELS   
        if(i_z==0):
            plot_label (subplot,'label',xlim,ylim,x_percentage=0.15,y_percentage=0.9, 
                        color='black',xlog=0,ylog=0,label=prefix_this_model,fontsize=13,fontweight='normal') 
            plot_label (subplot,'line',xlim,ylim,x_percentage=0.04, y_percentage=0.92,
                        color='red',x2_percentage=0.13,xlog=0,ylog=0,linestyle='-',linewidth=2)
   
            
    plt.tight_layout()
    plt.savefig('./fig/plots_sfrf.pdf')
    plt.savefig('./fig/HYW17_plots_sfrf.pdf')
    pdf.savefig()
    plt.close()

#end SFRF


def plot_sfrf_obs(i_z,subplot,xlim,ylim):
        
        if i_z==0:
            file = Datadir + '/sfrf_gruppioni2015_z0.00_z0.30.txt'
            obs = Table.read(file, format='ascii')      
            obs_x = (obs['log_sfr_low']+obs['log_sfr_high'])/2.            
            subplot.errorbar(obs_x, obs['log_phi'], yerr=obs['err'], fmt='o', markersize=5, 
                             markeredgecolor='blue', color='blue')
           
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.05, 
                        color='black', xlog=0, ylog=0, label='Gruppioni 2015-0.0<z<0.3', 
                        fontsize=11, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.065, 
                        color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.12) 
           
        if i_z==1:
            file = Datadir + '/sfrf_gruppioni2015_z0.80_z1.00.txt'
            obs = Table.read(file, format='ascii')      
            obs_x = (obs['log_sfr_low']+obs['log_sfr_high'])/2.            
            subplot.errorbar(obs_x, obs['log_phi'], yerr=obs['err'], fmt='o', markersize=5, 
                             color='blue', markeredgecolor='blue') 
            
            file = Datadir + '/sfrf_gruppioni2015_z1.00_z1.20.txt'
            obs = Table.read(file, format='ascii')      
            obs_x = (obs['log_sfr_low']+obs['log_sfr_high'])/2.            
            subplot.errorbar(obs_x, obs['log_phi'], yerr=obs['err'], fmt='o', markersize=5, 
                             mfc='white', markeredgecolor='blue', color='blue') 
                             
            file = Datadir + '/sfrf_Gruppioni2013_z0.70_1.00.txt'
            obs = Table.read(file, format='ascii')            
            subplot.errorbar(obs['log_sfr'], obs['phi'], yerr=[obs['err_down'],obs['err_up']], 
                             fmt='s', markersize=5, color='limegreen', markeredgecolor='limegreen') 
            
            file = Datadir + '/sfrf_Gruppioni2013_z1.00_1.20.txt'
            obs = Table.read(file, format='ascii')                     
            subplot.errorbar(obs['log_sfr'], obs['phi'], yerr=obs['err'], fmt='s', markersize=5, 
                             mfc='white', markeredgecolor='limegreen', color='limegreen') 
            
            file = Datadir + '/sfrf_Magnelli2011_z0.80_1.00.txt'
            obs = Table.read(file, format='ascii')          
            subplot.errorbar(obs['log_sfr'], obs['phi'], yerr=obs['err'], fmt='^', markersize=5, 
                             color='purple', markeredgecolor='purple') 
            
            file = Datadir + '/sfrf_Magnelli2011_z1.00_1.30.txt'
            obs = Table.read(file, format='ascii')                     
            subplot.errorbar(obs['log_sfr'], obs['phi'], yerr=[obs['err_down'],obs['err_up']], fmt='^', markersize=5, 
                             mfc='white', markeredgecolor='purple', color='purple') 
            
            file = Datadir + '/sfrf_Ly2011_z0.84.txt'
            obs = Table.read(file, format='ascii')            
            subplot.errorbar(obs['log_sfr'], obs['phi'], yerr=obs['err'], fmt='*', markersize=5, 
                             color='orange', markeredgecolor='orange')
            
            '''file = Datadir + '/sfrf_Sobral2013_z0.84.txt'
            obs = Table.read(file, format='ascii')            
            subplot.errorbar(obs['log_sfr'], obs['phi'], yerr=obs['err'], fmt='o', markersize=5, 
                             color='pink', markeredgecolor='pink')'''
                                    
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.35, 
                        color='black', xlog=0, ylog=0, label='Ly 2011-z=0.84', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.37, 
                        color='orange', xlog=0, ylog=0, sym='*', sym_size=5, err_size=0.12)
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.30, 
                        color='black', xlog=0, ylog=0, label='Magnelli 2011-0.7<z<1.0', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.315, 
                        color='purple', xlog=0, ylog=0, sym='^', sym_size=5, err_size=0.12)            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.25, 
                        color='black', xlog=0, ylog=0, label='Magnelli 2011-1.0<z<1.2', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.27, 
                        color='purple', xlog=0, ylog=0, sym='^', sym_size=5, err_size=0.12, mfc='white')
                          
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.20, 
                        color='black', xlog=0, ylog=0, label='Gruppioni 2015-0.8<z<1.0', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.22, 
                        color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.12)            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.15, 
                        color='black', xlog=0, ylog=0, label='Gruppioni 2015-1.0<z<1.2', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.17, 
                        color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.12, mfc='white')
                
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.10, 
                        color='black', xlog=0, ylog=0, label='Gruppioni 2013-0.7<z<1.0', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.115, 
                        color='limegreen', xlog=0, ylog=0, sym='s', sym_size=5, err_size=0.12)            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.05, 
                        color='black', xlog=0, ylog=0, label='Gruppioni 2013-1.0<z<1.2', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.065, 
                        color='limegreen', xlog=0, ylog=0, sym='s', sym_size=5, err_size=0.12, mfc='white')    
                
        if i_z==2:
            file = Datadir + '/sfrf_gruppioni2015_z1.70_z2.00.txt'
            obs = Table.read(file, format='ascii')      
            obs_x = (obs['log_sfr_low']+obs['log_sfr_high'])/2.            
            subplot.errorbar(obs_x, obs['log_phi'], yerr=obs['err'], fmt='o', markersize=5, 
                             markeredgecolor='blue', color='blue')
            
            file = Datadir + '/sfrf_gruppioni2015_z2.00_z2.50.txt'
            obs = Table.read(file, format='ascii')      
            obs_x = (obs['log_sfr_low']+obs['log_sfr_high'])/2.            
            subplot.errorbar(obs_x, obs['log_phi'], yerr=obs['err'], fmt='o', markersize=5, 
                             mfc='white', markeredgecolor='blue', color='blue') 
            
            file = Datadir + '/sfrf_Magnelli2011_z2.30.txt'
            obs = Table.read(file, format='ascii')                 
            subplot.errorbar(obs['log_sfr'], obs['phi'], yerr=[obs['err_down'],obs['err_up']], fmt='^', markersize=5, 
                             color='purple', markeredgecolor='purple')
            
            '''file = Datadir + '/sfrf_Parsa2015_z2.00.txt'
            obs = Table.read(file, format='ascii')                 
            subplot.errorbar(obs['log_sfr'], obs['phi'], yerr=obs['err'], fmt='o', markersize=5, 
                             color='orange', markeredgecolor='orange')'''
            
            file = Datadir + '/sfrf_Reddy2008_z1.80_2.30.txt'
            obs = Table.read(file, format='ascii')                 
            subplot.errorbar(obs['log_sfr'], obs['phi'], yerr=obs['err'], fmt='D', markersize=5, 
                             color='magenta', markeredgecolor='magenta')
            
            '''file = Datadir + '/sfrf_Sobral2013_z2.23.txt'
            obs = Table.read(file, format='ascii')            
            subplot.errorbar(obs['log_sfr'], obs['phi'], yerr=obs['err'], fmt='o', markersize=5, 
                             color='pink', markeredgecolor='pink')'''           
            
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.20, 
                        color='black', xlog=0, ylog=0, label='Gruppioni 2015-1.7<z<2.0', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.22, 
                        color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.12)            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.15, 
                        color='black', xlog=0, ylog=0, label='Gruppioni 2015-2.0<z<2.5', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.17, 
                        color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.12, mfc='white')
                
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.10, 
                        color='black', xlog=0, ylog=0, label='Magnelli 2011-z=2.3', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.115, 
                        color='purple', xlog=0, ylog=0, sym='^', sym_size=5, err_size=0.12) 
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.05, 
                        color='black', xlog=0, ylog=0, label='Reddy 2008-1.8<z<2.3', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.065, 
                        color='magenta', xlog=0, ylog=0, sym='D', sym_size=5, err_size=0.12, mfc='white')    
            
        if i_z==3:
            file = Datadir + '/sfrf_gruppioni2015_z2.50_z3.00.txt'
            obs = Table.read(file, format='ascii')      
            obs_x = (obs['log_sfr_low']+obs['log_sfr_high'])/2.            
            subplot.errorbar(obs_x, obs['log_phi'], yerr=obs['err'], fmt='o', markersize=5, 
                             markeredgecolor='blue', color='blue')
            
            file = Datadir + '/sfrf_gruppioni2015_z3.00_z4.20.txt'
            obs = Table.read(file, format='ascii')      
            obs_x = (obs['log_sfr_low']+obs['log_sfr_high'])/2.            
            subplot.errorbar(obs_x, obs['log_phi'], yerr=obs['err'], fmt='o', markersize=5, 
                             mfc='white', markeredgecolor='blue', color='blue')     
            
            file = Datadir + '/sfrf_Gruppioni2013_z2.50_3.00.txt'
            obs = Table.read(file, format='ascii')            
            subplot.errorbar(obs['log_sfr'], obs['phi'], yerr=obs['err'], fmt='s', markersize=5, 
                             color='limegreen', markeredgecolor='limegreen') 
            
            file = Datadir + '/sfrf_Reddy2008_z3.10.txt'
            obs = Table.read(file, format='ascii')                 
            subplot.errorbar(obs['log_sfr'], obs['phi'], yerr=obs['err'], fmt='D', markersize=5, 
                             color='magenta', markeredgecolor='magenta')
            
            '''file = Datadir + '/sfrf_Parsa2015_z3.00.txt'
            obs = Table.read(file, format='ascii')                 
            subplot.errorbar(obs['log_sfr'], obs['phi'], yerr=obs['err'], fmt='o', markersize=5, 
                             color='orange', markeredgecolor='orange')'''
            
            file = Datadir + '/sfrf_vanderBurg2010_z3.10.txt'
            obs = Table.read(file, format='ascii')            
            subplot.errorbar(obs['log_sfr'], obs['phi'], yerr=obs['err'], fmt='v', markersize=5, 
                             color='dimgrey', markeredgecolor='dimgrey') 
            
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.25, 
                        color='black', xlog=0, ylog=0, label='Reddy 2008-z=3.1', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.27, 
                        color='magenta', xlog=0, ylog=0, sym='D', sym_size=5, err_size=0.12)
                          
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.20, 
                        color='black', xlog=0, ylog=0, label='van der Burg 2010-z=3.1', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.22, 
                        color='dimgrey', xlog=0, ylog=0, sym='v', sym_size=5, err_size=0.12)       
                
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.15, 
                        color='black', xlog=0, ylog=0, label='Gruppioni 2015-0.8<z<1.0', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.17, 
                        color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.12)            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.10, 
                        color='black', xlog=0, ylog=0, label='Gruppioni 2015-1.0<z<1.2', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.12, 
                        color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.12, mfc='white')
          
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.05, 
                        color='black', xlog=0, ylog=0, label='Gruppioni 2013-2.5<z<3.0', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.065, 
                        color='limegreen', xlog=0, ylog=0, sym='s', sym_size=5, err_size=0.12) 
            
            

def gas_fraction(G_MR, ThisRedshiftList, pdf):
   
    for ii in range(0,len(ThisRedshiftList)):        
        
        char_redshift="%0.2f" % ThisRedshiftList[ii]
        
        xlim=[8.5,12.0]
        ylim=[-2.,1.0]
       
        bin=0.1
        plot_color=['red','purple']        
       
        fig = plt.figure(figsize=(5,4))
        subplot=plt.subplot()
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
            
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))
            
        xlab='$log_{10}(M_*[M_{\odot}])$'
        ylab='$M_{\mathrm{Cold}}/M_*$'     
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   
            
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]          
        G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (G0_MR['ColdGas']>0.)]
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
        #Fraction=np.log10(G0_MR['ColdGas']*1.e10*Hubble_h)-StellarMass 
        Fraction=G0_MR['ColdGas']/G0_MR['StellarMass']
   
        (x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles (bin, xlim[0], xlim[1], StellarMass, Fraction)          
        sel=(median!=0)        
        subplot.plot(x_binned[sel], np.log10(median[sel]),color=plot_color[ii], linewidth=2)     
        subplot.plot(x_binned[sel], np.log10(pc16[sel]),color=plot_color[ii], linewidth=2, linestyle='--')
        subplot.plot(x_binned[sel], np.log10(pc84[sel]),color=plot_color[ii], linewidth=2, linestyle='--')
    
    
        file = Datadir + 'peeples_2015.txt'
        obs = Table.read(file, format='ascii', data_start=0)       
        obs_x = (obs['col1']+obs['col2'])/2.-2*np.log10(Hubble_h_WMAP7)
        obs_y = np.log10(obs['col3'])       
        obs_y_err = [np.log10(obs['col3']/(obs['col3']-obs['col4'])),np.log10((obs['col3']+obs['col4'])/obs['col3'])]
               
        subplot.errorbar(obs_x, obs_y, yerr=obs_y_err,
                 fmt='o', markersize=5, ecolor='blue', color='blue')
              
        #MCMC sample
        if opt_plot_MCMC_sample==1:
            file = MCMCSampledir + 'mcmc_plus_obs_ColdGasFractionvsStellarMass_z'+char_redshift+'.txt' 
            if os.path.isfile(file):
                obs = Table.read(file, format='ascii')      
                subplot.plot(obs['col1']-2*np.log10(Hubble_h),np.log10(obs['col4']), color='black', linewidth=2)
                #subplot.errorbar(obs['col1'],obs['col2'], yerr=obs['col3'], 
                #                 fmt='o', markersize=5, ecolor='black', color='black')
            
                if ii==len(ThisRedshiftList)-1:
                    plot_label (subplot, 'label', xlim, ylim, 
                        x_percentage=0.55, y_percentage=0.85, color='black', xlog=0, ylog=0, 
                        label='MCMC sample', fontsize=13, fontweight='normal') 
                    plot_label (subplot, 'line', xlim, ylim,
                        x_percentage=0.44, y_percentage=0.87, color='black', x2_percentage=0.53, 
                        xlog=0, ylog=0, linestyle='-', linewidth=2)
            
        #LABELS        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.2, 
                    color='black', xlog=0, ylog=0, label='Peeples 2015', 
                    fontsize=13, fontweight='normal') 
        plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.225, 
                    color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15) 
           
        plot_label (subplot, 'label', xlim, ylim, 
                x_percentage=0.15, y_percentage=0.1, color='black', xlog=0, ylog=0, 
                label=prefix_this_model, fontsize=13, fontweight='normal') 
        plot_label (subplot, 'line', xlim, ylim,
                x_percentage=0.04, y_percentage=0.12, color='red', x2_percentage=0.13, 
                xlog=0, ylog=0, linestyle='-', linewidth=2)
        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.55, y_percentage=0.8, 
                    color='black', xlog=0, ylog=0, label='Gas Fraction', 
                    fontsize=13, fontweight='normal')   
            
    plt.tight_layout()
    plt.savefig('./fig/plots_gas_fraction.pdf')  
    pdf.savefig()
    plt.close()

#end gas fraction


def HI_over_Lr_vs_HI(G_MR, ThisRedshiftList, pdf):
   
    for ii in range(0,len(ThisRedshiftList)):        
        
        char_redshift="%0.2f" % ThisRedshiftList[ii]
        
        xlim=[7.5,11.0]
        ylim=[-1.,1.0]
        ylim=[-2.5,2.0]       
        bin=0.25
      
        fig = plt.figure(figsize=(5,4))
        subplot=plt.subplot()
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
            
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))
            
        xlab='$\mathrm{log_{10}}(M_{\mathrm{HI}}[M_{\odot}])$'
        ylab='$\mathrm{log_{10}}(M_{\mathrm{HI}}/L_{\mathrm{r}})$'     
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   
            
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]          
        G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (G0_MR['ColdGas']>0.)]        
        if(opt_rings==1):
            HI=(np.log10(G0_MR['ColdGas']*(1.-G0_MR['H2fraction'])*1.e10/Hubble_h))              
        else:            
            HI=(np.log10(G0_MR['ColdGas']*0.54*1.e10/Hubble_h))
           
        Lr=mag_to_lum(G0_MR['MagDust'][:,17])      
        #Mass=np.log10(G0_MR['StellarMass']*1.e10*Hubble_h)
        Fraction=np.log10((10**HI)/Lr)
        #Fraction=np.log10(10**HI/10**Mass)

        (x_binned, median, mean, pc16, pc84,rms)=median_and_percentiles (bin, xlim[0], xlim[1], HI, Fraction)
        #(x_binned, median, mean, pc16, pc84)=median_and_percentiles (bin, xlim[0], xlim[1], Mass, Fraction)
        sel=(median!=0)        
        subplot.plot(x_binned[sel], median[sel],color='red', linewidth=2)     
        subplot.plot(x_binned[sel], pc16[sel],color='red', linewidth=2, linestyle='--')
        subplot.plot(x_binned[sel], pc84[sel],color='red', linewidth=2, linestyle='--')
    
    
        #OBSERVATIONS
        file=Datadir+'Haynes_newfile'
        fits_table=fits.open(file)
        haynes = fits_table[1]
                  
        haynes_MHI=haynes.data['HI']
        haynes_Magr=haynes.data['mr']-5*(np.log10(haynes.data['distance']*1.e6)-1)
        haynes_Lr=mag_to_lum(haynes_Magr)
        fraction=np.log10(10**haynes_MHI/haynes_Lr)
        (x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles(bin, xlim[0], xlim[1]-0.25, haynes_MHI, fraction) 
        obs_x = x_binned
        obs_y = (pc84+pc16)/2       
        obs_y_err = (pc84-pc16)/2               
        subplot.errorbar(obs_x, obs_y, yerr=obs_y_err, fmt='o', markersize=5, ecolor='blue', color='blue')
                
        fits_table.close()
            
        #MCMC sample
        if opt_plot_MCMC_sample==1:
            file = MCMCSampledir + 'mcmc_plus_obs_HI_over_Lr_vs_HI_z'+char_redshift+'.txt' 
            if os.path.isfile(file):
                obs = Table.read(file, format='ascii')      
                subplot.plot(obs['col1']-2.*np.log10(Hubble_h),np.log10(obs['col4']), color='black', linewidth=2)
                #subplot.errorbar(obs['col1'],obs['col2'], yerr=obs['col3'], 
                #                 fmt='o', markersize=5, ecolor='black', color='black')
            
                if ii==len(ThisRedshiftList)-1:
                    plot_label (subplot, 'label', xlim, ylim, 
                        x_percentage=0.55, y_percentage=0.85, color='black', xlog=0, ylog=0, 
                        label='MCMC sample', fontsize=13, fontweight='normal') 
                    plot_label (subplot, 'line', xlim, ylim,
                        x_percentage=0.44, y_percentage=0.87, color='black', x2_percentage=0.53, 
                        xlog=0, ylog=0, linestyle='-', linewidth=2)
            
        #LABELS        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.2, 
                    color='black', xlog=0, ylog=0, label='Haynes 2011', 
                    fontsize=13, fontweight='normal') 
        plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.225, 
                    color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15) 
           
        plot_label (subplot, 'label', xlim, ylim, 
                x_percentage=0.15, y_percentage=0.1, color='black', xlog=0, ylog=0, 
                label=prefix_this_model, fontsize=13, fontweight='normal') 
        plot_label (subplot, 'line', xlim, ylim,
                x_percentage=0.04, y_percentage=0.12, color='red', x2_percentage=0.13, 
                xlog=0, ylog=0, linestyle='-', linewidth=2)
        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.55, y_percentage=0.9, 
                    color='black', xlog=0, ylog=0, label='HI/Lr', 
                    fontsize=13, fontweight='normal')   
            
    plt.tight_layout()
    plt.savefig('./fig/plots_HI_over_Lr_vs_HI.pdf')   
    pdf.savefig()
    plt.close()

#end gas fraction

def HI_over_Lr_vs_HI_bins(G_MR, G_MRII, ThisRedshiftList, pdf):
   
    for ii in range(0,len(ThisRedshiftList)):        
        
        char_redshift="%0.2f" % ThisRedshiftList[ii]
        
        xlim=[-1.4,1.4]
        ylim=[0.0,0.59]
      
        bin=0.25
       
        fig = plt.figure(figsize=(three_two_size_small[0],three_two_size_small[1]))
        grid = gridspec.GridSpec(3, 2)
        grid.update(wspace=0.0, hspace=0.0)

       
                    
        Nbins_HI=6
        bin_HI=0.5        
        log_MHI=np.arange(8.0,11.0,bin_HI)
        bin_hist=0.25
        
        for jj in range(0,Nbins_HI):
       
            subplot=plt.subplot(grid[jj])
            subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
            
            subplot.xaxis.set_major_locator(MultipleLocator(0.5))    
            subplot.xaxis.set_minor_locator(MultipleLocator(0.1))            
            subplot.yaxis.set_minor_locator(MultipleLocator(0.02))  
            subplot.yaxis.set_major_locator(MultipleLocator(0.1))
                        
            xlab='$\mathrm{log_{10}}(M_{\mathrm{HI}}/L_{\mathrm{r}})$'
            ylab='fraction'   
            if jj==1 or jj == 3 or jj==5:
                plt.tick_params(axis='y', which='both', left='on', labelleft='off')
                ylab=''              
            subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   

            
            #MODEL
            if(log_MHI[jj]>=9.5):
                (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
                G0=G_MR[sel]  
            else:
                (sel)=select_current_redshift(G_MRII, ThisRedshiftList, ii, FullSnapshotList_MRII)        
                G0=G_MRII[sel]
                       
            if(opt_rings==1):
                HI=(np.log10(G0['ColdGas']*(1.-G0['H2fraction'])*1.e10/Hubble_h))              
            else:            
                HI=(np.log10(G0['ColdGas']*0.54*1.e10/Hubble_h))
            Lr=mag_to_lum(G0['MagDust'][:,17])           
            HI_over_LR=np.log10((10**HI)/Lr)            
            
            sel=(HI>log_MHI[jj]-bin_HI/2.) & (HI<log_MHI[jj]+bin_HI/2)    
            bin_arr=np.arange(xlim[0],xlim[1]+bin_hist,bin_hist)
            hist=np.histogram(HI_over_LR[sel], bins=bin_arr, range=(xlim[0],xlim[1])) 
            
            x_axis=hist[1][0:len(hist[1][:])-1]+bin_hist/2.           
            y_axis=hist[0]/np.sum(hist[0])     
            subplot.plot(x_axis,y_axis, color='red', linewidth=2, linestyle='-') 
                   
        
            #OBSERVATIONS
            file=Datadir+'Haynes_newfile'
            fits_table=fits.open(file)
            haynes = fits_table[1]

            haynes_MHI=haynes.data['HI']
            haynes_Magr=haynes.data['mr']-5*(np.log10(haynes.data['distance']*1.e6)-1)
            haynes_Lr=mag_to_lum(haynes_Magr)
            haynes_HI_over_LR=np.log10(10**haynes_MHI/haynes_Lr)
            
            sel=( (haynes_MHI>log_MHI[jj]-bin_HI/2.) & (haynes_MHI<log_MHI[jj]+bin_HI/2.) )   
            bin_arr=np.arange(xlim[0],xlim[1]+bin_hist,bin_hist)
            hist=np.histogram(haynes_HI_over_LR[sel], bins=bin_arr, range=(xlim[0],xlim[1])) 
            
            obs_x=hist[1][0:len(hist[1][:])-1]+bin_hist/2.           
            obs_y=hist[0]/np.sum(hist[0])
            obs_y_err=np.sqrt(hist[0])/np.sum(hist[0])
            subplot.errorbar(obs_x,obs_y, obs_y_err, fmt='o', markersize=5, ecolor='blue', color='blue')
          
            #write observations
            '''file_name=('HI_over_LrvsHI_MassBins_' + "%0.2f" % (log_MHI[jj]-bin_HI/2.) + 
                       '_' + "%0.2f" % (log_MHI[jj]+bin_HI/2.) + '_z'+char_redshift+'.txt' )
            fa = open(Datadir+file_name, "w")
            fa.write("%d\n" % len(obs_x))
            for kk in range (0,len(obs_x)):
                if(kk<(len(obs_x)-1)):
                    bin=obs_x[kk+1]-obs_x[kk]
                else:
                    bin=obs_x[kk]-obs_x[kk-1]               
                fa.write("%0.2f " % ((obs_x[kk]+np.log10(0.7**2))-bin/2.0) + 
                         "%0.2f " % ((obs_x[kk]+np.log10(0.7**2))+bin/2.0) + 
                         "%0.2f " % (obs_y[kk]*0.7) + 
                         "%0.2f\n" % (obs_y_err[kk]))                 
            fa.close()  '''
        
            fits_table.close()

            #MCMC sample
            if opt_plot_MCMC_sample==1:
                file = (MCMCSampledir + 'mcmc_plus_obs_HI_over_LrvsHI_MassBins_' +  
                        "%0.2f" % (log_MHI[jj]-bin_HI/2.) + '_' + "%0.2f" % (log_MHI[jj]+bin_HI/2.) +
                        '_z'+char_redshift+'.txt') 
                if os.path.isfile(file):
                    obs = Table.read(file, format='ascii')      
                    subplot.plot(obs['col1'],obs['col4'], color='black', linewidth=2)           
                    if jj==0:
                        plot_label (subplot, 'label', xlim, ylim, 
                            x_percentage=0.12, y_percentage=0.55, color='black', xlog=0, ylog=0, 
                            label='MCMC sample', fontsize=13, fontweight='normal') 
                        plot_label (subplot, 'line', xlim, ylim,
                            x_percentage=0.04, y_percentage=0.57, color='black', x2_percentage=0.1, 
                            xlog=0, ylog=0, linestyle='-', linewidth=2)
            
            #LABELS              
            label="%0.2f" % (log_MHI[jj]-bin_HI/2.) + r'$<\mathrm{log_{10}}(M_{HI}[M_{\odot}])<$' + "%0.2f" % (log_MHI[jj]+bin_HI/2.)
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.85, 
                        color='black', xlog=0, ylog=0, label=label, 
                        fontsize=13, fontweight='normal') 
            
            if(jj==0):
                plot_label (subplot, 'label', xlim, ylim, 
                        x_percentage=0.12, y_percentage=0.75, color='black', xlog=0, ylog=0, 
                        label=prefix_this_model, fontsize=13, fontweight='normal') 
                plot_label (subplot, 'line', xlim, ylim,
                        x_percentage=0.04, y_percentage=0.78, color='red', x2_percentage=0.1, 
                        xlog=0, ylog=0, linestyle='-', linewidth=2)
        
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.12, y_percentage=0.65, 
                    color='black', xlog=0, ylog=0, label='Haynes (2011)', 
                    fontsize=13, fontweight='normal') 
                plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.07, y_percentage=0.675, 
                    color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.02) 
        
            
    plt.tight_layout()
    plt.savefig('./fig/plots_HI_over_Lr_vs_HI_bins.pdf')
    plt.savefig('./fig/HYW17_plots_HI_over_Lr_vs_HI_bins.pdf')
    pdf.savefig()
    plt.close()

#end gas fraction


def HI_MF(G_MR, Volume_MR, G_MRII, Volume_MRII, ThisRedshiftList, pdf):
    for ii in range(0,len(ThisRedshiftList)):        
        
        xlim=[7.0,11.5]
        ylim=[-6.0,0.0]
        bin=0.25
    
        plot_color=['red','purple']        
      
        fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
        subplot=plt.subplot()
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
            
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))
            
        xlab='$\mathrm{log_{10}}(\mathrm{M_{\mathrm{HI}}}[M_{\odot}])$'   
        ylab='$\mathrm{log_{10}}(\phi [\mathrm{Mpc^{-3}} \mathrm{log_{10}}(M_{\mathrm{HI}}^{-1})])$'     
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   
            
        #PREVIOUS MODELS       
        char_redshift="%0.2f" % ThisRedshiftList[ii]
        if do_previous_model1==1: 
            file = file_previous_model1+'_coldgas_MF.txt' 
            model = Table.read(file, format='ascii')
            subplot.plot(model['col1']-2.*np.log10(Hubble_h_WMAP7),model['col2']+3.*np.log10(Hubble_h_WMAP7),
                         color='red',linestyle=linestyle_previous_model1, linewidth=2)
      
        if do_previous_model2==1: 
            file = file_previous_model2+'_coldgas_MF.txt' 
            model = Table.read(file, format='ascii')
            subplot.plot(model['col1']-2.*np.log10(Hubble_h),model['col2']+3.*np.log10(Hubble_h),
                         color='red',linestyle=linestyle_previous_model2, linewidth=2) 
            
        #MODEL
        #MR
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)         
        G0_MR=G_MR[sel]   
        
        if(opt_rings==1):
            G0_MR=G0_MR[(G0_MR['ColdGas']>0.) & (G0_MR['H2fraction']<1.)]
            HI=(np.log10(G0_MR['ColdGas']*(1.-G0_MR['H2fraction'])*1.e10/Hubble_h))  
        else:
            G0_MR=G0_MR[G0_MR['ColdGas']>0.]
            HI=(np.log10(G0_MR['ColdGas']*0.54*1.e10/Hubble_h))  
        
        bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
        hist_MR=np.histogram(HI, bins=bin_arr, range=(xlim[0],xlim[1]))   
       
        #MRII
        if(MRII==1):
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, ii, FullSnapshotList_MRII)         
            G0_MRII=G_MRII[sel]             
            if(opt_rings==1):
                G0_MRII=G0_MRII[(G0_MRII['ColdGas']>0.) & (G0_MRII['H2fraction']<1.)]
                HI=(np.log10(G0_MRII['ColdGas']*(1.-G0_MRII['H2fraction'])*1.e10/Hubble_h))  
            else:
                G0_MRII=G0_MRII[G0_MRII['ColdGas']>0.]
                HI=(np.log10(G0_MRII['ColdGas']*0.54*1.e10/Hubble_h)) 
        
            bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
            hist_MRII=np.histogram(HI, bins=bin_arr, range=(xlim[0],xlim[1])) 
            
        #join MR+MRII & plot     
        if(MRII==1):
            cut_MR_MRII=9.5
            (x_axis,y_axis)=plot_joint_MR_MRII(hist_MR, hist_MRII, cut_MR_MRII, Volume_MR, Volume_MRII, 
                                               bin, subplot, color='red',linewidth=2, linestyle='-')
        else:
            x_axis=hist_MR[1][0:len(hist_MR[1][:])-1]+bin/2.           
            hist_MR=hist_MR[0]       
            y_axis=np.log10(hist_MR/(Volume_MR*bin))
            subplot.plot(x_axis,y_axis, color='red', linewidth=2, linestyle='-') 
        
        #OBSERVATIONS
        h=0.75
        file = Datadir + 'zwaan2005.txt'       
        obs = Table.read(file, format='ascii')      
        obs_x = obs['col1']-2.*np.log10(h) 
        obs_y = obs['col2']     
        obs_y_err = [-obs['col3'],obs['col4']]
       
        subplot.errorbar(obs_x, obs_y, yerr=obs_y_err,
                 fmt='o', markersize=5, ecolor='blue', color='blue')
        
        
        file = Datadir + 'haynes2011_gmf.txt'       
        obs = Table.read(file, format='ascii')      
        obs_x = obs['col1']
        obs_y = obs['col2']
        obs_y_err = [-obs['col3'],obs['col4']]
       
        subplot.errorbar(obs_x, obs_y, yerr=obs_y_err,
                 fmt='o', markersize=5, ecolor='limegreen', color='limegreen')
         
        #MCMC sample
        if opt_plot_MCMC_sample==1:
            file = MCMCSampledir + 'mcmc_plus_obs_ColdGasMassFunction_z'+char_redshift+'.txt' 
            if os.path.isfile(file):
                obs = Table.read(file, format='ascii')      
                subplot.plot(obs['col1']-2.*np.log10(Hubble_h),np.log10(obs['col4'])+3.*np.log10(Hubble_h), 
                             color='black', linewidth=2)
                           
                if ii==len(ThisRedshiftList)-1:
                    plot_label (subplot, 'label', xlim, ylim, 
                        x_percentage=0.55, y_percentage=0.85, color='black', xlog=0, ylog=0, 
                        label='MCMC sample', fontsize=13, fontweight='normal') 
                    plot_label (subplot, 'line', xlim, ylim,
                        x_percentage=0.44, y_percentage=0.87, color='black', x2_percentage=0.53, 
                        xlog=0, ylog=0, linestyle='-', linewidth=2)    
            
        #LABELS        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.5, 
                    color='black', xlog=0, ylog=0, label='Zwaan 2005', 
                    fontsize=13, fontweight='normal') 
        plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.525, 
                    color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15) 
        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.4, 
                    color='black', xlog=0, ylog=0, label='Haynes 2011', 
                    fontsize=13, fontweight='normal') 
        plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.425, 
                    color='limegreen', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15) 
        
        plot_label_three_models (subplot, xlim, ylim, position='bottom_left')
        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.75, y_percentage=0.8, 
                    color='black', xlog=0, ylog=0, label='HI MF', 
                    fontsize=13, fontweight='normal')   
        
    plt.tight_layout()
    plt.savefig('./fig/plots_HI_MF.pdf')
    plt.savefig('./fig/HYW17_plots_HI_MF.pdf')
    pdf.savefig()
    plt.close()
   

#end HI_MF






















def coldgas_vs_stellarmass(G_MR, ThisRedshiftList, pdf):
  
    plot_color=['red','purple']        
   
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    
    #OBSERVATIONS READ  
    file = Datadir+"/Saintonge2016_gasfraction.txt"   
    Saint16 = Table.read(file, format='ascii')  
    Saint16_mass=(Saint16['mass_bin_low']+Saint16['mass_bin_high'])/2.  
  
     
    
    for ii in range(0,len(ThisRedshiftList)):        
               
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]    
        #G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (G0_MR['ColdGas']>0.) & 
        #            (G0_MR['Vvir']>120.) & (G0_MR['BulgeMass']/G0_MR['StellarMass']<0.15)]
        G0_MR=G0_MR[(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)>10.) & (G0_MR['ColdGas']>0.) & 
                    (np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h))>-11.)] 
        #G0_MR=G0_MR[(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)>10.) & (G0_MR['ColdGas']>0.)]
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
        
        xlim=[9.5,11.5]
        ylim=[-2.0,0.5]
        bin=0.25
           
        subplot=plt.subplot()    
        subplot.set_ylim(ylim), subplot.set_xlim(xlim) 
        ylab='$\log_{10}(M_{\mathrm{gas}}/M_*)$' 
        subplot.set_ylabel(ylab, fontsize=14) 
        xlab='$\log_{10}(M_*[M_{\odot}])$'           
        subplot.set_xlabel(xlab, fontsize=14)   
          
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))
        
            
        #MODEL
        sel=G0_MR['ColdGas']>0.
        Fraction=np.log10(G0_MR['ColdGas']*1.e10/Hubble_h)-StellarMass 
         
           
        (x_binned, median,mean,pc16,pc84,rms)=median_and_percentiles(bin,xlim[0],xlim[1],StellarMass[sel],Fraction[sel])    
        sel=(median!=0)        
        subplot.plot(x_binned[sel],median[sel],color=plot_color[ii],linewidth=2)     
        subplot.plot(x_binned[sel],pc16[sel],color=plot_color[ii],linewidth=2,linestyle='--')
        subplot.plot(x_binned[sel],pc84[sel],color=plot_color[ii],linewidth=2,linestyle='--')
        
        
        #OBSERVATIONS PLOT  
        y_err=np.zeros(len(Saint16['fHI']),dtype=np.float32)
        Saint16_H2plusHI=Saint16['fH2']+Saint16['fHI']
        Saint16_H2plusHI_err=Saint16['fH2_err']+Saint16['fHI_err']
       
        y_err=[np.log10(Saint16_H2plusHI/(Saint16_H2plusHI-Saint16_H2plusHI_err)),
               np.log10((Saint16_H2plusHI+Saint16_H2plusHI_err)/Saint16_H2plusHI)]
        subplot.errorbar(Saint16_mass, np.log10(Saint16_H2plusHI),xerr=0.12,yerr=y_err,
                         fmt='o', markersize=5, ecolor='blue', color='blue')    
               
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.15, y_percentage=0.2, 
                    color='black', xlog=0, ylog=0, label=prefix_this_model, 
                    fontsize=13, fontweight='normal')             
        plot_label (subplot, 'line', xlim, ylim,x_percentage=0.05, y_percentage=0.2175, 
                    color='red', x2_percentage=0.12, xlog=0, ylog=0, linestyle='-', linewidth=2)

        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.15, y_percentage=0.12, 
                    color='black', xlog=0, ylog=0, label='Saintonge 2016', 
                    fontsize=13, fontweight='normal') 
        plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.12, y_percentage=0.14, 
                    color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.075)     

        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.6, y_percentage=0.85, 
                    color='black', xlog=0, ylog=0, label='$M_{\mathrm{cold}}/M_*$', 
                    fontsize=15, fontweight='normal') 


            
                
    plt.tight_layout()
    plt.savefig('./fig/plots_coldgas_vs_stellarmass_vs_stellarmass.pdf')   
    pdf.savefig()
    plt.close()
    
    
    
    
    
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
       
    for ii in range(0,len(ThisRedshiftList)):        
               
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]    
        #G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (G0_MR['ColdGas']>0.) & 
        #            (G0_MR['Vvir']>120.) & (G0_MR['BulgeMass']/G0_MR['StellarMass']<0.15)]
        G0_MR=G0_MR[(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)>10.) & (G0_MR['ColdGas']>0.) & 
                    (np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h))>-11.)] 
        #G0_MR=G0_MR[(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)>10.) & (G0_MR['ColdGas']>0.)]
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
        
        xlim=[9.5,11.5]
        ylim=[8.0,12.0]
        bin=0.25
           
        subplot=plt.subplot()    
        subplot.set_ylim(ylim), subplot.set_xlim(xlim) 
        ylab='$\log_{10}(M_{\mathrm{gas}}/M_*)$' 
        subplot.set_ylabel(ylab, fontsize=14) 
        xlab='$\log_{10}(M_*[M_{\odot}])$'           
        subplot.set_xlabel(xlab, fontsize=14)   
          
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))
        
            
        #MODEL
        sel=G0_MR['ColdGas']>0.
        Fraction=np.log10(G0_MR['ColdGas']*1.e10/Hubble_h)
         
           
        (x_binned, median,mean,pc16,pc84,rms)=median_and_percentiles(bin,xlim[0],xlim[1],StellarMass[sel],Fraction[sel])    
        sel=(median!=0)        
        subplot.plot(x_binned[sel],median[sel],color=plot_color[ii],linewidth=2)     
        subplot.plot(x_binned[sel],pc16[sel],color=plot_color[ii],linewidth=2,linestyle='--')
        subplot.plot(x_binned[sel],pc84[sel],color=plot_color[ii],linewidth=2,linestyle='--')
                   
                
    plt.tight_layout()   
    pdf.savefig()
    plt.close()


#end

















def main_sequence(G_MR, ThisRedshiftList, pdf):
           
    xlim=[8.5,11.5]
    ylim=[-2.5, 4]   
    bin=[0.25,0.25]
    #bin=[0.25,0.5]
    Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]

    fig = plt.figure(figsize=(one_five_size_large[0],one_five_size_large[1]))
    grid = gridspec.GridSpec(1,5)
    grid.update(wspace=0.0, hspace=0.0)

    for ii in range(0,len(ThisRedshiftList)):
        
        char_redshift="%0.1f" % ThisRedshiftList[ii]
        
        subplot=plt.subplot(grid[ii])
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        
        xlab='$\mathrm{log_{10}}(M_*[M_{\odot}])$'      
        if ii==0:
            ylab='$\mathrm{log_{10}}(\mathrm{SFR}[M_{\odot}yr^{-1}])$'
        else:
            ylab=''      
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)
        
        subplot.text(xlim[0]+2.,ylim[0]+.3,'z='+char_redshift, fontsize=14, fontweight='normal')
        
        #if ii==2:
        #    subplot.text(xlim[0]+0.3,ylim[0]+0.5,'MS evo', fontsize=16, fontweight='normal')
        
        
    
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
        subplot.yaxis.set_minor_locator(MultipleLocator(0.25))
    
        if ii>0:
            plt.tick_params(axis='y', which='both', left='on', labelleft='off')
           
        
        #MODEL
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)
        
        G0_MR=G_MR[sel]   
        G0_MR=G0_MR[(G0_MR['StellarMass']>0.)]
        log_StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])      
        SFR=G0_MR['Sfr']       
        log_SFR=np.log10(SFR)
       
        width=0.25  
        (slope,b)=get_slope(9.0,-2,11.0,-1.)
        sel=SFR < 0.001      
        log_SFR[sel]=(np.random.randn(len(log_SFR[sel]))*width*(slope*log_StellarMass[sel]+b) +
                      (slope*log_StellarMass[sel]+b))
      
        
        Ngals=len(G0_MR)
      
        H, xedges, yedges = np.histogram2d(log_StellarMass, log_SFR, bins=Nbins, 
                                           range=[[xlim[0]-bin[0]/2., xlim[1]+bin[0]/2.], 
                                                  [ylim[0]-bin[1]/2., ylim[1]+bin[1]/2.]])            
        extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
        plt.subplots_adjust(bottom=0.15, left=0.15)        
        mylevels = np.log10(np.logspace(-2.0, 0.0, num=10))  
        H = zoom(H, 20)        
        H=np.log10(H/np.amax(H))           
        cont=plt.contourf(H.transpose()[::], origin='lower', cmap='Greys_r', levels=mylevels, extent=extent)   
        
        if ii==len(ThisRedshiftList)-1:          
            plt.colorbar(format='%0.1f') 
            
        #OBSERVATIONS
        #values at all_z
        
        #ELBAZ2007
        obs_slope_elbaz2007 =[0.77, -99.0, 0.9, -99.0, -99.0]
        obs_offset_elbaz2007=[np.log10(8.7)-(0.77*11.), -99.0, np.log10(7.2)-9, -99.0, -99.0]
        obs_offset_low_elbaz2007=[np.log10(5.0)-(0.77*11.), -99.0, np.log10(3.6)-9, -99.0, -99.0]
        obs_offset_high_elbaz2007=[np.log10(16.1)-(0.77*11.), -99.0, np.log10(14.4)-9, -99.0, -99.0]
        
        #KARIM2011
        file = Datadir + 'karim2011_sfr_mass_sf.txt'       
        karim2011 = Table.read(file, format='ascii') 
        karim_low_z_limit        = karim2011['col4']
        karim_medium_mass        = karim2011['col3']
        karim_sfr                = karim2011['col19']
        karim_sfr_error_up   = karim2011['col20']
        karim_sfr_error_down = karim2011['col21']
        log_karim_sfr_error_up=np.log10((karim_sfr+karim_sfr_error_up)/karim_sfr)
        log_karim_sfr_error_down=np.log10(karim_sfr/(karim_sfr-karim_sfr_error_down))
        
        obs_x=np.arange(xlim[0], xlim[1], 0.01)
        
        if ThisRedshiftList[ii]==0.0:  
            #SDSS-DR7
            file=Datadir+'dr7_gal_final.fit'
            fits_table=fits.open(file)
            dr7_gal_final = fits_table[1]
            dr7_gal_final = dr7_gal_final.data[(dr7_gal_final.data['z'] > SDSS_min_z) & 
                                          (dr7_gal_final.data['z'] < SDSS_max_z)]
         
            bin=[0.25,0.25]
            hubble_h_WMAP7=0.7
            Nbinsx=(xlim[1]-xlim[0])/bin[0]
            Nbinsy=(ylim[1]-ylim[0])/bin[1]          
            obs_mstar = dr7_gal_final['jarle_median_mass'] 
            obs_sfr   = dr7_gal_final['median_sfr']
  
            mag=dr7_gal_final['jarle_mag'][:,2] #absolute mag
            zz=dr7_gal_final['z'] 
            max_d=10**((17.6-mag)/5.+1.)/1.e6 #apparent
            weight = 1./max_d**3

      
            frac=np.zeros([int(Nbinsx),int(Nbinsy)],dtype=np.float32) 
            mstar_c=np.arange(xlim[0],xlim[1],bin[0])
            sfr_c=np.arange(ylim[0],ylim[1],bin[1])
           
            xind=0
            for jj in range (0,int(Nbinsx)):
                yind=0  
                xid=xlim[0]+bin[0]*jj                 
                for kk in range (0,int(Nbinsy)):
                    yid=ylim[0]+bin[1]*kk
                    sel=((obs_mstar > xid) & (obs_mstar < xid+bin[0]) &
                         (obs_sfr > yid) & (obs_sfr < yid+bin[1]))                      
                    if len(sel)==0:
                        frac[xind,yind]=0.
                    else: 
                        frac[xind,yind]=np.sum(weight[sel])                             
                    yind+=1
                xind+=1
          
            H=frac
            extent = [xlim[0], xlim[1],ylim[0], ylim[1]]       
            plt.subplots_adjust(bottom=0.15, left=0.15)     
            mylevels = np.log10(np.logspace(-2.0, 0.0, num=10))  
            H = zoom(H, 20)        
            H=np.log10(H/np.amax(H))       
            cont=plt.contour(H.transpose()[::], origin='lower', colors='black', 
                             linewidths=0.5, linestyles='-', levels=mylevels, extent=extent)          
                  
            fits_table.close()
            
            #ELBAZ2007
            obs_y=obs_x*obs_slope_elbaz2007[ii] + obs_offset_elbaz2007[ii] 
            subplot.plot(obs_x, obs_y, color='firebrick', linewidth=2)            
            obs_y=obs_x*obs_slope_elbaz2007[ii] + obs_offset_low_elbaz2007[ii]
            subplot.plot(obs_x, obs_y, color='firebrick', linewidth=2, linestyle='--')            
            obs_y=obs_x*obs_slope_elbaz2007[ii] + obs_offset_high_elbaz2007[ii]
            subplot.plot(obs_x, obs_y, color='firebrick', linewidth=2, linestyle='--')
            
            plot_label (subplot, 'label', xlim, ylim, 
                    x_percentage=0.15, y_percentage=0.92, color='black', xlog=0, ylog=0, 
                    label='Elbaz 2007', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim,
                    x_percentage=0.04, y_percentage=0.94, color='firebrick', x2_percentage=0.13, 
                    xlog=0, ylog=0, linestyle='-', linewidth=2)
            plot_label (subplot, 'label', xlim, ylim, 
                    x_percentage=0.15, y_percentage=0.85, color='black', xlog=0, ylog=0, 
                    label='SDSS-DR7', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim,
                    x_percentage=0.04, y_percentage=0.87, color='black', x2_percentage=0.13, 
                    xlog=0, ylog=0, linestyle='-', linewidth=2)
            
            
        if ThisRedshiftList[ii]==0.4: 
            #KARIM2011
            sel=(karim_low_z_limit==0.2) & (karim_medium_mass>8.8)                
            subplot.errorbar(karim_medium_mass[sel], np.log10(karim_sfr[sel]), 
                             [log_karim_sfr_error_down[sel], log_karim_sfr_error_up[sel]], 
                             mfc='white', markeredgecolor='red', color='red', fmt='o', markersize=5)
            sel=(karim_low_z_limit==0.4) & (karim_medium_mass>8.9)                
            subplot.errorbar(karim_medium_mass[sel], np.log10(karim_sfr[sel]), 
                             [log_karim_sfr_error_down[sel], log_karim_sfr_error_up[sel]], 
                             markeredgecolor='red',color='red', fmt='o', markersize=5)
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.93, 
                        color='black', xlog=0, ylog=0, label='Karim 2011 - 0.2<z<0.4', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.95, 
                        color='red', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.1,mfc='white')
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.88, 
                        color='black', xlog=0, ylog=0, label='Karim 2011 - 0.4<z<0.6', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.90, 
                        color='red', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.1)
            
        if ThisRedshiftList[ii]==1.0:  
             #ELBAZ2007
            obs_y=obs_x*obs_slope_elbaz2007[ii] + obs_offset_elbaz2007[ii]
            subplot.plot(obs_x, obs_y, color='firebrick', linewidth=2)            
            obs_y=obs_x*obs_slope_elbaz2007[ii] + obs_offset_low_elbaz2007[ii]
            subplot.plot(obs_x, obs_y, color='firebrick', linewidth=2, linestyle='--')            
            obs_y=obs_x*obs_slope_elbaz2007[ii] + obs_offset_high_elbaz2007[ii]
            subplot.plot(obs_x, obs_y, color='firebrick', linewidth=2, linestyle='--')
            
            #KARIM2011
            sel=(karim_low_z_limit==0.8) & (karim_medium_mass>9.1)                
            subplot.errorbar(karim_medium_mass[sel], np.log10(karim_sfr[sel]), 
                             [log_karim_sfr_error_down[sel], log_karim_sfr_error_up[sel]],
                             mfc='white', markeredgecolor='red', color='red', fmt='o', markersize=5)
            sel=(karim_low_z_limit==1.0) & (karim_medium_mass>9.3)                
            subplot.errorbar(karim_medium_mass[sel], np.log10(karim_sfr[sel]), 
                             [log_karim_sfr_error_down[sel], log_karim_sfr_error_up[sel]],  
                             markeredgecolor='red',color='red', fmt='o', markersize=5)
            
            #Whitaker2013
            file = Datadir + 'whitaker2013_mass_vs_sfr_z0.5_1.0.txt'       
            obs = Table.read(file, format='ascii') 
            subplot.errorbar(obs['col1'], obs['col2'], obs['col3'], mfc='white', 
                             markeredgecolor='blue', color='blue', fmt='o', markersize=5)
            file = Datadir + 'whitaker2013_mass_vs_sfr_z1.0_1.5.txt'       
            obs = Table.read(file, format='ascii') 
            subplot.errorbar(obs['col1'], obs['col2'], obs['col3'], color='blue', fmt='o', markersize=5)
            
            
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.93, 
                        color='black', xlog=0, ylog=0, label='Karim 2011 - 0.8<z<1.0', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.95, 
                        color='red', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.1,mfc='white')
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.88, 
                        color='black', xlog=0, ylog=0, label='Karim 2011 - 1.0<z<1.2', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.90, 
                        color='red', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.1)
        
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.83, 
                        color='black', xlog=0, ylog=0, label='Whitaker 2013 - 0.5<z<1.0', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.85, 
                        color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.1,mfc='white')
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.78, 
                        color='black', xlog=0, ylog=0, label='Whitaker 2013 - 1.0<z<1.5', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.80, 
                        color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.1)
        
            
            
        if ThisRedshiftList[ii]==2.0:  
            #KARIM2011
            sel=(karim_low_z_limit==1.6) & (karim_medium_mass>9.6)                
            subplot.errorbar(karim_medium_mass[sel], np.log10(karim_sfr[sel]), 
                             [log_karim_sfr_error_down[sel], log_karim_sfr_error_up[sel]],
                             mfc='white', markeredgecolor='red', color='red', fmt='o', markersize=5)
            sel=(karim_low_z_limit==2.0) & (karim_medium_mass>9.8)                
            subplot.errorbar(karim_medium_mass[sel], np.log10(karim_sfr[sel]), 
                             [log_karim_sfr_error_down[sel], log_karim_sfr_error_up[sel]],  
                             markeredgecolor='red',color='red', fmt='o', markersize=5)
            
            #Whitaker2013
            file = Datadir + 'whitaker2013_mass_vs_sfr_z1.5_2.0.txt'       
            obs = Table.read(file, format='ascii') 
            subplot.errorbar(obs['col1'], obs['col2'], obs['col3'], mfc='white', 
                             markeredgecolor='blue', color='blue', fmt='o', markersize=5)
            file = Datadir + 'whitaker2013_mass_vs_sfr_z2.0_2.5.txt'       
            obs = Table.read(file, format='ascii') 
            subplot.errorbar(obs['col1'], obs['col2'], obs['col3'], color='blue', fmt='o', markersize=5)
            
            file = Datadir + '/Shivaei2016_sfr.txt'       
            obs = Table.read(file, format='ascii') 
            y_err= [np.log10(obs['S16']/(obs['S16']-obs['S16err'])),np.log10((obs['S16']+obs['S16err'])/obs['S16'])]
            subplot.errorbar(np.log10(obs['mass']), np.log10(obs['S16']),yerr=y_err, 
                             color='limegreen', markeredgecolor='limegreen',fmt='o', markersize=5)
                     
            '''plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.92, 
                        color='black', xlog=0, ylog=0, label='Karim 2011 - 1.6<z<2.0', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.945, 
                        color='red', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15,mfc='white')
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.85, 
                        color='black', xlog=0, ylog=0, label='Karim 2011 - 2.0<z<2.5', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.875, 
                        color='red', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15)
        
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.78, 
                        color='black', xlog=0, ylog=0, label='Whitaker 2013 - 1.5<z<2.0', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.805, 
                        color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15,mfc='white')
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.71, 
                        color='black', xlog=0, ylog=0, label='Whitaker 2013 - 2.0<z<2.5', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.735, 
                        color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15)
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.64, 
                        color='black', xlog=0, ylog=0, label='Shivaei 2016 - 1.4<z<2.6', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.665, 
                        color='limegreen', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15)'''
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.93, 
                        color='black', xlog=0, ylog=0, label='Karim 2011 - 1.6<z<2.0', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.95, 
                        color='red', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.1,mfc='white')
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.88, 
                        color='black', xlog=0, ylog=0, label='Karim 2011 - 2.0<z<2.5', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.90, 
                        color='red', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.1)
        
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.83, 
                        color='black', xlog=0, ylog=0, label='Whitaker 2013 - 1.5<z<2.0', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.85, 
                        color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.1,mfc='white')
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.78, 
                        color='black', xlog=0, ylog=0, label='Whitaker 2013 - 2.0<z<2.5', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.79, 
                        color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.1)
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.73, 
                        color='black', xlog=0, ylog=0, label='Shivaei 2016 - 1.4<z<2.6', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.75, 
                        color='limegreen', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.1)
            
            
        if ThisRedshiftList[ii]==3.0:
            #KARIM2011
            sel=(karim_low_z_limit==2.5) & (karim_medium_mass>10.0)                
            subplot.errorbar(karim_medium_mass[sel], np.log10(karim_sfr[sel]), 
                             [log_karim_sfr_error_down[sel], log_karim_sfr_error_up[sel]],
                             markeredgecolor='red',color='red', fmt='o', markersize=5)
         
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.93, 
                        color='black', xlog=0, ylog=0, label='Karim 2011 - 2.5<z<3.0', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.95, 
                        color='red', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.1)
    #endfor
        
    plt.tight_layout()
    plt.savefig('./fig/plots_main_sequence.pdf')
    plt.savefig('./fig/HYW17_plots_main_sequence.pdf')
    pdf.savefig()
    plt.close()
#endif stellar_mass_vs_sfr




def ur_vs_r(G_MR, ThisRedshiftList, pdf):
           
    xlim=[-23.,-14.]
    ylim=[1.0, 3.2]   
    bin=[0.5,0.025]
    Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]

    fig = plt.figure(figsize=(5,4))
    subplot=plt.subplot()
    
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)      
    subplot.set_xlabel('r', fontsize=14), subplot.set_ylabel('u-r', fontsize=14) 
    
    majorFormatter = FormatStrFormatter('%2d')
    subplot.xaxis.set_major_locator(MultipleLocator(2))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      

    for ii in range(0,len(ThisRedshiftList)):                     
             
        subplot.text(xlim[0]+0.2,ylim[0]+0.25,'z=0', fontsize=16, fontweight='normal')     
        subplot.text(xlim[0]+0.2,ylim[0]+0.05,'r vs u-r cut', fontsize=16, fontweight='normal')
                  
        #MODEL
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)
        
        G0_MR=G_MR[sel]   
        G0_MR=G0_MR[(G0_MR['MagDust'][:,15]<99.) & (G0_MR['MagDust'][:,17]<99.)]        
        color_ur=G0_MR['MagDust'][:,15]-G0_MR['MagDust'][:,17]  
        Magr=G0_MR['MagDust'][:,17]         
        Ngals=len(G0_MR)
      
        H, xedges, yedges = np.histogram2d(Magr, color_ur, bins=Nbins)            
        extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
        plt.subplots_adjust(bottom=0.15, left=0.15)        
        mylevels = np.linspace(1., Nbins[0]*4., Nbins[0]*4.)*Ngals/(Nbins[0]**2/0.2)        
        H = zoom(H, 20)        
        cont=plt.contourf(H.transpose()[::], origin='lower', cmap='Greys_r', levels=mylevels, extent=extent)        
        plt.colorbar(format='%d') 
        
        #BestFit Cut   
        bin=0.01
        x_arr=np.arange(xlim[0],xlim[1]+bin,bin)        
        y_arr=(offset_color_cut[0]-slope_color_cut[0]*np.tanh((x_arr+18.07)/1.09))
        subplot.plot(x_arr,y_arr,color='red', linestyle='-', linewidth=2)  
        
        #OBSERVATIONAL CUT 
        Nbin=0.01
        x_arr=np.arange(xlim[0],xlim[1]+Nbin,Nbin)       
        y_arr=2.06-0.244*np.tanh((x_arr+20.07)/1.09)        
        subplot.plot(x_arr,y_arr,color='blue', linestyle='--', linewidth=2) 
        
        #LABELS
        plot_label (subplot, 'label', xlim, ylim, 
                    x_percentage=0.15, y_percentage=0.9, color='black', xlog=0, ylog=0, 
                    label='Baldry 2004', fontsize=13, fontweight='normal') 
        plot_label (subplot, 'line', xlim, ylim,
                    x_percentage=0.04, y_percentage=0.92, color='blue', x2_percentage=0.13, 
                    xlog=0, ylog=0, linestyle='--', linewidth=2)
        
        plot_label (subplot, 'label', xlim, ylim, 
                    x_percentage=0.15, y_percentage=0.82, color='black', xlog=0, ylog=0, 
                    label='Best Fit cut', fontsize=13, fontweight='normal') 
        plot_label (subplot, 'line', xlim, ylim,
                    x_percentage=0.04, y_percentage=0.84, color='red', x2_percentage=0.13, 
                    xlog=0, ylog=0, linestyle='-', linewidth=2)
        
    #endfor
        
    plt.tight_layout()
    plt.savefig('./fig/plots_ur_vs_r.pdf')
    pdf.savefig()
    plt.close()
#endif ur_vs_r






def UVJ_colour(G_MR, ThisRedshiftList, pdf):
           
    xlim=[-0.5,2.5]
    ylim=[-0.5, 2.5]   
    bin=[0.05,0.05]
    Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]

    fig = plt.figure(figsize=(12,4))
    grid = gridspec.GridSpec(1, 4)
    grid.update(wspace=0.0, hspace=0.0)

    for ii in range(0,len(ThisRedshiftList)):
        
        char_redshift="%0.1f" % ThisRedshiftList[ii]
        
        subplot=plt.subplot(grid[ii])
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        if ii==0:
            ylabel='U-V'
        else:
             ylabel=''
        subplot.set_xlabel('V-J', fontsize=14), subplot.set_ylabel(ylabel, fontsize=14)
        
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
        subplot.yaxis.set_minor_locator(MultipleLocator(0.25)) 
        
        if ii>0:
            plt.tick_params(axis='y', which='both', left='on', labelleft='off')
        
        subplot.text(xlim[0]+2.,ylim[0]+.2,'z='+char_redshift, fontsize=16, fontweight='normal')
        
        if ii==2:
            subplot.text(xlim[0]+0.2,ylim[0]+0.2,'UVJ cut', fontsize=16, fontweight='normal')
                  
        #MODEL
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)
        
        G0_MR=G_MR[sel]   
        G0_MR=G0_MR[(G0_MR['MagDust'][:,0]<99.) & (G0_MR['MagDust'][:,2]<99.) & (G0_MR['MagDust'][:,7]<99.)]        
        color_UV=G0_MR['MagDust'][:,0]-G0_MR['MagDust'][:,2]  
        color_VJ=G0_MR['MagDust'][:,2]-G0_MR['MagDust'][:,7]       
        Ngals=len(G0_MR)
      
        H, xedges, yedges = np.histogram2d(color_VJ, color_UV, bins=Nbins)            
        extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
        plt.subplots_adjust(bottom=0.15, left=0.15)        
        mylevels = np.linspace(1., Nbins[0], Nbins[0])*Ngals/(Nbins[0]**2/2.5)        
        H = zoom(H, 20)        
        cont=plt.contourf(H.transpose()[::], origin='lower', cmap='Greys_r', levels=mylevels, extent=extent)        
        
        if ii==len(ThisRedshiftList)-1:
            plt.colorbar(format='%d') 
        
        #BestFit Cut   
        bin=0.01
        slope=slope_color_cut[ii+1]
        offset=offset_color_cut[ii+1]
        minimum_y=minimum_y_color_cut[ii+1]
        
        x_arr=np.arange(xlim[0],xlim[1]+bin,bin)
        cut1=np.zeros(len(x_arr),dtype=np.float32)+minimum_y      
        cut2=x_arr*slope+offset
  
        sel1=x_arr<((minimum_y-offset)/slope)
        subplot.plot(x_arr[sel1],cut1[sel1],color='red', linestyle='-', linewidth=2)  
         
        sel2=x_arr > ((minimum_y-offset)/slope)
        subplot.plot(x_arr[sel2],cut2[sel2],color='red', linestyle='-', linewidth=2)  
        
        #OBSERVATIONAL CUT    
        Nbin=0.01
        slope=0.88
        if(ThisRedshiftList[ii]<1.):
            offset=0.69 
        else: 
            offset=0.59

        x_arr=np.arange(xlim[0],xlim[1]+Nbin,Nbin)
        cut1=np.zeros(len(x_arr),dtype=np.float32)+1.3   
        cut2=x_arr*slope+offset
  
        sel1=x_arr<((1.3-offset)/slope)
        subplot.plot(x_arr[sel1],cut1[sel1],color='blue', linestyle='--', linewidth=2)  
  
        sel2=x_arr>((1.3-offset)/slope) 
        subplot.plot(x_arr[sel2],cut2[sel2],color='blue', linestyle='--', linewidth=2)  
        
        
        #LABELS
        if ii==0:
            plot_label (subplot, 'label', xlim, ylim, 
                        x_percentage=0.15, y_percentage=0.15, color='black', xlog=0, ylog=0, 
                        label='Muzzin 2013', fontsize=13, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim,
                        x_percentage=0.04, y_percentage=0.17, color='blue', x2_percentage=0.13, 
                        xlog=0, ylog=0, linestyle='--', linewidth=2)
        
            plot_label (subplot, 'label', xlim, ylim, 
                        x_percentage=0.15, y_percentage=0.07, color='black', xlog=0, ylog=0, 
                        label='Best Fit cut', fontsize=13, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim,
                        x_percentage=0.04, y_percentage=0.09, color='red', x2_percentage=0.13, 
                        xlog=0, ylog=0, linestyle='-', linewidth=2)
    #endfor
        
    plt.tight_layout()
    plt.savefig('./fig/plots_UVJ_colour.pdf')
    pdf.savefig()
    plt.close()
#endif UVJ_colour




def UVJ_grid(G_MR, ThisRedshiftList, pdf):
           
    xlim=[-0.5,2.5]
    ylim=[-0.5, 2.5]   
    bin=[0.05,0.05]
    Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]

    fig = plt.figure(figsize=(16,12))
    grid = gridspec.GridSpec(3, 4)
    grid.update(wspace=0.1, hspace=0.1)

    MassArray=[9.25,9.75,10.25,10.75]
    MassBin=0.5
         
    grid_num=0
    for ii in range(0,len(ThisRedshiftList)):
        
        char_redshift="%0.1f" % ThisRedshiftList[ii]
        
        
        for jj in range (0,len(MassArray)):
            
            subplot=plt.subplot(grid[grid_num])            
            subplot.set_ylim(ylim), subplot.set_xlim(xlim)   
            
            if ii<len(ThisRedshiftList)-1:
                x_label=''
            else:
                x_label='V-J'
            if jj==0:    
                y_label='U-V'
            else:
                y_label=''
            subplot.set_xlabel(x_label, fontsize=14), subplot.set_ylabel(y_label, fontsize=14)           
           
            majorFormatter = FormatStrFormatter('%d')
            subplot.xaxis.set_major_locator(MultipleLocator(1))    
            subplot.xaxis.set_minor_locator(MultipleLocator(0.5))   
            subplot.yaxis.set_major_locator(MultipleLocator(1)) 
            subplot.yaxis.set_minor_locator(MultipleLocator(0.5)) 
        
            if ii<len(ThisRedshiftList)-1:
                subplot.tick_params(axis='x', which='both', bottom='on', labelbottom='off')              
            if jj>0:
                plt.tick_params(axis='y', which='both', left='on', labelleft='off') 
                  
            mass_low="%0.1f" % (MassArray[jj]-MassBin/2.) 
            mass_high="%0.1f" % (MassArray[jj]+MassBin/2.) 
            label=mass_low+"$<log_{10}(M_{\star} [M_{\odot}])<$"+mass_high
            if ii==0:
                subplot.text(xlim[0],ylim[1]+.2,label,fontsize=16, fontweight='normal')
               
            if(jj==(len(MassArray)-1)):
                subplot.text(xlim[1]+.2,ylim[0]+1.7,'z='+char_redshift, fontsize=18, fontweight='normal', rotation=270)
        
            #MODEL
            (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)
        
            G0_MR=G_MR[sel]   
            Nrandom=1000
           
            G0_MR=G0_MR[(G0_MR['MagDust'][:,0]<99.) & (G0_MR['MagDust'][:,2]<99.) & (G0_MR['MagDust'][:,7]<99.) &
                        (np.log10(G0_MR['StellarMass']*1.e10/Hubble_h) > MassArray[jj]-MassBin/2.) &
                        (np.log10(G0_MR['StellarMass']*1.e10/Hubble_h) < MassArray[jj]+MassBin/2.) ] 
            G0_MR=np.random.choice(G0_MR, size=Nrandom)
            color_UV=G0_MR['MagDust'][:,0]-G0_MR['MagDust'][:,2]  
            color_VJ=G0_MR['MagDust'][:,2]-G0_MR['MagDust'][:,7]       
            Ngals=len(G0_MR)
            

            color=['purple','blue','green','yellow','orange','red']
            SSFR=np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h))            
            sel=((SSFR<-8.0) & (SSFR>-8.5))
            subplot.scatter(color_VJ[sel], color_UV[sel], s=20, color=color[0]) 
            sel=((SSFR<-8.5) & (SSFR>-9.0))
            subplot.scatter(color_VJ[sel], color_UV[sel], s=20, color=color[1]) 
            sel=((SSFR<-9.0) & (SSFR>-9.5))
            subplot.scatter(color_VJ[sel], color_UV[sel], s=20, color=color[2]) 
            sel=((SSFR<-9.5) & (SSFR>-10.0))
            subplot.scatter(color_VJ[sel], color_UV[sel], s=20, color=color[3]) 
            sel=((SSFR<-10.0) & (SSFR>-10.5))
            subplot.scatter(color_VJ[sel], color_UV[sel], s=20, color=color[4]) 
            sel=(SSFR<-10.5)
            subplot.scatter(color_VJ[sel], color_UV[sel], s=20, color=color[5]) 
        
            '''H, xedges, yedges = np.histogram2d(color_VJ, color_UV, bins=Nbins)            
            extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
            plt.subplots_adjust(bottom=0.15, left=0.15)        
            mylevels = np.linspace(1., Nbins[0], Nbins[0])*Ngals/(Nbins[0]**2/2.5)        
            H = zoom(H, 20)        
            cont=plt.contourf(H.transpose()[::], origin='lower', cmap='Greys_r', levels=mylevels, extent=extent)'''        
         
            #if ii==len(ThisRedshiftList)-1:
            #    plt.colorbar(format='%d') 
        
            #OBSERVATIONAL CUT    
            Nbin=0.01
            slope=0.88
            if(ThisRedshiftList[ii]<1.):
                offset=0.69 
            else: 
                offset=0.59

            x_arr=np.arange(xlim[0],xlim[1]+Nbin,Nbin)
            cut1=np.zeros(len(x_arr),dtype=np.float32)+1.3   
            cut2=x_arr*slope+offset
  
            sel1=x_arr<((1.3-offset)/slope)
            subplot.plot(x_arr[sel1],cut1[sel1],color='blue', linestyle='--', linewidth=2)  
  
            sel2=x_arr>((1.3-offset)/slope) 
            subplot.plot(x_arr[sel2],cut2[sel2],color='blue', linestyle='--', linewidth=2)  
        
        
            #LABELS
            if ii==0:
                plot_label (subplot, 'label', xlim, ylim, 
                            x_percentage=0.15, y_percentage=0.15, color='black', xlog=0, ylog=0, 
                            label='Muzzin 2013', fontsize=13, fontweight='normal') 
                plot_label (subplot, 'line', xlim, ylim,
                            x_percentage=0.04, y_percentage=0.17, color='blue', x2_percentage=0.13, 
                            xlog=0, ylog=0, linestyle='--', linewidth=2)
        
            grid_num+=1
        #endfor
        
    #endfor
    #plt.tight_layout()
    plt.savefig('./fig/plots_UVJ_colour_grid.pdf')
    pdf.savefig()
    plt.close()
#endif UVJ_grid
    


def morphology_vs_stellarmass(G_MR, G_MRII, ThisRedshiftList, pdf):
    
    for ii in range(0,len(ThisRedshiftList)):        
        
        char_redshift="%0.2f" % ThisRedshiftList[ii]
        
        xlim=[8.0,12.]
        ylim=[0.0,1.3]
        bin=0.25
    
        plot_color=['red','purple']        
       
        fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
        subplot=plt.subplot()
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
            
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))
            
        xlab='$\mathrm{log_{10}}(\mathrm{M_{\star}}[M_{\odot}])$'   
        ylab='Fraction'     
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   
       
        #MODEL
        #MR
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)         
        G0_MR=G_MR[sel]                   
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
        BulgeMassRatio=G0_MR['BulgeMass']/G0_MR['StellarMass']
        
        Mass_arr=np.arange(xlim[0],np.amax(StellarMass)+bin/2.,bin)
        BulgeFraction=np.zeros(len(Mass_arr),dtype=np.float32)
        DiskFraction=np.zeros(len(Mass_arr),dtype=np.float32)
        IrrFraction=np.zeros(len(Mass_arr),dtype=np.float32)
                  
        for ll in range(0,len(Mass_arr)):
                sel_bulge=G0_MR[(BulgeMassRatio>0.7) &
                                (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)]
                sel_disk=G0_MR[(BulgeMassRatio<0.7) & (BulgeMassRatio>0.01) &
                               (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)]
                sel_irr=G0_MR[(BulgeMassRatio<0.01) & 
                              (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)]
                #print(len(sel_bulge),len(sel_disk),len(sel_irr))
                if(float(len(sel_bulge))+float(len(sel_disk))+float(len(sel_irr))>0):
                    
                    BulgeFraction[ll]=float(len(sel_bulge))/(float(len(sel_bulge))+
                                                             float(len(sel_disk))+float(len(sel_irr)))
                    DiskFraction[ll]=float(len(sel_disk))/(float(len(sel_bulge))+float(len(sel_disk))+float(len(sel_irr)))
                    IrrFraction[ll]=float(len(sel_irr))/(float(len(sel_bulge))+float(len(sel_disk))+float(len(sel_irr)))
        
        subplot.plot(Mass_arr, BulgeFraction, color='red', linestyle='--', linewidth=2)
        subplot.plot(Mass_arr, DiskFraction, color='blue', linestyle='--', linewidth=2)
        subplot.plot(Mass_arr, IrrFraction, color='green', linestyle='--', linewidth=2)
        
       
        #MRII
        if(MRII==1):
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, ii, FullSnapshotList_MRII)         
            G0_MRII=G_MRII[sel]                   
            StellarMass=stellar_mass_with_err(G0_MRII, Hubble_h, ThisRedshiftList[ii])
            BulgeMassRatio=G0_MRII['BulgeMass']/G0_MRII['StellarMass']
        
            Mass_arr=np.arange(xlim[0],np.amax(StellarMass)+bin/2.,bin)
            BulgeFraction=np.zeros(len(Mass_arr),dtype=np.float32)
            DiskFraction=np.zeros(len(Mass_arr),dtype=np.float32)
            IrrFraction=np.zeros(len(Mass_arr),dtype=np.float32)
                  
            for ll in range(0,len(Mass_arr)):
                sel_bulge=G0_MRII[(BulgeMassRatio>0.7) &
                                (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)]
                sel_disk=G0_MRII[(BulgeMassRatio<0.7) & (BulgeMassRatio>0.01) &
                               (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)]
                sel_irr=G0_MRII[(BulgeMassRatio<0.01) & 
                                (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)]
                #print(len(sel_bulge),len(sel_disk),len(sel_irr))
                if(float(len(sel_bulge))+float(len(sel_disk))+float(len(sel_irr))>0):                     
                    BulgeFraction[ll]=float(len(sel_bulge))/(float(len(sel_bulge))
                                                             +float(len(sel_disk))+float(len(sel_irr)))
                    DiskFraction[ll]=float(len(sel_disk))/(float(len(sel_bulge))+float(len(sel_disk))+float(len(sel_irr)))
                    IrrFraction[ll]=float(len(sel_irr))/(float(len(sel_bulge))+float(len(sel_disk))+float(len(sel_irr)))

            subplot.plot(Mass_arr, BulgeFraction, color='red', linestyle='-', linewidth=2)
            subplot.plot(Mass_arr, DiskFraction, color='blue', linestyle='-', linewidth=2)
            subplot.plot(Mass_arr, IrrFraction, color='green', linestyle='-', linewidth=2)
              
        #OBSERVATIONS
        h=0.7
        file = Datadir + 'conselice2006_bulge_fract.txt'       
        obs = Table.read(file, format='ascii')       
        subplot.errorbar(obs['col1'], obs['col2'], obs['col3'],
                 fmt='o', markersize=5, ecolor='red', color='red')
        file = Datadir + 'conselice2006_disk_fract.txt'       
        obs = Table.read(file, format='ascii')       
        subplot.errorbar(obs['col1'], obs['col2'], obs['col3'],
                 fmt='o', markersize=5, ecolor='blue', color='blue')
        file = Datadir + 'conselice2006_irr_fract.txt'       
        obs = Table.read(file, format='ascii')       
        subplot.errorbar(obs['col1'], obs['col2'], obs['col3'],
                 fmt='o', markersize=5, ecolor='green', color='green')
      
        #MCMC sample
        if opt_plot_MCMC_sample==1:
            file = MCMCSampledir + 'mcmc_plus_obs_BulgeFraction_z'+char_redshift+'.txt'          
            if os.path.isfile(file):
                obs = Table.read(file, format='ascii')      
                subplot.plot(obs['col1']-2.*np.log10(Hubble_h),obs['col4'], color='black', linewidth=2)
                          
                if ii==len(ThisRedshiftList)-1:
                    plot_label (subplot, 'label', xlim, ylim, 
                        x_percentage=0.55, y_percentage=0.85, color='black', xlog=0, ylog=0, 
                        label='MCMC sample', fontsize=13, fontweight='normal') 
                    plot_label (subplot, 'line', xlim, ylim,
                        x_percentage=0.44, y_percentage=0.87, color='black', x2_percentage=0.53, 
                        xlog=0, ylog=0, linestyle='-', linewidth=2)    
            
        #LABELS  
        plot_label (subplot, 'label', xlim, ylim, 
                    x_percentage=0.11, y_percentage=0.91, color='black', xlog=0, ylog=0, 
                    label=prefix_this_model+' - MR', fontsize=13, fontweight='normal')             
        plot_label (subplot,'line',xlim,ylim,x_percentage=0.02,y_percentage=0.93,
                    color='red',x2_percentage=0.09,xlog=0,ylog=0,linestyle='--',linewidth=2)

        plot_label (subplot, 'label', xlim, ylim, 
                    x_percentage=0.11, y_percentage=0.85, color='black', xlog=0, ylog=0, 
                    label=prefix_this_model+' - MRII', fontsize=13, fontweight='normal')             
        plot_label (subplot,'line',xlim,ylim,x_percentage=0.02,y_percentage=0.87,
                    color='red',x2_percentage=0.09,xlog=0,ylog=0,linestyle='-',linewidth=2)
    
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.78, 
                    color='black', xlog=0, ylog=0, label='Conselice 2006', 
                    fontsize=13, fontweight='normal') 
        plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.805, 
                    color='red', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.025) 
     
      
    
    plt.tight_layout()
    plt.savefig('./fig/plots_morphology.pdf')
    plt.savefig('./fig/HYW17_plots_morphology.pdf')
    pdf.savefig()
    plt.close()
   

#end morphology_vs_stellarmass   
    
    
def old_sizes_vs_stellarmass(G_MR, ThisRedshiftList, pdf):
     
    fig = plt.figure(figsize=(two_one_size_small[0],two_one_size_small[1]))
    grid = gridspec.GridSpec(2, 1)
    grid.update(wspace=0.0, hspace=0.0)   
           
    for i_z in range(0,len(ThisRedshiftList)):            
        
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]          
               
        char_redshift="%0.2f" % ThisRedshiftList[i_z] 
        
        xlim=[9.25,11.75]       
        ylim=[-0.5,1.4]     
        bin=0.2    
            
        #i_gal_type==0 -> discs, i_gal_type==0 -> bulges
        for i_gal_type in range (0,2): 
            subplot=plt.subplot(grid[i_gal_type])    
            subplot.set_ylim(ylim), subplot.set_xlim(xlim)            
            #format axis
            majorFormatter = FormatStrFormatter('%d')
            subplot.xaxis.set_major_locator(MultipleLocator(0.5))    
            subplot.xaxis.set_minor_locator(MultipleLocator(.25))  
            subplot.yaxis.set_major_locator(MultipleLocator(0.5))    
            subplot.yaxis.set_minor_locator(MultipleLocator(0.1))  
            xlab='$\mathrm{log_{10}}(M_*[M_{\odot}])$'       
            ylab='$\mathrm{R_{50}}(\mathrm{Kpc})$'            
            if(i_gal_type==0):                 
                plt.tick_params(axis='x', which='both', bottom='on', labelbottom='off')
            else:
                subplot.set_xlabel(xlab, fontsize=14)               
            subplot.set_ylabel(ylab, fontsize=14)
            #subplot.set_yscale('log')
                          
             
            #OBSERVATIONS
            '''if(i_gal_type==0):
                fa=Datadir+'shen2003_discs.txt' 
            else:               
                fa=Datadir+'shen2003_bulges.txt'       
      
            (obs_x,obs_y,obs_y_err_down,obs_y_err_up)=read_data_with_err(fa) 
            log_obs_y_err=[np.log10((obs_y+obs_y_err_up)/obs_y),np.log10(obs_y/(obs_y-obs_y_err_down))]        
            #subplot.errorbar(obs_x,np.log10(obs_y),yerr=log_obs_y_err,fmt='o',markersize=5,ecolor='blue',color='blue')
            #subplot.errorbar(obs_x,obs_y,yerr=[obs_y_err_down,obs_y_err_up],fmt='o',markersize=5,ecolor='blue', color='blue')
            #subplot.plot(obs_x,np.log10(obs_y),color='blue',linewidth=2,linestyle='-')'''
         
            if(i_gal_type==0):
                fa=Datadir+'shen2003_discs_mypoints.txt' 
            else:               
                fa=Datadir+'shen2003_bulges_mypoints.txt'       
            Shen = Table.read(fa, format='ascii')  
           
            subplot.plot(Shen['col1'],np.log10(Shen['col2']),color='blue',linewidth=2,linestyle='-')           
            sigma1,sigma2,M0=0.47,0.34,3.98*1.e10
            scatter=sigma2+(sigma1-sigma2)/(1+(10**Shen['col1']/M0)**2)
            subplot.plot(Shen['col1'],np.log10(Shen['col2'])+scatter/2.,color='blue',linewidth=2,linestyle='--')
            subplot.plot(Shen['col1'],np.log10(Shen['col2'])-scatter/2,color='blue',linewidth=2,linestyle='--')
            
                      
            '''if(i_gal_type==0):
                fa = open(Datadir+"shen2003_discs_with_h.txt", "w")
            else:
                fa = open(Datadir+"shen2003_bulges_with_h.txt", "w")
            fa.write("%d\n" % len(obs_x))
            for kk in range (0,len(obs_x)):
                if(kk<(len(obs_x)-1)):
                    bin=obs_x[kk+1]-obs_x[kk]
                else:
                    bin=obs_x[kk]-obs_x[kk-1]               
                fa.write("%0.2f " % ((obs_x[kk]+np.log10(0.7**2))-bin/2.0) + 
                         "%0.2f " % ((obs_x[kk]+np.log10(0.7**2))+bin/2.0) + 
                         "%0.2f " % (obs_y[kk]*0.7) + 
                         "%0.2f\n" % ((obs_y_err_up[kk]+obs_y_err_down[kk])/2.*0.7))                  
            fa.close()  
            
            if(i_gal_type==0):
                file=Datadir+'shen2003_discs_with_h.txt'  
            else:
                fa = open(Datadir+"shen2003_bulges_with_h.txt", "w")
            f = open(file, 'r')     
            line = int(f.readline())     
            obs = Table.read(file, format='ascii', data_start=1, data_end=line+1)
            obs_xbin=obs['col1']+(obs['col2']-obs['col1'])/2.
            asy_yerror = [np.log10(obs['col3']/(obs['col3']-obs['col4'])), 
                          np.log10((obs['col3']+obs['col4'])/obs['col3'])]
            subplot.errorbar(obs_xbin, np.log10(obs['col3']),yerr=asy_yerror,
                             fmt='o', markersize=5, ecolor='blue', color='blue')'''
        
        
            #MODEL   
            if(i_gal_type==0):
                Gal=G0_MR[ (G0_MR['DiskMass']/G0_MR['StellarMass']>0.8)]  
            else:
                Gal=G0_MR[(G0_MR['BulgeMass']/G0_MR['StellarMass']>0.2)]  
            
            StellarMass=stellar_mass_with_err(Gal, Hubble_h, ThisRedshiftList[i_z])
            #StellarMass-=np.log10(Hubble_h**2)      
            StellarDiskRadius=Gal['StellarHalfMassRadius']*1000./Hubble_h #(from Mpc/h -> Kpc)    

            (x_binned, median, mean, pc16, pc84, rms) = median_and_percentiles (bin, xlim[0], xlim[1], 
                                                                                StellarMass, np.log10(StellarDiskRadius))      
            subplot.plot(x_binned, median,color='red', linewidth=2)           
            subplot.plot(x_binned, median+rms,color='red', linewidth=2, linestyle='--')
            subplot.plot(x_binned, median-rms,color='red', linewidth=2, linestyle='--')
    
            #MCMC sample
            if opt_plot_MCMC_sample==1:
                if(i_gal_type==0):
                    file = MCMCSampledir + 'mcmc_plus_obs_SizevsStellarMass_Discs_z'+char_redshift+'.txt' 
                else:
                    file = MCMCSampledir + 'mcmc_plus_obs_SizevsStellarMass_Bulges_z'+char_redshift+'.txt' 
                if os.path.isfile(file):
                    obs = Table.read(file, format='ascii')      
                    subplot.plot(obs['col1']-2.*np.log10(Hubble_h),np.log10(obs['col4'])-np.log10(Hubble_h), 
                                 color='black', linewidth=2)           
                    if i_z==len(ThisRedshiftList)-1:
                        plot_label (subplot,'label',xlim,ylim,x_percentage=0.55,y_percentage=0.85, 
                                    color='black',xlog=0,ylog=0,label='MCMC sample',fontsize=13,fontweight='normal') 
                        plot_label (subplot, 'line', xlim, ylim,
                                    x_percentage=0.44, y_percentage=0.87, color='black', x2_percentage=0.53, 
                                    xlog=0, ylog=0, linestyle='-', linewidth=2)

            if(i_gal_type==0):                           
                plot_label (subplot, 'label', xlim, ylim, 
                            x_percentage=0.12, y_percentage=0.9, color='black', xlog=0, ylog=0, 
                            label=prefix_this_model, fontsize=13, fontweight='normal')             
                plot_label (subplot,'line',xlim,ylim,x_percentage=0.02,y_percentage=0.92,
                            color='red',x2_percentage=0.10,xlog=0,ylog=0,linestyle='-',linewidth=2)

                plot_label (subplot, 'label', xlim, ylim, 
                            x_percentage=0.12, y_percentage=0.82, color='black', xlog=0, ylog=0, 
                            label='Shen 2003', fontsize=13, fontweight='normal')             
                plot_label (subplot,'line',xlim,ylim,x_percentage=0.02,y_percentage=0.84,
                            color='blue',x2_percentage=0.10,xlog=0,ylog=0,linestyle='-',linewidth=2)
                
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.35, y_percentage=0.1, 
                            color='black', xlog=0, ylog=0, label='disc dominated', 
                            fontsize=15, fontweight='normal')   
            else:
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.35, y_percentage=0.1, 
                            color='black', xlog=0, ylog=0, label='bulge dominated', 
                            fontsize=15, fontweight='normal')   
    
         
    
    plt.tight_layout()
    plt.savefig('./fig/plots_sizes_vs_stellarmass.pdf')
    plt.savefig('./fig/HYW17_plots_sizes_vs_stellarmass.pdf')
    pdf.savefig()
    plt.close() 
        
#end   sizes_vs_stellarmass 


def sizes_vs_stellarmass(G_MR, ThisRedshiftList, pdf):
     
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()    
    
    plot_color=['blue','red']
    
    for i_z in range(0,len(ThisRedshiftList)):            
        
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]          
               
        char_redshift="%0.2f" % ThisRedshiftList[i_z] 
        
        xlim=[9.25,11.75]       
        ylim=[-0.5,1.4]     
        bin=0.2    
            
        #i_gal_type==0 -> discs, i_gal_type==0 -> bulges
        for i_gal_type in range (0,2):            
            subplot.set_ylim(ylim), subplot.set_xlim(xlim)            
            #format axis
            majorFormatter = FormatStrFormatter('%d')
            subplot.xaxis.set_major_locator(MultipleLocator(0.5))    
            subplot.xaxis.set_minor_locator(MultipleLocator(.25))  
            subplot.yaxis.set_major_locator(MultipleLocator(0.5))    
            subplot.yaxis.set_minor_locator(MultipleLocator(0.1))  
            xlab='$\mathrm{log_{10}}(M_*[M_{\odot}])$'       
            ylab='$\mathrm{R_{50}}(\mathrm{Kpc})$'            
            subplot.set_xlabel(xlab, fontsize=14),subplot.set_ylabel(ylab, fontsize=14)
            #subplot.set_yscale('log')
                          
             
            #OBSERVATIONS
            '''if(i_gal_type==0):
                fa=Datadir+'shen2003_discs.txt' 
            else:               
                fa=Datadir+'shen2003_bulges.txt'       
      
            (obs_x,obs_y,obs_y_err_down,obs_y_err_up)=read_data_with_err(fa) 
            log_obs_y_err=[np.log10((obs_y+obs_y_err_up)/obs_y),np.log10(obs_y/(obs_y-obs_y_err_down))]        
            #subplot.errorbar(obs_x,np.log10(obs_y),yerr=log_obs_y_err,fmt='o',markersize=5,ecolor='blue',color='blue')
            #subplot.errorbar(obs_x,obs_y,yerr=[obs_y_err_down,obs_y_err_up],fmt='o',markersize=5,ecolor='blue', color='blue')
            #subplot.plot(obs_x,np.log10(obs_y),color='blue',linewidth=2,linestyle='-')'''
         
            if(i_gal_type==0):
                fa=Datadir+'shen2003_discs_mypoints.txt' 
            else:               
                fa=Datadir+'shen2003_bulges_mypoints.txt'       
            Shen = Table.read(fa, format='ascii')  
           
            #subplot.plot(Shen['col1'],np.log10(Shen['col2']),color=plot_color[i_gal_type],linewidth=2,linestyle=':')  
            sigma1,sigma2,M0=0.47,0.34,3.98*1.e10
            scatter=sigma2+(sigma1-sigma2)/(1+(10**Shen['col1']/M0)**2)
            #subplot.plot(Shen['col1'],np.log10(Shen['col2'])+scatter/2.,color='blue',linewidth=2,linestyle='--')
            #subplot.plot(Shen['col1'],np.log10(Shen['col2'])-scatter/2,color='blue',linewidth=2,linestyle='--')
           
            subplot.fill_between(Shen['col1'],np.log10(Shen['col2'])-scatter/2,
                                 np.log10(Shen['col2'])+scatter/2., facecolor=plot_color[i_gal_type], 
                                 interpolate=True, alpha=0.4, edgecolor=plot_color[i_gal_type])   
            
            '''if(i_gal_type==0):
                fa = open(Datadir+"shen2003_discs_with_h.txt", "w")
            else:
                fa = open(Datadir+"shen2003_bulges_with_h.txt", "w")
            fa.write("%d\n" % len(obs_x))
            for kk in range (0,len(obs_x)):
                if(kk<(len(obs_x)-1)):
                    bin=obs_x[kk+1]-obs_x[kk]
                else:
                    bin=obs_x[kk]-obs_x[kk-1]               
                fa.write("%0.2f " % ((obs_x[kk]+np.log10(0.7**2))-bin/2.0) + 
                         "%0.2f " % ((obs_x[kk]+np.log10(0.7**2))+bin/2.0) + 
                         "%0.2f " % (obs_y[kk]*0.7) + 
                         "%0.2f\n" % ((obs_y_err_up[kk]+obs_y_err_down[kk])/2.*0.7))                  
            fa.close()  
            
            if(i_gal_type==0):
                file=Datadir+'shen2003_discs_with_h.txt'  
            else:
                fa = open(Datadir+"shen2003_bulges_with_h.txt", "w")
            f = open(file, 'r')     
            line = int(f.readline())     
            obs = Table.read(file, format='ascii', data_start=1, data_end=line+1)
            obs_xbin=obs['col1']+(obs['col2']-obs['col1'])/2.
            asy_yerror = [np.log10(obs['col3']/(obs['col3']-obs['col4'])), 
                          np.log10((obs['col3']+obs['col4'])/obs['col3'])]
            subplot.errorbar(obs_xbin, np.log10(obs['col3']),yerr=asy_yerror,
                             fmt='o', markersize=5, ecolor='blue', color='blue')'''
        
        
            #MODEL   
            if(i_gal_type==0):
                Gal=G0_MR[ (G0_MR['DiskMass']/G0_MR['StellarMass']>0.8)]  
            else:
                Gal=G0_MR[(G0_MR['BulgeMass']/G0_MR['StellarMass']>0.2)]  
            
            StellarMass=stellar_mass_with_err(Gal, Hubble_h, ThisRedshiftList[i_z])
            #StellarMass-=np.log10(Hubble_h**2)      
            StellarDiskRadius=Gal['StellarHalfMassRadius']*1000./Hubble_h #(from Mpc/h -> Kpc)    
            #StellarDiskRadius=Gal['StellarDiskRadius']*1000./Hubble_h #(from Mpc/h -> Kpc)    

            (x_binned, median, mean, pc16, pc84, rms) = median_and_percentiles (bin, xlim[0], xlim[1], 
                                                                                StellarMass, np.log10(StellarDiskRadius))      
            subplot.plot(x_binned, median,color=plot_color[i_gal_type], linewidth=2)           
            subplot.plot(x_binned, median+rms,color=plot_color[i_gal_type], linewidth=2, linestyle='--')
            subplot.plot(x_binned, median-rms,color=plot_color[i_gal_type], linewidth=2, linestyle='--')
    
            #MCMC sample
            if opt_plot_MCMC_sample==1:
                if(i_gal_type==0):
                    file = MCMCSampledir + 'mcmc_plus_obs_SizevsStellarMass_Discs_z'+char_redshift+'.txt' 
                else:
                    file = MCMCSampledir + 'mcmc_plus_obs_SizevsStellarMass_Bulges_z'+char_redshift+'.txt' 
                if os.path.isfile(file):
                    obs = Table.read(file, format='ascii')      
                    subplot.plot(obs['col1']-2.*np.log10(Hubble_h),np.log10(obs['col4'])-np.log10(Hubble_h), 
                                 color='black', linewidth=2)           
                    if i_z==len(ThisRedshiftList)-1:
                        plot_label (subplot,'label',xlim,ylim,x_percentage=0.55,y_percentage=0.85, 
                                    color='black',xlog=0,ylog=0,label='MCMC sample',fontsize=13,fontweight='normal') 
                        plot_label (subplot, 'line', xlim, ylim,
                                    x_percentage=0.44, y_percentage=0.87, color='black', x2_percentage=0.53, 
                                    xlog=0, ylog=0, linestyle='-', linewidth=2)

            if(i_gal_type==0):                           
                plot_label (subplot, 'label', xlim, ylim, 
                            x_percentage=0.12, y_percentage=0.9, color='black', xlog=0, ylog=0, 
                            label='Bulge Dominated', fontsize=13, fontweight='normal')             
                plot_label (subplot,'line',xlim,ylim,x_percentage=0.02,y_percentage=0.92,
                            color='red',x2_percentage=0.10,xlog=0,ylog=0,linestyle='-',linewidth=2)

                plot_label (subplot, 'label', xlim, ylim, 
                            x_percentage=0.12, y_percentage=0.82, color='black', xlog=0, ylog=0, 
                            label='Disk Dominated', fontsize=13, fontweight='normal')             
                plot_label (subplot,'line',xlim,ylim,x_percentage=0.02,y_percentage=0.84,
                            color='blue',x2_percentage=0.10,xlog=0,ylog=0,linestyle='-',linewidth=2)
                
                ''' plot_label (subplot, 'label', xlim, ylim, x_percentage=0.35, y_percentage=0.1, 
                            color='black', xlog=0, ylog=0, label='disc dominated', 
                            fontsize=15, fontweight='normal')  ''' 
            '''else:
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.35, y_percentage=0.1, 
                            color='black', xlog=0, ylog=0, label='bulge dominated', 
                            fontsize=15, fontweight='normal') '''  
    
         
    
    plt.tight_layout()
    plt.savefig('./fig/plots_sizes_vs_stellarmass.pdf')
    plt.savefig('./fig/HYW17_plots_sizes_vs_stellarmass.pdf')
    pdf.savefig()
    plt.close() 
        
#end   sizes_vs_stellarmass 
    
    
    
    
    
def ssfr_hist(G_MR, G_MRII, ThisRedshiftList, pdf):
     
    xlim=[-13.5,-8.5]
    ylim=[0.0, 0.35]
    bin_hist=0.1
        
    mass_bin=0.5
    mass_limits=[8.25,11.75]
    mass_bin_arr=np.arange(mass_limits[0],mass_limits[1]+mass_bin,mass_bin)    
        
        
    fig = plt.figure(figsize=(two_four_size_large[0],two_four_size_large[1]))
    grid = gridspec.GridSpec(2, 4)
    grid.update(wspace=0.0, hspace=0.0)
   
    for ii in range(0,len(ThisRedshiftList)):        
         
        #DEFINE MODEL VARIABLES
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]      
        log_StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])               
        log_SSFR=np.log10((G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h)))
        
        #MRII
        if(MRII==1):
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, ii, FullSnapshotList_MRII)
            G0_MRII=G_MRII[sel]  
            log_StellarMass_MRII=stellar_mass_with_err(G0_MRII, Hubble_h, ThisRedshiftList[ii])                       
            log_SSFR_MRII=np.log10((G0_MRII['Sfr']/(G0_MRII['StellarMass']*1.e10/Hubble_h)))
        
        #DEFINE Observational VARIABLES
        file=Datadir+'dr7_gal_final.fit'
        fits_table=fits.open(file)
        dr7_gal_final = fits_table[1]
           
        #SALIM    
        '''fa = open(Datadir+"salim_GSWLC-X1.dat", "r") 
        index=0
        for line in fa:
            index+=1 
        fa.close()
        salim_z=np.zeros(int(index),dtype=np.float32) #8
        salim_log_stellarmass=np.zeros(int(index),dtype=np.float32)#10 
        salim_log_SFR=np.zeros(int(index),dtype=np.float32)#12
        salim_flag_mgs=np.zeros(int(index),dtype=np.float32)#26      
    
        fa = open(Datadir+"salim_GSWLC-X1.dat", "r") 
        index=0
        for line in fa:             
            fields = line.strip().split()             
            salim_z[index]=float(fields[7])
            salim_log_stellarmass[index]=float(fields[9])   
            salim_log_SFR[index]=float(fields[11]) 
            salim_flag_mgs[index]=float(fields[25]) 
            index+=1
        fa.close()'''
      
        '''MUCH SLOWER        
        file = Datadir+"salim_GSWLC-X1.dat"         
        Salim = Table.read(file, format='ascii')
        salim_z=Salim['col8']
        salim_log_stellarmass=Salim['col10']
        salim_log_SFR=Salim['col12']
        salim_flag_mgs=Salim['col26'] '''
            
        for ii in range(0,len(mass_bin_arr)):       
         
            subplot=plt.subplot(grid[ii])                
            subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
            #format axis
            majorFormatter = FormatStrFormatter('%d')
            subplot.xaxis.set_major_locator(MultipleLocator(1))    
            subplot.xaxis.set_minor_locator(MultipleLocator(0.25)) 
            subplot.yaxis.set_major_locator(MultipleLocator(0.1))    
            subplot.yaxis.set_minor_locator(MultipleLocator(0.025))
            xlab,ylab='',''                
            if (ii!=0) & (ii != 4):                
                plt.tick_params(axis='y', which='both', left='on', labelleft='off') 
            else:
                ylab='fraction'  
            if ii<4:
                plt.tick_params(axis='x', which='both', bottom='on', labelbottom='off')      
            else:
                xlab='$\log_{10}(\mathrm{SSFR})[\mathrm{yr}^{-1}]$'
            subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)     
                       
  
            #PLOT MODEL
            #MRII           
            if((mass_bin_arr[ii]<9.5) & MRII==1):
                sel=((log_StellarMass_MRII>mass_bin_arr[ii]-mass_bin/2.) 
                     &(log_StellarMass_MRII<mass_bin_arr[ii]+mass_bin/2.))
                log_SSFR_this_bin=log_SSFR_MRII[sel]
                log_StellarMass_this_bin=log_StellarMass_MRII[sel]
        
                slope,b,width=-0.3,-8.6,0.5            
                sel=log_SSFR_this_bin < -12.      
                log_SSFR_this_bin[sel]=np.log10(np.random.randn(len(log_SSFR_this_bin[sel]))*
                                       width*10**(slope*log_StellarMass_this_bin[sel]+b) + 
                                                10**(slope*log_StellarMass_this_bin[sel]+b))
                            
                bin_arr=np.arange(xlim[0],xlim[1]+bin_hist,bin_hist)
                hist=np.histogram(log_SSFR_this_bin, bins=bin_arr, range=(xlim[0],xlim[1]))            
                x_axis=hist[1][0:len(hist[1][:])-1]+bin_hist/2.           
                y_axis=hist[0]/np.sum(hist[0])     
                subplot.plot(x_axis,y_axis, color='red', linewidth=2, linestyle='-') 
                
            else:                
                sel=(log_StellarMass>mass_bin_arr[ii]-mass_bin/2.) & (log_StellarMass<mass_bin_arr[ii]+mass_bin/2.)
                log_SSFR_this_bin=log_SSFR[sel]
                log_StellarMass_this_bin=log_StellarMass[sel]
        
                slope,b,width=-0.3,-8.6,0.5                       
                sel=log_SSFR_this_bin < -12.0      
                log_SSFR_this_bin[sel]=np.log10(np.random.randn(len(log_SSFR_this_bin[sel]))*
                                       width*10**(slope*log_StellarMass_this_bin[sel]+b) + 
                                                10**(slope*log_StellarMass_this_bin[sel]+b))
                  
                bin_arr=np.arange(xlim[0],xlim[1]+bin_hist,bin_hist)
                hist=np.histogram(log_SSFR_this_bin, bins=bin_arr, range=(xlim[0],xlim[1]))            
                x_axis=hist[1][0:len(hist[1][:])-1]+bin_hist/2.           
                y_axis=hist[0]/np.sum(hist[0])     
                subplot.plot(x_axis,y_axis, color='red', linewidth=2, linestyle='-') 
        
        
       
            #PLOT OBSERVATIONS
            min_z=0.005
            max_z=0.2
            sel=((dr7_gal_final.data['jarle_median_mass']>mass_bin_arr[ii]-mass_bin/2.) & 
                 (dr7_gal_final.data['jarle_median_mass']<mass_bin_arr[ii]+mass_bin/2.) &
                 (dr7_gal_final.data['z'] > min_z) & (dr7_gal_final.data['z'] < max_z))                
            dr7_log_SSFR=dr7_gal_final.data['median_ssfr'][sel]        
            dr7_mag=dr7_gal_final.data['jarle_mag'][sel,2]
          
            (x_arr,frac)=plot_fraction_vmax_weighted(xlim[0],xlim[1],bin_hist,dr7_log_SSFR,dr7_mag)           
            subplot.plot(x_arr,frac, color='black', linewidth=2, linestyle='-') 
        
            #Salim
            '''min_z=0.005
            max_z=0.2
            sel=((salim_log_stellarmass>mass_bin_arr[ii]-mass_bin/2.) & 
                 (salim_log_stellarmass<mass_bin_arr[ii]+mass_bin/2.) &
                 (salim_z > min_z) & (salim_z < max_z) & (salim_flag_mgs==1))                
            salim_log_SSFR_this_bin=salim_log_SFR[sel]-salim_log_stellarmass[sel]
            
            bin_arr=np.arange(xlim[0],xlim[1]+bin_hist,bin_hist)
            hist=np.histogram(salim_log_SSFR_this_bin, bins=bin_arr, range=(xlim[0],xlim[1]))            
            x_axis=hist[1][0:len(hist[1][:])-1]+bin_hist/2.           
            y_axis=hist[0]/np.sum(hist[0])     
            subplot.plot(x_axis,y_axis, color='blue', linewidth=2, linestyle='-') '''
        
            if(ii==0):
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.12, y_percentage=0.79, 
                            color='black', xlog=0, ylog=0, label=prefix_this_model, fontsize=12, 
                            fontweight='normal')             
                plot_label (subplot, 'line', xlim, ylim, x_percentage=0.025, y_percentage=0.81, 
                            color='red', x2_percentage=0.1, xlog=0, ylog=0, linestyle='-', linewidth=2)
                
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.12, y_percentage=0.69, 
                            color='black', xlog=0, ylog=0, label='SDSS/DR7', fontsize=12, 
                            fontweight='normal')             
                plot_label (subplot, 'line', xlim, ylim, x_percentage=0.025, y_percentage=0.71, 
                            color='black', x2_percentage=0.1, xlog=0, ylog=0, linestyle='-', linewidth=2)
        
            #LABELS
            mass_low, mass_high=mass_bin_arr[ii]-mass_bin/2.,mass_bin_arr[ii]+mass_bin/2.
            label="%0.2f" % mass_low + r'$<\mathrm{log_{10}}(M_*[M_{\odot}])<$' + "%0.2f" % mass_high
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.88, 
                        color='black', xlog=0, ylog=0, label=label, fontsize=11, fontweight='normal') 
        
        
        #endfor - loop on mass bins
        
        
        fits_table.close()
        
    plt.tight_layout()
    plt.savefig('./fig/plots_ssfr_hist.pdf')
    plt.savefig('./fig/HYW17_plots_ssfr_hist.pdf')
    pdf.savefig()
    plt.close()

    
#end ssfr_hist
   
    
    
def SFH(G_MR, ThisRedshiftList, pdf):
     
    xlim=[0.,12.]
    ylim=[0.0, 1.5]
    bin_hist=0.1
        
    mass_bin=0.5
    mass_limits=[9.0,11.5]
    mass_bin_arr=np.arange(mass_limits[0],mass_limits[1]+mass_bin,mass_bin)    
    
    obs_z=['0.35','0.60','0.85','1.15','1.50','1.90']
    obs_LookBackTime=[2.42, 5.02, 6.29, 7.7, 8.81, 9.69]
    #obs_z=['0.60','0.85','1.15','1.50','1.90']
    #obs_LookBackTime=[5.02, 6.29, 7.7, 8.81, 9.69]
    
    color_plot=['purple','blue','green','yellow','orange','red']     
   
    fig = plt.figure(figsize=(7,6))
   
   
    for ii in range(0,len(ThisRedshiftList)):        
         
        #DEFINE MODEL VARIABLES
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]             
       
        subplot=plt.subplot()                
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
        ylab='fraction'          
        xlab='$Gyr$'
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   
        
        #Read              
        fa = open(DirName_MR+"SFH_Bins","rb")                
        nbins =  np.fromfile(fa,np.int32,1)
        template = np.dtype([('SnapNum',np.int32,1),
                             ('Bin',np.int32,1),
                             ('Lookbacktime',np.float64,1),                           
                             ('dt',np.float64,1),
                             ('nbins',np.int32,1)
                            ])
        SFH = np.fromfile(fa,template,nbins)    
        fa.close()            
        #we only need the SFH strucutre from the current snap
        SFH=SFH[SFH['SnapNum']==G0_MR['SnapNum'][0]]
        SFH['Lookbacktime']=SFH['Lookbacktime']/1.e9+obs_LookBackTime[ii]
        
        G0_MR_SFH=(G0_MR['sfh_DiskMass'][:,0:len(SFH)]+G0_MR['sfh_BulgeMass'][:,0:len(SFH)])
        log_StellarMass=(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)) 
        log_SSFR=(np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h))) 
        color_UV=G0_MR['MagDust'][:,0]-G0_MR['MagDust'][:,2]  
        color_VJ=G0_MR['MagDust'][:,2]-G0_MR['MagDust'][:,7]       
       
        for ll in range(0,len(mass_bin_arr)): 
            
            #MODEL
            sel=((log_StellarMass>mass_bin_arr[ll]-mass_bin/2.) & 
                 (log_StellarMass<mass_bin_arr[ll]+mass_bin/2.) &                             
                 (((color_VJ < (minimum_y_color_cut[1]-offset_color_cut[1])/slope_color_cut[1]) &
                               (color_UV > minimum_y_color_cut[1])) |
                  ((color_VJ > (minimum_y_color_cut[1]-offset_color_cut[1])/slope_color_cut[1]) &
                               (color_UV > (color_VJ*slope_color_cut[1] + offset_color_cut[1]))))              
                 ) 
                     
                        
                              
                
            
            G0_MR_SFH_this_bin=G0_MR_SFH[sel]
            y_arr=np.zeros(len(SFH),dtype=np.float32)
           
            for jj in range(0,len(SFH)):             
                y_arr[jj]=np.sum(G0_MR_SFH_this_bin[:,jj])/len(G0_MR_SFH_this_bin[:,jj])/SFH['dt'][jj]
            y_arr=y_arr/np.sum(y_arr)    
            #subplot.plot(SFH['Lookbacktime']/1.e9,y_arr/np.max(y_arr), color=color_plot[ll], linewidth=2, linestyle=':')
            
            for jj in range (0,2):
                for kk in range(1, len(y_arr)-1):
                    y_arr[kk]=(y_arr[kk-1]+y_arr[kk]+y_arr[kk+1])/3.
            subplot.plot(SFH['Lookbacktime'],y_arr/np.max(y_arr),color=color_plot[ll], linewidth=2, linestyle='-') 
          
            #OBSERVATIONS
            #Read 
            char_mass="%0.1f" % mass_bin_arr[ll]
            fa = open(Datadir+"/pacifici_sfh/sfhmed_"+char_mass+"_z"+obs_z[ii]+".dat", "r") 
            nbins=0
            for line in fa: 
                nbins+=1           
            fa.close()
               
            x_axis_obs=np.zeros(int(nbins),dtype=np.float32)
            y_axis_p25_obs=np.zeros(int(nbins),dtype=np.float32)  
            y_axis_p50_obs=np.zeros(int(nbins),dtype=np.float32)
            y_axis_p75_obs=np.zeros(int(nbins),dtype=np.float32)
            fa = open(Datadir+"/pacifici_sfh/sfhmed_"+char_mass+"_z"+obs_z[ii]+".dat", "r") 
            index=0
            for line in fa:            
                fields = line.strip().split()               
                x_axis_obs[index]=float(fields[0])
                y_axis_p25_obs[index]=float(fields[1]) 
                y_axis_p50_obs[index]=float(fields[2])
                y_axis_p75_obs[index]=float(fields[3])                
                index+=1 
                
            x_axis_obs+=obs_LookBackTime[ii]
            n_obs_bins=len(y_axis_p50_obs[y_axis_p50_obs>0.])
            #subplot.plot(x_axis_obs,y_axis_p50_obs/np.sum(y_axis_p50_obs)/n_obs_bins, color=color_plot[ll], linewidth=2, linestyle='-') 
            subplot.plot(x_axis_obs,y_axis_p50_obs/np.amax(y_axis_p50_obs), color=color_plot[ll], linewidth=2, linestyle='--') 
            fa.close()
                
        '''#AVERAGE FOR ALL MASSES    
        sel=(log_StellarMass>9.0) & (log_StellarMass<13.0)  
        G0_MR_SFH_this_bin=G0_MR_SFH[sel]
        y_arr=np.zeros(len(SFH),dtype=np.float32)
        for jj in range(0,len(log_StellarMass[sel])):
            y_arr+=G0_MR_SFH_this_bin[jj,:]
            
        y_arr=y_arr/len(log_StellarMass[sel]) 
        y_arr=y_arr/np.sum(y_arr)               
        subplot.plot(SFH['Lookbacktime']/1.e9,y_arr/np.max(y_arr), color='black', linewidth=2, linestyle='-')'''
      

        
    plt.tight_layout()
    plt.savefig('./fig/plots_SFH.pdf')
    pdf.savefig()
    plt.close()

    
    plt.rcParams.update({'font.size': 18, 'xtick.labelsize': 18, 'ytick.labelsize': 18})
#end SFH
    
    
    
    
def cooling_heating(G_MR, ThisRedshiftList, pdf):
          
    
      
    '''   
    fig = plt.figure(figsize=(7,7))
    xlim=[9.0,11.5]
    ylim=[4.0,9.0]
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim) 
    
    for ii in range(0,len(ThisRedshiftList)):        
           
        #SCATTER 
        if(ii==1): 
            (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
            G0_MR_unsel=G_MR[sel]           
            G0_MR=G0_MR_unsel[(np.log10(G0_MR_unsel['BlackHoleMass']*1.e10) > 0.) & (G0_MR_unsel['Type'] == 0)]   
        
            x_axis=np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)
            CoolingGas=G0_MR['CoolingGas']*1.e10/Hubble_h
            CoolingGas[CoolingGas==0.]=1.e5
            y_axis=np.log10(CoolingGas)
          
            subplot.scatter(x_axis, y_axis, s=5, color='blue')      
        
        
    plt.tight_layout()  
    pdf.savefig()
    plt.close()'''    
        
        
    SSFR_cut=[-11.,-11.0, -10.5,-10.]
   
    xlim=[9.0,12.0]
    ylim=[-2.0, 3.]
    #ylim=[-0.5, 0.5]
           
    #plot_color=['gold','orange','limegreen','green']        
    #cmap = plt.get_cmap('Greens')  
    #plot_color = [cmap(i) for i in np.linspace(0.5, 0.8, 3)]
    linewidth=[4,2,2,2]
    linestyle=['-','-','-','-']
    
    fig = plt.figure(figsize=(7,6))
  
    #constants
    UnitLength_in_cm=3.08568e+24
    UnitVelocity_in_cm_per_s=100000.
    UnitTime_in_s=UnitLength_in_cm / UnitVelocity_in_cm_per_s
    SOLAR_MASS=1.989e33
    UnitMass_in_g=1.989e+43
    SEC_PER_YEAR=3.155e7
    UnitEnergy_in_cgs = UnitMass_in_g * UnitLength_in_cm**2. / UnitTime_in_s**2.

    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)  

    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 

    xlab='$\mathrm{log_{10}}(M_{\star}[M_{\odot}])$'       
    ylab='$log_{10}(\mathrm{AGN\,Heating\,Rate}/\mathrm{Cooling\,Rate})$' 
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   

    #heating/cooling balance line    
    x_arr=np.arange(xlim[0],xlim[1]+0.05,0.05)     
    y_arr=x_arr*0.       
    subplot.plot(x_arr,y_arr,color='black', linestyle='--', linewidth=2)     
    
    for ii in range(0,len(ThisRedshiftList)):             
        #SCATTER 
        if(ii==0): 
            (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
            G0_MR_unsel=G_MR[sel]   
            #G0_MR=G0_MR_unsel
            '''G0_MR=G0_MR_unsel[(np.log10(G0_MR_unsel['StellarMass']*1.e10*Hubble_h) > 8.) & 
                              (np.log10(G0_MR_unsel['BlackHoleMass']*1.e10) > 4.) &
                              (G0_MR_unsel['Type'] == 0)]'''
            G0_MR=G0_MR_unsel[(np.log10(G0_MR_unsel['BlackHoleMass']*1.e10) > 0.) & (G0_MR_unsel['Type'] == 0)]

            
            #sel=G0_MR['HaloID']==5000000000190
            #G0_MR=G0_MR[sel]
            
            G0_MR=np.random.choice(G0_MR, size=2000)   
            log_StellarMass=(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h))               
            log_SSFR=np.log10((G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h)))

            Cooling=G0_MR['CoolingRate_beforeAGN'] * (UnitTime_in_s/SEC_PER_YEAR)*(SOLAR_MASS/UnitMass_in_g)
            AgnEfficiency=5.3e-3*(UnitTime_in_s*SOLAR_MASS)/(UnitMass_in_g*SEC_PER_YEAR)
            AGNcoeff = (1.34e5 / G0_MR['Vvir']) * (1.34e5 / G0_MR['Vvir'])           
            AGNrate= AgnEfficiency * G0_MR['BlackHoleMass']/Hubble_h* (G0_MR['HotGas']/Hubble_h) * 10.         
            EDDrate = 1.3e48 * G0_MR['BlackHoleMass'] / (UnitEnergy_in_cgs / UnitTime_in_s) / 9e10;           
            sel=AGNrate>EDDrate
            AGNrate[sel]=EDDrate[sel]
            #dt=1.616664e-05
            AGNaccreted=AGNrate          
            AGNheating = AGNcoeff * AGNaccreted 
            sel=AGNheating > G0_MR['CoolingRate_beforeAGN']          
            #AGNheating[sel] = G0_MR[sel]['CoolingRate_beforeAGN']
         
            #Cooling[Cooling==0.]=0.0000001
            x_axis=log_StellarMass
            y_axis=np.log10(AGNheating/Cooling)
            log_SSFR=log_SSFR
            
            '''Cooling=G0_MR['CoolingRate_beforeAGN']
            Heating=G0_MR['CoolingRate_beforeAGN']-G0_MR['CoolingRate']
         
            sel=Cooling>0.
            x_axis=log_StellarMass[sel]
            y_axis=np.log10(Heating[sel]/Cooling[sel])
            log_SSFR=log_SSFR[sel]'''
                
            subplot.scatter(x_axis, y_axis, s=5, color='blue')           
            sel=log_SSFR<SSFR_cut[ii]
            subplot.scatter(x_axis[sel], y_axis[sel], s=5, color='red') 

            
            #SECOND AXIS WITH MVIR  
            G0_MR=G0_MR_unsel[(np.log10(G0_MR_unsel['StellarMass']*1.e10) > 0.) & (G0_MR_unsel['Type'] == 0)]
            StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
            Mvir=np.log10(G0_MR['Mvir']*1.e10/Hubble_h)
            (x_binned,median_MR,mean,pc16,pc84,rms)=median_and_percentiles(1.0,xlim[0],xlim[1],StellarMass,Mvir) 
            y_axis=median_MR
       
            for jj in range(0,len(y_axis)):
                y_axis[jj] = int((y_axis[jj] * 10) + 0.5) / 10.0            
             
            ax2 = subplot.twiny()        
            ax2.set_xbound(subplot.get_xbound())        
            ax2.set_xticks(x_binned)
            ax2.set_xticklabels(y_axis)

            xlab='$\mathrm{log_{10}}(<M_{\mathrm{vir}}>[M_{\odot}])$'  
            ax2.set_xlabel(xlab, fontsize=14)
            ax2.xaxis.set_label_position('top')  
            
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.065, y_percentage=0.84, 
                        color='red', xlog=0, ylog=0, label='Passive (z=0)', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.85, 
                        color='red', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.) 

            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.065, y_percentage=0.8, 
                        color='blue', xlog=0, ylog=0, label='Star Forming (z=0)', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.81, 
                        color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.)


        #MEDIAN
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR_unsel=G_MR[sel]   
        #G0_MR=G0_MR_unsel[(np.log10(G0_MR_unsel['BlackHoleMass']*1.e10) > 0.) & (G0_MR_unsel['Type'] == 0)]
        G0_MR=G0_MR_unsel[(np.log10(G0_MR_unsel['StellarMass']*1.e10/Hubble_h) > xlim[0]-2.) & 
                          (np.log10(G0_MR_unsel['BlackHoleMass']*1.e10/Hubble_h) > ylim[0]-2.) &
                          (G0_MR_unsel['Type'] == 0)]
                         
        log_StellarMass=(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)) 
        Cooling=G0_MR['CoolingRate_beforeAGN'] * (UnitTime_in_s/SEC_PER_YEAR)*(SOLAR_MASS/UnitMass_in_g)
        AgnEfficiency=5.3e-3*(UnitTime_in_s*SOLAR_MASS)/(UnitMass_in_g*SEC_PER_YEAR)
        AGNcoeff = (1.34e5 / G0_MR['Vvir']) * (1.34e5 / G0_MR['Vvir'])           
        AGNrate= AgnEfficiency * G0_MR['BlackHoleMass']/Hubble_h* (G0_MR['HotGas']/Hubble_h) * 10.         
        EDDrate = 1.3e48 * G0_MR['BlackHoleMass'] / (UnitEnergy_in_cgs / UnitTime_in_s) / 9e10;   
        
        sel=AGNrate>EDDrate
        AGNrate[sel]=EDDrate[sel]              
        AGNheating = AGNcoeff * AGNrate  
           
        sel=((Cooling>0.) & (AGNheating>0.))        
        x_axis=log_StellarMass[sel]
        y_axis=np.log10(AGNheating[sel]/Cooling[sel])
            
            
        '''bin=0.15
        (x_binned_1, median_1, mean, pc16, pc84, rms)=median_and_percentiles (bin, ylim[0], ylim[1], y_axis, x_axis)        
        (x_binned_2, median_2, mean, pc16, pc84, rms)=median_and_percentiles (bin, xlim[0], xlim[1], x_axis, y_axis)  
        x=np.append(x_binned_2[x_binned_2<10.4],median_1[median_1>10.4])
        y=np.append(median_2[x_binned_2<10.4],x_binned_1[median_1>10.4]) 
            
        #smooth the bining  
        for jj in range (0,3):
            for kk in range(1, len(y)-1):
                y[kk]=(y[kk-1]+y[kk]+y[kk+1])/3.

        for jj in range (0,3):
            for kk in range(1, len(x)-1):
                x[kk]=(x[kk-1]+x[kk]+x[kk+1])/3.

        if(ii==0):
            sel=x<11.
        else:
            if (ii==1):
                sel=x<10.9
            else:
                sel=x<10.7
                  
        subplot.plot(x[sel], y[sel], color='darkorange', linewidth=linewidth[ii],linestyle=linestyle[ii])          
        subplot.plot(x_binned_2[sel], median_2[sel], color='darkorange', linewidth=linewidth[ii],linestyle=linestyle[ii])'''
        
        
        
        bin=[0.1,0.1]
        Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]
        H, xedges, yedges = np.histogram2d(x_axis, y_axis, bins=Nbins, range=[[xlim[0]-bin[0]/2., xlim[1]+bin[0]/2.], 
                                                                              [ylim[0]-bin[1]/2., ylim[1]+bin[1]/2.]]) 
        extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
        plt.subplots_adjust(bottom=0.15, left=0.15)     
        mylevels = np.log10(np.logspace(-2.0, 0.0, num=10)) 
        HH=H
        H = zoom(H, 20)        
        H=np.log10(H/np.amax(H))       
        #cont=plt.contourf(H.transpose()[::], origin='lower', cmap='Greys_r', levels=mylevels, extent=extent)  
                           
        xx=np.arange(xlim[0],xlim[1],bin[0])   
        yy=np.arange(ylim[0],ylim[1],bin[1])          
        mean2d_arr=np.zeros(len(HH[:,0]))
           
        for jj in range(0,len(HH[:,0])): 
            sel=(HH[jj,:]==np.amax(HH[jj,:]))               
            yy_aux=yy[sel]
            mean2d_arr[jj]=np.median(yy_aux)
                
        for jj in range (0,5):
            for kk in range(1, len(mean2d_arr)-1):
                mean2d_arr[kk]=(mean2d_arr[kk-1]+mean2d_arr[kk]+mean2d_arr[kk+1])/3.
         
        #sel=((xx<11.) & (xx>9.5))
        sel=((xx<11.) & (xx>9.0))
        subplot.plot(xx[sel], mean2d_arr[sel], color='darkorange', linewidth=linewidth[ii],linestyle=linestyle[ii])
         
               

        
        #x_per=0.75
        #y_per=[0.69,0.75,0.81]
        #x_per=0.75
        #y_per=[0.59,0.65,0.71]
        x_per=0.8
        y_per=[0.43,0.5,0.6,0.7]
        #y_per=[0.6,0.67,0.75,0.82]
        
        if ii==0:
            plot_label (subplot, 'label', xlim, ylim, x_percentage=x_per, y_percentage=y_per[0], 
                        color='darkorange', xlog=0, ylog=0, label='z=3', 
                        fontsize=13, fontweight='bold',backgroundcolor='white',back_alpha=0.9)            
        if ii==1:
            plot_label (subplot, 'label', xlim, ylim, x_percentage=x_per, y_percentage=y_per[1], 
                        color='darkorange', xlog=0, ylog=0, label='z=2', 
                        fontsize=13, fontweight='bold',backgroundcolor='white',back_alpha=0.9)            
        if ii==2:
            plot_label (subplot, 'label', xlim, ylim, x_percentage=x_per, y_percentage=y_per[2], 
                        color='darkorange', xlog=0, ylog=0, label='z=1', 
                        fontsize=13, fontweight='bold',backgroundcolor='white',back_alpha=0.9)    
                 
        if ii==3:
            plot_label (subplot, 'label', xlim, ylim, x_percentage=x_per, y_percentage=y_per[3], 
                        color='darkorange', xlog=0, ylog=0, label='z=0', 
                        fontsize=13, fontweight='bold',backgroundcolor='white',back_alpha=0.9)    
                
                
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.06, y_percentage=0.055, 
                    color='black', xlog=0, ylog=0, label='Median at different redshifts', 
                    fontsize=10, fontweight='normal') 
        
        plot_label (subplot, 'line', xlim, ylim, x_percentage=0.05, y_percentage=0.065, 
                    color='darkorange', x2_percentage=0.01, xlog=0, ylog=0, linestyle='-', linewidth=2)  
        
                
    plt.tight_layout()
    plt.savefig('./fig/plots_cooling_heating.pdf')
    pdf.savefig()
    plt.close()

    
    plt.rcParams.update({'font.size': 18, 'xtick.labelsize': 18, 'ytick.labelsize': 18})
#end cooling_heating

    
    
    
    
    
    
    
    

    
    
    
    
    
    
    
    
#def BHBM_by_sfr(G_MR, G_MRII, ThisRedshiftList, pdf):    
def BHBM_by_sfr(G_MR, ThisRedshiftList, pdf):
       
    plot_inset=1
    other_models=0
    
    '''xlim=[9.0,11.5]
    ylim=[-13.0, -7.]   
         
    fig = plt.figure(figsize=(10,8)) 
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
   
    xlab='$\mathrm{log_{10}}(M_{\star}[M_{\odot}])$'       
    ylab='$\mathrm{log_{10}}(M_{\mathrm{BH}}[M_{\odot}])$' 
    subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)   
    
    plot_color=['red','green','blue']    
        
    for ii in range(0,len(ThisRedshiftList)):        
                    
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]  
        log_StellarMass=np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)
        log_SSFR=np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h))     
        subplot.scatter(log_StellarMass, log_SSFR,color=plot_color[ii])
    
    
    plt.tight_layout()  
    pdf.savefig()
    plt.close()'''
    
    
    '''xlim=[9.0,11.5]
    ylim=[-2.0, 3.]   
         
    fig = plt.figure(figsize=(10,8)) 
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
   
    xlab='$\mathrm{log_{10}}(M_{\star}[M_{\odot}])$'       
    ylab='$\mathrm{log_{10}}(M_{\mathrm{BH}}[M_{\odot}])$' 
    subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)   
    
    plot_color=['red','green','blue']    
        
    for ii in range(0,1):        
                    
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]       
        G0_MR=np.random.choice(G0_MR, size=1000.)   
        
        log_StellarMass=np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)      
        a=np.log10(G0_MR['CoolingRate_beforeAGN'])    
        a[a==0.]=0.00001
        subplot.scatter(log_StellarMass, a,color='blue')        
        b=np.log10(G0_MR['CoolingRate'])  
        b[b==0.]=0.00001
        subplot.scatter(log_StellarMass, b,color=plot_color[ii])
      
    
    plt.tight_layout()  
    pdf.savefig()
    plt.close()'''
    
    
    
    #local_slope_color_cut=[0.075,0.3, 0.32]
    #local_offset_color_cut=[1.85,1.18,0.99]
    #NOAGN
    #local_slope_color_cut=[0.5,0.3, 0.32]
    #local_offset_color_cut=[1.085,1.18,0.99]
    #local_slope_color_cut=[0.075,0.48, 0.38]
    #local_offset_color_cut=[2.05,1.0,1.0]
    
    SSFR_cut=[-11.,-10.5,-10.]
   
    xlim=[9.0,11.5]
    ylim=[4.5, 10.]
    if(plot_inset==0):
        ylim=[4.5, 9.]
        
    inset_xlim=[9.0,11.5]
    inset_ylim=[-6.0, -1.]    
        
    plot_color=['red','orange','green','blue']        
 
    fig = plt.figure(figsize=(7,6)) 
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
      
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
    
    xlab='$\mathrm{log_{10}}(M_{\star}[M_{\odot}])$'       
    ylab='$\mathrm{log_{10}}(M_{\mathrm{BH}}[M_{\odot}])$' 
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   

    
    if(plot_inset==1):
        # this is an inset axes over the main axes
        plt.rcParams.update({'font.size': 10, 'xtick.labelsize': 10, 'ytick.labelsize': 10,
                             'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 
                             'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0})
        #inset_xlim=[4.,9.]
        #inset_ylim=[.0, .5]   
        inset = inset_axes(subplot, width="100%", height="100%", 
                           bbox_to_anchor=(0.16, 0.65, 0.38, 0.3), 
                           bbox_transform=subplot.figure.transFigure)
        #fig.subplots_adjust(hspace=0.4)
        inset.set_ylim(inset_ylim), inset.set_xlim(inset_xlim)    

        inset.xaxis.set_major_locator(MultipleLocator(1))    
        inset.xaxis.set_minor_locator(MultipleLocator(0.5))    
        inset.yaxis.set_label_position("right")
        xlab='$\mathrm{log_{10}}(M_{\star}[M_{\odot}])$'       
        ylab='$\mathrm{log_{10}}(M_{\mathrm{BH}}/M_{\mathrm{vir}}^{0.097}\mathrm{x}(1+z)^{-1.5})$' 
        inset.set_xlabel(xlab, fontsize=10), inset.set_ylabel(ylab, fontsize=10, rotation=270, labelpad=20)   
        inset.tick_params(axis='y', which='both', left='on', labelleft='off', right='on',labelright='on')
       
    for ii in range(0,len(ThisRedshiftList)):        
                    
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR_unsel=G_MR[sel]        
        G0_MR=G0_MR_unsel[(np.log10(G0_MR_unsel['StellarMass']*1.e10*Hubble_h) > xlim[0]-1.) & 
                          (np.log10(G0_MR_unsel['BlackHoleMass']*1.e10) > ylim[0]-1.) &
                          (G0_MR_unsel['Type'] == 0)]
       
        #median
        bin=0.25 
        log_StellarMass=(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)) 
        log_BHMass=(np.log10(G0_MR['BlackHoleMass']*1.e10/Hubble_h))     
        inset_y_axis=np.log10(G0_MR['BlackHoleMass']/(G0_MR['Mvir']**0.097)/((1+ThisRedshiftList[ii])**1.5))
           
        (x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles (bin, ylim[0], ylim[1],log_BHMass, log_StellarMass)
        #(x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles (bin, xlim[0], xlim[1], log_StellarMass,log_BHMass)
        (inset_x_binned, inset_median, inset_mean, inset_pc16,inset_pc84, rms) = median_and_percentiles(bin,inset_ylim[0],inset_ylim[1],inset_y_axis,log_StellarMass) 
        
        G0_MR=np.random.choice(G0_MR, size=5000)        
        log_StellarMass=(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)) 
        log_BHMass=(np.log10(G0_MR['BlackHoleMass']*1.e10/Hubble_h)) 
        inset_y_axis=np.log10(G0_MR['BlackHoleMass']/(G0_MR['Mvir']**0.097)/((1+ThisRedshiftList[ii])**1.5))
        log_SSFR=np.log10((G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h)))
        
      
        color_ur=G0_MR['MagDust'][:,15]-G0_MR['MagDust'][:,17]  
        Magr=G0_MR['MagDust'][:,17]-5.*np.log10(Hubble_h)
        color_UV=G0_MR['MagDust'][:,0]-G0_MR['MagDust'][:,2]  
        color_VJ=G0_MR['MagDust'][:,2]-G0_MR['MagDust'][:,7]   
        
                        
         
       
        if(ii==2):    
            '''bin=[0.25,0.25]
            Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]
            (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)
            aux_G0_MR_unsel=G_MR[sel]           
            aux_G0_MR=aux_G0_MR_unsel[(aux_G0_MR_unsel['StellarMass'] > 0.) & (aux_G0_MR_unsel['BlackHoleMass'] > 0.) &
                                      (aux_G0_MR_unsel['Type'] == 0)]
     
            aux_log_StellarMass=(np.log10(aux_G0_MR['StellarMass']*1.e10*Hubble_h)) 
            aux_log_BHMass=(np.log10(aux_G0_MR['BlackHoleMass']*1.e10)) 
            H, xedges, yedges = np.histogram2d(aux_log_StellarMass, aux_log_BHMass, bins=Nbins,
                                               range=[[xlim[0]-bin[0]/2., xlim[1]+bin[0]/2.], 
                                                      [ylim[0]-bin[1]/2., ylim[1]+bin[1]/2.]])            
            extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
            plt.subplots_adjust(bottom=0.15, left=0.15)     
            mylevels = np.log10(np.logspace(-2.0, 0.0, num=10))  
            H = zoom(H, 20)        
            H=np.log10(H/np.amax(H))       
            cont=subplot.contourf(H.transpose()[::], origin='lower', cmap='Greys_r', levels=mylevels, extent=extent)'''
            
            subplot.scatter(log_StellarMass, log_BHMass, s=5, color='blue') 
            #sel=(color_ur>(local_offset_color_cut[ii]-local_slope_color_cut[ii]*np.tanh((Magr+20.07)/1.09)))
            #sel=(color_ur>(local_offset_color_cut[ii]-local_slope_color_cut[ii]*np.tanh((Magr+18.07)/1.09)))
            #sel=((color_UV > minimum_y_color_cut[ii]) & 
            #     (color_UV > (color_VJ*local_slope_color_cut[ii] + local_offset_color_cut[ii])))
            sel=log_SSFR<SSFR_cut[ii]
            subplot.scatter(log_StellarMass[sel], log_BHMass[sel], s=5, color='red') 
            
        if(ii==0): 
            x_arr=np.arange(9.78,10.02,0.005)  
            #x_arr=np.arange(9.74,9.95,0.005)  
            #z=0
            (slope0,b0)=get_slope(x1=9.,y1=5.8,x2=11.5,y2=6.25)
            #z=2
            (slope2,b2)=get_slope(x1=9.,y1=6.6,x2=11.5,y2=7.05)
            if(other_models==0):
                subplot.fill_between(x_arr-2.*np.log10(Hubble_h),x_arr*slope0+b0,x_arr*slope2+b2, facecolor='lightgrey', 
                                     interpolate=True, alpha=0.4, edgecolor='black')   
                      
        if(ii==0):                    
            (slope,b)=get_slope(x1=9.,y1=5.8,x2=11.5,y2=6.5)
        else:
            if(ii==1):
                (slope,b)=get_slope(x1=9.,y1=6.25,x2=11.5,y2=6.95)
            else:
                (slope,b)=get_slope(x1=9.,y1=6.6,x2=11.5,y2=7.3)
            
            
        x_arr=np.arange(xlim[0]+2.*np.log10(Hubble_h),xlim[1]+0.05,0.05)     
        y_arr=x_arr*slope+b    
        if(other_models==0):
            subplot.plot(x_arr-2.*np.log10(Hubble_h),y_arr,color=plot_color[ii], linestyle='--', linewidth=2)         
                    
        #median        
        #sel=(x_binned>4.7) & (x_binned<8.5) &(median>xlim[0])  
        sel=(x_binned>5.2) & (x_binned<8.5) &(median>xlim[0])  
        subplot.plot(median[sel],x_binned[sel], color=plot_color[ii], linewidth=2,linestyle='-')
        #subplot.plot(x_binned, median,color=plot_color[ii], linewidth=2,linestyle='-')
      
        #inset 
        if(plot_inset==1):
            if(ii==0):       
                sel=log_SSFR<SSFR_cut[ii]
                inset.scatter(log_StellarMass[sel], inset_y_axis[sel], s=2.5, color='red') 
                sel=log_SSFR>SSFR_cut[ii]
                inset.scatter(log_StellarMass[sel], inset_y_axis[sel], s=2.5, color='blue')   
                
            if(ii==0): 
                x_arr=np.arange(9.8,10.22,0.005)  
                #z=0
                (slope0,b0)=get_slope(x1=9.,y1=-4.15,x2=11.,y2=-4.15)
                #z=2
                (slope2,b2)=get_slope(x1=9.,y1=-4.25,x2=11.,y2=-4.25)
                inset.fill_between(x_arr-2.*np.log10(Hubble_h),x_arr*slope0+b0,x_arr*slope2+b2, facecolor='lightgrey', 
                                   interpolate=True, alpha=0.4, edgecolor='black') 

                (slope,b)=get_slope(x1=9.,y1=-4.2,x2=11.,y2=-4.2)
                x_arr=np.arange(xlim[0],xlim[1]+0.05,0.05)     
                y_arr=x_arr*slope+b          
                inset.plot(x_arr-2.*np.log10(Hubble_h),y_arr,color='black', linestyle='--', linewidth=2)  
            
            #median    
            sel=(inset_x_binned>-5.) & (inset_x_binned<-2.) &(inset_median>inset_xlim[0])      
            inset.plot(inset_median[sel],inset_x_binned[sel], color=plot_color[ii], linewidth=2,linestyle='-')
      
        #LABELS   
        label="median z=%0.1f" % ThisRedshiftList[ii]
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.7, y_percentage=0.20-ii*0.034, 
                    color='black', xlog=0, ylog=0, label=label, 
                    fontsize=10, fontweight='normal', backgroundcolor='none') 
        
        plot_label (subplot, 'line', xlim, ylim,
                    x_percentage=0.66, y_percentage=0.207-ii*0.034, color=plot_color[ii], x2_percentage=0.69, 
                    xlog=0, ylog=0, linestyle='-', linewidth=2)
                      
        if(other_models==0):
            label="Quenched threshold z=%0.1f" % ThisRedshiftList[ii]
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.65, y_percentage=0.095-ii*0.034, 
                        color='black', xlog=0, ylog=0, label=label, 
                        fontsize=10, fontweight='normal', backgroundcolor='white') 
        
            plot_label (subplot, 'line', xlim, ylim,
                        x_percentage=0.59, y_percentage=0.1065-ii*0.034, color=plot_color[ii], x2_percentage=0.64, 
                        xlog=0, ylog=0, linestyle='--', linewidth=2)
        
        
        if ii==0:
            if(other_models==0):
                label='average quenching mass (0<z<2)'
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.52, y_percentage=0.34, 
                            color='black', xlog=0, ylog=0, label=label, fontsize=12, fontweight='normal',
                            rotation=3, backgroundcolor='none') 
                '''plot_label (subplot, 'label', xlim, ylim, x_percentage=0.48, y_percentage=0.43, 
                            color='black', xlog=0, ylog=0, label=label, fontsize=12, fontweight='normal',
                            rotation=4.5, backgroundcolor='none') '''
            
            if(plot_inset==1):                         
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.045, y_percentage=0.46, 
                            color='red', xlog=0, ylog=0, label='Passive (z=0)', 
                            fontsize=10, fontweight='normal') 
                plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.03, y_percentage=0.47, 
                            color='red', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.) 

                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.045, y_percentage=0.42, 
                            color='blue', xlog=0, ylog=0, label='Star Forming (z=0)', 
                            fontsize=10, fontweight='normal') 
                plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.03, y_percentage=0.43, 
                            color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.) 
            else:
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.065, y_percentage=0.84, 
                            color='red', xlog=0, ylog=0, label='Passive (z=0)', 
                            fontsize=10, fontweight='normal') 
                plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.85, 
                            color='red', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.) 

                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.065, y_percentage=0.8, 
                            color='blue', xlog=0, ylog=0, label='Star Forming (z=0)', 
                            fontsize=10, fontweight='normal') 
                plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.81, 
                            color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.)
            
            #inset
            if(plot_inset==1):
                plot_label (inset, 'label', inset_xlim, inset_ylim, x_percentage=0.18, y_percentage=0.875, 
                            color='black', xlog=0, ylog=0, label='Cooling-Heating balance', 
                            fontsize=10, fontweight='normal') 

                plot_label (inset, 'label', inset_xlim, inset_ylim, x_percentage=0.075, y_percentage=0.75, 
                            color='red', xlog=0, ylog=0, label='Passive (z=0)', 
                            fontsize=10, fontweight='normal') 
                plot_label (inset, 'symbol', inset_xlim, inset_ylim, x_percentage=0.05, y_percentage=0.775, 
                            color='red', xlog=0, ylog=0, sym='o', sym_size=5., err_size=0.) 

                plot_label (inset, 'label', inset_xlim, inset_ylim, x_percentage=0.075, y_percentage=0.65, 
                            color='blue', xlog=0, ylog=0, label='Star Forming (z=0)', 
                            fontsize=10, fontweight='normal') 
                plot_label (inset, 'symbol', inset_xlim, inset_ylim, x_percentage=0.05, y_percentage=0.675, 
                            color='blue', xlog=0, ylog=0, sym='o', sym_size=5., err_size=0.) 
                     
    plt.tight_layout()
    plt.savefig('./fig/plots_bhbm_by_sfr.pdf')
    pdf.savefig()
    plt.close()

    
    plt.rcParams.update({'font.size': 14, 'xtick.labelsize': 14, 'ytick.labelsize': 14})
#end BHBM_by_sfr





def AGN_quenching(G_MR, ThisRedshiftList, pdf):
    
    local_slope_color_cut=[0.075,0.3, 0.32]
    local_offset_color_cut=[1.85,1.18,0.99]
    #local_slope_color_cut=[0.075,0.48, 0.38]
    #local_offset_color_cut=[2.0,1.0,1.0]
    
    SSFR_cut=[-11.,-10.5,-10.]
    
    xlim=[9.0,11.5]
    ylim=[4., 10.5]   
        
    plot_color=['red','orange','green','blue']        
   
    fig = plt.figure(figsize=(7,6)) 
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
      
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
    
    xlab='$\mathrm{log_{10}}(M_{\star}[M_{\odot}])$'       
    ylab='$\mathrm{log_{10}}(M_{\mathrm{BH}}[M_{\odot}])$' 
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   

    
    for ii in range(0,len(ThisRedshiftList)):        
                    
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR_unsel=G_MR[sel]        
        G0_MR=G0_MR_unsel[(np.log10(G0_MR_unsel['StellarMass']*1.e10/Hubble_h) > xlim[0]-1.) & 
                          (np.log10(G0_MR_unsel['BlackHoleMass']*1.e10/Hubble_h) > ylim[0]-1.) &
                          (G0_MR_unsel['Type'] == 0)]
                
        x_axis=(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h))        
        y_axis=(np.log10(G0_MR['BlackHoleMass']*1.e10/Hubble_h)) 
             
        bin=0.1
        (x_binned_1, median_1, mean, pc16, pc84, rms)=median_and_percentiles (bin, ylim[0], ylim[1], y_axis, x_axis)        
        (x_binned_2, median_2, mean, pc16, pc84, rms)=median_and_percentiles (bin, xlim[0], xlim[1], x_axis, y_axis)  
        x=np.append(x_binned_2[x_binned_2<10.4],median_1[median_1>10.4])
        y=np.append(median_2[x_binned_2<10.4],x_binned_1[median_1>10.4])   
      
        #smooth the bining  
        for jj in range (0,3):
            for kk in range(1, len(y)-1):
                y[kk]=(y[kk-1]+y[kk]+y[kk+1])/3.

        for jj in range (0,3):
            for kk in range(1, len(x)-1):
                x[kk]=(x[kk-1]+x[kk]+x[kk+1])/3.

        sel=x<11.2
        #subplot.plot(x[sel], y[sel], color=plot_color[ii], linewidth=2,linestyle='-')    
            
        
        
        bin=[0.1,0.1]
        Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]
        H, xedges, yedges = np.histogram2d(x_axis, y_axis, bins=Nbins, range=[[xlim[0]-bin[0]/2., xlim[1]+bin[0]/2.], 
                                                                              [ylim[0]-bin[1]/2., ylim[1]+bin[1]/2.]]) 
        extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
        plt.subplots_adjust(bottom=0.15, left=0.15)     
        mylevels = np.log10(np.logspace(-2.0, 0.0, num=10)) 
        HH=H
        H = zoom(H, 20)        
        H=np.log10(H/np.amax(H))       
        #cont=plt.contourf(H.transpose()[::], origin='lower', cmap='Greys_r', levels=mylevels, extent=extent)  
                           
        xx=np.arange(xlim[0],xlim[1],bin[0])   
        yy=np.arange(ylim[0],ylim[1],bin[1])          
        mean2d_arr=np.zeros(len(HH[:,0]))
           
        for jj in range(0,len(HH[:,0])): 
            sel=(HH[jj,:]==np.amax(HH[jj,:]))               
            yy_aux=yy[sel]
            mean2d_arr[jj]=np.median(yy_aux)
                
        for jj in range (0,5):
            for kk in range(1, len(mean2d_arr)-1):
                mean2d_arr[kk]=(mean2d_arr[kk-1]+mean2d_arr[kk]+mean2d_arr[kk+1])/3.
         
        #sel=((xx<11.) & (xx>9.5))
        sel=((xx<11.) & (xx>9.0))
        subplot.plot(xx[sel], mean2d_arr[sel], color=plot_color[ii], linewidth=2,linestyle='-')
        
        
        
        
        
        #inset 
        ellipse = Ellipse(xy=(9.2, 5.2), width=3.0, height=1.5, angle=30,
                          fc='lightblue', edgecolor='lightblue', lw=2, alpha=0.3)
        subplot.add_patch(ellipse)   
        ellipse = Ellipse(xy=(10.85, 8.), width=4.2, height=1.2, angle=70,
                          fc='pink', edgecolor='pink', lw=2, alpha=0.3)
        subplot.add_patch(ellipse)   
        #subplot.plot(median[sel],x_binned[sel], color=plot_color[ii], linewidth=2,alpha=0.3)
        
        (slope,b)=get_slope(x1=9.,y1=7.5,x2=11.5,y2=4.5)
        x_arr=np.arange(xlim[0],xlim[1]+0.05,0.05)     
        y_arr=x_arr*slope+b          
        subplot.plot(x_arr,y_arr,color='black', linestyle=':', linewidth=2) 
           
        if(ii>0):
            x1=9.05
            x2=10.1    
            y1=8.05
            y2=8.77
                   
            x_arr=np.arange(x1,x2,0.001)     
            subplot.fill_between(x_arr,x_arr*0.+y1,x_arr*0.+y2, facecolor='lightgrey', 
                                 interpolate=True, edgecolor='black')     
            
        
        #QUENCHING THRESHOLD LINES
        if(ii==0):                    
            (slope,b)=get_slope(x1=9.,y1=5.8,x2=11.5,y2=6.5)
        else:
            if(ii==1):
                (slope,b)=get_slope(x1=9.,y1=6.25,x2=11.5,y2=6.95)
            else:
                (slope,b)=get_slope(x1=9.,y1=6.6,x2=11.5,y2=7.3)
            
            
        x_arr=np.arange(xlim[0],xlim[1]+0.05,0.05)     
        y_arr=x_arr*slope+b      
        subplot.plot(x_arr,y_arr,color=plot_color[ii], linestyle='--', linewidth=2)      
            
        if(ii==0):
            #SECOND AXIS WITH MVIR           
            StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
            Mvir=np.log10(G0_MR['Mvir']*1.e10/Hubble_h)
            (x_binned,median_MR,mean,pc16,pc84,rms)=median_and_percentiles(1.0,xlim[0],xlim[1],StellarMass,Mvir) 
            y_axis=median_MR
       
            for jj in range(0,len(y_axis)):
                y_axis[jj] = int((y_axis[jj] * 10) + 0.5) / 10.0            
             
            ax2 = subplot.twiny()        
            ax2.set_xbound(subplot.get_xbound())        
            ax2.set_xticks(x_binned)
            ax2.set_xticklabels(y_axis)

            xlab='$\mathrm{log_{10}}(<M_{\mathrm{vir}}>[M_{\odot}])$'  
            ax2.set_xlabel(xlab, fontsize=14)
            ax2.xaxis.set_label_position('top') 
            
            
            label='SN feedback'
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.12, y_percentage=0.69, 
                        color='black', xlog=0, ylog=0, label=label, fontsize=15, fontweight='normal')
            label='efficiency threshold'
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.64, 
                        color='black', xlog=0, ylog=0, label=label, fontsize=15, fontweight='normal')
                     
                
            label='SN feedback regulated'
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.02, y_percentage=0.2, 
                        color='blue', xlog=0, ylog=0, label=label, fontsize=13, fontweight='normal',rotation=10)
            label='cold-mode accretion'
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.15, 
                        color='blue', xlog=0, ylog=0, label=label, fontsize=13, fontweight='normal',rotation=10)            
                     
            label='in-situ star formation'
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.44, y_percentage=0.12, 
                        color='blue', xlog=0, ylog=0, label=label, fontsize=13, fontweight='normal')            
            label='disk formation'
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.48, y_percentage=0.07, 
                        color='blue', xlog=0, ylog=0, label=label, fontsize=13, fontweight='normal')
            label='weak BH growth'
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.47, y_percentage=0.02, 
                        color='blue', xlog=0, ylog=0, label=label, fontsize=13, fontweight='normal')
                   
        
            label='AGN feedback regulated'
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.48, y_percentage=0.82, 
                        color='red', xlog=0, ylog=0, label=label, fontsize=13, fontweight='normal',rotation=37) 
            label='hot-mode accretion'
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.54, y_percentage=0.75, 
                        color='red', xlog=0, ylog=0, label=label, fontsize=13, fontweight='normal',rotation=37) 
                        
            label='significant merger-driven growth'
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.28, y_percentage=0.9, 
                        color='red', xlog=0, ylog=0, label=label, fontsize=13, fontweight='normal') 
            label='bulge formation'
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.4, y_percentage=0.85, 
                        color='red', xlog=0, ylog=0, label=label, fontsize=13, fontweight='normal') 
            label='strong BH growth'
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.39, y_percentage=0.8, 
                        color='red', xlog=0, ylog=0, label=label, fontsize=13, fontweight='normal') 
     
            
    
            label='quenching threshold z=0'
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.62, y_percentage=0.4, 
                        color=plot_color[0], xlog=0, ylog=0, label=label, fontsize=12, fontweight='normal',
                        rotation=5, backgroundcolor='none')         
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.93, y_percentage=0.465, 
                        color=plot_color[1], xlog=0, ylog=0, label='z=1', fontsize=12, fontweight='normal',
                        rotation=5, backgroundcolor='none')         
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.93, y_percentage=0.52, 
                        color=plot_color[2], xlog=0, ylog=0, label='z=2', fontsize=12, fontweight='normal',
                        rotation=5, backgroundcolor='none')
            
            
    plt.tight_layout()
    plt.savefig('./fig/plots_AGN_quenching.pdf')
    pdf.savefig()
    plt.close()

    
    plt.rcParams.update({'font.size': 18, 'xtick.labelsize': 18, 'ytick.labelsize': 18})
#end AGN_quenching


def growth_channels(G_MR, G_MRII, ThisRedshiftList, pdf):
    
    
    xlim=[9.5,12.]
    ylim=[0.0,1.2]
    bin=0.25
    
    model_to_print='Hen15_no_feedback'
    
    for ii in range(0,len(ThisRedshiftList)):        
        
        #MODEL
        #MR
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)         
        G0_MR=G_MR[sel]                   
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
        #StellarMass=np.log10(G0_MR['Mvir']*1.e10/Hubble_h)
        InSituRatio=G0_MR['MassFromInSitu']/G0_MR['StellarMass']
        MergersRatio=G0_MR['MassFromMergers']/G0_MR['StellarMass']
                        
        (x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles (bin, xlim[0], xlim[1], StellarMass, InSituRatio) 
        sel=median!=0.
        #subplot.plot(x_binned[sel], median[sel],color='blue', linewidth=2)     
        
        x_axis=x_binned
        y_axis=median
        fa = open(Datadir+"growth_channels_"+model_to_print+"_insitu_median_z0.0.txt", "w")
        fa.write("%d\n" % len(x_axis))
        for kk in range (0,len(x_axis)):               
            fa.write("%0.2f " % x_axis[kk] + "%0.2f\n" % y_axis[kk])     
        fa.close() 
        
        
        (x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles (bin, xlim[0], xlim[1], StellarMass, MergersRatio) 
        sel=median!=0.
        #subplot.plot(x_binned[sel], median[sel],color='red', linewidth=2)     
      
        x_axis=x_binned
        y_axis=median
        fa = open(Datadir+"growth_channels_"+model_to_print+"_mergers_median_z0.0.txt", "w")
        fa.write("%d\n" % len(x_axis))
        for kk in range (0,len(x_axis)):               
            fa.write("%0.2f " % x_axis[kk] + "%0.2f\n" % y_axis[kk])     
        fa.close() 
    
          
    
    #PLOTS   
  
    fig = plt.figure(figsize=(two_two_size_small[0],two_two_size_small[1]))
    grid = gridspec.GridSpec(2, 2)
    grid.update(wspace=0.0, hspace=0.0)
    
    model_names=['Hen15_no_feedback','Hen15_only_AGN','Hen15_only_SN','Hen15']
    model_label=['no feedback','only AGN','only SN','Hen15']
        
    for ii in range(0,4):
        
        subplot=plt.subplot(grid[ii])

        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        if ii==2 or ii== 3: 
            xlab='$\mathrm{log_{10}}(\mathrm{M_{\star}}[M_{\odot}])$'    
            subplot.set_xlabel(xlab, fontsize=14)
      
        if ii==0 or ii== 2:
            ylab='Fraction'  
            subplot.set_ylabel(ylab, fontsize=14)
                       
    
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))
        subplot.yaxis.set_major_locator(MultipleLocator(0.5))
        subplot.yaxis.set_minor_locator(MultipleLocator(0.25))   
        
        if ii==1 or ii==3:
            plt.tick_params(axis='y', which='both', left='on', labelleft='off')
        if ii==0 or ii==1:    
            plt.tick_params(axis='x', which='both', bottom='on', labelbottom='off')
    
    
        fa = Datadir+"growth_channels_"+model_names[ii]+"_insitu_median_z0.0.txt"
        (x_axis,y_axis)=read_file(fa)
        subplot.plot(x_axis, y_axis,color='blue', linewidth=2)   
        
        fa = Datadir+"growth_channels_"+model_names[ii]+"_mergers_median_z0.0.txt"
        (x_axis,y_axis)=read_file(fa)
        subplot.plot(x_axis, y_axis,color='red', linewidth=2)
    
          
        #SECOND AXIS WITH MVIR   
        if ii==0 or ii==1:  
            StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[0])
            Mvir=np.log10(G0_MR['Mvir']*1.e10/Hubble_h)
            (x_binned,median_MR,mean,pc16,pc84,rms)=median_and_percentiles(1.0,xlim[0],xlim[1],StellarMass,Mvir) 
            y_axis=median_MR
       
            for jj in range(0,len(y_axis)):
                y_axis[jj] = int((y_axis[jj] * 10) + 0.5) / 10.0            
             
            ax2 = subplot.twiny()        
            ax2.set_xbound(subplot.get_xbound())        
            ax2.set_xticks(x_binned)
            ax2.set_xticklabels(y_axis)

            xlab='$\mathrm{log_{10}}(<M_{\mathrm{vir}}>[M_{\odot}])$'  
            ax2.set_xlabel(xlab, fontsize=14)
            ax2.xaxis.set_label_position('top') 
             
            subplot.tick_params(axis='x', which='both', bottom='off', labelbottom='off')
            
            
        #LABELS
        if ii==0:
                
            plot_label (subplot, 'label', xlim, ylim, 
                x_percentage=0.05, y_percentage=0.6, color='black', xlog=0, ylog=0, 
                label='Growth from:', fontsize=13, fontweight='normal') 
            
            plot_label (subplot, 'label', xlim, ylim, 
                    x_percentage=0.15, y_percentage=0.5, color='black', xlog=0, ylog=0, 
                    label='in-situ star formation', fontsize=13, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim,
                    x_percentage=0.05, y_percentage=0.53, color='blue', x2_percentage=0.12, 
                    xlog=0, ylog=0, linestyle='-', linewidth=2)
        
            plot_label (subplot, 'label', xlim, ylim, 
                    x_percentage=0.15, y_percentage=0.4, color='black', xlog=0, ylog=0, 
                    label='mergers', fontsize=13, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim,
                    x_percentage=0.05, y_percentage=0.43, color='red', x2_percentage=0.12, 
                    xlog=0, ylog=0, linestyle='-', linewidth=2)
        
        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.1, 
                    color='black', xlog=0, ylog=0, label=model_label[ii], fontsize=13, fontweight='normal') 
            
    plt.tight_layout()
    plt.savefig('./fig/plots_growth_channels.pdf')
    pdf.savefig()
    plt.close()
   

#end growth_channels   

def bluck_red_fractions(G_MR, ThisRedshiftList, pdf):
       
    ii=0
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)
                 
    G0_MR=G_MR[sel]   
               
    xmin=8.0
    xmax=13.0
    ymin=-12.
    ymax=-8.0
    bin=0.1

    fig = plt.figure(figsize=(9,9))
    grid = gridspec.GridSpec(1, 2)
    #subplotbluck=plt.subplot(grid[ii]) 
    subplot=plt.subplot()
               
        
        
    #FIND PASSIVE CUT
    subplot.set_ylim([ymin, ymax]), subplot.set_xlim([xmin, xmax])         
    xlab='$log_{10}(M_*[h^{-2}M_{\odot}])$'           
    ylab='$log_{10}(\phi [h^3 Mpc^{-3} log_{10}(M^{-1})])$'               
    subplot.set_xlabel(xlab, fontsize=14)
    subplot.set_ylabel(ylab, fontsize=14)     
        
        
    BHMass=np.log10(G0_MR['BlackHoleMass']*1.e10/Hubble_h) 
    StellarMass=np.log10(G0_MR['StellarMass']*1.e10/Hubble_h) 
    SFR=np.log10(G0_MR['Sfr'] ) 
    SSFR=np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h)) 
            
    sel=np.logical_and(SSFR>-100.,StellarMass>-10.)    
    StellarMass=StellarMass[sel]
    SSFR=SSFR[sel]
    #subplotbluck.plot(StellarMass, SSFR, 'o', markersize=2, color='blue')
             
    #subplot.hist2d(StellarMass, SSFR, bins=30, norm=LogNorm())  
    #subplot.hexbin(StellarMass, SSFR, gridsize=200)
    #plt.colorbar()
        
    #plt.tight_layout()
    pdf.savefig()
    plt.close()
        
    #PLOT RED FRACTIONS WITHOUT BINNING    
     
    G0_MR=G0_MR[G0_MR['Type']==0]           
    BHMass=np.log10(G0_MR['BlackHoleMass']*1.e10/Hubble_h)    
    SSFR=np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h)) 
    HaloMass=np.log10(G0_MR['Mvir']*1.e10/Hubble_h)   
    
    
    
    #BH MASS
    xmin=5.0
    xmax=9.0
    ymin=0.0
    ymax=1.0
    bin=0.5
    
    fig = plt.figure(figsize=(12,5))      
    subplot=plt.subplot(grid[0])
    subplot.set_ylim([ymin, ymax]), subplot.set_xlim([xmin, xmax])         
    xlab='$log_{10}(M_{BH}[M_{\odot}])$'           
    ylab='$f_{Quench}$'               
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
    subplot.yaxis.set_minor_locator(MultipleLocator(0.25))           
    #plt.tick_params(axis='y', which='both', left='on', labelleft='off')
    
    Nbins=int((xmax-xmin)/bin+1)  
    hist=np.array([],dtype=np.float64)
    x_array=np.array([],dtype=np.float64)
    Aindex=np.array(len(G0_MR),dtype=np.float64)
    Aindex=SSFR
    for ii in range(0,Nbins):              
        n_passive  = len(Aindex[(BHMass > (xmin+ii*bin)) & (BHMass < (xmin+(ii+1)*bin)) & (SSFR<-10.5)])
        n_total  = len(Aindex[(BHMass > (xmin+ii*bin)) & (BHMass < (xmin+(ii+1)*bin))])
        if n_blucktotal>0. :
            hist=np.append(hist,n_passive/n_total)
        else:
            hist=np.append(hist,0.)
        x_array=np.append(x_array,xmin+ii*bin+bin/2.)
        #endfor             
    subplot.plot(x_array,hist,color='red',linestyle='-', linewidth=2)
    
    
    #HALOMASS
    xmin=10.0
    xmax=16.0
    ymin=0.0
    ymax=1.0
    bin=0.5
    
    subplot=plt.subplot(grid[1])
    subplot.set_ylim([ymin, ymax]), subplot.set_xlim([xmin, xmax])         
    xlab='$log_{10}(M_{200c}[M_{\odot}])$'           
    ylab='$f_{Quench}$'               
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
    subplot.yaxis.set_minor_locator(MultipleLocator(0.25))           
    #plt.tick_params(axis='y', which='both', left='on', labelleft='off')
    
    Nbins=int((xmax-xmin)/bin+1)  
    hist=np.array([],dtype=np.float64)
    x_array=np.array([],dtype=np.float64)
    Aindex=np.array(len(G0_MR),dtype=np.float64)
    Aindex=SSFR
    for ii in range(0,Nbins):              
        n_passive  = len(Aindex[(HaloMass > (xmin+ii*bin)) & (HaloMass < (xmin+(ii+1)*bin)) & (SSFR<-10.5)])
        n_total  = len(Aindex[(HaloMass > (xmin+ii*bin)) & (HaloMass < (xmin+(ii+1)*bin))])
        if n_total>0. :
            hist=np.append(hist,n_passive/n_total)
        else:
            hist=np.append(hist,0.)
        x_array=np.append(x_array,xmin+ii*bin+bin/2.)
        #endfor             
    subplot.plot(x_array,hist,color='red',linestyle='-', linewidth=2)
            
        
    plt.tight_layout()
    plt.savefig('./fig/plots_bluck_red_fractions_no_binning.pdf')
    pdf.savefig()
    plt.close()
        
    #PLOT RED FRACTIONS
    xmin=5.0
    xmax=9.0
    ymin=0.0
    ymax=1.0
    bin=0.5
        
    fig = plt.figure(figsize=(12,5))      
    subplot=plt.subplot(grid[0])
    subplot.set_ylim([ymin, ymax]), subplot.set_xlim([xmin, xmax])         
    xlab='$log_{10}(M_{BH}[M_{\odot}])$'           
    ylab='$f_{Quench}$'               
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
    subplot.yaxis.set_minor_locator(MultipleLocator(0.25))           
    #plt.tick_params(axis='y', which='both', left='on', labelleft='off')
         
    BHMass=np.log10(G0_MR['BlackHoleMass']*1.e10/Hubble_h) 
    #BHMass=np.log10(G0_MR['BulgeMass']*1.e10/Hubble_h) 
    #StellarMass=np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)       
    SSFR=np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h)) 
    HaloMass=np.log10(G0_MR['Mvir']*1.e10/Hubble_h) 
        
    sel=G0_MR['Type']==0
    BHMass=BHMass[sel]
    SSFR=SSFR[sel]
    HaloMass=HaloMass[sel]
    Index_arr=HaloMass
        
    plot_colors=['blue','green','yellow','red']
    halo_min_mass=11.
    halo_bin=1.0
    for jj in range (0,4):
        print((halo_min_mass+jj*halo_bin),(halo_min_mass+(jj+1)*halo_bin))
        sel=np.logical_and(HaloMass > (halo_min_mass+jj*halo_bin), HaloMass < (halo_min_mass+(jj+1)*halo_bin))
        Nbins=int((xmax-xmin)/bin+1)  
        hist=np.array([],dtype=np.float64)
        x_array=np.array([],dtype=np.float64)
        ABHMass=BHMass[sel]
        ASSFR=SSFR[sel]
        Aindex=Index_arr[sel]
        for ii in range(0,Nbins):              
            n_passive  = len(Aindex[np.logical_and(np.logical_and(ABHMass > (xmin+ii*bin), ABHMass < (xmin+(ii+1)*bin)), 
                                           ASSFR<-10.5)])
            n_total  = len(Aindex[np.logical_and(ABHMass > (xmin+ii*bin), ABHMass < (xmin+(ii+1)*bin))])
            if n_total>0. :
                hist=np.append(hist,n_passive/n_total)
            else:
                hist=np.append(hist,0.)
            x_array=np.append(x_array,xmin+ii*bin+bin/2.)
        #endfor             
        subplot.plot(x_array,hist,color=plot_colors[jj],linestyle='-', linewidth=2)
            
    #endfor             
        
    #plt.tight_layout()
    #plt.savefig('./fig/plots.pdf')
    #pdf.savefig()
    #plt.close()
            
    #halo mass bins
    xmin=11.0
    xmax=15.0
    ymin=0.0
    ymax=1.0
    bin=0.5
        
    #fig = plt.figure(figsize=(9,9))      
    subplot=plt.subplot(grid[1]) 
    subplot.set_ylim([ymin, ymax]), subplot.set_xlim([xmin, xmax])         
    xlab='$log_{10}(M_{200c}[M_{\odot}])$'           
    ylab='$f_{Quench}$'               
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
    subplot.yaxis.set_minor_locator(MultipleLocator(0.25))           
    #plt.tick_params(axis='y', which='both', left='on', labelleft='off')
          
    BHMass=np.log10(G0_MR['BlackHoleMass']*1.e10/Hubble_h) 
    #BHMass=np.log10(G0_MR['BulgeMass']*1.e10/Hubble_h) 
    #StellarMass=np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)       
    SSFR=np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h)) 
    HaloMass=np.log10(G0_MR['Mvir']*1.e10/Hubble_h) 
        
    sel=G0_MR['Type']==0
    BHMass=BHMass[sel]
    SSFR=SSFR[sel]
    HaloMass=HaloMass[sel]
    Index_arr=HaloMass
    
    plot_colors=['blue','green','yellow','red']
    bh_min_mass=5.
    bh_bin=1.0
        
    for jj in range (0,4):
            
        print((bh_min_mass+jj*bh_bin),(bh_min_mass+(jj+1)*bh_bin))
        sel=np.logical_and(BHMass > (bh_min_mass+jj*bh_bin), BHMass < (bh_min_mass+(jj+1)*bh_bin))
            
        Nbins=int((xmax-xmin)/bin+1)  
        hist=np.array([],dtype=np.float64)
        x_array=np.array([],dtype=np.float64)
        AHaloMass=HaloMass[sel]
        ASSFR=SSFR[sel]
        Aindex=Index_arr[sel]
        for ii in range(0,Nbins):             
            n_passive  = len(Aindex[np.logical_and(np.logical_and(AHaloMass > (xmin+ii*bin), AHaloMass < (xmin+(ii+1)*bin)),
                                          ASSFR<-10.5)])           
            n_total  = len(Aindex[np.logical_and(AHaloMass > (xmin+ii*bin), AHaloMass < (xmin+(ii+1)*bin))])
            if n_total>0. :
                hist=np.append(hist,n_passive/n_total)
            else:
                hist=np.append(hist,0.)
            x_array=np.append(x_array,xmin+ii*bin+bin/2.)
        #endfor             
        subplot.plot(x_array,hist,color=plot_colors[jj],linestyle='-', linewidth=2)
            
    #endfor             
      
    plt.tight_layout()
    plt.savefig('./fig/plots_bluck_red_fractions.pdf')
    pdf.savefig()
    plt.close()
             
#end bluck_red_fractions



def satellite_quench_times(G_MR, ThisRedshiftList, pdf):
     
    plt.rcParams.update({'font.size': 10, 'xtick.labelsize': 10, 'ytick.labelsize': 10, 'axes.linewidth': 1, 
                     'xtick.major.size': 4, 'xtick.major.width': 1., 
                     'ytick.major.size': 4, 'ytick.major.width': 1., 
                     'xtick.minor.size': 2, 'xtick.minor.width': 1.,                   
                     'ytick.minor.size': 2, 'ytick.minor.width': 1.})         
        
    xlim=[-3.,8.]
    #ylim=[-13.0, -9.0]
    ylim=[0.0, 1.0]
    bin_hist=0.1
        
    mass_bin=0.5
    #mass_limits=[9.0,11.5]
    #mass_bin_arr=np.arange(mass_limits[0],mass_limits[1]+mass_bin,mass_bin)    
    #stellar_mass_low=9.4
    #stellar_mass_high=10.2
    Stellarmass_bin_arr=[9.0,9.5,10.0]   
    #Stellarmass_bin_arr=[10.0]   
    Halomass_bin_arr=[12.5,13.5,14.5]
    #Halomass_bin_arr=[14.5]
    
    color_plot=['blue','limegreen','red']     
    plot_line_style=[':','--','-']
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
 
    for ii in range(0,len(ThisRedshiftList)):        
         
        #DEFINE MODEL VARIABLES
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]             
       
        subplot=plt.subplot()                
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
        #ylab='$\mathrm{log_{10}}(\mathrm{SSFR}[yr^{-1}])$'       
        ylab='Quenched Fraction'         
        xlab='$\mathrm{Time \; Since \; Infall \;} [\mathrm{Gyr}]$'
        subplot.set_xlabel(xlab, fontsize=12), subplot.set_ylabel(ylab, fontsize=12)   
        
        #Read              
        fa = open(DirName_MR+"SFH_Bins","rb")                
        nbins =  np.fromfile(fa,np.int32,1)      
        template = np.dtype([('SnapNum',np.int32,1),
                             ('Bin',np.int32,1),
                             ('Lookbacktime',np.float64,1),                           
                             ('dt',np.float64,1),
                             ('nbins',np.int32,1)
                            ])
        SFH = np.fromfile(fa,template,int(nbins))         
        fa.close()            
        #we only need the SFH strucutre from the current snap
        SFH=SFH[SFH['SnapNum']==G0_MR['SnapNum'][0]]     
        SFH['Lookbacktime']=SFH['Lookbacktime']/1.e9
                
        G0_MR_SFH=(G0_MR['sfh_DiskMass'][:,0:len(SFH)]+G0_MR['sfh_BulgeMass'][:,0:len(SFH)])
        log_StellarMass=(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)) 
        StellarMass=(G0_MR['StellarMass']*1.e10/Hubble_h)
        log_SSFR=(np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h))) 
        color_UV=G0_MR['MagDust'][:,0]-G0_MR['MagDust'][:,2]  
        color_VJ=G0_MR['MagDust'][:,2]-G0_MR['MagDust'][:,7]       
        log_CentralMvir=(np.log10(G0_MR['CentralMvir']*1.e10/Hubble_h)) 
       
        file = "/net/bootes/export/data1/Workspace/LGal_Development_Branch/input/MRPlancksnaplist.txt"   
        zlist = Table.read(file, format='ascii')    
        #print(zlist['t/yr'][58]/1.e9)
        Age_Universe=1.38e+10/1.e9
        InfallSnapValue=40
        x=Age_Universe-SFH['Lookbacktime']-zlist['t/yr'][InfallSnapValue]/1.e9
        
        #sel will give and index to the infall time in the SFH array
        sel=(x<0.31) & (x>-0.1)  
        model_SFH=np.squeeze(G0_MR_SFH[:,sel],1)
        model_SSFR=(model_SFH*1.e10/Hubble_h)/SFH['dt'][sel]/StellarMass
        #print(model_SSFR.shape,StellarMass.shape) 
        
        for ii_stellarmass in range(0,len(Stellarmass_bin_arr)):
            for ii_halomass in range(0,len(Halomass_bin_arr)): 
                                            
                #MODEL
                sel=((log_StellarMass>Stellarmass_bin_arr[ii_stellarmass]-mass_bin/2.) & 
                     (log_StellarMass<Stellarmass_bin_arr[ii_stellarmass]+mass_bin/2.) & 
                     (log_CentralMvir>Halomass_bin_arr[ii_halomass]-mass_bin/2.) & 
                     (log_CentralMvir<Halomass_bin_arr[ii_halomass]+mass_bin/2.) & 
                     #(log_SSFR<-11.) & (G0_MR['InfallSnap']>44) & (G0_MR['InfallSnap']<50) & 
                     (G0_MR['InfallSnap']==InfallSnapValue) & 
                     (G0_MR['Type']>0))# &
                     #(np.log10(model_SSFR) > -11.0) )                    
                       
                aux_G0_MR_SFH_this_bin=G0_MR_SFH[sel,:]
                InfallSnap=G0_MR['InfallSnap'][sel]
                aux_Stellarmass=StellarMass[sel]
                               
                #print(aux_G0_MR_SFH_this_bin.shape)
                #NGals=len(G0_MR_SFH_this_bin)
                #y_arr=np.zeros(len(SFH)*len(G0_MR_SFH_this_bin),dtype=np.float32)         
                #x_arr=np.zeros(len(SFH)*len(G0_MR_SFH_this_bin),dtype=np.float32)
            
                #for jj in range(0,len(SFH)):             
                #    #y_arr[jj]=np.sum(G0_MR_SFH_this_bin[:,jj]*1.e10/Hubble_h)/len(G0_MR_SFH_this_bin[:,jj])/(SFH['dt'][jj])    
                #     y_arr[NGals*jj:NGals*(jj+1)]=(G0_MR_SFH_this_bin[:,jj]*1.e10/Hubble_h)/(SFH['dt'][jj])/StellarMass[sel]
                #     x_arr[NGals*jj:NGals*(jj+1)]=Age_Universe-SFH['Lookbacktime'][jj]-zlist['t/yr'][InfallSnap]/1.e9
            
                #bin=1.0
                #(x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles (bin, xlim[0], xlim[1], x_arr,y_arr)         
                #subplot.plot(x_binned, np.log10(median),color=color_plot[ll], linewidth=2, linestyle='-')
            
                NGals=len(aux_G0_MR_SFH_this_bin)
                bin=1.0
                x_new=np.zeros(len(np.arange(xlim[0],xlim[1],bin))*NGals,dtype=np.float32)   
                y_new=np.zeros(len(np.arange(xlim[0],xlim[1],bin))*NGals,dtype=np.float32) 
                y_new_2=np.zeros(len(np.arange(xlim[0],xlim[1],bin))*NGals,dtype=np.float32) 
                quenching_time=np.zeros(NGals,dtype=np.float32)   
                #print(NGals, Stellarmass_bin_arr[ii_stellarmass], Halomass_bin_arr[ii_halomass])
                
                for ii in range (0, len(aux_G0_MR_SFH_this_bin)):                                
                   
                    y=(aux_G0_MR_SFH_this_bin[ii,:]*1.e10/Hubble_h)/(SFH['dt'])/aux_Stellarmass[ii]                       
                    Nbins=len(np.arange(xlim[0],xlim[1],bin))
                    
                    x_new[Nbins*ii:Nbins*(ii+1)]=np.arange(xlim[0],xlim[1],bin)
                    y_new[Nbins*ii:Nbins*(ii+1)]=np.interp(x_new[Nbins*ii:Nbins*(ii+1)], x, y)         
                   
                    #tck = interpolate.splrep(x, y, s=0)
                    #y_new = interpolate.splev(x_new, tck, der=0)
                    #subplot.plot(x_new[Nbins*ii:Nbins*(ii+1)],np.log10(y_new[Nbins*ii:Nbins*(ii+1)]),
                    #             color=color_plot[ll], linewidth=2, linestyle=':') 
                   
                   
                    
                    #print(x)
                    sel=( (np.log10((aux_G0_MR_SFH_this_bin[ii,:]*1.e10/Hubble_h)/(SFH['dt'])/aux_Stellarmass[ii]) < -11.0) &
                        ((aux_G0_MR_SFH_this_bin[ii,:]*1.e10/Hubble_h)/(SFH['dt'])/aux_Stellarmass[ii] > 0.0))
                    aux_x=x[sel]   
                    
                    '''if(len(aux_x[aux_x>-3.]>0.)):                           
                        quenching_time[ii]=np.amin(aux_x[aux_x>-3.])                       
                    else:                        
                        quenching_time[ii]=10.'''
                    if(len(aux_x>0.)):                           
                        quenching_time[ii]=np.amin(aux_x)                       
                    else:                        
                        quenching_time[ii]=10.
                    
                #sel=(y_new>0.)    
                #(x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles (bin,xlim[0],xlim[1],
                #                                                                  x_new[sel],y_new[sel])   
                #subplot.plot(x_binned, np.log10(median),color=color_plot[ii_halomass], 
                #             linewidth=1, linestyle=plot_line_style[ii_stellarmass])
                
                #print(quenching_time)
              
            
            
                bin=1.0
                bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)      
                           
                hist=np.histogram(quenching_time, bins=bin_arr, range=(xlim[0],xlim[1])) 
                x_axis=hist[1][0:len(hist[1][:])-1]+bin/2.
                y_axis=hist[0]/np.sum(hist[0])            
         
                for ii in range(1,len(y_axis)):
                    y_axis[ii]+=y_axis[ii-1]
        
                subplot.plot(x_axis,y_axis, color=color_plot[ii_halomass], 
                             linewidth=1, linestyle=plot_line_style[ii_stellarmass]) 
         
                
                #x_offset=0.2*ii_stellarmass
                #x_offset=0.0
                #y_offset=0.05*ii_halomass+0.2*ii_stellarmass
                #plot_label (subplot,'line',xlim,ylim,x_percentage=0.02+x_offset,y_percentage=0.06+y_offset,
                #            color=color_plot[ii_halomass], x2_percentage=0.09+x_offset, xlog=0, ylog=0 ,
                #            linestyle=plot_line_style[ii_stellarmass],linewidth=2)
                #bin_low=Halomass_bin_arr[ii_halomass]-mass_bin/2.
                #bin_high=Halomass_bin_arr[ii_halomass]+mass_bin/2.
                #label="$M_{\mathrm{Halo}}=[$" + "%0.2f" % bin_low + "," + "%0.2f" % bin_high + "]"   
                #plot_label (subplot, 'label', xlim, ylim, x_percentage=0.1+x_offset, y_percentage=0.05+y_offset, 
                #            color='black', xlog=0, ylog=0, label=label, fontsize=8, fontweight='normal')        
            
            #bin_low=Stellarmass_bin_arr[ii_halomass]-mass_bin/2.
            #bin_high=Stellarmass_bin_arr[ii_halomass]+mass_bin/2.
            #label="$M_{\mathrm{sat}}=[$" + "%0.2f" % bin_low + "," + "%0.2f" % bin_high + "]" 
            #plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.1+y_offset, 
            #        color='black', xlog=0, ylog=0, label=label, fontsize=8, fontweight='normal') 
            
        label="$\mathrm{log_{10}}(M_{\mathrm{Halo}}[M_{\odot}])=[12.5,13.5,14.5]$"      
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.8, 
                    color='black', xlog=0, ylog=0, label=label, fontsize=10, fontweight='normal') 
        label="$\mathrm{log_{10}}(M_{\mathrm{sat}}[M_{\odot}])=[9.0,9.5,10.0]$"      
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.87, 
                    color='black', xlog=0, ylog=0, label=label, fontsize=10, fontweight='normal') 
             
        #y=np.arange(-10.5,-9.4,0.1)
        #subplot.plot(y/y-1.,y,linestyle=':', color='black')
        #plot_label (subplot, 'label', xlim, ylim, x_percentage=0.16, y_percentage=0.65, 
        #            color='black', xlog=0, ylog=0, label='Infall', fontsize=10, fontweight='normal')             
        x=np.arange(xlim[0],xlim[1]+1.0,0.1)
        subplot.plot(x,x/x-1.+0.5,linestyle=':', color='black')
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.78, y_percentage=0.57, 
                    color='black', xlog=0, ylog=0, label='More than', fontsize=9, fontweight='normal')
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.74, y_percentage=0.52, 
                    color='black', xlog=0, ylog=0, label='50% Quenched', fontsize=9, fontweight='normal')
       
        
        
            
    plt.tight_layout()
    plt.savefig('./fig/plots_satellite_quench_times.pdf')
    pdf.savefig()
    plt.close()

    
    
    
    
    
    
    
    
    
    
    
    
    
    plt.rcParams.update({'font.size': 14, 'xtick.labelsize': 14, 'ytick.labelsize': 14, 'axes.linewidth': 2, 
                     'xtick.major.size': 6, 'xtick.major.width': 1.5, 
                     'ytick.major.size': 6, 'ytick.major.width': 1.5, 
                     'xtick.minor.size': 3, 'xtick.minor.width': 1.,                   
                     'ytick.minor.size': 3, 'ytick.minor.width': 1.})  
        
    plot_xy=0
    
    if(plot_xy==1):
        plt.rcParams.update({'font.size': 18, 'xtick.labelsize': 18, 'ytick.labelsize': 18})   
        xlim=[75.0,90.]
        ylim=[85.0,100.0]
        fig = plt.figure(figsize=(10,10))      
        subplot=plt.subplot()  
    
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        xlab='$\mathrm{log_{10}}(M_*[M_{\odot}])$'      
        ylab='Satellite Fraction'
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)
        subplot.xaxis.set_major_locator(MultipleLocator(1.))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
        subplot.yaxis.set_major_locator(MultipleLocator(1.)) 
        subplot.yaxis.set_minor_locator(MultipleLocator(0.5))    
        
        sel=((G0_MR['Type']==0) & (G0_MR['Pos'][:,0]>75) & (G0_MR['Pos'][:,0]<90) &
            (G0_MR['Pos'][:,1]>85) & (G0_MR['Pos'][:,1]<100) & (G0_MR['Pos'][:,2]>90) & (G0_MR['Pos'][:,2]<92))
        G=G0_MR[sel]
    
        patches = []    
        for ii in range (0, len(G)):        
            subplot.add_artist(Circle((G['Pos'][ii,0],G['Pos'][ii,1]),radius=G['Rvir'][ii],color='black',fill=False)) 
       

        pdf.savefig()
        plt.close()

#end satellite_quench_times




def sat_fraction(G_MR, ThisRedshiftList, pdf):
    
    xlim=[8.5,11.5]
    ylim=[0., 1.] 
           
    plot_color=['red','purple']        
   
    fig = plt.figure(figsize=(15,4))
    grid = gridspec.GridSpec(1, 5)
    grid.update(wspace=0.0, hspace=0.0)
    
    for ii in range (0,len(ThisRedshiftList)):
           
        char_redshift="%0.2f" % ThisRedshiftList[ii]
            
        subplot=plt.subplot(grid[ii])
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        
        xlab='$\mathrm{log_{10}}(M_*[M_{\odot}])$'      
        if ii==0:
            ylab='Satellite Fraction'
        else:
            ylab=''      
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)
        
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
        
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
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
    

#end sat_fraction




def BHmass_in_radio(G_MR, ThisRedshiftList, pdf):
    
    xlim=[8.0,12.]
    ylim=[-5.0,0.0]
    bin=0.25
    
    plot_color=['red','purple']        
  
    fig = plt.figure(figsize=(5,4))
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
            
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.25))
           
    xlab='$\mathrm{log_{10}}(M_*[h^{-2}M_{\odot}])$'   
    ylab='$\mathrm{log_{10}}(\mathrm{M_{BHHot}/M_{BH})}$'     
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   
            
    
    
    for ii in range (0,len(ThisRedshiftList)):
           
        char_redshift="%0.2f" % ThisRedshiftList[ii]
              
        
        bin=0.25        
        Mass_arr=np.arange(xlim[0],xlim[1],bin)
        #MassInRadio=np.zeros(len(Mass_arr),dtype=np.float32)
              
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]   
        G0_MR=G0_MR[G0_MR['StellarMass']>0.]
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
        FractioninRadio=G0_MR['MassRadio']/G0_MR['StellarMass']
               
        (x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles (bin, xlim[0], xlim[1], StellarMass, FractioninRadio)    
        subplot.plot(x_binned, np.log10(median),color='red', linewidth=2)
        print(x_binned, median)        
                      
        #subplot.plot(Mass_arr, SatFraction, color='red', linestyle='-', linewidth=2) 
        
        #labels
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.55, 
                    color='black', xlog=0, ylog=0, label='z='+char_redshift[:-1], 
                    fontsize=14, fontweight='normal') 
    
        #if(ii==0):
        #    plot_label (subplot, 'label', xlim, ylim, 
        #            x_percentage=0.15, y_percentage=0.90, color='black', xlog=0, ylog=0, 
        #            label='galaxies', fontsize=13, fontweight='normal') 
        #    plot_label (subplot, 'line', xlim, ylim,
        #            x_percentage=0.04, y_percentage=0.92, color='red', x2_percentage=0.13, 
        #            xlog=0, ylog=0, linestyle='-', linewidth=2)
          
            
    plt.tight_layout()
    plt.savefig('./fig/plots_BH_radio.pdf')
    pdf.savefig()
    plt.close()

#end BHmass_in_radio



def sfr_massive_galaxies(G_MR, pdf):   
    plot_color=['red','purple']        
   
    fig = plt.figure(figsize=(15,15))
    grid = gridspec.GridSpec(3, 3)
    #grid.update(wspace=0.3, hspace=0.0)

    xlim=[0.0,4.0]
    ylim=[-13.0,-8.0]
    subplot1=plt.subplot(grid[0])
    subplot1.set_ylim(ylim), subplot1.set_xlim(xlim) 
    
    xlim=[0.0,4.0]
    ylim=[10.0,11.5]
    subplot2=plt.subplot(grid[1])    
    subplot2.set_ylim(ylim), subplot2.set_xlim(xlim)  
    
    xlim=[0.0,4.0]
    ylim=[5.0,9.0]    
    subplot3=plt.subplot(grid[2])    
    subplot3.set_ylim(ylim), subplot3.set_xlim(xlim)  
    
    xlim=[0.0,4.0]
    ylim=[10.5,14.5]
    subplot4=plt.subplot(grid[3])
    subplot4.set_ylim(ylim), subplot4.set_xlim(xlim)  
     
        
    xlim=[0.0,4.0]
    ylim=[8.0,11.0]
    subplot5=plt.subplot(grid[4])
    subplot5.set_ylim(ylim), subplot5.set_xlim(xlim)  
    
    xlim=[0.0,4.0]
    #ylim=[-3.0,3.0]
    ylim=[2.0,10.0]
    subplot6=plt.subplot(grid[5])
    subplot6.set_ylim(ylim), subplot6.set_xlim(xlim)  
    
    xlim=[0.0,4.0]
    ylim=[8.0,15.0]
    subplot7=plt.subplot(grid[6])
    subplot7.set_ylim(ylim), subplot7.set_xlim(xlim) 
    
    xlim=[0.0,4.0]
    #ylim=[-7.0,-5.0]
    ylim=[-7.0,-2.0]
    subplot8=plt.subplot(grid[7])
    subplot8.set_ylim(ylim), subplot8.set_xlim(xlim)
    
    xlim=[0.0,4.0]  
    ylim=[-3.0,0.1]
    subplot9=plt.subplot(grid[8])
    subplot9.set_ylim(ylim), subplot9.set_xlim(xlim)
    
    #format axis
    #majorFormatter = FormatStrFormatter('%d')
    #subplot1.xaxis.set_major_locator(MultipleLocator(1))    
    #subplot1.xaxis.set_minor_locator(MultipleLocator(0.25))
           
    xlab='$z$'   
    ylab='$log_{10}(\mathrm{SSFR})$'    
    subplot1.set_xlabel(xlab, fontsize=14), subplot1.set_ylabel(ylab, fontsize=14)
    ylab='$log_{10}(M_{*})$' 
    subplot2.set_xlabel(xlab, fontsize=14), subplot2.set_ylabel(ylab, fontsize=14)
    ylab='$log_{10}(M_{\mathrm{BH}})$' 
    subplot3.set_xlabel(xlab, fontsize=14), subplot3.set_ylabel(ylab, fontsize=14)
    ylab='$log_{10}(M_{\mathrm{200c}})$' 
    subplot4.set_xlabel(xlab, fontsize=14), subplot4.set_ylabel(ylab, fontsize=14)    
    ylab='$log_{10}(M_{\mathrm{Cold}})$' 
    subplot5.set_xlabel(xlab, fontsize=14), subplot5.set_ylabel(ylab, fontsize=14)
    ylab='$log_{10}(\mathrm{Cooling rate})$' 
    subplot6.set_xlabel(xlab, fontsize=14), subplot6.set_ylabel(ylab, fontsize=14)
    ylab='$log_{10}(M_{\mathrm{Hot}})$' 
    subplot7.set_xlabel(xlab, fontsize=14), subplot7.set_ylabel(ylab, fontsize=14)
    ylab='$log_{10}(\mathrm{rate ratio})$' 
    subplot8.set_xlabel(xlab, fontsize=14), subplot8.set_ylabel(ylab, fontsize=14)
    ylab='$Bulge Fraction$' 
    subplot9.set_xlabel(xlab, fontsize=14), subplot9.set_ylabel(ylab, fontsize=14)
            
    z0_snap=58          
    N_z=40  
    
    #STAR FORMING
    G0_MR=G_MR[G_MR['SnapNum']==z0_snap]      
    G0_MR=G0_MR[(G0_MR['StellarMass']*1.e10/Hubble_h>1.e10) & (G0_MR['StellarMass']*1.e10/Hubble_h<1.e11) & 
                ((G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h))>1.e-11) & (G0_MR['BlackHoleMass']>0.) &
                (G0_MR['Type']==0)]
    G0_MR=np.random.choice(G0_MR, size=30)
        
    z0_IDs=G0_MR['GalID']
    Ngals=len(z0_IDs)
   
    log_SSFR_array=np.zeros(N_z,dtype=np.float32)
    log_StellarMass_array=np.zeros(N_z,dtype=np.float32)
    log_BHMass_array=np.zeros(N_z,dtype=np.float32)    
    log_Mvir_array=np.zeros(N_z,dtype=np.float32)
    log_ColdGas_array=np.zeros(N_z,dtype=np.float32)
    log_CoolingRate_array=np.zeros(N_z,dtype=np.float32)
    log_HotGas_array=np.zeros(N_z,dtype=np.float32)
    log_RateRatio_array=np.zeros(N_z,dtype=np.float32)
    log_BulgeFract_array=np.zeros(N_z,dtype=np.float32)
    redshift=np.zeros(N_z,dtype=np.float32)
        
    #cmap = plt.get_cmap('rainbow')
    cmap = plt.get_cmap('winter')
    colors = [cmap(i) for i in np.linspace(0, 1, Ngals)]
    
    for ii in range (0, Ngals): 
    #for ii in range (0, 1):       
        Gal_allz=G_MR[(G_MR['GalID']>=G0_MR[ii]['GalID']) & (G_MR['GalID']<=G0_MR[ii]['MainLeafId'])]       
        kk=0
        min_snap=np.amin(Gal_allz['SnapNum'])
        if min_snap<(z0_snap-N_z)+1:
            min_snap=(z0_snap-N_z)+1
        #for jj in range (z0_snap-N_z+1,z0_snap+1):       
        for jj in range (min_snap,z0_snap+1):     
            sel=(Gal_allz['SnapNum']==jj)    
            if(len(Gal_allz[sel]['Sfr'])>0.):
                log_SSFR_array[kk]=np.log10(Gal_allz[sel]['Sfr']/(Gal_allz[sel]['StellarMass']*1.e10/Hubble_h))
                log_StellarMass_array[kk]=np.log10(Gal_allz[sel]['StellarMass']*1.e10/Hubble_h)
                log_BHMass_array[kk]=np.log10(Gal_allz[sel]['BlackHoleMass']*1.e10/Hubble_h)               
                log_Mvir_array[kk]=np.log10(Gal_allz[sel]['Mvir']*1.e10/Hubble_h)
                log_ColdGas_array[kk]=np.log10(Gal_allz[sel]['ColdGas']*1.e10/Hubble_h)
                log_CoolingRate_array[kk]=np.log10(Gal_allz[sel]['CoolingRate'])
                ##log_CoolingRate_array[kk]=np.log10(Gal_allz[sel]['CoolingGas']*1.e10/Hubble_h)
                log_HotGas_array[kk]=np.log10(Gal_allz[sel]['HotGas']*1.e10/Hubble_h)
                log_RateRatio_array[kk]=np.log10(Gal_allz[sel]['RadioAccretionRate']/Gal_allz[sel]['CoolingRate_beforeAGN'])
                log_BulgeFract_array[kk]=np.log10(Gal_allz[sel]['BulgeMass']/Gal_allz[sel]['StellarMass'])
                
                if(ii==0):
                    redshift[kk]=Gal_allz[sel]['Redshift']
                kk+=1
        #subplot.scatter(redshift, log_SSFR_array, s=5, color='black')   
        subplot1.plot(redshift, log_SSFR_array,color=colors[ii], linewidth=2)
        subplot2.plot(redshift, log_StellarMass_array,color=colors[ii], linewidth=2)
        subplot3.plot(redshift, log_BHMass_array,color=colors[ii], linewidth=2)       
        subplot4.plot(redshift, log_Mvir_array,color=colors[ii], linewidth=2)
        subplot5.plot(redshift, log_ColdGas_array,color=colors[ii], linewidth=2)
        subplot6.plot(redshift, log_CoolingRate_array,color=colors[ii], linewidth=2)
        subplot7.plot(redshift, log_HotGas_array,color=colors[ii], linewidth=2)
        subplot8.plot(redshift, log_RateRatio_array,color=colors[ii], linewidth=2)
        subplot9.plot(redshift, log_BulgeFract_array,color=colors[ii], linewidth=2)
    
    #PASSIVE
    G0_MR=G_MR[G_MR['SnapNum']==z0_snap]      
    G0_MR=G0_MR[(G0_MR['StellarMass']*1.e10/Hubble_h>1.e10) & (G0_MR['StellarMass']*1.e10/Hubble_h<1.e11) &
                ((G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h))<1.e-11) & (G0_MR['BlackHoleMass']>0.) &
                (G0_MR['Type']==0)]
    G0_MR=np.random.choice(G0_MR, size=200)  
          
    z0_IDs=G0_MR['GalID']
    Ngals=len(z0_IDs)
   
    log_SSFR_array=np.zeros(N_z,dtype=np.float32)
    redshift=np.zeros(N_z,dtype=np.float32)
           
    cmap = plt.get_cmap('autumn')
    #cmap = plt.get_cmap('rainbow')
    colors = [cmap(i) for i in np.linspace(0, 1, Ngals)]

    for ii in range (0, Ngals): 
    #for ii in range (0, 1): 
        Gal_allz=G_MR[(G_MR['GalID']>=G0_MR[ii]['GalID']) & (G_MR['GalID']<=G0_MR[ii]['MainLeafId'])]         
        kk=0
        
        min_snap=np.amin(Gal_allz['SnapNum'])
        if min_snap<(z0_snap-N_z)+1:
            min_snap=(z0_snap-N_z)+1
        #for jj in range (z0_snap-N_z+1,z0_snap+1):       
        for jj in range (min_snap,z0_snap+1):         
            sel=(Gal_allz['SnapNum']==jj)     
            if(len(Gal_allz[sel]['Sfr'])>0.):
                log_SSFR_array[kk]=np.log10(Gal_allz[sel]['Sfr']/(Gal_allz[sel]['StellarMass']*1.e10/Hubble_h))
                log_StellarMass_array[kk]=np.log10(Gal_allz[sel]['StellarMass']*1.e10/Hubble_h)
                log_BHMass_array[kk]=np.log10(Gal_allz[sel]['BlackHoleMass']*1.e10/Hubble_h)
                log_Mvir_array[kk]=np.log10(Gal_allz[sel]['Mvir']*1.e10/Hubble_h)
                log_ColdGas_array[kk]=np.log10(Gal_allz[sel]['ColdGas']*1.e10/Hubble_h)
                log_CoolingRate_array[kk]=np.log10(Gal_allz[sel]['CoolingRate'])
                ##log_CoolingRate_array[kk]=np.log10(Gal_allz[sel]['CoolingGas']*1.e10/Hubble_h)
                log_HotGas_array[kk]=np.log10(Gal_allz[sel]['HotGas']*1.e10/Hubble_h)
                log_RateRatio_array[kk]=np.log10(Gal_allz[sel]['RadioAccretionRate']/Gal_allz[sel]['CoolingRate_beforeAGN'])
                log_BulgeFract_array[kk]=np.log10(Gal_allz[sel]['BulgeMass']/Gal_allz[sel]['StellarMass'])
                if(ii==0):
                    redshift[kk]=Gal_allz[sel]['Redshift']
                kk+=1
        #subplot.scatter(redshift, log_SSFR_array, s=20, color=colors[ii])   
        subplot1.plot(redshift, log_SSFR_array,color=colors[ii], linewidth=2)
        subplot2.plot(redshift, log_StellarMass_array,color=colors[ii], linewidth=2)     
        subplot3.plot(redshift, log_BHMass_array,color=colors[ii], linewidth=2)
        subplot4.plot(redshift, log_Mvir_array,color=colors[ii], linewidth=2)
        subplot5.plot(redshift, log_ColdGas_array,color=colors[ii], linewidth=2)
        subplot6.plot(redshift, log_CoolingRate_array,color=colors[ii], linewidth=2)
        subplot7.plot(redshift, log_HotGas_array,color=colors[ii], linewidth=2)
        subplot8.plot(redshift, log_RateRatio_array,color=colors[ii], linewidth=2)
        subplot9.plot(redshift, log_BulgeFract_array,color=colors[ii], linewidth=2)
                                                  
    plt.tight_layout()
    plt.savefig('./fig/plots_sfr_massive_galaxies.pdf')
    pdf.savefig()
    plt.close()

#end sfr_massive_galaxies


def HotGas_fraction(G_MR, ThisRedshiftList, pdf):
    
    xlim=[10.0,15.]
    ylim=[0.0,0.2]    
    bin=[0.25,0.01]
    Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]
    plot_color=['red','purple']  
    
    fig = plt.figure(figsize=(8,6))
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
            
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.25))
           
    xlab='$\mathrm{log_{10}}(M_vir[h^{-1}M_{\odot}])$'   
    ylab='$\mathrm{log_{10}}(\mathrm{M_{Hot}/M_{vir})}$'     
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   
            
    
    
    for ii in range (0,len(ThisRedshiftList)):
           
        char_redshift="%0.2f" % ThisRedshiftList[ii]
           
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel] 
        
        Ngals=len(G0_MR)
        NN=10000.       
        sel= np.random.uniform(0.0,1.0,Ngals) < NN/Ngals 
        G00_MR=G0_MR[sel]
        plt.scatter(np.log10(G00_MR['Mvir']*1.e10), G00_MR['HotGas']/G00_MR['Mvir'], s=5, color='black')   
    
        G0_MR=G0_MR[(G0_MR['Mvir']>0.) & (G0_MR['Type']==0)]
        #G0_MR=G0_MR[(G0_MR['Mvir']>0.)]       
        Mvir=np.log10(G0_MR['Mvir']*1.e10)
        HotFraction=G0_MR['HotGas']/G0_MR['Mvir']
                  
        Ngals=len(G0_MR)    
        H, xedges, yedges = np.histogram2d(Mvir, HotFraction, bins=Nbins)            
        extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
        plt.subplots_adjust(bottom=0.15, left=0.15)        
        mylevels = np.linspace(1., Nbins[0], Nbins[0])*Ngals/(Nbins[0]**2/3.)        
        H = zoom(H, 20)        
        cont=plt.contourf(H.transpose()[::], origin='lower', cmap='Greys_r', levels=mylevels, extent=extent)        
        plt.colorbar(format='%d')    
        
            
        (x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles (bin[0], xlim[0], xlim[1], Mvir, HotFraction)    
        subplot.plot(x_binned, median,color='red', linewidth=2)
        print(x_binned, median)        
                      
        #subplot.plot(Mass_arr, SatFraction, color='red', linestyle='-', linewidth=2) 
        
        #labels
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.55, 
                    color='black', xlog=0, ylog=0, label='z='+char_redshift[:-1], 
                    fontsize=14, fontweight='normal') 
    
        #if(ii==0):
        #    plot_label (subplot, 'label', xlim, ylim, 
        #            x_percentage=0.15, y_percentage=0.90, color='black', xlog=0, ylog=0, 
        #            label='galaxies', fontsize=13, fontweight='normal') 
        #    plot_label (subplot, 'line', xlim, ylim,
        #            x_percentage=0.04, y_percentage=0.92, color='red', x2_percentage=0.13, 
        #            xlog=0, ylog=0, linestyle='-', linewidth=2)
          
            
    plt.tight_layout()
    plt.savefig('./fig/plots_HotGas_fraction.pdf')
    pdf.savefig()
    plt.close()

#end HotGas_fraction









def fabian_fb(G_MR, Volume_MR, ThisRedshiftList, pdf):
           
    xlim=[9.0,12.5]
    ylim=[-7., -0.8]
    ylim_2=[-0.3, 0.3]
    bin=0.25

    plot_colors=['red','orange','yellow','green','blue','purple']    
    
    fig = plt.figure(figsize=(7,8))
   
    gs1 = gridspec.GridSpec(1, 1)    
    gs1.update(left=0.15, right=0.95, top=0.95, bottom=0.4, wspace=0.0)
      
    gs2 = gridspec.GridSpec(2, 1)  
    gs2.update(left=0.15, right=0.95, top=0.4, bottom=-0.15, wspace=0.0, hspace=0.0)
           
           
    
    
    subplot_1=plt.subplot(gs1[0])         
    subplot_1.set_ylim(ylim), subplot_1.set_xlim(xlim)       
    xlab=''      
    ylab='$\Sigma \phi(>M_{\star})$'          
    subplot_1.set_xlabel(xlab, fontsize=14), subplot_1.set_ylabel(ylab, fontsize=14)
    plt.tick_params(axis='x', labelbottom='off')
   
    subplot_2=plt.subplot(gs2[0]) 
    subplot_2.set_ylim(ylim_2), subplot_2.set_xlim(xlim)       
    xlab='$\mathrm{log_{10}}(M_*[M_{\odot}])$'      
    ylab='$\Delta \Sigma \phi$'          
    subplot_2.set_xlabel(xlab, fontsize=14), subplot_2.set_ylabel(ylab, fontsize=14)
    majorFormatter = FormatStrFormatter('%d')
    subplot_2.yaxis.set_major_locator(MultipleLocator(0.2))      
    subplot_2.yaxis.set_minor_locator(MultipleLocator(0.1))          
    
    for ii in range(0,len(ThisRedshiftList)):
        
        char_redshift="%0.2f" % ThisRedshiftList[ii]
        
        #MODEL
        #MR
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)
        
        G0_MR=G_MR[sel]   
        G0_MR=G0_MR[G0_MR['StellarMass']>0.]
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
        StellarMass=StellarMass-2.*np.log10(Hubble_h)
        bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
        hist_MR=np.histogram(StellarMass, bins=bin_arr, range=(xlim[0],xlim[1]))   
       
        x_axis=hist_MR[1][0:len(hist_MR[1][:])-1]+bin/2.           
        hist_MR=hist_MR[0]/(Volume_MR*bin/Hubble_h**3)
        y_axis=np.zeros(len(hist_MR),dtype=np.float32)
        
        #Write
        fa = open(Datadir+"new_fabian_cumulative_ND_Hen15_"+char_redshift+".txt", "w")
        fa.write("%d\n" % len(hist_MR))
        for kk in range (0,len(hist_MR)):        
            y_axis[kk]=np.sum(hist_MR[kk:len(hist_MR)-1])             
            fa.write("%0.2f " % x_axis[kk] + "%0.2f\n" % np.log10(y_axis[kk]))         
        fa.close()
                
        #Read    
        fa = open(Datadir+"fabian_cumulative_ND_Hen15_"+char_redshift+".txt", "r") 
        index=0
        for line in fa:
            if(index==0):                
                fields = line.strip().split()    
                x_axis_hen15=np.zeros(int(fields[0]),dtype=np.float32)
                y_axis_hen15=np.zeros(int(fields[0]),dtype=np.float32)               
            else:
                fields = line.strip().split()               
                x_axis_hen15[index-1]=float(fields[0])
                y_axis_hen15[index-1]=float(fields[1])               
            index+=1    
        subplot_1.plot(x_axis_hen15,y_axis_hen15, color=plot_colors[ii], linewidth=2, linestyle='-') 
        
        fa = open(Datadir+"fabian_cumulative_ND_Hen15_fb_plus_001_"+char_redshift+".txt", "r") 
        index=0
        for line in fa:
            if(index==0):                
                fields = line.strip().split()    
                x_axis_hen15_fb_plus_001=np.zeros(int(fields[0]),dtype=np.float32)
                y_axis_hen15_fb_plus_001=np.zeros(int(fields[0]),dtype=np.float32)               
            else:
                fields = line.strip().split()               
                x_axis_hen15_fb_plus_001[index-1]=float(fields[0])
                y_axis_hen15_fb_plus_001[index-1]=float(fields[1])               
            index+=1    
        subplot_1.plot(x_axis_hen15_fb_plus_001,y_axis_hen15_fb_plus_001, color=plot_colors[ii], linewidth=2, linestyle='--') 
        
        fa = open(Datadir+"fabian_cumulative_ND_Hen15_fb_minus_001_"+char_redshift+".txt", "r") 
        index=0
        for line in fa:
            if(index==0):                
                fields = line.strip().split()    
                x_axis_hen15_fb_minus_001=np.zeros(int(fields[0]),dtype=np.float32)
                y_axis_hen15_fb_minus_001=np.zeros(int(fields[0]),dtype=np.float32)               
            else:
                fields = line.strip().split()               
                x_axis_hen15_fb_minus_001[index-1]=float(fields[0])
                y_axis_hen15_fb_minus_001[index-1]=float(fields[1])               
            index+=1    
        subplot_1.plot(x_axis_hen15_fb_minus_001,y_axis_hen15_fb_minus_001, color=plot_colors[ii], linewidth=2, linestyle=':') 
        
        subplot_2.plot(x_axis_hen15,y_axis_hen15-y_axis_hen15, 
                       color=plot_colors[ii], linewidth=2, linestyle='-') 
        subplot_2.plot(x_axis_hen15_fb_plus_001,y_axis_hen15_fb_plus_001-y_axis_hen15, 
                       color=plot_colors[ii], linewidth=2, linestyle='--') 
        subplot_2.plot(x_axis_hen15_fb_minus_001,y_axis_hen15_fb_minus_001-y_axis_hen15, 
                       color=plot_colors[ii], linewidth=2, linestyle=':') 
    
        #LABELS
        if ii==0:
            plot_label (subplot_1, 'label', xlim, ylim, x_percentage=0.81, y_percentage=0.95, 
                        color='black', xlog=0, ylog=0, label='fb=0.165', 
                        fontsize=13, fontweight='normal') 
            plot_label (subplot_1, 'line', xlim, ylim, x_percentage=0.75, y_percentage=0.96, 
                        color='black', x2_percentage=0.79, 
                        xlog=0, ylog=0, linestyle='--', linewidth=2)
            
            plot_label (subplot_1, 'label', xlim, ylim, x_percentage=0.81, y_percentage=0.91, 
                        color='black', xlog=0, ylog=0, label='fb=0.155', 
                        fontsize=13, fontweight='normal') 
            plot_label (subplot_1, 'line', xlim, ylim, x_percentage=0.75, y_percentage=0.92, 
                        color='black', x2_percentage=0.79, 
                        xlog=0, ylog=0, linestyle='-', linewidth=2)
             
            plot_label (subplot_1, 'label', xlim, ylim, x_percentage=0.81, y_percentage=0.87, 
                        color='black', xlog=0, ylog=0, label='fb=0.145', 
                        fontsize=13, fontweight='normal') 
            plot_label (subplot_1, 'line', xlim, ylim, x_percentage=0.75, y_percentage=0.88, 
                        color='black', x2_percentage=0.79, 
                        xlog=0, ylog=0, linestyle=':', linewidth=2)
      
        y_pos=[0.93,0.9,0.87,0.82,0.75,0.64]        
        plot_label (subplot_1, 'label', xlim, ylim, x_percentage=0.05, y_percentage=y_pos[ii], 
                    color='black', xlog=0, ylog=0, label='z='+char_redshift[:-1], 
                    fontsize=14, fontweight='normal')      
           
    #endfor


    #plt.tight_layout()
    plt.savefig('./fig/plots_fabian.pdf')
    pdf.savefig()
    plt.close()
#endif fabian_fb








def misc_plots(G_MR, FullSnapshotList, pdf):
       
    fig = plt.figure(figsize=(10,10))
        
    sel= np.logical_and(np.logical_and(np.logical_and(G_MR['SnapNum']==FullSnapshotList[0], G_MR['StellarMass']>0.),
          G_MR['BlackHoleMass']>0.), G_MR['DiskMass']>0.) 
    G0_MR=G_MR[sel]   
    Ngals=len(G0_MR['StellarMass']) 
  
    d = {'StellarMass' : pd.Series(np.log10((G0_MR['StellarMass'])*1.e10),index=np.zeros(Ngals,dtype=np.int32)),
         'BulgeMass' : pd.Series(np.log10((G0_MR['BulgeMass'])*1.e10),index=np.zeros(Ngals,dtype=np.int32)),
         'DiskMass' : pd.Series(np.log10((G0_MR['DiskMass'])*1.e10),index=np.zeros(Ngals,dtype=np.int32)),
         'BlackHoleMass' : pd.Series(np.log10((G0_MR['BlackHoleMass'])*1.e10),index=np.zeros(Ngals,dtype=np.int32)),
         'SSFR' : pd.Series(np.log10(G0_MR['Sfr']/G0_MR['StellarMass']),index=np.zeros(Ngals,dtype=np.int32)),
         'Type' : pd.Series(G0_MR['Type'],index=np.zeros(Ngals,dtype=np.int32)),
         'Activity':pd.Series(np.zeros(Ngals,dtype=np.str_),index=np.zeros(Ngals,dtype=np.int32))}       
       
    df = pd.DataFrame(d)    
    #df['SSFR'] = (df.weight/2000).astype(int) 
    
    NN=10000.       
    sel= np.random.uniform(0.0,1.0,Ngals) < NN/Ngals 
    df=df[sel]
    
    df.Activity[:]='Active'
    sel=df.SSFR<-10.5
    df.Activity[sel]='Passive'
    
    #df.Passive.map({1: 'Passive', 0: 'Active'})
       
    #histogram TODO
    #g = sns.FacetGrid(df, col="Activity", size=6, aspect=1)
    #g.map(sns.distplot, "StellarMass")    
          
    #scatter plot
    ##g = sns.FacetGrid(df, col='Activity', size=6, aspect=1)  
    ##g.map(plt.scatter, 'StellarMass', 'BlackHoleMass')     
   
    #linear regretion
    ##g = sns.FacetGrid(df, col='Activity', size=6, aspect=1)  
    ##g.map(sns.regplot, 'StellarMass', 'BlackHoleMass')      
        
    #contour plot    
    #g = sns.FacetGrid(df, col="Activity", row='Type', size=6, aspect=1)
    #g.map(sns.kdeplot, 'StellarMass', 'BlackHoleMass') 
    #plt.xlim(7.0, 12.0), plt.ylim(4.0, 12.0) 
    
    ##g = sns.FacetGrid(df, col="Activity", size=6, aspect=1)
    ##g.map(sns.kdeplot, 'StellarMass', 'BlackHoleMass') 
    ##plt.xlim(7.0, 12.0), plt.ylim(4.0, 10.0)  
    
    #multiple variables
    ##this one just produces scater plots (pairplot) with the option of a hist in the diagonal
    ##g = sns.pairplot(df[["StellarMass", "BlackHoleMass", "Activity"]], hue="Activity", diag_kind="hist")  
    ##for ax in g.axes.flat:  
    ##    plt.setp(ax.get_xticklabels(), rotation=45)
    
    #We were able to control three regions (the diagonal, the lower-left triangle, and the 
    #upper-right triangle) separately. Again, you can pipe in any plotting function that 
    #understands the data it's given.    
    ##g = sns.PairGrid(df[["StellarMass", "BlackHoleMass", "Activity"]], hue="Activity")  
    ##g.map_upper(sns.regplot)  
    ##g.map_lower(sns.residplot)  
    ##g.map_diag(plt.hist)  
    ##for ax in g.axes.flat:  
    ##    plt.setp(ax.get_xticklabels(), rotation=45)
    ##g.add_legend()  
    ##g.set(alpha=0.5)  
    
    ax = sns.kdeplot(df.StellarMass, df.BlackHoleMass,
                  cmap="Reds", shade=True, shade_lowest=False)
    ax = sns.kdeplot(df.StellarMass, df.DiskMass,
                  cmap="Blues", shade=True, shade_lowest=False)
    
    
    #joint plots
    sns.jointplot("StellarMass", "BlackHoleMass", data=df, kind='kde')  
    
    g = sns.JointGrid(x="StellarMass", y="BlackHoleMass", data=df)  
    g.plot_joint(sns.kdeplot, shade=True, n_levels=20, cut=10., bw=0.2, shade_lowest=False)  
    g.plot_marginals(sns.distplot)  
    
    
    #g = sns.JointGrid(x="StellarMass", y="BlackHoleMass", data=df)  
    #g.plot_joint(sns.regplot, order=2)  
    #g.plot_marginals(sns.distplot)  
    
    #pdf.savefig()
    #plt.close()
 
#end misc_plots
     
    
    
def test_resolution(G_MR, Volume_MR, G_MRII, Volume_MRII, ThisRedshiftList, pdf):
   
    fig = plt.figure(figsize=(two_two_size_small[0],two_two_size_small[1]))
    
    
    grid = gridspec.GridSpec(2, 2)
    grid.update(wspace=0.0, hspace=0.0)
       
        
    #Cold Gas    
    xlim=[7.0,11.5]
    ylim=[-5.5,-0.5]
    bin=0.25
    subplot=plt.subplot(grid[0])
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.25))
    subplot.yaxis.set_major_locator(MultipleLocator(1))    
    subplot.yaxis.set_minor_locator(MultipleLocator(0.25))   
    plt.tick_params(axis='x', which='both', top='on', labeltop='on', bottom='off', labelbottom='off')
    
    xlab='$\log_{10}(M_{\mathrm{Cold}}[M_{\odot}])$'   
    ylab='$\log_{10}(\phi [\mathrm{Mpc^{-3}} \mathrm{log_{10}}(M^{-1})])$'   
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)
    subplot.xaxis.set_label_position('top')  
        
    #MR  
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, 0, FullSnapshotList_MR)         
    G0_MR=G_MR[sel]   
    G0_MR=G0_MR[(G0_MR['ColdGas']>0.) & (G0_MR['Type']==0)]
    mass=(np.log10(G0_MR['ColdGas']*1.e10/Hubble_h))

    bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
    hist=np.histogram(mass, bins=bin_arr, range=(xlim[0],xlim[1]))   
    subplot.plot(hist[1][0:len(hist[1][:])-1]+bin/2.,np.log10(hist[0][:]/(Volume_MR*bin)),
                 color='red', linewidth=2, linestyle='--')

    #MRII
    (sel)=select_current_redshift(G_MRII, ThisRedshiftList, 0, FullSnapshotList_MRII)         
    G0_MRII=G_MRII[sel]   
    G0_MRII=G0_MRII[(G0_MRII['ColdGas']>0.) & (G0_MRII['Type']==0)]
    mass=(np.log10(G0_MRII['ColdGas']*1.e10/Hubble_h))

    bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
    hist=np.histogram(mass, bins=bin_arr, range=(xlim[0],xlim[1]))   
    subplot.plot(hist[1][0:len(hist[1][:])-1]+bin/2.,np.log10(hist[0][:]/(Volume_MRII*bin)),
                 color='red', linewidth=2)
        
    plot_label (subplot, 'label', xlim, ylim, 
                x_percentage=0.15, y_percentage=0.20, color='black', xlog=0, ylog=0, 
                label='ColdGas', fontsize=15, fontweight='normal')  
    
    
        
    #SFR  
    xlim=[-3.5,2.]    
    bin=0.25
    subplot=plt.subplot(grid[1])
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.25))
    subplot.yaxis.set_major_locator(MultipleLocator(1))    
    subplot.yaxis.set_minor_locator(MultipleLocator(0.25))   
    plt.tick_params(axis='x', which='both', top='on', labeltop='on', bottom='off', labelbottom='off')   
    plt.tick_params(axis='y', which='both', left='on', labelleft='off')
    
    xlab='$\log_{10}(\mathrm{SFR}[M_{\odot}yr^{-1}])$'  
    subplot.set_xlabel(xlab, fontsize=14)
    subplot.xaxis.set_label_position('top')  
    
    #MR
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, 0, FullSnapshotList_MR)         
    G0_MR=G_MR[sel]   
    G0_MR=G0_MR[G0_MR['Sfr']>0.]
    mass=(np.log10(G0_MR['Sfr']/Hubble_h))

    bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
    hist=np.histogram(mass, bins=bin_arr, range=(xlim[0],xlim[1]))   
    subplot.plot(hist[1][0:len(hist[1][:])-1]+bin/2.,np.log10(hist[0][:]/(Volume_MR*bin)),
                 color='red', linewidth=2, linestyle='--')

    #MRII
    (sel)=select_current_redshift(G_MRII, ThisRedshiftList, 0, FullSnapshotList_MRII)         
    G0_MRII=G_MRII[sel]   
    G0_MRII=G0_MRII[G0_MRII['Sfr']>0.]
    mass=(np.log10(G0_MRII['Sfr']/Hubble_h))

    bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
    hist=np.histogram(mass, bins=bin_arr, range=(xlim[0],xlim[1]))   
    subplot.plot(hist[1][0:len(hist[1][:])-1]+bin/2.,np.log10(hist[0][:]/(Volume_MRII*bin)),
                 color='red', linewidth=2)
    
    plot_label (subplot, 'label', xlim, ylim, 
                x_percentage=0.15, y_percentage=0.20, color='black', xlog=0, ylog=0, 
                label='SFR', fontsize=15, fontweight='normal') 
    
    
    #Stellar Mass   
    xlim=[7.0,12.5]   
    bin=0.25
    subplot=plt.subplot(grid[2])
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.25))
    subplot.yaxis.set_major_locator(MultipleLocator(1))    
    subplot.yaxis.set_minor_locator(MultipleLocator(0.25))   
    plt.tick_params(axis='x', which='both', top='off', labeltop='off')  
    
    xlab='$\log_{10}(M_*[M_{\odot}])$'   
    ylab='$\log_{10}(\phi [\mathrm{Mpc^{-3}} \mathrm{log_{10}}(M^{-1})])$'   
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)
   
    
    #MR
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, 0, FullSnapshotList_MR)         
    G0_MR=G_MR[sel]   
    G0_MR=G0_MR[G0_MR['StellarMass']>0.]
    mass=(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h))

    bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
    hist=np.histogram(mass, bins=bin_arr, range=(xlim[0],xlim[1]))   
    subplot.plot(hist[1][0:len(hist[1][:])-1]+bin/2.,np.log10(hist[0][:]/(Volume_MR*bin)),
                 color='red', linewidth=2, linestyle='--')

    #MRII
    (sel)=select_current_redshift(G_MRII, ThisRedshiftList, 0, FullSnapshotList_MRII)         
    G0_MRII=G_MRII[sel]   
    G0_MRII=G0_MRII[G0_MRII['StellarMass']>0.]
    mass=(np.log10(G0_MRII['StellarMass']*1.e10/Hubble_h))

    bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
    hist=np.histogram(mass, bins=bin_arr, range=(xlim[0],xlim[1]))   
    subplot.plot(hist[1][0:len(hist[1][:])-1]+bin/2.,np.log10(hist[0][:]/(Volume_MRII*bin)),
                 color='red', linewidth=2)
    
    plot_label (subplot, 'label', xlim, ylim, 
                x_percentage=0.15, y_percentage=0.20, color='black', xlog=0, ylog=0, 
                label='StellarMass', fontsize=15, fontweight='normal') 
           
        
    plot_label (subplot, 'label', xlim, ylim, 
                x_percentage=0.15, y_percentage=0.5, color='black', xlog=0, ylog=0, 
                label='MR', fontsize=13, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim,
                x_percentage=0.04, y_percentage=0.53, color='red', x2_percentage=0.13, 
                xlog=0, ylog=0, linestyle='--', linewidth=2)
    plot_label (subplot, 'label', xlim, ylim, 
                x_percentage=0.15, y_percentage=0.4, color='black', xlog=0, ylog=0, 
                label='MRII', fontsize=13, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim,
                x_percentage=0.04, y_percentage=0.43, color='red', x2_percentage=0.13, 
                xlog=0, ylog=0, linestyle='-', linewidth=2)    
        
    #BH Mass    
    xlim=[3.5,10.]   
    bin=0.25
    subplot=plt.subplot(grid[3])
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.25))
    subplot.yaxis.set_major_locator(MultipleLocator(1))    
    subplot.yaxis.set_minor_locator(MultipleLocator(0.25))   
    plt.tick_params(axis='x', which='both', top='off', labeltop='off') 
    plt.tick_params(axis='y', which='both', left='on', labelleft='off')
    
    xlab='$\log_{10}(M_{\mathrm{BH}}[M_{\odot}])$'         
    subplot.set_xlabel(xlab, fontsize=14)
    
          
    #MR  
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, 0, FullSnapshotList_MR)         
    G0_MR=G_MR[sel]   
    G0_MR=G0_MR[(G0_MR['BlackHoleMass']>0.) & (G0_MR['Type']==0)]
    mass=(np.log10(G0_MR['BlackHoleMass']*1.e10/Hubble_h))

    bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
    hist=np.histogram(mass, bins=bin_arr, range=(xlim[0],xlim[1]))   
    subplot.plot(hist[1][0:len(hist[1][:])-1]+bin/2.,np.log10(hist[0][:]/(Volume_MR*bin)),
                 color='red', linewidth=2, linestyle='--')

    #MRII
    (sel)=select_current_redshift(G_MRII, ThisRedshiftList, 0, FullSnapshotList_MRII)         
    G0_MRII=G_MRII[sel]   
    G0_MRII=G0_MRII[(G0_MRII['BlackHoleMass']>0.) & (G0_MRII['Type']==0)]
    mass=(np.log10(G0_MRII['BlackHoleMass']*1.e10/Hubble_h))

    bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
    hist=np.histogram(mass, bins=bin_arr, range=(xlim[0],xlim[1]))   
    subplot.plot(hist[1][0:len(hist[1][:])-1]+bin/2.,np.log10(hist[0][:]/(Volume_MRII*bin)),
                 color='red', linewidth=2)
        
    plot_label (subplot, 'label', xlim, ylim, 
                x_percentage=0.15, y_percentage=0.20, color='black', xlog=0, ylog=0, 
                label='BlackHoleMass', fontsize=15, fontweight='normal') 
    
    
    
    plt.tight_layout()
    plt.savefig('./fig/plots_test_resolution_rings.pdf')
    plt.savefig('./fig/HYW17_plots_test_resolution_rings.pdf')
    pdf.savefig()
    plt.close()
    
#end test_resolution_rings
    
    
def simple_tree_map(G_MR, pdf):
    
    fig = plt.figure(figsize=(15,15))
    subplot=plt.subplot()
    subplot.set_xlim([0.0,300.]), subplot.set_ylim([-5.0, 150.]) 
    subplot.set_xlabel('x', fontsize=14), subplot.set_ylabel('y', fontsize=14)
               
    G0_MR=G_MR[(G_MR['SnapNum']==255) & (G_MR['StellarMass']>1.0) & (G_MR['Type']==0)]      
    G0_MR=G_MR[(G_MR['GalID']>=G0_MR['GalID']) & (G_MR['GalID']<=G0_MR['LastProgGal'])]    
      
    print('Ngals in tree=',len(G0_MR))
        
    for jj in range (150,255):
        G_aux=G0_MR[G0_MR['SnapNum']==jj]
        G_desc=G0_MR[G0_MR['SnapNum']==jj+1]
     
        for ii in range (0,len(G_aux)):           
            if(G_aux['Type'][ii]==0):
                plot_color='blue'
            else:
                if(G_aux['Type'][ii]==1):
                    plot_color='green'
                else:
                    plot_color='red'
            mass=G_aux['StellarMass'][ii]*100.
            if(mass==0.):
                mass=1.                
            subplot.scatter(2+ii*10.0,255-jj,s=mass, color=plot_color)
            for kk in range (0,len(G_desc)):
                if(G_desc['GalID'][kk]==G_aux['DescendantGal'][ii]):
                    subplot.plot([2+ii*10.0,2+kk*10.0],[255-jj,255-jj-1],color='black')
            
#simple_tree_map


  
def full_tree_map(G_MR, pdf, object_type):       
    
    print('Doing tree for',object_type)
    
    #G0_MR=G_MR[(G_MR['SnapNum']==255) & (G_MR['StellarMass']>1.0) & (G_MR['Type']==0)]
    #G0_MR=G_MR[(G_MR['SnapNum']==255) & (G_MR['StellarMass']>.35) & (G_MR['Type']==0)]
    #G0_MR=G_MR[(G_MR['SnapNum']==255) & (G_MR['StellarMass']<0.03) & (G_MR['StellarMass']>0.01) & (G_MR['Type']==0)] 
    #print(len(G0_MR))
    G0_MR=G_MR[(G_MR['GalID']==4000000000000)]
    G0_MR=G_MR[(G_MR['GalID']>=G0_MR['GalID']) & (G_MR['GalID']<=G0_MR['LastProgGal'])]    
        
    print('Ngals in tree=',len(G0_MR['GalID']))
    
    done=np.zeros(len(G0_MR),dtype=np.float32)
    pos=np.zeros(len(G0_MR),dtype=np.float32)
    plot_snap=np.zeros(len(G0_MR),dtype=np.float32)
    marker_size=np.zeros(len(G0_MR),dtype=np.float32)
    plot_type=np.zeros(len(G0_MR),dtype=np.float32)
    
    gal=G0_MR[G0_MR['SnapNum']==255]
    prog=gal['FirstProgGal']
 
    print('ID of final gal=',gal['GalID'])
    ii=0
  
    #for ll in range (0,100):
    while(prog>0):     
        sel=G0_MR['GalID']==prog
        G_sel=G0_MR[sel]
        #if gal not done and there is a progenitor, plot galaxy and move to progenitor
        if((done[sel]==0) & (G_sel['FirstProgGal']>0)):           
            plot_type[sel]=G_sel['Type']            
            pos[sel]=ii
            plot_snap[sel]=255-G_sel['SnapNum'] 
            if(object_type=='galaxies'):
                marker_size[sel]=G_sel['StellarMass']      
            else:
                marker_size[sel]=G_sel['Mvir']      
            done[sel]=1
            prog=G_sel['FirstProgGal']
        else:
            #if galaxy not done but there is no progenitor, plot galaxy and move
            #either to the next prog (if it exists) or descendant (if nextprog=-1)
            if(done[sel]==0):
                plot_type[sel]=G_sel['Type']                
                pos[sel]=ii
                plot_snap[sel]=255-G_sel['SnapNum'] 
                if(object_type=='galaxies'):
                    marker_size[sel]=G_sel['StellarMass']      
                else:
                    marker_size[sel]=G_sel['Mvir']                  
                done[sel]=1
                
                if(G_sel['NextProgGal']>0):
                    prog=G_sel['NextProgGal']                     
                    ii+=1                     
                else:
                    prog=G_sel['DescendantGal']
            #if galaxy done move either to the next prog (if it exists) 
            #or descendant (if nextprog=-1)
            else:                                       
                if(G_sel['NextProgGal']>0):
                    prog=G_sel['NextProgGal']                   
                    ii+=1                    
                else:
                    prog=G_sel['DescendantGal']
        
    fig = plt.figure(figsize=(15,15))
    subplot=plt.subplot()
    subplot.set_xlim([-20.0,20.+np.amax(pos)*5.+50.]), subplot.set_ylim([-5.0, 255.]) 
    #subplot.set_xlim([-20.0,3700.]), subplot.set_ylim([40.0, 80.]) 
    subplot.set_xlabel('', fontsize=14), subplot.set_ylabel('Snap', fontsize=14)         
    subplot.xaxis.set_ticklabels([])
    
    if(object_type=='galaxies'):
        marker_scale=100.     
    else:
        marker_scale=5.   
                
    sel=((plot_type==0) & (marker_size>0.))    
    subplot.scatter(20+pos[sel]*5.0,plot_snap[sel],s=marker_size[sel]*marker_scale, color='blue')
    sel=((plot_type==1) & (marker_size>0.))   
    subplot.scatter(20+pos[sel]*5.0,plot_snap[sel],s=marker_size[sel]*marker_scale, color='green')
    sel=((plot_type==2) & (marker_size>0.))   
    subplot.scatter(20+pos[sel]*5.0,plot_snap[sel],s=marker_size[sel]*marker_scale, color='red')    
    sel=(marker_size==0.)  
    subplot.scatter(20+pos[sel]*5.0,plot_snap[sel],s=marker_size[sel]+marker_scale/100., color='black')
    
    print('done tree walk')    
                
    for jjj in range (0,len(G0_MR)):
        prog=G0_MR['FirstProgGal'][jjj] 
        sel=G0_MR['GalID']==prog
        G_sel=G0_MR[sel]
        if(prog>0):          
            subplot.plot([20+pos[jjj]*5.0,20+pos[sel]*5.0],
                         [255-G0_MR['SnapNum'][jjj],255-G_sel['SnapNum']],color='black')
            prog=G_sel['NextProgGal']
            sel=G0_MR['GalID']==prog
            G_sel=G0_MR[sel]
            while(prog>0):
                subplot.plot([20+pos[jjj]*5.0,20+pos[sel]*5.0],
                             [255-G0_MR['SnapNum'][jjj],255-G_sel['SnapNum']],color='black')
                prog=G_sel['NextProgGal']
                sel=G0_MR['GalID']==prog
                G_sel=G0_MR[sel]
               
    print('done second tree walk')      
             
    plt.tight_layout()
    plt.savefig('./fig/plots_full_tree_map.pdf')
    pdf.savefig()
    plt.close()

#full_tree_map   
