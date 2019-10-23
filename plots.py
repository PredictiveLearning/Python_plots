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
'''
import numpy as np
import pandas as pd
#import seaborn as sns
#sns.set_style('darkgrid')
import time
import inspect   
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.table import Table
from astropy.io import fits
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
from scipy import interpolate
from scipy.stats import binned_statistic_2d
from scipy.optimize import curve_fit
from scipy.integrate import simps
from scipy.integrate import newton_cotes
import sys
import csv
from scipy.ndimage import zoom
import os.path
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Ellipse, Circle
from matplotlib.collections import PatchCollection
from importlib import reload


import MillenniumSQL 
reload (MillenniumSQL)
from MillenniumSQL import *

import procedures
reload (procedures)
from procedures import *
import plots_input
reload (plots_input)
from plots_input import *


def stellar_mass_vs_halo_mass_vandokkum(ThisRedshiftList):
   
    ylim=[2.0,12.5]
    xlim=[2.5, 15.]
    
    fig = plt.figure(figsize=(10,10))
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim) 

    ylab='$\mathrm{log_{10}}(M_*/$'+Msun+'$)$'       
    xlab='$\mathrm{log_{10}}(M_{200c}/$'+Msun+'$)$'
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)
                  
    for ii in range (0,len(ThisRedshiftList)):
  
        #MODEL
        (sel)=select_current_redshift(G_MRII, ThisRedshiftList, ii, FullSnapshotList_MRII)        
        G0_MR=G_MRII[sel] 
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
    #pdf.savefig()
    #plt.close()
    
    
    xlim=[5.0, 11.]
    ylim=[7.0,13.]
    
    
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim) 

    xlab='$\mathrm{log_{10}}(M_*/$'+Msun+'$)$'       
    ylab='$\mathrm{log_{10}}(M_{200c}/$'+Msun+'$)$'
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)
    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.1)) 
    subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
    
    for ii in range (0,len(ThisRedshiftList)):
  
        #MODEL
        (sel)=select_current_redshift(G_MRII, ThisRedshiftList, ii, FullSnapshotList_MRII)        
        G0_MR=G_MRII[sel] 
        G0_MR=np.random.choice(G0_MR, size=100000)        
        G0_MR=G0_MR[G0_MR['StellarMass']>0.]
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])          
        HaloMass=np.log10((G0_MR['Len']*0.000768884)*1.e10/Hubble_h)       
        sel=((G0_MR['Type']==1) | (G0_MR['Type']==2))
        subplot.scatter(StellarMass[sel],HaloMass[sel],s=5, color='red')        
        sel=G0_MR['Type']==0
        subplot.scatter(StellarMass[sel],HaloMass[sel],s=5, color='blue')
  
         #MODEL
        (sel)=select_current_redshift(G_MRII, ThisRedshiftList, ii, FullSnapshotList_MRII)   
        G0_MR=G_MRII[sel] 
        G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & ((G0_MR['StellarMass']/(G0_MR['Len']*0.000768884))>0.1)]
        
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])          
        HaloMass=np.log10((G0_MR['Len']*0.000768884)*1.e10/Hubble_h)       
        #subplot.scatter(StellarMass,HaloMass,s=5, color='green') 
        
        print(G0_MR['StellarMass'][0:9]/(G0_MR['Len'][0:9]*0.000768884))
        
        x = [2*1e8]
        y = [1.5*1e8]
        subplot.scatter(np.log10(x),np.log10(y),s=50, marker='s', color='black')
        
        xx = np.arange(xlim[0],xlim[1],0.1)
        yy = xx
        subplot.plot(xx,yy,linestyle=':', color='black')
        yy = xx+1
        subplot.plot(xx,yy,linestyle=':', color='black')
        yy = xx+2
        subplot.plot(xx,yy,linestyle=':', color='black')
        yy = xx+3
        subplot.plot(xx,yy,linestyle=':', color='black')
    #endfor
    
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')   
    plt.close()
       
    return 

def stellar_mass_vs_halo_mass_fractional(ThisRedshiftList):

    #model_names=['Hen15_no_feedback','Hen15_only_AGN','Hen15_only_SN','Hen15']
    
    model_to_print='Hen15_other'    
    #TO WRITE OUT and plot FOR A GIVEN MODEL
    #PLOT for all models IS ACTUALLY MADE FROM FILES PREVIOUSLY READ
    
        
    xlim=[11., 14.5]
    ylim=[-3.,0.]    
    
    
    #PLOT current model
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1.0))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
    subplot.yaxis.set_minor_locator(MultipleLocator(0.5)) 
    
    xlab='$\log_{10}(M_{200c}/$'+Msun+'$)$' 
    ylab='$\log_{10}(M_*/M_{200c})$'
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14) 

    for ii in range (0,len(ThisRedshiftList)):
  
        #MODEL
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]             
        G0_MR=G0_MR[(G0_MR['StellarMass']>0.0) & (G0_MR['Type']==0) & 
                    (G0_MR['Mvir']*1e10/Hubble_h>xlim[0])]        
        G0_MR=np.random.choice(G0_MR, size=3000, replace=False)   
        log_StellarMass=stellar_mass_with_err(G0_MR, Hubble_h,ThisRedshiftList[ii])
        
        #sel=G0_MR['ICM']>0.
        #log_StellarMass[sel]+=np.log10(G0_MR['ICM'][sel]*1.e10/Hubble_h)              
        log_HaloMass=np.log10(G0_MR['Mvir']*1.e10/Hubble_h)         
        log_fraction=log_StellarMass-log_HaloMass       
        SSFR=G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h)
        
        '''log_StellarMass=np.zeros(len(G0_MR),dtype=np.float32)
        log_HaloMass=np.zeros(len(G0_MR),dtype=np.float32)
        log_fraction=np.zeros(len(G0_MR),dtype=np.float32)
        
        for jj in range (0,len(G0_MR)):
            dist=np.sqrt( (G0_MR['Pos'][:,0]-G0_MR['Pos'][jj,0])**2 +
                       (G0_MR['Pos'][:,1]-G0_MR['Pos'][jj,1])**2 + 
                       (G0_MR['Pos'][:,1]-G0_MR['Pos'][jj,1])**2 )
            sel=dist<G0_MR['Rvir'][jj]
            log_StellarMass[jj]=np.log10(np.sum(G0_MR['StellarMass'][sel]*1.e10/Hubble_h, axis=0) + 
                                         G0_MR['ICM'][jj]*1.e10/Hubble_h)
            #log_StellarMass[jj]=np.log10(np.sum(G0_MR['StellarMass'][sel]*1.e10/Hubble_h, axis=0))
            log_HaloMass[jj]=np.log10(G0_MR['Mvir'][jj]*1.e10/Hubble_h)           
            log_fraction[jj]=log_StellarMass[jj]-log_HaloMass[jj]'''
           
        #WRITE OUTPUT            
        #SF
        #redshift = 0.
        #t_h = 13.812*1.e9 
        #sel=SSFR>((1.+redshift)/(2.*t_h))
        sel = np.log10(SSFR)>-11.0
        subplot.scatter(log_HaloMass[sel],log_fraction[sel],s=5, color='blue') 
        x_axis=log_HaloMass[sel]
        y_axis=log_fraction[sel]
        fa = open(Datadir+"SMHM_"+model_to_print+"_SF_z0.0.txt", "w")
        fa.write("%d\n" % len(x_axis))
        for kk in range (0,len(x_axis)):               
            fa.write("%0.2f " % x_axis[kk] + "%0.2f\n" % y_axis[kk])         
        fa.close() 
            
        #PASSIVE
        #sel=SSFR<((1.+redshift)/(2.*t_h))
        sel = np.log10(SSFR)<-11.0
        subplot.scatter(log_HaloMass[sel],log_fraction[sel],s=5, color='red')         
        x_axis=log_HaloMass[sel]
        y_axis=log_fraction[sel]
        fa = open(Datadir+"SMHM_"+model_to_print+"_passive_z0.0.txt", "w")
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
        subplot.plot(x_binned, median,color='black', linewidth=2, linestyle='-') 
        
       
        '''x_axis=x_binned
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
        fa.close()'''
       
        
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)  
        G0_MR=G_MR[sel]  
        G0_MR_centrals = G0_MR[(G0_MR['StellarMass']>0.1) & (G0_MR['Mvir']*1e10/Hubble_h>1e12) & (G0_MR['Type']==0)]     
           
        G0_MR=G_MR[sel]             
        G0_MR=G0_MR[(G0_MR['StellarMass']>0.1) & (G0_MR['Mvir']*1e10/Hubble_h>1e11)] 
       
        
        log_StellarMass=np.zeros(len(G0_MR_centrals),dtype=np.float32)
        log_HaloMass=np.zeros(len(G0_MR_centrals),dtype=np.float32)
        log_fraction=np.zeros(len(G0_MR_centrals),dtype=np.float32)
        
        
        ''' print('N_centrals=',len(G0_MR_centrals))
        
        G0_MR_centrals = G0_MR_centrals[0:1]
        
        n_sats=len(G0_MR)
        n_centrals=len(G0_MR_centrals)
        dist=np.zeros(n_centrals*n_sats,dtype=np.float32)
        Rvir=np.zeros(n_centrals*n_sats,dtype=np.float32)
        
        print('replicating arrays')
        centrals_sel_x=np.repeat(G0_MR_centrals['Pos'][:,0],n_sats)
        centrals_sel_y=np.repeat(G0_MR_centrals['Pos'][:,1],n_sats)
        centrals_sel_z=np.repeat(G0_MR_centrals['Pos'][:,2],n_sats)
        sat_sel_x=np.repeat(G0_MR['Pos'][:,0],n_centrals)
        sat_sel_y=np.repeat(G0_MR['Pos'][:,1],n_centrals)
        sat_sel_z=np.repeat(G0_MR['Pos'][:,2],n_centrals)
        
            
        print('computing dist') 
        dist= np.sqrt((centrals_sel_x-sat_sel_x)*(centrals_sel_x-sat_sel_x) + 
                      (centrals_sel_y-sat_sel_y)*(centrals_sel_y-sat_sel_y) + 
                      (centrals_sel_z-sat_sel_z)*(centrals_sel_z-sat_sel_z))
        Rvir = np.repeat(G0_MR_centrals['Rvir'],n_sats)
        Sat_Mass =  np.repeat(G0_MR['StellarMass']*1.e10/Hubble_h,n_centrals)
        
        print(np.shape(dist))
        print(np.shape(Rvir))
        sel = dist < Rvir
        print(Sat_Mass[sel])
        #sel = G0_MR['FOFCentralGal'] == G0_MR_centrals['FOFCentralGal'][jj] 
        log_StellarMass[jj]=np.log10(np.sum(G0_MR['StellarMass'][sel]*1.e10/Hubble_h, axis=0) + 
                                         G0_MR['ICM'][jj]*1.e10/Hubble_h)            
        log_HaloMass[jj]=np.log10(G0_MR_centrals['Mvir'][jj]*1.e10/Hubble_h)             
        log_fraction[jj]=log_StellarMass[jj]-log_HaloMass[jj]'''
        
        #centrals+ satellites and ICM        
        '''for jj in range (0,len(G0_MR_centrals)):
            
            if(jj%int(len(G0_MR_centrals)/100)==0):
                print(int(jj/len(G0_MR_centrals)*100.), '% done')
                
            #sel = G0_MR['CentralMvir']==G0_MR_centrals['Mvir'][jj]
            dist=np.sqrt( (G0_MR['Pos'][:,0]-G0_MR_centrals['Pos'][jj,0])**2 +
                          (G0_MR['Pos'][:,1]-G0_MR_centrals['Pos'][jj,1])**2 + 
                          (G0_MR['Pos'][:,2]-G0_MR_centrals['Pos'][jj,2])**2 )
            sel = dist < G0_MR_centrals['Rvir'][jj] 
            #sel = G0_MR['FOFCentralGal'] == G0_MR_centrals['FOFCentralGal'][jj] 
            log_StellarMass[jj]=np.log10(np.sum(G0_MR['StellarMass'][sel]*1.e10/Hubble_h, axis=0) + 
                                         G0_MR['ICM'][jj]*1.e10/Hubble_h)            
            log_HaloMass[jj]=np.log10(G0_MR_centrals['Mvir'][jj]*1.e10/Hubble_h)             
            log_fraction[jj]=log_StellarMass[jj]-log_HaloMass[jj]
        
       
        bin=0.25
        (x_binned, median, mean, pc16, pc84,rms)= median_and_percentiles(bin, xlim[0], xlim[1]+1.,log_HaloMass,log_fraction)  
        subplot.plot(x_binned, median,color='black', linewidth=2, linestyle=':') 
               
        x_axis=x_binned
        y_axis=median
        fa = open(Datadir+"SMHM_"+model_to_print+"_median_with_sats_z0.0.txt", "w")
        fa.write("%d\n" % len(x_axis))
        for kk in range (0,len(x_axis)):               
            fa.write("%0.2f " % x_axis[kk] + "%0.2f\n" % y_axis[kk])     
        fa.close()'''
        
        
        
        
        
        #only satellites and ICM
        '''for jj in range (0,len(G0_MR)):
            dist=np.sqrt( (G0_MR['Pos'][:,0]-G0_MR['Pos'][jj,0])**2 +
                       (G0_MR['Pos'][:,1]-G0_MR['Pos'][jj,1])**2 + 
                       (G0_MR['Pos'][:,1]-G0_MR['Pos'][jj,1])**2 )
            sel=(dist<G0_MR['Rvir'][jj]) & (dist>0.)           
            log_StellarMass[jj]=np.log10(np.sum(G0_MR['StellarMass'][sel]*1.e10/Hubble_h, axis=0) + 
                                         G0_MR['ICM'][jj]*1.e10/Hubble_h)          
            #log_StellarMass[jj]=np.log10(np.sum(G0_MR['StellarMass'][sel]*1.e10/Hubble_h, axis=0))       
            log_HaloMass[jj]=np.log10(G0_MR['Mvir'][jj]*1.e10/Hubble_h)             
            log_fraction[jj]=log_StellarMass[jj]-log_HaloMass[jj]       
        (x_binned, median, mean, pc16, pc84,rms)= median_and_percentiles(bin, xlim[0], xlim[1],log_HaloMass,log_fraction)  
        subplot.plot(x_binned, median,color='black', linewidth=2, linestyle=':') '''
        
        
    
    #Behroozi2013
    z=ThisRedshiftList[0]   
    a=1/(1+z)
    neu=np.exp(-4*a*a)
    log_epsilon=-1.777+(-0.006*(a-1)-0.*z)*neu-0.119*(a-1)
    log_M1=11.514+(-1.793*(a-1)-0.251*z)*neu
    alpha=-1.412+0.731*(a-1.)*neu
    delta=3.508+(2.608*(a-1)-0.043*z)*neu
    gamma=0.316+(1.319*(a-1)+0.279*z)*neu 
           
    x=0.
    f_0=-np.log10(10**(alpha*x) + 1) + delta * ((np.log10(1 + np.exp(x)))**gamma)/(1 + np.exp(10**(-x)))             
    log_Mh=np.arange(8.0,16.0,0.01) 
    x=log_Mh-log_M1
    f_log_Mh=-np.log10(10**(alpha*x) + 1) + delta * ((np.log10(1 + np.exp(x)))**gamma)/(1 + np.exp(10**(-x)))      
    log_mstar=log_epsilon+log_M1+f_log_Mh-f_0     
    subplot.plot(log_Mh, log_mstar-log_Mh,color='limegreen', linewidth=3)
    
    plt.tight_layout()
    plt.savefig('./fig/plots_smhm_fractional_single_plot.pdf')
    #pdf.savefig()
    #plt.close()
    
    
    
    #PLOTS     

    #All models
    fig = plt.figure(figsize=(one_four_size_large[0],one_four_size_large[1]))
    grid = gridspec.GridSpec(1, 4)
    grid.update(wspace=0.0, hspace=0.0)

    
    model_names=['Hen15_no_feedback','Hen15_only_AGN','Hen15_only_SN','Hen15']
    #model_names=['Hen15_other','Hen15_only_AGN','Hen15_only_SN','Hen15']
    model_label=['no feedback','only AGN','only SN','Hen15']
    
    for ii in range(0,4):
        
        subplot=plt.subplot(grid[ii])

        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        xlab='$\log_{10}(M_{\mathrm{200c}}/$'+Msun+'$)$' 
        subplot.set_xlabel(xlab, fontsize=14 ,fontweight='bold')
      
        if ii==0:
            ylab='$\log_{10}(M_*/M_{\mathrm{200c}})$'
            subplot.set_ylabel(ylab, fontsize=14)
                       
    
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
        subplot.yaxis.set_minor_locator(MultipleLocator(0.25))
  
        if ii>0:
            plt.tick_params(axis='y', which='both', left='on', labelleft='off')
      
    
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
        
        #MODELS       
        fa = Datadir+"SMHM_"+model_names[ii]+"_SF_z0.0.txt"
        (x_axis,y_axis)=read_file(fa)    
        subplot.scatter(x_axis,y_axis,s=5, color='blue') 
        
        fa = Datadir+"SMHM_"+model_names[ii]+"_passive_z0.0.txt"
        (x_axis,y_axis)=read_file(fa)    
        subplot.scatter(x_axis,y_axis,s=5, alpha=0.8, color='red') 
            
        fa = Datadir+"SMHM_"+model_names[ii]+"_median_z0.0.txt"
        (x_axis,y_axis)=read_file(fa)
        subplot.plot(x_axis, y_axis,color='black', linewidth=2)   
         
        fa = Datadir+"SMHM_"+model_names[ii]+"_pc16_z0.0.txt"
        (x_axis,y_axis)=read_file(fa)
        subplot.plot(x_axis, y_axis,color='black', linewidth=2, linestyle='--') 
         
        fa = Datadir+"SMHM_"+model_names[ii]+"_pc84_z0.0.txt"
        (x_axis,y_axis)=read_file(fa)
        subplot.plot(x_axis, y_axis,color='black', linewidth=2, linestyle='--')  
    
        x_arr=np.arange(xlim[0],xlim[1]+0.05,0.05)           
        y_arr=x_arr*0.+np.log10(0.155)      
        subplot.plot(x_arr,y_arr,color='black', linestyle='--', linewidth=3)  
        
        fa = Datadir+"SMHM_"+model_names[ii]+"_median_with_sats_z0.0.txt"
        (x_axis,y_axis)=read_file(fa)       
        subplot.plot(x_axis,y_axis,color='black', linewidth=2, linestyle=':') 
         
        if(ii==0):
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.15, y_percentage=0.24, 
                        color='black', xlog=0, ylog=0, label='Model Median', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.06, y_percentage=0.26, x2_percentage=0.12, 
                        color='black', xlog=0, ylog=0, linestyle='-', linewidth=2)  
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.15, y_percentage=0.17, 
                        color='black', xlog=0, ylog=0, label='Centrals + sats + ICL', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.06, y_percentage=0.19, x2_percentage=0.12, 
                        color='black', xlog=0, ylog=0, linestyle=':', linewidth=2)  
            
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
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()
    
    return 




def stellar_mass_vs_halo_mass_fractional_allz(ThisRedshiftList):
 
    model_name = 'Hen15'
       
    plot_color = ['black', 'red', 'orange', 'green']   
    
        
    xlim=[11., 15.0]
    ylim=[-3.5,0.]    
    
    
    #PLOT current model
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1.0))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
    subplot.yaxis.set_minor_locator(MultipleLocator(0.5)) 
    
    xlab='$\log_{10}(M_{200c}/$'+Msun+'$)$' 
    ylab='$\log_{10}(M_*/M_{200c})$'
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14) 

    for ii in range (0,len(ThisRedshiftList)):
  
        char_redshift="%0.1f" % ThisRedshiftList[ii]

        #MODEL
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]             
        G0_MR=G0_MR[(G0_MR['StellarMass']>0.0) & (G0_MR['Type']==0)]  
        
        log_StellarMass = stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])  
        #log_StellarMass = np.log10(G0_MR['StellarMass']*1.e10/Hubble_h) 
        log_HaloMass=np.log10(G0_MR['Mvir']*1.e10/Hubble_h)  
        log_fraction=log_StellarMass-log_HaloMass
        
        bin=0.25
        (x_binned, median, mean, pc16, pc84,rms)= median_and_percentiles(bin, xlim[0], xlim[1],log_HaloMass,log_fraction) 
        sel = median<0.
        x_axis=x_binned[sel]
        y_axis=median[sel] 
        
        #subplot.plot(x_axis, y_axis,color='black', linewidth=2, linestyle='--') 
        
        fa = open(Datadir+"SMHM_allz_"+model_name+"_median_and_per_z"+char_redshift+".txt", "w")
        fa.write("%d\n" % len(x_axis))
        for kk in range (0,len(x_axis)):               
            fa.write("%0.2f " % x_axis[kk] + "%0.2f " % y_axis[kk]  + "%0.2f " % pc16[sel][kk]  + "%0.2f\n" % pc84[sel][kk])     
        fa.close() 
        
        fa = open(Datadir+"SMHM_allz_"+model_name+"_median_z"+char_redshift+".txt", "w")
        fa.write("%d\n" % len(x_axis))
        for kk in range (0,len(x_axis)):               
            fa.write("%0.2f " % x_axis[kk] + "%0.2f\n" % y_axis[kk])     
        fa.close() 
        
        y_axis = pc16[sel]
        fa = open(Datadir+"SMHM_allz_"+model_name+"_pc16_z"+char_redshift+".txt", "w")
        fa.write("%d\n" % len(x_axis))
        for kk in range (0,len(x_axis)):               
            fa.write("%0.2f " % x_axis[kk] + "%0.2f\n" % y_axis[kk])     
        fa.close() 
        
        y_axis = pc84[sel]
        fa = open(Datadir+"SMHM_allz_"+model_name+"_pc84_z"+char_redshift+".txt", "w")
        fa.write("%d\n" % len(x_axis))
        for kk in range (0,len(x_axis)):               
            fa.write("%0.2f " % x_axis[kk] + "%0.2f\n" % y_axis[kk])     
        fa.close() 
             
        
        
    
        #Behroozi2013
        z=ThisRedshiftList[ii]   
        a=1/(1+z)
        neu=np.exp(-4*a*a)
        log_epsilon=-1.777+(-0.006*(a-1)-0.*z)*neu-0.119*(a-1)
        log_M1=11.514+(-1.793*(a-1)-0.251*z)*neu
        alpha=-1.412+0.731*(a-1.)*neu
        delta=3.508+(2.608*(a-1)-0.043*z)*neu
        gamma=0.316+(1.319*(a-1)+0.279*z)*neu 
           
        x=0.
        f_0=-np.log10(10**(alpha*x) + 1) + delta * ((np.log10(1 + np.exp(x)))**gamma)/(1 + np.exp(10**(-x)))             
        log_Mh=np.arange(8.0,16.0,0.01) 
        x=log_Mh-log_M1
        f_log_Mh=-np.log10(10**(alpha*x) + 1) + delta * ((np.log10(1 + np.exp(x)))**gamma)/(1 + np.exp(10**(-x)))      
        log_mstar=log_epsilon+log_M1+f_log_Mh-f_0     
        subplot.plot(log_Mh, log_mstar-log_Mh,color=plot_color[ii], linewidth=2, linestyle='--')
    
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
      
        
        #MODELS   
        fa = Datadir+"SMHM_allz_"+model_name+"_median_z"+char_redshift+".txt"
        (x_axis,y_axis)=read_file(fa)
        fa = Datadir+"SMHM_allz_"+model_name+"_pc16_z"+char_redshift+".txt"
        (x_axis,pc16)=read_file(fa)
        fa = Datadir+"SMHM_allz_"+model_name+"_pc84_z"+char_redshift+".txt"
        (x_axis,pc84)=read_file(fa)
        
        #subplot.plot(x_axis, y_axis, color=plot_color[ii], linewidth=2)   
        #subplot.plot(x_axis, pc16, color=plot_color[ii], linewidth=2)   
        #subplot.plot(x_axis, pc84, color=plot_color[ii], linewidth=2) 
        error_down = y_axis-pc16
        error_up = pc84-y_axis
        asy_yerror = [error_down,error_up]
        subplot.errorbar(x_axis+ii*0.05, y_axis, yerr=asy_yerror, color=plot_color[ii], linewidth=2)   
           
        x_arr=np.arange(xlim[0],xlim[1]+0.05,0.05)           
        y_arr=x_arr*0.+np.log10(0.155)      
        subplot.plot(x_arr,y_arr,color='black', linestyle='--', linewidth=2)  
        
    '''fa = Datadir+"SMHM_"+model_names[ii]+"_median_with_sats_z0.0.txt"
    (x_axis,y_axis)=read_file(fa)       
    subplot.plot(x_axis,y_axis,color='black', linewidth=2, linestyle=':') '''
         
    '''if(ii==0):
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.15, y_percentage=0.24, 
                    color='black', xlog=0, ylog=0, label='Model Median', fontsize=10, fontweight='normal') 
        plot_label (subplot, 'line', xlim, ylim, x_percentage=0.06, y_percentage=0.26, x2_percentage=0.12, 
                    color='black', xlog=0, ylog=0, linestyle='-', linewidth=2)  
            
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.15, y_percentage=0.17, 
                    color='black', xlog=0, ylog=0, label='Centrals + sats + ICL', fontsize=10, fontweight='normal') 
        plot_label (subplot, 'line', xlim, ylim, x_percentage=0.06, y_percentage=0.19, x2_percentage=0.12, 
                    color='black', xlog=0, ylog=0, linestyle=':', linewidth=2)  
            
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.15, y_percentage=0.10, 
                    color='black', xlog=0, ylog=0, label='Behroozi et al. (2013)', fontsize=10, fontweight='normal') 
        plot_label (subplot, 'line', xlim, ylim, x_percentage=0.06, y_percentage=0.12, x2_percentage=0.12, 
                    color='limegreen', xlog=0, ylog=0, linestyle='-', linewidth=2) 
        
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.85, y_percentage=0.6, 
                color='black', xlog=0, ylog=0, label='$f_b$', 
                fontsize=15, fontweight='normal')
        
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.1, y_percentage=0.85, 
                color='black', xlog=0, ylog=0, label=model_name[ii], 
                fontsize=12, fontweight='normal')'''
        
    #endfor
    
    
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()
    
    return 
    
#end stellar_mass_vs_halo_mass_fractional_allz    
    
    
    
def all_masses_vs_halo_mass_fractional_allz(ThisRedshiftList):

    xlim=[11.0, 14.5]
    ylim=[-3.5,0.]    
        
        
    model_to_print='Hen15_other'    
    #TO WRITE OUT and plot FOR A GIVEN MODEL
    #PLOT for all models IS ACTUALLY MADE FROM FILES PREVIOUSLY READ
    
    
    #PLOT current model
    fig = plt.figure(figsize=(two_two_size_small[0],two_two_size_small[1]))
    grid = gridspec.GridSpec(2, 2)
    grid.update(wspace=0.0, hspace=0.0)

    for ii in range (0,len(ThisRedshiftList)):   
  
        char_redshift="%0.1f" % ThisRedshiftList[ii]

        subplot=plt.subplot(grid[ii])
        #subplot=plt.subplot()
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1.0))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.1)) 
        subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
    
       
        plt.tick_params(axis='x', which='both', bottom='on', labelbottom='on', top='on', labeltop='off')
        plt.tick_params(axis='y', which='both', left='on', labelleft='on', right='on', labelright='off')

        if (ii==0) or (ii==1):
            plt.tick_params(axis='x', which='both', bottom='on', labelbottom='off')
            xlab='' 
            ylab='$\log_{10}(M_*/M_{200c})$'
        if (ii==1) or (ii==3):
            plt.tick_params(axis='y', which='both', left='on', labelleft='off')
            xlab='$\log_{10}(M_{200c}/$'+Msun+'$)$' 
            ylab=''
            
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14) 
        
        #MEDIAN
        bin=0.25
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]             
        G0_MR=G0_MR[(G0_MR['StellarMass']>0.0) & (G0_MR['Type']==0)]       
        #G0_MR=G0_MR[(G0_MR['StellarMass']>0.0)]       
        log_StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii]) 
        log_BHMass=np.log10(G0_MR['BlackHoleMass']*1.e10/Hubble_h)  
        ColdGas=G0_MR['ColdGas']*1.e10/Hubble_h
        HotGas=G0_MR['HotGas']*1.e10/Hubble_h 
        EjectedGas=G0_MR['EjectedMass']*1.e10/Hubble_h  
        log_HaloMass=np.log10(G0_MR['Mvir']*1.e10/Hubble_h)  
        log_totBaryons=np.log10((G0_MR['BlackHoleMass']+G0_MR['ColdGas']+G0_MR['StellarMass']+
                                G0_MR['HotGas']+G0_MR['EjectedMass']+G0_MR['ICM'])*1.e10/Hubble_h)  
        
             
        
        log_fraction=log_BHMass-log_HaloMass      
        (x_binned, median, mean, pc16, pc84,rms)= median_and_percentiles(bin,xlim[0],xlim[1],log_HaloMass,log_fraction)  
        sel=median<0
        subplot.plot(x_binned[sel], median[sel],color='brown', linewidth=2, linestyle='--')
        
        log_fraction=log_StellarMass-log_HaloMass      
        (x_binned, median, mean, pc16, pc84,rms)= median_and_percentiles(bin,xlim[0],xlim[1],log_HaloMass,log_fraction) 
        sel=median<0
        subplot.plot(x_binned[sel], median[sel],color='orange', linewidth=2)        
                
        fraction=ColdGas/10**log_HaloMass        
        (x_binned, median, mean, pc16, pc84,rms)= median_and_percentiles(bin,xlim[0],xlim[1],log_HaloMass,fraction)        
        subplot.plot(x_binned, np.log10(median), color='blue', linewidth=2)
        
        fraction=HotGas/10**log_HaloMass        
        (x_binned, median, mean, pc16, pc84,rms)= median_and_percentiles(bin,xlim[0],xlim[1],log_HaloMass,fraction)       
        subplot.plot(x_binned, np.log10(median), color='red', linewidth=2)
                
        fraction=EjectedGas/10**log_HaloMass      
        (x_binned, median, mean, pc16, pc84,rms)= median_and_percentiles(bin,xlim[0],xlim[1],log_HaloMass,fraction)        
        subplot.plot(x_binned, np.log10(median),color='purple', linewidth=2)
        
        log_fraction=log_totBaryons-log_HaloMass      
        (x_binned, median, mean, pc16, pc84,rms)= median_and_percentiles(bin,xlim[0],xlim[1],log_HaloMass,log_fraction)  
        sel=median<0
        subplot.plot(x_binned[sel], median[sel],color='black', linewidth=2)  
    
        x_arr=np.arange(xlim[0],xlim[1]+0.05,0.05)           
        y_arr=x_arr*0.+np.log10(0.155)      
        subplot.plot(x_arr,y_arr,color='black', linestyle='--', linewidth=3)    
        
        
        
        #SATELLITES
        bin=0.25
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]             
        G0_MR_centrals=G0_MR[(G0_MR['StellarMass']>0.0) & (G0_MR['Type']==0)]       
        G0_MR_centrals=np.random.choice(G0_MR, size=10000)   
        
        log_ColdandStars=np.zeros(len(G0_MR_centrals),dtype=np.float32)
        log_HaloMass=np.zeros(len(G0_MR_centrals),dtype=np.float32)
        log_fraction=np.zeros(len(G0_MR_centrals),dtype=np.float32)
        
        '''for jj in range (0,len(G0_MR_centrals)):
            dist=np.sqrt( (G0_MR['Pos'][:,0]-G0_MR_centrals['Pos'][jj,0])**2 +
                       (G0_MR['Pos'][:,1]-G0_MR_centrals['Pos'][jj,1])**2 + 
                       (G0_MR['Pos'][:,2]-G0_MR_centrals['Pos'][jj,2])**2 )
            sel=(dist<G0_MR_centrals['Rvir'][jj]) & (dist>0.) & (G0_MR['StellarMass']>0.) & (G0_MR['ColdGas']>0.)         
            if(len(G0_MR['StellarMass'][sel])>0.):     
                #print(np.log10(G0_MR_centrals['Mvir'][jj]*1.e10/Hubble_h), 
                #      np.log10(G0_MR_centrals['StellarMass'][jj]*1.e10/Hubble_h))
                #print(np.log10(G0_MR['StellarMass'][sel]*1.e10/Hubble_h))      
                
                log_ColdandStars[jj]=np.log10(np.sum((G0_MR['StellarMass'][sel]+
                                                      G0_MR['ColdGas'][sel])*1.e10/Hubble_h, axis=0))            
                log_HaloMass[jj]=np.log10(G0_MR_centrals['Mvir'][jj]*1.e10/Hubble_h)             
                log_fraction[jj]=log_ColdandStars[jj]-log_HaloMass[jj]             
            else:
                log_HaloMass[jj]=np.log10(G0_MR_centrals['Mvir'][jj]*1.e10/Hubble_h)             
                log_fraction[jj]=0.0       
                
            #if((log_HaloMass[jj]<11.5) & (log_HaloMass[jj]>11.0)):
            #    print(log_ColdandStars[jj])
        
        (x_binned, median, mean, pc16, pc84,rms)= median_and_percentiles(bin,xlim[0],xlim[1],log_HaloMass,log_fraction) 
        sel=median<0
        subplot.plot(x_binned[sel], median[sel],color='blue', linestyle='--', linewidth=2)'''  
               
        label='z='+char_redshift
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.86, 
                    color='black', xlog=0, ylog=0, label=label, fontsize=13, fontweight='normal')
        
        if(ii==3):
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.72, y_percentage=0.10, 
                        color='black', xlog=0, ylog=0, label='BH', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.6, y_percentage=0.12, x2_percentage=0.7, 
                        color='brown', xlog=0, ylog=0, linestyle='--', linewidth=2)
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.72, y_percentage=0.20, 
                        color='black', xlog=0, ylog=0, label='Stars', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.6, y_percentage=0.22, x2_percentage=0.7, 
                        color='orange', xlog=0, ylog=0, linestyle='-', linewidth=2)
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.72, y_percentage=0.30, 
                        color='black', xlog=0, ylog=0, label='Cold', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.6, y_percentage=0.32, x2_percentage=0.7, 
                        color='blue', xlog=0, ylog=0, linestyle='-', linewidth=2)
            
            '''plot_label (subplot, 'label', xlim, ylim, x_percentage=0.72, y_percentage=0.40, 
                        color='black', xlog=0, ylog=0, label='Sats (C+S)', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.6, y_percentage=0.42, x2_percentage=0.7, 
                        color='blue', xlog=0, ylog=0, linestyle='--', linewidth=2)'''
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.72, y_percentage=0.50, 
                        color='black', xlog=0, ylog=0, label='Hot', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.6, y_percentage=0.52, x2_percentage=0.7, 
                        color='red', xlog=0, ylog=0, linestyle='-', linewidth=2)
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.72, y_percentage=0.60, 
                        color='black', xlog=0, ylog=0, label='Eject', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.6, y_percentage=0.62, x2_percentage=0.7, 
                        color='purple', xlog=0, ylog=0, linestyle='-', linewidth=2)
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.72, y_percentage=0.70, 
                        color='black', xlog=0, ylog=0, label='Tot_bar', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.6, y_percentage=0.72, x2_percentage=0.7, 
                        color='black', xlog=0, ylog=0, linestyle='-', linewidth=2)
        
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')   
    plt.close()    
    
    return   
#end all_masses_vs_halo_mass_fractional_allz






def eject_hot_masses_vs_halo_mass_fractional_allz(ThisRedshiftList):

    xlim=[11.0, 14.5]
    ylim=[-4.,0.]    
            
    #PLOT current model
    fig = plt.figure(figsize=(two_two_size_small[0],two_two_size_small[1]))
  
    plt.tick_params(axis='x', which='both', reset=1) 
    plt.tick_params(axis='y', which='both', reset=1) 

    
    for ii in range (0,len(ThisRedshiftList)):
  
        char_redshift="%0.1f" % ThisRedshiftList[ii]

        #subplot=plt.subplot(grid[ii])
        subplot=plt.subplot()
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(0.5))  
        subplot.yaxis.set_major_locator(MultipleLocator(1.0)) 
        subplot.xaxis.set_minor_locator(MultipleLocator(0.1)) 
        subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
              
        #plt.tick_params(axis='x', which='both', bottom='on', labelbottom='on')
        #plt.tick_params(axis='x', which='both', top='on', labeltop='on')

        if (ii==0) or (ii==1):
            #plt.tick_params(axis='x', which='both', bottom='on', labelbottom='off')
            xlab='' 
            ylab='$\log_{10}(M_*/M_{200c})$'
        if (ii==1) or (ii==3):
            #plt.tick_params(axis='y', which='both', left='on', labelleft='off')
            xlab='$\log_{10}(M_{200c}/$'+Msun+'$)$' 
            ylab=''
            
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14) 
        
        plt_line=['-', '--', '-.', ':']
        
        #MEDIAN
        bin=0.25
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]      
        HotGas=G0_MR['HotGas']*1.e10/Hubble_h 
        EjectedGas=G0_MR['EjectedMass']*1.e10/Hubble_h  
        log_HaloMass=np.log10(G0_MR['Mvir']*1.e10/Hubble_h)  
       
        fraction=HotGas/10**log_HaloMass        
        (x_binned, median, mean, pc16, pc84,rms)= median_and_percentiles(bin,xlim[0],xlim[1],log_HaloMass,fraction)       
        subplot.plot(x_binned, np.log10(median), color='red', linestyle=plt_line[ii],linewidth=2)
               
        fraction=EjectedGas/10**log_HaloMass      
        (x_binned, median, mean, pc16, pc84,rms)= median_and_percentiles(bin,xlim[0],xlim[1],log_HaloMass,fraction)        
        subplot.plot(x_binned, np.log10(median),color='purple', linestyle=plt_line[ii], linewidth=2)
        
        x_arr=np.arange(xlim[0],xlim[1]+0.05,0.05)           
        y_arr=x_arr*0.+np.log10(0.155)      
        subplot.plot(x_arr,y_arr,color='black', linestyle='--', linewidth=3)    
        
        
     
        label='z='+char_redshift
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.1, y_percentage=0.85+ii*0.05, 
                    color='black', xlog=0, ylog=0, label=label, fontsize=13, fontweight='normal')
        plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, y_percentage=0.86+ii*0.05, x2_percentage=0.09, 
                    color='red', xlog=0, ylog=0, linestyle=plt_line[ii], linewidth=2)
        
        
        if(ii==0):           
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.76, y_percentage=0.50, 
                        color='black', xlog=0, ylog=0, label='Hot', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.7, y_percentage=0.51, x2_percentage=0.75, 
                        color='red', xlog=0, ylog=0, linestyle='-', linewidth=2)
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.76, y_percentage=0.55, 
                        color='black', xlog=0, ylog=0, label='Eject', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.7, y_percentage=0.56, x2_percentage=0.75, 
                        color='purple', xlog=0, ylog=0, linestyle='-', linewidth=2)
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.76, y_percentage=0.60, 
                        color='black', xlog=0, ylog=0, label='Tot_bar', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.7, y_percentage=0.61, x2_percentage=0.75, 
                        color='black', xlog=0, ylog=0, linestyle='--', linewidth=2)
        
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close() 
    
    return   
    
    
    '''xlim=[0.,3.]
    ylim=[11.0, 14.5]    
        
        
    model_to_print='Hen15_other'    
    #TO WRITE OUT and plot FOR A GIVEN MODEL
    #PLOT for all models IS ACTUALLY MADE FROM FILES PREVIOUSLY READ
    
    
    #PLOT current model
    fig = plt.figure(figsize=(two_two_size_small[0],two_two_size_small[1]))
    
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1.0))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.1)) 
    subplot.yaxis.set_minor_locator(MultipleLocator(0.5)) 
    
       
    plt.tick_params(axis='x', which='both', bottom='on', labelbottom='on')
    plt.tick_params(axis='x', which='both', top='on', labeltop='on')
    ylab='$\log_{10}(M_*/M_{200c})$'
    xlab='$\log_{10}(M_{200c}/$'+Msun+'$)$' 
      
    time = np.array([13.721, 5.903, 3.316, 2.171])
    mvir = np.array([13.721, 5.903, 3.316, 2.171])
    mvir = 10.*3.*10**20/(mvir*1e9)
    redshift=np.array([0.0,1.0,2.0,3.0])
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)       
    subplot.plot(redshift, np.log10(mvir),color='brown', linewidth=2, linestyle='--')
    
    mvir = np.array([12.3,12.,11.9,11.7])
       
    subplot.plot(redshift, mvir,color='brown', linewidth=2, linestyle='--')
        
    pdf.savefig()
    plt.close()   '''
#end eject_hot_masses_vs_halo_mass_fractional_allz






    
def all_masses_vs_halo_mass_fractional_z0(ThisRedshiftList):

    xlim=[11.0, 14.5]
    ylim=[-6.,0.]    
      
    ii=0    
        
    #PLOT current model
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
  
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1.0))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
    subplot.yaxis.set_minor_locator(MultipleLocator(0.5)) 
    
    xlab='$\log_{10}(M_{200c}/$'+Msun+'$)$' 
    ylab='$\log_{10}(M/M_{200c})$'
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14) 
        
    #MEDIAN
    bin=0.1
    (sel)=select_current_redshift(G_MRII, ThisRedshiftList, ii, FullSnapshotList_MRII)        
    G0_MR=G_MRII[sel]             
    G0_MR=G0_MR[(G0_MR['StellarMass']>0.0) & (G0_MR['Type']==0)]       
    #G0_MR=G0_MR[(G0_MR['StellarMass']>0.0)]       
    log_HaloMass=np.log10((G0_MR['Mvir']*1.e10/Hubble_h))  
    StellarMass=10**stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii]) 
    BHMass=(G0_MR['BlackHoleMass']*1.e10/Hubble_h)  
    ColdGas=(G0_MR['ColdGas']*1.e10/Hubble_h)  
    HotGas=(G0_MR['HotGas']*1.e10/Hubble_h)  
    EjectedGas=(G0_MR['EjectedMass']*1.e10/Hubble_h)     
    totBaryons=((G0_MR['BlackHoleMass']+G0_MR['ColdGas']+G0_MR['StellarMass']+
                 G0_MR['HotGas']+G0_MR['EjectedMass']+G0_MR['ICM'])*1.e10/Hubble_h)  
    SFR=(G0_MR['Sfr'])                
    
           
    #fraction=totBaryons/10**log_HaloMass      
    #(x_binned, median, mean, pc16, pc84,rms)= median_and_percentiles(bin,xlim[0],xlim[1],log_HaloMass,fraction)   
    #subplot.plot(x_binned, np.log10(median),color='black', linewidth=2)  
    
    x_arr=np.arange(xlim[0],xlim[1]+0.05,0.05)           
    y_arr=x_arr*0.+np.log10(0.155)      
    subplot.plot(x_arr,y_arr,color='black', linestyle='--', linewidth=2) 
    
    fraction=EjectedGas/10**log_HaloMass      
    (x_binned, median, mean, pc16, pc84,rms)= median_and_percentiles(bin,xlim[0],xlim[1],log_HaloMass,fraction)
    subplot.fill_between(x_binned, np.log10(pc84), np.log10(pc16), interpolate=True, 
                         alpha=0.4, facecolor='purple', edgecolor='purple')  
    subplot.plot(x_binned, np.log10(median),color='purple', linewidth=2)  
    median_1 = np.log10(median)
    
    fraction=HotGas/10**log_HaloMass      
    (x_binned, median, mean, pc16, pc84,rms)= median_and_percentiles(bin,xlim[0],xlim[1],log_HaloMass,fraction)
    subplot.fill_between(x_binned, np.log10(pc84), np.log10(pc16), interpolate=True,
                         alpha=0.4,facecolor='red',edgecolor='red')  
    subplot.plot(x_binned, np.log10(median),color='red', linewidth=2)
    median_2 = np.log10(median)
    
    #SN switch off
    yy=np.arange(ylim[0]-0.5,ylim[1]+0.5,0.1)
    xx_1=np.zeros(len(yy))
    #xx_1+=np.amin(x_binned[median_1<median_2])  
    xx_1 += 11.665
    
    sel=G0_MR['ColdGas']>0.
    fraction=ColdGas/10**log_HaloMass      
    (x_binned, median, mean, pc16, pc84,rms)= median_and_percentiles(bin,xlim[0],xlim[1],log_HaloMass[sel],fraction[sel]) 
    subplot.fill_between(x_binned, np.log10(pc84), np.log10(pc16), interpolate=True,
                         alpha=0.4,facecolor='blue',edgecolor='blue')  
    subplot.plot(x_binned, np.log10(median),color='blue', linewidth=2)
    
    fraction=StellarMass/10**log_HaloMass      
    (x_binned, median, mean, pc16, pc84,rms)= median_and_percentiles(bin,xlim[0],xlim[1],log_HaloMass,fraction) 
    subplot.fill_between(x_binned, np.log10(pc84), np.log10(pc16), interpolate=True, 
                         alpha=0.4, facecolor='orange', edgecolor='orange') 
    subplot.plot(x_binned, np.log10(median),color='chocolate', linewidth=2)
    
    fraction=BHMass/10**log_HaloMass      
    (x_binned, median, mean, pc16, pc84,rms)= median_and_percentiles(bin,xlim[0],xlim[1],log_HaloMass,fraction)  
    sel = pc16<10**ylim[0]
    pc16[sel]=10**ylim[0]
    subplot.fill_between(x_binned, np.log10(pc84), np.log10(pc16), interpolate=True,
                         alpha=0.4, facecolor='brown', edgecolor='brown')  
    subplot.plot(x_binned, np.log10(median),color='brown', linewidth=2)
   
    #Mvir for quenching
    #sel = np.log10(SSFR)>-11.0
    SSFR=SFR/StellarMass     
    (x_binned, median, mean, pc16, pc84,rms)= median_and_percentiles(bin,xlim[0],xlim[1],log_HaloMass,SSFR)    
    yy=np.arange(ylim[0]-0.5,ylim[1]+0.5,0.1)
    xx_2=np.zeros(len(yy))
    Redshift=0.
    #xx_2+=np.amax(x_binned[np.log10(median)>np.log10((1+Redshift)/(2.*1.37e10))])  
    xx_2+=np.amax(x_binned[np.log10(median)>-11.0])  
    
    
    
    #region between SN switch off and quenching
    xx = np.arange(xlim[0]-0.5,xlim[1]+0.5,0.01)
    xx = xx[(xx>xx_1[0]) & (xx<xx_2[0])]
    yy_1 = np.zeros(len(xx)) + ylim[0]
    yy_2 = np.zeros(len(xx)) + ylim[1]     
    #subplot.fill_between(xx, yy_1, yy_2, interpolate=True,alpha=0.2, hatch='/', color='grey',
    #                     facecolor='grey', edgecolor='grey', linestyle = '--') 
    subplot.fill_between(xx, yy_1, yy_2, interpolate=True,alpha=0.4, color='grey', facecolor='grey') 
    yy = np.arange(ylim[0],ylim[1],0.01)
    xx = yy/yy-1.+ xx_1[0]
    #subplot.plot(xx,yy, linestyle=':', color='black', linewidth=1) 
    
    y = np.arange(ylim[0],ylim[1],0.01)
    xx = yy/yy-1.+ xx_2[0]
    #subplot.plot(xx,yy, linestyle=':', color='black', linewidth=1) 
    
    
    #peak SFR LINE
    (x_binned, median, mean, pc16, pc84,rms)= median_and_percentiles(bin,xlim[0],xlim[1],log_HaloMass,SFR)    
    yy=np.arange(ylim[0]-0.5,ylim[1]+0.5,0.01)
    xx=np.zeros(len(yy))    
    xx+=x_binned[median==np.max(median)]
    subplot.plot(xx,yy, linestyle='--', color='black', alpha=0.5) 
    
    #Mvir=12 line
    #yy=np.arange(ylim[0]-0.5,ylim[1]+0.5,0.01)
    #xx=np.zeros(len(yy))    
    #xx+=12.0    
    #subplot.plot(xx,yy, linestyle='--', color='black', alpha=0.5)   
    
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.11, y_percentage=0.31, 
                color='black', xlog=0, ylog=0, label='Tot_bar', fontsize=10, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.03, y_percentage=0.33, x2_percentage=0.09, 
                color='black', xlog=0, ylog=0, linestyle='--', linewidth=2)
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.11, y_percentage=0.26, 
                color='black', xlog=0, ylog=0, label='Eject', fontsize=10, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.03, y_percentage=0.28, x2_percentage=0.09, 
                color='purple', xlog=0, ylog=0, linestyle='-', linewidth=2)
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.11, y_percentage=0.21, 
                color='black', xlog=0, ylog=0, label='Hot', fontsize=10, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.03, y_percentage=0.23, x2_percentage=0.09, 
                color='red', xlog=0, ylog=0, linestyle='-', linewidth=2) 
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.11, y_percentage=0.16, 
                color='black', xlog=0, ylog=0, label='Cold', fontsize=10, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.03, y_percentage=0.18, x2_percentage=0.09, 
                color='blue', xlog=0, ylog=0, linestyle='-', linewidth=2) 
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.11, y_percentage=0.11, 
                color='black', xlog=0, ylog=0, label='Stars', fontsize=10, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.03, y_percentage=0.13, x2_percentage=0.09, 
                color='orange', xlog=0, ylog=0, linestyle='-', linewidth=2) 
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.11, y_percentage=0.06, 
                color='black', xlog=0, ylog=0, label='BH', fontsize=10, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.03, y_percentage=0.08, x2_percentage=0.09, 
                color='brown', xlog=0, ylog=0, linestyle='-', linewidth=2) 

    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HWL18_AMHM_fractional_z0.pdf')
    plt.close()  
   
    return 
#end all_masses_vs_halo_mass_fractional_z0    
    
    
    
    
def stellar_mass_vs_halo_mass_fractional_models_overplot(ThisRedshiftList):

    xlim=[11.0, 14.5]
    ylim=[-4.,0.]    
        
        
    #model_to_print='Hen15_ejection_cut_very_high_mass'  
    model_to_print='Hen15_other'
    model_names=['Hen15_ejection_cut_very_low_mass','Hen15_ejection_cut_low_mass',
                 'Hen15_ejection_cut_high_mass','Hen15_ejection_cut_very_high_mass','Hen15']
    #TO WRITE OUT and plot FOR A GIVEN MODEL
    #PLOT for all models IS ACTUALLY MADE FROM FILES PREVIOUSLY READ
    
    
    #PLOT current model
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
   
    for ii in range (0,len(ThisRedshiftList)):
  
        char_redshift="%0.1f" % ThisRedshiftList[ii]

        subplot=plt.subplot()
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1.0))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
        subplot.yaxis.set_minor_locator(MultipleLocator(0.5)) 
    
        ylab='$\log_{10}(M_*/M_{200c})$'
        xlab='$\log_{10}(M_{200c}/$'+Msun+'$)$' 
     
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14) 
        
        #MEDIAN
        bin=0.25
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]             
        G0_MR=G0_MR[(G0_MR['StellarMass']>0.0) & (G0_MR['Type']==0)]                
        log_StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii]) 
       
        log_HaloMass=np.log10(G0_MR['Mvir']*1.e10/Hubble_h)  
                   
        
        log_fraction=log_StellarMass-log_HaloMass      
        (x_binned, median, mean, pc16, pc84,rms)= median_and_percentiles(bin,xlim[0],xlim[1],log_HaloMass,log_fraction)
        sel=median<0
        subplot.plot(x_binned[sel], median[sel],color='orange', linewidth=2,linestyle='--')        
                
      
        x_axis=x_binned[sel]
        y_axis=median[sel]
        fa = open(Datadir+"SMHM_models_overplot_"+model_to_print+"_median_z0.0.txt", "w")
        fa.write("%d\n" % len(x_axis))
        for kk in range (0,len(x_axis)):               
            fa.write("%0.2f " % x_axis[kk] + "%0.2f\n" % y_axis[kk])     
        fa.close() 
        
        #MODELS       
        fa = Datadir+"SMHM_models_overplot_"+model_names[0]+"_median_z0.0.txt"
        (x_axis,y_axis)=read_file(fa)            
        subplot.plot(x_axis,y_axis,color='orange', linewidth=2)    
        
        fa = Datadir+"SMHM_models_overplot_"+model_names[1]+"_median_z0.0.txt"
        (x_axis,y_axis)=read_file(fa)            
        subplot.plot(x_axis,y_axis,color='orange', linewidth=2)    
        
        fa = Datadir+"SMHM_models_overplot_"+model_names[2]+"_median_z0.0.txt"
        (x_axis,y_axis)=read_file(fa)            
        subplot.plot(x_axis,y_axis,color='orange', linewidth=2)   
        
        fa = Datadir+"SMHM_models_overplot_"+model_names[3]+"_median_z0.0.txt"
        (x_axis,y_axis)=read_file(fa)            
        subplot.plot(x_axis,y_axis,color='orange', linewidth=2)   
        
        fa = Datadir+"SMHM_models_overplot_"+model_names[4]+"_median_z0.0.txt"
        (x_axis,y_axis)=read_file(fa)            
        subplot.plot(x_axis,y_axis,color='orange', linewidth=2)   
    
    
        x_arr=np.arange(xlim[0],xlim[1]+0.05,0.05)           
        y_arr=x_arr*0.+np.log10(0.155)      
        subplot.plot(x_arr,y_arr,color='black', linestyle='--', linewidth=3)    
        
       
       
               
        #label='z='+char_redshift
        #plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.86, 
        #            color='black', xlog=0, ylog=0, label=label, fontsize=13, fontweight='normal')
        
        if(ii==3):
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.72, y_percentage=0.10, 
                        color='black', xlog=0, ylog=0, label='BH', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.6, y_percentage=0.12, x2_percentage=0.7, 
                        color='brown', xlog=0, ylog=0, linestyle='--', linewidth=2)
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.72, y_percentage=0.20, 
                        color='black', xlog=0, ylog=0, label='Stars', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.6, y_percentage=0.22, x2_percentage=0.7, 
                        color='orange', xlog=0, ylog=0, linestyle='-', linewidth=2)
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.72, y_percentage=0.30, 
                        color='black', xlog=0, ylog=0, label='Cold', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.6, y_percentage=0.32, x2_percentage=0.7, 
                        color='blue', xlog=0, ylog=0, linestyle='-', linewidth=2)
            
            '''plot_label (subplot, 'label', xlim, ylim, x_percentage=0.72, y_percentage=0.40, 
                        color='black', xlog=0, ylog=0, label='Sats (C+S)', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.6, y_percentage=0.42, x2_percentage=0.7, 
                        color='blue', xlog=0, ylog=0, linestyle='--', linewidth=2)'''
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.72, y_percentage=0.50, 
                        color='black', xlog=0, ylog=0, label='Hot', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.6, y_percentage=0.52, x2_percentage=0.7, 
                        color='red', xlog=0, ylog=0, linestyle='-', linewidth=2)
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.72, y_percentage=0.60, 
                        color='black', xlog=0, ylog=0, label='Eject', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.6, y_percentage=0.62, x2_percentage=0.7, 
                        color='brown', xlog=0, ylog=0, linestyle='-', linewidth=2)
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.72, y_percentage=0.70, 
                        color='black', xlog=0, ylog=0, label='Tot_bar', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.6, y_percentage=0.72, x2_percentage=0.7, 
                        color='black', xlog=0, ylog=0, linestyle='-', linewidth=2)
        
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HWL18_SMHM_fractional_models_overplot.pdf')
    plt.close()    
    
    return  
    
 



def cumulative_ND(ThisRedshiftList):
 
    model_to_print = 'Hen15'
       
    plot_color = ['black', 'red', 'orange', 'green']   
    
        
    xlim=[0.0, 3.5]
    ylim=[-6.5,-2.5]    
           
    #PLOT current model
    #fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    #subplot=plt.subplot()
    #subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    fig = plt.figure(figsize=(two_one_size_small[0],two_one_size_small[1]))
    grid = gridspec.GridSpec(2, 1)
    grid.update(wspace=0.0, hspace=0.0)   
    
    Hen15_redshift_list = [0.1,0.4,1.0,2.0,3.0]
    #log10(M*[h^-2 M_sun]) = [9.00, 10.25, 10.75]
    #log10(M*[M_sun]) = [9.3, 10.5, 11.]
    Hen15_ND_red_1 = np.array([-1.96059, -2.06783, -2.38238, -2.84467, -3.74531])  
    Hen15_ND_red_2 = np.array([-2.46364, -2.60637, -2.94759, -3.62935, -4.75871])
    Hen15_ND_red_3 = np.array([-3.29765, -3.48007, -3.79894, -4.46123, -5.54333])
    
    Hen15_ND_blue_1 = np.array([-1.81627, -1.80221, -1.82272, -2.03774, -2.31911])  
    Hen15_ND_blue_2 = np.array([-3.11823, -2.90210, -2.87914, -3.19085, -3.60414])
    Hen15_ND_blue_3 = np.array([-4.37175, -4.07072, -4.13094, -4.28870, -4.66282])
        
    OBS_ND_red_1      = np.array([-1.80343, -2.10752, -2.38171])   
    OBS_ND_red_high_1 = np.array([0.102433, 0.101904, 0.173990])   
    OBS_ND_red_low_1  = np.array([0.134404, 0.134143, 0.296018])  
    OBS_red_z = [0.1,0.4,1.0,2.0]
    OBS_ND_red_2      = np.array([-2.51148, -2.68871, -2.86200, -3.62600])
    OBS_ND_red_high_2 = np.array([0.104875, 0.156997, 0.166824, 0.226057])
    OBS_ND_red_low_2  = np.array([0.138750, 0.249925, 0.275766, 0.504285])
    OBS_ND_red_3      = np.array([-3.50930, -3.55753, -3.61682, -4.59453])
    OBS_ND_red_high_3 = np.array([0.097801, 0.190485, 0.201206, 0.262409])
    OBS_ND_red_low_3  = np.array([0.126428, 0.351342, 0.393114, 0.812666])
   
    OBS_ND_blue_1      = np.array([-1.81167, -1.86209, -1.97233, -2.17243]) 
    OBS_ND_blue_high_1 = np.array([0.0619133, 0.0756772, 0.126469, 0.137743]) 
    OBS_ND_blue_low_1  = np.array([0.0724885, 0.0923479, 0.179968, 0.203716]) 
    OBS_ND_blue_2      = np.array([-3.00867, -3.01968, -3.18696, -3.28282, -3.65009])
    OBS_ND_blue_high_2 = np.array([0.119301, 0.178726, 0.132090, 0.193525, 0.224108])
    OBS_ND_blue_low_2  = np.array([0.165320, 0.317151, 0.196590, 0.365221, 0.492930])
    OBS_blue_z_3 = [0.4,1.0,2.0,3.0]
    OBS_ND_blue_3      = np.array([-4.26583, -4.48239, -4.46429, -4.48222])
    OBS_ND_blue_high_3 = np.array([0.263206, 0.243436, 0.222771, 0.252961])
    OBS_ND_blue_low_3  = np.array([0.834312, 0.702773, 0.566281, 0.686049])
    
    eagle_z = [0.1, 1.0, 2.0, 3.0]
    eagle_passive = [2.470e-04, 4.400e-05, 3.000e-06, 0.0]
    eagle_SF = [3.700e-05, 1.100e-04, 2.700e-05, 1.000e-06]
    

    
    #subplot.plot(Hen15_redshift_list, Hen15_ND_red_1, color='red', linestyle=':')   
    #subplot.plot(Hen15_redshift_list, Hen15_ND_red_2+np.log10(Hubble_h**3), color='red', linestyle=':')
    #subplot.plot(Hen15_redshift_list, Hen15_ND_red_3+np.log10(Hubble_h**3), color='red', linestyle=':')
    
    #subplot.plot(Hen15_redshift_list, Hen15_ND_blue_1, color='blue', linestyle=':')
    #subplot.plot(Hen15_redshift_list, Hen15_ND_blue_2+np.log10(Hubble_h**3), color='blue', linestyle=':')
    #subplot.plot(Hen15_redshift_list, Hen15_ND_blue_3+np.log10(Hubble_h**3), color='blue', linestyle=':')
      
    
    Mass_limits = [10.56, 11.06]
    #Mass_limits = [11.06]
    Max_Mass = 11.66
    
    #Mass_limits = [10.25, 10.75]
    #Max_Mass = 11.375
    
    kk = 0
    for i_mass in Mass_limits:
    
        subplot=plt.subplot(grid[kk])    
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)            
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1.0))    
        subplot.xaxis.set_minor_locator(MultipleLocator(.1))  
        subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
        subplot.yaxis.set_minor_locator(MultipleLocator(0.1))  
            
        xlab='z' 
        ylab='$\log_{10}(n[\mathrm{Mpc}^{-3}])$'       
        if(kk==0):                 
            plt.tick_params(axis='x', which='both', bottom='on', labelbottom='off')
        else:
            subplot.set_xlabel(xlab, fontsize=14)               
        subplot.set_ylabel(ylab, fontsize=14)
        #subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)
        
        if kk == 0:
            subplot.fill_between(Hen15_redshift_list, OBS_ND_blue_2+np.log10(Hubble_h**3)-OBS_ND_blue_low_2, 
                                 OBS_ND_blue_2+np.log10(Hubble_h**3)+OBS_ND_blue_high_2,
                                 interpolate=True,alpha=0.3,facecolor='blue')
            subplot.fill_between(OBS_red_z, OBS_ND_red_2+np.log10(Hubble_h**3)-OBS_ND_red_low_2, 
                                 OBS_ND_red_2+np.log10(Hubble_h**3)+OBS_ND_red_high_2,
                                 interpolate=True,alpha=0.3,facecolor='red')
        else:                
            subplot.fill_between(OBS_blue_z_3, OBS_ND_blue_3+np.log10(Hubble_h**3)-OBS_ND_blue_low_3, 
                                 OBS_ND_blue_3+np.log10(Hubble_h**3)+OBS_ND_blue_high_3,
                                 interpolate=True,alpha=0.3,facecolor='blue')
            subplot.fill_between(OBS_red_z, OBS_ND_red_3+np.log10(Hubble_h**3)-OBS_ND_red_low_3, 
                                 OBS_ND_red_3+np.log10(Hubble_h**3)+OBS_ND_red_high_3,
                                 interpolate=True,alpha=0.3,facecolor='red')
        kk+=1
                    
        cumulative_ND_red = np.zeros(len(ThisRedshiftList),dtype=np.float32)
        cumulative_ND_blue = np.zeros(len(ThisRedshiftList),dtype=np.float32)
    
    
        for ii in range (0,len(ThisRedshiftList)):
  
            char_redshift="%0.1f" % ThisRedshiftList[ii]
            char_Mass="%0.1f" % i_mass

            #MODEL
            (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
            G0_MR=G_MR[sel]             
            G0_MR=G0_MR[(G0_MR['StellarMass']>0.0) & (G0_MR['Type']==0)]  
            log_StellarMass = stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])#+np.log10(Hubble_h**2)          
            log_SFR=np.log10(G0_MR['Sfr'])
        
            #z=0.          
            if ii==0:
                color_ur=G0_MR['MagDust'][:,15]-G0_MR['MagDust'][:,17]                  
                Magr=G0_MR['MagDust'][:,17]                
                sel_red=( (color_ur>(offset_color_cut[ii]-slope_color_cut[ii]*np.tanh((Magr+18.07)/1.09))))  
                sel_blue=( (color_ur<(offset_color_cut[ii]-slope_color_cut[ii]*np.tanh((Magr+18.07)/1.09))))  
            #z>0.
            else:
                color_UV=G0_MR['MagDust'][:,0]-G0_MR['MagDust'][:,2]  
                color_VJ=G0_MR['MagDust'][:,2]-G0_MR['MagDust'][:,7]               
                sel_red = np.logical_and(color_UV > minimum_y_color_cut[ii], 
                                         color_UV > color_VJ*slope_color_cut[ii] + offset_color_cut[ii])
                sel_blue = np.logical_or(color_UV < minimum_y_color_cut[ii], 
                                         color_UV < color_VJ*slope_color_cut[ii] + offset_color_cut[ii])
              
              
                  
    
            sel_red = ( (log_SFR-log_StellarMass) < np.log10((1+ThisRedshiftList[ii])/(2.*1.37e10)) ) 
            model_MR = log_StellarMass[sel_red]
            model_MR = model_MR[(model_MR>i_mass) & (model_MR<Max_Mass)]
            weight_MR=model_MR/model_MR/(Volume_MR)
        
            bin=0.1      
            x_bin=np.arange(7.0, 12.0, bin)
            hist = np.zeros(len(x_bin),dtype=np.float32)
                     
            for jj in range(0, len(x_bin)):
                sel=(model_MR >= x_bin[jj]-bin/2.) & (model_MR <= x_bin[jj]+bin/2.)
                if(len(weight_MR[sel])>0.):
                    hist[jj] = np.sum(weight_MR[sel])/bin
          
            #print(np.log10(simps(hist, x=x_bin)) )
            cumulative_ND_red[ii] = np.log10(simps(hist, x=x_bin)) 
                    
            sel_blue = ( (log_SFR-log_StellarMass) > np.log10((1+ThisRedshiftList[ii])/(2.*1.37e10)) ) 
            model_MR=log_StellarMass[sel_blue]
            model_MR = model_MR[(model_MR>i_mass) & (model_MR<Max_Mass)]
            weight_MR=model_MR/model_MR/(Volume_MR)
        
            bin=0.1      
            x_bin=np.arange(7.0, 12.0, bin)
            hist = np.zeros(len(x_bin),dtype=np.float32)
                         
            for jj in range(0, len(x_bin)):
                sel=(model_MR >= x_bin[jj]-bin/2.) & (model_MR <= x_bin[jj]+bin/2.)
                if(len(weight_MR[sel])>0.):
                    hist[jj] = np.sum(weight_MR[sel])/bin
                    
          
            #print(np.log10(simps(hist, x=x_bin)) )
            cumulative_ND_blue[ii] = np.log10(simps(hist, x=x_bin)) 
           

        #Lgalaxies
        subplot.plot(ThisRedshiftList, cumulative_ND_red,color='red', linewidth=4, linestyle=':') 
        subplot.plot(ThisRedshiftList, cumulative_ND_blue,color='blue', linewidth=4, linestyle=':') 
      
        #EAGLE
        subplot.plot(eagle_z, np.log10(eagle_passive),color='red', linewidth=2, linestyle=':') 
        subplot.plot(eagle_z, np.log10(eagle_SF),color='blue', linewidth=2, linestyle=':') 
       
    
        #WRITE OUTPUT
        file = Datadir+"cumulative_ND_red_"+model_to_print+"M_"+ char_Mass +".txt"
        write_to_file(ThisRedshiftList, cumulative_ND_red, file)
        file = Datadir+"cumulative_ND_blue_"+model_to_print+"M_"+ char_Mass +".txt"
        write_to_file(ThisRedshiftList, cumulative_ND_blue, file)
            
    #LABELS    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.90, color='black', 
                xlog=0, ylog=0, label='$\log_{10}(M_*/$'+Msun+'$)>11.0$', fontsize=12, fontweight='normal') 
    
    plot_label (subplot,'label',xlim,ylim,x_percentage=0.12,y_percentage=0.83, 
                color='red',xlog=0,ylog=0,label='Passive', fontsize=10, fontweight='normal')             
    plot_label (subplot, 'line', xlim, ylim,x_percentage=0.05,y_percentage=0.845,
                color='red', x2_percentage=0.1,xlog=0, ylog=0, linestyle='-', linewidth=2)
            
    plot_label (subplot,'label',xlim,ylim,x_percentage=0.12,y_percentage=0.77, 
                color='blue',xlog=0,ylog=0,label='Star Forming', fontsize=10, fontweight='normal')             
    plot_label (subplot, 'line', xlim, ylim,x_percentage=0.05,y_percentage=0.785,
                color='blue', x2_percentage=0.1,xlog=0, ylog=0, linestyle='-', linewidth=2)
    
        
    plot_label (subplot,'label',xlim,ylim,x_percentage=0.61,y_percentage=0.89, 
                color='black',xlog=0,ylog=0,label='Illustris', fontsize=10, fontweight='normal')             
    plot_label (subplot, 'line', xlim, ylim,x_percentage=0.52,y_percentage=0.905,
                color='black', x2_percentage=0.59,xlog=0, ylog=0, linestyle='-', linewidth=2)
    
    plot_label (subplot,'label',xlim,ylim,x_percentage=0.61,y_percentage=0.83, 
                color='black',xlog=0,ylog=0,label='TNG300', fontsize=10, fontweight='normal')             
    plot_label (subplot, 'line', xlim, ylim,x_percentage=0.52,y_percentage=0.845,
                color='black', x2_percentage=0.59,xlog=0, ylog=0, linestyle='--', linewidth=2)
            
    plot_label (subplot,'label',xlim,ylim,x_percentage=0.61,y_percentage=0.77, 
                color='black',xlog=0,ylog=0,label='Eagle', fontsize=10, fontweight='normal')             
    plot_label (subplot, 'line', xlim, ylim,x_percentage=0.52,y_percentage=0.785,
                color='black', x2_percentage=0.59,xlog=0, ylog=0, linestyle=':', linewidth=2)
            
    plot_label (subplot,'label',xlim,ylim,x_percentage=0.61,y_percentage=0.71, 
                color='black',xlog=0,ylog=0,label='Magneticum', fontsize=10, fontweight='normal')             
    plot_label (subplot, 'line', xlim, ylim,x_percentage=0.52,y_percentage=0.725,
                color='black', x2_percentage=0.59,xlog=0, ylog=0, linestyle='-.', linewidth=2)
            
    plot_label (subplot,'label',xlim,ylim,x_percentage=0.61,y_percentage=0.65, 
                color='black',xlog=0,ylog=0,label='L-Galaxies (SAM)', fontsize=10, fontweight='normal')            
    plot_label (subplot, 'line', xlim, ylim,x_percentage=0.52,y_percentage=0.665,
                color='black', x2_percentage=0.59,xlog=0, ylog=0, linestyle=':', linewidth=4)    
                         
        
        
        
    '''#MODELS   
    fa = Datadir+"SMHM_allz_"+model_name+"_median_z"+char_redshift+".txt"
    (x_axis,y_axis)=read_file(fa)
    fa = Datadir+"SMHM_allz_"+model_name+"_pc16_z"+char_redshift+".txt"
    (x_axis,pc16)=read_file(fa)
    fa = Datadir+"SMHM_allz_"+model_name+"_pc84_z"+char_redshift+".txt"
    (x_axis,pc84)=read_file(fa)
        
    #subplot.plot(x_axis, y_axis, color=plot_color[ii], linewidth=2)   
    #subplot.plot(x_axis, pc16, color=plot_color[ii], linewidth=2)   
    #subplot.plot(x_axis, pc84, color=plot_color[ii], linewidth=2) 
    error_down = y_axis-pc16
    error_up = pc84-y_axis
    asy_yerror = [error_down,error_up]
    subplot.errorbar(x_axis+ii*0.05, y_axis, yerr=asy_yerror, color=plot_color[ii], linewidth=2)   
           
    x_arr=np.arange(xlim[0],xlim[1]+0.05,0.05)           
    y_arr=x_arr*0.+np.log10(0.155)      
    subplot.plot(x_arr,y_arr,color='black', linestyle='--', linewidth=2) ''' 
        
        
    #endfor
    
    
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()
    
    return 
    
#end cumulative_ND

def cooling_radius(ThisRedshiftList):

    xlim=[10., 14.5]
    ylim=[-2.,1.]    
        
        
    
    #PLOT current model
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1.0))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
    subplot.yaxis.set_minor_locator(MultipleLocator(0.5)) 
    
    xlab='$\log_{10}(M_{\mathrm{vir}}/$'+Msun+'$)$' 
    ylab='$R_{\mathrm{cool}}/R_{\mathrm{vir}}$'
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14) 

    plot_color=['red','blue']
    
    for ii in range (0,len(ThisRedshiftList)):
  
        #MODEL
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)
        
        #scater plot on subset
        G0_MR=G_MR[sel]   
        G0_MR=G0_MR[(G0_MR['StellarMass']>0.0) & (G0_MR['Type']==0)]        
        G0_MR=np.random.choice(G0_MR, size=2000)          
        log_HaloMass=np.log10(G0_MR['Mvir']*1.e10/Hubble_h)   
        log_fraction=np.log10(G0_MR['CoolingRadius']/G0_MR['Rvir'])      
        subplot.scatter(log_HaloMass,log_fraction,s=5, color='blue') 
      
        #median plot on full sample
        G0_MR=G_MR[sel]   
        G0_MR=G0_MR[(G0_MR['StellarMass']>0.0) & (G0_MR['Type']==0)]       
        log_HaloMass=np.log10(G0_MR['Mvir']*1.e10/Hubble_h)   
        log_fraction=np.log10(G0_MR['CoolingRadius']/G0_MR['Rvir'])
      
        bin=0.5        
        (x_binned, median, mean, pc16, pc84,rms)= median_and_percentiles(bin, xlim[0], xlim[1],log_HaloMass,log_fraction) 
        #subplot.plot(x_binned, median,color=plot_color[ii], linewidth=2,linestyle='-') 
        subplot.plot(x_binned, median,color='black', linewidth=2,linestyle='-') 
           
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.25, y_percentage=0.17, 
                    color='black', xlog=0, ylog=0, label='Model Median', fontsize=10, fontweight='normal') 
        plot_label (subplot, 'line', xlim, ylim, x_percentage=0.16, y_percentage=0.19, x2_percentage=0.22, 
                    color='black', xlog=0, ylog=0, linestyle='-', linewidth=2)   
        
        #Behroozi2013
        z=ThisRedshiftList[ii]   
        a=1/(1+z)
        neu=np.exp(-4*a*a)
        log_epsilon=-1.777+(-0.006*(a-1)-0.*z)*neu-0.119*(a-1)
        log_M1=11.514+(-1.793*(a-1)-0.251*z)*neu
        alpha=-1.412+0.731*(a-1.)*neu
        delta=3.508+(2.608*(a-1)-0.043*z)*neu
        gamma=0.316+(1.319*(a-1)+0.279*z)*neu 
        
        x=0.
        f_0=-np.log10(10**(alpha*x) + 1) + delta * ((np.log10(1 + np.exp(x)))**gamma)/(1 + np.exp(10**(-x)))             
        log_Mh=np.arange(8.0,16.0,0.01) 
        x=log_Mh-log_M1
        f_log_Mh=-np.log10(10**(alpha*x) + 1) + delta * ((np.log10(1 + np.exp(x)))**gamma)/(1 + np.exp(10**(-x)))      
        log_mstar=log_epsilon+log_M1+f_log_Mh-f_0     
        #subplot.plot(log_Mh, log_mstar-log_Mh+1.35,color=plot_color[ii], linewidth=2,linestyle='--')  
        subplot.plot(log_Mh, log_mstar-log_Mh+1.35,color='limegreen', linewidth=2,linestyle='-')  
        #subplot.plot(log_Mh, log_mstar-log_Mh+1.1,color='limegreen', linewidth=2,linestyle='-')  
      
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.25, y_percentage=0.10, 
                    color='black', xlog=0, ylog=0, label='Behroozi et al. (2013)', fontsize=10, fontweight='normal') 
        plot_label (subplot, 'line', xlim, ylim, x_percentage=0.16, y_percentage=0.12, x2_percentage=0.22, 
                    color='limegreen', xlog=0, ylog=0, linestyle='-', linewidth=2)
        
        
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()

    return 
#end cooling_radius    
    
    
def halo_mass_function(ThisRedshiftList):
        
    xlim=[10.5,15.0]
    ylim=[-8.5, 0.5]
    bin=0.25

    model_to_print='Hen15_other'
    
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()

    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    xlab='$\mathrm{log_{10}}(M_*/$'+Msun+'$)$'
    subplot.set_xlabel(xlab, fontsize=14)
    ylab='$\mathrm{log_{10}}(\phi / (\mathrm{Mpc^{-3}} \mathrm{dex}))$'
    subplot.set_ylabel(ylab, fontsize=14)
    
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
    subplot.yaxis.set_minor_locator(MultipleLocator(0.25))
  
    for i_z in range(0,len(ThisRedshiftList)):
        
        char_redshift="%0.2f" % ThisRedshiftList[i_z]
        
    
        #MR
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]   
        G0_MR=G0_MR[(G0_MR['Mvir']>0.) & (G0_MR['Type']==0)]
        Mass_MR=np.log10(G0_MR['Mvir']*1.e10/Hubble_h)
      
        (x_axis,y_axis) = plot_mass_function(subplot, Mass_MR, Volume_MR, xlim, bin, MRII==0, color='red')
        
        '''#MRII
        if(MRII==1):
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, i_z, FullSnapshotList_MRII)
        
            G0_MRII=G_MRII[sel]   
            G0_MRII=G0_MRII[G0_MRII['StellarMass']>0.]
            Mass_MRII=stellar_mass_with_err(G0_MRII, Hubble_h, ThisRedshiftList[i_z])
            
        if(MRII==0):
            (x_axis,y_axis) = plot_mass_function(subplot, Mass_MR, Volume_MR, xlim, bin, MRII, color='red')
        else:
            (x_axis,y_axis) = plot_mass_function(subplot, Mass_MR, Volume_MR, xlim, bin, MRII, 
                                                 Mass_MRII=Mass_MRII, Volume_MRII=Volume_MRII, color='red')'''
        
        #print(x_axis,y_axis)         
        
        #WRITE OUTPUT
        fa = open(Datadir+"HMF_"+model_to_print+"_z"+char_redshift+".txt", "w")
        fa.write("%d\n" % len(x_axis))
        for kk in range (0,len(x_axis)):               
            fa.write("%0.2f " % x_axis[kk] + "%0.2f\n" % y_axis[kk])         
        fa.close()    
    
    
    #endfor    
    subplot.fill_between(np.arange(11.5,12.6,0.1),np.arange(11.5,12.6,0.1)/np.arange(11.5,12.6,0.1)-1.+ylim[0],
                         np.arange(11.5,12.6,0.1)/np.arange(11.5,12.6,0.1)-1.+ylim[1],
                         facecolor='grey', interpolate=True, alpha=0.3)  
      
    
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')   
    plt.close()
        
    return 
    
#endif halo_mass_function==1:    


def halo_mass_density_above_mass(ThisRedshiftList):
    
    Time = np.array([0.771,0.942,1.186,1.558,2.171,3.316,5.903,6.516,7.289,8.143,
                     8.556,9.398,10.656,12.411, 13.721])
  
    #xlim=[0.0,9.0]
    xlim=[0.0,1.0]
    ylim=[-5.0, -2.0]
    bin=0.25

    model_to_print='Hen15'
    
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()

    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    xlab='$\log_{10}(1+z)$'
    subplot.set_xlabel(xlab, fontsize=14)
    ylab='$\mathrm{log_{10}}(\phi [\mathrm{Mpc^{-3}}])$'
    subplot.set_ylabel(ylab, fontsize=14)
    
    majorFormatter = FormatStrFormatter('%d')
    #subplot.xaxis.set_major_locator(MultipleLocator(1))   
    subplot.xaxis.set_major_locator(MultipleLocator(0.5))   
    #subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
    subplot.xaxis.set_minor_locator(MultipleLocator(0.1))      
    subplot.yaxis.set_minor_locator(MultipleLocator(0.25))
  
    ND_115_160 = np.zeros(len(ThisRedshiftList), dtype=np.float32)
    ND_115_125 = np.zeros(len(ThisRedshiftList), dtype = np.float32)
    
    bin=0.2
    hist_limits = [10.0,16.0]
    bin_arr=np.arange(hist_limits[0],hist_limits[1]+bin,bin)
    
    for i_z in range(0,len(ThisRedshiftList)):
        
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]   
        G0_MR=G0_MR[(G0_MR['Mvir']>0.) & (G0_MR['Type']==0)]
        Mass_MR=np.log10(G0_MR['Mvir']*1.e10/Hubble_h)
              
        hist_MR=np.histogram(Mass_MR, bins=bin_arr, range=(xlim[0],xlim[1]))      
        x_axis=hist_MR[1][0:len(hist_MR[1][:])-1]+bin/2.    
        hist_MR=hist_MR[0]       
        
        sel = (x_axis>=11.5)
        ND_115_160[i_z] = np.log10(sum(hist_MR[sel])/(Volume_MR))
        sel = (x_axis>=11.5) & (x_axis<=12.5)
        ND_115_125[i_z] = np.log10(sum(hist_MR[sel])/(Volume_MR))
        
        
    #WRITE OUTPUT
    fa = open(Datadir+"ND_115_125_"+model_to_print+".txt", "w")
    fa.write("%d\n" % len(ThisRedshiftList))
    for kk in range (0,len(ThisRedshiftList)):               
        fa.write("%0.2f " % ThisRedshiftList[kk] + "%0.2f\n" % ND_115_125[kk])         
    fa.close()
    
    fa = open(Datadir+"ND_115_160_"+model_to_print+".txt", "w")
    fa.write("%d\n" % len(ThisRedshiftList))
    for kk in range (0,len(ThisRedshiftList)):               
        fa.write("%0.2f " % ThisRedshiftList[kk] + "%0.2f\n" % ND_115_160[kk])         
    fa.close()
    
    
    #endfor
   
    subplot.plot(np.log10(1.+np.array(ThisRedshiftList)), ND_115_160, color='red',linestyle='--', linewidth=2)
    subplot.plot(np.log10(1.+np.array(ThisRedshiftList)), ND_115_125, color='red',linestyle='-', linewidth=2)
    
    
    
    for i_z in range(0,len(ThisRedshiftList)):
        
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]   
        G0_MR=G0_MR[(G0_MR['Mvir']>0.) & (G0_MR['Type']<2)]
        #InfallVmax
        Mass_MR=np.log10(G0_MR['Mvir']*1.e10/Hubble_h)
        
        Mvir = Vvir_to_Mvir(G0_MR['InfallVmax'], ThisRedshiftList[i_z], Omega=0.315 , OmegaLambda=0.685) 
        sel = (G0_MR['Type'] == 1)
        Mass_MR[sel]=np.log10(Mvir[sel]*1.e10/Hubble_h)
        
        hist_MR=np.histogram(Mass_MR, bins=bin_arr, range=(xlim[0],xlim[1]))      
        x_axis=hist_MR[1][0:len(hist_MR[1][:])-1]+bin/2.    
        hist_MR=hist_MR[0]      
      
        sel = (x_axis>=11.5) & (x_axis<=12.5)
        ND_115_125[i_z] = np.log10(sum(hist_MR[sel])/(Volume_MR))
    
    #endfor  
    subplot.plot(np.log10(1.+np.array(ThisRedshiftList)), ND_115_125, color='red',linestyle=':', linewidth=2)
        
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.1, y_percentage=0.45, color='black', 
                xlog=0, ylog=0, label='$11.5<\log_{10}(M_{\mathrm{vir}})<12.5$', fontsize=10, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.08, y_percentage=0.47, color='red', 
                x2_percentage=0.015, xlog=0, ylog=0, linestyle='-', linewidth=2)
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.1, y_percentage=0.37, color='black', 
                xlog=0, ylog=0, label='$\log_{10}(M_{\mathrm{vir}})>11.5$', fontsize=10, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.08, y_percentage=0.39, color='red', 
                x2_percentage=0.015, xlog=0, ylog=0, linestyle='--', linewidth=2)
                
    
    
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')   
    plt.close()
    
    
    
    
    #xlim=[0.0,9.0]
    xlim=[0.0,1.0]
    ylim=[7.0, 10.5]
    bin=0.25
    
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()

    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    xlab='$\log_{10}(1+z)$'
    subplot.set_xlabel(xlab, fontsize=14)
    ylab='$\mathrm{log_{10}}(\phi [\mathrm{Mpc^{-3}}])$'
    subplot.set_ylabel(ylab, fontsize=14)
    
    #subplot.xaxis.set_major_locator(MultipleLocator(1))   
    subplot.xaxis.set_major_locator(MultipleLocator(0.5))   
    #subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
    subplot.xaxis.set_minor_locator(MultipleLocator(0.1))      
    subplot.yaxis.set_minor_locator(MultipleLocator(0.25))
  
    ND_115_160 = np.zeros(len(ThisRedshiftList), dtype=np.float32)
    ND_115_125 = np.zeros(len(ThisRedshiftList), dtype = np.float32)
    ND_115_125_subhalos = np.zeros(len(ThisRedshiftList), dtype = np.float32)
    
    bin=0.2
    hist_limits = [10.0,16.0]
    bin_arr=np.arange(hist_limits[0],hist_limits[1]+bin,bin)
    
    for i_z in range(0,len(ThisRedshiftList)):
        
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]   
        G0_MR=G0_MR[(G0_MR['Mvir']>0.) & (G0_MR['Type']==0)]
        Mass_MR=np.log10(G0_MR['Mvir']*1.e10/Hubble_h)
        
        Mvir = Vvir_to_Mvir(G0_MR['InfallVmax'], ThisRedshiftList[i_z], Omega=0.315 , OmegaLambda=0.685) 
        sel = (G0_MR['Type'] == 1)
        Mass_MR[sel]=np.log10(Mvir[sel]*1.e10/Hubble_h)
        
        hist_MR=np.histogram(Mass_MR, bins=bin_arr, range=(xlim[0],xlim[1]))      
        x_axis=hist_MR[1][0:len(hist_MR[1][:])-1]+bin/2.    
        hist_MR=hist_MR[0]       
        
        sel = (x_axis>=11.5)
        ND_115_160[i_z] = np.log10(sum(hist_MR[sel]*(10**x_axis[sel]))/(Volume_MR))
        
        sel = (x_axis>=11.5) & (x_axis<=12.5)
        ND_115_125[i_z] = np.log10(sum(hist_MR[sel]*(10**x_axis[sel]))/(Volume_MR))
        
        
    #WRITE OUTPUT
    fa = open(Datadir+"ND_mass_density_115_125_"+model_to_print+".txt", "w")
    fa.write("%d\n" % len(ThisRedshiftList))
    for kk in range (0,len(ThisRedshiftList)):               
        fa.write("%0.2f " % ThisRedshiftList[kk] + "%0.2f\n" % ND_115_125[kk])         
    fa.close()
    
    fa = open(Datadir+"ND_mass_density_115_160_"+model_to_print+".txt", "w")
    fa.write("%d\n" % len(ThisRedshiftList))
    for kk in range (0,len(ThisRedshiftList)):               
        fa.write("%0.2f " % ThisRedshiftList[kk] + "%0.2f\n" % ND_115_160[kk])         
    fa.close()
    
    
    #endfor
    
    subplot.plot(np.log10(1.+np.array(ThisRedshiftList)), ND_115_160, color='red',linestyle='--', linewidth=2)
    subplot.plot(np.log10(1.+np.array(ThisRedshiftList)), ND_115_125, color='red',linestyle='-', linewidth=2)
    
    
    for i_z in range(0,len(ThisRedshiftList)):
        
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]   
        G0_MR=G0_MR[(G0_MR['Mvir']>0.) & (G0_MR['Type']<2)]
        Mass_MR=np.log10(G0_MR['Mvir']*1.e10/Hubble_h)
        
        Mvir = Vvir_to_Mvir(G0_MR['InfallVmax'], ThisRedshiftList[i_z], Omega=0.315 , OmegaLambda=0.685) 
        sel = (G0_MR['Type'] == 1)
        Mass_MR[sel]=np.log10(Mvir[sel]*1.e10/Hubble_h)
        
        hist_MR=np.histogram(Mass_MR, bins=bin_arr, range=(xlim[0],xlim[1]))      
        x_axis=hist_MR[1][0:len(hist_MR[1][:])-1]+bin/2.    
        hist_MR=hist_MR[0]       
     
        sel = (x_axis>=11.5) & (x_axis<=12.5)
        ND_115_125_subhalos[i_z] = np.log10(sum(hist_MR[sel]*(10**x_axis[sel]))/(Volume_MR))
    
    #endfor
 
    subplot.plot(np.log10(1.+np.array(ThisRedshiftList)), ND_115_125_subhalos, color='red',linestyle=':', linewidth=2)
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.1, y_percentage=0.45, color='black', 
                xlog=0, ylog=0, label='$11.5<\log_{10}(M_{\mathrm{vir}})<12.5$', fontsize=10, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.08, y_percentage=0.47, color='red', 
                x2_percentage=0.015, xlog=0, ylog=0, linestyle='-', linewidth=2)
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.1, y_percentage=0.37, color='black', 
                xlog=0, ylog=0, label='$\log_{10}(M_{\mathrm{vir}})>11.5$', fontsize=10, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.08, y_percentage=0.39, color='red', 
                x2_percentage=0.015, xlog=0, ylog=0, linestyle='--', linewidth=2)
                
    
    
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_mass_density_'+current_function+'.pdf')   
    plt.close()
    
    
    DM_Dt_115_160 = []
    #Dt = []
    for ii in range(0, len(ND_115_160)):
        if ii>0:
            DM_Dt_115_160.append((10**ND_115_160[ii]-10**ND_115_160[ii-1])/(Time[ii]-Time[ii-1])/(1e9))
            #Dt.append(Time[ii]-Time[ii-1])
    print(ND_115_160)
    print(DM_Dt_115_160)
    
    DM_Dt_115_125 = []   
    for ii in range(0, len(ND_115_125)):
        if ii>0:
            DM_Dt_115_125.append((10**ND_115_125[ii]-10**ND_115_125[ii-1])/(Time[ii]-Time[ii-1])/(1e9))
    
    DM_Dt_115_125_subhalos = []   
    for ii in range(0, len(ND_115_125)):
        if ii>0:
            DM_Dt_115_125_subhalos.append((10**ND_115_125_subhalos[ii]-10**ND_115_125_subhalos[ii-1])/(Time[ii]-Time[ii-1])/(1e9))
         
    
    #xlim=[0.0,9.0]
    xlim=[0.0,1.0]
    ylim=[-1.5, 1.0]
    bin=0.25
    
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()

    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    xlab='$\log_{10}(1+z)$'
    subplot.set_xlabel(xlab, fontsize=14)
    ylab='$\mathrm{log_{10}}(\Delta M/\Delta t / ($'+Msun+'$ yr^{-1} \mathrm{Mpc^{-3}}))$'
    subplot.set_ylabel(ylab, fontsize=14)
    
    #subplot.xaxis.set_major_locator(MultipleLocator(1))   
    subplot.xaxis.set_major_locator(MultipleLocator(0.5))   
    #subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
    subplot.xaxis.set_minor_locator(MultipleLocator(0.1))      
    subplot.yaxis.set_minor_locator(MultipleLocator(0.25))
    
    z_list = []
    for ii in range(0, len(ThisRedshiftList)):
        if ii>0:
            z_list.append((ThisRedshiftList[ii]+ThisRedshiftList[ii-1])/2.0) 
    
    
    subplot.plot(np.log10(1.+np.array(z_list)), smooth(np.log10(DM_Dt_115_160), 1), color='red',linestyle='--', linewidth=2)
    subplot.plot(np.log10(1.+np.array(z_list)), smooth(np.log10(DM_Dt_115_125), 1), color='red',linestyle='-', linewidth=2)
    subplot.plot(np.log10(1.+np.array(z_list)), smooth(np.log10(DM_Dt_115_125_subhalos), 1), color='red',linestyle=':', linewidth=2)
    
    #subplot.plot(np.log10(1.+np.array(z_list)), np.log10(DM_Dt_115_160), color='red',linestyle='--', linewidth=2)
    #subplot.plot(np.log10(1.+np.array(z_list)), np.log10(DM_Dt_115_125), color='red',linestyle='-', linewidth=2)
    #subplot.plot(np.log10(1.+np.array(z_list)), np.log10(DM_Dt_115_125_subhalos), color='red',linestyle=':', linewidth=2)


    #DM_Dt_115_125[-1]=0.08
    #subplot.plot(np.log10(1.+np.array(z_list)), np.log10(DM_Dt_115_125), color='red',linestyle=':', linewidth=2)
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.1, y_percentage=0.85, color='black', 
                xlog=0, ylog=0, label='$11.5<\log_{10}(M_{\mathrm{vir}})<12.5$', fontsize=10, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.08, y_percentage=0.87, color='red', 
                x2_percentage=0.015, xlog=0, ylog=0, linestyle='-', linewidth=2)
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.1, y_percentage=0.77, color='black', 
                xlog=0, ylog=0, label='$\log_{10}(M_{\mathrm{vir}})>11.5$', fontsize=10, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.08, y_percentage=0.79, color='red', 
                x2_percentage=0.015, xlog=0, ylog=0, linestyle='--', linewidth=2)
         
    
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_DM_Dt_'+current_function+'.pdf')   
    plt.close()
    
    return 
    
#endif halo_mass_density_above_mass    

def stellar_mass_function(ThisRedshiftList):
        
    xlim=[8.0,12.5]
    ylim=[-6.5, 0.5]
    bin=0.25
    
    fig = plt.figure(figsize=(two_two_size_small[0],two_two_size_small[1]))
    grid = gridspec.GridSpec(2, 2)
    grid.update(wspace=0.0, hspace=0.0)
    
    for i_z in range(0,len(ThisRedshiftList)):
        
        char_redshift="%0.2f" % ThisRedshiftList[i_z]
        
        subplot=plt.subplot(grid[i_z])

        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        if i_z==2 or i_z == 3: 
            xlab='$\mathrm{log_{10}}(M_*/$'+Msun+'$)$'
            subplot.set_xlabel(xlab, fontsize=16)
      
        if i_z==0 or i_z == 2:
            ylab='$\mathrm{log_{10}}(\phi /(\mathrm{Mpc^{-3}} \mathrm{dex}))$'
            subplot.set_ylabel(ylab, fontsize=16)
                       
    
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
        subplot.yaxis.set_minor_locator(MultipleLocator(0.25))
  
        if i_z==1 or i_z==3:
            subplot.tick_params(axis='y', which='both', left='on', labelleft='off')
        if i_z==0 or i_z==1:    
            subplot.tick_params(axis='x', which='both', bottom='on', labelbottom='off')
           
    
        #OBSERVATIONS             
        file = MCMCdir + '/ObsConstraints/StellarMassFunction_z'+char_redshift+'.txt'        
        f = open(file, 'r')     
        line = int(f.readline())     
        obs = Table.read(file, format='ascii', data_start=1, data_end=line+1)
        
        obs_xbin=obs['col1']+(obs['col2']-obs['col1'])/2.
        asy_yerror = [np.log10(obs['col3']/(obs['col3']-obs['col4'])), 
                      np.log10((obs['col3']+obs['col4'])/obs['col3'])]       
        subplot.errorbar(obs_xbin-2.*np.log10(Hubble_h), np.log10(obs['col3'])+3.*np.log10(Hubble_h),
                         yerr=asy_yerror,fmt='o', markersize=5, ecolor='blue', color='blue',zorder=+3,capsize=2)
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
      
    
        #MR
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]   
        G0_MR=G0_MR[G0_MR['StellarMass']>0.]
        Mass_MR=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[i_z])
      
        #MRII
        if(MRII==1):
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, i_z, FullSnapshotList_MRII)
        
            G0_MRII=G_MRII[sel]   
            G0_MRII=G0_MRII[G0_MRII['StellarMass']>0.]
            Mass_MRII = np.log10(G0_MRII['StellarMass']*1.e10/Hubble_h)
            #Mass_MRII=stellar_mass_with_err(G0_MRII, Hubble_h, ThisRedshiftList[i_z])
            
            #conn = connect("bhenriques", password="fPMitZK")
            #query = 'select stellarmass from import_test..SA_GALTREE_40 where snapnum='+str(G0_MRII['SnapNum'][0])
            #result = execute_query(conn, query)
            #Mass_MRII = np.log10(result['stellarmass']*1.e10/Hubble_h)
            
        #(x_axis,y_axis) = plot_mass_function(subplot, Mass_MR, Volume_MRII*100., xlim, bin, 0, color='red')
            
        if(MRII==0):
            (x_axis,y_axis) = plot_mass_function(subplot, Mass_MR, Volume_MR, xlim, bin, MRII, color='red')
        else:
            (x_axis,y_axis) = plot_mass_function(subplot, Mass_MR, Volume_MR, xlim, bin, MRII, 
                                                 Mass_MRII=Mass_MRII, Volume_MRII=Volume_MRII, color='red')
        
        #print(x_axis,y_axis)         
        
        #WRITE OUTPUT      
        if(write_to_file==1):
            df = pd.DataFrame({'log10_M':x_axis, 'log10_phi':y_axis})
            file = Datadir+file_to_write+'SMF'+str(f'_z{ThisRedshiftList[i_z]:0.2f}')+'.csv'
            df.to_csv(file,index=False)
            #df = pd.read_csv(file)
            #subplot.plot(df['log10_M'],df['log10_phi'], color='black')
        
    
        #MCMC sample
        if opt_plot_MCMC_sample==1:
            file = MCMCSampledir + 'mcmc_plus_obs_StellarMassFunction_z'+char_redshift+'.txt' 
            if os.path.isfile(file):
                obs = Table.read(file, format='ascii')      
                subplot.plot(obs['col1']-2.*np.log10(Hubble_h),np.log10(obs['col4'])+3.*np.log10(Hubble_h), 
                             color='black', linewidth=2) 
                #subplot.plot(obs['col1']-2.*np.log10(Hubble_h),np.log10(obs['col2'])+3.*np.log10(Hubble_h), 
                #             color='black', linewidth=2, linestyle='--')
                if i_z==len(ThisRedshiftList)-1:
                    plot_label (subplot, 'label', xlim, ylim, 
                        x_percentage=0.55, y_percentage=0.85, color='black', xlog=0, ylog=0, 
                        label='MCMC sample', fontsize=13, fontweight='normal') 
                    plot_label (subplot, 'line', xlim, ylim,
                        x_percentage=0.44, y_percentage=0.87, color='black', x2_percentage=0.53, 
                        xlog=0, ylog=0, linestyle='-', linewidth=2)
        
        #conn = connect("bhenriques", password="fPMitZK")
        #query = '''select .5*(.5+floor(log10(stellarmass*1.e10/0.673)/.5)) as mass,
        #          count(*) as count
        #     from Henriques2015a_new..MRIIscPlanck1
        #    where snapnum=62 and stellarmass>0.001
        # group by floor(log10(stellarmass*1.e10/0.673)/.5)
        # order by 1'''
        #result = execute_query(conn, query)
        #subplot.plot((result['mass']),np.log10(result['count']/(Volume_MRII*511*0.5)))
                
        #LABELS
        if i_z==0:
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.91, 
                        color='black', xlog=0, ylog=0, label='Observations used in MCMC', 
                        fontsize=13, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.935, 
                        color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15) 
            
        if i_z==len(ThisRedshiftList)-1:
            plot_label_two_models (subplot, xlim, ylim, position='bottom_left')
                                  
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.35, 
                    color='black', xlog=0, ylog=0, label='z='+char_redshift[:-1], 
                    fontsize=14, fontweight='normal')      
       
           
    #endfor

 
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYF19_'+current_function+'.pdf')
    plt.close()
        
    return 
    
#endif stellar_mass_function==1:


def stellar_mass_function_z0_overplot(ThisRedshiftList):
      
    plt.rcParams.update({'font.size': 14, 'xtick.labelsize': 14, 'ytick.labelsize': 14})    
        
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
    
    xlab='$\log_{10}(M_*/$'+Msun+'$)$'      
    ylab='$\log_{10}(\phi [\mathrm{Mpc^{-3}} \mathrm{log_{10}}(M^{-1})])$'         
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
            #G0_MRII=G0_MRII[(G0_MRII['StellarMass']>0.) & (G0_MRII['Type']<2)]
            Mvir=np.log10(G0_MRII['Mvir']*1.e10/Hubble_h*0.155)
            
            bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
            hist_MRII=np.histogram(Mvir, bins=bin_arr, range=(xlim[0],xlim[1]))   
           
        
        #join MR+MRII & plot        
        if(MRII==1):
            cut_MR_MRII=16.0
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
     
        xlab='$\log_{10}(<M_{\mathrm{200c, Hen15}}>/$'+Msun+'$)$'  
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
                color='black', xlog=0, ylog=0, label=r'$f_b \times M_{\mathrm{200c}}$', fontsize=13, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, y_percentage=0.185, 
                color='black', x2_percentage=0.1, xlog=0, ylog=0, linestyle=':', linewidth=2)
    
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HWL18_smf_z0.pdf')
    plt.close()

    return    
#end stellar_mass_function_z0_overplot

  
def stellar_mass_function_allz_overplot(ThisRedshiftList):
       
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
    
    xlab='$\mathrm{log_{10}}(M_*/$'+Msun+'$)$'
    ylab='$\mathrm{log_{10}}(\phi /(\mathrm{Mpc^{-3}} \mathrm{dex}))$'
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
    #pdf.savefig()
    #plt.close()
   
    
    
    
    
    #Active    
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
    subplot.yaxis.set_minor_locator(MultipleLocator(0.5)) 
    
    xlab='$\mathrm{log_{10}}(M_*/$'+Msun+'$)$'
    ylab='$\mathrm{log_{10}}(\phi/(\mathrm{Mpc^{-3}} \mathrm{dex}))$'
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
    #pdf.savefig()
    #plt.close()
    
    #PASSIVE    
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
    subplot.yaxis.set_minor_locator(MultipleLocator(0.5)) 
    
    xlab='$\mathrm{log_{10}}(M_*/$'+Msun+'$)$'
    ylab='$\mathrm{log_{10}}(\phi/(\mathrm{Mpc^{-3}} \mathrm{dex}))$'
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
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()
   
    return 
#endif stellar_mass_function_allz_overplot==1:
           
      

def stellar_mass_function_feedback_overplot(ThisRedshiftList):
        
    plt.rcParams.update({'font.size': 14, 'xtick.labelsize': 14, 'ytick.labelsize': 14})        
        
    xlim=[8.5,13.5]
    ylim=[-6., 1.5]
        
    plot_color=['red','purple']        
   
    fig = plt.figure(figsize=(one_four_size_large[0],one_four_size_large[1]))
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
            ylab='$\log_{10}(\phi [\mathrm{Mpc^{-3}} \mathrm{log_{10}}(M^{-1})])$'  
        else:
            ylab=''
        xlab='$\log_{10}(M_*/$'+Msun+'$)$'      
               
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14) 
        
        if jj>0:
            subplot.tick_params(axis='y', which='both', left='on', labelleft='off')
        
        #HALO MASS FUNCTION 
        bin=0.1
        for ii in range(0,len(ThisRedshiftList)):
            #MR
            (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)
        
            G0_MR=G_MR[sel]   
            G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (G0_MR['Type']<3)]   
            Mvir=np.log10(G0_MR['Mvir']*1.e10/Hubble_h*0.155)
            
            bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
            hist_MR=np.histogram(Mvir, bins=bin_arr, range=(xlim[0],xlim[1]))   
        
            #MRII
            if(MRII==1):
                (sel)=select_current_redshift(G_MRII, ThisRedshiftList, ii, FullSnapshotList_MRII)
        
                G0_MRII=G_MRII[sel]   
                G0_MRII=G0_MRII[(G0_MRII['StellarMass']>0.) & (G0_MRII['Type']<3)]
                Mvir=np.log10(G0_MRII['Mvir']*1.e10/Hubble_h*0.155)
        
                bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
                hist_MRII=np.histogram(Mvir, bins=bin_arr, range=(xlim[0],xlim[1]))   
           
        
            #join MR+MRII & plot        
            if(MRII==1):
                cut_MR_MRII=10.5
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
       
            if((jj==0)):
                #def fit_func(x, a, b, c):
                #    return ( ( a*(10**(x-c))**(b+1.)) * np.exp(-10**(x-c)) * np.log(10) )     
                def fit_func(x, a, b, c, d, e):
                    return ( ( a*(10**(x-e))**(b+1.) + c*(10**(x-e))**(d+1.) ) * np.exp(-10**(x-e)) * np.log(10) )  
                
                sel=(y_axis>ylim[0]) & (x_axis>xlim[0]+2.) #& (x_axis<xlim[1]-2.)
                p_init=[0.0019, -0.77, 0.00074, -1.6, 10.80]
                #p_init=[0.0019, -2.07, 10.80]
                par, err = curve_fit(fit_func, x_axis[sel],10**y_axis[sel],p_init)                
                Mass = np.arange(xlim[0],xlim[1], 0.1)                 
                #MF = ( (par[0]*(10**(Mass-par[2]))**(par[1]+1.)) * np.exp(-10**(Mass-par[2])) * np.log(10) ) 
                MF = ( (par[0]*(10**(Mass-par[4]))**(par[1]+1.) + par[2]*(10**(Mass-par[4]))**(par[3]+1.)) * 
                      np.exp(-10**(Mass-par[4])) * np.log(10) )      
                #subplot.plot(Mass, np.log10(MF), color='red', linewidth=2, linestyle=':')  
                
            else:
                if((jj==1)):
                    #def fit_func(x, a, b, c):
                    #    return ( ( a*(10**(x-c))**(b+1.)) * np.exp(-10**(x-c)) * np.log(10) )     
                    def fit_func(x, a, b, c, d, e):
                        return ( ( a*(10**(x-e))**(b+1.) + c*(10**(x-e))**(d+1.) ) * np.exp(-10**(x-e)) * np.log(10) )  
                
                    sel=(y_axis>ylim[0]) & (x_axis>xlim[0]) #& (x_axis<xlim[1]-2.)  
                    
                    ''' if(ii==0):
                        p_init=[0.00019, -0.77, 0.0008, -2.13, 10.80]
                    else:
                        if(ii==1):
                            p_init=[0.00019, -0.77, 0.001, -2.13, 11.00]
                        else:
                            if(ii==2):
                                p_init=[0.00019, -0.77, 0.003, -1.83, 11.00]
                            else:
                                p_init=[0.00019, -0.77, 0.003, -1.83, 11.10]
                                
                    #p_init=[0.0019, -2.07, 10.80]
                    par, err = curve_fit(fit_func, x_axis[sel],10**y_axis[sel],p_init)                
                    Mass = np.arange(xlim[0],xlim[1], 0.1)                 
                    #MF = ( (par[0]*(10**(Mass-par[2]))**(par[1]+1.)) * np.exp(-10**(Mass-par[2])) * np.log(10) ) 
                    MF = ( (par[0]*(10**(Mass-par[4]))**(par[1]+1.) + par[2]*(10**(Mass-par[4]))**(par[3]+1.)) * 
                          np.exp(-10**(Mass-par[4])) * np.log(10) )      
                    subplot.plot(Mass, np.log10(MF), color='red', linewidth=2, linestyle='-') ''' 
                         
                    if(ii==0):
                        p_init=[0.00019, -0.77, 0.003, -1.83, 11.2]
                        
                    else:
                        if(ii==1):
                            p_init=[0.00019, -0.77, 0.003, -1.83, 11.00]                            
                        else:
                            if(ii==2):
                                p_init=[0.00019, -0.77, 0.001, -2.13, 11.00]
                            else:
                                p_init=[0.00019, -0.77, 0.0008, -2.13, 10.80]
                    par = p_init        
                    #MF = ( (par[0]*(10**(Mass-par[4]))**(par[1]+1.) + par[2]*(10**(Mass-par[4]))**(par[3]+1.)) * 
                    #      np.exp(-10**(Mass-par[4])) * np.log(10) )      
                    #subplot.plot(Mass, np.log10(MF), color='red', linewidth=2, linestyle='-')  
                else:
                    def fit_func(x, a, b, c, d, e):
                        return ( ( a*(10**(x-e))**(b+1.) + c*(10**(x-e))**(d+1.) ) * np.exp(-10**(x-e)) * np.log(10) )     

                    sel=y_axis>ylim[0]
                    p_init=[0.0019, -0.77, 0.00074, -1.73, 10.80]
                    par, err = curve_fit(fit_func, x_axis[sel],10**y_axis[sel],p_init)               
                    Mass = np.arange(xlim[0],xlim[1], 0.1)        
                    MF = ( (par[0]*(10**(Mass-par[4]))**(par[1]+1.) + par[2]*(10**(Mass-par[4]))**(par[3]+1.)) * 
                          np.exp(-10**(Mass-par[4])) * np.log(10) )            
                    #subplot.plot(Mass, np.log10(MF), color='red', linewidth=2, linestyle=':')        
           
            #sel = x_axis == np.amin(x_axis[x_axis>par[4]])
            #subplot.scatter(x_axis[sel], y_axis[sel], marker='o', color='blue', s=50)
            
            '''if((jj==0) & (ii==0)):
                print(x_axis,(10**y_axis)*500.**3)
                aux_x=np.zeros(len(y_axis)-1)
                aux_y=np.zeros(len(y_axis)-1)
                integral=np.zeros(len(y_axis)-1)
                for iii in range(0,len(y_axis)-1):
                    aux_x[iii]=((10**x_axis[iii])+(10**x_axis[iii+1]))/2.
                    aux_y[iii]=((10**y_axis[iii]*500.**3)+(10**y_axis[iii+1]*500.**3))/2.
                    integral[iii]=((10**x_axis[iii+1])-(10**x_axis[iii]))*aux_y[iii]
                print(aux_x,np.log10(integral))
                #print(np.log10(np.sum(integral[np.log10(aux_x)<11.0],axis=0)))
                #print(np.log10(np.sum(integral[np.log10(aux_x)>11.0],axis=0)))
                print(np.sum(integral[(np.log10(aux_x)>8.0) & (np.log10(aux_x)<9.0)],axis=0)/
                                       np.sum(integral,axis=0))
                print(np.sum(integral[(np.log10(aux_x)>9.0) & (np.log10(aux_x)<10.0)],axis=0)/
                                       np.sum(integral,axis=0))
                print(np.sum(integral[(np.log10(aux_x)>10.0) & (np.log10(aux_x)<11.0)],axis=0)/
                                       np.sum(integral,axis=0))
                print(np.sum(integral[(np.log10(aux_x)>11.0) & (np.log10(aux_x)<12.0)],axis=0)/
                                       np.sum(integral,axis=0))
                print(np.sum(integral[(np.log10(aux_x)>12.0) & (np.log10(aux_x)<13.0)],axis=0)/
                                       np.sum(integral,axis=0))
                print(np.sum(integral[(np.log10(aux_x)>13.0) & (np.log10(aux_x)<14.0)],axis=0)/
                                       np.sum(integral,axis=0))
                print(np.sum(integral[(np.log10(aux_x)>14.0) & (np.log10(aux_x)<15.0)],axis=0)/
                                       np.sum(integral,axis=0))'''
                
            #Hen15    
            if(jj<3):                   
                fa = Datadir+model_name[3]+char_redshift+".txt"
                (x_axis,y_axis)=read_file(fa)    
                subplot.plot(x_axis, y_axis, color=colors[ii], linewidth=2, linestyle=':') 
           
                
            #LABELS       
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.9, 
                        color='black', xlog=0, ylog=0, label=label_name[jj], fontsize=12, fontweight='normal') 
           
            if jj==0:
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.15, y_percentage=0.2, 
                            color='blue', xlog=0, ylog=0, label=r'$f_b \times M_{\mathrm{200c}}$', 
                            fontsize=12, fontweight='normal') 
                plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, y_percentage=0.215, 
                            color='blue', x2_percentage=0.13, xlog=0, ylog=0, linestyle='--', linewidth=2)
                
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.15, y_percentage=0.135,
                            color='blue', xlog=0, ylog=0, label=r'$M_{*}$', 
                            fontsize=12, fontweight='normal') 
                plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, y_percentage=0.155, 
                            color='blue', x2_percentage=0.13, xlog=0, ylog=0, linestyle='-', linewidth=2)
                
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.15, y_percentage=0.07, 
                            color='blue', xlog=0, ylog=0, label='Hen15', 
                            fontsize=12, fontweight='normal') 
                plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, y_percentage=0.085, 
                            color='blue', x2_percentage=0.13, xlog=0, ylog=0, linestyle=':', linewidth=2)
    
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.41, 
                            color=colors[0], xlog=0, ylog=0, label='z=0',fontsize=12, fontweight='normal') 
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.36, 
                            color=colors[1], xlog=0, ylog=0, label='z=1',fontsize=12, fontweight='normal') 
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.31, 
                            color=colors[2], xlog=0, ylog=0, label='z=2',fontsize=12, fontweight='normal') 
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.26, 
                            color=colors[3], xlog=0, ylog=0, label='z=3',fontsize=12, fontweight='normal') 
        
        
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HWL18_smf_feedback_overplot.pdf')
    plt.close()
    
    return 

#end stellar_mass_function_feedback_overplot



def redfraction_color_cut_cuts(ThisRedshiftList):
  
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
            xlim=[-25.,-13.]
            ylim=[0.7, 2.5]     
            #ylim=[0.1, 3.5]     
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
            y_arr=offset_color_cut[ii]-slope_color_cut[ii]*np.tanh((x_arr+18.07)/1.09)
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
           
            sel=(color_UV < minimum_y_color_cut[ii]) | (color_UV < (color_VJ*slope_color_cut[ii])+offset_color_cut[ii])
           
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
        
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()
    
    return  
#endif redfraction_color_cut_cuts


def redfraction_color_cut(ThisRedshiftList):
       
    xlim=[8.0,11.5]
    ylim=[0., 1.0]   
    bin=[0.1,0.05]
    Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]

    fig = plt.figure(figsize=(one_five_size_large[0],one_five_size_large[1]))
    grid = gridspec.GridSpec(1, 5)
    grid.update(wspace=0.0, hspace=0.0)

    for ii in range(0,len(ThisRedshiftList)):
        
        char_redshift="%0.2f" % ThisRedshiftList[ii]
        
        subplot = plt.subplot(grid[ii])
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        
        xlab='$\mathrm{log_{10}}(M_*/$'+Msun+'$)$'      
        if ii==0:
            ylab='Red Fraction'
        else:
            ylab=''      
        subplot.set_xlabel(xlab, fontsize=15), subplot.set_ylabel(ylab, fontsize=15)
        
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
                 fmt='o', markersize=5, ecolor='blue', color='blue',zorder=+3,capsize=2)
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
                sel_red=G0_MR[(color_ur>(offset_color_cut[ii]-slope_color_cut[ii]*np.tanh((Magr+18.07)/1.09))) &
                    (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)]
                sel_all=G0_MR[(StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)]     
                if len(sel_all)>0.:
                    RedFraction[ll]=float(len(sel_red))/float(len(sel_all))  
                else:
                    RedFraction[ll]=0.
                    
                if((MRII==1) & (Mass_arr[ll]<9.5)):  
                    sel_red=G0_MRII[(color_ur_MRII>(offset_color_cut[ii]-slope_color_cut[ii] *
                                                    np.tanh((Magr_MRII+18.07)/1.09))) &
                                    (StellarMass_MRII>Mass_arr[ll]-bin/2.) & (StellarMass_MRII<Mass_arr[ll]+bin/2.)]
                    sel_all=G0_MRII[(StellarMass_MRII>Mass_arr[ll]-bin/2.) & (StellarMass_MRII<Mass_arr[ll]+bin/2.)]     
                    if len(sel_all)>0.:
                        RedFraction[ll]=float(len(sel_red))/float(len(sel_all))  
                    else:
                        RedFraction[ll]=0.
                     
                    
        
        #z>0.
        else:                      
            Mass = np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)
            SSFR = np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h))   
            color_UV=G0_MR['MagDust'][:,0]-G0_MR['MagDust'][:,2]  
            color_VJ=G0_MR['MagDust'][:,2]-G0_MR['MagDust'][:,7]       
                  
            for ll in range(0,len(Mass_arr)):               
                #sel_red=G0_MR[(color_UV > minimum_y_color_cut[ii]) & 
                #              (color_UV > (color_VJ*slope_color_cut[ii])+offset_color_cut[ii]) &            
                #              (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.) & (G0_MR['Type']==0)]  
                sel_red=G0_MR[(SSFR<np.log10(2.*(1+ThisRedshiftList[ii])**2./(1.37e10))-1.0) &
                              (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.) ] 
                #sel_red=G0_MR[(SSFR<np.log10((1+ThisRedshiftList[ii])/(2.*1.37e10))) &
                #              (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.) ] 
                #sel_red=G0_MR[(SSFR<-11.0) &
                #              (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.) ] 
                
                sel_all=G0_MR[(StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.) ]
                if len(sel_all)>0.:
                    RedFraction[ll]=float(len(sel_red))/float(len(sel_all))
                else:
                    RedFraction[ll]=0.
                    
        subplot.plot(Mass_arr, RedFraction, color='red', linestyle='-', linewidth=2) 
        
        
        #WRITE OUTPUT      
        if(write_to_file==1):
            df = pd.DataFrame({'log10_M':Mass_arr, 'RedFraction':RedFraction})
            file = Datadir+file_to_write+'RedFraction'+str(f'_z{ThisRedshiftList[ii]:0.2f}')+'.csv'
            df.to_csv(file,index=False)
            #df = pd.read_csv(file)
            #subplot.plot(df['log10_M'],df['RedFraction'], color='black')
    
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
    
    
        #Previous Model        
        if do_previous_model2==1: 
            file = file_previous_model2+'_redfrac_colorcut_z'+char_redshift+'.txt' 
            model = Table.read(file, format='ascii')
            subplot.plot(model['col1']-2.*np.log10(Hubble_h), model['col2'],
                         color='red',linestyle=linestyle_previous_model2, linewidth=2)           
        # prefix_previous_model2
        #LABELS    
        if ii==4:
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.095, y_percentage=0.91, 
                        color='black', xlog=0, ylog=0, label='Observations used in MCMC', 
                        fontsize=11, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.06, y_percentage=0.925, 
                        color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.025) 
        
            plot_label (subplot,'label',xlim,ylim,x_percentage=0.16,y_percentage=0.82, 
                        color='black',xlog=0,ylog=0,label=prefix_this_model, fontsize=14, fontweight='normal')             
            plot_label (subplot, 'line', xlim, ylim,x_percentage=0.04,y_percentage=0.84,
                        color='red', x2_percentage=0.14,xlog=0, ylog=0, linestyle='-', linewidth=2)
            
            plot_label (subplot,'label',xlim,ylim,x_percentage=0.16,y_percentage=0.74, 
                        color='black',xlog=0,ylog=0,label=prefix_previous_model2, fontsize=14, fontweight='normal')     
            plot_label (subplot, 'line', xlim, ylim,x_percentage=0.04,y_percentage=0.76,
                        color='red', x2_percentage=0.14,xlog=0, ylog=0, linestyle='--', linewidth=2)
        
        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.5, 
                    color='black', xlog=0, ylog=0, label='z='+char_redshift[:-1], 
                    fontsize=15, fontweight='normal') 
            
            
            
        #CHANGE THIS*****************    
        #if ii==len(ThisRedshiftList)-1:
        #    plot_label_three_models (subplot, xlim, ylim, position='top_left')
                        
    #endfor
        
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYF19_redfraction_color_cut.pdf')  
    plt.close() 
    
    return 
#endif redfraction_color_cut


def read_model_0(file):
    with open(Datadir + file + '.dat', 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=' ')
        xx =[]
        yy = []
        yy_2 = []
        for row in reader:
            xx.append(float(row[1]))
            yy.append(float(row[2]))
            yy_2.append(float(row[3]))
            
    return np.array(xx), np.array(yy), np.array(yy_2)

def redfraction_SFR_cut(ThisRedshiftList):
    
    model_to_print = 'Hen15'
    
    xlim=[9.0,11.5]
    ylim=[0., 1.0]   
    bin=[0.1,0.05]
    Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]

    fig = plt.figure(figsize=(one_three_size_large[0],one_three_size_large[1]))
    grid = gridspec.GridSpec(1, 3)
    grid.update(wspace=0.0, hspace=0.0)

    for ii in range(0,len(ThisRedshiftList)):
        
        char_redshift="%0.2f" % ThisRedshiftList[ii]
        
        subplot = plt.subplot(grid[ii])
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        
        xlab='$\mathrm{log_{10}}(M_*/$'+Msun+'$)$'      
        if ii==0:
            ylab='Passive Fraction'
        else:
            ylab=''      
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)
              
        subplot.text(xlim[0]+0.98,ylim[1]+.025,'z='+char_redshift[:-1], fontsize=14, fontweight='normal')
                        
    
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
        subplot.yaxis.set_minor_locator(MultipleLocator(0.1))
    
        if ii>0:
            plt.tick_params(axis='y', which='both', left='on', labelleft='off')
                                 
        #OBSERVATIONS 
        if ii==0:
            #wetzel data for centrals at z=0
            xx = np.array([9.900, 10.300, 10.700, 11.100])
            yy = np.array([0.210, 0.418, 0.706, 0.933])
            error_down = np.array([0.005, 0.004, 0.005, 0.006])
            error_up = np.array([0.006, 0.004, 0.005, 0.005])
            asy_err = [error_down, error_up]
            subplot.errorbar(xx, yy, asy_err, fmt='o', markersize=5, ecolor='black', color='black')
        else:
            #MCMC for higher z
            file = MCMCdir + '/ObsConstraints/RedFraction_z'+char_redshift+'.txt'        
            f = open(file, 'r')     
            line = int(f.readline())     
            obs = Table.read(file, format='ascii', data_start=1, data_end=line+1)        
            obs_xbin=obs['col1']+(obs['col2']-obs['col1'])/2.             
            plot_color = 'grey'
            subplot.errorbar(obs_xbin-2*np.log10(Hubble_h), obs['col3'],obs['col4'], fmt='o', 
                             markersize=5, ecolor=plot_color, color=plot_color)
       
       
        #model 0 - Illustris
        if ii == 0: 
            file = 'centralsILL_135'
            (xx, yy, yy_2) = read_model_0(file)
            subplot.plot(xx, yy/yy_2, color=Illustris_color, linestyle='-', linewidth=2) 
        if ii == 1:
            file = 'centralsILL_85'
            (xx, yy, yy_2) = read_model_0(file)
            subplot.plot(xx, yy/yy_2, color=Illustris_color, linestyle='-', linewidth=2) 
        if ii == 2:        
            file = 'centralsILL_68'
            (xx, yy, yy_2) = read_model_0(file)
            subplot.plot(xx, yy/yy_2, color=Illustris_color, linestyle='-', linewidth=2) 
        
        #model 1 - TNG300
        if ii == 0: 
            file = 'centrals300_99'
            (xx, yy, yy_2) = read_model_0(file)
            subplot.plot(xx, yy/yy_2, color=TNG300_color, linestyle='-', linewidth=2) 
        if ii == 1:
            file = 'centrals300_50'
            (xx, yy, yy_2) = read_model_0(file)
            subplot.plot(xx, yy/yy_2, color=TNG300_color, linestyle='-', linewidth=2) 
        if ii == 2:        
            file = 'centrals300_33'
            (xx, yy, yy_2) = read_model_0(file)
            subplot.plot(xx, yy/yy_2, color=TNG300_color, linestyle='-', linewidth=2) 
       
            
        #MODEL
        bin=0.25        
        Mass_arr=np.arange(xlim[0],xlim[1],bin)
        RedFraction=np.zeros(len(Mass_arr),dtype=np.float32)
        #if ii == 0:
        #    SFR_cut = -11.0
        #else:
        #    SFR_cut = np.log10((1+ThisRedshiftList[ii])/(2.*1.37e10))
        SFR_cut = np.log10((1+ThisRedshiftList[ii])/(2.*1.37e10))
                
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]          
        log_StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
        log_SFR=np.log10(G0_MR['Sfr'])
        
        if(MRII==1):
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, ii, FullSnapshotList_MRII)        
            G0_MRII=G_MRII[sel]   
            log_StellarMass_MRII=stellar_mass_with_err(G0_MRII, Hubble_h, ThisRedshiftList[ii])
            log_SFR_MRII=np.log10(G0_MRII['Sfr'])
           
       
        for ll in range(0,len(Mass_arr)):                
            '''sel_red=G0_MR[((log_SFR-log_StellarMass)<SFR_cut) & (log_StellarMass>Mass_arr[ll]-bin/2.) & 
                          (log_StellarMass<Mass_arr[ll]+bin/2.)]
            sel_all=G0_MR[(log_StellarMass>Mass_arr[ll]-bin/2.) & (log_StellarMass<Mass_arr[ll]+bin/2.)]  
             
            if len(sel_all)>0.:
                RedFraction[ll]=float(len(sel_red))/float(len(sel_all))  
            else:
                RedFraction[ll]=0.
                    
            if((MRII==1) & (Mass_arr[ll]<9.5)):  
                sel_red=G0_MRII[((log_SFR_MRII-log_StellarMass_MRII)<SFR_cut) &
                                (log_StellarMass_MRII>Mass_arr[ll]-bin/2.) & (log_StellarMass_MRII<Mass_arr[ll]+bin/2.)]
                sel_all=G0_MRII[(log_StellarMass_MRII>Mass_arr[ll]-bin/2.) & (log_StellarMass_MRII<Mass_arr[ll]+bin/2.)]  
                if len(sel_all)>0.:
                    RedFraction[ll]=float(len(sel_red))/float(len(sel_all))  
                else:
                    RedFraction[ll]=0.'''
            sel_red=G0_MRII[((log_SFR_MRII-log_StellarMass_MRII)<SFR_cut) &
                                (log_StellarMass_MRII>Mass_arr[ll]-bin/2.) & (log_StellarMass_MRII<Mass_arr[ll]+bin/2.)]
            sel_all=G0_MRII[(log_StellarMass_MRII>Mass_arr[ll]-bin/2.) & (log_StellarMass_MRII<Mass_arr[ll]+bin/2.)]  
            if len(sel_all)>0.:
                RedFraction[ll]=float(len(sel_red))/float(len(sel_all))  
            else:
                RedFraction[ll]=0.
                        
      
        subplot.plot(Mass_arr, RedFraction, color=LGalaxies_color, linestyle='-', linewidth=2) 
      
        #WRITE OUTPUT
        file = Datadir+"RedFraction_SFR_cut_"+model_to_print+"_z"+char_redshift+".txt"
        write_to_file(Mass_arr, RedFraction, file)
            
        #LABELS    
        if ii==0:
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.91, 
                        color='black', xlog=0, ylog=0, label='Wetzel et al. 2012', 
                        fontsize=12, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.935, 
                        color='black', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.03) 
        
            plot_label (subplot,'label',xlim,ylim,x_percentage=0.12,y_percentage=0.83, 
                        color=Illustris_color,xlog=0,ylog=0,label='Illustris', fontsize=10, fontweight='normal')             
            plot_label (subplot, 'line', xlim, ylim,x_percentage=0.02,y_percentage=0.845,
                        color=Illustris_color, x2_percentage=0.1,xlog=0, ylog=0, linestyle='-', linewidth=2)
            
            plot_label (subplot,'label',xlim,ylim,x_percentage=0.12,y_percentage=0.77, 
                        color=TNG300_color,xlog=0,ylog=0,label='TNG300', fontsize=10, fontweight='normal')             
            plot_label (subplot, 'line', xlim, ylim,x_percentage=0.02,y_percentage=0.785,
                        color=TNG300_color, x2_percentage=0.1,xlog=0, ylog=0, linestyle='-', linewidth=2)
            
            plot_label (subplot,'label',xlim,ylim,x_percentage=0.12,y_percentage=0.71, 
                        color=EAGLE_color,xlog=0,ylog=0,label='Eagle', fontsize=10, fontweight='normal')             
            plot_label (subplot, 'line', xlim, ylim,x_percentage=0.02,y_percentage=0.725,
                        color=EAGLE_color, x2_percentage=0.1,xlog=0, ylog=0, linestyle='-', linewidth=2)
            
            plot_label (subplot,'label',xlim,ylim,x_percentage=0.12,y_percentage=0.65, color=Magneticum_Box2_color,
                        xlog=0,ylog=0,label='Magneticum', fontsize=10, fontweight='normal')             
            plot_label (subplot, 'line', xlim, ylim,x_percentage=0.02,y_percentage=0.665, color=Magneticum_Box2_color, 
                        x2_percentage=0.1,xlog=0, ylog=0, linestyle='-', linewidth=2)
            
            plot_label (subplot,'label',xlim,ylim,x_percentage=0.12,y_percentage=0.59, color=LGalaxies_color, 
                        xlog=0,ylog=0,label='L-Galaxies (Hen15)', fontsize=10, fontweight='normal')        
            plot_label (subplot, 'line', xlim, ylim,x_percentage=0.02,y_percentage=0.605, color=LGalaxies_color, 
                        x2_percentage=0.1,xlog=0, ylog=0, linestyle='-', linewidth=2)
        
  
        
        #ay2 = subplot.twiny() 
        #ylab = 'z='+char_redshift[:-1]
        #ay2.set_ylabel(ylab, fontsize=14)
        #ax2.set_ybound(subplot.get_ybound())        
        #ax2.set_yticks(y_binned)
        #ax2.set_yticklabels(y_axis_2)
        
            
            
        #CHANGE THIS*****************    
        #if ii==len(ThisRedshiftList)-1:
        #    plot_label_three_models (subplot, xlim, ylim, position='top_left')
                        
    #endfor
        
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')   
    plt.close() 
    
    return 
#endif redfraction_SFR_cut


def wetzel_passive_fraction_vs_stellar_mass(ThisRedshiftList):

    model_to_print = 'Hen15'
    
    xlim=[9.5, 11.5]
    ylim=[0.0, 1.0]   
       
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.25)) 
    
    xlab='$\mathrm{log_{10}}(M_*/$'+Msun+'$)$'       
    ylab='Passive Fraction'         
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14) 
        
    halo_mass = np.array([13.25,13.75,14.25])    
    plot_color = ['green', 'yellow', 'orange']    
    for ii in range (0,len(ThisRedshiftList)):    
        
        #Model
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]             
       
        log_StellarMass=np.log10(G0_MR['StellarMass']*1e10/Hubble_h)
        log_SFR=np.log10(G0_MR['Sfr'])
        log_HaloMass = np.log10(G0_MR['CentralMvir']*1e10/Hubble_h)
        
        for jj in range(0, len(halo_mass)):
            char_redshift="%0.2f" % ThisRedshiftList[ii]
            char_halo_mass="%0.2f" % halo_mass[jj]
            char_min_mass="%0.1f" % (halo_mass[jj]-0.25)
            char_max_mass="%0.1f" % (halo_mass[jj]+0.25)
            
            #OBSERVATIONS
            file=Datadir+'wetzel1_red_fraction_vs_stellar_mass_' + char_min_mass + '_' + char_max_mass + '.txt'            
            with open(file, 'r') as csvfile:
                reader = csv.reader(csvfile, delimiter=' ')
                kk=0
                for row in reader:
                    if kk==0:
                        xx = row
                    if kk==1:    
                        yy = row
                    if kk==2:
                        y_err_down = row
                    if kk==3:    
                        y_err_up = row
                        
                    kk+=1
            
            xx = [float(x) for x in xx]
            yy = [float(x) for x in yy]
            y_err_down = [float(x) for x in y_err_down]
            y_err_up = [float(x) for x in y_err_up]
            
            #print(np.array(xx), np.array(yy), np.array(y_err_down), np.array(y_err_up))        
            subplot.errorbar(np.array(xx), np.array(yy), [np.array(y_err_down), np.array(y_err_up)],  
                            fmt='o', markersize=5, ecolor=plot_color[jj], color=plot_color[jj])
               
            max_obs_mass=max(np.array(xx)[np.array(yy) > 0.])
            min_obs_mass=min(np.array(xx)[np.array(yy) > 0.])
            
            #MODEL         
            bin=0.1
            Mass_arr=np.arange(min_obs_mass-0.5, max_obs_mass+bin,bin)
            RedFraction=np.zeros(len(Mass_arr),dtype=np.float32)
            
            for ll in range(0,len(Mass_arr)):                  
                sel_red = (((log_SFR - log_StellarMass)<-11.0) & (G0_MR['Type']>0) & 
                           (log_HaloMass>halo_mass[jj]-0.25) & (log_HaloMass<halo_mass[jj]+0.25) & 
                           (log_StellarMass>Mass_arr[ll]-bin/2.) & (log_StellarMass<Mass_arr[ll]+bin/2.))
                
                sel_all = ((log_HaloMass>halo_mass[jj]-0.25) & (log_HaloMass<halo_mass[jj]+0.25) & 
                           (G0_MR['Type']>0) & (log_StellarMass>Mass_arr[ll]-bin/2.) & (log_StellarMass<Mass_arr[ll]+bin/2.))
            
                if len(G0_MR[sel_all])>0.:
                    RedFraction[ll]=float(len(G0_MR[sel_red]))/float(len(G0_MR[sel_all]))  
            
            sel = (Mass_arr>=min_obs_mass-bin/2.) & (Mass_arr<=max_obs_mass+bin/2.)            
            subplot.plot(Mass_arr[sel], RedFraction[sel],color=plot_color[jj], linewidth=2)
      
    
            #WRITE OUTPUT
            file = Datadir+"Wetzel_passive_fration"+model_to_print+"_M"+ char_halo_mass +"_z"+char_redshift+".txt"
            write_to_file(Mass_arr, RedFraction, file)
     
   
    #endfor
        
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')  
    plt.close()

    return   
  
#end wetzel

def metals_vs_stellarmass(ThisRedshiftList):
    
    xlim=[9.0,12.0]
    ylim=[-1.3, 0.7]   
    bin=0.2
        
    plot_color=['red','purple']        
   
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.25))     
    subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
    
    xlab='$\log_{10}(M_*/$'+Msun+'$)$'       
    ylab='$\log_{10}(Z_*/\mathrm{Z_{\odot}})$'         
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14) 
        
    #for ii in range (0,len(ThisRedshiftList)):
    for ii in range (0,1):        
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
            subplot.plot(xmass,obsp50+np.log10(0.02)-np.log10(0.0142),color='blue',linestyle='-',linewidth=2)
            subplot.plot(xmass,obsp16+np.log10(0.02)-np.log10(0.0142),color='blue',linestyle='--',linewidth=2)
            subplot.plot(xmass,obsp84+np.log10(0.02)-np.log10(0.0142),color='blue',linestyle='--',linewidth=2) 
            
            
            file=Datadir+'Obs_MZsR_Zahid+17_Zsun00142.fits'    
            fits_table=fits.open(file)   
            cols = fits_table[1].columns
            #cols.info()            
            Zahid = fits_table[1]           
            #['LOGMSTAR', 'LOGZSTAR_LW', 'LOGZSTAR_MW']
            subplot.scatter(Zahid.data['LOGMSTAR'],Zahid.data['LOGZSTAR_LW'],color='green',s=20, marker='o', alpha=0.5)
            subplot.scatter(Zahid.data['LOGMSTAR'],Zahid.data['LOGZSTAR_MW'],color='purple',s=20, marker='o', alpha=0.5)
            
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
            label=prefix_this_model
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.11, y_percentage=0.92, color='black', 
                        xlog=0, ylog=0, label=label, fontsize=12, fontweight='normal')             
            plot_label (subplot,'line',xlim,ylim,x_percentage=0.04, x2_percentage=0.09, y_percentage=0.94,
                        color=plot_color[0],xlog=0,ylog=0,linestyle='-',linewidth=2)
            
                 
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.11, y_percentage=0.85, 
                        color='black', xlog=0, ylog=0, label=prefix_previous_model2, fontsize=12, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, x2_percentage=0.09, y_percentage=0.87, 
                        color='red', xlog=0, ylog=0, linestyle=':', linewidth=2)
        
            #z=3
            #label=prefix_this_model+', z=3'
            #plot_label (subplot, 'label', xlim, ylim, 
            #            x_percentage=0.12, y_percentage=0.82, color='black', xlog=0, ylog=0, 
            #            label=label, fontsize=13, fontweight='normal')             
            #plot_label (subplot,'line',xlim,ylim,x_percentage=0.02,y_percentage=0.84,
            #            color=plot_color[1],x2_percentage=0.10,xlog=0,ylog=0,linestyle='-',linewidth=2)
            
            #LABELS
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.4, y_percentage=0.19, color='black', 
                        xlog=0, ylog=0, label='SDSS-DR2 (LW, Gallazzi2005)', fontsize=10, fontweight='normal')     
            plot_label (subplot,'line',xlim,ylim,x_percentage=0.33,y_percentage=0.205,
                        color='blue',x2_percentage=0.38,xlog=0,ylog=0,linestyle='-',linewidth=2)
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.4, y_percentage=0.12, color='black', 
                        xlog=0, ylog=0,label='SDSS-DR7 (LW, Zahid2017)', fontsize=10, fontweight='normal')       
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.37, y_percentage=0.135, 
                        color='green', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.001, alpha=0.5) 
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.4, y_percentage=0.05, color='black', 
                        xlog=0, ylog=0,label='SDSS-DR7 (MW, Zahid2017)', fontsize=10, fontweight='normal')       
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.37, y_percentage=0.065, 
                        color='purple', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.001, alpha=0.5 )
            #plot_label_three_models (subplot, xlim, ylim, position='bottom_left')
                  
        #Model
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0=G_MR[sel]             
        '''The selection on type<2 is needed because these galaxies have their cooling artifically 
        shut down and therefore erroneously enriched by the SNIA and AGB channels without 
        cooling to balance it out'''
        #G0=G0[(G0['StellarMass']>0.) & (G0['Type']<2) & (np.log10(G0['Sfr'])>-2.0) &
        #      (np.log10(G0['Sfr']/(G0['StellarMass']*1.e10/Hubble_h))>-11.0)]  
        
        SSFR = np.log10(G0['Sfr']/(G0['StellarMass']*1.e10/Hubble_h))       
        G0=G0[(G0['StellarMass']>0.) & (G0['Type']<2) & (SSFR>np.log10(2.*(1+ThisRedshiftList[ii])**2./(1.37e10))-1.0)]  
        
      
        
        if(opt_detailed_enrichment==0): 
            G0=G0[G0['MetalsStellarMass']>0.]
        else:
            G0=G0[(G0['MetalsStellarMass'][:,0]+G0['MetalsStellarMass'][:,1]+G0['MetalsStellarMass'][:,2]>0.)]
                 
        StellarMass=stellar_mass_with_err(G0, Hubble_h, ThisRedshiftList[ii])        
        (Metallicity, StellarMass)=get_metallicity(G0, 'StellarMass', StellarMass)       
        (x_binned_MR, median_MR, mean_MR, pc16_MR, pc84_MR,rms_MR)=median_and_percentiles (bin, xlim[0]-1.0, xlim[1]+1.0, 
                                                                               StellarMass, Metallicity)    
       
        if(MRII==1):
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, ii, FullSnapshotList_MRII)        
            G0=G_MRII[sel]    
            
            G0=G0[(G0['StellarMass']>0.) & (G0['Type']<2) & (np.log10(G0['Sfr'])>-2.0) &
              (np.log10(G0['Sfr']/(G0['StellarMass']*1.e10/Hubble_h))>-11.0)]        
            if(opt_detailed_enrichment==0): 
                G0=G0[G0['MetalsStellarMass']>0.]
            else:
                G0=G0[(G0['MetalsStellarMass'][:,0]+G0['MetalsStellarMass'][:,1]+G0['MetalsStellarMass'][:,2]>0.)]
       
            StellarMass=stellar_mass_with_err(G0, Hubble_h, ThisRedshiftList[ii])  
            (Metallicity, StellarMass)=get_metallicity(G0, 'StellarMass', StellarMass)           
            (x_binned_MRII,median_MRII,mean_MRII,pc16_MRII,pc84_MRII,rms_MRII)=median_and_percentiles(bin,xlim[0]-1.0, 
                                                                                                      xlim[1]+1.0,
                                                                                       StellarMass,Metallicity) 
                              
            cut_MR_MRII=8.5    
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
        
        sel=((pc16>ylim[0]) & (median!=0.))
        subplot.plot(x_binned[sel], median[sel],color=plot_color[ii], linewidth=2)
        if (ii==0):
            subplot.plot(x_binned[sel], pc16[sel],color=plot_color[ii], linewidth=2, linestyle='--')
            subplot.plot(x_binned[sel], pc84[sel],color=plot_color[ii], linewidth=2, linestyle='--')
    
    
        #get metals from database
        #conn = connect("bhenriques", password="fPMitZK")
        #query = '''select top 100 stellarmass as stellarmass, metalsstellarmass as metals 
        #             from Henriques2015a_new..MRIIscPlanck1 
        #            where snapnum=62 and stellarmass>1.0'''
        #result = execute_query(conn, query)        
        #Mass_MRII = np.log10(result['stellarmass']*1.e10/Hubble_h)
        #metals = result['metals']/result['stellarmass']/0.0142
        #subplot.scatter(Mass_MRII,np.log10(metals))
       
        
        #WRITE OUTPUT      
        if(write_to_file==1):
            df = pd.DataFrame({'log10_M':x_binned[sel], 'Metallicity':median[sel], 'pc16':pc16[sel], 'pc84':pc84[sel]})
            file = Datadir+file_to_write+'StellarMetals_vs_StellarMass'+str(f'_z{ThisRedshiftList[ii]:0.2f}')+'.csv'
            df.to_csv(file,index=False)
            #df = pd.read_csv(file)
            #subplot.plot(df['log10_M'],df['Metallicity'], color='black')
        
        
        #Previous model
        if do_previous_model2==1:         
            file = file_previous_model2+'_StellarMetals_vs_StellarMass'+str(f'_z{ThisRedshiftList[ii]:0.2f}')+'.csv'
            df = pd.read_csv(file)
            subplot.plot(df['log10_M'],df['Metallicity'], color='red', linestyle=':')
    
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
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYF19_stellarmetals_vs_stellarmass.pdf')
    plt.close()

    return   
#end metals_vs_stellarmass


def gasmetals_vs_stellarmass(ThisRedshiftList):
      
    #xlim=[7.5,12.0]
    #ylim=[8.0, 9.5] 
    xlim=[9.,11.]
    ylim=[7.8, 9.5] 
    bin=0.2
        
    plot_color=['red','purple']        
   
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.1)) 
    #subplot.yaxis.set_major_locator(MultipleLocator(0.5))    
    subplot.yaxis.set_major_locator(MultipleLocator(0.5))    
    subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
    
    xlab='$\mathrm{log_{10}}(M_*/$'+Msun+'$)$'       
    ylab='$12 + \log_{10}$(O/H)$_{\mathrm{gas}}$'         
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)    
    
   
    
    #for ii in range (0,len(ThisRedshiftList)):
    for ii in range (0,1):    
           
        char_redshift="%0.2f" % ThisRedshiftList[ii]    
        
        
        #Observations
        #KewleyEllison08
        file=Datadir+'Obs_MZgR_KewleyEllison08.fits'
        fits_table=fits.open(file)   
        cols = fits_table[1].columns
        #cols.info() 
        obs = fits_table[1]      
        
        #cols = ['Z94', 'KK04', 'KD02', 'M91', 'D02', 'PP04_O3N2', 'PP04_N2', 'P01', 'P05'] 
        cols = ['T04', 'Z94', 'KK04', 'KD02', 'M91', 'D02', 'PP04_O3N2', 'PP04_N2']         
        plt_colors = [plt.cm.Blues(i) for i in np.linspace(0.8, 0.5,len(cols))]
        
        '''xx = np.arange(8.5, 11.0, 0.05)      
        for jj in range(0, len(cols)):
            sel = obs.data['DIAGNOSTICS'] == cols[jj]
            yy = obs.data['A'][sel] + obs.data['B'][sel]*xx + obs.data['C'][sel]*xx**2 + obs.data['D'][sel]*xx**3
            subplot.plot(xx,yy, color = plt_colors[jj])
            
            if jj < 4:
                xp = 0.05
                yp = 0.85-jj*0.05
            else:
                xp = 0.17
                yp = 0.85-(jj-4)*0.05
            plot_label (subplot, 'label', xlim, ylim,  x_percentage=xp, y_percentage=yp, color=plt_colors[jj], 
                        xlog=0, ylog=0, label=cols[jj], fontsize=10, fontweight='normal')
            
        plot_label (subplot, 'label', xlim, ylim,  x_percentage=0.05, y_percentage=0.92, color='black', 
                    xlog=0, ylog=0, label='Kewley&Ellison2008', fontsize=12, fontweight='normal') ''' 
               
            
        xx = np.arange(8.5, 11.0, 0.05) 
        max_yy = []
        min_yy=[]
        for jj in xx:
            yy = obs.data['A'][:-2] + obs.data['B'][:-2]*jj + obs.data['C'][:-2]*jj**2 + obs.data['D'][:-2]*jj**3  
            max_yy.append(max(yy))
            min_yy.append(min(yy))            
        subplot.fill_between(xx,min_yy,max_yy, facecolor='lightblue', interpolate=True, alpha=0.4)  
        
        for jj in range(0, len(cols)):
            sel = obs.data['DIAGNOSTICS'] == cols[jj]
            yy = obs.data['A'][sel] + obs.data['B'][sel]*xx + obs.data['C'][sel]*xx**2 + obs.data['D'][sel]*xx**3            
            subplot.plot(xx,yy, color = 'cornflowerblue', alpha=0.4)
            
            if jj < 4:
                xp = 0.05
                yp = 0.85-jj*0.05
            else:
                xp = 0.17
                yp = 0.85-(jj-4)*0.05
            plot_label (subplot, 'label', xlim, ylim,  x_percentage=xp, y_percentage=yp, color='cornflowerblue', 
                        xlog=0, ylog=0, label=cols[jj], fontsize=10, fontweight='normal')
        plot_label (subplot, 'label', xlim, ylim,  x_percentage=0.05, y_percentage=0.92, color='black', 
                    xlog=0, ylog=0, label='Kewley&Ellison2008', fontsize=12, fontweight='normal')  
            
        file=Datadir+'Obs_MZgR_Yates+19_fit.fits'
        fits_table=fits.open(file)   
        cols = fits_table[1].columns
        #cols.info()    
        obs = fits_table[1]           
        subplot.plot(obs.data['LOGMSTAR'],obs.data['Z'], color='purple')
        subplot.plot(obs.data['LOGMSTAR'],obs.data['Z_1SIG_UPPER'], color='purple', linestyle='--')
        subplot.plot(obs.data['LOGMSTAR'],obs.data['Z_1SIG_LOWER'], color='purple', linestyle='--')
               
        #file=Datadir+'Obs_MZgR_Yates+12_fit.fits'
        #fits_table=fits.open(file)   
        #cols = fits_table[1].columns
        ##cols.info()    
        #obs = fits_table[1]           
        #subplot.plot(obs.data['LOGMSTAR'],obs.data['Z'], color='purple')
        #subplot.plot(obs.data['LOGMSTAR'],obs.data['Z_1SIG_UPPER'], color='purple', linestyle='--')
        #subplot.plot(obs.data['LOGMSTAR'],obs.data['Z_1SIG_LOWER'], color='purple', linestyle='--')
            
        
        '''#Tremonti2004
        obs_bin=0.5        
        log_mstar=np.arange(xlim[0], xlim[1],obs_bin)
        log_gas_metallicity=-1.492+1.847*log_mstar-0.08026*(log_mstar**2)
        #subplot.plot(log_mstar, log_gas_metallicity,color='black', linewidth=2, linestyle=':')
        subplot.errorbar(log_mstar, log_gas_metallicity, yerr=0.1, color='blue', markeredgecolor='blue', fmt='o')'''
        
        '''#Maiolino
        obs_bin=0.1  
        log_mstar=np.arange(9.5, 11.0,obs_bin)            
        #z=0.1 and 3
        #logM0=np.array([11.18,12.87])
        #K0=np.array([9.04,8.90])
        #z=0.1 and 2.2
        logM0=np.array([11.18,12.38])
        K0=np.array([9.04,8.99])
        gas_metallicity=-0.0864* ((log_mstar-logM0[ii])**2)+K0[ii]       
        subplot.plot(log_mstar, gas_metallicity,color=plot_color[ii], linewidth=2, linestyle=':')'''           
    
        '''#Andrews2012
        file = Datadir + '/brett_gas_metallicity2012.txt' 
        brett12 = Table.read(file, format='ascii')
        err=[-1.*brett12['err_down'],brett12['err_up']]
        subplot.errorbar(brett12['log_mass'], brett12['log_gas_metallicity'], yerr=err, color='black', fmt='o')'''
        
        '''#Zahid2012
        file = Datadir + '/Zahid2012.txt' 
        zahid12 = Table.read(file, format='ascii')
        err=[-1.*zahid12['err'],zahid12['err']]
        subplot.errorbar(zahid12['log_mass'], zahid12['log_gas_metallicity'], yerr=err, color='black', fmt='o')'''
        
        '''#Zahid2014
        log_mstar=np.arange(xlim[0], xlim[1],obs_bin)
        logM0=np.array([9.219,10.06])
        Z0=np.array([9.102,9.08])
        gama=np.array([0.513,0.61])
        gas_metallicity=Z0[ii]+np.log10(1.-exp(-(log_mstar/logM0)**gama[ii]))       
        subplot.plot(log_mstar, gas_metallicity,color='black', linewidth=2, linestyle=':')'''
             
        #create OBSconstraint file for MCMC 
        '''if (ii==0 or ii==1):
            obs_bin=0.5
            log_mstar=np.arange(xlim[0], xlim[1],obs_bin)          
            gas_metallicity=-0.0864* ((log_mstar-logM0[ii])**2)+K0[ii] 
            
           
            obs_y=np.zeros(len(log_mstar),dtype=np.float64)
            obs_y_err=np.zeros(len(log_mstar),dtype=np.float64)           
            for kk in range (0,len(log_mstar)):
                print('%0.2f %0.2f %0.2f %0.2f' % ((log_mstar[kk]+2*np.log10(Hubble_h_WMAP7)-obs_bin/2.)
                                                   ,(log_mstar[kk]+2*np.log10(Hubble_h_WMAP7)+obs_bin/2),
                                                   10**(gas_metallicity[kk]),10**(gas_metallicity[kk])*0.3))  
                      
                obs_y[kk]=10**(gas_metallicity[kk])
                obs_y_err[kk]=10**(gas_metallicity[kk])*0.3          
            obs_y_err = [np.log10(obs_y/(obs_y-obs_y_err)),np.log10((obs_y+obs_y_err)/obs_y)]         
            subplot.errorbar(log_mstar, np.log10(obs_y), yerr=obs_y_err, color='black', fmt='o')'''
        
        
        
        
        #Model
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0=G_MR[sel] 
        
        SSFR = np.log10(G0['Sfr']/(G0['StellarMass']*1.e10/Hubble_h))       
        G0=G0[(G0['StellarMass']>0.) & (G0['Type']<2) & (SSFR>np.log10(2.*(1+ThisRedshiftList[ii])**2./(1.37e10))-1.0)]        
      
        if(opt_detailed_enrichment==0): 
            G0=G0[G0['MetalsColdGas']>0.]
        else:
            G0=G0[(G0['MetalsColdGas'][:,0]+G0['MetalsColdGas'][:,1]+G0['MetalsColdGas'][:,2]>0.)]
        
        if(opt_rings==0):
            G0=G0[np.log10(G0['Sfr'])>-2.0]        
        else:
            G0=G0[np.log10(np.sum(G0['SfrRings'], axis=1))>-2.0]        
        
        
        StellarMass=stellar_mass_with_err(G0, Hubble_h, ThisRedshiftList[ii])                 
        (Metallicity, StellarMass)=get_metallicity(G0, 'ColdGas', StellarMass)        
       
        (x_binned_MR, median_MR, mean_MR, pc16_MR, pc84_MR,rms_MR)=median_and_percentiles(bin, xlim[0], xlim[1], 
                                                                                          StellarMass, Metallicity) 
       
        if(MRII==1):
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, ii, FullSnapshotList_MRII)        
            G0=G_MRII[sel]  
            
            G0=G0[(G0['ColdGas']>0.) & (G0['Type']<2) & (np.log10(G0['Sfr']/(G0['StellarMass']*1.e10/Hubble_h))>-11.0)]        
            if(opt_detailed_enrichment==0): 
                G0=G0[G0['MetalsColdGas']>0.]
            else:
                G0=G0[(G0['MetalsColdGas'][:,0]+G0['MetalsColdGas'][:,1]+G0['MetalsColdGas'][:,2]>0.)]
        
            if(opt_rings==0):
                G0=G0[np.log10(G0['SfrRings'])>-2.0]        
            else:
                G0=G0[np.log10(np.sum(G0['SfrRings'], axis=1))>-2.0]  
            
            StellarMass=stellar_mass_with_err(G0, Hubble_h, ThisRedshiftList[ii])            
            (Metallicity, StellarMass)=get_metallicity(G0, 'ColdGas', StellarMass)    
            sel = Metallicity>0.
            (x_binned_MRII,median_MRII,mean_MRII,pc16_MRII,pc84_MRII,rms_MRII)=median_and_percentiles(bin,xlim[0],xlim[1],
                                                                                       StellarMass[sel],Metallicity[sel]) 
            
                      
            cut_MR_MRII=10.0   
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
        
        #WRITE OUTPUT      
        if(write_to_file==1):
            df = pd.DataFrame({'log10_M':x_binned[sel], 'Metallicity':median[sel], 'pc16':pc16[sel], 'pc84':pc84[sel]})
            file = Datadir+file_to_write+'GasMetals_vs_StellarMass'+str(f'_z{ThisRedshiftList[ii]:0.2f}')+'.csv'
            df.to_csv(file,index=False)
            #df = pd.read_csv(file)
            #subplot.plot(df['log10_M'],df['Metallicity'], color='black')
        
        
        #Previous model
        if do_previous_model2==1:         
            file = file_previous_model2+'_GasMetals_vs_StellarMass'+str(f'_z{ThisRedshiftList[ii]:0.2f}')+'.csv'
            df = pd.read_csv(file)
            subplot.plot(df['log10_M'],df['Metallicity'], color='red', linestyle=':')
            
        
        
        '''#guo10    
        if (ii==0):
            x=[9.25,9.75,10.25,10.75,11.25,11.75,12.25]
            y=[8.71,8.89,9.07,9.15,9.12,9.24,9.40] 
            subplot.plot(x,y,color='blue', linewidth=2)'''
            
        
        #MCMC sample
        if (ii==0):
            if opt_plot_MCMC_sample==1:
                file = MCMCSampledir + 'mcmc_plus_obs_ColdGasMetallicityvsStellarMass_z0.00.txt'            
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
                  
        if (ii==1):
            if opt_plot_MCMC_sample==1:
                file = MCMCSampledir + 'mcmc_plus_obs_ColdGasMetallicityvsStellarMass_z2.00.txt'            
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
        label=prefix_this_model
        plot_label (subplot, 'label', xlim, ylim,  x_percentage=0.49, y_percentage=0.21, color='black', 
                    xlog=0, ylog=0, label=label, fontsize=12, fontweight='normal')             
        plot_label (subplot,'line',xlim,ylim,x_percentage=0.42,y_percentage=0.23,
                    color=plot_color[0],x2_percentage=0.47,xlog=0,ylog=0,linestyle='-',linewidth=2)        
                  
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.49, y_percentage=0.13, 
                    color='black', xlog=0, ylog=0, label=prefix_previous_model2, fontsize=12, fontweight='normal') 
        plot_label (subplot, 'line', xlim, ylim, x_percentage=0.42, y_percentage=0.15, 
                    color='red', x2_percentage=0.47, xlog=0, ylog=0, linestyle=':', linewidth=2)
        
        plot_label (subplot, 'label', xlim, ylim,  x_percentage=0.49, y_percentage=0.05, color='black', 
                    xlog=0, ylog=0, label='T$_\mathrm{e}$-based (Yates2019)', fontsize=12, fontweight='normal')             
        plot_label (subplot,'line',xlim,ylim,x_percentage=0.42,y_percentage=0.07,
                    color='purple',x2_percentage=0.47,xlog=0,ylog=0,linestyle='-',linewidth=2)
        
       
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
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYF19_gasmetals_vs_stellarmass.pdf')
    plt.close()

    return 
#end gasmetals_vs_stellarmass 


def get_metallicity(G0, property_name, StellarMass):
        
    #IF RINGS ARE ON WE CAN MEASURE THINGS WITH THE RIGHT APERTURE (SAME AS SDSS, 30kpc) 
    if(opt_rings==1):
        #SDSS_fibre_size=30
        SDSS_fibre_size=6.1
        sel_Rings=RingRadius<SDSS_fibre_size  
    
    '''The selection on type<2 is needed because these galaxies have their cooling artifically 
    shut down and therefore erroneously enriched by the SNIA and AGB channels without 
    cooling to balance it out'''
       
    #IF there is no Oxigen, we use total metals
    if(property_name=='ColdGas'):
        if(opt_individual_elements==0):              
            if(opt_rings==0):   
                if(opt_detailed_enrichment==0): 
                    MetalsMass=G0['MetalsColdGas']
                else:                    
                    MetalsMass=G0['MetalsColdGas'][:,0]+G0['MetalsColdGas'][:,1]+G0['MetalsColdGas'][:,2]  
                Mass=G0['ColdGas']
            else:  
                sel_Rings=(RingRadius < SDSS_fibre_size) 
                ColdGasRings = G0['ColdGasRings'][:,sel_Rings]
                H2fractionRings = G0['H2fractionRings'][:,sel_Rings]               
                MetalsColdGasRings = (G0['MetalsColdGasRings'][:,sel_Rings,0] + 
                                      G0['MetalsColdGasRings'][:,sel_Rings,1] + 
                                      G0['MetalsColdGasRings'][:,sel_Rings,2])
                
                ColdGasRings[H2fractionRings==0.] = 0.              
                MetalsColdGasRings[H2fractionRings==0.] = 0.
                
                MetalsMass = np.sum(MetalsColdGasRings, axis=1)
                Mass = np.sum(ColdGasRings,axis=1) 
                   
            #Metallicity=np.log10(MetalsMass[sel]/Mass[sel]/0.0134)+8.69  
            #the form below agrees better with individual elements
            Metallicity=np.log10(MetalsMass/Mass/0.0134)+8.69           
            StellarMass=StellarMass         
            
        else:
            
            '''Hydrogen = G0['ColdGasRings_elements'][:,:,0]
            Oxygen = G0['ColdGasRings_elements'][:,:,4]
            
            #SFR weighting
            #SFRfrac[ii,*] = SFRRings[ii,*]/np.max(SFRRings[ii,*],axis=1))
            SFRfrac = SFRRings/np.max(SFRRings[ii,*],axis=1)
            #OH[ii] = TOTAL((oxygen[ii,*]/hydrogen[ii,*])*(1.008/16.0)*SFRfrac[ii,*]) / TOTAL(SFRfrac[ii])'''
         
            
            if(opt_rings==0):  
                N_H=G0['ColdGas_elements'][:,0]/1.
                N_O=G0['ColdGas_elements'][:,4]/16.  
                sel=((N_O>0.) & (N_H>0.))
                N_H = N_H[sel]
                N_O = N_O[sel]
            else:        
                #TO ONLY SELECT SF RINGS (DOENS"T SEEM TO MAKE A DIFFERENCE)               
                OH = np.zeros(len(G0))
                Hydrogen = G0['ColdGasRings_elements'][:,:,0]
                Oxygen = G0['ColdGasRings_elements'][:,:,4]
                SFRRings = G0['SfrRings']
                
                for idx in range (0, len(G0)):  
                    SFRfrac = np.zeros(RNUM)                    
                    SFRfrac = SFRRings[idx,:]/np.max(SFRRings[idx,:])
                    
                    sel = (SFRRings[idx,:]>0.) & (RingRadius<SDSS_fibre_size) & (Hydrogen[idx,:]>0.) & (Oxygen[idx,:]>0.)
                    if(np.sum(SFRfrac[sel])>0):
                        OH[idx] = np.log10(np.sum((Oxygen[idx,sel]/Hydrogen[idx,sel])*(1.008/16.0)*SFRfrac[sel]) /
                                           np.sum(SFRfrac[sel]))
                    else:
                        OH[idx] = -12.
                   
                '''_H=np.zeros(len(G0))
                N_O=np.zeros(len(G0))
                for jj in range (0, len(G0)):                    
                    sel=((G0['H2fractionRings'][jj,:]>0.) & (RingRadius<SDSS_fibre_size))   
                    N_H[jj]=np.sum(G0['ColdGasRings_elements'][jj,sel,0]/1. , axis=0)  
                    N_O[jj]=np.sum(G0['ColdGasRings_elements'][jj,sel,4]/16., axis=0) 
                    
                    
                sel=((N_O>0.) & (N_H>0.))
                N_H = N_H[sel]
                N_O = N_O[sel] ''' 
                                 
            #Metallicity=12.+np.log10(N_O/N_H)           
            Metallicity=12.+OH   
            StellarMass=StellarMass
            
            '''MetalsStellarMassRings = G0['MetalsDiskMassRings'] + G0['MetalsBulgeMassRings']
                    StellarMassRings = G0['DiskMassRings'] + G0['BulgeMassRings']
                    MaxMass = np.max(StellarMassRings[:,sel_Rings],axis=1)
                    
                    MetalsMass=np.sum(MetalsStellarMassRings[:,sel_Rings] /
                                      StellarMassRings[:,sel_Rings] * 
                                      StellarMassRings[:,sel_Rings], axis=1) /  
                                      MaxMass/(np.sum(StellarMassRings[:,sel_Rings])/MaxMass)'''
    else:
        if(property_name=='StellarMass'):
            if(opt_rings==0):   
                if(opt_detailed_enrichment==0): 
                    MetalsMass=G0['MetalsStellarMass']
                else:
                    MetalsMass=G0['MetalsStellarMass'][:,0]+G0['MetalsStellarMass'][:,1]+G0['MetalsStellarMass'][:,2]  
                Mass=G0['StellarMass']
            else:
                if(opt_detailed_enrichment==0): 
                    MetalsMass=np.sum( (G0['MetalsDiskMassRings'][:,sel_Rings] +
                                        G0['MetalsBulgeMassRings'][:,sel_Rings]), axis=1)
                else:
                    MetalsMass=np.sum( (G0['MetalsDiskMassRings'][:,sel_Rings,0] +
                                        G0['MetalsDiskMassRings'][:,sel_Rings,1] +
                                        G0['MetalsDiskMassRings'][:,sel_Rings,2] + 
                                        G0['MetalsBulgeMassRings'][:,sel_Rings,0] +
                                        G0['MetalsBulgeMassRings'][:,sel_Rings,1] +
                                        G0['MetalsBulgeMassRings'][:,sel_Rings,2]), axis=1)
                Mass=np.sum(G0['DiskMassRings'][:,sel_Rings]+G0['BulgeMassRings'][:,sel_Rings],axis=1)
                                   
            sel=((Mass>0.) & (MetalsMass>0.))
            Metallicity=np.log10(MetalsMass[sel]/Mass[sel]/0.0142)            
            StellarMass=StellarMass[sel]            
    #MetalsMassRings
    #MassRings
             
    
    return (Metallicity, StellarMass)




def effective_yield(ThisRedshiftList):
        
    #xlim=[7.5,12.0]
    #ylim=[8.0, 9.5] 
    xlim=[7.,12.0]
    ylim=[-4.0, 0.0] 
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
    
    xlab='$\mathrm{log_{10}}(M_{\mathrm{bar}}/$'+Msun+'$)$'       
    ylab='$\log_{10}(y_{\mathrm{eff}})$'         
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)    
        
    for ii in range (0,len(ThisRedshiftList)):
           
        char_redshift="%0.2f" % ThisRedshiftList[ii]    
            
        #Model
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0=G_MR[sel]             
        '''The selection on type<2 is needed because these galaxies have their cooling artifically 
        shut down and therefore erroneously enriched by the SNIA and AGB channels without 
        cooling to balance it out'''
        if(opt_detailed_enrichment==0): 
            G0=G0[(np.log10(G0['Sfr']/(G0['StellarMass']*1.e10/Hubble_h)) >-11.0) & (G0['ColdGas']>0.) &
                  (G0['MetalsColdGas']>0.) & (G0['Type']<2)]
        else:
            G0=G0[(np.log10(G0['Sfr']/(G0['StellarMass']*1.e10/Hubble_h)) >-11.0) & (G0['ColdGas']>0.) &
                  (G0['MetalsColdGas'][:,0]+G0['MetalsColdGas'][:,1]+G0['MetalsColdGas'][:,2]>0.) & 
                  (G0['Type']<2)]
                
       
        StellarMass=stellar_mass_with_err(G0, Hubble_h, ThisRedshiftList[ii])        
           
        if(opt_rings==0):   
            if(opt_detailed_enrichment==0): 
                MetalsMass=G0['MetalsColdGas']
            else:
                MetalsMass=G0['MetalsColdGas'][:,0]+G0['MetalsColdGas'][:,1]+G0['MetalsColdGas'][:,2]  
            Mass_in_Gas=G0['ColdGas']
            Mass_in_Stars=G0['DiskMass']
        else:
            SDSS_fibre_size=30
            sel_Rings=RingRadius<SDSS_fibre_size  
            if(opt_detailed_enrichment==0): 
                MetalsMass=np.sum(G0['MetalsColdGasRings'][:,sel_Rings], axis=1)
            else:
                MetalsMass=np.sum( (G0['MetalsColdGasRings'][:,sel_Rings,0] +
                                    G0['MetalsColdGasRings'][:,sel_Rings,1] +
                                    G0['MetalsColdGasRings'][:,sel_Rings,2]), axis=1)
            Mass_in_Gas=np.sum(G0['ColdGasRings'][:,sel_Rings],axis=1)  
            Mass_in_Stars=np.sum(G0['DiskMassRings'][:,sel_Rings],axis=1)  
            
        sel=((Mass_in_Gas>0.) & (MetalsMass>0.))
        eff_yield=np.log10((MetalsMass/Mass_in_Gas)/np.log(1./(Mass_in_Gas/(Mass_in_Gas+Mass_in_Stars))))          
        Mbar=np.log10((Mass_in_Gas+Mass_in_Stars)*1.e10/Hubble_h)           
        (x_binned_MR, median_MR, mean_MR, pc16_MR, pc84_MR,rms_MR)=median_and_percentiles (bin, xlim[0], xlim[1], 
                                                                               Mbar[sel], eff_yield[sel])    
       
        if(MRII==1):
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, ii, FullSnapshotList_MRII)        
            G0=G_MRII[sel]    
            if(opt_detailed_enrichment==0): 
                G0=G0[(np.log10(G0['Sfr']/(G0['StellarMass']*1.e10/Hubble_h)) >-11.0) & (G0['ColdGas']>0.) &
                      (G0['MetalsColdGas']>0.) & (G0['Type']<2)]
            else:
                G0=G0[(np.log10(G0['Sfr']/(G0['StellarMass']*1.e10/Hubble_h)) >-11.0) & (G0['ColdGas']>0.) &
                      (G0['MetalsColdGas'][:,0]+G0['MetalsColdGas'][:,1]+G0['MetalsColdGas'][:,2]>0.) & 
                      (G0['Type']<2)]
       
            StellarMass=stellar_mass_with_err(G0, Hubble_h, ThisRedshiftList[ii])  
        
            if(opt_rings==0):   
                if(opt_detailed_enrichment==0): 
                    MetalsMass=G0['MetalsColdGas']
                else:
                    MetalsMass=G0['MetalsColdGas'][:,0]+G0['MetalsColdGas'][:,1]+G0['MetalsColdGas'][:,2]  
                Mass_in_Gas=G0['ColdGas']
                Mass_in_Stars=G0['StellarMass']
            else:            
                if(opt_detailed_enrichment==0): 
                    MetalsMass=np.sum(G0['MetalsColdGasRings'][:,sel_Rings], axis=1)
                else:
                    MetalsMass=np.sum( (G0['MetalsColdGasRings'][:,sel_Rings,0] +
                                        G0['MetalsColdGasRings'][:,sel_Rings,1] +
                                        G0['MetalsColdGasRings'][:,sel_Rings,2]), axis=1)
                Mass_in_Gas=np.sum(G0['ColdGasRings'][:,sel_Rings],axis=1)  
                Mass_in_Stars=np.sum(G0['DiskMassRings'][:,sel_Rings],axis=1)  
            
            sel=((Mass_in_Gas>0.) & (MetalsMass>0.))
            eff_yield=np.log10((MetalsMass/Mass_in_Gas)/ np.log(1./(Mass_in_Gas/(Mass_in_Gas+Mass_in_Stars)))) 
            Mbar=np.log10((Mass_in_Gas+Mass_in_Stars)*1.e10/Hubble_h)           
            (x_binned_MRII,median_MRII,mean_MRII,pc16_MRII,pc84_MRII,rms_MRII)=median_and_percentiles(bin,xlim[0],xlim[1],
                                                                                       Mbar[sel],eff_yield[sel])
                              
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
        
    
      
        #LABELS        
        #z=0
        label=prefix_this_model+', z=0'
        plot_label (subplot, 'label', xlim, ylim, 
                    x_percentage=0.09, y_percentage=0.9, color='black', xlog=0, ylog=0, 
                    label=label, fontsize=13, fontweight='normal')             
        plot_label (subplot,'line',xlim,ylim,x_percentage=0.02,y_percentage=0.92,
                    color=plot_color[0],x2_percentage=0.07,xlog=0,ylog=0,linestyle='-',linewidth=2)
          
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
        
     
    #endfor
    
    plt.tight_layout()
    plt.savefig('./fig/plots_effective_yield.pdf')   
    #pdf.savefig()
    #plt.close()
   
    
    #metallicity vs gas fraction
    
    #xlim=[7.5,12.0]
    #ylim=[8.0, 9.5] 
    xlim=[0.0,0.8]
    ylim=[7.0, 10.0] 
    bin=0.2
        
    plot_color=['red','purple']        
   
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(0.2))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.1)) 
    #subplot.yaxis.set_major_locator(MultipleLocator(0.5))    
    subplot.yaxis.set_major_locator(MultipleLocator(0.5))    
    subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
    
    xlab='$f_{\mathrm{SF,gas}}$'       
    ylab='$12 + \log_{10}$(O/H)$_{gas}$'          
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)    
        
    for ii in range (0,len(ThisRedshiftList)):
           
        char_redshift="%0.2f" % ThisRedshiftList[ii]    
            
        #Model
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0=G_MR[sel]             
        '''The selection on type<2 is needed because these galaxies have their cooling artifically 
        shut down and therefore erroneously enriched by the SNIA and AGB channels without 
        cooling to balance it out'''
        if(opt_detailed_enrichment==0): 
            G0=G0[(np.log10(G0['Sfr']/(G0['StellarMass']*1.e10/Hubble_h)) >-11.0) & (G0['ColdGas']>0.) &
                  (G0['MetalsColdGas']>0.) & (G0['Type']<2)]
        else:
            G0=G0[(np.log10(G0['Sfr']/(G0['StellarMass']*1.e10/Hubble_h)) >-11.0) & (G0['ColdGas']>0.) &
                  (G0['MetalsColdGas'][:,0]+G0['MetalsColdGas'][:,1]+G0['MetalsColdGas'][:,2]>0.) & 
                  (G0['Type']<2)]
                
       
        if(opt_rings==0):   
            if(opt_detailed_enrichment==0): 
                MetalsMass=G0['MetalsColdGas']
            else:
                MetalsMass=G0['MetalsColdGas'][:,0]+G0['MetalsColdGas'][:,1]+G0['MetalsColdGas'][:,2]  
            Mass_in_Gas=G0['ColdGas']
            Mass_in_Stars=G0['DiskMass']
        else:
            SDSS_fibre_size=30
            sel_Rings=RingRadius<SDSS_fibre_size  
            if(opt_detailed_enrichment==0): 
                MetalsMass=np.sum(G0['MetalsColdGasRings'][:,sel_Rings], axis=1)
            else:
                MetalsMass=np.sum( (G0['MetalsColdGasRings'][:,sel_Rings,0] +
                                    G0['MetalsColdGasRings'][:,sel_Rings,1] +
                                    G0['MetalsColdGasRings'][:,sel_Rings,2]), axis=1)
            Mass_in_Gas=np.sum(G0['ColdGasRings'][:,sel_Rings],axis=1)  
            Mass_in_Stars=np.sum(G0['DiskMassRings'][:,sel_Rings],axis=1)  
        
        
        sel=((Mass_in_Gas>0.) & (MetalsMass>0.))
        Metallicity=np.log10(MetalsMass/Mass_in_Gas/0.0134)+8.69  
        GasFraction=Mass_in_Gas/(Mass_in_Gas+Mass_in_Stars)         
                   
        (x_binned_MR, median_MR, mean_MR, pc16_MR, pc84_MR,rms_MR)=median_and_percentiles (bin, xlim[0], xlim[1], 
                                                                               GasFraction[sel], Metallicity[sel])    
       
        if(MRII==1):
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, ii, FullSnapshotList_MRII)        
            G0=G_MRII[sel]    
            if(opt_detailed_enrichment==0): 
                G0=G0[(np.log10(G0['Sfr']/(G0['StellarMass']*1.e10/Hubble_h)) >-11.0) & (G0['ColdGas']>0.) &
                      (G0['MetalsColdGas']>0.) & (G0['Type']<2)]
            else:
                G0=G0[(np.log10(G0['Sfr']/(G0['StellarMass']*1.e10/Hubble_h)) >-11.0) & (G0['ColdGas']>0.) &
                      (G0['MetalsColdGas'][:,0]+G0['MetalsColdGas'][:,1]+G0['MetalsColdGas'][:,2]>0.) & 
                      (G0['Type']<2)]
       
            StellarMass=stellar_mass_with_err(G0, Hubble_h, ThisRedshiftList[ii])  
        
            if(opt_rings==0):   
                if(opt_detailed_enrichment==0): 
                    MetalsMass=G0['MetalsColdGas']
                else:
                    MetalsMass=G0['MetalsColdGas'][:,0]+G0['MetalsColdGas'][:,1]+G0['MetalsColdGas'][:,2]  
                Mass_in_Gas=G0['ColdGas']
                Mass_in_Stars=G0['StellarMass']
            else:            
                if(opt_detailed_enrichment==0): 
                    MetalsMass=np.sum(G0['MetalsColdGasRings'][:,sel_Rings], axis=1)
                else:
                    MetalsMass=np.sum( (G0['MetalsColdGasRings'][:,sel_Rings,0] +
                                        G0['MetalsColdGasRings'][:,sel_Rings,1] +
                                        G0['MetalsColdGasRings'][:,sel_Rings,2]), axis=1)
                Mass_in_Gas=np.sum(G0['ColdGasRings'][:,sel_Rings],axis=1)  
                Mass_in_Stars=np.sum(G0['DiskMassRings'][:,sel_Rings],axis=1)  
            
            sel=((Mass_in_Gas>0.) & (MetalsMass>0.))
            Metallicity=np.log10(MetalsMass/Mass_in_Gas/0.0134)+8.69  
            GasFraction=Mass_in_Gas/(Mass_in_Gas+Mass_in_Stars)         
                   
            (x_binned_MRII,median_MRII,mean_MRII,pc16_MRII,pc84_MRII,rms_MRII)=median_and_percentiles(bin,xlim[0],xlim[1],
                                                                                       GasFraction[sel], Metallicity[sel]) 
                              
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
        
    
      
        #LABELS        
        #z=0
        label=prefix_this_model+', z=0'
        plot_label (subplot, 'label', xlim, ylim, 
                    x_percentage=0.09, y_percentage=0.9, color='black', xlog=0, ylog=0, 
                    label=label, fontsize=13, fontweight='normal')             
        plot_label (subplot,'line',xlim,ylim,x_percentage=0.02,y_percentage=0.92,
                    color=plot_color[0],x2_percentage=0.07,xlog=0,ylog=0,linestyle='-',linewidth=2)
          
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
        
     
    #endfor
     
  
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()

    return   
#end effective_yield 


def iron_clusters(ThisRedshiftList):
      
    
    xlim=[-0.4,1.2]
    ylim=[-1.5, 0.2] 
    bin=0.2
        
    plot_color=['red','purple']        
   
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(0.5))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.1)) 
    #subplot.yaxis.set_major_locator(MultipleLocator(0.5))    
    subplot.yaxis.set_major_locator(MultipleLocator(0.5))    
    subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
    
    xlab='$\log_{10}(\mathrm{kT}_{500}/\mathrm{keV})$'       
    ylab='$\log_{10} (\overline{\mathrm{Z}}_{\mathrm{Fe,500}})$'         
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)    
    
   
    
    #for ii in range (0,len(ThisRedshiftList)):
    for ii in range (0,1):    
           
        char_redshift="%0.2f" % ThisRedshiftList[ii]    
        
        #Model        
        file=Datadir+'GCE_TZR_data_wholeMillI_200619.fits'
        fits_table=fits.open(file)   
        cols = fits_table[1].columns       
        model = fits_table[1]       
        err =[np.log10(model.data['MEANZFE500'])-np.log10(model.data['MEANZFE500']-model.data['MEANZFE500_ERR']),
              np.log10(model.data['MEANZFE500']+model.data['MEANZFE500_ERR'])-np.log10(model.data['MEANZFE500'])]
       
        x = np.log10(model.data['T500'])
        y = np.log10(model.data['MEANZFE500'])
        sel = (y<x*0.09-0.872) | (x>0.65) | ((y>x*-0.41-0.42) & (x>-0.4))
        print(len(x[sel]))
        subplot.errorbar(x[sel], y[sel], fmt='o', markersize=3, ecolor='black', color='black', 
                         elinewidth=1,capsize=2,zorder=-3,alpha=1.0)
        
                
        bin=[0.1,0.1]
        Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]        
        H, xedges, yedges = np.histogram2d(np.log10(model.data['T500']), np.log10(model.data['MEANZFE500']),
                                           bins=Nbins,range=[[xlim[0]-bin[0]/2., xlim[1]+bin[0]/2.], 
                                                  [ylim[0]-bin[1]/2., ylim[1]+bin[1]/2.]])        
        extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
        plt.subplots_adjust(bottom=0.15, left=0.15)     
        mylevels = np.log10(np.logspace(-2.0, 0.0, num=10))  
        H = zoom(H, 20)        
        H=np.log10(H/np.amax(H))       
        cont=plt.contourf(H.transpose()[::], origin='lower', cmap='Greys_r', levels=mylevels, extent=extent)  
            
        plt.colorbar(format='%0.1f') 
        
        
        #Observations        
        file=Datadir+'Obs_TZR_Yates+17.fits'
        fits_table=fits.open(file)   
        cols = fits_table[1].columns
        #print(cols.info())    
        obs = fits_table[1]       
        obs_x_err = [-obs.data['LOGT500_LOWER_ERR'], obs.data['LOGT500_UPPER_ERR']]
        obs_y_err = [-obs.data['LOGFEH500_LOWER_ERR'],obs.data['LOGFEH500_UPPER_ERR']]
        subplot.errorbar(obs.data['LOGT500'], obs.data['LOGFEH500'],xerr= obs_x_err, yerr=obs_y_err,
                         fmt='o', markersize=3, ecolor='darkorange', color='darkorange',  elinewidth=1,capsize=2,alpha=0.9)
        
        file=Datadir+'Obs_TZR_Mernier+18.fits'
        fits_table=fits.open(file)   
        cols = fits_table[1].columns
        #print(cols.info())    
        obs = fits_table[1]       
        obs_x_err = [-obs.data['LOGT500_LOWER_ERR'], obs.data['LOGT500_UPPER_ERR']]
        obs_y_err = [-obs.data['LOGFEH500_LOWER_ERR'],obs.data['LOGFEH500_UPPER_ERR']]
        subplot.errorbar(obs.data['LOGT500'], obs.data['LOGFEH500'],xerr= obs_x_err, yerr=obs_y_err,
                         fmt='o', markersize=3, ecolor='steelblue', color='royalblue',  elinewidth=1,capsize=2,alpha=0.9)    
        
        
       
        
        
        #LABELS        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.575, y_percentage=0.9, 
                    color='black', xlog=0, ylog=0, label='Mernier 2018', 
                    fontsize=12, fontweight='normal') 
        plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.55, y_percentage=0.925, 
                    color='steelblue', xlog=0, ylog=0, sym='o', sym_size=3, err_size=0.03,alpha=0.9)
        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.575, y_percentage=0.84, 
                    color='black', xlog=0, ylog=0, label='Yates 2017', 
                    fontsize=12, fontweight='normal') 
        plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.55, y_percentage=0.865, 
                    color='darkorange', xlog=0, ylog=0, sym='o', sym_size=3, err_size=0.03,alpha=0.9)
        
        #plot_label (subplot, 'label', xlim, ylim, x_percentage=0.575, y_percentage=0.78, 
        #            color='black', xlog=0, ylog=0, label=prefix_this_model, 
        #            fontsize=12, fontweight='normal') 
        #plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.55, y_percentage=0.80, 
        #            color='black', xlog=0, ylog=0, sym='o', sym_size=3, err_size=0.001,alpha=1.0)
        
       
    #endfor
     
    
    
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYF19_iron_clusters.pdf')
    plt.close()

    return 
#end iron_clusters 


def metals_mass(ThisRedshiftList):
     
    #xlim=[7.5,12.0]
    #ylim=[8.0, 9.5] 
    xlim=[7.,12.0]
    ylim=[5., 12.] 
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
    
    xlab='$\mathrm{log_{10}}(M_*/$'+Msun+'$)$'       
    ylab='$\mathrm{log_{10}}(M/$'+Msun+'$)$'         
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)    
        
    for ii in range (0,len(ThisRedshiftList)):
           
        char_redshift="%0.2f" % ThisRedshiftList[ii]    
            
        #Model
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0=G_MR[sel]             
        '''The selection on type<2 is needed because these galaxies have their cooling artifically 
        shut down and therefore erroneously enriched by the SNIA and AGB channels without 
        cooling to balance it out'''
        if(opt_detailed_enrichment==0): 
            G0=G0[(np.log10(G0['Sfr']/(G0['StellarMass']*1.e10/Hubble_h)) >-11.0) & (G0['ColdGas']>0.) &
                  (G0['MetalsColdGas']>0.) & (G0['Type']<2)]
        else:
            G0=G0[(np.log10(G0['Sfr']/(G0['StellarMass']*1.e10/Hubble_h)) >-11.0) & (G0['ColdGas']>0.) &
                  (G0['MetalsColdGas'][:,0]+G0['MetalsColdGas'][:,1]+G0['MetalsColdGas'][:,2]>0.) & 
                  (G0['Type']<2)]
                
       
        StellarMass=stellar_mass_with_err(G0, Hubble_h, ThisRedshiftList[ii])  
        #MetalsMass=G0['MetalsColdGas'][:,0]+G0['MetalsColdGas'][:,1]+G0['MetalsColdGas'][:,2]  
        #Mass=G0['ColdGas']
        #Metallicity=np.log10(MetalsMass/Mass/0.0134)+8.69         
        #(x_binned_MR, median_MR, mean_MR, pc16_MR, pc84_MR,rms_MR)=median_and_percentiles (bin, xlim[0], xlim[1], 
        #                                                                       StellarMass, Metallicity)    
       
        StellarMass=stellar_mass_with_err(G0, Hubble_h, ThisRedshiftList[ii])           
        Mass=np.log10(G0['ColdGas']*1.e10/Hubble_h) 
        Metals=np.log10((G0['MetalsColdGas'][:,0]+G0['MetalsColdGas'][:,1]+G0['MetalsColdGas'][:,2])*1.e10/Hubble_h)
        (x_binned_MR, median_MR, mean, pc16, pc84,rms)=median_and_percentiles (bin,xlim[0],xlim[1],StellarMass,Mass) 
        (x_binned2_MR, median2_MR, mean, pc16, pc84,rms)=median_and_percentiles (bin,xlim[0],xlim[1],StellarMass,Metals) 
        
        (sel)=select_current_redshift(G_MRII, ThisRedshiftList, ii, FullSnapshotList_MRII)        
        G0=G_MRII[sel]  
        if(opt_detailed_enrichment==0): 
            G0=G0[(np.log10(G0['Sfr']/(G0['StellarMass']*1.e10/Hubble_h)) >-11.0) & (G0['ColdGas']>0.) &
                  (G0['MetalsColdGas']>0.) & (G0['Type']<2)]
        else:
            G0=G0[(np.log10(G0['Sfr']/(G0['StellarMass']*1.e10/Hubble_h)) >-11.0) & (G0['ColdGas']>0.) &
                  (G0['MetalsColdGas'][:,0]+G0['MetalsColdGas'][:,1]+G0['MetalsColdGas'][:,2]>0.) & 
                  (G0['Type']<2)]                
       
        StellarMass=stellar_mass_with_err(G0, Hubble_h, ThisRedshiftList[ii])           
        Mass=np.log10(G0['ColdGas']*1.e10/Hubble_h) 
        Metals=np.log10((G0['MetalsColdGas'][:,0]+G0['MetalsColdGas'][:,1]+G0['MetalsColdGas'][:,2])*1.e10/Hubble_h)
        (x_binned_MRII, median_MRII,mean,pc16,pc84,rms)=median_and_percentiles (bin,xlim[0], xlim[1], StellarMass, Mass) 
        (x_binned2_MRII, median2_MRII,mean,pc16,pc84,rms)=median_and_percentiles(bin,xlim[0],xlim[1],StellarMass,Metals) 
       
        
        cut_MR_MRII=9.0    
        sel=x_binned_MRII<cut_MR_MRII
        x_binned_MRII=x_binned_MRII[sel]
        median_MRII=median_MRII[sel] 
        x_binned2_MRII=x_binned2_MRII[sel]
        median2_MRII=median2_MRII[sel]
           
        sel=x_binned_MR>cut_MR_MRII
        x_binned_MR=x_binned_MR[sel]
        median_MR=median_MR[sel] 
        x_binned2_MR=x_binned2_MR[sel]
        median2_MR=median2_MR[sel] 
          
        x_binned=np.concatenate((x_binned_MRII,x_binned_MR), axis=0)
        median=np.concatenate((median_MRII,median_MR), axis=0)
        x_binned2=np.concatenate((x_binned2_MRII,x_binned2_MR), axis=0)
        median2=np.concatenate((median2_MRII,median2_MR), axis=0)
                  
        subplot.plot(x_binned, median,color=plot_color[ii], linewidth=2)             
        subplot.plot(x_binned2, median2,color=plot_color[ii], linewidth=2, linestyle='--')
       
        
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
     
        
        
      
     
    
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()

    return  
#end metals_vs_mass 



def test_metal_evo(ThisRedshiftList):
    
    xlim=[8.,12.0]
    ylim=[7.5, 9.5] 
    bin=0.2
        
    plot_color=['darkred', 'red', 'orange', 'yellow']        
    
    fig = plt.figure(figsize=(two_two_size_large[0],two_two_size_large[1]))
    grid = gridspec.GridSpec(2, 2)
  
    #plot different colours
    #colormap = plt.cm.gist_ncar
    colormap = plt.cm.Blues
    plot_color = [colormap(i) for i in np.linspace(0.9, 0.3, len(ThisRedshiftList))]
               
  
    for ii in range (0,len(ThisRedshiftList)):
           
        char_redshift="%0.2f" % ThisRedshiftList[ii]    
            
        #Model
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0=G_MR[sel]             
        
        if(opt_detailed_enrichment==0): 
            G0=G0[(np.log10(G0['Sfr']/(G0['StellarMass']*1.e10/Hubble_h)) >-11.0) & (G0['ColdGas']>0.) &
                  (G0['MetalsColdGas']>0.) & (G0['Type']<2)]
        else:
            G0=G0[(np.log10(G0['Sfr']/(G0['StellarMass']*1.e10/Hubble_h)) >-11.0) & (G0['ColdGas']>0.) &
                  (G0['MetalsColdGas'][:,0]+G0['MetalsColdGas'][:,1]+G0['MetalsColdGas'][:,2]>0.) & 
                  (G0['Type']<2)]
      
       
        G0=G0[(np.log10(G0['Sfr']/(G0['StellarMass']*1.e10/Hubble_h)) >-11.0) & (G0['ColdGas']>0.)]
              
        for jj in range(0, 4):
            
            StellarMass=stellar_mass_with_err(G0, Hubble_h, ThisRedshiftList[ii])  
                      
            if(jj==0):
                ylim=[7.8, 9.2] 
                ylab='$12 + \log_{10}$(O/H)$_{gas}$'         
                (Metallicity, StellarMass)=get_metallicity(G0, 'ColdGas', StellarMass)
                property_yy =  Metallicity
                
            if(jj==1):
                ylim=[-1.5, 0.5]   
                ylab='$Z_*$'         
                (Metallicity, StellarMass)=get_metallicity(G0, 'StellarMass', StellarMass)
                property_yy =  Metallicity
        
            if(jj==2): 
                ylim=[-2.0, 2.5]  
                ylab='Sfr'
                property_yy =  np.log10(G0['Sfr'])
                
            if(jj==3):  
                ylim=[0.0, 1.]  
                ylab='Gas Fraction'
                property_yy =  G0['ColdGas']/(G0['ColdGas']+G0['StellarMass'])     
                
            subplot=plt.subplot(grid[jj])            
            subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
            xlab='$\mathrm{log_{10}}(M_*/$'+Msun+'$)$'             
            subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)    
                           
            (x_binned, median, mean, pc16, pc84,rms)=median_and_percentiles (bin, xlim[0], xlim[1], StellarMass, property_yy) 
            if(jj==0):
                sel = median>0.  
            else:
                sel = median>-100.
        
            #subplot.plot(x_binned[sel], median[sel], color=plot_color[ii], linewidth=2)      
            #subplot.plot(x_binned[sel], pc16[sel], color=plot_color[ii], linewidth=2, linestyle='--')
            #subplot.plot(x_binned[sel], pc84[sel], color=plot_color[ii], linewidth=2, linestyle='--')
            subplot.plot(x_binned[sel], median[sel], color=plot_color[ii],linewidth=2)      
            #subplot.plot(x_binned[sel], pc16[sel], linewidth=2, linestyle='--')
            #subplot.plot(x_binned[sel], pc84[sel], linewidth=2, linestyle='--')
         
            if(jj==0):
                #Observations
                obs_bin=0.1
                #Maiolino
                if(ii<2):    
                #if(ii==0):    
                    log_mstar=np.arange(9.5, 11.0,obs_bin)            
                    #z=0.1 and 3
                    #logM0=np.array([11.18,12.87])
                    #K0=np.array([9.04,8.90])
                    #z=0.1 and 2.2
                    logM0=np.array([11.18,12.38])
                    K0=np.array([9.04,8.99])
                    gas_metallicity=-0.0864* ((log_mstar-logM0[ii])**2)+K0[ii]       
                    subplot.plot(log_mstar, gas_metallicity,color='red', linewidth=2, linestyle=':')
        
    
    
    
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()

    return 
#end test_metal_evo 


def BHBM(ThisRedshiftList):
   
    for ii in range(0,len(ThisRedshiftList)): 
        
        xlim=[8.5,12.5]
        #xlim=[11.,14.]
        ylim=[5.0, 10.5]
        #ylim=[4.0, 10.5]
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
            
        xlab='$\mathrm{log_{10}}(M_{\mathrm{Bulge}}/$'+Msun+'$)$'   
        #xlab='$\mathrm{log_{10}}(M_{*}/$'+Msun+'$)$'   
        #xlab='$\mathrm{log_{10}}(M_{vir}/$'+Msun+'$)$'   
        ylab='$\mathrm{log_{10}}(M_{\mathrm{BH}}/$'+Msun+'$)$' 
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   
            
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR_unsel=G_MR[sel]   
        G0_MR=G0_MR_unsel[(G0_MR_unsel['BulgeMass'] > 0.) & (G0_MR_unsel['BlackHoleMass'] > 0.)]
        Ngals=len(G0_MR) 
       
        #G0_MR_random = np.random.choice(G0_MR,100000)    
        #BulgeMass=(np.log10(G0_MR_random['BulgeMass']*1.e10/Hubble_h))        
        #BHMass=(np.log10(G0_MR_random['BlackHoleMass']*1.e10/Hubble_h)) 
        #plt.scatter(BulgeMass, BHMass, s=5, color='black')   
        
        
        BulgeMass=(np.log10(G0_MR['BulgeMass']*1.e10/Hubble_h)) 
        #BulgeMass=(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)) 
        #BulgeMass=(np.log10(G0_MR['Mvir']*1.e10/Hubble_h)) 
        BHMass=(np.log10(G0_MR['BlackHoleMass']*1.e10/Hubble_h))                     
        #plt.scatter(BulgeMass, BHMass, s=5, color='black')          
        
        #plot all point at high mass
        #sel = (BulgeMass>11.9) | ((BulgeMass<11.9) & ((BHMass<2.67*BulgeMass-23.0) | (BHMass>0.71*BulgeMass+1.13)))      
        sel = ( (BulgeMass>11.9) | 
                ((BulgeMass<11.9) & (BulgeMass>11.2) & (BHMass>6.0) & (BHMass<2.67*BulgeMass-23.0)) |
                ((BulgeMass>8.5) & (BulgeMass<11.9) &(BHMass>0.71*BulgeMass+1.13)) |
                ((BulgeMass<11.5) & (BHMass<1.6*BulgeMass-12.0)) 
              )
      
        print(len(BulgeMass[sel]))        
        plt.scatter(BulgeMass[sel], BHMass[sel], s=2, color='black')  
               
        H, xedges, yedges = np.histogram2d(BulgeMass, BHMass, bins=Nbins,
                                           range=[[xlim[0]-bin[0]/2., xlim[1]+bin[0]/2.], 
                                                  [ylim[0]-bin[1]/2., ylim[1]+bin[1]/2.]])            
        extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
        plt.subplots_adjust(bottom=0.15, left=0.15)     
        mylevels = np.log10(np.logspace(-3.5, 0.0, num=10))  
        H = zoom(H, 20)        
        H=np.log10(H/np.amax(H))       
        cont=plt.contourf(H.transpose()[::], origin='lower', cmap='Greys_r', levels=mylevels, extent=extent)  
        #cont=plt.contourf(H.transpose()[::], origin='lower', cmap='rainbow', levels=mylevels, extent=extent)     
        plt.colorbar(format='%0.1f') 
        #(ax, cmap=None, norm=None, alpha=None, values=None, boundaries=None, orientation='vertical', 
        #ticklocation='auto', extend='neither', spacing='uniform', ticks=None, format=None, 
        #drawedges=False, filled=True, extendfrac=None, extendrect=False, label='')
        #plt.scatter(BulgeMass[sel], BHMass[sel], s=2, color='red',alpha=0.4)    
     
        
        file = Datadir + 'mcconnel2012.dat'
        obs = Table.read(file, format='ascii', data_start=20)     
        obs_x = np.log10(obs['col14'])
        obs_y = np.log10(obs['col3'])        
        obs_x_err=np.zeros(len(obs_x),dtype=np.float64)+0.24 
        obs_y_err = [np.log10(obs['col3']/obs['col4']),np.log10(obs['col5']/obs['col3'])]
       
        subplot.errorbar(obs_x, obs_y,xerr= obs_x_err, yerr=obs_y_err,
                 fmt='o', markersize=5, ecolor='blue', color='blue',capsize=2)
       
        '''x_arr=np.arange(xlim[0],xlim[1],0.01)
        y_arr=-0.952381*x_arr+17.2        
        subplot.plot(x_arr,y_arr)
        y_arr=-0.952381*x_arr+18.88
        subplot.plot(x_arr,y_arr)
        y_arr=1.05*x_arr-2.91961
        subplot.plot(x_arr,y_arr)'''
        
        #LABELS        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.9, 
                    color='black', xlog=0, ylog=0, label='McConnell2012', 
                    fontsize=13, fontweight='normal') 
        plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.925, 
                    color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15) 
       
        
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    #plt.savefig('./fig/HWL17_bhbm.pdf')
    plt.savefig('./fig/HYF19_bhbm.pdf')
    plt.close()

    return   
#end BHBM


def BHMvir(ThisRedshiftList):
   
    for ii in range(0,len(ThisRedshiftList)):        
        
        xlim=[8.5,15.5]
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
            
        xlab='$\mathrm{log_{10}}(M_{\mathrm{Bulge}}/$'+Msun+'$)$'       
        ylab='$\mathrm{log_{10}}(M_{\mathrm{BH}}/$'+Msun+'$)$' 
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   
            
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR_unsel=G_MR[sel]   
        G0_MR=G0_MR_unsel[(G0_MR_unsel['BulgeMass'] > 0.) & (G0_MR_unsel['BlackHoleMass'] > 0.)]
        Ngals=len(G0_MR) 
       
        BulgeMass=(np.log10(G0_MR['Mvir']*1.e10/Hubble_h)) 
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
          
        #G0_MR=np.random.choice(G0_MR, size=2000)           
        #BulgeMass=(np.log10(G0_MR['BulgeMass']*1.e10/Hubble_h)) 
        #BHMass=(np.log10(G0_MR['BlackHoleMass']*1.e10/Hubble_h))     
        #plt.scatter(BulgeMass, BHMass, s=5, color='black')          
        
        
        (slope,b)=get_slope(11.,5.5,13.,7.8)
        print(slope,b)
        x_arr=np.arange(xlim[0],xlim[1],0.1)
        subplot.plot(x_arr,x_arr*slope+b, color='red', linewidth=2, linestyle='-')           
        
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
                    color='black', xlog=0, ylog=0, label='BHMvir', 
                    fontsize=13, fontweight='normal')
        
    plt.tight_layout()   
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()

    return  
#end BHMvir



def SFRF(ThisRedshiftList):
    
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
            xlab='$\mathrm{log_{10}}(\mathrm{SFR}/($'+Msun+'$yr^{-1}))$' 
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
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYJ18_sfrf.pdf')
    plt.close()

    return 
#end SFRF

def SFRD(ThisRedshiftList):
       
    xlim=[-1.0,1.0]
    ylim=[-2.5, 0.0]
    bin=[0.25,0.25]
    Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]
    
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
            
            
    xlab='$z$'       
    ylab=r'$\log_{10}(\rho_{\mathrm{SFR}}/($'+Msun+'$\mathrm{yr}^{-1}\mathrm{Mpc}^{-3}))$'     
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   
    #subplot.set_xscale('symlog')
    
    #format axis
    #majorFormatter = FormatStrFormatter('%d')
    #majorFormatter = FormatStrFormatter('%d')
    
    #from matplotlib.ticker import StrMethodFormatter, NullFormatter
    #subplot.xaxis.set_major_formatter(StrMethodFormatter('{x:0.0f}'))
    #subplot.xaxis.set_minor_formatter(StrMethodFormatter('{x:0.0f}'))
    
    subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
    subplot.yaxis.set_minor_locator(MultipleLocator(0.1))    
    #subplot.xaxis.set_major_locator(MultipleLocator(0.5))    
    #subplot.xaxis.set_minor_locator(MultipleLocator(1.0)) 
    
    arr = [0.1, 0.2, 0.5, 1, 2, 5, 10]     
    subplot.set_xticks(np.log10(np.array(arr)))
    subplot.set_xticklabels([str(x) for x in arr])    
    arr = [0.09,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0] 
    subplot.set_xticks(np.log10(np.array(arr)), minor=True)
        
    #OBSERVATIONS    
    obs_z=np.arange(10**xlim[0]-1.0+0.001,10**xlim[1],0.1)    
    A=10**(-0.997*(obs_z-1.243))
    B=10**(0.241*(obs_z-1.243))
    obs_sfrd=0.180/(A+B)
    
    obs_z_err=obs_z/obs_z-1.0   
    obs_z_err[obs_z < 0.9]+=0.13   
    obs_z_err[(obs_z > 0.9) & (obs_z < 1.7)]+=0.17    
    obs_z_err[(obs_z > 1.7) & (obs_z < 3.)]+=0.19    
    obs_z_err[(obs_z > 3.) & (obs_z < 12.)]+=0.27

    #subplot.plot((obs_z+1), np.log10(obs_sfrd),color='gray', linewidth=2)    
    obs_z_err=smooth(obs_z_err,1000)   
    subplot.fill_between(np.log10(obs_z),np.log10(obs_sfrd)-obs_z_err,np.log10(obs_sfrd)+obs_z_err,
                         facecolor='grey', interpolate=True, alpha=0.3)  
       
        
        
    #Bouwens2012   
    obs_sfrd=[-1.12,-1.51,-1.71,-1.90,-2.13]
    obs_sfrd_err=[0.05,0.06,0.08,0.10,0.11]
    obs_z=[3.8,5.,5.9,6.8,8.0]
    subplot.errorbar(np.log10(obs_z),obs_sfrd, yerr=obs_sfrd_err, fmt='o', markersize=5, mfc='darkviolet',
                     markeredgecolor='darkviolet', color='darkviolet',capsize=2)
    
    #driver2018
    obs_redshift = np.array([0.05, 0.10, 0.17, 0.24, 0.32, 0.40, 0.50, 0.62, 0.75, 0.91, 
                             1.10, 1.32, 1.60, 1.98, 2.40, 2.93, 3.50, 4.00, 4.63])
    obs_sfrd     = np.array([-1.95, -1.82, -1.90, -1.77, -1.75, -1.79, -1.73, -1.56, -1.42, -1.29, 
                             -1.31, -1.27, -1.17, -1.30, -1.29, -1.28, -1.33, -1.42, -1.45])

    #already subtracted from SFRD
    error_eddington = np.array([0.03, 0.03, 0.02, 0.01, 0.01, 0.01, 0.04, 0.05, 0.06, 0.05,
                                0.04, 0.03, 0.02, 0.04, 0.04, 0.04, 0.03, 0.04, 0.03])
    error_poisson   = np.array([0.00, 0.01, 0.00, 0.00, 0.00, 0.01, 0.01, 0.00, 0.01, 0.00, 
                                0.00, 0.00, 0.00, 0.01, 0.01, 0.01, 0.01, 0.04, 0.04]  )
    error_cv        = np.array([0.07, 0.05, 0.04, 0.05, 0.06, 0.06, 0.09, 0.07, 0.06, 0.07,
                                0.05, 0.06, 0.06, 0.07, 0.04, 0.04, 0.03, 0.05, 0.04])
    error_agn       = np.array([0.00, 0.01, 0.00, 0.00, 0.01, 0.01, 0.03, 0.02, 0.04, 0.01,
                                0.01, 0.02, 0.03, 0.06, 0.09, 0.11, 0.08, 0.02, 0.04])

    subplot.errorbar(np.log10(obs_redshift), obs_sfrd, yerr=(error_poisson+error_cv+error_agn), fmt='o', markersize=5, 
                             markeredgecolor='blue', color='blue',zorder=+3,capsize=2)
    
    
    model_z=np.zeros(len(ThisRedshiftList),dtype=np.float32)
    model_sfrd=np.zeros(len(ThisRedshiftList),dtype=np.float32)
    
    for ii in range(0,len(ThisRedshiftList)):       
        
        if(MRII==1):
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, ii, FullSnapshotList_MRII)        
            G0_MR_unsel=G_MRII[sel]  
        else:
            (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
            G0_MR_unsel=G_MR[sel]  
            
        G0_MR=G0_MR_unsel[(G0_MR_unsel['StellarMass'] > 0.00001) & ~np.isnan(G0_MR_unsel['Sfr'])]
        
        model_sfrd[ii]=np.sum(G0_MR['Sfr'], axis=0)/Volume_MRII        
        model_z[ii]=ThisRedshiftList[ii]
              
    subplot.plot(np.log10(model_z), np.log10(model_sfrd),color='red', linewidth=2, linestyle='-')
    #subplot.plot(model_z, np.log10(model_sfrd),color='red', linewidth=2, linestyle='-')
    
    #WRITE OUTPUT      
    if(write_to_file==1):
        df = pd.DataFrame({'log10_z':np.log10(model_z), 'log10_SFRD':np.log10(model_sfrd)})
        file = Datadir+file_to_write+'SFRD_MRII.csv'
        df.to_csv(file,index=False)
        #df = pd.read_csv(file)
        #subplot.plot(df['log10_z'],df['log10_SFRD'], color='black')
        
    for ii in range(0,len(ThisRedshiftList)):            
            
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR_unsel=G_MR[sel]   
        G0_MR=G0_MR_unsel[(G0_MR_unsel['StellarMass'] > 0.00001) & ~np.isnan(G0_MR_unsel['Sfr'])]
        Ngals=len(G0_MR) 
       
        model_sfrd[ii]=np.sum(G0_MR['Sfr']*10**(np.random.randn(len(G0_MR['Sfr']))*0.1), axis=0)/Volume_MR
        model_z[ii]=ThisRedshiftList[ii]
       
    subplot.plot(np.log10(model_z), np.log10(model_sfrd),color='red', linewidth=2,linestyle='--')
    
    
    #WRITE OUTPUT      
    if(write_to_file==1):
        df = pd.DataFrame({'log10_z':np.log10(model_z), 'log10_SFRD':np.log10(model_sfrd)})
        file = Datadir+file_to_write+'SFRD_MR.csv'
        df.to_csv(file,index=False)
        #df = pd.read_csv(file)
        #subplot.plot(df['log10_z'],df['log10_SFRD'], color='black')
                
    '''fa = open(Datadir+"SFRD_noSN.txt", "w")    
    for kk in range (0,len(model_sfrd)):                     
        fa.write("%0.2f " % model_z[kk] + "%0.2f\n" % np.log10(model_sfrd[kk]))         
    fa.close()'''
    
    '''file=Datadir+"SFRD_noAGN.txt"
    if os.path.isfile(file):
        model = Table.read(file, format='ascii')    
    model = Table.read(file, format='ascii')      
    subplot.plot(model['col1']+1,model['col2'], color='blue', linewidth=2,linestyle='-')
    
    file=Datadir+"SFRD_noSN.txt"
    if os.path.isfile(file):
        model = Table.read(file, format='ascii')    
    model = Table.read(file, format='ascii')      
    subplot.plot(model['col1']+1,model['col2'], color='brown', linewidth=2,linestyle='-')
    
    file=Datadir+"SFRD_nofeedback.txt"
    if os.path.isfile(file):
        model = Table.read(file, format='ascii')    
    model = Table.read(file, format='ascii')      
    subplot.plot(model['col1']+1,model['col2'], color='black', linewidth=2,linestyle='-')'''
    
    
    '''plot_label (subplot, 'label', xlim, ylim, x_percentage=0.12, y_percentage=0.2, 
                color='black', xlog=0, ylog=0, label=prefix_previous_model2, fontsize=10, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, y_percentage=0.22, 
                color='red', x2_percentage=0.1, xlog=0, ylog=0, linestyle=linestyle_previous_model2, linewidth=2)'''
    
    #PREVIOUS MODELS           
    if do_previous_model2==1: 
        file = file_previous_model2+'_sfrd.txt' 
        model = Table.read(file, format='ascii')
        subplot.plot(np.log10(model['col1']-1),np.log10(model['col2']*(Hubble_h)),color='red',
                     linestyle=':', linewidth=2) 
    
    #LABELS        
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.11, y_percentage=0.9, 
                color='black', xlog=0, ylog=0, label='Driver2018', fontsize=12, fontweight='normal') 
    plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.07, y_percentage=0.92, 
                color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.05) 
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.11, y_percentage=0.83, 
                color='black', xlog=0, ylog=0, label='Bouwens2012 ', fontsize=12, fontweight='normal') 
    plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.07, y_percentage=0.85, 
                color='darkviolet', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.05) 
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.11, y_percentage=0.76, color='black', 
                xlog=0, ylog=0, label='Behroozi2013', fontsize=12, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.01, x2_percentage=0.08, y_percentage=0.78, 
                color='grey', alpha=0.3, xlog=0, ylog=0, linestyle='-', linewidth=6)
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.11, y_percentage=0.69, color='black', 
                xlog=0, ylog=0, label=prefix_this_model+' MSII', fontsize=12, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.01, x2_percentage=0.08, y_percentage=0.71, 
                color='red', xlog=0, ylog=0, linestyle='-', linewidth=2)
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.11, y_percentage=0.62, color='black', 
                xlog=0, ylog=0, label=prefix_this_model+' MS', fontsize=12, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.01, x2_percentage=0.08, y_percentage=0.64, 
                color='red', xlog=0, ylog=0, linestyle='--', linewidth=2)
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.11, y_percentage=0.55, color='black', 
                xlog=0, ylog=0, label=prefix_previous_model2, fontsize=12, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.01, x2_percentage=0.08, y_percentage=0.57, 
                color='red', xlog=0, ylog=0, linestyle=':', linewidth=2)
    
    
    
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYF19_SFRD.pdf')
    plt.close()

    return     
#end SFRD




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
            
            

def gas_fraction(ThisRedshiftList):
   
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
            
        xlab='$log_{10}(M_*/$'+Msun+'$)$'
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
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()

    return    
#end gas fraction


def HI_over_Lr_vs_HI(ThisRedshiftList):
   
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
            
        xlab='$\mathrm{log_{10}}(M_{\mathrm{HI}}/$'+Msun+'$)$'
        ylab='$\mathrm{log_{10}}(M_{\mathrm{HI}}/L_{\mathrm{r}})$'     
        subplot.set_xlabel(xlab, fontsize=15), subplot.set_ylabel(ylab, fontsize=15)   
            
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
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()

    return 
#end gas fraction

def HI_over_Lr_vs_HI_bins(ThisRedshiftList):
   
    for ii in range(0,len(ThisRedshiftList)):        
        
        char_redshift="%0.2f" % ThisRedshiftList[ii]
        
        xlim=[-1.4,1.4]
        ylim=[0.0,0.7]
      
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
            subplot.yaxis.set_minor_locator(MultipleLocator(0.05))  
            subplot.yaxis.set_major_locator(MultipleLocator(0.2))
                        
            xlab='$\mathrm{log_{10}}((M_{\mathrm{HI}}/L_{\mathrm{r}})/($'+Msun+'$/$'+Lsun+'$))$'
            ylab='fraction'   
            if jj==1 or jj == 3 or jj==5:
                plt.tick_params(axis='y', which='both', left='on', labelleft='off')
                ylab=''              
            subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   

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
            obs_y=hist[0]*1./np.sum(hist[0])
            obs_y_err=np.sqrt(hist[0])/np.sum(hist[0])
            subplot.errorbar(obs_x,obs_y, obs_y_err, fmt='o', markersize=5, ecolor='blue', color='blue',zorder=+3,capsize=2)
          
        
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
            
            
            #MODEL
            if((log_MHI[jj]<9.5) & (MRII==1)):
                (sel)=select_current_redshift(G_MRII, ThisRedshiftList, ii, FullSnapshotList_MRII)        
                G0=G_MRII[sel]  
            else:
                (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
                G0=G_MR[sel]
                       
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
            y_axis=hist[0]*1./np.sum(hist[0])     
            subplot.plot(x_axis,y_axis, color='red', linewidth=2, linestyle='-') 
                   
               
            #WRITE OUTPUT      
            if(write_to_file==1):
                df = pd.DataFrame({'log10_HI_r':x_axis, 'Fraction':y_axis})
                file = Datadir + file_to_write + 'HI_over_r_bins' + str(f'_MHI{log_MHI[jj]-bin_HI/2:0.2f}') + \
                       str(f'_{log_MHI[jj]+bin_HI/2:0.2f}') + str(f'_z{ThisRedshiftList[ii]:0.2f}')+'.csv'
                df.to_csv(file,index=False)
                #df = pd.read_csv(file)
                #subplot.plot(df['log10_HI_r'],df['Fraction'], color='black')
            
            
            #Previous Model
            if do_previous_model2==1:                 
                file = file_previous_model2 + '_HI_over_r_bins' + str(f'_MHI{log_MHI[jj]-bin_HI/2:0.2f}') + \
                       str(f'_{log_MHI[jj]+bin_HI/2:0.2f}') + str(f'_z{ThisRedshiftList[ii]:0.2f}')+'.csv'
                df = pd.read_csv(file)               
                subplot.plot(df['log10_HI_r'],df['Fraction'], color='red', linestyle=linestyle_previous_model2)
               
            

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
            label="%0.2f" % (log_MHI[jj]-bin_HI/2.) + r'$<\mathrm{log_{10}}(M_{HI}/$'+Msun+'$)<$' + "%0.2f" % (log_MHI[jj]+bin_HI/2.)
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.85, 
                        color='black', xlog=0, ylog=0, label=label, 
                        fontsize=13, fontweight='normal') 
            
            if(jj==0):
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.14, y_percentage=0.75, color='black', 
                            xlog=0, ylog=0, label=prefix_this_model, fontsize=13, fontweight='normal') 
                plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, y_percentage=0.78, color='red',
                            x2_percentage=0.12, xlog=0, ylog=0, linestyle='-', linewidth=2)
        
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.14, y_percentage=0.65, color='black', 
                            xlog=0, ylog=0, label=prefix_previous_model2, fontsize=12, fontweight='normal') 
                plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, y_percentage=0.68, color='red',
                            x2_percentage=0.12, xlog=0, ylog=0, linestyle=linestyle_previous_model2, linewidth=2)
            
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.12, y_percentage=0.55, color='black', 
                            xlog=0, ylog=0, label='Haynes2011',  fontsize=13, fontweight='normal') 
                plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.07, y_percentage=0.575, 
                    color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.02) 
        
            
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYF19_HI_over_Lr_vs_HI_bins.pdf')
    plt.close()

    return   
#end gas fraction


def HI_MF(ThisRedshiftList):

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
            
        xlab='$\mathrm{log_{10}}(M_{\mathrm{HI}}/$'+Msun+'$)$'   
        ylab='$\mathrm{log_{10}}(\phi/(\mathrm{Mpc^{-3}} \mathrm{dex}))$'     
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
                    
        #WRITE OUTPUT      
        if(write_to_file==1):
            df = pd.DataFrame({'log10_HI':x_axis, 'log10_phi':y_axis})
            file = Datadir+file_to_write+'HI_MF'+str(f'_z{ThisRedshiftList[ii]:0.2f}')+'.csv'
            df.to_csv(file,index=False)
            #df = pd.read_csv(file)
            #subplot.plot(df['log10_HI'],df['log10_phi'], color='black')
        
        #OBSERVATIONS
        h=0.75
        file = Datadir + 'zwaan2005.txt'       
        obs = Table.read(file, format='ascii')      
        obs_x = obs['col1']-2.*np.log10(h) 
        obs_y = obs['col2']     
        obs_y_err = [-obs['col3'],obs['col4']]
       
        subplot.errorbar(obs_x, obs_y, yerr=obs_y_err,fmt='o', markersize=5, ecolor='blue', color='blue',zorder=+3,capsize=2)
        
        
        file = Datadir + 'haynes2011_gmf.txt'       
        obs = Table.read(file, format='ascii')      
        obs_x = obs['col1']
        obs_y = obs['col2']
        obs_y_err = [-obs['col3'],obs['col4']]       
        subplot.errorbar(obs_x, obs_y, yerr=obs_y_err,fmt='o', markersize=5, ecolor='purple',
                         color='purple',zorder=+3,capsize=2)
        #subplot.errorbar(obs_x, obs_y, yerr=obs_y_err,fmt='o', markersize=5, ecolor='green', color='green',zorder=+3)
         
        #JONES 2018    
        df = pd.read_csv(Datadir+'Jones2018_HIMF_err.csv')      
        subplot.errorbar(df['x'], df['y'], yerr=[-df['err_down'],df['err_up']],fmt='o', markersize=5, 
                         ecolor='limegreen', color='limegreen',capsize=2)
           
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
                    color='black', xlog=0, ylog=0, label='Zwaan2005', 
                    fontsize=12, fontweight='normal') 
        plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.52, 
                    color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15) 
        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.43, 
                    color='black', xlog=0, ylog=0, label='Haynes2011', 
                    fontsize=12, fontweight='normal') 
        plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.45, 
                    color='purple', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15) 
        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.36, 
                    color='black', xlog=0, ylog=0, label='Jones2018', 
                    fontsize=12, fontweight='normal') 
        plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.38, 
                    color='limegreen', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15) 
        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.12, y_percentage=0.2, 
                    color='black', xlog=0, ylog=0, label=prefix_previous_model2, fontsize=10, fontweight='normal') 
        plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, y_percentage=0.22, 
                    color='red', x2_percentage=0.1, xlog=0, ylog=0, linestyle=linestyle_previous_model2, linewidth=2)
                    
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.12, y_percentage=0.13, 
                    color='black', xlog=0, ylog=0, label=prefix_this_model, fontsize=10, fontweight='normal') 
        plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, y_percentage=0.15, 
                    color='red', x2_percentage=0.1, xlog=0, ylog=0, linestyle='-', linewidth=2) 
        
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYF19_HI_MF.pdf')
    plt.close()
   
    return 
    '''xlim=[7.0,11.5]
    ylim=[4.0,12.0]
    bin=0.25
    
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
            
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.25))
            
    xlab='$\mathrm{log_{10}}(\mathrm{M_*}/$'+Msun+'$)$'   
    ylab='$\mathrm{log_{10}}(\mathrm{M_{\mathrm{HI}}}/$'+Msun+'$)$'      
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   
          
    #MODEL
    #MR
    (sel)=select_current_redshift(G_MRII, ThisRedshiftList, ii, FullSnapshotList_MRII)         
    G0_MR=G_MRII[sel]   
        
    G0_MR=G0_MR[(G0_MR['ColdGas']>0.) & (G0_MR['H2fraction']<1.)]
    HI=(np.log10(G0_MR['ColdGas']*(1.-G0_MR['H2fraction'])*1.e10/Hubble_h))  
    H2=(np.log10(G0_MR['ColdGas']*(G0_MR['H2fraction'])*1.e10/Hubble_h))  
  
        
    subplot.scatter(np.log10(G0_MR['StellarMass']*1e10/Hubble_h),HI,color='blue',marker='o',s=5)
    subplot.scatter(np.log10(G0_MR['StellarMass']*1e10/Hubble_h),H2,color='red',marker='o',s=5)
       

    
    pdf.savefig()
    plt.close()'''
   
#end HI_MF






def H2D(ThisRedshiftList):

    
        
    xlim=[-0.25,5.0]
    ylim=[6.7,8.5]
    bin=0.25

    plot_color=['red','purple']        

    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    

    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.25))
    subplot.yaxis.set_minor_locator(MultipleLocator(0.1))

    xlab='redshift'   
    ylab='$\\rho_{\mathrm{H_2}}/($'+Msun+'$ \mathrm{Mpc}^{-3})$'     
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   
    
    H2D_MR = np.zeros(len(ThisRedshiftList),dtype=np.float32)
    H2D_MRII = np.zeros(len(ThisRedshiftList),dtype=np.float32)
    for ii in range(0,len(ThisRedshiftList)):             
        #MODEL
        #MR
        '''(sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)         
        G0_MR=G_MR[sel]
        G0_MR = G0_MR[G0_MR['H2fraction']>0.]
        H2D_MR[ii] = np.sum(G0_MR['ColdGas']*1.e10/Hubble_h*G0_MR['H2fraction'])/Volume_MR'''
        
        if(MRII==1):
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, ii, FullSnapshotList_MRII)         
            G0_MRII=G_MRII[sel]  
            G0_MRII = G0_MRII[G0_MRII['H2fraction']>0.]
            H2D_MRII[ii] = np.sum(G0_MRII['ColdGas']*1.e10/Hubble_h*G0_MRII['H2fraction'])/Volume_MRII
        #else:
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)         
        G0_MR=G_MR[sel]  
        G0_MR = G0_MR[G0_MR['H2fraction']>0.]
        H2D_MR[ii] = np.sum(G0_MR['ColdGas']*1.e10/Hubble_h*G0_MR['H2fraction'])/Volume_MR
  
    if(MRII==1):
        subplot.plot(ThisRedshiftList,np.log10(H2D_MRII),color='red',linestyle='-')
    #else:
    subplot.plot(ThisRedshiftList,np.log10(H2D_MR),color='red',linestyle='--')    
    
    if(write_to_file==1):
        df = pd.DataFrame({'redshift':ThisRedshiftList, 'log10_H2D':np.log10(H2D_MR)})
        file = Datadir+file_to_write+'H2D_MR'+'.csv'
        df.to_csv(file,index=False)
        
        df = pd.DataFrame({'redshift':ThisRedshiftList, 'log10_H2D':np.log10(H2D_MRII)})
        file = Datadir+file_to_write+'H2D_MRII'+'.csv'
        df.to_csv(file,index=False)
   

    #OBS Decarli2019
    x_down = np.array([1.0,2.0,3.0])
    x_up = np.array([1.7,3.1,4.475])
    y_down = np.array([7.63,7.26,6.97])
    y_up = np.array([8.05,8.10,7.77])
    
    x = (x_down+x_up)/2.
    y = (y_down+y_up)/2.
    x_err = (x_up-x_down)/2.
    y_err = (y_up-y_down)/2.  
    subplot.errorbar(x, y, xerr=x_err,yerr=y_err, fmt='o', markersize=5, ecolor='blue', color='blue',zorder=+3) 
    
    #Saintonge
    #x = 0.0
    #y = np.log10(1.1e7)
    #y_err=[[0.23],[0.26]]
    #subplot.errorbar(x,y, yerr=y_err, fmt='o', markersize=5, ecolor='purple', color='purple',zorder=+3) 
   
    #Keres2003
    x = 0.0
    y = np.log10(3.1e7*Hubble_h)
    y_err_up = np.log10(3.1e7*Hubble_h+1.2e7*Hubble_h) - y
    y_err_down = y - np.log10(3.1e7*Hubble_h-1.2e7*Hubble_h)    
    y_err=[[y_err_down, y_err_up]]  
    subplot.errorbar(x,y, yerr=y_err, fmt='o', markersize=5, ecolor='purple', color='purple',zorder=+3) 

    #LAbels
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.12, y_percentage=0.92, color='black', 
                xlog=0, ylog=0, label=prefix_this_model+' MSII', fontsize=12, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, y_percentage=0.94, color='red',
                x2_percentage=0.1, xlog=0, ylog=0, linestyle='-', linewidth=2)
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.12, y_percentage=0.86, color='black', 
                xlog=0, ylog=0, label=prefix_this_model+' MS', fontsize=12, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, y_percentage=0.88, color='red',
                x2_percentage=0.1, xlog=0, ylog=0, linestyle='--', linewidth=2)
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.54, y_percentage=0.92, 
                color='black', xlog=0, ylog=0, label='Decarli2019',  fontsize=12, fontweight='normal') 
    plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.52, y_percentage=0.94, 
                color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.04) 
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.54, y_percentage=0.86, 
                color='black', xlog=0, ylog=0, label='Keres2003',  fontsize=12, fontweight='normal') 
    plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.52, y_percentage=0.88, 
                color='purple', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.04) 
    
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYF19_H2D.pdf')
    plt.close()
    
    return 

#end H2D











def coldgas_vs_stellarmass(ThisRedshiftList):
  
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
        SSFR = np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h))   
        G0_MR=G0_MR[(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)>10.) & (G0_MR['ColdGas']>0.) & 
                    (SSFR>np.log10(2.*(1+ThisRedshiftList[ii])**2./(1.37e10))-1.0)] 
                   
        #G0_MR=G0_MR[(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)>10.) & (G0_MR['ColdGas']>0.)]
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
        
        xlim=[9.5,11.5]
        ylim=[-2.0,0.5]
        bin=0.25
           
        subplot=plt.subplot()    
        subplot.set_ylim(ylim), subplot.set_xlim(xlim) 
        ylab='$\log_{10}(M_{\mathrm{gas}}/M_*)$' 
        subplot.set_ylabel(ylab, fontsize=14) 
        xlab='$\log_{10}(M_*/$'+Msun+'$)$'           
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
                         fmt='o', markersize=5, ecolor='blue', color='blue',zorder=+3)    
               
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
    #pdf.savefig()
    #plt.close()
   
    
    
    
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
        xlab='$\log_{10}(M_*/$'+Msun+'$)$'           
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
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()

    return 
#end

















def main_sequence(ThisRedshiftList):
       
    xlim=[8.5,11.5]
    #xlim=[8.5,15.5]
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
        
        xlab='$\mathrm{log_{10}}(M_*/$'+Msun+'$)$'      
        if ii==0:
            ylab='$\mathrm{log_{10}}(\mathrm{SFR}/($'+Msun+'$yr^{-1}))$'
        else:
            ylab=''      
        subplot.set_xlabel(xlab, fontsize=15), subplot.set_ylabel(ylab, fontsize=15)
        
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
        #log_StellarMass=np.log10(G0_MR['Mvir']*1e10/Hubble_h)
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
            
        #(slope,b)=get_slope(10.5,-2.0,12.,0.5)
        #print(slope,b)
        #x_arr=np.arange(xlim[0],xlim[1],0.1)
        #subplot.plot(x_arr,x_arr*slope+b, color='red', linewidth=2, linestyle='-')            
            
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
                    label='Elbaz2007', fontsize=10, fontweight='normal') 
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
                        color='black', xlog=0, ylog=0, label='Karim2011 - 0.2<z<0.4', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.95, 
                        color='red', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.1,mfc='white')
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.88, 
                        color='black', xlog=0, ylog=0, label='Karim2011 - 0.4<z<0.6', 
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
                        color='black', xlog=0, ylog=0, label='Karim2011 - 0.8<z<1.0', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.95, 
                        color='red', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.1,mfc='white')
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.88, 
                        color='black', xlog=0, ylog=0, label='Karim2011 - 1.0<z<1.2', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.90, 
                        color='red', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.1)
        
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.83, 
                        color='black', xlog=0, ylog=0, label='Whitaker2013 - 0.5<z<1.0', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.85, 
                        color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.1,mfc='white')
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.78, 
                        color='black', xlog=0, ylog=0, label='Whitaker2013 - 1.0<z<1.5', 
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
                        color='black', xlog=0, ylog=0, label='Karim2011 - 1.6<z<2.0', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.95, 
                        color='red', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.1,mfc='white')
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.88, 
                        color='black', xlog=0, ylog=0, label='Karim2011 - 2.0<z<2.5', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.90, 
                        color='red', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.1)
        
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.83, 
                        color='black', xlog=0, ylog=0, label='Whitaker2013 - 1.5<z<2.0', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.85, 
                        color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.1,mfc='white')
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.78, 
                        color='black', xlog=0, ylog=0, label='Whitaker2013 - 2.0<z<2.5', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.79, 
                        color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.1)
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.73, 
                        color='black', xlog=0, ylog=0, label='Shivaei2016 - 1.4<z<2.6', 
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
                        color='black', xlog=0, ylog=0, label='Karim2011 - 2.5<z<3.0', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.95, 
                        color='red', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.1)
    #endfor
        
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYF19_main_sequence.pdf')
    plt.close()
    
    return 
#endif stellar_mass_vs_sfr












def SSFR_mass(ThisRedshiftList):
    
    xlim=[8.5,11.5]   
    ylim=[-3.5, 0.5]   
    bin=[0.25,0.25]
    #bin=[0.25,0.5]
    Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]

    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
   
    for ii in range(0,len(ThisRedshiftList)):
        
        char_redshift="%0.1f" % ThisRedshiftList[ii]
        
        subplot=plt.subplot()
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        
        xlab='$\mathrm{log_{10}}(M_*/$'+Msun+'$)$'        
        ylab='$\mathrm{log_{10}}(\mathrm{SSFR}/($'+Msun+'$yr^{-1}))$'
         
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)
        
        subplot.text(xlim[0]+2.,ylim[0]+.3,'z='+char_redshift, fontsize=14, fontweight='normal')
        
    
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
        subplot.yaxis.set_minor_locator(MultipleLocator(0.25))
     
        
        #MODEL
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)
        
        G0_MR=G_MR[sel]   
        G0_MR=G0_MR[(G0_MR['StellarMass']>0.)]
        log_StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii]) 
        log_StellarMass=np.log10(G0_MR['StellarMass']*1e10/Hubble_h)
        #log_StellarMass=np.log10(G0_MR['Mvir']*1e10/Hubble_h)
        SFR=G0_MR['Sfr']       
        log_SFR=np.log10(SFR)
       
        width=0.25  
        (slope,b)=get_slope(9.0,-2,11.0,-1.)
        sel=SFR < 0.001      
        log_SFR[sel]=(np.random.randn(len(log_SFR[sel]))*width*(slope*log_StellarMass[sel]+b) +
                      (slope*log_StellarMass[sel]+b))
      
        
        Ngals=len(G0_MR)
      
        H, xedges, yedges = np.histogram2d(log_StellarMass, np.log10(1.e9*10**log_SFR/10**log_StellarMass), bins=Nbins, 
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
            log_obs_mstar = dr7_gal_final['jarle_median_mass'] 
            log_obs_ssfr   = np.log10(1.e9*10**dr7_gal_final['median_sfr']/10**dr7_gal_final['jarle_median_mass'])
  
            mag=dr7_gal_final['jarle_mag'][:,2] #absolute mag
            zz=dr7_gal_final['z'] 
            max_d=10**((17.6-mag)/5.+1.)/1.e6 #apparent
            weight = 1./max_d**3

      
            frac=np.zeros([int(Nbinsx),int(Nbinsy)],dtype=np.float32) 
                      
            xind=0
            for jj in range (0,int(Nbinsx)):
                yind=0  
                xid=xlim[0]+bin[0]*jj                 
                for kk in range (0,int(Nbinsy)):
                    yid=ylim[0]+bin[1]*kk
                    sel=((log_obs_mstar > xid) & (log_obs_mstar < xid+bin[0]) &
                         (log_obs_ssfr > yid) & (log_obs_ssfr < yid+bin[1]))                      
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
            
          
            plot_label (subplot, 'label', xlim, ylim, 
                    x_percentage=0.15, y_percentage=0.85, color='black', xlog=0, ylog=0, 
                    label='SDSS-DR7', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim,
                    x_percentage=0.04, y_percentage=0.87, color='black', x2_percentage=0.13, 
                    xlog=0, ylog=0, linestyle='-', linewidth=2)
            
           
          
        
    #endfor
        
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()
    
    return   
#endif opt_SSFR_mass


















def ur_vs_r(ThisRedshiftList):
        
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
        #(sel)=select_current_redshift(G_MRII, ThisRedshiftList, ii, FullSnapshotList_MRII)        
        #G0_MR=G_MRII[sel]
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
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYF19_'+current_function+'.pdf')
    plt.close()
    
    return     
#endif ur_vs_r



def SSFR_evo(ThisRedshiftList):
    
    xlim=[9.5,11.5]
    ylim=[-13.0, -7.5]   
    bin=[0.25,0.25]
    Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]

    fig = plt.figure(figsize=(12,4))
    grid = gridspec.GridSpec(1, 4)
    grid.update(left=0.075, right = 0.93, wspace=0.0, hspace=0.0)

    for ii in range(0,len(ThisRedshiftList)):
        
        char_redshift="%0.1f" % ThisRedshiftList[ii]
        
        subplot=plt.subplot(grid[ii])
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        if ii==0:
            ylabel='$\log_{10}(\mathrm{sSFR}/\mathrm{yr}^{-1})$'
        else:
             ylabel=''
        xlab='$\mathrm{log_{10}}(M_*/$'+Msun+'$)$'        
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylabel, fontsize=14)
        
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
        subplot.yaxis.set_minor_locator(MultipleLocator(0.25)) 
        
                  
        #MODEL
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)
        
        G0_MR=G_MR[sel] 
        #G0_MR=G0_MR[G0_MR['Type']==0]
        '''NGals=len(G0_MR)
        Nrandom = 10000
        if(Nrandom>NGals):
            Nrandom=NGals
        G0_MR=np.random.choice(G0_MR, size=Nrandom) 
        log_StellarMass = np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)     
        SFR=G0_MR['Sfr']       
        log_SFR = np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h))       
        #subplot.scatter(Mass, SSFR, marker='o',s=1) 
        
        
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
            plt.colorbar(format='%0.1f') '''
        
        
                   
        log_StellarMass=np.log10(G0_MR['StellarMass']*1e10/Hubble_h)       
        SFR=G0_MR['Sfr']       
        log_SFR=np.log10(SFR)
       
        width=0.25  
        (slope,b)=get_slope(9.0,-2,11.0,-1.)
        sel=SFR < 0.001      
        log_SFR[sel]=(np.random.randn(len(log_SFR[sel]))*width*(slope*log_StellarMass[sel]+b) +
                      (slope*log_StellarMass[sel]+b))
      
        
        Ngals=len(G0_MR)
      
        H, xedges, yedges = np.histogram2d(log_StellarMass, np.log10(10**log_SFR/10**log_StellarMass), bins=Nbins, 
                                           range=[[xlim[0]-bin[0]/2., xlim[1]+bin[0]/2.], 
                                                  [ylim[0]-bin[1]/2., ylim[1]+bin[1]/2.]])            
        extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
        plt.subplots_adjust(bottom=0.15, left=0.15)        
        mylevels = np.log10(np.logspace(-2.0, 0.0, num=10))  
        H = zoom(H, 20)        
        H=np.log10(H/np.amax(H))           
        cont=plt.contourf(H.transpose()[::], origin='lower', cmap='Greys_r', levels=mylevels, extent=extent)   
        
      
        if ii==0:         
            cbar_ax = fig.add_axes([0.94, 0.15, 0.01, 0.73])            
            plt.colorbar(cax=cbar_ax, format='%0.1f')          
            cbar_ax.xaxis.set_major_locator(MultipleLocator(100))    
            
        #BestFit Cut 
        xx = np.arange(xlim[0],xlim[1]+0.5,0.01)
        yy = xx/xx-1.0 + np.log10(2.*(1+ThisRedshiftList[ii])**2/(1.37e10)) 
        subplot.plot(xx,yy,color='red', linestyle='-', linewidth=2)  
        
        xx = np.arange(xlim[0],xlim[1]+0.5,0.01)
        yy = xx/xx-1.0 + np.log10(2.*(1+ThisRedshiftList[ii])**2/(1.37e10)) +0.3
        subplot.plot(xx,yy,color='red', linestyle=':', linewidth=2) 
        
        xx = np.arange(xlim[0],xlim[1]+0.5,0.01)
        yy = xx/xx-1.0 + np.log10(2.*(1+ThisRedshiftList[ii])**2/(1.37e10)) -0.3
        subplot.plot(xx,yy,color='red', linestyle=':', linewidth=2)
        
        xx = np.arange(xlim[0],xlim[1]+0.5,0.01)
        yy = xx/xx-1.0 + np.log10(2.*(1+ThisRedshiftList[ii])**2/(1.37e10)) -1.0
        subplot.plot(xx,yy,color='red', linestyle='--', linewidth=2)
        
         
        
        #LABELS
        if ii>0:
            plt.tick_params(axis='y', which='both', left='on', labelleft='off')
        
        subplot.text(xlim[1]-0.75,ylim[1]-.5,'z='+char_redshift, fontsize=12, fontweight='normal')
        
        if ii==len(ThisRedshiftList)-1:
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.155, y_percentage=0.23, 
                        color='black', xlog=0, ylog=0, label='MS', fontsize=13, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim,x_percentage=0.04, y_percentage=0.25, 
                        color='red', x2_percentage=0.135, xlog=0, ylog=0, linestyle='-', linewidth=2)
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.155, y_percentage=0.15, 
                        color='black', xlog=0, ylog=0, label='MS $\pm$ 0.3 dex', fontsize=13, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim,x_percentage=0.04, y_percentage=0.17, 
                        color='red', x2_percentage=0.135, xlog=0, ylog=0, linestyle=':', linewidth=2)
            
        
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.155, y_percentage=0.07, 
                        color='black', xlog=0, ylog=0,  label='Passive cut', fontsize=13, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, y_percentage=0.09, 
                        color='red', x2_percentage=0.135,  xlog=0, ylog=0, linestyle='--', linewidth=2)
    #endfor
        
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYF19_SSFR_evo.pdf')
    plt.close()
    
    return   
#endif SSFR_evo


def UVJ_colour(ThisRedshiftList):
    
    xlim=[0.0,1.5]
    ylim=[0.0, 2.]   
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
        NGals=len(G0_MR)
      
        '''H, xedges, yedges = np.histogram2d(color_VJ, color_UV, bins=Nbins)            
        extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
        plt.subplots_adjust(bottom=0.15, left=0.15)        
        mylevels = np.linspace(1., Nbins[0], Nbins[0])*NGals/(Nbins[0]**2/2.5)        
        H = zoom(H, 20)        
        cont=plt.contourf(H.transpose()[::], origin='lower', cmap='Greys_r', levels=mylevels, extent=extent)        
        
        if ii==len(ThisRedshiftList)-1:
            plt.colorbar(format='%d')''' 
        
        
        H, xedges, yedges = np.histogram2d(color_VJ, color_UV, bins=Nbins,
                                           range=[[xlim[0]-bin[0]/2., xlim[1]+bin[0]/2.], 
                                                  [ylim[0]-bin[1]/2., ylim[1]+bin[1]/2.]])            
        extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
        plt.subplots_adjust(bottom=0.15, left=0.15)     
        mylevels = np.log10(np.logspace(-1.0, 0.0, num=10))  
        H = zoom(H, 20)        
        H=np.log10(H/np.amax(H))       
        cont=plt.contourf(H.transpose()[::], origin='lower', cmap='Greys_r', levels=mylevels, extent=extent)  
        #cont=plt.contourf(H.transpose()[::], origin='lower', cmap='rainbow', levels=mylevels, extent=extent)
        if ii==len(ThisRedshiftList)-1:
            plt.colorbar(format='%0.1f')
        
        
        Nrandom = 10000
        if(Nrandom>NGals):
            Nrandom=NGals
        G0_MR=np.random.choice(G0_MR, size=Nrandom) 
        Mass = np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)
        SSFR = np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h))  
        color_UV=G0_MR['MagDust'][:,0]-G0_MR['MagDust'][:,2]  
        color_VJ=G0_MR['MagDust'][:,2]-G0_MR['MagDust'][:,7]        
        
        sel = SSFR < np.log10((1+ThisRedshiftList[ii])/(2.*1.37e10))       
        #subplot.scatter(color_VJ, color_UV, marker='o',s=1, color='blue') 
        subplot.scatter(color_VJ[sel], color_UV[sel], marker='o',s=1, color='red') 
        
        
        
        
        #BestFit Cut         
        slope=slope_color_cut[ii+1]
        offset=offset_color_cut[ii+1]
        minimum_y=minimum_y_color_cut[ii+1]
        
        x_arr=np.arange(xlim[0],xlim[1]+0.5,0.01)
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
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()
    
    return   
#endif UVJ_colour




def UVJ_grid(ThisRedshiftList):
         
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
            label=mass_low+"$<log_{10}(M_{\star} /$'+Msun+'$)<$"+mass_high
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
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()
    
    return     
#endif UVJ_grid
    


def morphology_vs_stellarmass(ThisRedshiftList):
    
    for i_z in range(0,len(ThisRedshiftList)):        
        
        char_redshift="%0.2f" % ThisRedshiftList[i_z]
        
        xlim=[9.0,12.]
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
            
        xlab='$\mathrm{log_{10}}(\mathrm{M_{\star}}/$'+Msun+'$)$'   
        ylab='Fraction'     
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   
       
        #MODEL
        #MR
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)         
        G0_MR=G_MR[sel]                   
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[i_z])
        BulgeMassRatio=G0_MR['BulgeMass']/G0_MR['StellarMass']
        
        Mass_arr=np.arange(xlim[0],np.amax(StellarMass)+bin/2.,bin)
        BulgeFraction=np.zeros(len(Mass_arr),dtype=np.float32)
        CompFraction=np.zeros(len(Mass_arr),dtype=np.float32)
        DiskFraction=np.zeros(len(Mass_arr),dtype=np.float32)
                  
        for ll in range(0,len(Mass_arr)):
                sel_bulge=G0_MR[(BulgeMassRatio>0.7) &
                                (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)]
                sel_comp=G0_MR[(BulgeMassRatio<0.7) & (BulgeMassRatio>0.3) &
                               (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)]
                sel_disk=G0_MR[(BulgeMassRatio<0.3) & 
                                (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)]
                '''sel_bulge=G0_MR[(BulgeMassRatio>0.5) &
                                (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)]    
                sel_comp=G0_MR[(BulgeMassRatio<0.7) & (BulgeMassRatio>0.3) &
                               (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)]
                sel_disk=G0_MR[(BulgeMassRatio<0.5) & 
                                (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)]
                #print(len(sel_bulge),len(sel_disk),len(sel_irr))'''
                if(float(len(sel_bulge))+float(len(sel_disk))+float(len(sel_comp))>0):                     
                    BulgeFraction[ll]=float(len(sel_bulge))/(float(len(sel_bulge))+float(len(sel_disk))+float(len(sel_comp)))
                    CompFraction[ll]=float(len(sel_comp))/(float(len(sel_bulge))+float(len(sel_disk))+float(len(sel_comp)))
                    DiskFraction[ll]=float(len(sel_disk))/(float(len(sel_bulge))+float(len(sel_disk))+float(len(sel_comp)))
                    
        subplot.plot(Mass_arr, BulgeFraction, color='red', linestyle='--', linewidth=2)
        subplot.plot(Mass_arr, CompFraction, color='green', linestyle='--', linewidth=2)
        subplot.plot(Mass_arr, DiskFraction, color='blue', linestyle='--', linewidth=2)
        
        #WRITE OUTPUT      
        if(write_to_file==1):
            df = pd.DataFrame({'log10_M':Mass_arr, 'BulgeFraction':BulgeFraction, 
                               'CompFraction':CompFraction, 'DiskFraction':DiskFraction})
            file = Datadir+file_to_write+'Morphology_vs_StellarMass_MR'+str(f'_z{ThisRedshiftList[i_z]:0.2f}') + '.csv'
            df.to_csv(file,index=False)
            #df = pd.read_csv(file)
            #subplot.plot(df['log10_M'],df['BulgeFraction'], color='black')
                
        #MRII
        if(MRII==1):
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, i_z, FullSnapshotList_MRII)         
            G0_MRII=G_MRII[sel]                   
            StellarMass=stellar_mass_with_err(G0_MRII, Hubble_h, ThisRedshiftList[i_z])
            BulgeMassRatio=G0_MRII['BulgeMass']/G0_MRII['StellarMass']
        
            Mass_arr=np.arange(xlim[0],np.amax(StellarMass)+bin/2.,bin)
            BulgeFraction=np.zeros(len(Mass_arr),dtype=np.float32)
            CompFraction=np.zeros(len(Mass_arr),dtype=np.float32)
            DiskFraction=np.zeros(len(Mass_arr),dtype=np.float32)
            
                  
            for ll in range(0,len(Mass_arr)):
                sel_bulge=G0_MRII[(BulgeMassRatio>0.7) &
                                (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)]
                sel_comp=G0_MRII[(BulgeMassRatio<0.7) & (BulgeMassRatio>0.3) &
                               (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)]
                sel_disk=G0_MRII[(BulgeMassRatio<0.3) & 
                                (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)]
                '''sel_bulge=G0_MRII[(BulgeMassRatio>0.5) &
                                (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)] 
                sel_comp=G0_MRII[(BulgeMassRatio<0.7) & (BulgeMassRatio>0.3) &
                               (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)]
                sel_disk=G0_MRII[(BulgeMassRatio<0.5) & 
                                (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)]'''
                #print(len(sel_bulge),len(sel_disk),len(sel_irr))
                if(float(len(sel_bulge))+float(len(sel_disk))+float(len(sel_comp))>0):                     
                    BulgeFraction[ll]=float(len(sel_bulge))/(float(len(sel_bulge))+float(len(sel_disk))+float(len(sel_comp)))
                    CompFraction[ll]=float(len(sel_comp))/(float(len(sel_bulge))+float(len(sel_disk))+float(len(sel_comp)))
                    DiskFraction[ll]=float(len(sel_disk))/(float(len(sel_bulge))+float(len(sel_disk))+float(len(sel_comp)))
                    
            subplot.plot(Mass_arr, BulgeFraction, color='red', linestyle='-', linewidth=2)
            subplot.plot(Mass_arr, CompFraction, color='green', linestyle='-', linewidth=2)
            subplot.plot(Mass_arr, DiskFraction, color='blue', linestyle='-', linewidth=2)
            
            
            if(write_to_file==1):
                df = pd.DataFrame({'log10_M':Mass_arr, 'BulgeFraction':BulgeFraction, 
                               'CompFraction':CompFraction, 'DiskFraction':DiskFraction})
                file = Datadir+file_to_write+'Morphology_vs_StellarMass_MRII'+str(f'_z{ThisRedshiftList[i_z]:0.2f}')+'.csv'
                df.to_csv(file,index=False)
                #df = pd.read_csv(file)
                #subplot.plot(df['log10_M'],df['BulgeFraction'], color='black')
            
        #OBSERVATIONS
        h=0.7
        file = Datadir + 'conselice2006_bulge_fract.txt'       
        obs = Table.read(file, format='ascii')       
        #subplot.errorbar(obs['col1'], obs['col2'], obs['col3'], fmt='o', markersize=5, ecolor='red', color='red',capsize=2)
        file = Datadir + 'conselice2006_disk_fract.txt'       
        obs = Table.read(file, format='ascii')       
        #subplot.errorbar(obs['col1'], obs['col2'], obs['col3'], fmt='o', markersize=5, ecolor='blue', color='blue',capsize=2)
        file = Datadir + 'conselice2006_irr_fract.txt'       
        obs = Table.read(file, format='ascii')       
        #subplot.errorbar(obs['col1'], obs['col2'], obs['col3'],fmt='o', markersize=5, ecolor='green',color='green',capsize=2)
      
        #MCMC sample
        if opt_plot_MCMC_sample==1:
            file = MCMCSampledir + 'mcmc_plus_obs_BulgeFraction_z'+char_redshift+'.txt'          
            if os.path.isfile(file):
                obs = Table.read(file, format='ascii')      
                subplot.plot(obs['col1']-2.*np.log10(Hubble_h),obs['col4'], color='black', linewidth=2)
                          
                if i_z==len(ThisRedshiftList)-1:
                    plot_label (subplot, 'label', xlim, ylim, 
                        x_percentage=0.55, y_percentage=0.85, color='black', xlog=0, ylog=0, 
                        label='MCMC sample', fontsize=13, fontweight='normal') 
                    plot_label (subplot, 'line', xlim, ylim,
                        x_percentage=0.44, y_percentage=0.87, color='black', x2_percentage=0.53, 
                        xlog=0, ylog=0, linestyle='-', linewidth=2)    
          
        
        #bluck2019_frac
        obs_df = pd.read_csv(Datadir+'bluck2019_frac.dat', delimiter =' ')       
        obs_df = obs_df.iloc[::2,:]
        subplot.errorbar(obs_df['log_stellar_mass'], obs_df['frac_sph'], obs_df['frac_sph_err'], fmt='o', 
                         markersize=5, ecolor='red',color='red',capsize=2)  
        subplot.errorbar(obs_df['log_stellar_mass'], obs_df['frac_comp'], obs_df['frac_comp_err'], fmt='o', 
                         markersize=5, ecolor='green',color='green',capsize=2) 
        subplot.errorbar(obs_df['log_stellar_mass'], obs_df['frac_disc'], obs_df['frac_disc_err'], fmt='o', 
                         markersize=5, ecolor='blue',color='blue',capsize=2)
         
        '''obs_df = pd.read_csv(Datadir+'bluck2019_frac_05.dat', delimiter =' ')        
        obs_df = obs_df.iloc[::2,:]
        subplot.errorbar(obs_df['log_stellar_mass'], obs_df['frac_BD'], obs_df['frac_BD_err'], fmt='o', 
                         markersize=5, ecolor='red',color='red',capsize=2)  
        subplot.errorbar(obs_df['log_stellar_mass'], obs_df['frac_DD'], obs_df['frac_DD_err'], fmt='o', 
                         markersize=5, ecolor='blue',color='blue',capsize=2)''' 
                    
        #LABELS  
        plot_label (subplot, 'label', xlim, ylim,  x_percentage=0.11, y_percentage=0.91, 
                    color='black', xlog=0, ylog=0, label=prefix_this_model+' - MS', fontsize=12, fontweight='normal')          
        plot_label (subplot,'line',xlim,ylim,x_percentage=0.03,y_percentage=0.93,
                    color='red',x2_percentage=0.09,xlog=0,ylog=0,linestyle='--',linewidth=2)

        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.11, y_percentage=0.85, 
                    color='black', xlog=0, ylog=0, label=prefix_this_model+' - MSII', fontsize=12, fontweight='normal')  
        plot_label (subplot,'line',xlim,ylim,x_percentage=0.03,y_percentage=0.87,
                    color='red',x2_percentage=0.09,xlog=0,ylog=0,linestyle='-',linewidth=2)
    
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.11, y_percentage=0.79, 
                    color='black', xlog=0, ylog=0, label='Bluck2019',  fontsize=12, fontweight='normal') 
        plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.08, y_percentage=0.81, 
                    color='red', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.025) 
     
      
    
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYF19_morphology.pdf')
    plt.close()
   
    return 
#end morphology_vs_stellarmass   
    

    
def morphology_contour(ThisRedshiftList):
    
    for i_z in range(0,len(ThisRedshiftList)):        
        
        char_redshift="%0.2f" % ThisRedshiftList[i_z]
        
        xlim=[9.5,12.]
        ylim=[0.0,1.]
        bin=[0.25,0.25]
        Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]
        
        fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
        subplot=plt.subplot()
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
            
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))
            
        xlab='$\mathrm{log_{10}}(\mathrm{M_{\star}}/$'+Msun+'$)$'   
        ylab='B/T'     
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   
       
        #MRII
        if(MRII==1):
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, i_z, FullSnapshotList_MRII)         
            G0_MRII=G_MRII[sel]                   
            StellarMass=stellar_mass_with_err(G0_MRII, Hubble_h, ThisRedshiftList[i_z])
            BulgeMassRatio=G0_MRII['BulgeMass']/G0_MRII['StellarMass']
        
                       
            #plt.scatter(BulgeMass, BHMass, s=5, color='black')          
            H, xedges, yedges = np.histogram2d(StellarMass, BulgeMassRatio, bins=Nbins,
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
            
    
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    #plt.savefig('./fig/HYF19_morphology.pdf')
    plt.close()
   
    return 
#end morphology_contour

def morphology_SMF(ThisRedshiftList):
    
    B_T_ranges=np.array([0.0, 0.2, 0.5, 0.8, 1.0])
    
    
    obs_mstar0   = np.array([10.74, 10.998, 10.926,11.023])
    obs_alpha0   = np.array([-1.524, -1.270, -0.764, -0.580])
    obs_phistar0 = np.array([0.332, 0.489, 1.016, 1.203])
    for i_z in range(0,len(ThisRedshiftList)):        
        
        if(MRII==1):
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, i_z, FullSnapshotList_MRII)         
            G0_MRII=G_MRII[sel] 
        if(MR==1):    
            (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)         
            G0_MR=G_MR[sel]       
        #G0_MRII=G0_MR
        
        char_redshift="%0.2f" % ThisRedshiftList[i_z]
        
        xlim=[9.,12.]
        ylim=[-6.0,-1.5]
        bin=0.25
       
        
        fig = plt.figure(figsize=(two_two_size_large[0],two_two_size_large[1]))        
        grid = gridspec.GridSpec(2, 2)
        grid.update(wspace=0.0, hspace=0.0)      
       
        for jj in range(0,4):    
            subplot=plt.subplot(grid[jj])
            subplot.set_ylim(ylim), subplot.set_xlim(xlim)    

            #format axis
            majorFormatter = FormatStrFormatter('%d')
            subplot.xaxis.set_major_locator(MultipleLocator(1))    
            subplot.xaxis.set_minor_locator(MultipleLocator(0.25))

            if(jj==1 or jj==3):                 
                plt.tick_params(axis='y', which='both', left='on', labelleft='off')
            if(jj==0 or jj==1):                 
                plt.tick_params(axis='x', which='both', bottom='on', labelbottom='off')
            
            xlab='$\mathrm{log_{10}}(\mathrm{M_{\star}}/$'+Msun+'$)$'   
            if(jj==0 or jj==2):
                ylab='$\mathrm{log_{10}}(\phi /(\mathrm{Mpc^{-3}} \mathrm{dex}))$'     
            else:
                ylab=''
            subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   
 
               
            #MR     
            if(MR==1):           
                G0_MR=G0_MR[G0_MR['StellarMass']>0.]
                #Mass_MR=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[i_z])
                Mass_MR = np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)
                BulgeMassRatio=G0_MR['BulgeMass']/G0_MR['StellarMass']
                sel = (BulgeMassRatio>=B_T_ranges[jj]) & (BulgeMassRatio<=B_T_ranges[jj+1])
                #(x_axis,y_axis) = plot_mass_function(subplot, Mass_MR[sel], Volume_MR, xlim, bin, MRII, color='red')
                bin_arr=np.arange(xlim[0]-bin/2.0,xlim[1]+bin/2.0,bin)
                hist_MR=np.histogram(Mass_MR[sel], bins=bin_arr, range=(xlim[0],xlim[1]))     
                x_axis=hist_MR[1][0:len(hist_MR[1][:])-1]+bin/2.           
                hist_MR=hist_MR[0]       
                y_axis=np.log10(hist_MR/(Volume_MR*bin))
                subplot.plot(x_axis,y_axis, color='red', linewidth=2, linestyle='--') 
                
            #MRII
            if(MRII==1):           
                G0_MRII=G0_MRII[G0_MRII['StellarMass']>0.]
                Mass_MRII = np.log10(G0_MRII['StellarMass']*1.e10/Hubble_h)
                BulgeMassRatio=G0_MRII['BulgeMass']/G0_MRII['StellarMass']
                sel = (BulgeMassRatio>=B_T_ranges[jj]) & (BulgeMassRatio<=B_T_ranges[jj+1])
                bin_arr=np.arange(xlim[0]-bin/2.0,xlim[1]+bin/2.0,bin)
                hist_MR=np.histogram(Mass_MRII[sel], bins=bin_arr, range=(xlim[0],xlim[1]))     
                x_axis=hist_MR[1][0:len(hist_MR[1][:])-1]+bin/2.           
                hist_MR=hist_MR[0]       
                y_axis=np.log10(hist_MR/(Volume_MRII*bin))
                subplot.plot(x_axis,y_axis, color='red', linewidth=2, linestyle='-')

                
            #Observations
            MM=np.arange(8.0,12.0,0.1)
            Mstar = 10**obs_mstar0[jj]
            alpha1 = obs_alpha0[jj]+1.0
            phistar1 = obs_phistar0[jj]*1.e-3
            M=(10**MM)
            LLF1 = (phistar1*(M/Mstar)**(alpha1)*np.exp(-M/Mstar))*np.log(10)
            subplot.plot(np.log10(M),np.log10(LLF1), color='blue', linewidth=2, linestyle='-') 
            
            if(jj==0):
                label = str(f"{B_T_ranges[jj]:0.1f}<=B/T<{B_T_ranges[jj+1]:0.1f}")
            elif (jj==3):    
                 label = str(f"{B_T_ranges[jj]:0.1f}<B/T<={B_T_ranges[jj+1]:0.1f}")
            else:
                 label = str(f"{B_T_ranges[jj]:0.1f}<B/T<{B_T_ranges[jj+1]:0.1f}")
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.5, y_percentage=0.9, 
                        color='black', xlog=0, ylog=0, label=label,  fontsize=12, fontweight='normal')   
            
            if(jj==0):
                plot_label (subplot,'label',xlim,ylim,x_percentage=0.13,y_percentage=0.45, color='black',
                            xlog=0,ylog=0,label='MR',fontsize=12,fontweight='normal') 
                plot_label (subplot, 'line', xlim, ylim,x_percentage=0.04, y_percentage=0.47, color='red',
                            x2_percentage=0.12, xlog=0, ylog=0, linestyle='--', linewidth=2)
                plot_label (subplot,'label',xlim,ylim,x_percentage=0.13,y_percentage=0.38, color='black',
                            xlog=0,ylog=0,label='MRII',fontsize=12,fontweight='normal') 
                plot_label (subplot, 'line', xlim, ylim,x_percentage=0.04, y_percentage=0.40, color='red',
                            x2_percentage=0.12, xlog=0, ylog=0, linestyle='-', linewidth=2)
                plot_label (subplot,'label',xlim,ylim,x_percentage=0.13,y_percentage=0.31, color='black',
                            xlog=0,ylog=0,label='Thanjavur2016',fontsize=12,fontweight='normal') 
                plot_label (subplot, 'line', xlim, ylim,x_percentage=0.04, y_percentage=0.33, color='blue',
                            x2_percentage=0.12, xlog=0, ylog=0, linestyle='-', linewidth=2)
        
        #if(MRII==0):
        #    (x_axis,y_axis) = plot_mass_function(subplot, Mass_MR, Volume_MR, xlim, bin, MRII, color='red')
        #else:
        #    (x_axis,y_axis) = plot_mass_function(subplot, Mass_MR, Volume_MR, xlim, bin, MRII, 
        #                                         Mass_MRII=Mass_MRII, Volume_MRII=Volume_MRII, color='red')    
    
    
  
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYF19_morphology_SMF.pdf')
    plt.close()
   
    return 
#end morphology_SMF

def morphology_hists(ThisRedshiftList):
    
    '''for i_z in range(0,len(ThisRedshiftList)):        
        
        char_redshift="%0.2f" % ThisRedshiftList[i_z]
        
        xlim=[0.,1.]
        ylim=[0.0,1.]
        bin_hist=0.1
        bin=0.4
        
        mass = np.arange(9.2,11.6,bin)
        print(mass)
        plt_color=['purple','blue','green','yellow','orange','red','brown','black','black','black','black']
        cmap = plt.get_cmap('rainbow')  
        plt_color = [cmap(i) for i in np.linspace(0., 1., len(mass))]
        
        fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
        subplot=plt.subplot()
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
            
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(0.25))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.1))
            
        xlab='$B/T$'   
        ylab='fraction'     
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   
       
        #MRII
        if(MRII==1):
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, i_z, FullSnapshotList_MRII)         
            G0_MRII=G_MRII[sel] 
            
            StellarMass=stellar_mass_with_err(G0_MRII, Hubble_h, ThisRedshiftList[i_z])
            BulgeMassRatio=G0_MRII['BulgeMass']/G0_MRII['StellarMass']
            for idx in range(0, len(mass)):
                sel = (StellarMass>mass[idx]-bin/2.) & (StellarMass<mass[idx]+bin/2.)
            
                bin_arr=np.arange(xlim[0],xlim[1]+bin_hist,bin_hist)
                hist=np.histogram(BulgeMassRatio[sel], bins=bin_arr, range=(xlim[0],xlim[1]))            
                x_axis=hist[1][0:len(hist[1][:])-1]+bin_hist/2.           
                y_axis=hist[0]/np.sum(hist[0])     
                subplot.plot(x_axis,y_axis, color=plt_color[idx], linewidth=2, linestyle='-') 
            
        
            
    
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'_1.pdf')
    #plt.savefig('./fig/HYF19_morphology.pdf')'''
       


    for i_z in range(0,len(ThisRedshiftList)):        
        
        char_redshift="%0.2f" % ThisRedshiftList[i_z]
        
        xlim=[0.,1.]
        ylim=[-6.0,-1.]
        bin_hist=0.2
        bin=0.4
        
        mass = np.arange(9.2,11.6,bin)
        print(mass)
        plt_color=['purple','blue','green','yellow','orange','red','brown','black','black','black','black']
        cmap = plt.get_cmap('rainbow')  
        plt_color = [cmap(i) for i in np.linspace(0., 1., len(mass))]
        
        fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
        subplot=plt.subplot()
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
            
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(0.2))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.1))
            
        xlab='$B/T$'   
        ylab='$\mathrm{log_{10}}(\phi /(\mathrm{Mpc^{-3}} \mathrm{dex}))$'     
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   
           
        mass_switch = 100.4
        if(MR==1):
            (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)         
            G0_MR=G_MR[sel] 
            
            StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[i_z])
            BulgeMassRatio=G0_MR['BulgeMass']/G0_MR['StellarMass']
            for idx in range(0, len(mass)):
                if(mass[idx]<mass_switch):
                    sel = (StellarMass>mass[idx]-bin/2.) & (StellarMass<mass[idx]+bin/2.)

                    bin_arr=np.arange(xlim[0]-bin_hist/2.0,xlim[1]+bin_hist/2.0,bin_hist)
                    hist=np.histogram(BulgeMassRatio[sel], bins=bin_arr, range=(xlim[0],xlim[1]))            
                    x_axis=hist[1][0:len(hist[1][:])-1]+bin_hist/2.           
                    y_axis=hist[0]/(Volume_MR*bin_hist)
                    #y_axis=hist[0]
                    subplot.plot(x_axis,np.log10(y_axis), color=plt_color[idx], linewidth=2, linestyle='--') 
                
        #MRII
        if(MRII==1):
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, i_z, FullSnapshotList_MRII)         
            G0_MRII=G_MRII[sel] 
            
            StellarMass=stellar_mass_with_err(G0_MRII, Hubble_h, ThisRedshiftList[i_z])
            BulgeMassRatio=G0_MRII['BulgeMass']/G0_MRII['StellarMass']
            for idx in range(0, len(mass)):
                if(mass[idx]<mass_switch):
                    sel = (StellarMass>mass[idx]-bin/2.) & (StellarMass<mass[idx]+bin/2.)
                    bin_arr=np.arange(xlim[0]-bin_hist/2.0,xlim[1]+bin_hist/2.0,bin_hist)
                    hist=np.histogram(BulgeMassRatio[sel], bins=bin_arr, range=(xlim[0],xlim[1]))            
                    x_axis=hist[1][0:len(hist[1][:])-1]+bin_hist/2.           
                    y_axis=hist[0]/(Volume_MRII*bin_hist)
                    #y_axis=hist[0]

                    subplot.plot(x_axis,np.log10(y_axis), color=plt_color[idx], linewidth=2, linestyle='-') 
    
    xx = [0.05,0.15,0.25,0.37,0.49,0.61,0.73,0.85,1.0,1.0,1.0]    
    for idx in range(0, len(mass)):
        plot_label (subplot, 'label', xlim, ylim, x_percentage=xx[idx], y_percentage=0.9, color=plt_color[idx], 
                    xlog=0, ylog=0, label=str(f"{mass[idx]:0.1f}"),  fontsize=12, fontweight='normal')     
            
    
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    #plt.savefig('./fig/HYF19_morphology.pdf')
    plt.close()


    return 
#end morphology_hists       





def old_sizes_vs_stellarmass(ThisRedshiftList):
     
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
            xlab='$\mathrm{log_{10}}(M_*/$'+Msun+'$)$'       
            ylab='$R_{50}(\mathrm{Kpc})$'            
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
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYW17_plots_sizes_vs_stellarmass.pdf')
    plt.close() 

    return 
#end   sizes_vs_stellarmass 


def sizes_vs_stellarmass_shen(ThisRedshiftList):

    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()    
    
    plot_color=['blue','red']
    
    gal_types = ['discs', 'bulges']
       
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
            xlab='$\mathrm{log_{10}}(M_*/$'+Msun+'$)$'       
            ylab='$\mathrm{log_{10}}(R_{50}/\mathrm{kpc})$'            
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
            
            #subplot.plot(Shen['col1'],np.log10(Shen['col2']), color=plot_color[i_gal_type], linestyle=':', alpha=0.7)   
            
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
            HalfMassRadius=Gal['StellarHalfMassRadius']*1000./Hubble_h #(from Mpc/h -> Kpc)      
            #StellarDiskRadius=Gal['StellarDiskRadius']*1000./Hubble_h #(from Mpc/h -> Kpc)    

            (x_binned, median, mean, pc16, pc84, rms) = median_and_percentiles (bin, xlim[0], xlim[1], 
                                                                                StellarMass, np.log10(HalfMassRadius))      
            subplot.plot(x_binned, median,color=plot_color[i_gal_type], linewidth=2)           
            subplot.plot(x_binned, median+rms,color=plot_color[i_gal_type], linewidth=2, linestyle='--')
            subplot.plot(x_binned, median-rms,color=plot_color[i_gal_type], linewidth=2, linestyle='--')
    
    
            #WRITE OUTPUT      
            if(write_to_file==1):
                df = pd.DataFrame({'log10_M':x_binned, 'HalfMassRadius':median, 'pc16':median-rms, 'pc84':median+rms})
                file = Datadir+file_to_write+'Sizes_vs_StellarMass_'+gal_types[i_gal_type] + \
                       str(f'_z{ThisRedshiftList[i_z]:0.2f}') + '.csv'
                df.to_csv(file,index=False)
                #df = pd.read_csv(file)
                #subplot.plot(df['log10_M'],df['HalfMassRadius'], color='black')
    
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
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYF19_sizes_vs_stellarmass.pdf')
    plt.close() 
  
    return
#end   sizes_vs_stellarmass 



def sizes_vs_stellarmass_allz(ThisRedshiftList):

    model_to_print = 'Hen15'
    
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()    
    
    plot_color=['black','red', 'orange', 'green']
    plot_line_style=['-','--','-.',':']
    
    xlim=[9.,11.5]       
    ylim=[-0.1,0.9]   
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)            
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(0.5))    
    subplot.xaxis.set_minor_locator(MultipleLocator(.25))  
    subplot.yaxis.set_major_locator(MultipleLocator(0.5))    
    subplot.yaxis.set_minor_locator(MultipleLocator(0.1))  
    xlab='$\log_{10}(M_*/$'+Msun+'$)$'       
    ylab='$\log_{10}(\mathrm{R_{50}}/\mathrm{Kpc})$'            
    subplot.set_xlabel(xlab, fontsize=14),subplot.set_ylabel(ylab, fontsize=14)
    #subplot.set_yscale('log')
    
    obs_mass= np.array([9.25,9.75,10.25,10.75,11.25])
    obs_df = pd.read_csv(Datadir+'vanderwel2014.csv')
    obs_z025 = obs_df.iloc[0:3,:]
    obs_z125 = obs_df.iloc[6:9,:]
    obs_z225 = obs_df.iloc[12:15,:]
   
    #late_types = 
    #print(obs_z125.iloc[0,[11,13,15,17,19]])
    #print(obs_z125.iloc[1,[11,13,15,17,19]])
    #print(obs_z125.iloc[2,[11,13,15,17,19]])
   
    for i_z in range(0,len(ThisRedshiftList)):            
        
        if(MRII==1):
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, i_z, FullSnapshotList_MRII)        
            G0_MR=G_MRII[sel]   
        else:
            (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)        
            G0_MR=G_MR[sel]   
               
        char_redshift="%0.2f" % ThisRedshiftList[i_z] 
        
        #log_StellarMass = stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])    
        log_StellarMass=np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)
        log_SFR=np.log10(G0_MR['Sfr'])
        #BulgeSize=G0_MR['BulgeSize']*1000./Hubble_h #(from Mpc/h -> Kpc)        
        #StellarDiskRadius=G0_MR['StellarDiskRadius']*1000./Hubble_h #(from Mpc/h -> Kpc)         
        Size = G0_MR['StellarHalfMassRadius']*1000./Hubble_h #(from Mpc/h -> Kpc) 
        #Size = G0_MR['GasDiskRadius']/3.*1000./Hubble_h #(from Mpc/h -> Kpc) 
        
        if(i_z==0): 
            
            #observations
            
            #Mosleh
            file=Datadir+'mosleh2013.fits'
            hdul = fits.open(file)
            obs_data = hdul[1].data
            #print(hdul[1].header)
            #print(obs_data['re_np_kpc_corr'], obs_data['Mass_final'], obs_data['SSFR_MEDIAN'])
          
            #mosleh2013 passive
            bin=0.2
            mass_bins = np.arange(9.1,11.3,0.2)
            re_median = np.zeros(len(mass_bins), dtype=np.float32)
            re_pc16 = np.zeros(len(mass_bins), dtype=np.float32)
            re_pc84 = np.zeros(len(mass_bins), dtype=np.float32)
            for jj in range(0, len(mass_bins)):
                sel = ((obs_data['SSFR_MEDIAN']<-11.0) & (obs_data['Mass_final']>mass_bins[jj]-bin/2.0) &
                       (obs_data['Mass_final']<mass_bins[jj]+bin/2.0))                
                obs_re = np.log10(obs_data['re_np_kpc_corr'][sel])
                N_values = len(obs_re)
                N_bootstrap = 1000
                medians = np.zeros(N_bootstrap, dtype=np.float32)
                for ii in range(0, N_bootstrap):
                    bootstrap_sample = np.random.choice(obs_re, size=N_values, replace=True)
                    medians[ii] = np.median(bootstrap_sample)
                    
                re_median[jj] = np.median(medians)
                re_pc16[jj] = np.sort(medians)[int(16*len(medians)/100)]      
                re_pc84 [jj] = np.sort(medians)[int(84*len(medians)/100)]      
            
            err = (re_pc84-re_pc16)/2.
            subplot.errorbar(mass_bins,re_median, err,fmt='o',markersize=5, ecolor='red', color='red',capsize=2) 
            
            #mosleh2013 active
            mass_bins = np.arange(9.1,11.1,0.2)
            re_median = np.zeros(len(mass_bins), dtype=np.float32)
            re_pc16 = np.zeros(len(mass_bins), dtype=np.float32)
            re_pc84 = np.zeros(len(mass_bins), dtype=np.float32)
            for jj in range(0, len(mass_bins)):
                sel = ((obs_data['SSFR_MEDIAN']>-11.0) & (obs_data['Mass_final']>mass_bins[jj]-bin/2.0) &
                       (obs_data['Mass_final']<mass_bins[jj]+bin/2.0))                
                obs_re = np.log10(obs_data['re_np_kpc_corr'][sel])
                N_values = len(obs_re)
                N_bootstrap = 1000
                medians = np.zeros(N_bootstrap, dtype=np.float32)
                for ii in range(0, N_bootstrap):
                    bootstrap_sample = np.random.choice(obs_re, size=N_values, replace=True)
                    medians[ii] = np.median(bootstrap_sample)
                    
                re_median[jj] = np.median(medians)
                re_pc16[jj] = np.sort(medians)[int(16*len(medians)/100)]      
                re_pc84 [jj] = np.sort(medians)[int(84*len(medians)/100)]      
            
            err = (re_pc84-re_pc16)/2.            
            subplot.errorbar(mass_bins,re_median, err,fmt='o',markersize=5, ecolor='blue', color='blue',capsize=2) 
            
            #vanderwel
            #subplot.plot(obs_mass, obs_z025.iloc[1,[1,3,5,7,9]]) 
            #subplot.plot(obs_mass, obs_z025.iloc[1,[11,13,15,17,19]], color='black')  
            #subplot.plot(obs_mass, obs_z125.iloc[1,[11,13,15,17,19]], color='black')  
            #subplot.plot(obs_mass, obs_z225.iloc[1,[11,13,15,17,19]], color='black')
            #subplot.text(9.7,0.65,'z=0.25',fontsize=10)
            #subplot.text(9.7,0.55,'z=1.25',fontsize=10)
            #subplot.text(9.7,0.4,'z=2.25',fontsize=10)
            #subplot.text(9.1,0.72,'vanderwel2014',fontsize=10)
            
            sel = ((log_SFR-log_StellarMass)<np.log10((1+ThisRedshiftList[i_z])**2/(1.37e10/2.)) -1.0)
            bin=0.25
            (x_binned, median, mean, pc16, pc84, rms) = median_and_percentiles (bin, xlim[0], xlim[1], log_StellarMass[sel],
                                                                                np.log10(Size[sel])) 
            subplot.plot(x_binned, median,color='red', linewidth=2,linestyle=plot_line_style[i_z])           
            if(write_to_file==1):
                df = pd.DataFrame({'log10_M':x_binned, 'log10_R':median})
                file = Datadir+file_to_write+'Sizes_vs_StellarMass_Passive'+str(f'_z{ThisRedshiftList[i_z]:0.2f}')+'.csv'
                df.to_csv(file,index=False)        
                
                
        sel = ((log_SFR-log_StellarMass)>np.log10((1+ThisRedshiftList[i_z])**2/(1.37e10/2.)) -1.0)
        bin=0.25
        (x_binned, median, mean, pc16, pc84, rms) = median_and_percentiles (bin, xlim[0], xlim[1], log_StellarMass[sel],
                                                                            np.log10(Size[sel]))   
        sel = (x_binned<11.0) |  ((median>0.) & (x_binned>11.0))
        subplot.plot(x_binned[sel], median[sel],color='blue', linewidth=2,linestyle=plot_line_style[i_z])           
        if(write_to_file==1):
            df = pd.DataFrame({'log10_M':x_binned, 'log10_R':median})
            file = Datadir+file_to_write+'Sizes_vs_StellarMass_Active'+str(f'_z{ThisRedshiftList[i_z]:0.2f}')+'.csv'
            df.to_csv(file,index=False)      
        
        #WRITE OUTPUT        
        #file = Datadir+"size_vs_stellarmass_disks_"+model_to_print+"_z"+char_redshift+".txt"
        #write_to_file(x_binned, median, file) 
    
    
    #labels
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.1, y_percentage=0.9, 
                color='black', xlog=0, ylog=0, label='Mosleh2013 (z=0)', 
                fontsize=10, fontweight='normal') 
    plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.04, y_percentage=0.915, 
                color='red', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.035)
    plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.07, y_percentage=0.915, 
                color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.035)
   
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.04, y_percentage=0.79, color='black', 
                xlog=0, ylog=0, label=prefix_this_model, fontsize=10,  fontweight='normal')             
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.1, y_percentage=0.74, color='black', 
                xlog=0, ylog=0, label='z=0', fontsize=10,  fontweight='normal')             
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, y_percentage=0.755,  color='blue',
                x2_percentage=0.09, xlog=0, ylog=0, linestyle='-', linewidth=1.5)
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.1, y_percentage=0.69, color='black', 
                xlog=0, ylog=0, label='z=1', fontsize=10,  fontweight='normal')             
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, y_percentage=0.705,  color='blue',
                x2_percentage=0.09, xlog=0, ylog=0, linestyle='--', linewidth=1.5)
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.1, y_percentage=0.64, color='black', 
                xlog=0, ylog=0, label='z=2', fontsize=10,  fontweight='normal')             
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, y_percentage=0.655,  color='blue',
                x2_percentage=0.09, xlog=0, ylog=0, linestyle='-.', linewidth=1.5)
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.1, y_percentage=0.59, color='black', 
                xlog=0, ylog=0, label='z=3', fontsize=10,  fontweight='normal')             
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, y_percentage=0.605,  color='blue',
                x2_percentage=0.09, xlog=0, ylog=0, linestyle=':', linewidth=1.5)

    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYF19_sizes_vs_stellarmass.pdf')
    plt.close() 
  
    return
#end   sizes_vs_stellarmass_allz
    
    
    
    
    
def ssfr_hist(ThisRedshiftList):
    
    xlim=[-13.5,-8.5]
    ylim=[0.0, 0.35]
    bin_hist=0.1
        
    mass_bin=0.5
    mass_limits=[8.25,11.75]
    mass_bin_arr=np.arange(mass_limits[0],mass_limits[1]+mass_bin,mass_bin)    
        
        
    fig = plt.figure(figsize=(two_four_size_large[0],two_four_size_large[1]))
    grid = gridspec.GridSpec(2, 4)
    grid.update(wspace=0.0, hspace=0.0)
   
    for i_z in range(0,len(ThisRedshiftList)):        
         
        #DEFINE MODEL VARIABLES
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]         
        log_StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[i_z])      
        log_SSFR = np.zeros(len(G0_MR),dtype=np.float32)
        if(opt_rings==1):    
            for i_gal in range(0, len(G0_MR)):   
                sel = RingRadius<G0_MR['StellarHalfLightRadius'][i_gal]*1000.
                #sel = RingRadius<6.
                SFR = np.sum(G0_MR['SfrRings'][i_gal,sel])
                Mass = np.sum(G0_MR['DiskMassRings'][i_gal,sel])+np.sum(G0_MR['BulgeMassRings'][i_gal,sel])
                log_SSFR[i_gal]=np.log10(SFR/(Mass*1.e10/Hubble_h))            
        else:    
            log_SSFR=np.log10((G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h)))       
        #N_H[jj]=np.sum(G0['ColdGasRings_elements'][jj,sel,0]/1. , axis=0)  
        #MRII
        if(MRII==1):
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, i_z, FullSnapshotList_MRII)
            G0_MRII=G_MRII[sel]  
            log_StellarMass_MRII=stellar_mass_with_err(G0_MRII, Hubble_h, ThisRedshiftList[i_z])             
            log_SSFR_MRII = np.zeros(len(G0_MRII),dtype=np.float32)
            if(opt_rings==1):   
                for i_gal in range(0, len(G0_MRII)):   
                    sel = RingRadius<G0_MRII['StellarHalfLightRadius'][i_gal]*1000.
                    #sel = RingRadius<6.
                    SFR = np.sum(G0_MRII['SfrRings'][i_gal,sel])
                    Mass = np.sum(G0_MRII['DiskMassRings'][i_gal,sel])+np.sum(G0_MRII['BulgeMassRings'][i_gal,sel])
                    log_SSFR_MRII[i_gal]=np.log10(SFR/(Mass*1.e10/Hubble_h))                    
                else:    
                    log_SSFR_MRII=np.log10((G0_MRII['Sfr']/(G0_MRII['StellarMass']*1.e10/Hubble_h)))      
        
        #DEFINE Observational VARIABLES
        file=Datadir+'dr7_gal_final.fit'
        fits_table=fits.open(file)
        dr7_gal_final = fits_table[1]        
        #print(fits_table[1].columns)
        
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
            mass_low = mass_bin_arr[ii]-mass_bin/2.
            mass_high = mass_bin_arr[ii]+mass_bin/2.
            
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
                xlab='$\log_{10}(\mathrm{sSFR}/\mathrm{yr}^{-1})$'
            subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)     
                       
  
            #PLOT MODEL
            #MRII           
            if((mass_bin_arr[ii]<9.5) & MRII==1):
                sel=((log_StellarMass_MRII>mass_low) & (log_StellarMass_MRII<mass_high))
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
                sel=(log_StellarMass>mass_low) & (log_StellarMass<mass_high)
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
        
        
            #WRITE OUTPUT      
            if(write_to_file==1):
                df = pd.DataFrame({'log10_SSFR':x_axis, 'Fraction':y_axis})                
                file = Datadir + file_to_write + 'SSFR_hist' + str(f'_M{mass_low:0.2f}') + str(f'_{mass_high:0.2f}') + \
                       str(f'_z{ThisRedshiftList[i_z]:0.2f}')+'.csv'
                df.to_csv(file,index=False)
                #df = pd.read_csv(file)
                #subplot.plot(df['log10_SSFR'],df['Fraction'], color='black')
            
            
            #Previous models
            if do_previous_model2==1:           
                file = file_previous_model2 + '_ssfr_hist_' + str(f'{mass_low:0.1f}') + str(f'_{mass_high:0.1f}') + \
                       '_z0.10.txt' 
                model = Table.read(file, format='ascii')
                subplot.plot(model['col1'],model['col2'], color='red',linestyle=linestyle_previous_model2, linewidth=2)
            
            
            #PLOT OBSERVATIONS
            min_z=0.005
            max_z=0.2
            
            sel=((dr7_gal_final.data['jarle_median_mass']>mass_low) & (dr7_gal_final.data['jarle_median_mass']<mass_high) &
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
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.15, y_percentage=0.79, color='black', 
                            xlog=0, ylog=0, label=prefix_this_model, fontsize=12,  fontweight='normal')             
                plot_label (subplot, 'line', xlim, ylim, x_percentage=0.03, y_percentage=0.81,  color='red',
                            x2_percentage=0.13, xlog=0, ylog=0, linestyle='-', linewidth=2)
                
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.15, y_percentage=0.69, color='black', 
                            xlog=0, ylog=0, label=prefix_previous_model2, fontsize=12, fontweight='normal') 
                plot_label (subplot, 'line', xlim, ylim, x_percentage=0.03, y_percentage=0.71,  color='red',
                            x2_percentage=0.13, xlog=0, ylog=0, linestyle=linestyle_previous_model2, linewidth=2)
                
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.15, y_percentage=0.59, color='black', 
                            xlog=0, ylog=0, label='SDSS/DR7', fontsize=12, fontweight='normal')             
                plot_label (subplot, 'line', xlim, ylim, x_percentage=0.03, y_percentage=0.61, color='black',
                            x2_percentage=0.13, xlog=0, ylog=0, linestyle='-', linewidth=2)
        
            #LABELS
            mass_low, mass_high=mass_bin_arr[ii]-mass_bin/2.,mass_bin_arr[ii]+mass_bin/2.
            label="%0.2f" % mass_low + r'$<\mathrm{log_{10}}(M_*/$'+Msun+'$)<$' + "%0.2f" % mass_high
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.88, 
                        color='black', xlog=0, ylog=0, label=label, fontsize=11, fontweight='normal') 
        
        
        #endfor - loop on mass bins
        
        
        fits_table.close()
        
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYF19_ssfr_hist.pdf')
    plt.close()

    return   
#end ssfr_hist
   
    
    
def SFH(ThisRedshiftList):
     
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
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()

    plt.rcParams.update({'font.size': 14, 'xtick.labelsize': 14, 'ytick.labelsize': 14})
    return  
    
#end SFH
    
    
    
    
def cooling_heating(ThisRedshiftList):
       
    plt.rcParams.update({'font.size': 14, 'xtick.labelsize': 14, 'ytick.labelsize': 14})
      
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
               
    SSFR_cut=[-11.,-11., -10.5,-10.]
    plot_color=['black','dimgrey', 'darkgrey','lightgrey']
    xlim=[8.0,12.0]
    ylim=[-3.0, 3.]
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

    xlab='$\log_{10}(M_{\star}/$'+Msun+'$)$'       
    ylab='$\log_{10}(\mathrm{AGN\,Heating\,Rate}/\mathrm{Cooling\,Rate})$' 
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   

    #heating/cooling balance line    
    x_arr=np.arange(xlim[0],xlim[1]+0.05,0.05)     
    y_arr=x_arr*0.       
    subplot.plot(x_arr,y_arr,color='black', linestyle='--', linewidth=2)     
    
    for ii in range(0,len(ThisRedshiftList)):         
        #SCATTER 
        if(ii==0): 
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, ii, FullSnapshotList_MRII)        
            G0_MR_unsel=G_MRII[sel]            
            #G0_MR=G0_MR_unsel
            '''G0_MR=G0_MR_unsel[(np.log10(G0_MR_unsel['StellarMass']*1.e10*Hubble_h) > 8.) & 
                              (np.log10(G0_MR_unsel['BlackHoleMass']*1.e10) > 4.) &
                              (G0_MR_unsel['Type'] == 0)]'''
            G0_MR=G0_MR_unsel[(np.log10(G0_MR_unsel['BlackHoleMass']*1.e10) > 0.) & (G0_MR_unsel['Type'] == 0) &
                              (np.log10(G0_MR_unsel['StellarMass']/Hubble_h*1.e10) > xlim[0])]           
            #G0_MR=G0_MR_unsel[(G0_MR_unsel['Type'] == 0)]

            
            #sel=G0_MR['HaloID']==5000000000190
            #G0_MR=G0_MR[sel]
            
            Ngals = 2000
            if Ngals > len(G0_MR):
                Ngals = len(G0_MR) 
            G0_MR_rand=np.random.choice(G0_MR, size=Ngals, replace=False)    
            G0_MR=G0_MR_rand
            
            log_StellarMass=(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h))  
            #log_StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
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
            AGNheating = AGNcoeff * AGNrate
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
             
              
            sel=log_SSFR>SSFR_cut[ii]
            subplot.scatter(x_axis[sel], y_axis[sel], s=5, color='blue') 
            sel=log_SSFR<SSFR_cut[ii]           
            subplot.scatter(x_axis[sel], y_axis[sel], s=5, color='red') 
                        
            
            #SECOND XAXIS WITH MVIR  
            G0_MR=G0_MR_unsel[(np.log10(G0_MR_unsel['StellarMass']*1.e10) > 0.) & (G0_MR_unsel['Type'] == 0)]
            StellarMass=np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)
            Mvir=np.log10(G0_MR['Mvir']*1.e10/Hubble_h)
            (x_binned,median_MR,mean,pc16,pc84,rms)=median_and_percentiles(1.0,xlim[0],xlim[1],StellarMass,Mvir) 
            y_axis=median_MR            
            for jj in range(0,len(y_axis)):
                y_axis[jj] = int((y_axis[jj] * 10) + 0.5) / 10.0            
            #there are no 10^12 galaxies in MRII so their mean Mvir is taken from the database for MR
            y_axis[jj] = float('{:0.1f}'.format(np.log10(22221.01754115783*1e10)))
            #y_axis[jj-1] = float('{:0.1f}'.format(np.log10(5000.*1e10)))
            
            ax2 = subplot.twiny()             
            ax2.set_xbound(subplot.get_xbound())        
            ax2.set_xticks(x_binned)
            ax2.set_xticklabels(y_axis)
            
            xlab='$\log_{10}(<M_{\mathrm{200c}}>/$'+Msun+'$)$'  
            ax2.set_xlabel(xlab, fontsize=14)
            ax2.xaxis.set_label_position('top')  
                        
                    
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.065, y_percentage=0.84, 
                        color='red', xlog=0, ylog=0, label='Passive (z=0)', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.85, 
                        facecolors='red', xlog=0, ylog=0, sym='o', sym_size=20, err_size=0.) 

            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.065, y_percentage=0.8, 
                        color='blue', xlog=0, ylog=0, label='Star Forming (z=0)', 
                        fontsize=10, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.81, 
                        facecolors='blue', xlog=0, ylog=0, sym='o', sym_size=20, err_size=0.)


        #MEDIAN
        (sel)=select_current_redshift(G_MRII, ThisRedshiftList, ii, FullSnapshotList_MRII)        
        G0_MR_unsel=G_MRII[sel]   
        #G0_MR=G0_MR_unsel[(np.log10(G0_MR_unsel['BlackHoleMass']*1.e10) > 0.) & (G0_MR_unsel['Type'] == 0)]
        G0_MR=G0_MR_unsel[(np.log10(G0_MR_unsel['StellarMass']*1.e10/Hubble_h) > xlim[0]-1.) & 
                          (np.log10(G0_MR_unsel['BlackHoleMass']*1.e10/Hubble_h) > ylim[0]-1.) &
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
        #sel=((Cooling>0.) & (AGNheating>0.) & (log_StellarMass>9.0))    
        x_axis=log_StellarMass[sel]
        y_axis=np.log10(AGNheating[sel]/Cooling[sel])
        BlackHoleMass=np.log10(G0_MR['BlackHoleMass']*1.e10/Hubble_h)
        BlackHoleMass=BlackHoleMass[sel]
        
      
        sel = (y_axis>0.0) & (y_axis<0.2) & (x_axis>xlim[0]) & (x_axis<xlim[1])
        transition_mass = np.median(x_axis[sel])
        #print(transition_mass)
       
        #bin=0.1
        #sel = y_axis>0.        
        #(x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles (bin, xlim[0], xlim[1], x_axis, y_axis)    
        #sel=((x_binned<11.2) & (x_binned>8.5))
        #subplot.plot(x_binned[sel], median[sel], color='darkorange', linewidth=linewidth[ii],linestyle=linestyle[ii])
        
        #y-axis binning
        #bin=0.1
        #(x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles (bin, ylim[0], ylim[1], y_axis, x_axis)    
        #sel=(median>0.)
        #subplot.plot(median[sel], x_binned[sel], color='darkorange', linewidth=linewidth[ii],linestyle=linestyle[ii])
        #print(x_binned, median)
        
        
        #FIT LINE        
        '''if(ii==0):
            xx=np.arange(xxlim[0],xxlim[1],0.1)             
            #yy = 1./(2.5e13)*20**(xx*0.9)-0.5          
            #subplot.plot(xx, yy, color='darkorange', linewidth=linewidth[ii],linestyle=linestyle[ii]) 
            #yy = 1.9**(xx)/300.-0.5               
            #subplot.plot(xx, yy, color='red', linewidth=linewidth[ii],linestyle=linestyle[ii])  
            #yy = 3.7**(xx)/800000.-0.5       #increase base to increase slope, increse denominator to move to the right   
            #subplot.plot(xx, yy, color='blue', linewidth=linewidth[ii],linestyle=linestyle[ii])  
        
            def fit_func(x, a, b, c, d):
                return ((a**(x*d))/b+c)  
              
                
            sel=(x_axis<xlim[1]) & (x_axis>xlim[0]) #& (x_axis<xlim[1]-2.)                 
            p_init=[3.7,800000.,-0.5,1.0]
            par, err = curve_fit(fit_func, x_axis[sel],y_axis[sel],p_init, maxfev=10000)  
           
            xx=np.arange(xlim[0],xlim[1],0.1)  
            #par=[3.7,800000.,-0.5,1.0]
            yy = (par[0]**(xx*par[3]))/par[1]+par[2]
            print(par)
            subplot.plot(xx, yy, color='darkorange', linewidth=linewidth[ii],linestyle=linestyle[ii]) '''
                        
        
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
        
        
        #2-D region 
        bin=[0.4,0.4]   
        sel = (x_axis > xlim[0]) &  (x_axis < xlim[1]) &  (y_axis > ylim[0]) &  (y_axis < ylim[1])
        Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]
        H, xedges, yedges = np.histogram2d(x_axis[sel], y_axis[sel], bins=Nbins, range=[[xlim[0]-bin[0]/2., xlim[1]+bin[0]/2.], 
                                                                              [ylim[0]-bin[1]/2., ylim[1]+bin[1]/2.]]) 
        #H, xedges, yedges = np.histogram2d(x_axis, y_axis, bins=Nbins, range=[[xlim[0], xlim[1]], [ylim[0], ylim[1]]]) 
        extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
        plt.subplots_adjust(bottom=0.15, left=0.15)     
        mylevels = np.log10(np.logspace(-1.5, 0.0, num=5))[[1,3]]
        #mylevels = np.log10(np.logspace(-1.5, 0.0, num=6)) 
        HH=H
        H = zoom(H, 20)        
        H = np.log10(H/np.amax(H))        
        
        cont=subplot.contour(H.transpose()[::], origin='lower', colors=plot_color[ii], levels=mylevels, 
                             extent=extent, linestyles='solid')  
               
        '''
        #bin=[0.1,0.1]  
        #FOR HEN15
        #bin=[0.1,0.2]   
        #FOR AGN
        bin=[0.1,0.1]   
        Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]
        H, xedges, yedges = np.histogram2d(x_axis, y_axis, bins=Nbins, range=[[xlim[0]-bin[0]/2., xlim[1]+bin[0]/2.], 
                                                                              [ylim[0]-bin[1]/2., ylim[1]+bin[1]/2.]]) 
        #H, xedges, yedges = np.histogram2d(x_axis, y_axis, bins=Nbins, range=[[xlim[0], xlim[1]], [ylim[0], ylim[1]]]) 
        extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
        plt.subplots_adjust(bottom=0.15, left=0.15)     
        mylevels = np.log10(np.logspace(-2.0, 0.0, num=10)) 
        #mylevels = np.log10(np.logspace(-1.5, 0.0, num=6)) 
        HH=H
        H = zoom(H, 20)        
        H = np.log10(H/np.amax(H))   
        #if(ii==0):
        #    cont=subplot.contourf(H.transpose()[::], origin='lower', cmap='Greys_r', levels=mylevels, extent=extent)  
        
                   
        xx=np.arange(xlim[0],xlim[1],bin[0])   
        yy=np.arange(ylim[0],ylim[1],bin[1])          
        mean2d_arr=np.zeros(len(HH[:,0]))       
        for jj in range(0,len(HH[:,0])): 
            sel=(HH[jj,:]==np.amax(HH[jj,:]))               
            yy_aux=yy[sel]
            mean2d_arr[jj]=np.mean(yy_aux)
        
        
        #FOR HEN15
        if(ii==1):
            sel=xx>10.5
            mean2d_arr[sel]=mean2d_arr[sel]+0.75
        
        if(ii==0 or ii==1):
            sel=((xx<10.9) & (xx>8.0))
        else:
            if(ii==2):                
                sel=((xx<10.8) & (xx>8.0))
            else:
                sel=((xx<10.8) & (xx>8.0))
                
        #sel=((xx<11.) & (xx>9.0))
        xx=xx[sel]
        mean2d_arr=mean2d_arr[sel]
        
        if(ii==2):            
            sel=(xx>10.79) & (xx<10.81)
            mean2d_arr[sel]=mean2d_arr[sel]+0.5
        if(ii==3):            
            sel=(xx>10.79) & (xx<10.81)
            mean2d_arr[sel]=mean2d_arr[sel]+0.5
            
        
        for jj in range (0,5):
            for kk in range(1, len(mean2d_arr)-1):
                mean2d_arr[kk]=(mean2d_arr[kk-1]+mean2d_arr[kk]+mean2d_arr[kk+1])/3.
            
        #FOR AGN only
        ANG = 0
        if AGN==1:
            for jj in range (0,5):
                for kk in range(1, len(mean2d_arr)-1):
                    mean2d_arr[kk]=(mean2d_arr[kk-1]+mean2d_arr[kk]+mean2d_arr[kk+1])/3.
        
            if(ii<2):
                sel=((xx<10.9) & (xx>9.0))
            else:
                if(ii==2):
                    sel=((xx<10.7) & (xx>9.0))
                else:
                    sel=((xx<10.5) & (xx>9.0))
            xx=xx[sel]
            mean2d_arr=mean2d_arr[sel]
      
        #if(ii==0):        
        #    subplot.plot(xx-bin[0]/2., mean2d_arr, color='darkorange', linewidth=linewidth[ii],linestyle=linestyle[ii])
        ''' 
               
        #SECOND YAXIS with Black Hole Mass  
        '''if(ii==0):            
            (y_binned,median_MR,mean,pc16,pc84,rms)=median_and_percentiles(1.0,ylim[0],ylim[1],y_axis,BlackHoleMass) 
            y_axis_2=median_MR
       
            for jj in range(0,len(y_axis_2)):
                y_axis_2[jj] = int((y_axis_2[jj] * 10) + 0.5) / 10.0            
             
            ax2 = subplot.twinx()        
            ax2.set_ybound(subplot.get_ybound())        
            ax2.set_yticks(y_binned)
            ax2.set_yticklabels(y_axis_2)

            ylab='$\log_{10}(<M_{\mathrm{BH}}>/$'+Msun+'$)$'  
            ax2.set_ylabel(ylab, fontsize=14)
            ax2.yaxis.set_label_position('right')'''  
        
        #x_per=0.75
        #y_per=[0.69,0.75,0.81]
        #x_per=0.75
        #y_per=[0.59,0.65,0.71]
        #x_per=0.8
        #y_per=[0.33,0.45,0.56,0.71]
        #y_per=[0.54,0.65,0.75,0.84]
        
        #Hen15
        #x_per=0.65       
        #y_per=[0.50,0.6,0.73,0.8]
        
        #AGN
        x_per=0.65    
        y_per=[0.50,0.65,0.73,0.8]
        #y_per=[0.5,0.63,0.71,0.8]
        
        '''if ii==0:
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
                
                
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.26, y_percentage=0.055, 
                    color='black', xlog=0, ylog=0, label='Ridge line of the 2D distribution at different redshifts', 
                    fontsize=10, fontweight='normal') 
        
        plot_label (subplot, 'line', xlim, ylim, x_percentage=0.25, y_percentage=0.065, 
                    color='darkorange', x2_percentage=0.21, xlog=0, ylog=0, linestyle='-', linewidth=2)'''  
        
        snap_list=['61','42','34','29']
        model = 'Hen15_MRII'
        fa = Datadir+"cooling_heating_evo_transition_mass_at_snap_"+model+'_snap_'+snap_list[ii]+".txt"
        fa = open(fa, "r")                            
        fields = list(fa)        
        median = float(fields[0].strip())
        pc16 = float(fields[1].strip())
        pc84 = float(fields[2].strip())
        print(median,pc16,pc84)
        #subplot.plot([median,median], [-0.1,0.1], color=plot_color[ii], linewidth=4,linestyle='-')
        #subplot.plot([transition_mass,transition_mass], [0.05,0.2], color=plot_color[ii], linewidth=5,linestyle='-')
        print(transition_mass)
        
        y_value = [-2.45,-2.6,-2.75,-2.9]       
        subplot.fill_between([pc16,median,pc84],[y_value[ii],y_value[ii],y_value[ii]],
                             [y_value[ii]+0.15,y_value[ii]+0.15,y_value[ii]+0.15], 
                              facecolor=plot_color[ii], interpolate=True, alpha=0.4, edgecolor='black')  
                
        subplot.plot([median,median], [y_value[ii]+0.03,y_value[ii]+0.12], color=plot_color[ii], linewidth=2,linestyle='-') 
     
        
        
        
        #transition halo mass
        '''snap_list=['57','38','30','25']
        model = 'Hen15'
        fa = Datadir+"cooling_heating_evo_halo_mass_transition_mass_at_snap_"+model+'_snap_'+snap_list[ii]+".txt"
        fa = open(fa, "r")                            
        fields = list(fa)        
        median = float(fields[0].strip())
        pc16 = float(fields[1].strip())
        pc84 = float(fields[2].strip())
        print(median,pc16,pc84)
        #subplot.plot([median,median], [-0.1,0.1], color=plot_color[ii], linewidth=4,linestyle='-')
        #subplot.plot([transition_mass,transition_mass], [0.05,0.2], color=plot_color[ii], linewidth=5,linestyle='-')
        #print(transition_mass)
        
        #values need to be converted to stellar mass to put on plot
        if(ii==0):
            median=10.08
            pc16=9.5
            pc84=10.38
            #11.8134 11.52719 12.05525 
        if(ii==1):
            median=10.6
            pc16=10.5
            pc84=10.95
        #12.25494 12.09569 12.62711
        if(ii==2):
            median=10.75
            pc16=10.5
            pc84=11.05
        #12.43689 12.17839 12.83333
        if(ii==3):
            median=10.95
            pc16=10.6
            pc84=11.1
        #12.58611 12.26601 12.91644
        
        y_value = [-.45,-.6,-.75,-.9]       
        subplot.fill_between([pc16,median,pc84],[y_value[ii],y_value[ii],y_value[ii]],
                             [y_value[ii]+0.15,y_value[ii]+0.15,y_value[ii]+0.15], 
                              facecolor=plot_color[ii], interpolate=True, alpha=0.4, edgecolor='black')  

        subplot.plot([median,median], [y_value[ii]+0.03,y_value[ii]+0.12], color=plot_color[ii], linewidth=2,linestyle='-') 
           
        #11.8134 11.52719 12.05525        
        #12.25494 12.09569 12.62711
        #12.43689 12.17839 12.83333
        #12.58611 12.26601 12.91644
        '''
        
        #Hen15
        x_per=[0.67,0.67,0.55,0.43]
        y_per=[0.7,0.4,0.27,0.18]
        
        if ii==3:
            plot_label (subplot, 'label', xlim, ylim, x_percentage=x_per[3], y_percentage=y_per[3], 
                        color=plot_color[ii], xlog=0, ylog=0, label='z=3', fontsize=13, fontweight='bold')            
        if ii==2:
            plot_label (subplot, 'label', xlim, ylim, x_percentage=x_per[2], y_percentage=y_per[2], 
                        color=plot_color[ii], xlog=0, ylog=0, label='z=2', fontsize=13, fontweight='bold')            
        if ii==1:
            plot_label (subplot, 'label', xlim, ylim, x_percentage=x_per[1], y_percentage=y_per[1], 
                        color=plot_color[ii], xlog=0, ylog=0, label='z=1', fontsize=13, fontweight='bold')    
                 
        if ii==0:
            plot_label (subplot, 'label', xlim, ylim, x_percentage=x_per[0], y_percentage=y_per[0], 
                        color=plot_color[ii], xlog=0, ylog=0, label='z=0', fontsize=13, fontweight='bold')  
        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.08, y_percentage=0.925, 
                    color='black', xlog=0, ylog=0, label='2D distribution at different redshifts', 
                    fontsize=10, fontweight='normal') 
        
        plot_label (subplot, 'line', xlim, ylim, x_percentage=0.07, y_percentage=0.935, 
                    color='black', x2_percentage=0.03, xlog=0, ylog=0, linestyle='-', linewidth=2) 
        
                
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    #plt.savefig('./fig/HWL18_cooling_heating_MRII_new.pdf')
    plt.savefig('./fig/HWL18_cooling_heating_only_AGN_MRII_new.pdf')
    plt.close()

    
    plt.rcParams.update({'font.size': 14, 'xtick.labelsize': 14, 'ytick.labelsize': 14})
    
            
    return 
#end cooling_heating

    
    
    
    
    
    
    
    

    
def schmidt_kenn(ThisRedshiftList):
    
    xlim=[-5.0, 3.3]
    ylim=[-10.0, 1.7]
   
    fig = plt.figure(figsize=(7,6)) 
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
      
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
    
    xlab='$\mathrm{log_{10}}(\Sigma_{H_2})$'
    ylab='$\mathrm{log_{10}}(\Sigma_{Sfr})$'      
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   

    for ii in range(0,len(ThisRedshiftList)):        
                    
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]      
        G0_MR = G0_MR[G0_MR['Sfr']>10.]
        G0_MR=np.random.choice(G0_MR, size=1000) 
        
        sigma_SFR = np.zeros([len(G0_MR),RNUM])
        sigma_H2 = np.zeros([len(G0_MR),RNUM])
        area = np.zeros(RNUM)
        
        xx = 10**np.arange(-1.0,5.0,0.1)
        
        for jj in range(0, RNUM):
            if(jj==0):
                area = RingRadius[jj]*RingRadius[jj]*1e6
            else:
                area = (RingRadius[jj]*RingRadius[jj]-RingRadius[jj-1]*RingRadius[jj-1])*1e6
                
        for jj in range(0, len(G0_MR)):
            sigma_SFR[jj,:]=G0_MR['SfrRings'][jj,:]/(area/1e6)
            sigma_H2[jj,:] = G0_MR['ColdGasRings'][jj,:]*1e10/Hubble_h*G0_MR['H2fractionRings'][jj,:]/area
        
        subplot.scatter(np.log10(sigma_H2), np.log10(sigma_SFR), s=5, color='blue')
        
        a_sf = 0.004
        N_sf = 1.0
        sigma_crit = 70
        yy = a_sf*(xx/10.)*(1+xx/sigma_crit)**N_sf
        
        subplot.plot(np.log10(xx), np.log10(yy), color='black', linestyle='-')
        #print(np.log10(xx), np.log10(yy))
        
        #print(np.log10(sigma_H2), np.log10(sigma_SFR))
        
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()

    
    plt.rcParams.update({'font.size': 14, 'xtick.labelsize': 14, 'ytick.labelsize': 14})
    
    return    
#end BHBM_by_sfr    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#def BHBM_by_sfr(G_MR, G_MRII, ThisRedshiftList):    
def BHBM_by_sfr(ThisRedshiftList):
      
    plot_inset=1
    other_models=0
    
    '''xlim=[9.0,11.5]
    ylim=[-13.0, -7.]   
         
    fig = plt.figure(figsize=(10,8)) 
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
   
    xlab='$\mathrm{log_{10}}(M_{\star}/$'+Msun+'$)$'       
    ylab='$\mathrm{log_{10}}(M_{\mathrm{BH}}/$'+Msun+'$)$' 
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
   
    xlab='$\mathrm{log_{10}}(M_{\star}/$'+Msun+'$)$'       
    ylab='$\mathrm{log_{10}}(M_{\mathrm{BH}}/$'+Msun+'$)$' 
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
    
    xlab='$\mathrm{log_{10}}(M_{\star}/$'+Msun+'$)$'       
    ylab='$\mathrm{log_{10}}(M_{\mathrm{BH}}/$'+Msun+'$)$' 
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
        xlab='$\mathrm{log_{10}}(M_{\star}/$'+Msun+'$)$'       
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
           
        #(x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles (bin, ylim[0], ylim[1],log_BHMass, log_StellarMass)        
        #(inset_x_binned, inset_median, inset_mean, inset_pc16,inset_pc84, rms) = median_and_percentiles(bin,inset_ylim[0],inset_ylim[1],inset_y_axis,log_StellarMass) 
      
        
        #2d median BHBM 
        (x_binned, median)=median_2d(log_StellarMass,log_BHMass,xlim,ylim,[0.1,0.05])
      
    
        #2d median cooling heating     
        (inset_x_binned, inset_median)=median_2d(log_StellarMass,inset_y_axis,xlim,inset_ylim,[0.1,0.05])
          
        #color_ur=G0_MR['MagDust'][:,15]-G0_MR['MagDust'][:,17]  
        #Magr=G0_MR['MagDust'][:,17]-5.*np.log10(Hubble_h)
        #color_UV=G0_MR['MagDust'][:,0]-G0_MR['MagDust'][:,2]  
        #color_VJ=G0_MR['MagDust'][:,2]-G0_MR['MagDust'][:,7]   
        
                        
          
        G0_MR=np.random.choice(G0_MR, size=2000)        
        log_StellarMass=(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)) 
        log_BHMass=(np.log10(G0_MR['BlackHoleMass']*1.e10/Hubble_h)) 
        inset_y_axis=np.log10(G0_MR['BlackHoleMass']/(G0_MR['Mvir']**0.097)/((1+ThisRedshiftList[ii])**1.5))
        log_SSFR=np.log10((G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h)))
         
       
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
            x_arr=np.arange(9.97,10.63,0.005)  
            #x_arr=np.arange(9.74,9.95,0.005)  
            #z=0
            (slope0,b0)=get_slope(x1=9.,y1=5.85,x2=11.5,y2=6.57)
            #z=2
            (slope2,b2)=get_slope(x1=9.,y1=6.6,x2=11.5,y2=7.35)
            if(other_models==0):
                subplot.fill_between(x_arr-2.*np.log10(Hubble_h),x_arr*slope0+b0,x_arr*slope2+b2, facecolor='lightgrey', 
                                     interpolate=True, alpha=0.4, edgecolor='black')   
                      
        if(ii==0):                    
            (slope,b)=get_slope(x1=9.,y1=5.85,x2=11.5,y2=6.57)
        else:
            if(ii==1):
                (slope,b)=get_slope(x1=9.,y1=6.27,x2=11.5,y2=6.97)
            else:
                (slope,b)=get_slope(x1=9.,y1=6.6,x2=11.5,y2=7.35)
            
            
        x_arr=np.arange(xlim[0]+2.*np.log10(Hubble_h),xlim[1]+0.05,0.05)     
        y_arr=x_arr*slope+b    
        if(other_models==0):
            subplot.plot(x_arr-2.*np.log10(Hubble_h),y_arr,color=plot_color[ii], linestyle='--', linewidth=2)         
                    
        #median        
        #sel=(x_binned>4.7) & (x_binned<8.5) &(median>xlim[0])         
        #subplot.plot(median[sel],x_binned[sel], color=plot_color[ii], linewidth=2,linestyle='-')        
        subplot.plot(x_binned, median, color=plot_color[ii], linewidth=2,linestyle='-')
        
        
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
            #sel=(inset_x_binned>-5.) & (inset_x_binned<-2.) &(inset_median>inset_xlim[0])      
            #inset.plot(inset_median[sel],inset_x_binned[sel], color=plot_color[ii], linewidth=2,linestyle='-')            
            inset.plot(inset_x_binned, inset_median, color=plot_color[ii], linewidth=2,linestyle='-')
      
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
                #label='average quenching mass (0<z<2)'
                #plot_label (subplot, 'label', xlim, ylim, x_percentage=0.52, y_percentage=0.3, 
                #            color='black', xlog=0, ylog=0, label=label, fontsize=12, fontweight='normal',
                #            rotation=5, backgroundcolor='none') 
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
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HWL18_plots_bhbm_by_sfr.pdf')
    plt.close()
    
    plt.rcParams.update({'font.size': 14, 'xtick.labelsize': 14, 'ytick.labelsize': 14})
    
    return    
#end BHBM_by_sfr





def AGN_quenching(ThisRedshiftList):
    
    local_slope_color_cut=[0.075,0.3, 0.32]
    local_offset_color_cut=[1.85,1.18,0.99]
    #local_slope_color_cut=[0.075,0.48, 0.38]
    #local_offset_color_cut=[2.0,1.0,1.0]
    
    SSFR_cut=[-11.,-10.5,-10.]
    
    xlim=[9.0,12.0]
    ylim=[4., 10.5]   
        
    plot_color=['red','orange','green','blue']        
   
    fig = plt.figure(figsize=(7,6)) 
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
      
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
    
    xlab='$\mathrm{log_{10}}(M_{\star}/$'+Msun+'$)$'       
    ylab='$\mathrm{log_{10}}(M_{\mathrm{BH}}/$'+Msun+'$)$' 
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   

    
    for ii in range(0,len(ThisRedshiftList)):        
                    
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR_unsel=G_MR[sel]        
        G0_MR=G0_MR_unsel[(np.log10(G0_MR_unsel['StellarMass']*1.e10/Hubble_h) > xlim[0]-1.) & 
                          (np.log10(G0_MR_unsel['BlackHoleMass']*1.e10/Hubble_h) > ylim[0]-1.) &
                          (G0_MR_unsel['Type'] == 0)]
        
        '''G0_MR=np.random.choice(G0_MR, size=2000)           
        x_axis=(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h))        
        y_axis=(np.log10(G0_MR['BlackHoleMass']*1.e10/Hubble_h)) 
        log_SSFR=np.log10((G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h)))
        
        if(ii==1):       
                sel=log_SSFR<SSFR_cut[ii]
                subplot.scatter(x_axis[sel], y_axis[sel], s=2.5, color='red') 
                sel=log_SSFR>SSFR_cut[ii]
                subplot.scatter(x_axis[sel], y_axis[sel], s=2.5, color='blue')'''   
        
        
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR_unsel=G_MR[sel]        
        G0_MR=G0_MR_unsel[(np.log10(G0_MR_unsel['StellarMass']*1.e10/Hubble_h) > xlim[0]-1.) & 
                          (np.log10(G0_MR_unsel['BlackHoleMass']*1.e10/Hubble_h) > ylim[0]-1.) &
                          (G0_MR_unsel['Type'] == 0)]
        x_axis=(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h))        
        y_axis=(np.log10(G0_MR['BlackHoleMass']*1.e10/Hubble_h)) 
        
        #OLD MEDIAN
        '''bin=0.1
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
        subplot.plot(x[sel], y[sel], color=plot_color[ii], linewidth=2,linestyle='-') '''   
            
        
        #NEW MEDIAN
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
        sel=((xx<11.5) & (xx>9.0))
        subplot.plot(xx[sel], mean2d_arr[sel], color=plot_color[ii], linewidth=2,linestyle='-')
        
        
        
        
        
        #inset 
        ellipse = Ellipse(xy=(9.2, 5.2), width=3.0, height=1.5, angle=40,
                          fc='lightblue', edgecolor='lightblue', lw=2, alpha=0.3)
        subplot.add_patch(ellipse)   
        ellipse = Ellipse(xy=(11.05, 8.15), width=4.2, height=1.2, angle=72,
                          fc='pink', edgecolor='pink', lw=2, alpha=0.3)
        subplot.add_patch(ellipse)   
        ##subplot.plot(median[sel],x_binned[sel], color=plot_color[ii], linewidth=2,alpha=0.3)
        
        (slope,b)=get_slope(x1=9.,y1=7.5,x2=11.5,y2=4.5)
        x_arr=np.arange(xlim[0],xlim[1]+0.05,0.05)     
        y_arr=x_arr*slope+b          
        subplot.plot(x_arr,y_arr,color='black', linestyle=':', linewidth=2) 
           
        if(ii>0):
            x1=9.07
            x2=10.3   
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

            xlab='$\mathrm{log_{10}}(<M_{\mathrm{vir}}>/$'+Msun+'$)$'  
            ax2.set_xlabel(xlab, fontsize=14)
            ax2.xaxis.set_label_position('top') 
            
            
            label='SN feedback'
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.12, y_percentage=0.69, 
                        color='black', xlog=0, ylog=0, label=label, fontsize=15, fontweight='normal')
            label='efficiency threshold'
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.64, 
                        color='black', xlog=0, ylog=0, label=label, fontsize=15, fontweight='normal')
                     
                
            label='SN feedback regulated'
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.04, y_percentage=0.23, 
                        color='blue', xlog=0, ylog=0, label=label, fontsize=13, fontweight='normal',rotation=15)
            label='cold-mode accretion'
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.07, y_percentage=0.18, 
                        color='blue', xlog=0, ylog=0, label=label, fontsize=13, fontweight='normal',rotation=15)            
                     
            label='in-situ star formation'
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.37, y_percentage=0.12, 
                        color='blue', xlog=0, ylog=0, label=label, fontsize=13, fontweight='normal')            
            label='disk formation'
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.41, y_percentage=0.07, 
                        color='blue', xlog=0, ylog=0, label=label, fontsize=13, fontweight='normal')
            label='weak BH growth'
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.40, y_percentage=0.02, 
                        color='blue', xlog=0, ylog=0, label=label, fontsize=13, fontweight='normal')
                   
        
            label='AGN feedback regulated'
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.48, y_percentage=0.82, 
                        color='red', xlog=0, ylog=0, label=label, fontsize=13, fontweight='normal',rotation=37) 
            label='hot-mode accretion'
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.54, y_percentage=0.75, 
                        color='red', xlog=0, ylog=0, label=label, fontsize=13, fontweight='normal',rotation=37) 
                        
            label='significant merger-driven growth'
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.21, y_percentage=0.9, 
                        color='red', xlog=0, ylog=0, label=label, fontsize=13, fontweight='normal') 
            label='bulge formation'
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.33, y_percentage=0.85, 
                        color='red', xlog=0, ylog=0, label=label, fontsize=13, fontweight='normal') 
            label='strong BH growth'
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.32, y_percentage=0.8, 
                        color='red', xlog=0, ylog=0, label=label, fontsize=13, fontweight='normal') 
     
            
    
            label='quenching threshold z=0'
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.62, y_percentage=0.41, 
                        color=plot_color[0], xlog=0, ylog=0, label=label, fontsize=12, fontweight='normal',
                        rotation=5, backgroundcolor='none')         
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.93, y_percentage=0.48, 
                        color=plot_color[1], xlog=0, ylog=0, label='z=1', fontsize=12, fontweight='normal',
                        rotation=5, backgroundcolor='none')         
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.93, y_percentage=0.535, 
                        color=plot_color[2], xlog=0, ylog=0, label='z=2', fontsize=12, fontweight='normal',
                        rotation=5, backgroundcolor='none')
            
            
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/plots_AGN_quenching.pdf')
    plt.close()
    
    plt.rcParams.update({'font.size': 14, 'xtick.labelsize': 14, 'ytick.labelsize': 14})
    
    return    
#end AGN_quenching


def growth_channels(ThisRedshiftList):
    
    xlim=[9.5,12.5]
    ylim=[0.0,1.2]
    bin=0.2
    
    model_to_print='Hen15_other_'
    sim='MRII'
    
    for ii in range(0,len(ThisRedshiftList)):        
        
        #MODEL
        #MR
        (sel)=select_current_redshift(G_MRII, ThisRedshiftList, ii, FullSnapshotList_MRII)         
        G0_MR=G_MRII[sel]                   
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
        StellarMass = np.log10(G0_MR['StellarMass']*1e10/Hubble_h)       
        StellarMassICM = np.log10(G0_MR['StellarMass']*1e10/Hubble_h+G0_MR['ICM']*1e10/Hubble_h)
        HaloMass = np.log10(G0_MR['Mvir']*1e10/Hubble_h)  
        #StellarMass=np.log10(G0_MR['Mvir']*1.e10/Hubble_h)
        InSituRatio=G0_MR['MassFromInSitu']/G0_MR['StellarMass']
        MergersRatio=G0_MR['MassFromMergers']/G0_MR['StellarMass']
        
        InSituICMRatio=G0_MR['MassFromInSitu']/(G0_MR['StellarMass']+G0_MR['ICM'])
        MergersICMRatio=(G0_MR['MassFromMergers']+G0_MR['ICM'])/(G0_MR['StellarMass']+G0_MR['ICM'])
           
       
            
        '''(x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles (bin, xlim[0], xlim[1], StellarMass, InSituRatio) 
        #sel=median!=0.
        #subplot.plot(x_binned[sel], median[sel],color='blue', linewidth=2)     
        
        x_axis=x_binned
        y_axis=median
        fa = open(Datadir+"growth_channels_"+model_to_print+sim+"_insitu_median_z0.0.txt", "w")
        fa.write("%d\n" % len(x_axis))
        for kk in range (0,len(x_axis)):               
            fa.write("%0.2f " % x_axis[kk] + "%0.2f\n" % y_axis[kk])     
        fa.close() 
        
        (x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles (bin,xlim[0],xlim[1],StellarMass,InSituICMRatio) 
        
        x_axis=x_binned
        y_axis=median
        fa = open(Datadir+"growth_channels_"+model_to_print+sim+"_insitu_with_ICM_median_z0.0.txt", "w")
        fa.write("%d\n" % len(x_axis))
        for kk in range (0,len(x_axis)):               
            fa.write("%0.2f " % x_axis[kk] + "%0.2f\n" % y_axis[kk])     
        fa.close() 
        
        (x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles (bin,xlim[0],xlim[1],StellarMassICM,InSituICMRatio) 
        
        x_axis=x_binned
        y_axis=median
        fa = open(Datadir+"growth_channels_"+model_to_print+sim+"_insitu_with_ICM_median_ICM_mass_z0.0.txt", "w")
        fa.write("%d\n" % len(x_axis))
        for kk in range (0,len(x_axis)):               
            fa.write("%0.2f " % x_axis[kk] + "%0.2f\n" % y_axis[kk])     
        fa.close() '''
        
        
        
        x_variable = StellarMass    
        xx_variable = StellarMassICM             
        y_variable_1 = MergersRatio
        y_variable_2 = MergersICMRatio
        y_variable_3 = InSituRatio
        y_variable_4 = InSituICMRatio
        
        Nbins=int((xlim[1]-xlim[0])/bin+1)  
        median_xx=np.zeros(Nbins, np.float32) 
        median_x=np.zeros(Nbins, np.float32) 
        median_1=np.zeros(Nbins, np.float32) 
        median_2=np.zeros(Nbins, np.float32) 
        median_3=np.zeros(Nbins, np.float32) 
        median_4=np.zeros(Nbins, np.float32) 
        x_min=xlim[0]-bin/2.
  
        for ii in range(0,Nbins):
            x_variable_sel=x_variable[(x_variable > (x_min+(ii*bin))) & (x_variable < (x_min+(ii+1.0)*bin))] 
            xx_variable_sel=xx_variable[(x_variable > (x_min+(ii*bin))) & (x_variable < (x_min+(ii+1.0)*bin))] 
            y_variable_sel_1=y_variable_1[(x_variable > (x_min+(ii*bin))) & (x_variable < (x_min+(ii+1.0)*bin))] 
            y_variable_sel_2=y_variable_2[(x_variable > (x_min+(ii*bin))) & (x_variable < (x_min+(ii+1.0)*bin))] 
            y_variable_sel_3=y_variable_3[(x_variable > (x_min+(ii*bin))) & (x_variable < (x_min+(ii+1.0)*bin))] 
            y_variable_sel_4=y_variable_4[(x_variable > (x_min+(ii*bin))) & (x_variable < (x_min+(ii+1.0)*bin))] 
               
            if(len(x_variable_sel) > 0):           
                median_1[ii]=np.median(y_variable_sel_1)
                median_2[ii]=np.median(y_variable_sel_2)
                median_3[ii]=np.median(y_variable_sel_3)
                median_4[ii]=np.median(y_variable_sel_4)
                median_xx[ii]=np.median(xx_variable_sel)
                median_x[ii]=np.median(x_variable_sel)
           
        x_binned=np.arange(Nbins)*((x_min+(Nbins*bin))-(x_min+(Nbins*0.0)))/(Nbins*1.)+x_min+bin/2.    
        
        x_axis=x_binned
        y_axis=median_1
        fa = open(Datadir+"growth_channels_"+model_to_print+sim+"_mergers_median_z0.0.txt", "w")
        fa.write("%d\n" % len(x_axis))
        for kk in range (0,len(x_axis)):               
            fa.write("%0.2f " % x_axis[kk] + "%0.2f\n" % y_axis[kk])     
        fa.close() 
        
        x_axis=x_binned
        y_axis=median_2
        fa = open(Datadir+"growth_channels_"+model_to_print+sim+"_mergers_with_ICM_median_z0.0.txt", "w")
        fa.write("%d\n" % len(x_axis))
        for kk in range (0,len(x_axis)):               
            fa.write("%0.2f " % x_axis[kk] + "%0.2f\n" % y_axis[kk])     
        fa.close()
                
        x_axis = median_xx
        y_axis = median_2
        fa = open(Datadir+"growth_channels_"+model_to_print+sim+"_mergers_with_ICM_median_ICM_mass_z0.0.txt", "w")
        fa.write("%d\n" % len(x_axis))
        for kk in range (0,len(x_axis)):               
            fa.write("%0.2f " % x_axis[kk] + "%0.2f\n" % y_axis[kk])     
        fa.close()
              
        
        x_axis=x_binned
        y_axis=median_3
        fa = open(Datadir+"growth_channels_"+model_to_print+sim+"_insitu_median_z0.0.txt", "w")
        fa.write("%d\n" % len(x_axis))
        for kk in range (0,len(x_axis)):               
            fa.write("%0.2f " % x_axis[kk] + "%0.2f\n" % y_axis[kk])     
        fa.close() 
        
        x_axis=x_binned
        y_axis=median_4
        fa = open(Datadir+"growth_channels_"+model_to_print+sim+"_insitu_with_ICM_median_z0.0.txt", "w")
        fa.write("%d\n" % len(x_axis))
        for kk in range (0,len(x_axis)):               
            fa.write("%0.2f " % x_axis[kk] + "%0.2f\n" % y_axis[kk])     
        fa.close()
        
        x_axis = median_xx
        y_axis = median_4
        fa = open(Datadir+"growth_channels_"+model_to_print+sim+"_insitu_with_ICM_median_ICM_mass_z0.0.txt", "w")
        fa.write("%d\n" % len(x_axis))
        for kk in range (0,len(x_axis)):               
            fa.write("%0.2f " % x_axis[kk] + "%0.2f\n" % y_axis[kk])     
        fa.close()
        
    
        #SECOND AXIS WITH MVIR   
        G0_MR=G0_MR[G0_MR['Type']==0]
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[0])            
        Mvir=np.log10(G0_MR['Mvir']*1.e10/Hubble_h)      
        (x_binned,median_MR,mean,pc16,pc84,rms)=median_and_percentiles(1.0,xlim[0],xlim[1],StellarMass,Mvir) 
        y_axis=median_MR
       
        for jj in range(0,len(y_axis)):
            y_axis[jj] = int((y_axis[jj] * 10) + 0.5) / 10.0 
            
        fa = open(Datadir+"growth_channels_"+model_to_print+"_mvir_axis.txt", "w")
        fa.write("%d\n" % len(x_binned))
        for kk in range (0,len(x_binned)):               
            fa.write("%0.2f " % x_binned[kk] + "%0.2f\n" % y_axis[kk])     
        fa.close()     
    
    
    
        '''#halo mass
    
        x_variable = HaloMass    
        xx_variable = HaloMass             
        y_variable_1 = MergersRatio
        y_variable_2 = MergersICMRatio
        y_variable_3 = InSituRatio
        y_variable_4 = InSituICMRatio
        
        bin=0.5
        Nbins=int((15.-11.)/bin+1)  
        median_xx=np.zeros(Nbins, np.float32) 
        median_x=np.zeros(Nbins, np.float32) 
        median_1=np.zeros(Nbins, np.float32) 
        median_2=np.zeros(Nbins, np.float32) 
        median_3=np.zeros(Nbins, np.float32) 
        median_4=np.zeros(Nbins, np.float32) 
        x_min=11.-bin/2.
  
        for ii in range(0,Nbins):
            x_variable_sel=x_variable[(x_variable > (x_min+(ii*bin))) & (x_variable < (x_min+(ii+1.0)*bin))] 
            xx_variable_sel=xx_variable[(x_variable > (x_min+(ii*bin))) & (x_variable < (x_min+(ii+1.0)*bin))] 
            y_variable_sel_1=y_variable_1[(x_variable > (x_min+(ii*bin))) & (x_variable < (x_min+(ii+1.0)*bin))] 
            y_variable_sel_2=y_variable_2[(x_variable > (x_min+(ii*bin))) & (x_variable < (x_min+(ii+1.0)*bin))] 
            y_variable_sel_3=y_variable_3[(x_variable > (x_min+(ii*bin))) & (x_variable < (x_min+(ii+1.0)*bin))] 
            y_variable_sel_4=y_variable_4[(x_variable > (x_min+(ii*bin))) & (x_variable < (x_min+(ii+1.0)*bin))] 
               
            if(len(x_variable_sel) > 0):           
                median_1[ii]=np.median(y_variable_sel_1)
                median_2[ii]=np.median(y_variable_sel_2)
                median_3[ii]=np.median(y_variable_sel_3)
                median_4[ii]=np.median(y_variable_sel_4)
                median_xx[ii]=np.median(xx_variable_sel)
                median_x[ii]=np.median(x_variable_sel)
           
        x_binned=np.arange(Nbins)*((x_min+(Nbins*bin))-(x_min+(Nbins*0.0)))/(Nbins*1.)+x_min+bin/2.    
        
        save_x1 = median_xx
        save_y1 = median_1
        save_x3 = median_xx
        save_y3 = median_3
        
        print(save_x1, save_y1)
        print(save_x3, save_y3)'''
        
        
    #PLOTS   
  
    fig = plt.figure(figsize=(one_four_size_large[0],one_four_size_large[1]))
    grid = gridspec.GridSpec(1, 4)
    grid.update(wspace=0.0, hspace=0.0)
    
    model_names=['Hen15_no_feedback_','Hen15_only_AGN_','Hen15_only_SN_','Hen15_']
    model_label=['no feedback','only AGN','only SN','Hen15']
        
    for ii in range(0,4):
        
        subplot=plt.subplot(grid[ii])
  
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        xlab='$\log_{10}(M_{\star}/$'+Msun+'$)$'    
        subplot.set_xlabel(xlab, fontsize=14)
      
        if ii==0:
            ylab='Fraction'  
            subplot.set_ylabel(ylab, fontsize=14)
                       
    
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))
        subplot.yaxis.set_major_locator(MultipleLocator(0.5))
        subplot.yaxis.set_minor_locator(MultipleLocator(0.25))   
              
        
        if ii>0:
            plt.tick_params(axis='y', which='both', left='on', labelleft='off')
            
        
        fa = Datadir+"growth_channels_"+model_names[ii]+sim+"_insitu_median_z0.0.txt"
        (x_axis,y_axis)=read_file(fa)        
        subplot.plot(x_axis, y_axis,color='blue', linewidth=2)   
                
        fa = Datadir+"growth_channels_"+model_names[ii]+sim+"_mergers_median_z0.0.txt"       
        (x_axis,y_axis)=read_file(fa)
        subplot.plot(x_axis, y_axis,color='red', linewidth=2)
        
        fa = Datadir+"growth_channels_"+model_names[ii]+sim+"_insitu_median_z0.0.txt"
        (x_axis,y_axis)=read_file(fa)        
        subplot.scatter(x_axis, y_axis, marker='o', color='blue', s=20)
       
        fa = Datadir+"growth_channels_"+model_names[ii]+sim+"_mergers_median_z0.0.txt"
        (x_axis,y_axis)=read_file(fa)
        subplot.scatter(x_axis, y_axis, marker='o', color='red', s=20)   
        
        #with ICM   
        fa = Datadir+"growth_channels_"+model_names[ii]+sim+"_insitu_with_ICM_median_ICM_mass_z0.0.txt"
        (x_axis,y_axis)=read_file(fa)
        subplot.scatter(x_axis, y_axis, marker='o', facecolors='none', edgecolors='blue', s=20)   
        #subplot.plot(x_axis, y_axis,color='blue', linewidth=2, linestyle=':')
       
        fa = Datadir+"growth_channels_"+model_names[ii]+sim+"_mergers_with_ICM_median_ICM_mass_z0.0.txt"
        (x_axis,y_axis)=read_file(fa)
        subplot.scatter(x_axis, y_axis, marker='o', facecolors='none', edgecolors='red',s=20)
        #subplot.plot(x_axis, y_axis,color='red', linewidth=2, linestyle=':')
   
        fa = Datadir+"growth_channels_"+model_names[ii]+sim+"_insitu_with_ICM_median_ICM_mass_z0.0.txt"
        (x_axis,y_axis)=read_file(fa)
        subplot.plot(x_axis, y_axis,color='blue', linewidth=2, linestyle=':')
        #subplot.scatter(x_axis, y_axis, marker='o', facecolors='none', edgecolors='blue', s=20)   
          
        fa = Datadir+"growth_channels_"+model_names[ii]+sim+"_mergers_with_ICM_median_ICM_mass_z0.0.txt"
        (x_axis,y_axis)=read_file(fa)
        subplot.plot(x_axis, y_axis,color='red', linewidth=2, linestyle=':')
        #subplot.scatter(x_axis, y_axis, marker='o', facecolors='none', edgecolors='red',s=20)
    
    
    
        #Observations
        xx = [11.25]
        yy = [0.32]
        subplot.scatter(xx, yy, marker='x', facecolors='blue', edgecolors='blue', s=50)   
        yy = [0.68]
        subplot.scatter(xx, yy, marker='x', facecolors='red', edgecolors='red', s=50)   
        if ii==3:
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.02, y_percentage=0.45, color='black', 
                        xlog=0, ylog=0, label='Ownsworth (2014)', fontsize=9, fontweight='normal') 
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.11, y_percentage=0.38, color='black', 
                        xlog=0, ylog=0, label='in-situ', fontsize=9, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.08, y_percentage=0.395, 
                        xlog=0, ylog=0, sym='x', facecolors='blue', edgecolors='blue', sym_size=30, err_size=0)
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.11, y_percentage=0.32, color='black', 
                        xlog=0, ylog=0, label='accreted', fontsize=9, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.08, y_percentage=0.335, 
                        xlog=0, ylog=0, sym='x', facecolors='red', edgecolors='red', sym_size=30, err_size=0)   
            
        #SECOND AXIS WITH MVIR   
        #G0_MR=G0_MR[G0_MR['Type']==0]
        #StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[0])            
        #Mvir=np.log10(G0_MR['Mvir']*1.e10/Hubble_h)
        #(x_binned,median_MR,mean,pc16,pc84,rms)=median_and_percentiles(1.0,xlim[0],xlim[1],StellarMass,Mvir) 
        #y_axis=median_MR
       
        #for jj in range(0,len(y_axis)):
        #    y_axis[jj] = int((y_axis[jj] * 10) + 0.5) / 10.0            
             
        fa = Datadir+"growth_channels_"+model_names[ii]+sim+"_mvir_axis.txt"
        #fa = Datadir+"growth_channels_"+model_to_print+"_mvir_axis.txt"
        (x_axis,y_axis)=read_file(fa)    
            
        ax2 = subplot.twiny()        
        ax2.set_xbound(subplot.get_xbound())        
        ax2.set_xticks(x_axis)
        ax2.set_xticklabels(y_axis)

        xlab='$\log_{10}(<M_{\mathrm{200c}}>/$'+Msun+'$)$'  
        ax2.set_xlabel(xlab, fontsize=14)
        ax2.xaxis.set_label_position('top') 
             
        #subplot.tick_params(axis='x', which='both', bottom='off', labelbottom='off')
            
            
        #LABELS
        if ii==0:
                
            plot_label (subplot, 'label', xlim, ylim, 
                x_percentage=0.05, y_percentage=0.92, color='black', xlog=0, ylog=0, 
                label='Growth from:', fontsize=10, fontweight='normal') 
            
            plot_label (subplot, 'label', xlim, ylim, 
                    x_percentage=0.15, y_percentage=0.86, color='black', xlog=0, ylog=0, 
                    label='in-situ star formation', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim,
                    x_percentage=0.05, y_percentage=0.88, color='blue', x2_percentage=0.12, 
                    xlog=0, ylog=0, linestyle='-', linewidth=2)
        
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.085, y_percentage=0.88, 
                        xlog=0, ylog=0, sym='o', facecolors='blue', edgecolors='blue', sym_size=20, err_size=0)   
        
        
            plot_label (subplot, 'label', xlim, ylim, 
                    x_percentage=0.15, y_percentage=0.80, color='black', xlog=0, ylog=0, 
                    label='accreted', fontsize=10, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim,
                    x_percentage=0.05, y_percentage=0.82, color='red', x2_percentage=0.12, 
                    xlog=0, ylog=0, linestyle='-', linewidth=2)
            
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.085, y_percentage=0.82, 
                        xlog=0, ylog=0, sym='o', facecolors='red', edgecolors='red', sym_size=20, err_size=0)   
            
            
            plot_label (subplot, 'label', xlim, ylim, 
                    x_percentage=0.15, y_percentage=0.74, color='black', xlog=0, ylog=0, 
                    label='including ICM', fontsize=10, fontweight='normal') 
            
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.03, y_percentage=0.76, 
                        color='red', x2_percentage=0.085, xlog=0, ylog=0, linestyle=':', linewidth=2)
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.085, y_percentage=0.76, 
                        color='blue', x2_percentage=0.14, xlog=0, ylog=0, linestyle=':', linewidth=2)
                
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.06, y_percentage=0.76, 
                        xlog=0, ylog=0, sym='o', facecolors='none', edgecolors='red', sym_size=20, err_size=0)   
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.11, y_percentage=0.76, 
                        xlog=0, ylog=0, sym='o', facecolors='none', edgecolors='blue', sym_size=20, err_size=0)   
     
            #subplot.scatter(9.8,0.80, marker='o', s=20, facecolors='none', edgecolors='red')
            
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.05, 
                    color='black', xlog=0, ylog=0, label=model_label[ii], fontsize=10, fontweight='normal') 
            
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HWL18_growth_channels_MRII.pdf')
    plt.close()
   
    return
#end growth_channels   

def bluck_red_fractions(ThisRedshiftList):
        
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
    xlab='$log_{10}(M_*/(h^{-2}$'+Msun+'$))$'           
    ylab='$log_{10}(\phi /(h^3 Mpc^{-3} \mathrm{dex}))$'               
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
    xlab='$log_{10}(M_{BH}/$'+Msun+'$)$'           
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
    xlab='$log_{10}(M_{200c}/$'+Msun+'$)$'           
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
    xlab='$log_{10}(M_{BH}/$'+Msun+'$)$'           
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
    xlab='$log_{10}(M_{200c}/$'+Msun+'$)$'           
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
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/plots_bluck_red_fractions.pdf')
    plt.close()

    return   
#end bluck_red_fractions



def satellite_quench_times(ThisRedshiftList):
     
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
            
        label="$\mathrm{log_{10}}(M_{\mathrm{Halo}}/$'+Msun+'$)=[12.5,13.5,14.5]$"      
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.8, 
                    color='black', xlog=0, ylog=0, label=label, fontsize=10, fontweight='normal') 
        label="$\mathrm{log_{10}}(M_{\mathrm{sat}}/$'+Msun+'$)=[9.0,9.5,10.0]$"      
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
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ''' plt.rcParams.update({'font.size': 14, 'xtick.labelsize': 14, 'ytick.labelsize': 14, 'axes.linewidth': 2, 
                     'xtick.major.size': 6, 'xtick.major.width': 1.5, 
                     'ytick.major.size': 6, 'ytick.major.width': 1.5, 
                     'xtick.minor.size': 3, 'xtick.minor.width': 1.,                   
                     'ytick.minor.size': 3, 'ytick.minor.width': 1.})  
        
    plot_xy=0
    
    if(plot_xy==1):
        plt.rcParams.update({'font.size': 14, 'xtick.labelsize': 14, 'ytick.labelsize': 14})   
        xlim=[75.0,90.]
        ylim=[85.0,100.0]
        fig = plt.figure(figsize=(10,10))      
        subplot=plt.subplot()  
    
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        xlab='$\mathrm{log_{10}}(M_*/$'+Msun+'$)$'      
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
       

        #pdf.savefig()
        #plt.close()'''

    return       
#end satellite_quench_times




def sat_fraction(ThisRedshiftList):
    
    sat_numbers=0
    
    
    xlim=[8.5,11.5]
    ylim=[0., 1.] 
   
    fig = plt.figure(figsize=(15,4))
    grid = gridspec.GridSpec(1, 5)
    grid.update(wspace=0.0, hspace=0.0)
    
    for i_z in range (0,len(ThisRedshiftList)):
           
        char_redshift="%0.2f" % ThisRedshiftList[i_z]
            
        subplot=plt.subplot(grid[i_z])
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        
        xlab='$\mathrm{log_{10}}(M_*/$'+Msun+'$)$'      
        if i_z==0:
            ylab='Satellite Fraction'
        else:
            ylab=''      
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)
        
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
        subplot.yaxis.set_minor_locator(MultipleLocator(0.1))
        
        if i_z>0:
            plt.tick_params(axis='y', which='both', left='on', labelleft='off')
        
        bin=0.25        
        Mass_arr=np.arange(xlim[0],xlim[1],bin)
        SatFraction=np.zeros(len(Mass_arr),dtype=np.float32)
        HaloSatFraction=np.zeros(len(Mass_arr),dtype=np.float32)
        
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]   
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[i_z])
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
    
        if(i_z==0):
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
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    #pdf.savefig()
    plt.close()
   
    
    xlim=[8.5,11.5]
    ylim=[-7., 1.] 
    bin=0.25 
    
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
      
    linestyle=['-', '--', '-.', ':']    
        
    for i_z in range (0,len(ThisRedshiftList)-1):
           
        char_redshift="%0.2f" % ThisRedshiftList[i_z]
            
        subplot=plt.subplot()
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        
        xlab='$\mathrm{log_{10}}(M_*/$'+Msun+'$)$'      
        ylab=ylab='$\mathrm{log_{10}}(\phi [\mathrm{Mpc^{-3}} \mathrm{log_{10}}(M^{-1})])$'      
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)
        
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
        subplot.yaxis.set_minor_locator(MultipleLocator(0.1))
                
        #MR
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]   
        G0_MR=G0_MR[G0_MR['StellarMass']>0.]
        Mass_MR=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[i_z])
      
        #MRII
        if(MRII==1):
            (sel)=select_current_redshift(G_MRII, ThisRedshiftList, i_z, FullSnapshotList_MRII)
        
            G0_MRII=G_MRII[sel]   
            G0_MRII=G0_MRII[G0_MRII['StellarMass']>0.]
            Mass_MRII=stellar_mass_with_err(G0_MRII, Hubble_h, ThisRedshiftList[i_z])
            
        if(MRII==0):
            (x_axis,y_axis) = plot_mass_function(subplot, Mass_MR, Volume_MR, xlim, bin, MRII, 
                                                 linestyle=linestyle[i_z], color='red')
        else:
            (x_axis,y_axis) = plot_mass_function(subplot, Mass_MR, Volume_MR, xlim, bin, MRII, 
                                                 Mass_MRII=Mass_MRII, Volume_MRII=Volume_MRII, 
                                                 linestyle=linestyle[i_z], color='red')
    
        
        if(i_z==2):
            x_axis_z10=x_axis; y_axis_z10=y_axis
        if(i_z==3):
            x_axis_z20=x_axis; y_axis_z20=y_axis
    
        #labels
        if(i_z==0):
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.85, 
                        color='black', xlog=0, ylog=0, label='z=[0.0,0.4,1.0,2.0]', 
                        fontsize=14, fontweight='normal') 
    
    plt.tight_layout()
    plt.savefig('./fig/plots_sat_MF.pdf')
    #pdf.savefig()
    #plt.close()
    
    
    
    
    xlim=[8.5,11.5]
    ylim=[0.0, 4.] 
       
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
         
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        
    xlab='$\mathrm{log_{10}}(M_*/$'+Msun+'$)$'      
    ylab=ylab='$\phi_{z=1}/\phi_{z=2}$'      
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)
        
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
    subplot.yaxis.set_minor_locator(MultipleLocator(0.1))
     
    subplot.plot(x_axis_z10,10**y_axis_z10/10**y_axis_z20, color='red', linestyle='-')
       
    plt.tight_layout()
    plt.savefig('./fig/plots_sat_MF_dif.pdf')
    #pdf.savefig()
    #plt.close()
 
    
    
    
    
    
    
    if(sat_numbers==1):
    
        xlim=[12.5,16.0]
        ylim=[0., 3.] 
    
        fig = plt.figure(figsize=(15,4))
        grid = gridspec.GridSpec(1, 5)
        grid.update(wspace=0.0, hspace=0.0)
    
        for ii in range (0,len(ThisRedshiftList)):
           
            char_redshift="%0.2f" % ThisRedshiftList[ii]
            
            subplot=plt.subplot(grid[ii])
            subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        
            xlab='$\log_{10}(M_{\mathrm{vir}}/$'+Msun+'$)$'      
            if ii==0:
                ylab='$\log_{10}$(Satellite Number)'
            else:
                ylab=''      
            subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)
                
            if ii>0:
                plt.tick_params(axis='y', which='both', left='on', labelleft='off')
        
            subplot.xaxis.set_major_locator(MultipleLocator(1))    
            subplot.xaxis.set_minor_locator(MultipleLocator(0.25)) 
        
            bin=0.25      
            Mass_arr=np.arange(xlim[0],xlim[1],bin)
            SatNumber=np.zeros(len(Mass_arr),dtype=np.float32)
             
        
            (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
            G0_MR=G_MR[sel]   
            Mvir=np.log10(G0_MR['CentralMvir']*1e10/Hubble_h)
            Type=G0_MR['Type']
    
    
            #Mass > 9.0    
            StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
            G0_MR_sat=G0_MR[StellarMass>9.0]       
            Mvir_sat=np.log10(G0_MR_sat['CentralMvir']*1e10/Hubble_h)
            Type_sat=G0_MR_sat['Type']
        
            for ll in range(0,len(Mass_arr)):
                sel_sat=G0_MR_sat[(Type_sat>0) & (Mvir_sat>Mass_arr[ll]-bin/2.) & (Mvir_sat<Mass_arr[ll]+bin/2.)] 
                sel_centrals=G0_MR[(Type==0) & (Mvir>Mass_arr[ll]-bin/2.) & (Mvir<Mass_arr[ll]+bin/2.)]
                if len(sel_centrals)>0.:
                    SatNumber[ll]=float(len(sel_sat))/float(len(sel_centrals))               
                else:               
                    SatNumber[ll]=0.                                      
            subplot.plot(Mass_arr, np.log10(SatNumber), color='red', linestyle='-', linewidth=2) 
                
            
            #Mass > 10.0    
            StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
            G0_MR_sat=G0_MR[StellarMass>10.0]       
            Mvir_sat=np.log10(G0_MR_sat['CentralMvir']*1e10/Hubble_h)
            Type_sat=G0_MR_sat['Type']
        
            for ll in range(0,len(Mass_arr)):
                sel_sat=G0_MR_sat[(Type_sat>0) & (Mvir_sat>Mass_arr[ll]-bin/2.) & (Mvir_sat<Mass_arr[ll]+bin/2.)]
                sel_centrals=G0_MR[(Type==0) & (Mvir>Mass_arr[ll]-bin/2.) & (Mvir<Mass_arr[ll]+bin/2.)]
                if len(sel_centrals)>0.:
                    SatNumber[ll]=float(len(sel_sat))/float(len(sel_centrals))               
                else:               
                    SatNumber[ll]=0.                                      
            subplot.plot(Mass_arr, np.log10(SatNumber), color='red', linestyle='--', linewidth=2) 
        
        
        
        
            #MagK < 27.0          
            appMag=abs_to_app_mag(G0_MR['MagDust'][:,10], ThisRedshiftList[ii], Hubble_h, 0.315)
            G0_MR_sat=G0_MR[appMag<27.0]       
            Mvir_sat=np.log10(G0_MR_sat['CentralMvir']*1e10/Hubble_h)
            Type_sat=G0_MR_sat['Type']
        
            for ll in range(0,len(Mass_arr)):
                sel_sat=G0_MR_sat[(Type_sat>0) & (Mvir_sat>Mass_arr[ll]-bin/2.) & (Mvir_sat<Mass_arr[ll]+bin/2.)] 
                sel_centrals=G0_MR[(Type==0) & (Mvir>Mass_arr[ll]-bin/2.) & (Mvir<Mass_arr[ll]+bin/2.)]
                if len(sel_centrals)>0.:
                    SatNumber[ll]=float(len(sel_sat))/float(len(sel_centrals))               
                else:               
                    SatNumber[ll]=0.                                      
            subplot.plot(Mass_arr, np.log10(SatNumber), color='green', linestyle='-', linewidth=2) 
        
            #MagK < 25.0          
            appMag=abs_to_app_mag(G0_MR['MagDust'][:,10], ThisRedshiftList[ii], Hubble_h, 0.315)
            G0_MR_sat=G0_MR[appMag<25.0]       
            Mvir_sat=np.log10(G0_MR_sat['CentralMvir']*1e10/Hubble_h)
            Type_sat=G0_MR_sat['Type']
        
            for ll in range(0,len(Mass_arr)):
                sel_sat=G0_MR_sat[(Type_sat>0) & (Mvir_sat>Mass_arr[ll]-bin/2.) & (Mvir_sat<Mass_arr[ll]+bin/2.)] 
                sel_centrals=G0_MR[(Type==0) & (Mvir>Mass_arr[ll]-bin/2.) & (Mvir<Mass_arr[ll]+bin/2.)]
                if len(sel_centrals)>0.:
                    SatNumber[ll]=float(len(sel_sat))/float(len(sel_centrals))               
                else:               
                    SatNumber[ll]=0.                                      
            subplot.plot(Mass_arr, np.log10(SatNumber), color='green', linestyle='--', linewidth=2) 
        
            #MagK < 23.0          
            appMag=abs_to_app_mag(G0_MR['MagDust'][:,10], ThisRedshiftList[ii], Hubble_h, 0.315)
            G0_MR_sat=G0_MR[appMag<23.0]       
            Mvir_sat=np.log10(G0_MR_sat['CentralMvir']*1e10/Hubble_h)
            Type_sat=G0_MR_sat['Type']
        
            for ll in range(0,len(Mass_arr)):
                sel_sat=G0_MR_sat[(Type_sat>0) & (Mvir_sat>Mass_arr[ll]-bin/2.) & (Mvir_sat<Mass_arr[ll]+bin/2.)] 
                sel_centrals=G0_MR[(Type==0) & (Mvir>Mass_arr[ll]-bin/2.) & (Mvir<Mass_arr[ll]+bin/2.)]
                if len(sel_centrals)>0.:
                    SatNumber[ll]=float(len(sel_sat))/float(len(sel_centrals))               
                else:               
                    SatNumber[ll]=0.                                      
            subplot.plot(Mass_arr, np.log10(SatNumber), color='green', linestyle=':', linewidth=2) 
        
            #labels
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.7, y_percentage=0.55, 
                        color='black', xlog=0, ylog=0, label='z='+char_redshift[:-1], 
                        fontsize=14, fontweight='normal') 
        
            if(ii==4):
                plot_label (subplot, 'label', xlim, ylim, 
                            x_percentage=0.15, y_percentage=0.90, color='black', xlog=0, ylog=0, 
                            label='sat mass>10.0', fontsize=12, fontweight='normal') 
                plot_label (subplot, 'line', xlim, ylim,
                            x_percentage=0.04, y_percentage=0.92, color='red', x2_percentage=0.13, 
                            xlog=0, ylog=0, linestyle='-', linewidth=2)
                plot_label (subplot, 'label', xlim, ylim, 
                            x_percentage=0.15, y_percentage=0.84, color='black', xlog=0, ylog=0, 
                            label='sat mass>9.0', fontsize=12, fontweight='normal') 
                plot_label (subplot, 'line', xlim, ylim,
                            x_percentage=0.04, y_percentage=0.86, color='red', x2_percentage=0.13, 
                            xlog=0, ylog=0, linestyle='--', linewidth=2)

                plot_label (subplot, 'label', xlim, ylim, 
                            x_percentage=0.15, y_percentage=0.78, color='black', xlog=0, ylog=0, 
                            label='sat K<27', fontsize=12, fontweight='normal') 
                plot_label (subplot, 'line', xlim, ylim,
                            x_percentage=0.04, y_percentage=0.80, color='green', x2_percentage=0.13, 
                            xlog=0, ylog=0, linestyle='-', linewidth=2)
                plot_label (subplot, 'label', xlim, ylim, 
                            x_percentage=0.15, y_percentage=0.72, color='black', xlog=0, ylog=0, 
                            label='sat K<25', fontsize=12, fontweight='normal') 
                plot_label (subplot, 'line', xlim, ylim,
                            x_percentage=0.04, y_percentage=0.74, color='green', x2_percentage=0.13, 
                            xlog=0, ylog=0, linestyle='--', linewidth=2)
                plot_label (subplot, 'label', xlim, ylim, 
                            x_percentage=0.15, y_percentage=0.66, color='black', xlog=0, ylog=0, 
                            label='sat K<23', fontsize=12, fontweight='normal') 
                plot_label (subplot, 'line', xlim, ylim,
                            x_percentage=0.04, y_percentage=0.68, color='green', x2_percentage=0.13, 
                            xlog=0, ylog=0, linestyle=':', linewidth=2)
        
        plt.tight_layout()         
        plt.savefig('./fig/plots_sat_numbers.pdf')
        plt.close()
       
    
    
    
    return   
#end sat_fraction



def number_density_massive_gals(ThisRedshiftList):
   
    '''xlim=[10.,12.5]
    ylim=[-13.0, -7.0]
    bin=0.25

    fig = plt.figure(figsize=(two_two_size_small[0],two_two_size_small[1]))
    grid = gridspec.GridSpec(2, 2)
    grid.update(wspace=0.0, hspace=0.0)
    
    for i_z in range(0,len(ThisRedshiftList)):
        
        subplot=plt.subplot(grid[i_z])
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        if i_z==2 or i_z == 3: 
            xlab='Mass'
            subplot.set_xlabel(xlab, fontsize=14)      
        if i_z==0 or i_z == 2:
            ylab='SSFR'
            subplot.set_ylabel(ylab, fontsize=14)                       
    
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25)) 
        subplot.yaxis.set_major_locator(MultipleLocator(1))  
        subplot.yaxis.set_minor_locator(MultipleLocator(0.25))
  
        if i_z==1 or i_z==3:
            plt.tick_params(axis='y', which='both', left='on', labelleft='off')
        if i_z==0 or i_z==1:    
            plt.tick_params(axis='x', which='both', bottom='on', labelbottom='off')
            
        #MODEL
        #MR
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]   
        G0_MR=G0_MR[G0_MR['StellarMass']>0.]
        Mass=G0_MR['StellarMass']*1e10/Hubble_h   
        plt.scatter(np.log10(Mass), np.log10(G0_MR['Sfr']/Mass), s=5, color='black') 
                  
    plt.tight_layout()   
    pdf.savefig()
    plt.close()'''
    
    
    
    xlim=[10.5,12.5]
    ylim=[-8.5, -2.5]
    bin=0.25

    fig = plt.figure(figsize=(two_two_size_small[0],two_two_size_small[1]))
    grid = gridspec.GridSpec(2, 2)
    grid.update(wspace=0.0, hspace=0.0)
    
    for i_z in range(0,len(ThisRedshiftList)):
        
        char_redshift="%0.2f" % ThisRedshiftList[i_z]
        
        subplot=plt.subplot(grid[i_z])

        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        if i_z==2 or i_z == 3: 
            xlab='$\mathrm{log_{10}}(M_*/$'+Msun+'$)$'
            subplot.set_xlabel(xlab, fontsize=14)
      
        if i_z==0 or i_z == 2:
            ylab='$\mathrm{log_{10}}(\phi [\mathrm{Mpc^{-3}}])$'
            subplot.set_ylabel(ylab, fontsize=14)
                       
    
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25)) 
        subplot.yaxis.set_major_locator(MultipleLocator(1))  
        subplot.yaxis.set_minor_locator(MultipleLocator(0.25))
  
        if i_z==1 or i_z==3:
            plt.tick_params(axis='y', which='both', left='on', labelleft='off')
        if i_z==0 or i_z==1:    
            plt.tick_params(axis='x', which='both', bottom='on', labelbottom='off')
        
    
        #MODEL
        #MR
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)
        
        G0_MR=G_MR[sel]   
        G0_MR=G0_MR[G0_MR['StellarMass']>0.]
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[i_z])
        #StellarMass=np.log10((G0_MR['MassFromInSitu']+G0_MR['MassFromMergers']+G0_MR['MassFromBursts'])*1.e10*Hubble_h)
        
        bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
        hist_MR=np.histogram(StellarMass, bins=bin_arr, range=(xlim[0],xlim[1]))   
       
        x_axis=hist_MR[1][0:len(hist_MR[1][:])-1]+bin/2.           
        hist_MR=hist_MR[0]       
        y_axis=np.log10(hist_MR/(Volume_MR*bin))
        subplot.plot(x_axis,y_axis, color='black', linewidth=2, linestyle='-') 
         
        #small_error    
        np.random.seed(seed=10)
        StellarMass=(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)+
                     np.random.randn(len(G0_MR['StellarMass']))*0.04*(1+ThisRedshiftList[i_z]))
        bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
        hist_MR=np.histogram(StellarMass, bins=bin_arr, range=(xlim[0],xlim[1]))   
       
        x_axis=hist_MR[1][0:len(hist_MR[1][:])-1]+bin/2.           
        hist_MR=hist_MR[0]       
        y_axis=np.log10(hist_MR/(Volume_MR*bin))
        subplot.plot(x_axis,y_axis, color='black', linewidth=2, linestyle='--') 
        
        
        
        #MODEL SF
        #MR
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)
        
        G0_MR=G_MR[sel]  
        if(i_z<2):
            G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1e10/Hubble_h))>-10.)]
        else:
            G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1e10/Hubble_h))>-9.)]
            
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[i_z])
        #StellarMass=np.log10((G0_MR['MassFromInSitu']+G0_MR['MassFromMergers']+G0_MR['MassFromBursts'])*1.e10*Hubble_h)
        
        bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
        hist_MR=np.histogram(StellarMass, bins=bin_arr, range=(xlim[0],xlim[1]))   
       
        x_axis=hist_MR[1][0:len(hist_MR[1][:])-1]+bin/2.           
        hist_MR=hist_MR[0]       
        y_axis=np.log10(hist_MR/(Volume_MR*bin))
        subplot.plot(x_axis,y_axis, color='blue', linewidth=2, linestyle='-') 
         
        #small_error    
        np.random.seed(seed=10)
        StellarMass=(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)+
                     np.random.randn(len(G0_MR['StellarMass']))*0.04*(1+ThisRedshiftList[i_z]))
        bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
        hist_MR=np.histogram(StellarMass, bins=bin_arr, range=(xlim[0],xlim[1]))   
       
        x_axis=hist_MR[1][0:len(hist_MR[1][:])-1]+bin/2.           
        hist_MR=hist_MR[0]       
        y_axis=np.log10(hist_MR/(Volume_MR*bin))
        subplot.plot(x_axis,y_axis, color='blue', linewidth=2, linestyle='--') 
        
        
        
        #MODEL Quenched
        #MR
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)
        
        G0_MR=G_MR[sel]  
        if(i_z<2):
            G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1e10/Hubble_h))<-10.)]
        else:
            G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1e10/Hubble_h))<-9.)]
            
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[i_z])
        #StellarMass=np.log10((G0_MR['MassFromInSitu']+G0_MR['MassFromMergers']+G0_MR['MassFromBursts'])*1.e10*Hubble_h)
        
        bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
        hist_MR=np.histogram(StellarMass, bins=bin_arr, range=(xlim[0],xlim[1]))   
       
        x_axis=hist_MR[1][0:len(hist_MR[1][:])-1]+bin/2.           
        hist_MR=hist_MR[0]       
        y_axis=np.log10(hist_MR/(Volume_MR*bin))
        subplot.plot(x_axis,y_axis, color='red', linewidth=2, linestyle='-') 
         
        #small_error    
        np.random.seed(seed=10)
        StellarMass=(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)+
                     np.random.randn(len(G0_MR['StellarMass']))*0.04*(1+ThisRedshiftList[i_z]))
        bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
        hist_MR=np.histogram(StellarMass, bins=bin_arr, range=(xlim[0],xlim[1]))   
       
        x_axis=hist_MR[1][0:len(hist_MR[1][:])-1]+bin/2.           
        hist_MR=hist_MR[0]       
        y_axis=np.log10(hist_MR/(Volume_MR*bin))
        subplot.plot(x_axis,y_axis, color='red', linewidth=2, linestyle='--') 
        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.7, y_percentage=0.55, 
                    color='black', xlog=0, ylog=0, label='z='+char_redshift[:-1], 
                    fontsize=14, fontweight='normal') 
        
        #LABELS 
        if i_z==3:
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.15, y_percentage=0.90, 
                        color='black', xlog=0, ylog=0, label='total', fontsize=12, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, y_percentage=0.92, 
                        color='black', x2_percentage=0.13, xlog=0, ylog=0, linestyle='-', linewidth=2)
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.15, y_percentage=0.82, 
                        color='black', xlog=0, ylog=0, label='SF', fontsize=12, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, y_percentage=0.84, 
                        color='blue', x2_percentage=0.13, xlog=0, ylog=0, linestyle='-', linewidth=2)
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.15, y_percentage=0.76, 
                        color='black', xlog=0, ylog=0, label='quenched', fontsize=12, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.04, y_percentage=0.78, 
                        color='red', x2_percentage=0.13, xlog=0, ylog=0, linestyle='-', linewidth=2)
            
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()
   
    return 
#end    
    
    
    
    
def BHmass_in_radio(ThisRedshiftList):
    
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
           
    xlab='$\mathrm{log_{10}}(M_*[h^{-2}$'+Msun+'$])$'   
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
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()

    return    
#end BHmass_in_radio



def sfr_massive_galaxies():   

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
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()

    return    
#end sfr_massive_galaxies


def HotGas_fraction(ThisRedshiftList):
    
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
           
    xlab='$\mathrm{log_{10}}(M_vir[h^{-1}$'+Msun+'$])$'   
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
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()

    return 
#end HotGas_fraction









def fabian_fb(ThisRedshiftList):     
        
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
    xlab='$\mathrm{log_{10}}(M_*/$'+Msun+'$)$'      
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


    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()
    
    return 
#endif fabian_fb








def misc_plots(FullSnapshotList):
    
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
    
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()
 
    return 
#end misc_plots
     
def surface_density_vs_stellar_mass(ThisRedshiftList):
    
    
    #xlim=[-15.0,-9.0]
    #ylim=[5.0, 11.0]
    xlim=[9.0,12.5]
    #ylim=[5.0, 12.0]
    #ylim=[-13.0, -9.0]
    ylim=[-5.0, 5.0]
   
    bin=0.25
    
    ThisRedshiftList=[0.0]
    
    color=['black','grey','violet','aqua','green','lightgreen','yellow','orange','pink', 'red','darkred', 'brown']
    
    fig = plt.figure(figsize=(two_two_size_large[0],two_two_size_large[1]))    
    grid = gridspec.GridSpec(2, 2)
    grid.update(wspace=0.0, hspace=0.0)
    
    subplot=plt.subplot(grid[0])
    
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    xlab='$\mathrm{log_{10}}(M_*/$'+Msun+'$)$'
    #xlab='$\mathrm{log_{10}}(sSFR[yr^{-1}])$'
    subplot.set_xlabel(xlab, fontsize=14)
     
    #ylab='$\mathrm{log_{10}}(\Sigma_{\mathrm{H}_2} [$'+Msun+'$ \mathrm{Kpc}^{-2}])$'
    #ylab='$\mathrm{log_{10}}(sSFR[yr^{-1}])$'
    #ylab='$\mathrm{log_{10}}(SFR[$'+Msun+'$yr^{-1}])$'
    ylab='$\mathrm{log_{10}}(\Sigma_{SFR} [$'+Msun+'$ \mathrm{Kpc}^{-2}])$'
    subplot.set_ylabel(ylab, fontsize=14)
    
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
    subplot.yaxis.set_minor_locator(MultipleLocator(0.25))
    
    area = np.zeros(RNUM, dtype=np.float32) 
    for ii in range(0, RNUM):
        if ii==0:
            area[ii] = 3.14 * RingRadius[ii]*RingRadius[ii]
        else:
            area[ii] = 3.14 * (RingRadius[ii]*RingRadius[ii]-RingRadius[ii-1]*RingRadius[ii-1])
              
    for i_z in range(0,len(ThisRedshiftList)):
        
        char_redshift="%0.2f" % ThisRedshiftList[i_z]
        
        #all galaxies
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)        
        G0_MR=G_MR[sel] 
        StellarMass = np.log10(G0_MR['StellarMass']/Hubble_h*1.e10)
        G0_MR = G0_MR[(StellarMass>9.0)]
                
        total_ngals = len(G0_MR)       
        Ngals=2000
        if(Ngals>total_ngals):
            Ngals=total_ngals
        G0_MR=np.random.choice(G0_MR, size=Ngals, replace=False) 
        #print("Ngals = %d" % len(G0_MR))
            
        StellarMass = np.log10(G0_MR['StellarMass']/Hubble_h*1.e10)
        GasMass =  np.log10(G0_MR['ColdGas']/Hubble_h*1.e10)
        SSFR = np.log10(G0_MR['Sfr']/10**StellarMass)
        
        subplot.scatter(StellarMass,SSFR, s=1, color='black')
        
        #SurfaceDensity = np.log10((G0_MR['ColdGas']*G0_MR['H2fraction']/Hubble_h*1e10)/G0_MR['GasDiskRadius'])
        SurfaceDensity = np.log10((G0_MR['ColdGasRings'][:,0]/Hubble_h*1e10) / area[0])
               
        for ii in range(0,RNUM):
            SurfaceDensity = np.log10((G0_MR['ColdGasRings'][:,ii]*G0_MR['H2fractionRings'][:,ii]/Hubble_h*1e10) / area[ii]) 
            SfrRings = np.log10(G0_MR['SfrRings'][:,ii]/(G0_MR['StellarMass']/Hubble_h*1.e10))
            #SfrRings = np.log10(G0_MR['SfrRings'][:,ii])
            #subplot.scatter(StellarMass,SurfaceDensity, s=5, color=color[ii])
            #subplot.scatter(StellarMass,SfrRings, s=5, color=color[ii])
            #subplot.scatter(SSFR,SurfaceDensity, s=5, color=color[ii])
        
        #subplot.scatter(StellarMass, SurfaceDensity, s=5)
        #subplot.scatter(SSFR,GasMass, s=5, color='blue')
        
        
        #massive galaxies
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)        
        G0_MR=G_MR[sel] 
        StellarMass = np.log10(G0_MR['StellarMass']/Hubble_h*1.e10)
        G0_MR = G0_MR[(StellarMass>11.0) & (StellarMass<11.5)]
                
        total_ngals = len(G0_MR)       
        Ngals=1000
        if(Ngals>total_ngals):
            Ngals=total_ngals
        G0_MR=np.random.choice(G0_MR, size=Ngals, replace=False) 
        #print("Ngals = %d" % len(G0_MR))
            
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[i_z])
        SSFR = np.log10(G0_MR['Sfr']/10**StellarMass)
        #SurfaceDensity = np.log10((G0_MR['ColdGas']*G0_MR['H2fraction']/Hubble_h*1e10)/G0_MR['GasDiskRadius'])
        SurfaceDensity = np.log10((G0_MR['ColdGasRings'][:,0]*G0_MR['H2fractionRings'][:,0]/Hubble_h*1e10) / area[0])
        
       
        for ii in range(0,RNUM):
            SurfaceDensity = np.log10((G0_MR['ColdGasRings'][:,ii]*G0_MR['H2fractionRings'][:,ii]/Hubble_h*1e10) / area[ii]) 
            sSfrRings = np.log10(G0_MR['SfrRings'][:,ii]/(G0_MR['StellarMass']/Hubble_h*1.e10))
            SfrRings = np.log10(G0_MR['SfrRings'][:,ii])
            SFRDensity = np.log10(G0_MR['SfrRings'][:,ii] / area[ii]) 
            #subplot.scatter(SSFR,SurfaceDensity, s=5, color=color[ii])
            ##subplot.scatter(sSfrRings,SurfaceDensity, s=5, color=color[ii])
            
            #subplot.scatter(StellarMass,SurfaceDensity, s=5, color=color[ii])
            #subplot.scatter(StellarMass,sSfrRings, s=5, color=color[ii])
            subplot.scatter(StellarMass,SfrRings, s=5, color=color[ii])
            #subplot.scatter(StellarMass,SFRDensity, s=5, color=color[ii])
            
                        
        #subplot.scatter(StellarMass, SurfaceDensity, s=5)
        #subplot.scatter(SSFR,SurfaceDensity, s=5, color='red')
        
        
        
         #10^10 galaxies
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)        
        G0_MR=G_MR[sel] 
        StellarMass = np.log10(G0_MR['StellarMass']/Hubble_h*1.e10)
        G0_MR = G0_MR[(StellarMass>10.0) & (StellarMass<10.5)]
                
        total_ngals = len(G0_MR)       
        Ngals=1000
        if(Ngals>total_ngals):
            Ngals=total_ngals
        G0_MR=np.random.choice(G0_MR, size=Ngals, replace=False) 
        #print("Ngals = %d" % len(G0_MR))
            
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[i_z])
        SSFR = np.log10(G0_MR['Sfr']/10**StellarMass)
        #SurfaceDensity = np.log10((G0_MR['ColdGas']*G0_MR['H2fraction']/Hubble_h*1e10)/G0_MR['GasDiskRadius'])
        SurfaceDensity = np.log10((G0_MR['ColdGasRings'][:,0]*G0_MR['H2fractionRings'][:,0]/Hubble_h*1e10) / area[0])
        
        color=['black','grey','violet','aqua','green','lightgreen','yellow','orange','pink', 'red','darkred', 'brown']
        for ii in range(0,RNUM):
            SurfaceDensity = np.log10((G0_MR['ColdGasRings'][:,ii]*G0_MR['H2fractionRings'][:,ii]/Hubble_h*1e10) / area[ii]) 
            SfrRings = np.log10(G0_MR['SfrRings'][:,ii]/(G0_MR['StellarMass']/Hubble_h*1.e10))
            #subplot.scatter(SSFR,SurfaceDensity, s=5, color=color[ii])
            ##subplot.scatter(SfrRings,SurfaceDensity, s=5, color=color[ii])
                        
        #subplot.scatter(StellarMass, SurfaceDensity, s=5)
        #subplot.scatter(SSFR,SurfaceDensity, s=5, color='fuchsia')
       
      
    
    plt.tight_layout()    
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()
    
    
    
   
  
    return
#end surface_density_vs_stellar_mass

def test_resolution(ThisRedshiftList):
   
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
    plt.tick_params(axis='x', which='both', top=True, labeltop=True, bottom=False, labelbottom=False)
    
    xlab='$\log_{10}(M_{\mathrm{Cold}}/$'+Msun+'$)$'   
    ylab='$\log_{10}(\phi /(\mathrm{Mpc^{-3}} \mathrm{dex}))$'   
    subplot.set_xlabel(xlab, fontsize=15), subplot.set_ylabel(ylab, fontsize=15)
    subplot.xaxis.set_label_position('top')  
        
    #MR  
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, 0, FullSnapshotList_MR)         
    G0_MR=G_MR[sel]   
    G0_MR=G0_MR[(G0_MR['ColdGas']>0.) & (G0_MR['Type']==0)]
    mass=(np.log10(G0_MR['ColdGas']*1.e10/Hubble_h))

    bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
    hist=np.histogram(mass, bins=bin_arr, range=(xlim[0],xlim[1]))   
    xx = hist[1][0:len(hist[1][:])-1]+bin/2.
    yy = np.log10(hist[0][:]/(Volume_MR*bin))
    subplot.plot(xx,yy,color='red', linewidth=2, linestyle='--')

    #WRITE OUTPUT      
    if(write_to_file==1):
        df = pd.DataFrame({'log10_M':xx, 'log10_phi':yy})
        file = Datadir+file_to_write+'Resolution_ColdGasMF_MR'+str(f'_z{ThisRedshiftList[0]:0.2f}')+'.csv'
        df.to_csv(file,index=False)
        #df = pd.read_csv(file)
        #subplot.plot(df['log10_M'],df['log10_phi'], color='black')
        
    #MRII
    (sel)=select_current_redshift(G_MRII, ThisRedshiftList, 0, FullSnapshotList_MRII)         
    G0_MRII=G_MRII[sel]   
    G0_MRII=G0_MRII[(G0_MRII['ColdGas']>0.) & (G0_MRII['Type']==0)]
    mass=(np.log10(G0_MRII['ColdGas']*1.e10/Hubble_h))

    bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
    hist=np.histogram(mass, bins=bin_arr, range=(xlim[0],xlim[1]))   
    xx = hist[1][0:len(hist[1][:])-1]+bin/2.
    yy = np.log10(hist[0][:]/(Volume_MRII*bin))
    subplot.plot(xx,yy,color='red', linewidth=2, linestyle='-')
    
    #WRITE OUTPUT      
    if(write_to_file==1):
        df = pd.DataFrame({'log10_M':xx, 'log10_phi':yy})
        file = Datadir+file_to_write+'Resolution_ColdGasMF_MRII'+str(f'_z{ThisRedshiftList[0]:0.2f}')+'.csv'
        df.to_csv(file,index=False)
        #df = pd.read_csv(file)
        #subplot.plot(df['log10_M'],df['log10_phi'], color='black')
        
    
    plot_label (subplot, 'label', xlim, ylim, 
                x_percentage=0.15, y_percentage=0.15, color='black', xlog=0, ylog=0, 
                label='ColdGas', fontsize=12, fontweight='normal')  
    
    
        
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
    plt.tick_params(axis='x', which='both', top=True, labeltop=True, bottom=False, labelbottom=False)   
    plt.tick_params(axis='y', which='both', left=True, labelleft=False)
    
    xlab='$\log_{10}(\mathrm{SFR}/($'+Msun+'$yr^{-1}))$'  
    subplot.set_xlabel(xlab, fontsize=15)
    subplot.xaxis.set_label_position('top')  
    
    #MR
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, 0, FullSnapshotList_MR)         
    G0_MR=G_MR[sel]   
    G0_MR=G0_MR[G0_MR['Sfr']>0.]
    mass=(np.log10(G0_MR['Sfr']/Hubble_h))

    bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
    hist=np.histogram(mass, bins=bin_arr, range=(xlim[0],xlim[1]))   
    xx = hist[1][0:len(hist[1][:])-1]+bin/2.
    yy = np.log10(hist[0][:]/(Volume_MR*bin))
    subplot.plot(xx,yy,color='red', linewidth=2, linestyle='--')

    #WRITE OUTPUT      
    if(write_to_file==1):
        df = pd.DataFrame({'log10_M':xx, 'log10_phi':yy})
        file = Datadir+file_to_write+'Resolution_SFRF_MR'+str(f'_z{ThisRedshiftList[0]:0.2f}')+'.csv'
        df.to_csv(file,index=False)
        #df = pd.read_csv(file)
        #subplot.plot(df['log10_M'],df['log10_phi'], color='black')

    #MRII
    (sel)=select_current_redshift(G_MRII, ThisRedshiftList, 0, FullSnapshotList_MRII)         
    G0_MRII=G_MRII[sel]   
    G0_MRII=G0_MRII[G0_MRII['Sfr']>0.]
    mass=(np.log10(G0_MRII['Sfr']/Hubble_h))

    bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
    hist=np.histogram(mass, bins=bin_arr, range=(xlim[0],xlim[1]))   
    xx = hist[1][0:len(hist[1][:])-1]+bin/2.
    yy = np.log10(hist[0][:]/(Volume_MRII*bin))
    subplot.plot(xx,yy,color='red', linewidth=2, linestyle='-')
    
    #WRITE OUTPUT      
    if(write_to_file==1):
        df = pd.DataFrame({'log10_M':xx, 'log10_phi':yy})
        file = Datadir+file_to_write+'Resolution_SFRF_MRII'+str(f'_z{ThisRedshiftList[0]:0.2f}')+'.csv'
        df.to_csv(file,index=False)
        #df = pd.read_csv(file)
        #subplot.plot(df['log10_M'],df['log10_phi'], color='black')
    
    plot_label (subplot, 'label', xlim, ylim, 
                x_percentage=0.15, y_percentage=0.15, color='black', xlog=0, ylog=0, 
                label='SFR', fontsize=12, fontweight='normal') 
    
    
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
    plt.tick_params(axis='x', which='both', top=False, labeltop=False)  
    
    xlab='$\log_{10}(M_*/$'+Msun+'$)$'   
    ylab='$\log_{10}(\phi /(\mathrm{Mpc^{-3}} \mathrm{dex}))$'   
    subplot.set_xlabel(xlab, fontsize=15), subplot.set_ylabel(ylab, fontsize=15)
   
    
    #MR
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, 0, FullSnapshotList_MR)         
    G0_MR=G_MR[sel]   
    G0_MR=G0_MR[G0_MR['StellarMass']>0.]
    mass=(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h))

    bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
    hist=np.histogram(mass, bins=bin_arr, range=(xlim[0],xlim[1]))   
    xx = hist[1][0:len(hist[1][:])-1]+bin/2.
    yy = np.log10(hist[0][:]/(Volume_MR*bin))
    subplot.plot(xx,yy,color='red', linewidth=2, linestyle='--')

    #WRITE OUTPUT      
    if(write_to_file==1):
        df = pd.DataFrame({'log10_M':xx, 'log10_phi':yy})
        file = Datadir+file_to_write+'Resolution_SMF_MR'+str(f'_z{ThisRedshiftList[0]:0.2f}')+'.csv'
        df.to_csv(file,index=False)
        #df = pd.read_csv(file)
        #subplot.plot(df['log10_M'],df['log10_phi'], color='black')

    #MRII
    (sel)=select_current_redshift(G_MRII, ThisRedshiftList, 0, FullSnapshotList_MRII)         
    G0_MRII=G_MRII[sel]   
    G0_MRII=G0_MRII[G0_MRII['StellarMass']>0.]
    mass=(np.log10(G0_MRII['StellarMass']*1.e10/Hubble_h))

    bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
    hist=np.histogram(mass, bins=bin_arr, range=(xlim[0],xlim[1]))   
    xx = hist[1][0:len(hist[1][:])-1]+bin/2.
    yy = np.log10(hist[0][:]/(Volume_MRII*bin))
    subplot.plot(xx,yy,color='red', linewidth=2, linestyle='-')
    
    #WRITE OUTPUT      
    if(write_to_file==1):
        df = pd.DataFrame({'log10_M':xx, 'log10_phi':yy})
        file = Datadir+file_to_write+'Resolution_SMF_MRII'+str(f'_z{ThisRedshiftList[0]:0.2f}')+'.csv'
        df.to_csv(file,index=False)
        #df = pd.read_csv(file)
        #subplot.plot(df['log10_M'],df['log10_phi'], color='black')
    
    plot_label (subplot, 'label', xlim, ylim, 
                x_percentage=0.15, y_percentage=0.15, color='black', xlog=0, ylog=0, 
                label='StellarMass', fontsize=12, fontweight='normal') 
           
        
    plot_label (subplot, 'label', xlim, ylim, 
                x_percentage=0.15, y_percentage=0.5, color='black', xlog=0, ylog=0, 
                label='MS', fontsize=13, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim,
                x_percentage=0.04, y_percentage=0.53, color='red', x2_percentage=0.13, 
                xlog=0, ylog=0, linestyle='--', linewidth=2)
    plot_label (subplot, 'label', xlim, ylim, 
                x_percentage=0.15, y_percentage=0.4, color='black', xlog=0, ylog=0, 
                label='MSII', fontsize=13, fontweight='normal') 
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
    plt.tick_params(axis='x', which='both', top=False, labeltop=False) 
    plt.tick_params(axis='y', which='both', left=True, labelleft=False)
    
    xlab='$\log_{10}(M_{\mathrm{BH}}/$'+Msun+'$)$'         
    subplot.set_xlabel(xlab, fontsize=15)
    
          
    #MR  
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, 0, FullSnapshotList_MR)         
    G0_MR=G_MR[sel]   
    G0_MR=G0_MR[(G0_MR['BlackHoleMass']>0.) & (G0_MR['Type']==0)]
    mass=(np.log10(G0_MR['BlackHoleMass']*1.e10/Hubble_h))

    bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
    hist=np.histogram(mass, bins=bin_arr, range=(xlim[0],xlim[1]))   
    xx = hist[1][0:len(hist[1][:])-1]+bin/2.
    yy = np.log10(hist[0][:]/(Volume_MR*bin))
    subplot.plot(xx,yy,color='red', linewidth=2, linestyle='--')

    #WRITE OUTPUT      
    if(write_to_file==1):
        df = pd.DataFrame({'log10_M':xx, 'log10_phi':yy})
        file = Datadir+file_to_write+'Resolution_BHMF_MR'+str(f'_z{ThisRedshiftList[0]:0.2f}')+'.csv'
        df.to_csv(file,index=False)
        #df = pd.read_csv(file)
        #subplot.plot(df['log10_M'],df['log10_phi'], color='black')

    #MRII
    (sel)=select_current_redshift(G_MRII, ThisRedshiftList, 0, FullSnapshotList_MRII)         
    G0_MRII=G_MRII[sel]   
    G0_MRII=G0_MRII[(G0_MRII['BlackHoleMass']>0.) & (G0_MRII['Type']==0)]
    mass=(np.log10(G0_MRII['BlackHoleMass']*1.e10/Hubble_h))

    bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
    hist=np.histogram(mass, bins=bin_arr, range=(xlim[0],xlim[1]))   
    xx = hist[1][0:len(hist[1][:])-1]+bin/2.
    yy = np.log10(hist[0][:]/(Volume_MRII*bin))
    subplot.plot(xx,yy,color='red', linewidth=2, linestyle='-')
    
    #WRITE OUTPUT      
    if(write_to_file==1):
        df = pd.DataFrame({'log10_M':xx, 'log10_phi':yy})
        file = Datadir+file_to_write+'Resolution_BHMF_MRII'+str(f'_z{ThisRedshiftList[0]:0.2f}')+'.csv'
        df.to_csv(file,index=False)
        #df = pd.read_csv(file)
        #subplot.plot(df['log10_M'],df['log10_phi'], color='black')
        
    plot_label (subplot, 'label', xlim, ylim, 
                x_percentage=0.15, y_percentage=0.15, color='black', xlog=0, ylog=0, 
                label='BlackHoleMass', fontsize=12, fontweight='normal') 
    
    
    
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYF19_test_resolution_rings.pdf')
    plt.close()
    
    return 
#end test_resolution_rings
    
    
def H2fraction_fits():
    
    xlim=[1.,4.]    
    ylim=[-1.,0.] 
    
    xx = [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
          1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 
          2., 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 
          3., 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9]
    
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))   
    #subplot=fig.add_subplot(1,1,1)
    subplot=plt.subplot()
    
    
    metallicity = np.array([0.01, 0.03, 0.09, 0.3, 1., 3.0])    
    #metallicity = np.array([0.01, 0.02, 0.04,0.08, 0.16, 0.32, 0.64, 1.3, 2.6])    
    
    SigmaHRings = np.arange(0.0,4.0,0.1)
    print(SigmaHRings)
    #SFR  
   
    bin=0.25
    
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    #majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.25))
    subplot.yaxis.set_major_locator(MultipleLocator(1))    
    subplot.yaxis.set_minor_locator(MultipleLocator(0.25))   
    
    xlab='simgaH'  
    subplot.set_xlabel(xlab, fontsize=14)
    ylab='H2fract'  
    subplot.set_ylabel(ylab, fontsize=14)
    
    for element in metallicity:    
        print(element)
        khi = 3.1/4.1 * (1.+3.1 * (element**0.365) )
        tau = 0.066 * 10**SigmaHRings * element;
        s = np.log(1.+0.6*khi+0.01*khi*khi)/(0.6*tau);
        #H2fraction = np.max(4-2.*s,0.)/(4.+s)
        H2fraction = 1-0.75*s/(1+0.25*s)
        print(H2fraction)
        subplot.plot(SigmaHRings, np.log10(H2fraction))
    
    
    plt.tight_layout()   
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()
    
    return 
#end H2fraction_fits
    
    
    

def SFH_bins_1():
      
    
    xlim=[0,60]
    ylim=[0,14] 
   
    plot_color=['red','purple']        
   
    fig = plt.figure(figsize=(one_two_size_large[0],one_two_size_large[1]))
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis   
    subplot.xaxis.set_major_locator(MultipleLocator(10.0))    
    subplot.xaxis.set_minor_locator(MultipleLocator(1.0))    
    subplot.yaxis.set_major_locator(MultipleLocator(2.0))    
    subplot.yaxis.set_minor_locator(MultipleLocator(2.0)) 
    
    xlab='Snap'       
    ylab='Loockbacktime (Gyr)'         
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)    
   
    subplot.xaxis.set_label_position('top')
    
    for key, spine in subplot.spines.items():
        spine.set_visible(False)
    #subplot.spines['right'].set_visible(False)
    #subplot.spines['top'].set_visible(False)
    
    subplot.tick_params(which='both', right=False, bottom=False)
    subplot.tick_params(which='both', top=True)
    subplot.tick_params(axis='x',labelbottom=False,labeltop=True) 
    
    
    file=Datadir+'Database_snapshots_table.fits'
    fits_table=fits.open(file)
    snap_table = fits_table[1]       
    
    FullRedshiftList=snap_table.data['z'][::-1]    
    FullTimeList=snap_table.data['lookBackTime'][::-1]
    FullSnapList=snap_table.data['SnapNum'][::-1]
    
    RedshiftList=sorted(FullRedshiftList[(FullRedshiftList>0.) & (FullSnapList>0)], reverse=True)           
    TimeList=sorted(FullTimeList[(FullRedshiftList>0.) & (FullSnapList>0)], reverse=True)      
    SnapList=sorted(FullSnapList[(FullRedshiftList>0.) & (FullSnapList>0)])   
    
    
    
    #Read              
    fa = open(DirName_MR+"SFH_Bins","rb")                
    nbins =  np.fromfile(fa,np.int32,1)   
    template = np.dtype([('SnapNum',np.int32,1),
                         ('Bin',np.int32,1),
                         ('Lookbacktime',np.float64,1),                           
                         ('dt',np.float64,1),
                         ('nbins',np.int32,1)
                        ])
    SFH = np.fromfile(fa,template,nbins[0])    
    fa.close()            
    #we only need the SFH strucutre from the current snap  
    SFH=SFH[SFH['SnapNum']<=np.max(SnapList)]
    
    Snap_list = set(SFH['SnapNum'])
    
    plt_color = [plt.cm.GnBu(i) for i in np.linspace(0.3, 0.9, len(Snap_list))]   
    plt_color_edge = [plt.cm.GnBu(i) for i in np.linspace(0.4, 1.0, len(Snap_list))]   
    
    
    label_list = [x for x in range(1, max(Snap_list)+1,8)]
   
    for snap in range(2, max(Snap_list)+1,2):         
        SFH_isnap = SFH[SFH['SnapNum']==snap]        
        for ii in range(0, len(SFH_isnap['SnapNum'])):
            yerr = SFH_isnap['dt'][ii]/1.e9/2.
            y = SFH_isnap['Lookbacktime'][ii]/1.e9+(TimeList[snap-1])
           
        
            subplot.fill_between([snap-0.2, snap+0.2], [y-yerr, y-yerr], [y+yerr, y+yerr], interpolate=True, 
                                 edgecolor=plt_color_edge[snap-1], facecolor=plt_color[snap-1], alpha=0.7) 
            
            y = np.log10(SFH_isnap['Lookbacktime'][ii]) 
            subplot.text(snap+0.3,13.3-ii*0.6,str(f"{y:0.1f}"), fontsize=7, color='black') 
            
        subplot.text(snap+0.3,13.3-(ii+1)*0.6,str("z="+str(f"{RedshiftList[snap-1]:0.2f}  |")), 
                     fontsize=7, color='black', rotation=90) 
        
        
    subplot.text(15,1.0,str("Age Bins in $\log_{10}$(yr)"), fontsize=12, color='black') 
        
    '''if(y>1):
        subplot.text(snap+0.2,14-ii*0.7,str(f"{y:0.1f}"), fontsize=7, color='black')    
    elif(y>0.1):    
        subplot.text(snap+0.2,14-ii*0.7,str(f"{y:0.2f}"), fontsize=7, color='black')  
    else:     
        subplot.text(snap+0.2,14-ii*0.7,str(f"{y:0.1e}"), fontsize=7, color='black')'''
           
                
            
    
    '''inset_xlim=[0,60]
    inset_ylim=[0,14]    
   
    inset = inset_axes(subplot, width="100%", height="100%", 
                       bbox_to_anchor=(0.1, 0.05, 0.48, 0.4), 
                       bbox_transform=subplot.figure.transFigure)
    #fig.subplots_adjust(hspace=0.4)
    inset.set_ylim(inset_ylim), inset.set_xlim(inset_xlim)    
    
    #for key, spine in inset.spines.items():
    #    spine.set_visible(False)   
    inset.tick_params(which='both', right=False, left=False, bottom=False, top=False)
    inset.tick_params(axis='y', which='both', labelleft='off') 
    inset.tick_params(axis='x', which='both', labelbottom='off') 

    for snap in range(1, max(Snap_list)+1,4):         
        SFH_isnap = SFH[SFH['SnapNum']==snap]               
        for ii in range(0, len(SFH_isnap['SnapNum'])):
            yerr = SFH_isnap['dt'][ii]/1.e9/2.
            y = SFH_isnap['Lookbacktime'][ii]/1.e9+(TimeList[snap-1])
            inset.text(snap-1.0,y,SFH_isnap['nbins'][ii], fontsize=10, color='black')'''
    
    #LABELS        
    '''plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.9, 
                color='black', xlog=0, ylog=0, label='Mernier 2018', 
                fontsize=12, fontweight='normal') 
    plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.925, 
                color='steelblue', xlog=0, ylog=0, sym='o', sym_size=3, err_size=0.04)

    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.84, 
                color='black', xlog=0, ylog=0, label='Yates 2017', 
                fontsize=12, fontweight='normal') 
    plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.865, 
                color='darkorange', xlog=0, ylog=0, sym='o', sym_size=3, err_size=0.04)

    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.78, 
                color='black', xlog=0, ylog=0, label=prefix_this_model, 
                fontsize=12, fontweight='normal') 
    plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.805, 
                color='black', xlog=0, ylog=0, sym='o', sym_size=3, err_size=0.001)'''

       
    
    
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYF19_SFH_bins.pdf')
    plt.close()

    return 
#end SFH_bins 

    
    
    
    
   
        
def SFH_bins():
      
    
    xlim=[10,60]
    ylim=[6,10.5] 
   
    plot_color=['red','purple']        
   
    fig = plt.figure(figsize=(one_two_size_large[0],one_two_size_large[1]))
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis   
    subplot.xaxis.set_major_locator(MultipleLocator(10.0))    
    subplot.xaxis.set_minor_locator(MultipleLocator(5.0))    
    subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
    subplot.yaxis.set_minor_locator(MultipleLocator(1.0)) 
    
    xlab='Snap'       
    ylab='$\log_{10}(yr)$'         
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)    
   
    for key, spine in subplot.spines.items():
        spine.set_visible(False)
    
    subplot.tick_params(which='both', top=False)
   
    
    file=Datadir+'Database_snapshots_table.fits'
    fits_table=fits.open(file)
    snap_table = fits_table[1]       
    
    FullRedshiftList=snap_table.data['z'][::-1]    
    FullTimeList=snap_table.data['lookBackTime'][::-1]*1.e10
    FullSnapList=snap_table.data['SnapNum'][::-1]
    
    RedshiftList=sorted(FullRedshiftList[(FullRedshiftList>0.) & (FullSnapList>0)], reverse=True)           
    TimeList=sorted(FullTimeList[(FullRedshiftList>0.) & (FullSnapList>0)], reverse=True)      
    SnapList=sorted(FullSnapList[(FullRedshiftList>0.) & (FullSnapList>0)])   
    
    
    
    #Read   
    fa = open(DirName_MR+"SFH_Bins","rb")                
    nbins =  np.fromfile(fa,np.int32,1)   
    template = np.dtype([('SnapNum',np.int32,1),
                         ('Bin',np.int32,1),
                         ('Lookbacktime',np.float64,1),                           
                         ('dt',np.float64,1),
                         ('nbins',np.int32,1)
                        ])
    SFH = np.fromfile(fa,template,nbins[0])    
    fa.close()            
    #we only need the SFH strucutre from the current snap  
    SFH=SFH[SFH['SnapNum']<=np.max(SnapList)]
    
    Snap_list = set(SFH['SnapNum'])
    
    plt_color = [plt.cm.GnBu(i) for i in np.linspace(0.1, 0.7, len(Snap_list))]   
    #plt_color = [plt.cm.GnBu(i) for i in np.linspace(0.3, 0.9, len(Snap_list))]   
    plt_color_edge = [plt.cm.GnBu(i) for i in np.linspace(0.4, 1.0, len(Snap_list))]   
    
    
    label_list = [x for x in range(1, max(Snap_list)+1,8)]
    #print(TimeList[-1])
    for snap in range(11, max(Snap_list)+1):         
        SFH_isnap = SFH[SFH['SnapNum']==snap] 
       
        for ii in range(0, len(SFH_isnap['SnapNum'])-1): 
            yerr = (np.log10(SFH_isnap['Lookbacktime'][ii]+SFH_isnap['dt'][ii]/2.) - 
                    np.log10(SFH_isnap['Lookbacktime'][ii]-SFH_isnap['dt'][ii]/2.) )/2. 
            y = np.log10(SFH_isnap['Lookbacktime'][ii])
            
            subplot.fill_between([snap-0.2, snap+0.2], [y-yerr, y-yerr], [y+yerr, y+yerr], interpolate=True, 
                                 edgecolor='black', facecolor=plt_color[snap-1]) 
        
        ii+=1
        y_top = np.log10(SFH_isnap['Lookbacktime'][ii]+SFH_isnap['dt'][ii]/2.)+0.02
        y_bottom = 4.0
       
        #subplot.fill_between([snap-0.2, snap+0.2], [y_bottom, y_bottom], [y_top, y_top], interpolate=True, 
        #                     edgecolor='black', facecolor=plt_color[snap-1]) 
        subplot.fill_between([snap-0.2, snap+0.2], [y_bottom, y_bottom], [y_top, y_top], interpolate=True, 
                             edgecolor='black', facecolor=plt_color[snap-1]) 
          
        
    #subplot.set_yscale( "log" )    
    subplot.text(15,1.0,str("Age Bins in $\log_{10}$(yr)"), fontsize=12, color='black') 
     
    
    x_axis_2 =  RedshiftList[0:60:10]   
    
    for ii in range(0, len(x_axis_2)):
        x_axis_2[ii] = int((x_axis_2[ii] * 10) + 0.5) / 10.0            
    x_axis_2[0]=50   
    ax2 = subplot.twiny()        
  
    ax2.set_xbound(subplot.get_xbound())        
    ax2.set_xticks([x for x in range(10,60,10)])
    ax2.set_xticklabels(x_axis_2)
    
    #ax2.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
    ax2.set_xlabel('redshift', fontsize=14)
    ax2.xaxis.set_label_position('top')
                
    for key, spine in ax2.spines.items():
        spine.set_visible(False)
    
    ax2.tick_params(which='both', bottom=False)
       
    
    
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYF19_SFH_bins.pdf')
    plt.close()

    return 
#end SFH_bins 

    
    
    
    
   
        
    