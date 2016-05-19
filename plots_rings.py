'''
gasfractions_vs_stellarmass
milkyway_sfr_and_gas_profiles
evo_milkyway_gas_profile
gas_metallicity_gradients
SFR_gradients
'''

import numpy as np
import pandas as pd
#import seaborn as sns
#sns.set_style('darkgrid')

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.table import Table
from importlib import reload
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
import sys
from scipy.ndimage import zoom

import procedures
reload (procedures)
from procedures import *
import plots_input
reload (plots_input)
from plots_input import *

def milkyway_sfr_and_gas_profiles(G_MR, RingRadius, RNUM, ThisRedshiftList, pdf):
  
    plot_color=['red','purple']        
    plt.rcParams.update({'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 
                             'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0})
    fig = plt.figure(figsize=(10,10))
    grid = gridspec.GridSpec(2, 2)
    #grid.update(wspace=0.0, hspace=0.0)
        
    for ii in range(0,len(ThisRedshiftList)):        
        
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]          
        #G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (G0_MR['DiskMass']>0.) &       
        #            (G0_MR['Vvir']>200.) & (G0_MR['Vvir']<235.) & (G0_MR['Type']==0)                
        #            & (G0_MR['BulgeMass']/G0_MR['StellarMass']<0.15)]
        #            #& (G0_MR['ColdGas']/G0_MR['StellarMass']>0.01)]
        '''G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (G0_MR['DiskMass']>0.) &       
                    (G0_MR['Vvir']>200.) & (G0_MR['Vvir']<235.) & (G0_MR['Type']==0) & 
                    (G0_MR['BulgeMass']/G0_MR['StellarMass']<0.15)]'''
        G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (G0_MR['DiskMass']>0.) &       
                    (G0_MR['Vvir']>200.) & (G0_MR['Vvir']<235.)]
      
        xlim=[0.0,20.0]
        ylim=[0.,3.5]
        
        #Sigmastar
        bin=2.           
        subplot=plt.subplot(grid[0])    
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)   
            
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(5))    
        subplot.xaxis.set_minor_locator(MultipleLocator(1.))            
        xlab='$r[\mathrm{kpc}]$'
        ylab='$\Sigma_{\mathrm{star}}[M_{\odot}/\mathrm{pc^2}]$'     
        subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)   
                   
        Sigma=np.zeros(RNUM,dtype=np.float32)
        Sigma_mean=np.zeros(RNUM,dtype=np.float32)
        pc16=np.zeros(RNUM,dtype=np.float32) 
        pc84=np.zeros(RNUM,dtype=np.float32) 
        for ii in range(0,RNUM):
            y_variable=(G0_MR['DiskMassRings'][:,ii]*1e10/Hubble_h)/(3.14*RingRadius[ii]*RingRadius[ii]*1e6)
            Sigma[ii]=np.median(y_variable)
            Sigma_mean[ii]=np.mean(y_variable)
            y_sorted = np.sort(y_variable)
            pc16[ii] = y_sorted[16*len(y_variable)/100]      
            pc84[ii] = y_sorted[84*len(y_variable)/100]  
         
        subplot.plot(RingRadius, np.log10(Sigma_mean),color='red', linewidth=2, linestyle=':')
        subplot.plot(RingRadius, np.log10(Sigma),color='red', linewidth=2)     
        subplot.plot(RingRadius, np.log10(pc16),color='red', linestyle='--') 
        subplot.plot(RingRadius, np.log10(pc84),color='red', linestyle='--') 
        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.88, 
                    color='black', xlog=0, ylog=0, label='$200<V_{\mathrm{vir}}[km/s]<235$', 
                    fontsize=13, fontweight='normal') 
        
        
        ylim=[0.,3.]
        #SigmaHI       
        bin=2.           
        subplot=plt.subplot(grid[1])    
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)   
            
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(5))    
        subplot.xaxis.set_minor_locator(MultipleLocator(1.))            
        xlab='$r[\mathrm{kpc}]$'
        ylab='$\Sigma_{\mathrm{HI}}[M_{\odot}/\mathrm{pc^2}]$'     
        subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)   
                    
        Sigma=np.zeros(RNUM,dtype=np.float32)
        Sigma_mean=np.zeros(RNUM,dtype=np.float32)
        pc16=np.zeros(RNUM,dtype=np.float32) 
        pc84=np.zeros(RNUM,dtype=np.float32)         
        for ii in range(0,RNUM):           
            Mass=(G0_MR['ColdGasRings'][:,ii]*1e10/Hubble_h*(1.-G0_MR['H2fractionRings'][:,ii]))            
            y_variable=Mass/(3.14*RingRadius[ii]*RingRadius[ii]*1e6)
            Sigma[ii]=np.median(y_variable)
            Sigma_mean[ii]=np.mean(y_variable)
            y_sorted = np.sort(y_variable)
            pc16[ii] = y_sorted[16*len(y_variable)/100]      
            pc84[ii] = y_sorted[84*len(y_variable)/100]  
                 
        subplot.plot(RingRadius, np.log10(Sigma_mean),color='red', linewidth=2, linestyle=':')        
        subplot.plot(RingRadius, np.log10(Sigma),color='red', linewidth=2)     
        subplot.plot(RingRadius, np.log10(pc16),color='red', linestyle='--') 
        subplot.plot(RingRadius, np.log10(pc84),color='red', linestyle='--') 
                
        '''print('StellarMass, DiskMass, BulgeMass,ColdGas,BlackHoleMass')    
        print(np.log10(G0_MR['StellarMass'][0]*1.e10),np.log10(G0_MR['DiskMass'][0]*1.e10),
              np.log10(G0_MR['BulgeMass'][0]*1.e10), np.log10(G0_MR['ColdGas'][0]*1.e10),
              np.log10(G0_MR['BlackHoleMass'][0]*1.e10))
        print(np.log10(G0_MR['StellarMass'][10]*1.e10),np.log10(G0_MR['DiskMass'][10]*1.e10),
              np.log10(G0_MR['BulgeMass'][10]*1.e10),np.log10(G0_MR['ColdGas'][10]*1.e10),
              np.log10(G0_MR['BlackHoleMass'][10]*1.e10))'''
        
        '''print(len(G0_MR))
        for ii in range(0,len(G0_MR)): '''
        
        '''for ii in range(0,len(G0_MR)):          
            Mass=(G0_MR['ColdGasRings'][ii,:]*1e10/Hubble_h*(1-G0_MR['H2fractionRings'][ii,:]))
            y_variable=(3.14*RingRadius*RingRadius*1e6)
            y_variable=Mass/y_variable
            subplot.plot(RingRadius, np.log10(y_variable),color='black', linewidth=2) '''
            
        #SigmaH2
        bin=2.           
        subplot=plt.subplot(grid[2])    
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)   
            
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(5))    
        subplot.xaxis.set_minor_locator(MultipleLocator(1.))            
        xlab='$r[\mathrm{kpc}]$'
        ylab='$\Sigma_{\mathrm{H_2}}[M_{\odot}/\mathrm{pc^2}]$'     
        subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)   
            
               
        Sigma=np.zeros(RNUM,dtype=np.float32)
        Sigma_mean=np.zeros(RNUM,dtype=np.float32)   
        pc16=np.zeros(RNUM,dtype=np.float32) 
        pc84=np.zeros(RNUM,dtype=np.float32) 
        for ii in range(0,RNUM):
            Mass=(G0_MR['ColdGasRings'][:,ii]*1e10/Hubble_h*(G0_MR['H2fractionRings'][:,ii]))
            y_variable=Mass/(3.14*RingRadius[ii]*RingRadius[ii]*1e6)
            Sigma[ii]=np.median(y_variable)
            Sigma_mean[ii]=np.mean(y_variable)
            y_sorted = np.sort(y_variable)
            pc16[ii] = y_sorted[16*len(y_variable)/100]      
            pc84[ii] = y_sorted[84*len(y_variable)/100]  
                 
        subplot.plot(RingRadius, np.log10(Sigma_mean),color='red', linewidth=2, linestyle=':')        
        subplot.plot(RingRadius, np.log10(Sigma),color='red', linewidth=2)     
        subplot.plot(RingRadius, np.log10(pc16),color='red', linestyle='--') 
        subplot.plot(RingRadius, np.log10(pc84),color='red', linestyle='--') 
        
        '''for ii in range(0,len(G0_MR)):          
            Mass=(G0_MR['ColdGasRings'][ii,:]*1e10/Hubble_h*(G0_MR['H2fractionRings'][ii,:]))
            y_variable=(3.14*RingRadius*RingRadius*1e6)
            y_variable=Mass/y_variable
            subplot.plot(RingRadius, np.log10(y_variable),color='black', linewidth=2) '''
        
        #SigmaGas
        bin=2.           
        subplot=plt.subplot(grid[3])    
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)   
            
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(5))    
        subplot.xaxis.set_minor_locator(MultipleLocator(1.))            
        xlab='$r[\mathrm{kpc}]$'
        ylab='$\Sigma_{\mathrm{gas}}[M_{\odot}/\mathrm{pc^2}]$'     
        subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)   
            
               
        Sigma=np.zeros(RNUM,dtype=np.float32)
        Sigma_mean=np.zeros(RNUM,dtype=np.float32)
        pc16=np.zeros(RNUM,dtype=np.float32) 
        pc84=np.zeros(RNUM,dtype=np.float32) 
        for ii in range(0,RNUM):
            Mass=G0_MR['ColdGasRings'][:,ii]*1e10/Hubble_h
            y_variable=Mass/(3.14*RingRadius[ii]*RingRadius[ii]*1e6)
            Sigma[ii]=np.median(y_variable)
            Sigma_mean[ii]=np.mean(y_variable)
            y_sorted = np.sort(y_variable)
            pc16[ii] = y_sorted[16*len(y_variable)/100]      
            pc84[ii] = y_sorted[84*len(y_variable)/100]  
                 
        subplot.plot(RingRadius, np.log10(Sigma_mean),color='red', linewidth=2, linestyle=':')         
        subplot.plot(RingRadius, np.log10(Sigma),color='red', linewidth=2)     
        subplot.plot(RingRadius, np.log10(pc16),color='red', linestyle='--') 
        subplot.plot(RingRadius, np.log10(pc84),color='red', linestyle='--') 
        
        '''for ii in range(0,len(G0_MR)):             
            Mass=(G0_MR['ColdGasRings'][ii,:]*1e10/Hubble_h)
            y_variable=(3.14*RingRadius*RingRadius*1e6)
            y_variable=Mass/y_variable
            subplot.plot(RingRadius, np.log10(y_variable),color='black', linewidth=2)''' 
        
        #LABELS        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.2, 
                    color='black', xlog=0, ylog=0, label='Peeples 2015', 
                    fontsize=13, fontweight='normal') 
        plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.225, 
                    color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15) 
       
            
    plt.tight_layout()
    plt.savefig('./fig/milkyway_sfr_and_gas_profiles.pdf')
    pdf.savefig()
    plt.close()

#end  milkyway_sfr_and_gas_profiles














def gas_metallicity_gradients(G_MR, RingRadius, RNUM, ThisRedshiftList, pdf):
       
    ii=0   
            
    plot_color=['blue','green','red']        
    plt.rcParams.update({'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 
                             'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0})
    
    fig = plt.figure(figsize=(12,5))
    grid = gridspec.GridSpec(1, 2)
    
    #metallicity versus physical radius
    subplot_1=plt.subplot(grid[0])    
    xlim=[0.0,15.0]
    ylim=[8.0, 9.0]  
    subplot_1.set_ylim(ylim), subplot_1.set_xlim(xlim)      
    xlab='$r$[kpc]'           
    ylab='$12 + log_{10}$(O/H)$_{gas}$'               
    subplot_1.set_xlabel(xlab, fontsize=16), subplot_1.set_ylabel(ylab, fontsize=16)  
    
    #metallicity versus physical radius/gas disk size      
    subplot_2=plt.subplot(grid[1])
    xlim=[-1.0,0.1]
    ylim=[8.0, 9.0] 
    subplot_2.set_ylim(ylim), subplot_2.set_xlim(xlim)    
    xlab='$r/r_{d}$'           
    ylab='$12 + log_{10}$(O/H)$_{gas}$'               
    subplot_2.set_xlabel(xlab, fontsize=16), subplot_2.set_ylabel(ylab, fontsize=16)  
      
    median_metallicity=np.zeros(RNUM,dtype=np.float32)
    low_mass_limits=[9.5,10.0,10.5]
    #low_mass_limits=[10.0]
    massbin=0.5   
   

    
    #Model
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)                 
    G0_MR_unsel=G_MR[sel] 
    
    if(opt_detailed_enrichment==1):   
        G0_MR_unsel=G0_MR_unsel[(G0_MR_unsel['MetalsColdGas'][:,0] + 
                                 G0_MR_unsel['MetalsColdGas'][:,1] + 
                                 G0_MR_unsel['MetalsColdGas'][:,2])>.0]
    else:
        G0_MR_unsel=G0_MR_unsel[(G0_MR_unsel['MetalsColdGas']>.0)] #& 
                                #(G0_MR_unsel['BulgeMass']/G0_MR_unsel['StellarMass']<.3)]
        
    for kk in range(0,len(low_mass_limits)):
      
        G0_MR=G0_MR_unsel[(np.log10(G0_MR_unsel['StellarMass']*1.e10/Hubble_h) > low_mass_limits[kk]) & 
                          (np.log10(G0_MR_unsel['StellarMass']*1.e10/Hubble_h) < low_mass_limits[kk]+massbin)]
              
        for jj in range(0,RNUM):             
            if(opt_detailed_enrichment==1):                  
                MetalsColdGasRing=(G0_MR['MetalsColdGasRings'][:,jj,0] +
                                   G0_MR['MetalsColdGasRings'][:,jj,1] +
                                   G0_MR['MetalsColdGasRings'][:,jj,2])
            else:            
                MetalsColdGasRing=G0_MR['MetalsColdGasRings'][:,jj]
                
                                
            ColdGasRing=G0_MR['ColdGasRings'][:,jj]
                      
            metallicity=MetalsColdGasRing[MetalsColdGasRing>0.]/ColdGasRing[MetalsColdGasRing>0.]/0.016            
            median_metallicity[jj]=np.log10(np.median(metallicity))+8.69           
            #median_metallicity[jj]=np.median(np.log10(metallicity))+8.69
               
        #metallicity versus physical radius
        subplot_1.plot(RingRadius, median_metallicity,color=plot_color[kk], linewidth=2)  
        
        '''for jj in range(0,20): 
            metallicity=(G0_MR['MetalsColdGasRings'][jj,0,:] +
                         G0_MR['MetalsColdGasRings'][jj,1,:] +
                         G0_MR['MetalsColdGasRings'][jj,2,:])/G0_MR['ColdGasRings'][jj,:]/0.02
            subplot_1.plot(RingRadius, np.log10(metallicity)+8.69,color=plot_color[kk], linewidth=2)'''
        #print(RingRadius, median_metallicity)
        #metallicity versus physical radius/gas disk size  
        #Rings=np.log10(RingRadius/np.median(G0_MR['GasDiskRadius']*1000./Hubble_h))
        median_radius=np.median(G0_MR['StellarDiskRadius'][G0_MR['StellarDiskRadius']>0.]*1000./Hubble_h)
        Rings=np.log10(RingRadius/median_radius)
        subplot_2.plot(Rings, median_metallicity, color=plot_color[kk], linewidth=2)   
        
        label="%0.1f" % low_mass_limits[kk] + "<$M_{\star}[M_{\odot}]$<" + "%0.1f" % (low_mass_limits[kk]+massbin) 
        plot_label (subplot_2, 'label', xlim, ylim, 
                    x_percentage=0.15, y_percentage=0.21-(kk*0.075), color='black', xlog=0, ylog=0, 
                    label=label, fontsize=12, fontweight='normal')             
        plot_label (subplot_2, 'line', xlim, ylim,
                    x_percentage=0.05, y_percentage=0.22-(kk*0.075), color=plot_color[kk], x2_percentage=0.12, 
                    xlog=0, ylog=0, linestyle='-', linewidth=2)
    
        #endfor
    #endfor
              
    plt.tight_layout()
    plt.savefig('./fig/plots_gas_gradients_1.pdf')
    pdf.savefig()
    plt.close()
    
    
    
    
    
    
    
    #gradient for gas disks
    
    fig = plt.figure(figsize=(12,5))  
    grid = gridspec.GridSpec(1, 2)
    subplot=plt.subplot(grid[0])    
    xlim=[0.0,15.0]
    ylim=[-0.75, 0.75]  
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)      
    xlab='$r$[kpc]'           
    ylab='$\log_{10}$$(Z_{\mathrm{gas}}/Z_\odot)$'               
    subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)  
    
    
    mean_metallicity=np.zeros(RNUM,dtype=np.float32)
     
    #Model
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)                 
    G0_MR_unsel=G_MR[sel] 
   
    if(opt_detailed_enrichment==1): 
        G0_MR=G0_MR_unsel[(G0_MR_unsel['StellarMass']>0.) & (G0_MR_unsel['DiskMass']>0.) & 
                          ((G0_MR_unsel['MetalsDiskMass'][:,0] + G0_MR_unsel['MetalsDiskMass'][:,1] +
                            G0_MR_unsel['MetalsDiskMass'][:,2])>0.) & (G0_MR_unsel['Type']==0) &  
                          #(np.log10(G0_MR_unsel['StellarMass']*1.e10/Hubble_h)>10.5) &
                          #(np.log10(G0_MR_unsel['StellarMass']*1.e10/Hubble_h)>11.) &
                          (G0_MR_unsel['Vvir']>200.) & (G0_MR_unsel['Vvir']<235.)] 
    else:
        G0_MR=G0_MR_unsel[(G0_MR_unsel['StellarMass']>0.) & (G0_MR_unsel['DiskMass']>0.) & 
                          (G0_MR_unsel['MetalsDiskMass']>0.) & (G0_MR_unsel['Type']==0) &  
                          #(np.log10(G0_MR_unsel['StellarMass']*1.e10/Hubble_h)>10.5) &
                          #(np.log10(G0_MR_unsel['StellarMass']*1.e10/Hubble_h)>11.) &
                          (G0_MR_unsel['Vvir']>200.) & (G0_MR_unsel['Vvir']<235.)]   

        
    for jj in range(0,RNUM):             
        if(opt_detailed_enrichment==1):                  
            MetalsColdGasRing=(G0_MR['MetalsColdGasRings'][:,jj,0] +
                               G0_MR['MetalsColdGasRings'][:,jj,1] +
                               G0_MR['MetalsColdGasRings'][:,jj,2])
        else:         
            MetalsColdGasRing=G0_MR['MetalsColdGasRings'][:,jj]
                                               
        ColdGasRing=G0_MR['ColdGasRings'][:,jj]
                      
        metallicity=MetalsColdGasRing[MetalsColdGasRing>0.]/ColdGasRing[MetalsColdGasRing>0.]/0.016         
        mean_metallicity[jj]=np.log10(np.mean(metallicity))
    #endfor  
    
    subplot.plot(RingRadius, mean_metallicity,color='red', linewidth=2)  
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.88, 
                    color='black', xlog=0, ylog=0, label='$200<V_{\mathrm{vir}}[km/s]<235$ (central galaxies)', 
                    fontsize=13, fontweight='normal')   
        
        
        
        
    
    #gradient for stellar disks    
    subplot=plt.subplot(grid[1])    
    xlim=[0.0,15.0]
    ylim=[-0.75, 0.75]  
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)      
    xlab='$r$[kpc]'           
    ylab='$\log_{10}$$(Z_{\star}/Z_\odot)$'               
    subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)  
    
    
    mean_metallicity=np.zeros(RNUM,dtype=np.float32)
     
    #Model
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)                 
    G0_MR_unsel=G_MR[sel] 
  
    if(opt_detailed_enrichment==1): 
        G0_MR=G0_MR_unsel[(G0_MR_unsel['StellarMass']>0.) & (G0_MR_unsel['DiskMass']>0.) & 
                          ((G0_MR_unsel['MetalsDiskMass'][:,0] + G0_MR_unsel['MetalsDiskMass'][:,1] +
                            G0_MR_unsel['MetalsDiskMass'][:,2])>0.) & (G0_MR_unsel['Type']==0) &  
                          #(np.log10(G0_MR_unsel['StellarMass']*1.e10/Hubble_h)>10.5) &
                          #(np.log10(G0_MR_unsel['StellarMass']*1.e10/Hubble_h)>11.) &
                          (G0_MR_unsel['Vvir']>200.) & (G0_MR_unsel['Vvir']<235.)] 
    else:
        G0_MR=G0_MR_unsel[(G0_MR_unsel['StellarMass']>0.) & (G0_MR_unsel['DiskMass']>0.) & 
                          (G0_MR_unsel['MetalsDiskMass']>0.) & (G0_MR_unsel['Type']==0) &  
                          #(np.log10(G0_MR_unsel['StellarMass']*1.e10/Hubble_h)>10.5) &
                          #(np.log10(G0_MR_unsel['StellarMass']*1.e10/Hubble_h)>11.) &
                          (G0_MR_unsel['Vvir']>200.) & (G0_MR_unsel['Vvir']<235.)]   
   
    #print(len(G0_MR))
        
    for jj in range(0,RNUM):             
        if(opt_detailed_enrichment==1):                  
            MetalsColdGasRing=(G0_MR['MetalsDiskMassRings'][:,jj,0] +
                               G0_MR['MetalsDiskMassRings'][:,jj,1] +
                               G0_MR['MetalsDiskMassRings'][:,jj,2])
        else:         
            MetalsColdGasRing=G0_MR['MetalsDiskMassRings'][:,jj]
                                               
        ColdGasRing=G0_MR['DiskMassRings'][:,jj]
                      
        metallicity=MetalsColdGasRing[MetalsColdGasRing>0.]/ColdGasRing[MetalsColdGasRing>0.]/0.016  
        if(len(metallicity)>0.):
            mean_metallicity[jj]=np.log10(np.mean(metallicity))
        else:           
            mean_metallicity[jj]=0.
    #endfor  
    
    subplot.plot(RingRadius, mean_metallicity,color='red', linewidth=2)  
                         
        
    plt.tight_layout()
    plt.savefig('./fig/plots_gas_grandients_2.pdf')
    pdf.savefig()
    plt.close()
    
#end gas_metallicity_gradients











def SFR_gradients(G_MR, RingRadius, RNUM, ThisRedshiftList, pdf):
          
    plt.rcParams.update({'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 
                             'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0})
    fig = plt.figure(figsize=(10,10))
    grid = gridspec.GridSpec(2, 2)
          
    xlim=[0.0,6.0]
    ylim=[-5.0, 0.0]    
    xlab='$r/r_d$[kpc]'           
    ylab='$\Sigma_{\mathrm{SFR}}$'   
          
    low_mass_limits=[9.0,9.5,10.0,10.5]  
    massbin=0.5
           
    #Model
    ii=0
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)                 
    G0_MR_unsel=G_MR[sel] 
    G0_MR_unsel=G0_MR_unsel[G0_MR_unsel['BulgeMass']/G0_MR_unsel['StellarMass']<0.8]
        
    #metallicity versus physical radius
    for kk in range(0,len(low_mass_limits)):       
        mean_SFR=np.zeros(RNUM,dtype=np.float32)
        radius=np.zeros(RNUM,dtype=np.float32)
        G0_MR=G0_MR_unsel[(np.log10(G0_MR_unsel['StellarMass']*1.e10/Hubble_h) > low_mass_limits[kk]) & 
                          (np.log10(G0_MR_unsel['StellarMass']*1.e10/Hubble_h) < low_mass_limits[kk]+massbin)]       
        
        for jj in range(0,RNUM):    
            area=(3.14*RingRadius[jj]*RingRadius[jj])           
            SFR=G0_MR['SfrRings'][:,jj]/area
            mean_SFR[jj]=np.log10(np.mean(SFR))
        
        G=G0_MR[G0_MR['StellarDiskRadius']>0.]
        radius=RingRadius/np.median(G['StellarDiskRadius']/3.*1000./Hubble_h)
        
        subplot=plt.subplot(grid[kk]) 
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)      
        subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16) 
        subplot.plot(radius, mean_SFR,color='red', linewidth=2)      
                
        label="%0.1f" % low_mass_limits[kk] + "$<M_{\star}[M_{\odot}]<$" + "%0.1f"  % (low_mass_limits[kk]+massbin)
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.88, 
                    color='black', xlog=0, ylog=0, label=label, 
                    fontsize=13, fontweight='normal') 
        
        #endfor
    #endfor
   
        
    plt.tight_layout()
    plt.savefig('./fig/plots_metals_vs_stellarmass.pdf')
    pdf.savefig()
    plt.close()
    
#end SFR_gradients













def gasfractions_vs_stellarmass(G_MR, ThisRedshiftList, pdf):
  
    plot_color=['red','purple']        
    plt.rcParams.update({'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 
                             'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0})
    fig = plt.figure(figsize=(20,5))
    grid = gridspec.GridSpec(1, 4)

    for ii in range(0,len(ThisRedshiftList)):        
        
        xlim=[9.5,11.5]
        ylim=[-2.0,0.0]
        
        #ColdGas
        bin=0.25           
        subplot=plt.subplot(grid[0])    
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
            
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))
            
        xlab='$log_{10}(M_*[h^{-2}M_{\odot}])$'
        ylab='$M_{\mathrm{Cold}}/M_*$'     
        subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)   
            
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]          
        #G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (G0_MR['ColdGas']>0.) & 
        #            (G0_MR['Vvir']>120.) & (G0_MR['BulgeMass']/G0_MR['StellarMass']<0.15)]
        G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (G0_MR['ColdGas']>0.) & 
                    (G0_MR['BulgeMass']/G0_MR['StellarMass']<0.15)]
        
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
        Fraction=np.log10(G0_MR['ColdGas']*1.e10*Hubble_h)-StellarMass 
   
        (x_binned, median, mean, pc16, pc84)=median_and_percentiles (bin, xlim[0], xlim[1], StellarMass, Fraction)          
        sel=(median!=0)        
        subplot.plot(x_binned[sel], median[sel],color=plot_color[ii], linewidth=2)     
        subplot.plot(x_binned[sel], pc16[sel],color=plot_color[ii], linewidth=2, linestyle='--')
        subplot.plot(x_binned[sel], pc84[sel],color=plot_color[ii], linewidth=2, linestyle='--')
       
       
            
            
        #HI
        xlim=[9.5,11.5]
        ylim=[-2.0,0.0]
        bin=0.25           
        subplot=plt.subplot(grid[1])    
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
            
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))
            
        xlab='$log_{10}(M_*[h^{-2}M_{\odot}])$'
        ylab='$M_{\mathrm{HI}}/M_*$'     
        subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)   
            
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]          
        G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (G0_MR['ColdGas']>0.) & (G0_MR['H2fraction']<1.)
                    & (G0_MR['Vvir']>120.) & (G0_MR['BulgeMass']/G0_MR['StellarMass']<0.15)]        
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
        Fraction=np.log10(G0_MR['ColdGas']*(1.-G0_MR['H2fraction'])*1.e10*Hubble_h)-StellarMass 
 
        (x_binned, median, mean, pc16, pc84)=median_and_percentiles (bin, xlim[0], xlim[1], StellarMass, Fraction)          
        sel=(median!=0)        
        subplot.plot(x_binned[sel], median[sel],color=plot_color[ii], linewidth=2)     
        subplot.plot(x_binned[sel], pc16[sel],color=plot_color[ii], linewidth=2, linestyle='--')
        subplot.plot(x_binned[sel], pc84[sel],color=plot_color[ii], linewidth=2, linestyle='--')
    
            
        #HII
        xlim=[9.5,11.5]
        ylim=[-2.5,-0.5]
        bin=0.25           
        subplot=plt.subplot(grid[2])    
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
            
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))
            
        xlab='$log_{10}(M_*[h^{-2}M_{\odot}])$'
        ylab='$M_{\mathrm{H_2}}/M_*$'     
        subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)   
            
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]          
        G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (G0_MR['ColdGas']>0.) & (G0_MR['H2fraction']>0.) & 
                    (G0_MR['Vvir']>120.) & (G0_MR['BulgeMass']/G0_MR['StellarMass']<0.15)]
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
        Fraction=np.log10(G0_MR['ColdGas']*(G0_MR['H2fraction'])*1.e10*Hubble_h)-StellarMass 
 
        (x_binned, median, mean, pc16, pc84)=median_and_percentiles (bin, xlim[0], xlim[1], StellarMass, Fraction)          
        sel=(median!=0)        
        subplot.plot(x_binned[sel], median[sel],color=plot_color[ii], linewidth=2)     
        subplot.plot(x_binned[sel], pc16[sel],color=plot_color[ii], linewidth=2, linestyle='--')
        subplot.plot(x_binned[sel], pc84[sel],color=plot_color[ii], linewidth=2, linestyle='--')
        
        
        #HII/HI
        xlim=[9.5,11.5]
        ylim=[-1.5,0.5]
        bin=0.25           
        subplot=plt.subplot(grid[3])    
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
            
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))
            
        xlab='$log_{10}(M_*[h^{-2}M_{\odot}])$'
        ylab='$M_{\mathrm{H_2}}/M_{HI}$'     
        subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)   
            
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]          
        G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (G0_MR['ColdGas']>0.) & (G0_MR['H2fraction']>0.)& 
                    (G0_MR['Vvir']>120.) & (G0_MR['BulgeMass']/G0_MR['StellarMass']<0.15)]
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
        Fraction=np.log10(G0_MR['H2fraction']/(1.-G0_MR['H2fraction']))
 
        (x_binned, median, mean, pc16, pc84)=median_and_percentiles (bin, xlim[0], xlim[1], StellarMass, Fraction)          
        sel=(median!=0)        
        subplot.plot(x_binned[sel], median[sel],color=plot_color[ii], linewidth=2)     
        subplot.plot(x_binned[sel], pc16[sel],color=plot_color[ii], linewidth=2, linestyle='--')
        subplot.plot(x_binned[sel], pc84[sel],color=plot_color[ii], linewidth=2, linestyle='--')
        
    plt.tight_layout()
    plt.savefig('./fig/plots_gasfractions_vs_stellarmass.pdf')
    pdf.savefig()
    plt.close()

#end


def H2fraction_vs_stellarmass(G_MR, ThisRedshiftList, pdf):
  
    plot_color=['red','purple']        
    plt.rcParams.update({'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 
                             'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0})
    fig = plt.figure(figsize=(7,7))
   
    for ii in range(0,len(ThisRedshiftList)):        
                   
        #HII
        xlim=[9.5,11.5]
        ylim=[-2.5,1.0]
        bin=0.25           
        subplot=plt.subplot()    
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
            
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))
            
        xlab='$log_{10}(M_*[h^{-2}M_{\odot}])$'
        ylab='$f_{\mathrm{H_2}}$'     
        subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)   
            
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]          
        G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (G0_MR['ColdGas']>0.) & (G0_MR['H2fraction']>0.)]
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
        #Fraction=np.log10(G0_MR['H2fraction'])
        Fraction=G0_MR['H2fraction']/(1.-G0_MR['H2fraction'])
 
        (x_binned, median, mean, pc16, pc84)=median_and_percentiles (bin, xlim[0], xlim[1], StellarMass, Fraction)          
        sel=(median!=0)        
        subplot.plot(x_binned[sel], np.log10(median[sel]),color=plot_color[ii], linewidth=2)     
        subplot.plot(x_binned[sel], np.log10(pc16[sel]),color=plot_color[ii], linewidth=2, linestyle='--')
        subplot.plot(x_binned[sel], np.log10(pc84[sel]),color=plot_color[ii], linewidth=2, linestyle='--')
      
        
       
            
    plt.tight_layout()
    plt.savefig('./fig/plots_gasfractions_vs_stellarmass.pdf')
    pdf.savefig()
    plt.close()

#end H2fraction_vs_stellarmass











def evo_milkyway_gas_profile(G_MR, RingRadius, RNUM, ThisRedshiftList, pdf):
  
    plot_color=['red','purple']        
    plt.rcParams.update({'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 
                             'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0})
    fig = plt.figure(figsize=(15,5))  
    grid = gridspec.GridSpec(1, 3)
    #grid.update(wspace=0.0, hspace=0.0)
        
    for ii in range(0,len(ThisRedshiftList)):        
        
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]         
        selected_Gal=G0_MR[(G0_MR['StellarMass']>0.) & (G0_MR['DiskMass']>0.) &  
                           (G0_MR['Vvir']/Hubble_h>200.) & (G0_MR['Vvir']/Hubble_h<235.) & (G0_MR['Type']==0) 
                           & (G0_MR['BulgeMass']/G0_MR['StellarMass']<0.15)] 
        print(len(selected_Gal))    
        selected_Gal=selected_Gal[0]              
        MainBranch=G_MR[(G_MR['GalID']>=selected_Gal['GalID']) & (G_MR['GalID']<=selected_Gal['MainLeafId'])]
        print(selected_Gal['GalID'],selected_Gal['MainLeafId'])  
             
            
        print('diskmass=', np.log10(selected_Gal['DiskMass']*1.e10), 
              'disk fraction=', selected_Gal['DiskMass']/selected_Gal['StellarMass'],
              'stellarmass=', np.log10(selected_Gal['StellarMass']*1.e10))    
            
        # Have a look at the colormaps here and decide which one you'd like:
        # http://matplotlib.org/1.2.1/examples/pylab_examples/show_colormaps.html
        num_plots=len(MainBranch)
        colormap = plt.cm.gist_ncar_r
        plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0.2, 1.0, num_plots)])
        #colormap = plt.cm.gist_ncar
        #plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0., 0.8, num_plots)])
            
        xlim=[0.0,30.0]
        ylim=[0.,4.]
               
        #SigmaGas
        bin=2.           
        subplot=plt.subplot(grid[0])    
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)         
        xlab='$r[\mathrm{kpc}]$'
        ylab='$\Sigma_{\mathrm{gas}}[M_{\odot}/\mathrm{pc^2}]$'     
        subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)                       
        colormap = plt.cm.gist_ncar_r
        plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0.2, 1.0, num_plots)])
            
        for jj in range (0,np.amax(MainBranch['SnapNum'])+1): 
        #for jj in range (0,35): 
            Gal=MainBranch[MainBranch['SnapNum']==jj]   
            #print(Gal['SnapNum'],Gal['ColdGasRings'])
            if(len(Gal)>0):
                #print(Gal['SnapNum'])
                Sigma=np.zeros(RNUM,dtype=np.float32)          
                for ii in range(0,RNUM):
                    Mass=Gal['ColdGasRings'][0][ii]*1e10/Hubble_h 
                    Sigma[ii]= Mass/(3.14*RingRadius[ii]*RingRadius[ii]*1e6)                     
                subplot.plot(RingRadius, np.log10(Sigma), linewidth=1)            
          
        #H2
        bin=2.           
        subplot=plt.subplot(grid[1])    
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)      
        xlab='$r[\mathrm{kpc}]$'
        ylab='$\Sigma_{\mathrm{H_2}}[M_{\odot}/\mathrm{pc^2}]$'     
        subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)           
        colormap = plt.cm.gist_ncar_r
        plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0.2, 1.0, num_plots)])    
            
        for jj in range (0,np.amax(MainBranch['SnapNum'])+1):        
            Gal=MainBranch[MainBranch['SnapNum']==jj]          
            if(len(Gal)>0):
                #print(Gal['SnapNum'])
                Sigma=np.zeros(RNUM,dtype=np.float32)          
                for ii in range(0,RNUM):
                    Mass=Gal['ColdGasRings'][0][ii]*1e10/Hubble_h*Gal['H2fractionRings'][0][ii]/1.3                 
                    Sigma[ii]= Mass/(3.14*RingRadius[ii]*RingRadius[ii]*1e6)                     
                subplot.plot(RingRadius, np.log10(Sigma), linewidth=1)            
          

        #Stars
        bin=2.           
        subplot=plt.subplot(grid[2])    
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)        
        xlab='$r[\mathrm{kpc}]$'
        ylab='$\Sigma_{\mathrm{stars}}[M_{\odot}/\mathrm{pc^2}]$'     
        subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)        
        colormap = plt.cm.gist_ncar_r
        plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0.2, 1.0, num_plots)])
              
        for jj in range (0,np.amax(MainBranch['SnapNum'])+1):        
            Gal=MainBranch[MainBranch['SnapNum']==jj]          
            if(len(Gal)>0):
                #print(Gal['SnapNum'])
                Sigma=np.zeros(RNUM,dtype=np.float32)          
                for ii in range(0,RNUM):
                    Mass=Gal['DiskMassRings'][0][ii]*1e10/Hubble_h   
                    Sigma[ii]= Mass/(3.14*RingRadius[ii]*RingRadius[ii]*1e6)                     
                subplot.plot(RingRadius, np.log10(Sigma), linewidth=1)            
          
        #option1
        #color = (0, jj*1./num_plots , 0, 1)
        #option2
        #color=iter(cm.rainbow(np.linspace(0,1,n)))
        #for i in range(n):
        #c=next(color)
        #ax1.plot(x, y,c=c)
         
        #LABELS        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.2, 
                    color='black', xlog=0, ylog=0, label='Peeples 2015', 
                    fontsize=13, fontweight='normal') 
        plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.225, 
                    color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15) 
       
            
    plt.tight_layout()
    plt.savefig('./fig/plots_evo_milkyway_and_gas_profile.pdf')
    pdf.savefig()
    plt.close()

#end evo_milkyway_and_gas_profile






def test_H2_prescriptions(G_MR, RingRadius, RNUM, ThisRedshiftList, pdf):
  
    

    for ii in range(0,len(ThisRedshiftList)):        
                   
        #HII
        xlim=[0.0,4.0]        
        ylim=[-2.5,0.2]
                
        plot_color=['red','purple']        
        plt.rcParams.update({'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 
                             'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0})
        fig = plt.figure(figsize=(15,5))
                   
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]          
        G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (G0_MR['ColdGas']>0.) & (G0_MR['H2fraction']>0.)]
        SigmaGas=np.zeros(len(G0_MR)*RNUM,dtype=np.float64)
        Fraction=np.zeros(len(G0_MR)*RNUM,dtype=np.float64)
        Metallicity=np.zeros(len(G0_MR)*RNUM,dtype=np.float64)
        StellarDensity=np.zeros(len(G0_MR)*RNUM,dtype=np.float64)
        for ii in range(0,RNUM):
            area=(3.14*RingRadius[ii]*RingRadius[ii]*1e6)
            SigmaGas[ii*len(G0_MR):(ii+1)*len(G0_MR)]=np.log10((G0_MR['ColdGasRings'][:,ii]*1e10/Hubble_h)/area)
            Fraction[ii*len(G0_MR):(ii+1)*len(G0_MR)]=G0_MR['H2fractionRings'][:,ii]
            Metallicity[ii*len(G0_MR):(ii+1)*len(G0_MR)]=np.log10(G0_MR['MetalsColdGasRings'][:,ii]/
                                                                  G0_MR['ColdGasRings'][:,ii]/0.02)
            StellarDensity[ii*len(G0_MR):(ii+1)*len(G0_MR)]=np.log10((G0_MR['DiskMassRings'][:,ii]*1e10/Hubble_h)/area)
        #SigmaGas=np.log10((G0_MR['ColdGas']*1e10/Hubble_h)/(3.14*G0_MR['GasDiskRadius']*G0_MR['GasDiskRadius']*1e6*1e6))       
        #Fraction=np.log10(G0_MR['H2fraction'])
        
        #Bin by metallicity
        grid = gridspec.GridSpec(1, 2)
        subplot=plt.subplot(grid[0])
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)                    
        xlab='$log_{10}(\Sigma_{\mathrm{gas}}[M_{\odot}\mathrm{pc}^{-2}])$'
        ylab='$log_{10}(f_{\mathrm{H_2}})$'  
        subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)   
        
        #KMT10 model
        plot_color=['red','red','green','green','blue','blue','black']
        metallicityr=[-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0]
        metallicityr=np.power(10,metallicityr)        
        for ii in range(0,7):               
            SigmaH=np.arange(xlim[0],xlim[1]+0.1,0.001)
            khi=3.1*1.0*(1+3.1*metallicityr[ii]**0.365)/4.1         
            tau=0.066*pow(10,SigmaH)*metallicityr[ii];
            s=np.log(1+0.6*khi+0.01*khi*khi)/(0.6*tau);        
            fraction=np.log10(1-0.75*s/(1+0.25*s));
            subplot.plot(SigmaH, fraction,color=plot_color[ii], linewidth=2, linestyle='--') 
                
        plot_color=['red','green','blue','black']
        x_bins=[-1.25,-0.5,0.25,1.0]  
        bin=0.5
        sel_property=Metallicity        
        for ii in range(0,len(x_bins)):            
            sel=((sel_property>x_bins[ii]-bin) & (sel_property<x_bins[ii]+bin))
            if(len(SigmaGas[sel])>0):
                bin=0.1
                (x_binned, median, mean, pc16, pc84)=median_and_percentiles (bin, xlim[0], xlim[1], 
                                                                             SigmaGas[sel], Fraction[sel])               
                subplot.plot(x_binned, np.log10(median),color=plot_color[ii], linewidth=2)     
               
          
        
        #Bin by stellar surface density
        subplot=plt.subplot(grid[1])
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)              
        subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)  
       
        plot_color=['red','green','blue','black']
        x_bins=[-1.25,-0.5,0.25,1.0]  
        bin=0.5
        sel_property=Metallicity
        x_bins=[1.,10.,100.,1000.]   
        sel_property=StellarDensity
        
        for ii in range(0,len(x_bins)):
            sel=((sel_property>x_bins[ii]-x_bins[ii]*0.1) & (sel_property<x_bins[ii]+x_bins[ii]*0.1))
            #sel=((sel_property>x_bins[ii]-bin) & (sel_property<x_bins[ii]+bin))
            if(len(SigmaGas[sel])>0):
                bin=0.1
                (x_binned, median, mean, pc16, pc84)=median_and_percentiles (bin, xlim[0], xlim[1], 
                                                                             SigmaGas[sel], Fraction[sel])               
                subplot.plot(x_binned, np.log10(median),color=plot_color[ii], linewidth=2)     
        
        #LABELS          
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
    plt.savefig('./fig/plots_gasfractions_vs_stellarmass.pdf')
    pdf.savefig()
    plt.close()

#end test_H2_prescriptions









def test_rings(G_MR, RingRadius, RNUM, ThisRedshiftList, pdf):
        
    ii=0   
        
    plot_color=['blue','green','red']        
    plt.rcParams.update({'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 
                         'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0})
         
    #Gas Mass (Total vs Rings)
    fig = plt.figure(figsize=(7,7))
    subplot_1=plt.subplot()    
    xlim=[5.0,12.0]
    ylim=[5.0,12.0]        
    subplot_1.set_ylim(ylim), subplot_1.set_xlim(xlim)      
    xlab='$\mathrm{Total} - M_{\mathrm{Cold}}$'           
    ylab='$\mathrm{Rings Sum} - M_{\mathrm{Cold}}$'               
    subplot_1.set_xlabel(xlab, fontsize=16), subplot_1.set_ylabel(ylab, fontsize=16)  
    
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)                 
    G0_MR=G_MR[sel]     
    #G0_MR=G0_MR[np.log10(G0_MR['ColdGas']*1e10)>2.]    
    ColdGasRings=np.sum(G0_MR['ColdGasRings'],axis=1)           
    subplot_1.scatter(np.log10(G0_MR['ColdGas'][0:999]*1e10),np.log10(ColdGasRings[0:999]*1e10),s=5, color='black')  
   
    plt.tight_layout()       
    pdf.savefig()
    plt.close()       
            
        
        
    #Mass of Metals in Gas (Total vs Rings)
    fig = plt.figure(figsize=(7,7))
    subplot_1=plt.subplot()    
    xlim=[4.0,11.0]
    ylim=[4.0, 11.0]        
    subplot_1.set_ylim(ylim), subplot_1.set_xlim(xlim)      
    xlab='$\mathrm{Total} - M_{\mathrm{MetalsCold}}$'           
    ylab='$\mathrm{Rings Sum} - M_{\mathrm{MetalsCold}}$'               
    subplot_1.set_xlabel(xlab, fontsize=16), subplot_1.set_ylabel(ylab, fontsize=16)  

    (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)                 
    G0_MR=G_MR[sel]     
    
    if(opt_detailed_enrichment==1):
        G0_MR=G0_MR[np.log10(G0_MR['MetalsColdGas'][:,0]*1.e10)>2.]
    else:
        G0_MR=G0_MR[np.log10(G0_MR['MetalsColdGas']*1.e10)>2.]
        
    if(opt_detailed_enrichment==1):
        MetalscoldGas=G0_MR['MetalsColdGas'][:,0]+G0_MR['MetalsColdGas'][:,1]+G0_MR['MetalsColdGas'][:,2]      
        MetalscoldGasRings=G0_MR['MetalsColdGasRings'][:,:,0] + G0_MR['MetalsColdGasRings'][:,:,1] + \
        G0_MR['MetalsColdGasRings'][:,:,2]        
        MetalscoldGasRings=np.sum(MetalscoldGasRings,axis=1)        
    else:    
        MetalscoldGas=G0_MR['MetalsColdGas']
        MetalscoldGasRings=np.sum(G0_MR['MetalsColdGasRings'],axis=1)
    
    GG=G0_MR[0]   
    subplot_1.scatter(np.log10(MetalscoldGas[0:999]*1e10), np.log10(MetalscoldGasRings[0:999]*1e10),s=5, color='black')  

    plt.tight_layout()       
    pdf.savefig()
    plt.close()   

    
    #Disk Mass (Total vs Rings)
    fig = plt.figure(figsize=(7,7))
    subplot_1=plt.subplot()    
    xlim=[5.0,12.0]
    ylim=[5.0,12.0]        
    subplot_1.set_ylim(ylim), subplot_1.set_xlim(xlim)      
    xlab='$\mathrm{Total} - M_{\mathrm{DiskMass}}$'           
    ylab='$\mathrm{Rings Sum} - M_{\mathrm{DiskMass}}$'               
    subplot_1.set_xlabel(xlab, fontsize=16), subplot_1.set_ylabel(ylab, fontsize=16)  
    
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)                 
    G0_MR=G_MR[sel]     
    G0_MR=G0_MR[np.log10(G0_MR['DiskMass']*1.e10)>2.]    
    ColdGasRings=np.sum(G0_MR['DiskMassRings'],axis=1)           
    subplot_1.scatter(np.log10(G0_MR['DiskMass'][0:999]*1e10),np.log10(ColdGasRings[0:999]*1e10),s=5, color='black')  

    plt.tight_layout()       
    pdf.savefig()
    plt.close() 
    
    #DiskMass Metallicity (Total vs Rings)
    fig = plt.figure(figsize=(7,7))
    subplot_1=plt.subplot()    
    xlim=[6.0,10.0]
    ylim=[6.0,10.0]        
    subplot_1.set_ylim(ylim), subplot_1.set_xlim(xlim)      
    xlab='$\mathrm{Total} - M_{\mathrm{MetalsDiskMass}}$'           
    ylab='$\mathrm{Rings Sum} - M_{\mathrm{MetalsDiskMass}}$'               
    subplot_1.set_xlabel(xlab, fontsize=16), subplot_1.set_ylabel(ylab, fontsize=16)  

    #Model
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)                 
    G0_MR=G_MR[sel]     
    if(opt_detailed_enrichment==1):
        G0_MR=G0_MR[np.log10(G0_MR['MetalsDiskMass'][:,0]*1.e10)>2.]
    else:
        G0_MR=G0_MR[np.log10(G0_MR['MetalsDiskMass']*1.e10)>2.]
        
    if(opt_detailed_enrichment==1):
        MetalscoldGas=G0_MR['MetalsDiskMass'][:,0]+G0_MR['MetalsDiskMass'][:,1]+G0_MR['MetalsDiskMass'][:,2]      
        MetalscoldGasRings=G0_MR['MetalsDiskMassRings'][:,:,0] + G0_MR['MetalsDiskMassRings'][:,:,1] + \
                           G0_MR['MetalsDiskMassRings'][:,:,2] 
        #MetalscoldGas=G0_MR['MetalsDiskMass'][:,1]  
        #MetalscoldGasRings=G0_MR['MetalsDiskMassRings'][:,:,1]      
        MetalscoldGasRings=np.sum(MetalscoldGasRings,axis=1)
    else:
        MetalscoldGas=G0_MR['MetalsDiskMass']        
        MetalscoldGasRings=np.sum(G0_MR['MetalsDiskMassRings'],axis=1)
        
    subplot_1.scatter(np.log10(MetalscoldGas[0:999]*1e10), np.log10(MetalscoldGasRings[0:999]*1e10),s=5, color='black')  

    plt.tight_layout()       
    pdf.savefig()
    plt.close()       

    
    
    #Cold Gas H+He vs cold gas(Total vs Rings)
    if(opt_detailed_enrichment==1):
        fig = plt.figure(figsize=(7,7))
        subplot_1=plt.subplot()    
        xlim=[6.0,10.0]
        ylim=[6.0, 10.0]        
        subplot_1.set_ylim(ylim), subplot_1.set_xlim(xlim)      
        xlab='Total - H in ColdGas'           
        ylab='Rings Sum - ColdGas'               
        subplot_1.set_xlabel(xlab, fontsize=16), subplot_1.set_ylabel(ylab, fontsize=16)    

        #Model
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)                 
        G0_MR=G_MR[sel]    
        
        ColdGas_elements=G0_MR['ColdGas_elements'][:,0]
        ColdGas_elementsRings=G0_MR['ColdGasRings'][:,:]*1e10    
        ColdGas_elementsRings=np.sum(ColdGas_elementsRings,axis=1)
        
        subplot_1.scatter(np.log10(ColdGas_elements[0:2000]), np.log10(ColdGas_elementsRings[0:2000]),s=5, color='red')  
      
        plt.tight_layout()        
        pdf.savefig()
        plt.close()
    
    #Cold Gas Elements(Total vs Rings)
    if(opt_detailed_enrichment==1):
        fig = plt.figure(figsize=(7,7))
        subplot_1=plt.subplot()    
        xlim=[6.0,8.0]
        ylim=[6.0, 8.0]        
        subplot_1.set_ylim(ylim), subplot_1.set_xlim(xlim)      
        xlab='Total - ColdGas in Element Nr6'           
        ylab='Rings Sum - ColdGas in Element Nr6'               
        subplot_1.set_xlabel(xlab, fontsize=16), subplot_1.set_ylabel(ylab, fontsize=16)    

        #Model
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)                 
        G0_MR=G_MR[sel]    
        
        ColdGas_elements=G0_MR['ColdGas_elements'][:,6]     
        ColdGas_elementsRings=G0_MR['ColdGas_elementsRings'][:,:,6]    
        ColdGas_elementsRings=np.sum(ColdGas_elementsRings,axis=1)
        
        subplot_1.scatter(np.log10(ColdGas_elements[0:2000]), np.log10(ColdGas_elementsRings[0:2000]),s=5, color='red')  
      
        plt.tight_layout()        
        pdf.savefig()
        plt.close()
        
    #DiskMass Elements(Total vs Rings)
    if(opt_detailed_enrichment==1):
        fig = plt.figure(figsize=(7,7))
        subplot_1=plt.subplot()    
        xlim=[6.0,8.0]
        ylim=[6.0, 8.0]        
        subplot_1.set_ylim(ylim), subplot_1.set_xlim(xlim)      
        xlab='Total - DiskMass in Element Nr6'           
        ylab='Rings Sum - DiskMass in Element Nr6'               
        subplot_1.set_xlabel(xlab, fontsize=16), subplot_1.set_ylabel(ylab, fontsize=16)    

        #Model
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)                 
        G0_MR=G_MR[sel]    
        
        DiskMass_elements=G0_MR['DiskMass_elements'][:,6]     
        DiskMass_elementsRings=G0_MR['DiskMass_elementsRings'][:,:,6]    
        DiskMass_elementsRings=np.sum(DiskMass_elementsRings,axis=1)
        
        subplot_1.scatter(np.log10(DiskMass_elements[0:2000]), np.log10(DiskMass_elementsRings[0:2000]),s=5, color='red')  
      
        plt.tight_layout()        
        pdf.savefig()
        plt.close()     

        
#end test_rings        