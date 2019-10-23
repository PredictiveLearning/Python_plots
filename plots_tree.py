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
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from matplotlib.patches import Ellipse, Circle
from matplotlib.collections import PatchCollection
from matplotlib import colorbar
from importlib import reload
import random



import procedures
reload (procedures)
from procedures import *
import plots_input
reload (plots_input)
from plots_input import *

def simple_tree_map():
    
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


  
def full_tree_map(object_type):       
    
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
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()

    return 
#full_tree_map 



def all_masses_evo():
       
    xlim=[0.0,7.5]
    ylim=[6.0, 15]
    #ylim=[-4.0, 0.0]
    bin=0.25
    
    model = 'Hen15_test'
    
    compute_and_write = 0
    compute_and_write_SFR = 0
    
    #Mvir_bins = np.arange(12., 16., 1.0)
    #Mvir_bins = np.array([14.,15.])
    #Mvir_bins = np.array([14.1,14.5])
    Mvir_bins = np.array([11.,12.,13.,14.])
    #Mvir_bins = np.array([12.,14.5])
    
    if(len(Mvir_bins)==2):
        fig = plt.figure(figsize=(one_two_size_small[0],one_two_size_small[1])) 
        grid = gridspec.GridSpec(1, 2)
    elif(len(Mvir_bins)==3):
        fig = plt.figure(figsize=(one_three_size_large[0],one_three_size_large[1])) 
        grid = gridspec.GridSpec(1, 3)
    elif(len(Mvir_bins)==4):
        fig = plt.figure(figsize=(one_four_size_large[0],one_four_size_large[1])) 
        grid = gridspec.GridSpec(1, 4)
    grid.update(wspace=0.0, hspace=0.0)    
    
    for idx, element in enumerate(Mvir_bins):
        
        char_mvir="%0.1f" % element
        
        subplot=plt.subplot(grid[idx])    
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(2.0))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
        subplot.yaxis.set_major_locator(MultipleLocator(2.0))    
        subplot.yaxis.set_minor_locator(MultipleLocator(0.5))     
        xlab='$z$'; ylab='$\log_{10}(M[\mathrm{M}_{\odot}])$' 
        
        #the following three lines shouldnt be necessary but they were at some point for this
        #plots to have thick marks inside on all 4 axis 
        #wich refers to major or minor
        subplot.tick_params(axis='both', which='both', direction='in')
        subplot.tick_params(axis='y', which='both', right='on', labelright='off')
        subplot.tick_params(axis='x', which='both', top='on', labelright='off')
        
        if(len(Mvir_bins)==2):           
            if(idx==1):
                ylab=''
                subplot.tick_params(axis='y', which='both',left='on', labelleft='off')                
            subplot.set_xlabel(xlab, fontsize=16); subplot.set_ylabel(ylab, fontsize=16)    
        elif(len(Mvir_bins)==3):      
            if(idx>0):
                ylab=''
                subplot.tick_params(axis='y',labelleft='off')
            subplot.set_xlabel(xlab, fontsize=16); subplot.set_ylabel(ylab, fontsize=16)    
        elif(len(Mvir_bins)==4):          
            if(idx>0):
                ylab=''               
                subplot.tick_params(axis='y', which='both', left='on', labelleft='off')            
            subplot.set_xlabel(xlab, fontsize=16); subplot.set_ylabel(ylab, fontsize=16) 
        
        
        
        if(compute_and_write == 1):    
            #total_ngals=len(G_MR[(G_MR['SnapNum']==62) & (G_MR['Type']==0) & 
            #                     (G_MR['Mvir']>(10**(element-0.5))/1.e10) & (G_MR['Mvir']<(10**(element+0.5))/1.e10)])
            total_ngals=len(G_MR[(G_MR['SnapNum']==58) & (G_MR['Type']==0) & 
                                 (G_MR['Mvir']>(10**(element-0.5))/1.e10) & (G_MR['Mvir']<(10**(element+0.5))/1.e10)])
            Nsnaps=64
            #Ngals=2
            Ngals=2000
            if(Ngals>total_ngals):
                Ngals=total_ngals
            print('Ngals per mass bin=',Ngals)       
            #G0_MR=np.random.choice(G_MR[(G_MR['SnapNum']==62) & (G_MR['Type']==0) & 
            #                            (G_MR['Mvir']>(10**(element-0.5))/1.e10) & (G_MR['Mvir']<(10**(element+0.5))/1.e10)],
            #                       size=Ngals, replace=False) 
            G0_MR=np.random.choice(G_MR[(G_MR['SnapNum']==58) & (G_MR['Type']==0) & 
                                        (G_MR['Mvir']>(10**(element-0.5))/1.e10) & (G_MR['Mvir']<(10**(element+0.5))/1.e10)],
                                   size=Ngals, replace=False) 
              
            Redshift=np.zeros(Nsnaps*Ngals,dtype=np.float32)    
            Mvir=np.zeros(Nsnaps*Ngals,dtype=np.float32)
            ColdGas=np.zeros(Nsnaps*Ngals,dtype=np.float32)
            HotGas=np.zeros(Nsnaps*Ngals,dtype=np.float32)
            EjectedGas=np.zeros(Nsnaps*Ngals,dtype=np.float32)
            StellarMass=np.zeros(Nsnaps*Ngals,dtype=np.float32)
            BlackHoleMass=np.zeros(Nsnaps*Ngals,dtype=np.float32)
            SFR=np.zeros(Nsnaps*Ngals,dtype=np.float32)
          
            for ii in range (0,Ngals): 
                initial_gal=G0_MR['GalID'][ii]
                G0_MR_tree=G_MR[(G_MR['GalID']>=G0_MR['GalID'][ii]) & (G_MR['GalID']<=G0_MR['LastProgGal'][ii]) 
                               & (G_MR['Type']<2)]    
      
                snap=0  
                currentgal=initial_gal
                while True:    
                    Gal = G0_MR_tree[G0_MR_tree['GalID']==currentgal]        
                    Redshift[Nsnaps*ii+snap] = Gal['Redshift']        
                    Mvir[Nsnaps*ii+snap] = np.log10(Gal['Mvir']*1.e10/Hubble_h)            
                    ColdGas[Nsnaps*ii+snap]=np.log10(Gal['ColdGas']*1.e10/Hubble_h)
                    HotGas[Nsnaps*ii+snap]=np.log10(Gal['HotGas']*1.e10/Hubble_h)
                    EjectedGas[Nsnaps*ii+snap]=np.log10(Gal['EjectedMass']*1.e10/Hubble_h)
                    StellarMass[Nsnaps*ii+snap]=np.log10(Gal['StellarMass']*1.e10/Hubble_h)
                    BlackHoleMass[Nsnaps*ii+snap]=np.log10(Gal['BlackHoleMass']*1.e10/Hubble_h)
                    #SFR[Nsnaps*ii+snap]=np.log10(Gal['Sfr'])
           
                    prog=Gal['FirstProgGal']
                    currentgal=prog    
                    snap+=1
                    if prog==-1:
                        break
                    
                if(ii%(Ngals/10)==0):
                    print(ii, Ngals)
      
            #select random galaxy to overplot
            if(Ngals>1):
                random_gal = np.random.randint(1, Ngals, size=1)  
            else:
                random_gal=[0]
            ran_gal_z             = Redshift[Nsnaps*random_gal[0]:Nsnaps*(random_gal[0]+1)]
            ran_gal_Mvir          = Mvir[Nsnaps*random_gal[0]:Nsnaps*(random_gal[0]+1)]
            ran_gal_EjectedGas    = EjectedGas[Nsnaps*random_gal[0]:Nsnaps*(random_gal[0]+1)]
            ran_gal_HotGas        = HotGas[Nsnaps*random_gal[0]:Nsnaps*(random_gal[0]+1)]
            ran_gal_ColdGas       = ColdGas[Nsnaps*random_gal[0]:Nsnaps*(random_gal[0]+1)]
            ran_gal_StellarMass   = StellarMass[Nsnaps*random_gal[0]:Nsnaps*(random_gal[0]+1)]
            ran_gal_BlackHoleMass = BlackHoleMass[Nsnaps*random_gal[0]:Nsnaps*(random_gal[0]+1)]
        
            sel_ran               = ran_gal_Mvir > 0.        
            ran_gal_z             = ran_gal_z[sel_ran]
            ran_gal_Mvir          = ran_gal_Mvir[sel_ran]       
            ran_gal_EjectedGas    = ran_gal_EjectedGas[sel_ran]       
            ran_gal_HotGas        = ran_gal_HotGas[sel_ran]       
            ran_gal_ColdGas       = ran_gal_ColdGas[sel_ran]       
            ran_gal_StellarMass   = ran_gal_StellarMass[sel_ran]       
            ran_gal_BlackHoleMass = ran_gal_BlackHoleMass[sel_ran]       
        
            #write to file
            (x_binned, median, mean, pc16, pc84, rms)= median_and_percentiles_fixed_xx(Redshift, Mvir) 
            fa = open(Datadir+"all_masses_evo_"+model+'_mvir_'+char_mvir+"_mvir.txt", "w")
            fa.write("%d\n" % len(x_binned))
            for kk in range (0,len(x_binned)):               
                fa.write("%0.2f " % x_binned[kk] + "%0.2f " % median[kk]  + "%0.2f " % pc16[kk]  + "%0.2f\n" % pc84[kk])   
            fa.close() 
            
            (x_binned, median, mean, pc16, pc84, rms)= median_and_percentiles_fixed_xx(Redshift, EjectedGas)  
            fa = open(Datadir+"all_masses_evo_"+model+'_mvir_'+char_mvir+"_EjectedGas.txt", "w")
            fa.write("%d\n" % len(x_binned))
            for kk in range (0,len(x_binned)):               
                fa.write("%0.2f " % x_binned[kk] + "%0.2f " % median[kk]  + "%0.2f " % pc16[kk]  + "%0.2f\n" % pc84[kk])   
            fa.close() 
            
            (x_binned, median, mean, pc16, pc84, rms)= median_and_percentiles_fixed_xx(Redshift, HotGas)    
            fa = open(Datadir+"all_masses_evo_"+model+'_mvir_'+char_mvir+"_HotGas.txt", "w")
            fa.write("%d\n" % len(x_binned))
            for kk in range (0,len(x_binned)):               
                fa.write("%0.2f " % x_binned[kk] + "%0.2f " % median[kk]  + "%0.2f " % pc16[kk]  + "%0.2f\n" % pc84[kk])   
            fa.close() 
            
            (x_binned, median, mean, pc16, pc84, rms)= median_and_percentiles_fixed_xx(Redshift, ColdGas)        
            fa = open(Datadir+"all_masses_evo_"+model+'_mvir_'+char_mvir+"_ColdGas.txt", "w")
            fa.write("%d\n" % len(x_binned))
            for kk in range (0,len(x_binned)):               
                fa.write("%0.2f " % x_binned[kk] + "%0.2f " % median[kk]  + "%0.2f " % pc16[kk]  + "%0.2f\n" % pc84[kk])   
            fa.close() 
            
            (x_binned, median, mean, pc16, pc84, rms)= median_and_percentiles_fixed_xx(Redshift, StellarMass)    
            fa = open(Datadir+"all_masses_evo_"+model+'_mvir_'+char_mvir+"_StellarMass.txt", "w")
            fa.write("%d\n" % len(x_binned))
            for kk in range (0,len(x_binned)):               
                fa.write("%0.2f " % x_binned[kk] + "%0.2f " % median[kk]  + "%0.2f " % pc16[kk]  + "%0.2f\n" % pc84[kk])   
            fa.close() 
            
            (x_binned, median, mean, pc16, pc84, rms)= median_and_percentiles_fixed_xx(Redshift, BlackHoleMass)    
            fa = open(Datadir+"all_masses_evo_"+model+'_mvir_'+char_mvir+"_BlackHoleMass.txt", "w")
            fa.write("%d\n" % len(x_binned))
            for kk in range (0,len(x_binned)):               
                fa.write("%0.2f " % x_binned[kk] + "%0.2f " % median[kk]  + "%0.2f " % pc16[kk]  + "%0.2f\n" % pc84[kk])   
            fa.close() 
        
        
        
        if(compute_and_write_SFR == 1):            
            total_ngals=len(G_MR[(G_MR['SnapNum']==58) & (G_MR['Type']==0) & 
                                 (G_MR['Mvir']>(10**(element-0.5))/1.e10) & (G_MR['Mvir']<(10**(element+0.5))/1.e10)])
            Nsnaps=64          
            Ngals=1000
            if(Ngals>total_ngals):
                Ngals=total_ngals
            print('Ngals per mass bin=',Ngals)          
            G0_MR=np.random.choice(G_MR[(G_MR['SnapNum']==58) & (G_MR['Type']==0) & 
                                        (G_MR['Mvir']>(10**(element-0.5))/1.e10) & (G_MR['Mvir']<(10**(element+0.5))/1.e10)],
                                   size=Ngals, replace=False) 
              
            Redshift=np.zeros(Nsnaps*Ngals,dtype=np.float32)   
      
            SFR=np.zeros(Nsnaps*Ngals,dtype=np.float32)
          
            for ii in range (0,Ngals): 
                initial_gal=G0_MR['GalID'][ii]
                G0_MR_tree=G_MR[(G_MR['GalID']>=G0_MR['GalID'][ii]) & (G_MR['GalID']<=G0_MR['LastProgGal'][ii]) 
                               & (G_MR['Type']<2)]    
      
                snap=0  
                currentgal=initial_gal
                while True:    
                    Gal = G0_MR_tree[G0_MR_tree['GalID']==currentgal]        
                    Redshift[Nsnaps*ii+snap] = Gal['Redshift']        
                    SFR[Nsnaps*ii+snap]=np.log10(Gal['Sfr'])
           
                    prog=Gal['FirstProgGal']
                    currentgal=prog    
                    snap+=1
                    if prog==-1:
                        break
                    
                if(ii%(Ngals/10)==0):
                    print(ii, Ngals)
                
            #write to file
            (x_binned, median, mean, pc16, pc84, rms)= median_and_percentiles_fixed_xx(Redshift, SFR) 
            fa = open(Datadir+"all_masses_evo_"+model+'_mvir_'+char_mvir+"_SFR.txt", "w")
            fa.write("%d\n" % len(x_binned))
            for kk in range (0,len(x_binned)):               
                fa.write("%0.2f " % x_binned[kk] + "%0.2f " % median[kk]  + "%0.2f " % pc16[kk]  + "%0.2f\n" % pc84[kk])   
            fa.close() 
                    
        #end if(compute_and_write == 1):
        
        
           
        fa = Datadir+"all_masses_evo_"+model+'_mvir_'+char_mvir+"_mvir.txt"
        (x_binned,median, pc16, pc84)=read_data_with_err(fa)  
        subplot.fill_between(x_binned, pc84, pc16, interpolate=True,alpha=0.8, facecolor='grey', edgecolor='grey') 
        subplot.plot(x_binned, median, linestyle='-', linewidth=2, color='black', alpha=0.8)
        #subplot.plot(ran_gal_z, ran_gal_Mvir, linestyle='-', linewidth=2, color='black')
        
        
        
        #Mvir=12 line
        '''yy=np.arange(ylim[0]-0.5,ylim[1]+0.5,0.1)
        xx=np.zeros(len(yy))
        if(len(x_binned[median>=12.0])>0):
            xx+=np.amax(x_binned[median>=12.0])            
            subplot.plot(xx,yy, linestyle='--', color='black', alpha=0.5)'''   
                    
        fa = Datadir+"all_masses_evo_"+model+'_mvir_'+char_mvir+"_EjectedGas.txt"
        (x_binned,median, pc16, pc84)=read_data_with_err(fa)  
        sel = pc16<ylim[0]
        pc16[sel]=ylim[0]
        sel = (pc84!=pc16)
        subplot.fill_between(x_binned[sel], pc84[sel], pc16[sel], interpolate=True, 
                             alpha=0.4, facecolor='purple', edgecolor='purple') 
        subplot.plot(x_binned[sel], median[sel], linestyle='-', linewidth=2, color='purple', alpha=0.8)
        #subplot.plot(ran_gal_z, ran_gal_EjectedGas, linestyle='-', linewidth=2, color='grey')        
        eject = median
        
        fa = Datadir+"all_masses_evo_"+model+'_mvir_'+char_mvir+"_HotGas.txt"
        (x_binned,median, pc16, pc84)=read_data_with_err(fa)  
        sel = pc16<ylim[0]
        pc16[sel]=ylim[0]
        subplot.fill_between(x_binned, pc84, pc16, interpolate=True,alpha=0.4, facecolor='red', edgecolor='red') 
        subplot.plot(x_binned, median, linestyle='-', linewidth=2, color='red', alpha=0.8)
        #subplot.plot(ran_gal_z, ran_gal_HotGas, linestyle='-', linewidth=2, color='red')
        hot = median
        
        #hot greater than eject line
        yy=np.arange(ylim[0]-0.5,ylim[1]+0.5,0.1)
        xx_1=np.zeros(len(yy))
        if(len(x_binned[hot>=eject])>0):
            xx_1+=np.amax(x_binned[(hot>=eject) & (x_binned<xlim[1])])            
            #subplot.plot(xx_1,yy, linestyle='--', color='black',alpha=0.1) 
        #print(xx,yy)
                
        fa = Datadir+"all_masses_evo_"+model+'_mvir_'+char_mvir+"_ColdGas.txt"
        (x_binned,median, pc16, pc84)=read_data_with_err(fa)
        subplot.fill_between(x_binned, pc84, pc16, interpolate=True,alpha=0.4, facecolor='blue', edgecolor='blue') 
        subplot.plot(x_binned, median, linestyle='-', linewidth=2, color='blue', alpha=0.8)
        #subplot.plot(ran_gal_z, ran_gal_ColdGas, linestyle='-', linewidth='2', color='blue')
        
        fa = Datadir+"all_masses_evo_"+model+'_mvir_'+char_mvir+"_StellarMass.txt"
        (x_binned,median, pc16, pc84)=read_data_with_err(fa)
        subplot.fill_between(x_binned, pc84, pc16, interpolate=True,alpha=0.4, facecolor='orange', edgecolor='orange') 
        subplot.plot(x_binned, median, linestyle='-', linewidth=2, color='chocolate', alpha=0.8)
        #subplot.plot(ran_gal_z, ran_gal_StellarMass, linestyle='-', linewidth='2', color='orange')
        StellarMass = median
        x_binned_1 = x_binned
                
        fa = Datadir+"all_masses_evo_"+model+'_mvir_'+char_mvir+"_BlackHoleMass.txt"
        (x_binned,median, pc16, pc84)=read_data_with_err(fa)
        sel = pc16<ylim[0]
        pc16[sel]=ylim[0]
        sel = pc84<ylim[0]
        pc84[sel]=ylim[0]
        
        subplot.fill_between(x_binned, pc84, pc16, interpolate=True,alpha=0.4, facecolor='brown', edgecolor='brown') 
        subplot.plot(x_binned, median, linestyle='-', linewidth=2, color='brown', alpha=0.8)
        #subplot.plot(ran_gal_z, ran_gal_BlackHoleMass, linestyle='-', linewidth='2', color='brown')
       
    
        fa = Datadir+"all_masses_evo_"+model+'_mvir_'+char_mvir+"_SFR.txt"
        (x_binned,median, pc16, pc84)=read_data_with_err(fa)
        #print(x_binned,median)
        #subplot.fill_between(x_binned, pc84, pc16, interpolate=True,alpha=0.4, facecolor='blue', edgecolor='blue') 
        #subplot.plot(x_binned, median+9., linestyle='--', linewidth=2, color='black', alpha=0.8)
        SFR = median
        x_binned_2 = x_binned
               
        #SFR peak line
        yy=np.arange(ylim[0]-0.5,ylim[1]+0.5,0.1)
        xx=np.zeros(len(yy))  
        max_median = median==np.amax(median)        
        xx+=np.mean(x_binned[max_median])
        if(idx>0):
            subplot.plot(xx,yy, linestyle='--', color='black', alpha=0.5)   
        
        #quenching line
        yy=np.arange(ylim[0]-0.5,ylim[1]+0.5,0.1)
        xx=np.arange(xlim[0]-0.5,xlim[1]+0.5,0.1)        
        StellarMass=np.interp(xx, x_binned_1, StellarMass) 
        SFR=np.interp(xx, x_binned_2, SFR)   
        
        #print(x_binned[(SFR-StellarMass)<=np.log10((1+x_binned)/(2.*1.37e10))])
        if(len(xx[(SFR-StellarMass)<=np.log10((1+xx)/(2.*1.37e10))])>0):
            xx_2 = xx+np.amax(xx[(SFR-StellarMass)<=np.log10((1+xx)/(2.*1.37e10))])
        else:
            xx_2=xx
       
        #fill region between quenching and hot_over_ejected
        xx = np.arange(xlim[0]-0.5,xlim[1]+0.5,0.01)
        xx = xx[(xx<xx_1[0]) & (xx>xx_2[0])]
        yy_1 = np.zeros(len(xx)) + ylim[0]
        yy_2 = np.zeros(len(xx)) + ylim[1]     
        #subplot.fill_between(xx, yy_1, yy_2, interpolate=True,alpha=0.4, hatch='/', color='grey',
        #                     facecolor='grey', edgecolor='grey', linestyle = '--') 
        subplot.fill_between(xx, yy_1, yy_2, interpolate=True,alpha=0.3, color='grey', facecolor='grey') 
        
    
        if(idx==0):
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.15, y_percentage=0.93, color='black', xlog=0, ylog=0,
                        label='$M_{200\mathrm{c}}$', fontsize=10)    
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.05, y_percentage=0.945, 
                        color='grey', x2_percentage=0.12, xlog=0, ylog=0, linestyle='-', linewidth=2) 
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.15, y_percentage=0.88, color='black', xlog=0, ylog=0,
                        label='$M_{\mathrm{ejected}}$', fontsize=10)    
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.05, y_percentage=0.895, 
                        color='purple', x2_percentage=0.12, xlog=0, ylog=0, linestyle='-', linewidth=2) 
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.15, y_percentage=0.83, color='black', xlog=0, ylog=0,
                        label='$M_{\mathrm{hot}}$', fontsize=10)    
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.05, y_percentage=0.845, 
                        color='red', x2_percentage=0.12, xlog=0, ylog=0, linestyle='-', linewidth=2) 
            
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.45, y_percentage=0.93, color='black', xlog=0, ylog=0,
                        label='$M_{\mathrm{cold}}$', fontsize=10)    
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.35, y_percentage=0.945, 
                        color='blue', x2_percentage=0.42, xlog=0, ylog=0, linestyle='-', linewidth=2) 
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.45, y_percentage=0.88, color='black', xlog=0, ylog=0,
                        label='$M_*$', fontsize=10)    
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.35, y_percentage=0.895, 
                        color='orange', x2_percentage=0.42, xlog=0, ylog=0, linestyle='-', linewidth=2) 
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.45, y_percentage=0.83, color='black', xlog=0, ylog=0,
                        label='$M_{\mathrm{BH}}$', fontsize=10)    
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.35, y_percentage=0.845, 
                        color='brown', x2_percentage=0.42, xlog=0, ylog=0, linestyle='-', linewidth=2) 
                    
        label='$M_{\mathrm{200c,0}}\sim$%0.1f' % element
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.62, y_percentage=0.9, color='black', 
                    xlog=0, ylog=0, label=label, fontsize=10)    
        
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HWL18_mass_fractions_evo_all_masses.pdf')
    plt.close()
    
    return    
#end  mass_fractions_evo_all_masses   
    
    
    
    
    
    
def mass_fractions_evo_all_single_gal():    
    
    plt.rcParams.update({'font.size': 12, 'xtick.labelsize': 12, 'ytick.labelsize': 12})    
    
    print('\n\n')
    print('Doing single plot')
    
    #SINGLE PLOT
    xlim=[0.0,11.]
    ylim=[6.0, 15.]
    bin=0.25

    fig = plt.figure(figsize=(one_one_size_small[0]+1.0,one_one_size_small[1]))
    
    #num_gals=len(G_MR[(G_MR['SnapNum']==58) & (G_MR['Mvir']>10000.0) & (G_MR['Type']==0)]) 
    num_gals=len(G_MRII[(G_MRII['SnapNum']==62) & (G_MRII['Mvir']>1000.0) & (G_MRII['Type']==0)]) 
    print('Ngals=',num_gals)
       
    subplot=plt.subplot()    
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(2.0))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
    subplot.yaxis.set_major_locator(MultipleLocator(2.0))    
    subplot.yaxis.set_minor_locator(MultipleLocator(0.5)) 

    xlab='$z$'
    ylab='$\log_{10}(M[\mathrm{M}_{\odot}])$'          
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)    
            
    #select one random galaxy to plot    
    #idx = np.random.randint(1, num_gals, size=1)    
    
    #use file 0 and the standard Hen15 output catalog
    #idx = 4
    #G0_MR = G_MR[(G_MR['SnapNum']==58) & (G_MR['Mvir']>10000.0) & (G_MR['Type']==0)]
    
    #use file 41 and the standard Hen15 MRII output catalog
    idx=0   
    G0_MR = G_MRII[(G_MRII['SnapNum']==62) & (G_MRII['Mvir']>1000.0) & (G_MRII['Type']==0)]    
    
    #StellarMass = G_MRII['StellarMass']*1.e10/Hubble_h
    #Sfr = G_MRII['Sfr']
    #log_SSFR = np.log10(Sfr/StellarMass)
    #log_StellarMass= np.log10(StellarMass)
    #G0_MR = G_MRII[(G_MRII['SnapNum']==62) & (log_StellarMass>9.0) & (log_StellarMass<9.1) 
    #               & (log_SSFR<-11.0) & (G_MRII['Type']==0)]
    print(len(G0_MR))
    #idx=0
    initial_gal=G0_MR['GalID'][idx]
    G0_MR_tree=G_MRII[(G_MRII['GalID']>=G0_MR['GalID'][idx]) & (G_MRII['GalID']<=G0_MR['LastProgGal'][idx])]    
      
    print('Ngals in tree=',len(G0_MR_tree))
    
    print(np.log10(G0_MR_tree['StellarMass'][G0_MR_tree['GalID']==initial_gal]*1e10/Hubble_h))
    print(np.log10(G0_MR_tree['Mvir'][G0_MR_tree['GalID']==initial_gal]*1e10/Hubble_h))
    
    Nsnaps=0  
    currentgal=initial_gal
    while True:        
        prog=G0_MR_tree['FirstProgGal'][G0_MR_tree['GalID']==currentgal]
        currentgal=prog       
        Nsnaps+=1
        if prog==-1:
            break
    
    Redshift=np.zeros(Nsnaps,dtype=np.float32)    
    Mvir=np.zeros(Nsnaps,dtype=np.float32)
    EjectedMass=np.zeros(Nsnaps,dtype=np.float32)
    HotGas=np.zeros(Nsnaps,dtype=np.float32)
    ColdGas=np.zeros(Nsnaps,dtype=np.float32)
    StellarMass=np.zeros(Nsnaps,dtype=np.float32)
    StellarMassICM=np.zeros(Nsnaps,dtype=np.float32)
    BlackHoleMass=np.zeros(Nsnaps,dtype=np.float32)
    SFR=np.zeros(Nsnaps,dtype=np.float32)
    
    
    Nsnaps=0  
    currentgal=initial_gal
    while True:    
        Redshift[Nsnaps]=G0_MR_tree['Redshift'][G0_MR_tree['GalID']==currentgal]        
        Mvir[Nsnaps]=np.log10(G0_MR_tree['Mvir'][G0_MR_tree['GalID']==currentgal]*1.e10/Hubble_h)
        EjectedMass[Nsnaps]=np.log10(G0_MR_tree['EjectedMass'][G0_MR_tree['GalID']==currentgal]*1.e10/Hubble_h)
        HotGas[Nsnaps]=np.log10(G0_MR_tree['HotGas'][G0_MR_tree['GalID']==currentgal]*1.e10/Hubble_h)
        ColdGas[Nsnaps]=np.log10(G0_MR_tree['ColdGas'][G0_MR_tree['GalID']==currentgal]*1.e10/Hubble_h)
        StellarMass[Nsnaps]=np.log10(G0_MR_tree['StellarMass'][G0_MR_tree['GalID']==currentgal]*1.e10/Hubble_h)
        StellarMassICM[Nsnaps]=np.log10(G0_MR_tree['StellarMass'][G0_MR_tree['GalID']==currentgal]*1.e10/Hubble_h + 
                                        G0_MR_tree['ICM'][G0_MR_tree['GalID']==currentgal]*1.e10/Hubble_h)
        BlackHoleMass[Nsnaps]=np.log10(G0_MR_tree['BlackHoleMass'][G0_MR_tree['GalID']==currentgal]*1.e10/Hubble_h)
        SFR[Nsnaps]=np.log10(G0_MR_tree['Sfr'][G0_MR_tree['GalID']==currentgal])
           
        prog=G0_MR_tree['FirstProgGal'][G0_MR_tree['GalID']==currentgal]
        currentgal=prog    
        Nsnaps+=1
        if prog==-1:
            break
    
    
    subplot.plot(Redshift, Mvir, linestyle='-', linewidth=2, color='black')
    subplot.plot(Redshift, EjectedMass, linestyle='-', linewidth=2, color='purple')
    subplot.plot(Redshift, HotGas, linestyle='-', linewidth=2, color='red')
    subplot.plot(Redshift, ColdGas, linestyle='-', linewidth=2, color='blue')
    subplot.plot(Redshift, StellarMass, linestyle='-', linewidth=2, color='orange')
    subplot.plot(Redshift, StellarMassICM, linestyle='--', linewidth=2, color='orange')    
    subplot.plot(Redshift, BlackHoleMass, linestyle='-', linewidth=2, color='brown')
        
        
    #yy=np.arange(ylim[0]-0.5,ylim[1]+0.5,0.1)
    #xx=np.zeros(len(yy))
    #xx+=np.amax(Redshift[Mvir>=12.0])
    #subplot.plot(xx,yy, linestyle=':', linewidth=1, color='black')
    
    yy=np.arange(ylim[0]-0.5,ylim[1]+0.5,0.1)
    xx=np.zeros(len(yy))
    xx+=2.85
    subplot.plot(xx,yy, linestyle='--', linewidth=1, color='black', alpha=0.5)
  
    #plot_label (subplot, 'label', xlim, ylim, x_percentage=0.415, y_percentage=0.9, 
    #            color='black', xlog=0, ylog=0, label='$\log_{10}(M_{\mathrm{vir}})<12.0$', fontsize=10)
    
    #plot_label (subplot, 'label', xlim, ylim, x_percentage=0.08, y_percentage=0.9, 
    #            color='black', xlog=0, ylog=0, label='$\log_{10}(M_{\mathrm{vir}})>12.0$', fontsize=10)

    #plot_label (subplot, 'label', xlim, ylim, x_percentage=0.23, y_percentage=0.9, color='black', xlog=0, ylog=0,
    #            label='$\log_{10}(M_{\mathrm{200c}}/M_{\odot})>12.0$', fontsize=10,backgroundcolor='white',back_alpha=0.8)    
    #plot_label (subplot, 'line', xlim, ylim, x_percentage=0.31, y_percentage=0.85, 
    #            color='black', x2_percentage=0.41, xlog=0, ylog=0, linestyle='-', linewidth=2)    
    #plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.31, y_percentage=0.85, 
    #            color='black', xlog=0, ylog=0, sym='<', sym_size=5, err_size=0.0001,mfc='white')
    
    yy=np.arange(ylim[0]-0.5,ylim[1]+0.5,0.1)
    xx_1=np.zeros(len(yy))
    xx_1+=np.amax(Redshift[(SFR-StellarMass)<=np.log10((1+Redshift)/(2.*1.37e10))])
    #subplot.plot(xx_1,yy, linestyle=':', linewidth=1, color='black')
   
    yy=np.arange(ylim[0]-0.5,ylim[1]+0.5,0.1)
    xx_2=np.zeros(len(yy))
    xx_2+=np.amax(Redshift[EjectedMass<=HotGas])
    #subplot.plot(xx_2,yy, linestyle=':', linewidth=1, color='black')    

    xx = np.arange(xlim[0]-0.5,xlim[1]+0.5,0.01)
    xx = xx[(xx>xx_1[0]) & (xx<xx_2[0])]
    yy_1 = np.zeros(len(xx)) + ylim[0]
    yy_2 = np.zeros(len(xx)) + ylim[1]     
    #subplot.fill_between(xx, yy_1, yy_2, interpolate=True,alpha=0.1, hatch='/', color='grey',
    #                     facecolor='grey', edgecolor='grey', linestyle = '--') 
    subplot.fill_between(xx, yy_1, yy_2, interpolate=True,alpha=0.3, color='grey', facecolor='grey') 
        
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.08, y_percentage=0.915, color='black', xlog=0, ylog=0,
                label='Quenched', fontsize=10,backgroundcolor='white',back_alpha=0.6)    
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.21, y_percentage=0.87, 
                color='black', x2_percentage=0.11, xlog=0, ylog=0, linestyle='-', linewidth=2)    
    plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.11, y_percentage=0.87, 
                color='black', xlog=0, ylog=0, sym='<', sym_size=5, err_size=0.00001,mfc='white')
    
    
    #plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.31, color='black', xlog=0, ylog=0,
    #            label='~4.2 Gyr', fontsize=7,backgroundcolor='white',back_alpha=0.2) 
    #plot_label (subplot, 'label', xlim, ylim, x_percentage=0.2, y_percentage=0.31, color='black', xlog=0, ylog=0,
    #            label='~2.2 Gyr', fontsize=7,backgroundcolor='white',back_alpha=0.2) 
    
    #plot_label (subplot, 'line', xlim, ylim, x_percentage=0.055, y_percentage=0.36, 
    #            color='black', x2_percentage=0.355, xlog=0, ylog=0, linestyle='-', linewidth=1)    
    #subplot.plot([0.6,0.6],[9.1,9.4],color='black',linestyle='-', linewidth=1)
    #subplot.plot([1.75,1.75],[9.1,9.4],color='black',linestyle='-', linewidth=1)
    #subplot.plot([3.95,3.95],[9.1,9.4],color='black',linestyle='-', linewidth=1)
    
       
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.67, y_percentage=0.91, color='black', xlog=0, ylog=0,
                label='$M_{\mathrm{200c}}$', fontsize=10)    
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.6, y_percentage=0.925, 
                color='grey', x2_percentage=0.656, xlog=0, ylog=0, linestyle='-', linewidth=2) 
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.67, y_percentage=0.86, color='black', xlog=0, ylog=0,
                label='$M_{\mathrm{ejected}}$', fontsize=10)    
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.6, y_percentage=0.875, 
                color='purple', x2_percentage=0.656, xlog=0, ylog=0, linestyle='-', linewidth=2) 
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.67, y_percentage=0.81, color='black', xlog=0, ylog=0,
                label='$M_{\mathrm{hot}}$', fontsize=10)    
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.6, y_percentage=0.825, 
                color='red', x2_percentage=0.656, xlog=0, ylog=0, linestyle='-', linewidth=2)     
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.67, y_percentage=0.76, color='black', xlog=0, ylog=0,
                label='$M_{\mathrm{cold}}$', fontsize=10)    
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.6, y_percentage=0.775, 
                color='blue', x2_percentage=0.656, xlog=0, ylog=0, linestyle='-', linewidth=2) 
    
   
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.87, y_percentage=0.91, color='black', xlog=0, ylog=0,
                label='$M_*$', fontsize=10)    
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.8, y_percentage=0.925, 
                color='orange', x2_percentage=0.856, xlog=0, ylog=0, linestyle='-', linewidth=2) 
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.87, y_percentage=0.86, color='black', xlog=0, ylog=0,
                label='$M_{*+\mathrm{ICM}}$', fontsize=10)    
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.8, y_percentage=0.875, 
                color='orange', x2_percentage=0.856, xlog=0, ylog=0, linestyle='--', linewidth=2)     
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.87, y_percentage=0.81, color='black', xlog=0, ylog=0,
                label='$M_{\mathrm{BH}}$', fontsize=10)    
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.8, y_percentage=0.825, 
                color='brown', x2_percentage=0.856, xlog=0, ylog=0, linestyle='-', linewidth=2)     
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.87, y_percentage=0.75, color='black', xlog=0, ylog=0,
                label='$SFR$', fontsize=10)    
    plot_label (subplot, 'line', xlim, ylim, x_percentage=0.8, y_percentage=0.765, 
                color='black', x2_percentage=0.856, xlog=0, ylog=0, linestyle='--', linewidth=2) 

    #Second x-axis
    file=Datadir+'Database_snapshots_table.fits'
    fits_table=fits.open(file)
    snap_table = fits_table[1]
         
    xx=np.arange(xlim[0],xlim[1],1.0)   
    x_axis=np.interp(xx, snap_table.data['z'][::-1], snap_table.data['lookBackTime'][::-1]) 
       
    for jj in range(0,len(x_axis)):
        x_axis[jj] = int((x_axis[jj] * 10) + 0.5) / 10.0            
             
            
    plt.rcParams.update({'font.size': 10, 'xtick.labelsize': 10, 'ytick.labelsize': 10})            
    ax2 = subplot.twiny()        
    ax2.set_xbound(subplot.get_xbound())        
    ax2.set_xticks(xx)
    ax2.set_xticklabels(x_axis)
    plt.rcParams.update({'font.size': 12, 'xtick.labelsize': 12, 'ytick.labelsize': 12})    

    xlab='Age (Gyr)'  
    ax2.set_xlabel(xlab, fontsize=14)
    ax2.xaxis.set_label_position('top')  
    
    #Second Axis    
    ylim2=[-2.5, 5.5] 
    ylim2=[-1.7, 5.5] 
    subplot_2 = subplot.twinx()             
    subplot_2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    subplot_2.yaxis.set_minor_locator(MultipleLocator(0.5)) 
    subplot_2.yaxis.set_major_locator(MultipleLocator(2.0)) 
    subplot_2.set_ylim(ylim2) 
    subplot_2.set_xlim(xlim) 
        
    ylab='$\log_{10}(SFR[\mathrm{M}_{\odot}\mathrm{yr}^{-1}])$'  
    subplot_2.set_ylabel(ylab, fontsize=14)
    subplot_2.yaxis.set_label_position('right')  
    
    subplot_2.plot(Redshift, SFR, linestyle='--', linewidth=2, color='black')
       
    #plot_label (subplot, 'label', xlim, ylim, x_percentage=0.87, y_percentage=0.23, 
    #            color='black', xlog=0, ylog=0, label='$SFR$', fontsize=13, fontweight='normal')
        
       
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HWL18_mass_fractions_evo_all_single_gal.pdf')
    plt.close()
    
    return 
#end mass_fractions_evo_all_single_gal





def mass_fractions_evo_all_multiple_gals():    
    
    plt.rcParams.update({'font.size': 12, 'xtick.labelsize': 12, 'ytick.labelsize': 12})    
    
    print('\n\n')
    print('Doing single plot')
    
    #SINGLE PLOT
    xlim=[0.0,11.]
    ylim=[6.0, 15.]
    bin=0.25

    fig = plt.figure(figsize=(one_one_size_small[0]+1.0,one_one_size_small[1]))
    
    num_gals=len(G_MR[(G_MR['SnapNum']==58) & (G_MR['Mvir']>10000.0) & (G_MR['Type']==0)]) 
    print('Ngals=',num_gals)
       
    subplot=plt.subplot()    
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(2.0))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
    subplot.yaxis.set_major_locator(MultipleLocator(2.0))    
    subplot.yaxis.set_minor_locator(MultipleLocator(0.5)) 

    xlab='$z$'
    ylab='$\log_{10}(M[\mathrm{M}_{\odot}])$'          
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)    
            
    #select one random galaxy to plot    
    idx = np.random.randint(1, num_gals, size=1)    
    #use file 0 and the standard Hen15 output catalog
    
    
    #for ii in range(0,9):
    for ii in range(0,1):
        idx = 3
        G0_MR = G_MR[(G_MR['SnapNum']==58) & (G_MR['Mvir']>10000.0) & (G_MR['Type']==0)]
        initial_gal=G0_MR['GalID'][idx]
        G0_MR_tree=G_MR[(G_MR['GalID']>=G0_MR['GalID'][idx]) & (G_MR['GalID']<=G0_MR['LastProgGal'][idx])]    
      
        print('Ngals in tree=',len(G0_MR_tree))
    
        print(np.log10(G0_MR_tree['StellarMass'][G0_MR_tree['GalID']==initial_gal]*1e10/Hubble_h))
        print(np.log10(G0_MR_tree['Mvir'][G0_MR_tree['GalID']==initial_gal]*1e10/Hubble_h))
    
        Nsnaps=0  
        currentgal=initial_gal
        while True:        
            prog=G0_MR_tree['FirstProgGal'][G0_MR_tree['GalID']==currentgal]
            currentgal=prog       
            Nsnaps+=1
            if prog==-1:
                break
    
        Redshift=np.zeros(Nsnaps,dtype=np.float32)    
        Mvir=np.zeros(Nsnaps,dtype=np.float32)  
        EjectedMass=np.zeros(Nsnaps,dtype=np.float32)
        HotGas=np.zeros(Nsnaps,dtype=np.float32)
        ColdGas=np.zeros(Nsnaps,dtype=np.float32)
        StellarMass=np.zeros(Nsnaps,dtype=np.float32)    
        BlackHoleMass=np.zeros(Nsnaps,dtype=np.float32)
        SFR=np.zeros(Nsnaps,dtype=np.float32)
    
        Nsnaps=0  
        currentgal=initial_gal
        while True:    
            Redshift[Nsnaps]=G0_MR_tree['Redshift'][G0_MR_tree['GalID']==currentgal]        
            Mvir[Nsnaps]=np.log10(G0_MR_tree['Mvir'][G0_MR_tree['GalID']==currentgal]*1.e10/Hubble_h)  
            EjectedMass[Nsnaps]=np.log10(G0_MR_tree['EjectedMass'][G0_MR_tree['GalID']==currentgal]*1.e10/Hubble_h)
            HotGas[Nsnaps]=np.log10(G0_MR_tree['HotGas'][G0_MR_tree['GalID']==currentgal]*1.e10/Hubble_h)       
            ColdGas[Nsnaps]=np.log10(G0_MR_tree['ColdGas'][G0_MR_tree['GalID']==currentgal]*1.e10/Hubble_h)
            StellarMass[Nsnaps]=np.log10(G0_MR_tree['StellarMass'][G0_MR_tree['GalID']==currentgal]*1.e10/Hubble_h)     
            BlackHoleMass[Nsnaps]=np.log10(G0_MR_tree['BlackHoleMass'][G0_MR_tree['GalID']==currentgal]*1.e10/Hubble_h)
            SFR[Nsnaps]=np.log10(G0_MR_tree['Sfr'][G0_MR_tree['GalID']==currentgal])
           
            prog=G0_MR_tree['FirstProgGal'][G0_MR_tree['GalID']==currentgal]
            currentgal=prog    
            Nsnaps+=1
            if prog==-1:
                break
       
        subplot.plot(Redshift, EjectedMass, linestyle='-', linewidth=1, color='purple')
        subplot.plot(Redshift, HotGas, linestyle='-', linewidth=1, color='red')
        subplot.plot(Redshift, ColdGas, linestyle='-', linewidth=1, color='blue')
        subplot.plot(Redshift, StellarMass, linestyle='-', linewidth=1, color='orange')    
        subplot.plot(Redshift, BlackHoleMass, linestyle='-', linewidth=1, color='brown')
        
        
        yy=np.arange(ylim[0]-0.5,ylim[1]+0.5,0.1)
        xx=np.zeros(len(yy))
        xx+=np.amax(Redshift[Mvir>=12.0])
        subplot.plot(xx,yy, linestyle=':', linewidth=1, color='black')
   
        #Second x-axis
        file=Datadir+'Database_snapshots_table.fits'
        fits_table=fits.open(file)
        snap_table = fits_table[1]
         
        xx=np.arange(xlim[0],xlim[1],1.0)   
        x_axis=np.interp(xx, snap_table.data['z'][::-1], snap_table.data['lookBackTime'][::-1]) 
       
        for jj in range(0,len(x_axis)):
            x_axis[jj] = int((x_axis[jj] * 10) + 0.5) / 10.0            
             
            
        plt.rcParams.update({'font.size': 10, 'xtick.labelsize': 10, 'ytick.labelsize': 10})            
        ax2 = subplot.twiny()        
        ax2.set_xbound(subplot.get_xbound())        
        ax2.set_xticks(xx)
        ax2.set_xticklabels(x_axis)
        plt.rcParams.update({'font.size': 12, 'xtick.labelsize': 12, 'ytick.labelsize': 12})    

        xlab='Age (Gyr)'  
        ax2.set_xlabel(xlab, fontsize=14)
        ax2.xaxis.set_label_position('top')  
    
        #Second Axis    
        ylim2=[-2.5, 5.5]        
        subplot_2 = subplot.twinx()             
        subplot_2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        subplot_2.yaxis.set_minor_locator(MultipleLocator(0.5)) 
        subplot_2.yaxis.set_major_locator(MultipleLocator(2.0)) 
        subplot_2.set_ylim(ylim2) 
        subplot_2.set_xlim(xlim) 
        
        ylab='$\log_{10}(SFR[\mathrm{M}_{\odot}\mathrm{yr}^{-1}])$'  
        subplot_2.set_ylabel(ylab, fontsize=14)
        subplot_2.yaxis.set_label_position('right')  
    
        subplot_2.plot(Redshift, SFR, linestyle='-', linewidth=1, color='black')
        
       
    plt.tight_layout()    
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()
    
    return 
#end mass_fractions_evo_all_multiple_gals




def properties_evo():    
    
    plt.rcParams.update({'font.size': 12, 'xtick.labelsize': 12, 'ytick.labelsize': 12})    
    
    print('\n\n')
    print('Doing single plot')
    
    #SINGLE PLOT
    xlim=[-0.5,5.]
    ylim=[6.0, 15.]
    bin=0.25

    fig = plt.figure(figsize=(two_two_size_large[0]+1.0,two_two_size_large[1]))
    grid = gridspec.GridSpec(2, 2)
        
        
    num_gals=len(G_MR[(G_MR['SnapNum']==58)  & (G_MR['Type']==0) & 
                      (np.log10(G_MR['StellarMass']/Hubble_h*1.e10)>11.0) & 
                      (np.log10(G_MR['Sfr']/(G_MR['StellarMass']/Hubble_h*1.e10))>-11.0)]) 
    print('Ngals=',num_gals)
       
    subplot=plt.subplot(grid[0])    
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(2.0))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
    subplot.yaxis.set_major_locator(MultipleLocator(2.0))    
    subplot.yaxis.set_minor_locator(MultipleLocator(0.5)) 

    xlab='$z$'
    ylab='$\log_{10}(M[\mathrm{M}_{\odot}])$'          
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)    
            
    my_linestyle=['-','--','-.',':','-','-','-','-','-','-','-','-','-',
                  '-','-','-','-','-','-','-','-','-','-','-','-','-']
    
    IDs = [40000006000005, 40000011000005, 40000012000182, 40000016000005, 40008118003404,
           40008120000005, 40008126000073, 40008148000005, 40008165000005, 40008167000005]
    
    IDs = [40000016000005, 40008126000073, 40008148000005]
    #IDs = [40008120000005, 40008126000073, 40008148000005]
    
    IDs = [40008126000073]
    
    IDs = [40008118003404,40008120000005, 40008126000073]
    
    #num_gals = min(10, num_gals)
    num_gals = len(IDs)
    
    #for idx in range(0,num_gals): 
    for idx in range(0,num_gals):     
        #G0_MR = G_MR[(G_MR['SnapNum']==58) & (G_MR['Type']==0) &
        #             (np.log10(G_MR['StellarMass']/Hubble_h*1.e10)>11.0) & 
        #             (np.log10(G_MR['Sfr']/(G_MR['StellarMass']/Hubble_h*1.e10))>-11.0)]
        #initial_gal=G0_MR['GalID'][idx]
        #G0_MR_tree=G_MR[(G_MR['GalID']>=G0_MR['GalID'][idx]) & (G_MR['GalID']<=G0_MR['LastProgGal'][idx])]    
              
        G0_MR = G_MR[(G_MR['SnapNum']==58) & (G_MR['HaloID']==IDs[idx]) & (G_MR['Type']==0)]
        initial_gal=G0_MR['GalID']
        G0_MR_tree=G_MR[(G_MR['GalID']>=G0_MR['GalID']) & (G_MR['GalID']<=G0_MR['LastProgGal'])]    
      
    
        mass = np.log10(G0_MR_tree['StellarMass'][G0_MR_tree['GalID']==initial_gal]*1.e10/Hubble_h)
        halomass = np.log10(G0_MR_tree['Mvir'][G0_MR_tree['GalID']==initial_gal]*1.e10/Hubble_h)
        print('Mass=%0.2f, Mvir=%0.2f, Ngals in tree=%d' % (mass, halomass, len(G0_MR_tree)))
        #print(G0_MR_tree['HaloID'][G0_MR_tree['GalID']==initial_gal])
    
        Nsnaps=0  
        currentgal=initial_gal
        while True:        
            prog=G0_MR_tree['FirstProgGal'][G0_MR_tree['GalID']==currentgal]
            currentgal=prog       
            Nsnaps+=1
            if prog==-1:
                break
    
        Redshift=np.zeros(Nsnaps,dtype=np.float32)    
        Mvir=np.zeros(Nsnaps,dtype=np.float32)  
        EjectedMass=np.zeros(Nsnaps,dtype=np.float32)
        HotGas=np.zeros(Nsnaps,dtype=np.float32)
        ColdGas=np.zeros(Nsnaps,dtype=np.float32)
        StellarMass=np.zeros(Nsnaps,dtype=np.float32)    
        BlackHoleMass=np.zeros(Nsnaps,dtype=np.float32)
        SFR=np.zeros(Nsnaps,dtype=np.float32)
    
        Nsnaps=0  
        currentgal=initial_gal
        while True:    
            Redshift[Nsnaps]=G0_MR_tree['Redshift'][G0_MR_tree['GalID']==currentgal]        
            Mvir[Nsnaps]=np.log10(G0_MR_tree['Mvir'][G0_MR_tree['GalID']==currentgal]*1.e10/Hubble_h)  
            EjectedMass[Nsnaps]=np.log10(G0_MR_tree['EjectedMass'][G0_MR_tree['GalID']==currentgal]*1.e10/Hubble_h)
            HotGas[Nsnaps]=np.log10(G0_MR_tree['HotGas'][G0_MR_tree['GalID']==currentgal]*1.e10/Hubble_h)       
            ColdGas[Nsnaps]=np.log10(G0_MR_tree['ColdGas'][G0_MR_tree['GalID']==currentgal]*1.e10/Hubble_h)
            StellarMass[Nsnaps]=np.log10(G0_MR_tree['StellarMass'][G0_MR_tree['GalID']==currentgal]*1.e10/Hubble_h)     
            BlackHoleMass[Nsnaps]=np.log10(G0_MR_tree['BlackHoleMass'][G0_MR_tree['GalID']==currentgal]*1.e10/Hubble_h)
            SFR[Nsnaps]=np.log10(G0_MR_tree['Sfr'][G0_MR_tree['GalID']==currentgal])
           
            prog=G0_MR_tree['FirstProgGal'][G0_MR_tree['GalID']==currentgal]
            currentgal=prog    
            Nsnaps+=1
            if prog==-1:
                break
       
        #subplot.plot(Redshift, EjectedMass, linestyle=my_linestyle[idx], linewidth=1, color='purple')
        subplot.plot(Redshift, HotGas, linestyle=my_linestyle[idx], linewidth=1, color='red')
        subplot.plot(Redshift, ColdGas, linestyle=my_linestyle[idx], linewidth=1, color='blue')
        subplot.plot(Redshift, StellarMass, linestyle=my_linestyle[idx], linewidth=1, color='orange')    
        subplot.plot(Redshift, BlackHoleMass, linestyle=my_linestyle[idx], linewidth=1, color='brown')
        subplot.plot(Redshift, SFR-StellarMass+20., linestyle=my_linestyle[idx], linewidth=1, color='black')
        
       
    plt.tight_layout()    
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()
    
    return 
#end properties_evo










    
    


















def halo_growth():
         
    plt.rcParams.update({'font.size': 14, 'xtick.labelsize': 12, 'ytick.labelsize': 12})            
        
    file=Datadir+'Database_snapshots_table.fits'
    fits_table=fits.open(file)
    snap_table = fits_table[1]       
    FullRedshiftList=snap_table.data['z'][::-1]
    FullTimeList=snap_table.data['lookBackTime'][::-1]
    FullTimeList=FullTimeList[FullRedshiftList>0.]
    FullRedshiftList=FullRedshiftList[FullRedshiftList>0.]
      
    plot_color=['gray','blue','red','green']
    Property_Name=['$M_{\mathrm{vir}}$','$M_{\mathrm{cold}}$','$M_*$','$SFR$']
    plot_linestyle=[':','-.','--','-']
    
    #SINGLE PLOT
    xlim=[0.0,0.9]  
    
    fig = plt.figure(figsize=(two_two_size_small[0],two_two_size_small[1]))
    grid = gridspec.GridSpec(2, 2)
    grid.update(wspace=0.0, hspace=0.0)
    
    for i_plot in range(0,4):
        if(i_plot==0):
            ylim=[9.0, 15.]
            ylab='$\log_{10}(M_{\mathrm{vir}}[\mathrm{M}_{\odot}])$'
        else:
            if(i_plot==1):
                ylim=[7.0, 12.]
                ylab='$\log_{10}(M_{\mathrm{Cold}}[\mathrm{M}_{\odot}])$'
            else:
                if(i_plot==2):
                    ylim=[7.0, 12.]
                    ylab='$\log_{10}(M_*[\mathrm{M}_{\odot}])$'
                else:
                    ylim=[-1.0, 3.]
                    ylab='$\log_{10}(SFR[\mathrm{M}_{\odot}\mathrm{yr}^{-1}])$'
                    
        subplot=plt.subplot(grid[i_plot])  
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(0.5))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.1)) 
        subplot.yaxis.set_major_locator(MultipleLocator(2.0))    
        subplot.yaxis.set_minor_locator(MultipleLocator(0.5)) 
        
        if(i_plot==2 or i_plot==3):
            xlab='$\log_{10}(1+z)$'
            subplot.set_xlabel(xlab, fontsize=14)
            
        subplot.set_ylabel(ylab, fontsize=14)    
             
       
  
        if i_plot==1 or i_plot==3:
            subplot.tick_params(axis='y', which='both', left='on', labelleft='off', labelright='on')
            subplot.yaxis.set_label_position("right")

        if i_plot==0 or i_plot==1:    
            subplot.tick_params(axis='x', which='both', bottom='on', labelbottom='off')    
            
            
        HaloMassBins=np.arange(11.0,15.0,1.0)
    
        for ii in range(0, len(HaloMassBins)):
        
            G0_MR=G_MR[(G_MR['SnapNum']==58) & (G_MR['Type']==0) & 
                       (np.log10(G_MR['Mvir']*1.e10/Hubble_h)>HaloMassBins[ii]-0.1) & 
                       (np.log10(G_MR['Mvir']*1.e10/Hubble_h)<HaloMassBins[ii]+0.1)]  
        
            print('HaloMass=',HaloMassBins[ii],'Ntrees=',len(G0_MR))
       
            NTrees=min(100,len(G0_MR))
            G0_MR=G0_MR[0:NTrees]
            
            Property_Array=np.zeros([len(FullRedshiftList),NTrees],dtype=np.float32)
          
            #for jj in range(0,len(G0_MR)):
            for jj in range(0,NTrees):
                initial_gal=G0_MR['GalID'][jj]                 
                G0_MR_tree=G_MR[(G_MR['GalID']>=G0_MR['GalID'][jj]) & (G_MR['GalID']<=G0_MR['LastProgGal'][jj])]    
               
                Nsnaps=0  
                currentgal=initial_gal
                while True:
                    sel=G0_MR_tree['GalID']==currentgal
                    if(i_plot==0):
                        Property_Array[Nsnaps,jj]=G0_MR_tree['Mvir'][sel]*1.e10/Hubble_h
                    else:
                        if(i_plot==1):
                            Property_Array[Nsnaps,jj]=G0_MR_tree['ColdGas'][sel]*1.e10/Hubble_h
                        else:
                            if(i_plot==2):
                                Property_Array[Nsnaps,jj]=G0_MR_tree['StellarMass'][sel]*1.e10/Hubble_h
                            else:
                                Property_Array[Nsnaps,jj]=G0_MR_tree['Sfr'][sel]
                  
                    prog=G0_MR_tree['FirstProgGal'][G0_MR_tree['GalID']==currentgal]
                    currentgal=prog    
                    Nsnaps+=1
                    if prog==-1:
                        break
  
            (x_array,median,mean,pc16,pc84,rms)=convert_2d_array_into_region(FullRedshiftList, Property_Array)              
            #if(i_plot==0 or i_plot==2):
            #    subplot.fill_between(np.log10(x_array+1),np.log10(pc84),np.log10(pc16),interpolate=True,alpha=0.4,
            #                         facecolor=plot_color[i_plot],edgecolor=plot_color[i_plot])         
            subplot.plot(np.log10(x_array+1),np.log10(median), linestyle=plot_linestyle[ii], linewidth=2, color=plot_color[i_plot])
              
                
        #2nd axis        
        if(i_plot<2):
            #create an array with the dimension of the previous axis
            major_ticks=np.arange(xlim[0],xlim[1]+0.1,0.2) 
            minor_ticks=np.arange(xlim[0],xlim[1]+0.1,0.1) 
            #interpolate the values of the new axis
            x_array_new=np.interp(major_ticks, np.log10(FullRedshiftList+1), FullTimeList)                
            #one decimal place
            for iii in range(0,len(x_array_new)):
                x_array_new[iii] = int((x_array_new[iii] * 10) + 0.5) / 10.0            
              
            ax2 = subplot.twiny()        
            ax2.set_xbound(subplot.get_xbound())        
            ax2.set_xticks(major_ticks)
            #ax2.set_xticks(minor_ticks, minor = True)
            ax2.set_xticklabels(x_array_new)
                
            xlab='LookBackTime(Gyr)'  
            ax2.set_xlabel(xlab, fontsize=14)
            ax2.xaxis.set_label_position('top')  
            
    
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.7, y_percentage=0.8, 
                    color='black', xlog=0, ylog=0, label=Property_Name[i_plot], fontsize=13, fontweight='normal')     
   
    plt.rcParams.update({'font.size': 14, 'xtick.labelsize': 14, 'ytick.labelsize': 14})            
       
    plt.tight_layout()   
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()
  
    return 
#halo_growth    



def halo_growth_normalized():
       
    plt.rcParams.update({'font.size': 14, 'xtick.labelsize': 12, 'ytick.labelsize': 12})        
        
        
    file=Datadir+'Database_snapshots_table.fits'
    fits_table=fits.open(file)
    snap_table = fits_table[1]       
    FullRedshiftList=snap_table.data['z'][::-1]
    FullTimeList=snap_table.data['lookBackTime'][::-1]
    FullTimeList=FullTimeList[FullRedshiftList>0.]
    FullRedshiftList=FullRedshiftList[FullRedshiftList>0.]
      
    plot_color=['gray','blue','red','green']
    Property_Name=['$M_{\mathrm{vir}}$','$M_{\mathrm{cold}}$','$M_*$','$SFR$']
    plot_linestyle=[':','-.','--','-']
    
    #SINGLE PLOT
    xlim=[0.0,0.9]  
    
    fig = plt.figure(figsize=(two_two_size_small[0],two_two_size_small[1]))
    grid = gridspec.GridSpec(2, 2)
    grid.update(wspace=0.0, hspace=0.0)
    
    for i_plot in range(0,4):
        if(i_plot==0):
            ylim=[-1.9, 0.]
            ylab='$\log_{10}(M_{\mathrm{vir}}(z)/M_{\mathrm{vir}}(z_0))$'
        else:
            if(i_plot==1):
                ylim=[-1.9, 2.]
                ylab='$\log_{10}(M_{\mathrm{Cold}}(z)/M_{\mathrm{Cold}}(z_0))$'
            else:
                if(i_plot==2):
                    ylim=[-1.9, 0.]
                    ylab='$\log_{10}(M_*(z)/M_*(z_0))$'
                else:
                    ylim=[-1.0, 3.]
                    ylab='$\log_{10}(SFR[\mathrm{M}_{\odot}\mathrm{yr}^{-1}])$'
                    
        subplot=plt.subplot(grid[i_plot])  
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(0.5))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.1)) 
        subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
        subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
        
        if(i_plot==2 or i_plot==3):
            xlab='$\log_{10}(1+z)$'
            subplot.set_xlabel(xlab, fontsize=14)
            
        subplot.set_ylabel(ylab, fontsize=14)    
             
       
  
        if i_plot==1 or i_plot==3:
            subplot.tick_params(axis='y', which='both', left='on', labelleft='off', labelright='on')
            subplot.yaxis.set_label_position("right")

        if i_plot==0 or i_plot==1:    
            subplot.tick_params(axis='x', which='both', bottom='on', labelbottom='off')    
            
            
        HaloMassBins=np.arange(11.0,15.0,1.0)
    
        for ii in range(0, len(HaloMassBins)):
        
            G0_MR=G_MR[(G_MR['SnapNum']==58) & (G_MR['Type']==0) & 
                       (np.log10(G_MR['Mvir']*1.e10/Hubble_h)>HaloMassBins[ii]-0.1) & 
                       (np.log10(G_MR['Mvir']*1.e10/Hubble_h)<HaloMassBins[ii]+0.1)]  
        
            print('HaloMass=',HaloMassBins[ii],'Ntrees=',len(G0_MR))
       
            NTrees=min(100,len(G0_MR))
            G0_MR=G0_MR[0:NTrees]
            
            Property_Array=np.zeros([len(FullRedshiftList),NTrees],dtype=np.float32)
          
            #for jj in range(0,len(G0_MR)):
            for jj in range(0,NTrees):
                initial_gal=G0_MR['GalID'][jj]                 
                G0_MR_tree=G_MR[(G_MR['GalID']>=G0_MR['GalID'][jj]) & (G_MR['GalID']<=G0_MR['LastProgGal'][jj])]    
               
                Nsnaps=0  
                currentgal=initial_gal
                while True:
                    sel=G0_MR_tree['GalID']==currentgal
                    if(i_plot==0):
                        Property_Array[Nsnaps,jj]=G0_MR_tree['Mvir'][sel]/G0_MR['Mvir'][jj]
                    else:
                        if(i_plot==1):
                            Property_Array[Nsnaps,jj]=G0_MR_tree['ColdGas'][sel]/G0_MR['ColdGas'][jj]
                        else:
                            if(i_plot==2):
                                Property_Array[Nsnaps,jj]=G0_MR_tree['StellarMass'][sel]/G0_MR['StellarMass'][jj]
                            else:
                                Property_Array[Nsnaps,jj]=G0_MR_tree['Sfr'][sel]
                  
                    prog=G0_MR_tree['FirstProgGal'][G0_MR_tree['GalID']==currentgal]
                    currentgal=prog    
                    Nsnaps+=1
                    if prog==-1:
                        break
  
            (x_array,median,mean,pc16,pc84,rms)=convert_2d_array_into_region(FullRedshiftList, Property_Array)             
            subplot.plot(np.log10(x_array+1), np.log10(median), linestyle=plot_linestyle[ii], 
                         linewidth=2, color=plot_color[i_plot])
              
        if(i_plot<2):
            #create an array with the dimension of the previous axis
            major_ticks=np.arange(xlim[0],xlim[1]+0.1,0.2) 
            minor_ticks=np.arange(xlim[0],xlim[1]+0.1,0.1) 
            #interpolate the values of the new axis
            x_array_new=np.interp(major_ticks, np.log10(FullRedshiftList+1), FullTimeList)                
            #one decimal place
            for iii in range(0,len(x_array_new)):
                x_array_new[iii] = int((x_array_new[iii] * 10) + 0.5) / 10.0            
              
            ax2 = subplot.twiny()        
            ax2.set_xbound(subplot.get_xbound())        
            ax2.set_xticks(major_ticks)
            #ax2.set_xticks(minor_ticks, minor = True)
            ax2.set_xticklabels(x_array_new)
                
            xlab='LookBackTime(Gyr)'  
            ax2.set_xlabel(xlab, fontsize=14)
            ax2.xaxis.set_label_position('top')  
            
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.7, y_percentage=0.8, 
                    color='black', xlog=0, ylog=0, label=Property_Name[i_plot], fontsize=13, fontweight='normal')     
   
        if(i_plot==0):            
            plot_label(subplot, 'label', xlim, ylim, x_percentage=0.14, y_percentage=0.26, color='black', 
                       xlog=0, ylog=0, label='$M_{\mathrm{vir(z=0)}}=10^{14}M_{\odot}$', fontsize=10, fontweight='normal') 
            plot_label(subplot, 'line', xlim, ylim, x_percentage=0.12, y_percentage=0.28, 
                        color='black', x2_percentage=0.03, xlog=0, ylog=0, linestyle='-', linewidth=2)
            
            plot_label(subplot, 'label', xlim, ylim, x_percentage=0.14, y_percentage=0.19, color='black', 
                       xlog=0, ylog=0, label='$M_{\mathrm{vir(z=0)}}=10^{13}M_{\odot}$', fontsize=10, fontweight='normal') 
            plot_label(subplot, 'line', xlim, ylim, x_percentage=0.12, y_percentage=0.21, 
                        color='black', x2_percentage=0.03, xlog=0, ylog=0, linestyle='--', linewidth=2)
            
            plot_label(subplot, 'label', xlim, ylim, x_percentage=0.14, y_percentage=0.12, color='black', 
                       xlog=0, ylog=0, label='$M_{\mathrm{vir(z=0)}}=10^{12}M_{\odot}$', fontsize=10, fontweight='normal') 
            plot_label(subplot, 'line', xlim, ylim, x_percentage=0.12, y_percentage=0.14, 
                        color='black', x2_percentage=0.03, xlog=0, ylog=0, linestyle='-.', linewidth=2)
            
            plot_label(subplot, 'label', xlim, ylim, x_percentage=0.14, y_percentage=0.05, color='black', 
                       xlog=0, ylog=0, label='$M_{\mathrm{vir(z=0)}}=10^{11}M_{\odot}$', fontsize=10, fontweight='normal') 
            plot_label(subplot, 'line', xlim, ylim, x_percentage=0.12, y_percentage=0.07, 
                        color='black', x2_percentage=0.03, xlog=0, ylog=0, linestyle=':', linewidth=2)
            
        
       
    plt.tight_layout()   
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()
    
    plt.rcParams.update({'font.size': 14, 'xtick.labelsize': 14, 'ytick.labelsize': 14})    
    
    return 
#halo_growth_normalized    


def baryon_fraction():
       
    file=Datadir+'Database_snapshots_table.fits'
    fits_table=fits.open(file)
    snap_table = fits_table[1]       
    FullRedshiftList=snap_table.data['z'][::-1]
    FullRedshiftList=FullRedshiftList[FullRedshiftList>0.]
      
    #plot_color=['gray','blue','red','green']
    #Property_Name=['$M_{\mathrm{vir}}$','$M_{\mathrm{cold}}$','$M_*$','$SFR$']
    plot_linestyle=[':','-.','--','-']
    
    #SINGLE PLOT
    xlim=[0.0,11.]  
    ylim=[-3.5, -0.5]
            
        
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()  
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    #majorFormatter = FormatStrFormatter('%d')
    #subplot.xaxis.set_major_locator(MultipleLocator(2.0))    
    #subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
    #subplot.yaxis.set_major_locator(MultipleLocator(2.0))    
    #subplot.yaxis.set_minor_locator(MultipleLocator(0.5)) 
   
    ylab='$\log_{10}(M_*/M_{\mathrm{vir}})$'
    xlab='$z$'
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)    
             
    HaloMassBins=np.arange(11.0,15.0,1.0)
    
    for ii in range(0, len(HaloMassBins)):
        
        G0_MR=G_MR[(G_MR['SnapNum']==58) & (G_MR['Type']==0) & 
                   (np.log10(G_MR['Mvir']*1.e10/Hubble_h)>HaloMassBins[ii]-0.1) & 
                   (np.log10(G_MR['Mvir']*1.e10/Hubble_h)<HaloMassBins[ii]+0.1)]  
        
        print('HaloMass=',HaloMassBins[ii],'Ntrees=',len(G0_MR))
       
        Property_Array=np.zeros([len(FullRedshiftList),len(G0_MR)],dtype=np.float32)
      
        #for jj in range(0,len(G0_MR)):
        for jj in range(0,min(10,len(G0_MR))):
            initial_gal=G0_MR['GalID'][jj]   
            G0_MR_tree=G_MR[(G_MR['GalID']>=G0_MR['GalID'][jj]) & (G_MR['GalID']<=G0_MR['LastProgGal'][jj])]    
               
            Nsnaps=0  
            currentgal=initial_gal
            while True:
                sel=G0_MR_tree['GalID']==currentgal
                Property_Array[Nsnaps,jj]=G0_MR_tree['StellarMass'][sel]/G0_MR_tree['Mvir'][sel]
                   
                prog=G0_MR_tree['FirstProgGal'][G0_MR_tree['GalID']==currentgal]
                currentgal=prog    
                Nsnaps+=1
                if prog==-1:
                    break
  
        (x_array,median,mean,pc16,pc84,rms)=convert_2d_array_into_region(FullRedshiftList, Property_Array)
        #subplot.fill_between(x_array,np.log10(pc84),np.log10(pc16),interpolate=True,alpha=0.4,
        #                     facecolor='gray',edgecolor='gray')         
        subplot.plot(x_array,np.log10(median), linestyle=plot_linestyle[ii], linewidth=2, color='black')
            
            
        
    xx=np.arange(xlim[0]-0.5,xlim[1]+0.5,0.1)    
    subplot.plot(xx,xx/xx-1.0+np.log10(0.155), linestyle=':', linewidth=1, color='black')
    
    #plot_label (subplot, 'label', xlim, ylim, x_percentage=0.7, y_percentage=0.8,
    #            color='black', xlog=0, ylog=0, label=Property_Name[i_plot], fontsize=13, fontweight='normal')
   
        
       
    plt.tight_layout()   
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()
  
    return 
#baryon_fraction


def halo_growth_rate():
       
    file=Datadir+'Database_snapshots_table.fits'
    fits_table=fits.open(file)
    snap_table = fits_table[1]       
    FullRedshiftList=snap_table.data['z'][::-1]
    FullTimeList=snap_table.data['lookBackTime'][::-1]
    FullTimeList=FullTimeList[FullRedshiftList>0.]
    FullRedshiftList=FullRedshiftList[FullRedshiftList>0.]   
    
    plot_color=['gray','blue','red','green']
    Property_Name=['$M_{\mathrm{vir}}$','$M_{\mathrm{cold}}$','$M_*$','$SFR$']
    #plot_linestyle=[':','-.','--','-']
    plot_linestyle=[':','--','-']
    
    #SINGLE PLOT
    xlim=[0.0,15.]  
    
    fig = plt.figure(figsize=(two_two_size_small[0],two_two_size_small[1]))
    grid = gridspec.GridSpec(2, 2)
    grid.update(wspace=0.0, hspace=0.0)
    
    #for i_plot in range(0,1):   
    for i_plot in range(0,4):
        if(i_plot==0):
            ylim=[-2.5, 1.0]
            ylab='$\log_{10}(\Delta M_{\mathrm{vir}}/\Delta t[\mathrm{M}_{\odot} \mathrm{Gyr}^{-1]})$'
        else:
            if(i_plot==1):
                ylim=[-2.5, 1.0]
                ylab='$\log_{10}(\Delta M_{\mathrm{Cold}}/\Delta t[\mathrm{M}_{\odot}\mathrm{Gyr}^{-1]}])$'
            else:
                if(i_plot==2):
                    ylim=[-2.5, 1.0]
                    ylab='$\log_{010}(\Delta M_*/\Delta t[\mathrm{M}_{\odot}\mathrm{Gyr}^{-1]}])$'
                else:
                    ylim=[0.0, 2.]
                    ylab='$\log_{10}(SFR[\mathrm{M}_{\odot}\mathrm{yr}^{-1}])$'
                    
        subplot=plt.subplot(grid[i_plot])  
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(2.0))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
        subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
        subplot.yaxis.set_minor_locator(MultipleLocator(0.5)) 
        
        if(i_plot==2 or i_plot==3):
            xlab='Look Back Time (Gyr)'
            subplot.set_xlabel(xlab, fontsize=14)
            
        subplot.set_ylabel(ylab, fontsize=14)    
             
       
  
        if i_plot==1 or i_plot==3:
            subplot.tick_params(axis='y', which='both', left='on', labelleft='off', labelright='on')
            subplot.yaxis.set_label_position("right")

        if i_plot==0 or i_plot==1:    
            subplot.tick_params(axis='x', which='both', bottom='on', labelbottom='off')    
            
            
        #HaloMassBins=np.arange(11.0,15.0,1.0)
        HaloMassBins=np.arange(11.0,14.0,1.0)
    
        for ii in range(0, len(HaloMassBins)):
        
            G0_MR=G_MR[(G_MR['SnapNum']==58) & (G_MR['Type']==0) & 
                       (np.log10(G_MR['Mvir']*1.e10/Hubble_h)>HaloMassBins[ii]-0.1) & 
                       (np.log10(G_MR['Mvir']*1.e10/Hubble_h)<HaloMassBins[ii]+0.1)]  
        
            print('HaloMass=',HaloMassBins[ii],'Ntrees=',len(G0_MR))
       
            Ntrees=min(200,len(G0_MR))
            Property_Array=np.zeros([len(FullTimeList),Ntrees],dtype=np.float32)
            Property_Array_rate=np.zeros([len(FullTimeList),Ntrees],dtype=np.float32)
      
            #for jj in range(0,len(G0_MR)):
            for jj in range(0,Ntrees):
                initial_gal=G0_MR['GalID'][jj]   
                G0_MR_tree=G_MR[(G_MR['GalID']>=G0_MR['GalID'][jj]) & (G_MR['GalID']<=G0_MR['LastProgGal'][jj])]    
                   
                Nsnaps=0  
                currentgal=initial_gal
                while True:
                    #print(Nsnaps)
                    if(Nsnaps==0):
                        delta_t=FullTimeList[1]-FullTimeList[0]
                    else:
                        delta_t=FullTimeList[Nsnaps]-FullTimeList[Nsnaps-1]
                    
                    sel=G0_MR_tree['GalID']==currentgal
                    if(i_plot==0):
                        Property_Array[Nsnaps,jj]=G0_MR_tree['Mvir'][sel]/G0_MR['Mvir'][jj]
                        #print(Property_Array[Nsnaps,jj])
                    else:
                        if(i_plot==1):
                            Property_Array[Nsnaps,jj]=G0_MR_tree['ColdGas'][sel]/G0_MR['ColdGas'][jj]
                        else:
                            if(i_plot==2):
                                Property_Array[Nsnaps,jj]=G0_MR_tree['StellarMass'][sel]/G0_MR['StellarMass'][jj]
                            else:
                                Property_Array[Nsnaps,jj]=np.log10(G0_MR_tree['Sfr'][sel])
                                Property_Array_rate[Nsnaps,jj]=np.log10(G0_MR_tree['Sfr'][sel])
                  
                    if(i_plot<3):
                        if(Nsnaps==0):
                            Property_Array_rate[0,jj]=0.
                        else:
                            Property_Array_rate[Nsnaps,jj]=(Property_Array[Nsnaps-1,jj]-
                                                             Property_Array[Nsnaps,jj])/delta_t
                       
                        #print((Property_Array[Nsnaps-1,jj]-Property_Array[Nsnaps,jj]))
                        
                    prog=G0_MR_tree['FirstProgGal'][G0_MR_tree['GalID']==currentgal]
                    currentgal=prog    
                    Nsnaps+=1
                    if prog==-1:
                        break
  
            (x_array,median,mean,pc16,pc84,rms)=convert_2d_array_into_region(FullTimeList, Property_Array)
            #if(i_plot<2):
            #    print(median)
            '''if(i_plot==0 or i_plot==2):
                subplot.fill_between(x_array,pc84,pc16,interpolate=True,alpha=0.4,
                                     facecolor=plot_color[i_plot],edgecolor=plot_color[i_plot])   '''   
            if(i_plot<3):
                median=np.log10(median)
            subplot.plot(x_array,median, linestyle=plot_linestyle[ii], linewidth=2, color=plot_color[i_plot])
              
    
    
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.7, y_percentage=0.8, 
                    color='black', xlog=0, ylog=0, label=Property_Name[i_plot], fontsize=13, fontweight='normal')     
   
        
       
    plt.tight_layout()   
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()
    
    return 
#halo_growth_rate 



def accretion_history(G_MR, Volume_MR):
   
    file=Datadir+'Database_snapshots_table.fits'
    fits_table=fits.open(file)
    snap_table = fits_table[1]       
    FullRedshiftList=snap_table.data['z'][::-1]
    FullTimeList=snap_table.data['lookBackTime'][::-1]
    
    FullTimeList=FullTimeList[FullRedshiftList>0.]
    FullRedshiftList=FullRedshiftList[FullRedshiftList>0.]
   
    #FullRedshiftList=FullTimeList
    
    plot_linestyle=[':','-.','--','-']

    fig = plt.figure(figsize=(two_one_size_small[0]+1.0,two_one_size_small[1]))
    grid = gridspec.GridSpec(2, 2)
    grid.update(wspace=0.0, hspace=0.0)
    
    '''
    
    #HaloMassBins=np.arange(11.0,15.0,1.0)
    HaloMassBins=np.arange(14.0,15.0,1.0)
    for kk in range (0,len(HaloMassBins)):
    
        G0_MR=G_MR[(G_MR['SnapNum']==58) & (G_MR['Type']==0) & 
                   (np.log10(G_MR['Mvir']*1.e10/Hubble_h)>HaloMassBins[kk]-0.1) & 
                   (np.log10(G_MR['Mvir']*1.e10/Hubble_h)<HaloMassBins[kk]+0.1)]  
                 
        Ntrees=min(1,len(G0_MR))
        
        print('Mass=%0.1f Ntrees=%d selected=%d' % (HaloMassBins[kk], len(G0_MR), Ntrees))
    
        bin=0.5
        bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
        x_axis=bin_arr[0:len(bin_arr)-1]+bin/2.    
        Mvir_hist=np.zeros([len(bin_arr)-1,Ntrees],dtype=np.float32)    
        totalProgMass=np.zeros(Ntrees,dtype=np.float32)    
        FinalMass=np.zeros(Ntrees,dtype=np.float32)    
    
        for ii in range(0,Ntrees):
            initial_gal=G0_MR['GalID'][ii]   
            G0_MR_tree=G_MR[(G_MR['GalID']>=G0_MR['GalID'][ii]) & (G_MR['GalID']<=G0_MR['LastProgGal'][ii])]  
            FinalMass[ii]=np.log10(G0_MR['Mvir'][ii]*1.e10/Hubble_h)
            
            #count number of non main progenitors
            Nprogenitors=0  
            currentgal=initial_gal
            while True:    
                if(len(G0_MR_tree['GalID'][G0_MR_tree['DescendantGal']==currentgal])>1):
                    Nprogenitors+=len(G0_MR_tree['GalID'][G0_MR_tree['DescendantGal']==currentgal])-1                 
                prog=G0_MR_tree['FirstProgGal'][G0_MR_tree['GalID']==currentgal]
                currentgal=prog               
                if prog==-1:
                    break
                
            
            Redshift=np.zeros(Nprogenitors,dtype=np.float32)    
            Mvir=np.zeros(Nprogenitors,dtype=np.float32)
            ColdGas=np.zeros(Nprogenitors,dtype=np.float32)
            StellarMass=np.zeros(Nprogenitors,dtype=np.float32)
            BlackHoleMass=np.zeros(Nprogenitors,dtype=np.float32)
            SFR=np.zeros(Nprogenitors,dtype=np.float32)

            MainGalMvir=np.zeros(len(FullRedshiftList),dtype=np.float32)        
            MainGalColdGas=np.zeros(len(FullRedshiftList),dtype=np.float32)
            MainGalStellarMass=np.zeros(len(FullRedshiftList),dtype=np.float32)
            MainGalBlackHoleMass=np.zeros(len(FullRedshiftList),dtype=np.float32)
            MainGalSFR=np.zeros(len(FullRedshiftList),dtype=np.float32)
        
            Nprogenitors=0  
            Nsnaps=0
            currentgal=initial_gal
            while True:    
                sel=G0_MR_tree['GalID']==currentgal
                MainGalMvir[Nsnaps]=np.log10(G0_MR_tree['Mvir'][sel]*1.e10/Hubble_h)
                MainGalColdGas[Nsnaps]=np.log10(G0_MR_tree['ColdGas'][sel]*1.e10/Hubble_h)
                MainGalStellarMass[Nsnaps]=np.log10(G0_MR_tree['StellarMass'][sel]*1.e10/Hubble_h)            
                MainGalBlackHoleMass[Nsnaps]=np.log10(G0_MR_tree['BlackHoleMass'][sel]*1.e10/Hubble_h)
                #print(G0_MR_tree['Redshift'][sel])
                #print("MainGalMass=%0.5e" % 10**MainGalMvir[Nsnaps])
                if(len(G0_MR_tree['GalID'][G0_MR_tree['DescendantGal']==currentgal])>1):
                    #select all progenitor gals
                    Prog_Gals=G0_MR_tree[(G0_MR_tree['DescendantGal']==currentgal)] 
                    #print(np.log10(Prog_Gals['Mvir']*1.e10/Hubble_h))
                    #sel main progenitor
                    sel=Prog_Gals['GalID']==G0_MR_tree['FirstProgGal'][G0_MR_tree['GalID']==currentgal]
                
                    #sel all non main progenitor gals
                    sel=Prog_Gals['GalID']!=G0_MR_tree['FirstProgGal'][G0_MR_tree['GalID']==currentgal]                    
                    Prog_Gals=Prog_Gals[sel]
                    #print(np.log10(Prog_Gals['Mvir']*1.e10/Hubble_h))
                    progMass=0
                    for jj in range(0,len(Prog_Gals)):
                        Redshift[Nprogenitors+jj]=Prog_Gals['Redshift'][jj]        
                        Mvir[Nprogenitors+jj]=np.log10(Prog_Gals['Mvir'][jj]*1.e10/Hubble_h)
                        ColdGas[Nprogenitors+jj]=np.log10(Prog_Gals['ColdGas'][jj]*1.e10/Hubble_h)
                        StellarMass[Nprogenitors+jj]=np.log10(Prog_Gals['StellarMass'][jj]*1.e10/Hubble_h)
                        BlackHoleMass[Nprogenitors+jj]=np.log10(Prog_Gals['BlackHoleMass'][jj]*1.e10/Hubble_h)
                        SFR[Nprogenitors+jj]=np.log10(Prog_Gals['Sfr'][jj])
                        progMass+=Prog_Gals['Mvir'][jj]*1.e10/Hubble_h
                    
                    #Nprogenitors+=len(G0_MR_tree['GalID'][G0_MR_tree['DescendantGal']==currentgal])-1 
                    Nprogenitors+=len(Prog_Gals)
                #print("progMass=%0.5e" % progMass)
                #print("ratio=%0.5f\n" % (progMass/10**MainGalMvir[Nsnaps]))
                
                prog=G0_MR_tree['FirstProgGal'][G0_MR_tree['GalID']==currentgal]
                currentgal=prog    
                Nsnaps+=1  
                if prog==-1:
                    break
                    
        
            hist_MR=np.histogram(Mvir, bins=bin_arr, range=(xlim[0],xlim[1]))          
            hist_MR=hist_MR[0]         
            #Mvir_hist[:,ii]=hist_MR/len(Mvir)
            Mvir_hist[:,ii]=hist_MR/len(Mvir)#*(10**x_axis/1.e14)
            totalProgMass[ii]=np.sum(10**Mvir,axis=0)        
        
            #print("tot prog mass=%0.1f" % totalProgMass[ii])
        
        (x_array,median,mean,pc16,pc84,rms)=convert_2d_array_into_region(x_axis,Mvir_hist)
        subplot.fill_between(x_array,np.log10(pc84),np.log10(pc16),interpolate=True,alpha=0.4,
                             facecolor='gray',edgecolor='gray')         
        subplot.plot(x_array,np.log10(median), linestyle=plot_linestyle[kk], linewidth=2, color='gray')
    
        #print(np.median(totalProgMass/10**(FinalMass)))
        
        #subplot.plot(x_axis,Mvir_hist[:,0], color='red', linewidth=2, linestyle=':', drawstyle='steps') '''
         
    subplot=plt.subplot()  
    #xlim=[8.5,15.]
    xlim=[10.5,15.]
    xlim=[-5.0,0.]
    ylim=[-4.0, 0.0]
    
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    #majorFormatter = FormatStrFormatter('%d')
    #subplot.xaxis.set_major_locator(MultipleLocator(2.0))    
    #subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
    #subplot.yaxis.set_major_locator(MultipleLocator(2.0))    
    #subplot.yaxis.set_minor_locator(MultipleLocator(0.5)) 
   
    xlab='$\log_{10}(M_{\mathrm{vir}}[\mathrm{M}_{\odot}])$'
    ylab='Fraction'          
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)    
    
    #11.0,12.0,13.0,14.0
    HaloMassBins=np.arange(11.0,15.0,1.0)
    HaloMassBins=np.arange(12.0,15.0,1.0)    
    for kk in range (0,len(HaloMassBins)):   
        G0_MR=G_MR[(G_MR['SnapNum']==58) & (G_MR['Type']==0) & 
                   (np.log10(G_MR['Mvir']*1.e10/Hubble_h)>HaloMassBins[kk]-0.1) & 
                   (np.log10(G_MR['Mvir']*1.e10/Hubble_h)<HaloMassBins[kk]+0.1)]  
                 
        Ntrees=min(100,len(G0_MR))
        
       
        '''trees=[]
        #find trees that have more than a single halo
        for ii in range(0,Ntrees):
            initial_gal=G0_MR['GalID'][ii]   
            G0_MR_tree=G_MR[(G_MR['GalID']>=G0_MR['GalID'][ii]) & (G_MR['GalID']<=G0_MR['LastProgGal'][ii])]  
                      
            GalsInTree = len(G0_MR_tree)
            SnapsInTree = (1+G0_MR['SnapNum'][ii]-G_MR['SnapNum'][G_MR['GalID']==G0_MR['LastProgGal'][ii]])
            #print((GalsInTree-1),SnapsInTree)
            if((GalsInTree-1)<=SnapsInTree):
                continue
            else:                
                trees.append(ii)
                
        Ntrees = len(trees)       
        print('Mass=%0.1f Ntrees=%d selected=%d' % (HaloMassBins[kk], len(G0_MR), Ntrees))'''
        
        bin=0.2
        limits = [10.0,15.0]
        bin_arr = np.arange(limits[0],limits[1]+bin,bin)
        x_axis = bin_arr[0:len(bin_arr)-1]+bin/2.    
        Mvir_hist=np.zeros([len(bin_arr)-1,Ntrees],dtype=np.float32)    
        totalProgMass=np.zeros(Ntrees,dtype=np.float32)    
        FinalMass=np.zeros(Ntrees,dtype=np.float32)    
    
        for ll in range(0,Ntrees):
            #ii = trees[ll]
            ii = ll
            initial_gal=G0_MR['GalID'][ii]   
            G0_MR_tree=G_MR[(G_MR['GalID']>=G0_MR['GalID'][ii]) & (G_MR['GalID']<=G0_MR['LastProgGal'][ii])]  
            FinalMass[ll]=np.log10(G0_MR['Mvir'][ii]*1.e10/Hubble_h)
            
            #skip tree with a single galaxy
            GalsInTree = len(G0_MR_tree)
            SnapsInTree = (1+G0_MR['SnapNum'][ii]-G_MR['SnapNum'][G_MR['GalID']==G0_MR['LastProgGal'][ii]])            
            if((GalsInTree)-1<=SnapsInTree):
                continue
           
        
            #count number of non main progenitors
            Nprogenitors=0  
            currentgal=initial_gal
            while True:    
                descendants = G0_MR_tree['DescendantGal']==currentgal
                if(len(G0_MR_tree['GalID'][descendants])>1):
                    Nprogenitors+=(len(G0_MR_tree['GalID'][descendants])-1)                
                prog=G0_MR_tree['FirstProgGal'][G0_MR_tree['GalID']==currentgal]
                currentgal=prog               
                if prog==-1:
                    break
            #print("Nprogenitors: %d" %(Nprogenitors))    
            
            Redshift=np.zeros(Nprogenitors,dtype=np.float32)    
            Mvir=np.zeros(Nprogenitors,dtype=np.float32)           
            MainGalMvir=np.zeros(len(FullRedshiftList),dtype=np.float32)        
            
        
            Nprogenitors=0  
            Nsnaps=0
            currentgal=initial_gal
            while True:    
                sel=G0_MR_tree['GalID']==currentgal
                MainGalMvir[Nsnaps]=np.log10(G0_MR_tree['Mvir'][sel]*1.e10/Hubble_h)
             
                #select all progenitor gals
                G_Prog_Gals=G0_MR_tree[(G0_MR_tree['DescendantGal']==currentgal)] 
                #select main progenitor
                FirstProgID = G0_MR_tree['FirstProgGal'][G0_MR_tree['GalID']==currentgal]  
                G_FirstProg = G0_MR_tree[G0_MR_tree['GalID']==FirstProgID]  
                    
                #CHANGE HERE
                #LOOP THROUGH ALL THE PROGENITORS TO FIND WHEN TEHY MERGED
                           
                #sel all non main progenitor gals
                FirstProgID = G0_MR_tree['FirstProgGal'][G0_MR_tree['GalID']==currentgal]   
                sel=G_Prog_Gals['GalID']!= FirstProgID                
                G_Prog_Gals=G_Prog_Gals[sel]
                #print(np.log10(Prog_Gals['Mvir']*1.e10/Hubble_h))
                progMass=0
                for jj in range(0,len(G_Prog_Gals)):
                    Redshift[Nprogenitors+jj]=G_Prog_Gals['Redshift'][jj]  
                    Mvir[Nprogenitors+jj]=np.log10(G_Prog_Gals['Mvir'][jj]*1.e10/Hubble_h)                       
                     
                    progMass+=G_Prog_Gals['Mvir'][jj]*1.e10/Hubble_h
                    
                #Nprogenitors+=len(G0_MR_tree['GalID'][G0_MR_tree['DescendantGal']==currentgal])-1 
                Nprogenitors+=len(G_Prog_Gals)
               
                #stop if the main progenitor of current galaxy has Mvir<0.1 x final_Mvir
                if(np.log10(G_FirstProg['Mvir']*1.e10/Hubble_h) < (FinalMass[ll]-1.)):  
                    break  
            
                #stop if there are no more progenitors
                prog=G0_MR_tree['FirstProgGal'][G0_MR_tree['GalID']==currentgal]
                currentgal=prog 
                if(prog==-1):
                    break
                Nsnaps+=1  
                
                
                
                  
                
                #if prog==-1:
                #    break
                    
            #print(Mvir)
            hist_MR=np.histogram(Mvir, bins=bin_arr, range=(xlim[0],xlim[1]))          
            hist_MR=hist_MR[0]         
            #divide by bin value over (final_mass-0.1*final_mass) to get the percetange fraction of mass
            #in each bin
            Mass_dif = 10**FinalMass[ll]-10**(FinalMass[ll]-1.)
            Mvir_hist[:,ll]=hist_MR * (10**x_axis)/Mass_dif          
            totalProgMass[ll]=np.log10(np.sum(10**Mvir,axis=0))        
        
        #for ll in range(0,Ntrees):
            #print(Mvir_hist[:,ii])
            #print(10**totalProgMass[ll]/10**FinalMass[ll])
        
        #CHANGE HERE
        final_hist_mean=np.zeros(len(x_axis),dtype=np.float32)   
        final_hist_median=np.zeros(len(x_axis),dtype=np.float32)   
        for ii in range(0, len(x_axis)):
            final_hist_mean[ii] = np.mean(Mvir_hist[ii,:])
            final_hist_median[ii] = np.median(Mvir_hist[ii,:])
       
        subplot.plot(x_axis-HaloMassBins[kk], np.log10(final_hist_mean), linestyle=plot_linestyle[kk], linewidth=2, color='gray')
        #subplot.plot(x_axis, np.log10(final_hist_median), linestyle=plot_linestyle[kk], linewidth=2, color='gray')
        print("Mass in mergers=%0.2f" % sum(final_hist_mean))
        
       
        #subplot.plot(np.log10((10**x_axis/10**HaloMassBins[kk])), np.log10(final_hist_mean), 
        #             linestyle=plot_linestyle[kk], linewidth=2, color='gray')
         
    plt.tight_layout()   
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()        
    
    plot_prog_at_single_redshift=0
    
    if(plot_prog_at_single_redshift==1):
            
        subplot=plt.subplot()   
        xlim=[-5.0,0.]   
        ylim=[-3., 0.0]
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
        xlab='$\log_{10}(M_{\mathrm{vir}}[\mathrm{M}_{\odot}])$'
        ylab='Fraction'          
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)    
    
        #11.0,12.0,13.0,14.0   
        HaloMassBins=np.arange(12.0,15.0,1.0)    
        
        for kk in range (0,len(HaloMassBins)):
    
            G0_MR=G_MR[(G_MR['SnapNum']==58) & (G_MR['Type']==0) & 
                       (np.log10(G_MR['Mvir']*1.e10/Hubble_h)>HaloMassBins[kk]-0.1) & 
                       (np.log10(G_MR['Mvir']*1.e10/Hubble_h)<HaloMassBins[kk]+0.1)]  
                 
            Ntrees=min(10,len(G0_MR))
            print(Ntrees)
       
          
            totalProgMass=np.zeros(Ntrees,dtype=np.float32)  
            target_redshift=np.zeros(Ntrees,dtype=np.float32)    
            target_mass=np.zeros(Ntrees,dtype=np.float32)    
            FinalMass=np.zeros(Ntrees,dtype=np.float32)    
    
            for ll in range(0,Ntrees):           
                initial_galID = G0_MR['GalID'][ll]   
                FinalMass[ll] = np.log10(G0_MR['Mvir'][ll]*1.e10/Hubble_h)
                G_tree = G_MR[(G_MR['GalID']>=G0_MR['GalID'][ll]) & (G_MR['GalID']<=G0_MR['LastProgGal'][ll])]  
              
                currentgalID = initial_galID
                while True:    
                
                    FirstProgID = G_tree['FirstProgGal'][G_tree['GalID']==currentgalID]
                    G_FirstProg = G_tree[G_tree['GalID']==FirstProgID]
                
                    #if main progenitor Mvir = 0.1 x Mvir(z=0)   
                    if(FirstProgID>-1):                    
                        if(np.log10(G_FirstProg['Mvir']*1.e10/Hubble_h) < (FinalMass[ll]-1.)):                       
                            target_redshift[ll] = G_FirstProg['Redshift']
                            target_mass[ll] = np.log10(G_FirstProg['Mvir']*1.e10/Hubble_h)
                            break  
                    #if there are no more progenitors, break        
                    else:
                        break
               
                    currentgalID=FirstProgID    
                    Nsnaps+=1  
                
            bin=0.2
            hist_lim = [10.0,15.0]   
            bin_arr=np.arange(hist_lim[0],hist_lim[1]+bin,bin)
            Mvir_hist=np.zeros([len(bin_arr)-1,Ntrees],dtype=np.float32)       
        
            for ll in range(0,Ntrees):
                initial_galID = G0_MR['GalID'][ll]   
           
                G_target_snap = G_MR[(G_MR['GalID']>=G0_MR['GalID'][ll]) & 
                                     (G_MR['GalID']<=G0_MR['LastProgGal'][ll]) &
                                     (G_MR['Type']<2) & (G_MR['Redshift']==target_redshift[ll])]  
            
            
                hist_MR=np.histogram(np.log10(G_target_snap['Mvir']*1.e10/Hubble_h), bins=bin_arr)      
                x_axis = hist_MR[1][0:len(hist_MR[1][:])-1]+bin/2.               
                Mvir_hist[:,ll] = hist_MR[0]*((10**x_axis)/(10**FinalMass[ll]-10**target_mass[ll]))
           
        
            final_hist_mean=np.zeros(len(bin_arr)-1,dtype=np.float32)   
            final_hist_median=np.zeros(len(bin_arr)-1,dtype=np.float32)   
            for ii in range(0, len(x_axis)):
                final_hist_mean[ii] = np.mean(Mvir_hist[ii,:])
                final_hist_median[ii] = np.median(Mvir_hist[ii,:])
          
        
            subplot.plot(np.log10((10**x_axis/10**HaloMassBins[kk])),np.log10(final_hist_mean), linestyle=plot_linestyle[kk], linewidth=2, color='gray')
       
     
        plt.tight_layout()   
        current_function =  inspect.getframeinfo(inspect.currentframe()).function   
        plt.savefig('./fig/plots_'+current_function+'_2.pdf')
        plt.close()
    
    
    #QUERIES FOR HALOS
    
    '''select top 10 log10(M_crit200*1.e10/0.67), 
                  HaloID        
    from MPAHaloTrees..MrscPlanck1        
    where M_crit200>10000. and snapnum=58'''
    
    '''HaloID    LAstProg
    84000007005505	84000007008523'''
    
    return 
#accretion_history

def cooling_heating_evo():
   
    #constants
    UnitLength_in_cm=3.08568e+24
    UnitVelocity_in_cm_per_s=100000.
    UnitTime_in_s=UnitLength_in_cm / UnitVelocity_in_cm_per_s
    SOLAR_MASS=1.989e33
    UnitMass_in_g=1.989e+43
    SEC_PER_YEAR=3.155e7
    UnitEnergy_in_cgs = UnitMass_in_g * UnitLength_in_cm**2. / UnitTime_in_s**2.

    
    #sel=log_SSFR>SSFR_cut[ii]
    #subplot.scatter(x_axis[sel], y_axis[sel], s=5, color='blue')           
        
    xlim=[9.0,12.0]
    ylim=[-3.0, 3.0]
    #ylim=[-4.0, 0.0]
    bin=0.25
    
    model = 'Hen15_test'
    
    compute_and_write = 1
    
    #snap_list = np.array([62, 42, 34, 29])
    #snap_list = np.array([58, 38, 30, 25])
    snap_list = np.array([58, 38, 30])
    SSFR_cut=[-11.,-11.0, -10.5,-7.5]
    
    
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1])) 
    subplot=plt.subplot()    
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1.0))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
    subplot.yaxis.set_major_locator(MultipleLocator(2.0))    
    subplot.yaxis.set_minor_locator(MultipleLocator(0.5))     
    xlab='$\log_{10}(M_*)$'; ylab='$\log_{10}(Heating/Cooling)$' 
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)
    
    #Ngals=10000
    Ngals=100
    
    transition_mass = np.zeros(len(snap_list)*Ngals,dtype=np.float32)
    transition_snap = np.zeros(len(snap_list)*Ngals,dtype=np.float32)
    
    for idx, element in enumerate(snap_list):
        
        char_snap ="%d" % element
                
        if(compute_and_write == 1):            
            log_StellarMass=(np.log10(G_MR['StellarMass']*1.e10/Hubble_h))                    
            Cooling=G_MR['CoolingRate_beforeAGN'] * (UnitTime_in_s/SEC_PER_YEAR)*(SOLAR_MASS/UnitMass_in_g)
            AgnEfficiency=5.3e-3*(UnitTime_in_s*SOLAR_MASS)/(UnitMass_in_g*SEC_PER_YEAR)
            AGNcoeff = (1.34e5**2 / G_MR['Vvir']**2)           
            AGNrate= AgnEfficiency * G_MR['BlackHoleMass']/Hubble_h* (G_MR['HotGas']/Hubble_h) * 10.         
            EDDrate = 1.3e48 * G_MR['BlackHoleMass'] / (UnitEnergy_in_cgs / UnitTime_in_s) / 9e10; 
            
            AGNrate[AGNrate>EDDrate]=EDDrate[AGNrate>EDDrate]              
            AGNheating = AGNcoeff * AGNrate  
            
            sel = ( (G_MR['SnapNum']==element) & (G_MR['Type']==0) & (AGNheating>Cooling) & 
                    (np.log10(G_MR['StellarMass']*1.e10/Hubble_h)>9.5) )
            #(np.log10(G_MR['Sfr']/(G_MR['StellarMass']*1.e10/Hubble_h)) < SSFR_cut[idx]) )   
                    
            
            #log_SSFR=
            total_ngals=len(G_MR[sel])
            #Nsnaps=64
            Nsnaps=58
          
            if(Ngals>total_ngals):
                Ngals=total_ngals
            print('Snap %d, Ngals=%d' % (element, Ngals))       
          
            G0_MR=np.random.choice(G_MR[sel], size=Ngals, replace=False) 
                         
            Redshift=np.zeros(Nsnaps*Ngals,dtype=np.float32)    
            Mvir=np.zeros(Nsnaps*Ngals,dtype=np.float32)           
            HotGas=np.zeros(Nsnaps*Ngals,dtype=np.float32)          
            StellarMass=np.zeros(Nsnaps*Ngals,dtype=np.float32)
            BlackHoleMass=np.zeros(Nsnaps*Ngals,dtype=np.float32)
            cooling_heating = np.zeros(Nsnaps*Ngals,dtype=np.float32)
          
           
            
            for ii in range (0,Ngals): 
                initial_gal=G0_MR['GalID'][ii]
                G0_MR_tree=G_MR[(G_MR['GalID']>=G0_MR['GalID'][ii]) & (G_MR['GalID']<=G0_MR['LastProgGal'][ii]) 
                               & (G_MR['Type']<2)]    
      
                snap=0  
                currentgal=initial_gal
                transition_flag = 0
                while True:                     
                    Gal = G0_MR_tree[G0_MR_tree['GalID']==currentgal]                    
                    Redshift[Nsnaps*ii+snap] = Gal['Redshift']        
                    Mvir[Nsnaps*ii+snap] = np.log10(Gal['Mvir']*1.e10/Hubble_h)                     
                    HotGas[Nsnaps*ii+snap]=np.log10(Gal['HotGas']*1.e10/Hubble_h)                   
                    StellarMass[Nsnaps*ii+snap]=np.log10(Gal['StellarMass']*1.e10/Hubble_h)
                    BlackHoleMass[Nsnaps*ii+snap]=np.log10(Gal['BlackHoleMass']*1.e10/Hubble_h)
                   
           
                    log_StellarMass=(np.log10(Gal['StellarMass']*1.e10/Hubble_h))                    
                    Cooling=Gal['CoolingRate_beforeAGN'] * (UnitTime_in_s/SEC_PER_YEAR)*(SOLAR_MASS/UnitMass_in_g)
                    AgnEfficiency=5.3e-3*(UnitTime_in_s*SOLAR_MASS)/(UnitMass_in_g*SEC_PER_YEAR)
                    AGNcoeff = (1.34e5 / Gal['Vvir']) * (1.34e5 / Gal['Vvir'])           
                    AGNrate= AgnEfficiency * Gal['BlackHoleMass']/Hubble_h* (Gal['HotGas']/Hubble_h) * 10.         
                    EDDrate = 1.3e48 * Gal['BlackHoleMass'] / (UnitEnergy_in_cgs / UnitTime_in_s) / 9e10;           
                    if(AGNrate>EDDrate):
                        AGNrate=EDDrate              
                    AGNheating = AGNcoeff * AGNrate  
                    if(AGNheating/Cooling == 0.):
                        cooling_heating[Nsnaps*ii+snap]=-99.
                    else:   
                        cooling_heating[Nsnaps*ii+snap]=np.log10(AGNheating/Cooling)
                    
                    if((AGNheating<Cooling) & (transition_flag==0)):
                        transition_mass[idx*ii+ii] = log_StellarMass
                        transition_snap[idx*ii+ii] = Gal['SnapNum'] 
                        transition_flag = 1
                
        
                    prog=Gal['FirstProgGal']
                    currentgal=prog    
                    snap+=1
                    if prog==-1:                       
                        break
                    
                if(Ngals>10):
                    if(ii%(Ngals/10)==0):
                        print(ii, Ngals)           
        
        
            #subplot.scatter(StellarMass, cooling_heating, s=5, color='red') 
            #write to file            
            bin = 0.2
            sel = cooling_heating>-99.
            (x_binned, median, mean, pc16, pc84, rms)= median_and_percentiles(bin, xlim[0], max(StellarMass)+bin/2., 
                                                                              StellarMass[sel], cooling_heating[sel]) 
            sel = median!=0.
            median = median[sel]
            x_binned = x_binned[sel]
            pc16 = pc16[sel]
            pc84 = pc84[sel]
            #(x_binned, median, mean, pc16, pc84, rms)= median_and_percentiles_fixed_xx(StellarMass, cooling_heating) 
            fa = open(Datadir+"cooling_heating_evo_"+model+'_snap_'+char_snap+".txt", "w")
            fa.write("%d\n" % len(x_binned))
            for kk in range (0,len(x_binned)):               
                fa.write("%0.2f " % x_binned[kk] + "%0.2f " % median[kk]  + "%0.2f " % pc16[kk]  + "%0.2f\n" % pc84[kk])   
            fa.close() 
            
            #print(np.mean(transition_mass))
            median = np.median(transition_mass)
            y_sorted = np.sort(transition_mass)
            pc16 = y_sorted[int(16*len(transition_mass)/100)]      
            pc84 = y_sorted[int(84*len(transition_mass)/100)]  
            fa = open(Datadir+"cooling_heating_evo_transition_mass_"+model+'_snap_'+char_snap+".txt", "w")           
            fa.write("%0.5f\n" % median)
            fa.write("%0.5f\n" % pc16)
            fa.write("%0.5f\n" % pc84)
            fa.close() 
           
        
           
        fa = Datadir+"cooling_heating_evo_"+model+'_snap_'+char_snap+".txt"
        (x_binned,median, pc16, pc84)=read_data_with_err(fa)  
        subplot.fill_between(x_binned, pc84, pc16, interpolate=True,alpha=0.8, facecolor='grey', edgecolor='grey') 
        subplot.plot(x_binned, median, linestyle='-', linewidth=2, color='black', alpha=0.8)
        #subplot.plot(ran_gal_z, ran_gal_Mvir, linestyle='-', linewidth=2, color='black')
        
        
        if(idx==0):
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.15, y_percentage=0.93, color='black', xlog=0, ylog=0,
                        label='$M_{200\mathrm{c}}$', fontsize=10)    
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.05, y_percentage=0.945, 
                        color='grey', x2_percentage=0.12, xlog=0, ylog=0, linestyle='-', linewidth=2) 
           
        label='$M_{\mathrm{200c,0}}\sim$%0.1f' % element
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.62, y_percentage=0.9, color='black', 
                    xlog=0, ylog=0, label=label, fontsize=10)    
    
    fa = Datadir+"cooling_heating_evo_transition_mass_"+model+'_snap_'+char_snap+".txt"
    fa = open(fa, "r")                            
    fields = list(fa)        
    median = float(fields[0].strip())
    pc16 = float(fields[1].strip())
    pc84 = float(fields[2].strip())
    #print(median,pc16,pc84)
    
    snap_list = np.array([61, 51, 42, 34, 29])
    snap_list = np.array([57, 47, 38, 30, 25])
    #snap_list = np.array([42, 29])
    for idx, element in enumerate(snap_list):  
        char_snap ="%d" % element
            
        sel = ((transition_snap > element-2) & ((transition_snap < element+2)))
        print(element, len(transition_mass[sel]))    
        print(np.mean(transition_mass[sel]))      
        median = np.median(transition_mass[sel])
        y_sorted = np.sort(transition_mass[sel])
        pc16 = y_sorted[int(16*len(transition_mass[sel])/100)]      
        pc84 = y_sorted[int(84*len(transition_mass[sel])/100)]  
        fa = open(Datadir+"cooling_heating_evo_transition_mass_at_snap_"+model+'_snap_'+char_snap+".txt", "w")           
        fa.write("%0.5f\n" % median)
        fa.write("%0.5f\n" % pc16)
        fa.write("%0.5f\n" % pc84)
        fa.close() 
            
        fa = Datadir+"cooling_heating_evo_transition_mass_at_snap_"+model+'_snap_'+char_snap+".txt"
        fa = open(fa, "r")                            
        fields = list(fa)        
        median = float(fields[0].strip())
        pc16 = float(fields[1].strip())
        pc84 = float(fields[2].strip())
        print(median,pc16,pc84)
    
    '''(61, 62)  10.1122 (10.05583, 9.70537, 10.50799)
       (51, 129) 10.2605 (10.31141, 9.89888, 10.64233)
       (42, 162) 10.3772 (10.42888, 10.11866, 10.70417)
       (34, 204) 10.5206 (10.52698, 10.29112, 10.81783)
       (29, 27)  10.5599 (10.46732, 10.42009, 10.79474)

    '''
    


    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HWL18_cooling_heating_evo.pdf')
    plt.close()
    
    return    
#end  cooling_heating_evo  


def cooling_heating_evo_halo_mass():
   
    #constants
    UnitLength_in_cm=3.08568e+24
    UnitVelocity_in_cm_per_s=100000.
    UnitTime_in_s=UnitLength_in_cm / UnitVelocity_in_cm_per_s
    SOLAR_MASS=1.989e33
    UnitMass_in_g=1.989e+43
    SEC_PER_YEAR=3.155e7
    UnitEnergy_in_cgs = UnitMass_in_g * UnitLength_in_cm**2. / UnitTime_in_s**2.

    
    #sel=log_SSFR>SSFR_cut[ii]
    #subplot.scatter(x_axis[sel], y_axis[sel], s=5, color='blue')           
        
    xlim=[11.0,15.0]
    ylim=[-3.0, 3.0]
    #ylim=[-4.0, 0.0]
    bin=0.25
    
    model = 'Hen15'
    
    compute_and_write = 1
    
    #snap_list = np.array([62, 42, 34, 29])
    #snap_list = np.array([58, 38, 30, 25])
    snap_list = np.array([58, 38, 30])
    SSFR_cut=[-11.,-11.0, -10.5,-7.5]
    
    
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1])) 
    subplot=plt.subplot()    
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1.0))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
    subplot.yaxis.set_major_locator(MultipleLocator(2.0))    
    subplot.yaxis.set_minor_locator(MultipleLocator(0.5))     
    xlab='$\log_{10}(M_*)$'; ylab='$\log_{10}(Heating/Cooling)$' 
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)
    
    #Ngals=10000
    Ngals=100
    
    transition_mass = np.zeros(len(snap_list)*Ngals,dtype=np.float32)
    transition_snap = np.zeros(len(snap_list)*Ngals,dtype=np.float32)
    
    for idx, element in enumerate(snap_list):
        
        char_snap ="%d" % element
                
        if(compute_and_write == 1):            
            log_StellarMass=(np.log10(G_MR['StellarMass']*1.e10/Hubble_h))         
            log_Mvir=(np.log10(G_MR['Mvir']*1.e10/Hubble_h))           
            Cooling=G_MR['CoolingRate_beforeAGN'] * (UnitTime_in_s/SEC_PER_YEAR)*(SOLAR_MASS/UnitMass_in_g)
            AgnEfficiency=5.3e-3*(UnitTime_in_s*SOLAR_MASS)/(UnitMass_in_g*SEC_PER_YEAR)
            AGNcoeff = (1.34e5**2 / G_MR['Vvir']**2)           
            AGNrate= AgnEfficiency * G_MR['BlackHoleMass']/Hubble_h* (G_MR['HotGas']/Hubble_h) * 10.         
            EDDrate = 1.3e48 * G_MR['BlackHoleMass'] / (UnitEnergy_in_cgs / UnitTime_in_s) / 9e10; 
            
            AGNrate[AGNrate>EDDrate]=EDDrate[AGNrate>EDDrate]              
            AGNheating = AGNcoeff * AGNrate  
            
            sel = ( (G_MR['SnapNum']==element) & (G_MR['Type']==0) & (AGNheating>Cooling) & 
                    (np.log10(G_MR['Mvir']*1.e10/Hubble_h)>10.5) )
            #(np.log10(G_MR['Sfr']/(G_MR['StellarMass']*1.e10/Hubble_h)) < SSFR_cut[idx]) )   
                    
            
            #log_SSFR=
            total_ngals=len(G_MR[sel])
            #Nsnaps=64
            Nsnaps=58
          
            if(Ngals>total_ngals):
                Ngals=total_ngals
            print('Snap %d, Ngals=%d' % (element, Ngals))       
          
            G0_MR=np.random.choice(G_MR[sel], size=Ngals, replace=False) 
                         
            Redshift=np.zeros(Nsnaps*Ngals,dtype=np.float32)    
            Mvir=np.zeros(Nsnaps*Ngals,dtype=np.float32)           
            HotGas=np.zeros(Nsnaps*Ngals,dtype=np.float32)          
            StellarMass=np.zeros(Nsnaps*Ngals,dtype=np.float32)
            BlackHoleMass=np.zeros(Nsnaps*Ngals,dtype=np.float32)
            cooling_heating = np.zeros(Nsnaps*Ngals,dtype=np.float32)
          
           
            
            for ii in range (0,Ngals): 
                initial_gal=G0_MR['GalID'][ii]
                G0_MR_tree=G_MR[(G_MR['GalID']>=G0_MR['GalID'][ii]) & (G_MR['GalID']<=G0_MR['LastProgGal'][ii]) 
                               & (G_MR['Type']<2)]    
      
                snap=0  
                currentgal=initial_gal
                transition_flag = 0
                while True:                     
                    Gal = G0_MR_tree[G0_MR_tree['GalID']==currentgal]                     
                    Redshift[Nsnaps*ii+snap] = Gal['Redshift']        
                    Mvir[Nsnaps*ii+snap] = np.log10(Gal['Mvir']*1.e10/Hubble_h)  
                    
                    HotGas[Nsnaps*ii+snap]=np.log10(Gal['HotGas']*1.e10/Hubble_h)                   
                    StellarMass[Nsnaps*ii+snap]=np.log10(Gal['StellarMass']*1.e10/Hubble_h)
                    BlackHoleMass[Nsnaps*ii+snap]=np.log10(Gal['BlackHoleMass']*1.e10/Hubble_h)
                   
                    log_Mvir=(np.log10(Gal['Mvir']*1.e10/Hubble_h))           
                    log_StellarMass=(np.log10(Gal['StellarMass']*1.e10/Hubble_h))                    
                    Cooling=Gal['CoolingRate_beforeAGN'] * (UnitTime_in_s/SEC_PER_YEAR)*(SOLAR_MASS/UnitMass_in_g)
                    AgnEfficiency=5.3e-3*(UnitTime_in_s*SOLAR_MASS)/(UnitMass_in_g*SEC_PER_YEAR)
                    AGNcoeff = (1.34e5 / Gal['Vvir']) * (1.34e5 / Gal['Vvir'])           
                    AGNrate= AgnEfficiency * Gal['BlackHoleMass']/Hubble_h* (Gal['HotGas']/Hubble_h) * 10.         
                    EDDrate = 1.3e48 * Gal['BlackHoleMass'] / (UnitEnergy_in_cgs / UnitTime_in_s) / 9e10;           
                    if(AGNrate>EDDrate):
                        AGNrate=EDDrate              
                    AGNheating = AGNcoeff * AGNrate  
                    if(AGNheating/Cooling == 0.):
                        cooling_heating[Nsnaps*ii+snap]=-99.
                    else:   
                        cooling_heating[Nsnaps*ii+snap]=np.log10(AGNheating/Cooling)
                    
                    if((AGNheating<Cooling) & (transition_flag==0)):
                        transition_mass[idx*ii+ii] = log_Mvir
                        transition_snap[idx*ii+ii] = Gal['SnapNum'] 
                        transition_flag = 1
                
        
                    prog=Gal['FirstProgGal']
                    currentgal=prog    
                    snap+=1
                    if prog==-1:                       
                        break
                    
                if(Ngals>10):
                    if(ii%(Ngals/10)==0):
                        print(ii, Ngals)           
        
        
            #subplot.scatter(StellarMass, cooling_heating, s=5, color='red') 
            #write to file            
            bin = 0.2
            sel = cooling_heating>-99.
            (x_binned, median, mean, pc16, pc84, rms)= median_and_percentiles(bin, xlim[0], max(Mvir)+bin/2., 
                                                                              Mvir[sel], cooling_heating[sel]) 
            sel = median!=0.
            median = median[sel]
            x_binned = x_binned[sel]
            pc16 = pc16[sel]
            pc84 = pc84[sel]
            #(x_binned, median, mean, pc16, pc84, rms)= median_and_percentiles_fixed_xx(StellarMass, cooling_heating) 
            fa = open(Datadir+"cooling_heating_evo_halo_mass_"+model+'_snap_'+char_snap+".txt", "w")
            fa.write("%d\n" % len(x_binned))
            for kk in range (0,len(x_binned)):               
                fa.write("%0.2f " % x_binned[kk] + "%0.2f " % median[kk]  + "%0.2f " % pc16[kk]  + "%0.2f\n" % pc84[kk])   
            fa.close() 
            
            #print(np.mean(transition_mass))
            median = np.median(transition_mass)
            y_sorted = np.sort(transition_mass)
            pc16 = y_sorted[int(16*len(transition_mass)/100)]      
            pc84 = y_sorted[int(84*len(transition_mass)/100)]  
            fa = open(Datadir+"cooling_heating_evo_halo_mass_transition_mass_"+model+'_snap_'+char_snap+".txt", "w")           
            fa.write("%0.5f\n" % median)
            fa.write("%0.5f\n" % pc16)
            fa.write("%0.5f\n" % pc84)
            fa.close() 
           
        
           
        fa = Datadir+"cooling_heating_evo_halo_mass_"+model+'_snap_'+char_snap+".txt"
        (x_binned,median, pc16, pc84)=read_data_with_err(fa)  
        subplot.fill_between(x_binned, pc84, pc16, interpolate=True,alpha=0.8, facecolor='grey', edgecolor='grey') 
        subplot.plot(x_binned, median, linestyle='-', linewidth=2, color='black', alpha=0.8)
        #subplot.plot(ran_gal_z, ran_gal_Mvir, linestyle='-', linewidth=2, color='black')
        
        
        if(idx==0):
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.15, y_percentage=0.93, color='black', xlog=0, ylog=0,
                        label='$M_{200\mathrm{c}}$', fontsize=10)    
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.05, y_percentage=0.945, 
                        color='grey', x2_percentage=0.12, xlog=0, ylog=0, linestyle='-', linewidth=2) 
           
        label='$M_{\mathrm{200c,0}}\sim$%0.1f' % element
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.62, y_percentage=0.9, color='black', 
                    xlog=0, ylog=0, label=label, fontsize=10)    
    
    fa = Datadir+"cooling_heating_evo_halo_mass_transition_mass_"+model+'_snap_'+char_snap+".txt"
    fa = open(fa, "r")                            
    fields = list(fa)        
    median = float(fields[0].strip())
    pc16 = float(fields[1].strip())
    pc84 = float(fields[2].strip())
    #print(median,pc16,pc84)
    
    snap_list = np.array([61, 51, 42, 34, 29])
    snap_list = np.array([57, 47, 38, 30, 25])
    #snap_list = np.array([42, 29])
    for idx, element in enumerate(snap_list):  
        char_snap ="%d" % element
            
        sel = ((transition_snap > element-2) & ((transition_snap < element+2)))
        print(element, len(transition_mass[sel]))    
        print(np.mean(transition_mass[sel]))      
        median = np.median(transition_mass[sel])
        y_sorted = np.sort(transition_mass[sel])
        pc16 = y_sorted[int(16*len(transition_mass[sel])/100)]      
        pc84 = y_sorted[int(84*len(transition_mass[sel])/100)]  
        fa = open(Datadir+"cooling_heating_evo_halo_mass_transition_mass_at_snap_"+model+'_snap_'+char_snap+".txt", "w")           
        fa.write("%0.5f\n" % median)
        fa.write("%0.5f\n" % pc16)
        fa.write("%0.5f\n" % pc84)
        fa.close() 
            
        fa = Datadir+"cooling_heating_evo_halo_mass_transition_mass_at_snap_"+model+'_snap_'+char_snap+".txt"
        fa = open(fa, "r")                            
        fields = list(fa)        
        median = float(fields[0].strip())
        pc16 = float(fields[1].strip())
        pc84 = float(fields[2].strip())
        print(median,pc16,pc84)
    
    '''(61, 62)  10.1122 (10.05583, 9.70537, 10.50799)
       (51, 129) 10.2605 (10.31141, 9.89888, 10.64233)
       (42, 162) 10.3772 (10.42888, 10.11866, 10.70417)
       (34, 204) 10.5206 (10.52698, 10.29112, 10.81783)
       (29, 27)  10.5599 (10.46732, 10.42009, 10.79474)

    '''
    


    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HWL18_cooling_heating_evo.pdf')
    plt.close()
    
    return    
#end  cooling_heating_evo  
    
       
    
def cooling_heating_evo_old():
   
    #constants
    UnitLength_in_cm=3.08568e+24
    UnitVelocity_in_cm_per_s=100000.
    UnitTime_in_s=UnitLength_in_cm / UnitVelocity_in_cm_per_s
    SOLAR_MASS=1.989e33
    UnitMass_in_g=1.989e+43
    SEC_PER_YEAR=3.155e7
    UnitEnergy_in_cgs = UnitMass_in_g * UnitLength_in_cm**2. / UnitTime_in_s**2.

    
    #sel=log_SSFR>SSFR_cut[ii]
    #subplot.scatter(x_axis[sel], y_axis[sel], s=5, color='blue')           
        
    xlim=[9.0,12.0]
    ylim=[-3.0, 3.0]
    #ylim=[-4.0, 0.0]
    bin=0.25
    
    model = 'Hen15_test'
    
    compute_and_write = 1
    
    #snap_list = np.array([62, 42, 34, 29])
    snap_list = np.array([58, 38, 30, 25])
    snap_list = np.array([58])
    SSFR_cut=[-11.,-11.0, -10.5,-10.]
    
    
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1])) 
    subplot=plt.subplot()    
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1.0))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.5)) 
    subplot.yaxis.set_major_locator(MultipleLocator(2.0))    
    subplot.yaxis.set_minor_locator(MultipleLocator(0.5))     
    xlab='$\log_{10}(M_*)$'; ylab='$\log_{10}(Heating/Cooling)$' 
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)
    
    for idx, element in enumerate(snap_list):
        
        char_snap ="%d" % element
                
        if(compute_and_write == 1): 
            sel = ( (G_MR['SnapNum']==element) & (G_MR['Type']==0) & 
                    (np.log10(G_MR['StellarMass']*1.e10/Hubble_h) > 10.5) & 
                    (np.log10(G_MR['StellarMass']*1.e10/Hubble_h) < 11.2) & 
                    (np.log10(G_MR['Sfr']/(G_MR['StellarMass']*1.e10/Hubble_h)) < SSFR_cut[idx]) )
            #sel = ( (G_MR['SnapNum']==element) & (G_MR['Type']==0) & 
            #        (np.log10(G_MR['StellarMass']*1.e10/Hubble_h) > 10.5) & 
            #        (np.log10(G_MR['StellarMass']*1.e10/Hubble_h) < 11.2))
            
            #log_SSFR=
            total_ngals=len(G_MR[sel])
            #Nsnaps=64
            Nsnaps=58
            Ngals=100
            #Ngals=2000
            if(Ngals>total_ngals):
                Ngals=total_ngals
            print('Snap %d, Ngals=%d' % (element, Ngals))       
          
            G0_MR=np.random.choice(G_MR[sel], size=Ngals, replace=False) 
              
            Redshift=np.zeros(Nsnaps*Ngals,dtype=np.float32)    
            Mvir=np.zeros(Nsnaps*Ngals,dtype=np.float32)           
            HotGas=np.zeros(Nsnaps*Ngals,dtype=np.float32)          
            StellarMass=np.zeros(Nsnaps*Ngals,dtype=np.float32)
            BlackHoleMass=np.zeros(Nsnaps*Ngals,dtype=np.float32)
            cooling_heating = np.zeros(Nsnaps*Ngals,dtype=np.float32)
          
          
            for ii in range (0,Ngals): 
                initial_gal=G0_MR['GalID'][ii]
                G0_MR_tree=G_MR[(G_MR['GalID']>=G0_MR['GalID'][ii]) & (G_MR['GalID']<=G0_MR['LastProgGal'][ii]) 
                               & (G_MR['Type']<2)]    
      
                snap=0  
                currentgal=initial_gal
                while True:                     
                    Gal = G0_MR_tree[G0_MR_tree['GalID']==currentgal]                    
                    Redshift[Nsnaps*ii+snap] = Gal['Redshift']        
                    Mvir[Nsnaps*ii+snap] = np.log10(Gal['Mvir']*1.e10/Hubble_h)                     
                    HotGas[Nsnaps*ii+snap]=np.log10(Gal['HotGas']*1.e10/Hubble_h)                   
                    StellarMass[Nsnaps*ii+snap]=np.log10(Gal['StellarMass']*1.e10/Hubble_h)
                    BlackHoleMass[Nsnaps*ii+snap]=np.log10(Gal['BlackHoleMass']*1.e10/Hubble_h)
                   
           
                    log_StellarMass=(np.log10(Gal['StellarMass']*1.e10/Hubble_h))                    
                    Cooling=Gal['CoolingRate_beforeAGN'] * (UnitTime_in_s/SEC_PER_YEAR)*(SOLAR_MASS/UnitMass_in_g)
                    AgnEfficiency=5.3e-3*(UnitTime_in_s*SOLAR_MASS)/(UnitMass_in_g*SEC_PER_YEAR)
                    AGNcoeff = (1.34e5 / Gal['Vvir']) * (1.34e5 / Gal['Vvir'])           
                    AGNrate= AgnEfficiency * Gal['BlackHoleMass']/Hubble_h* (Gal['HotGas']/Hubble_h) * 10.         
                    EDDrate = 1.3e48 * Gal['BlackHoleMass'] / (UnitEnergy_in_cgs / UnitTime_in_s) / 9e10;           
                    if(AGNrate>EDDrate):
                        AGNrate=EDDrate              
                    AGNheating = AGNcoeff * AGNrate  
                    if(AGNheating/Cooling ==0.):
                        cooling_heating[Nsnaps*ii+snap]=-99.
                    else:   
                        cooling_heating[Nsnaps*ii+snap]=np.log10(AGNheating/Cooling)
                    
                        
                
        
                    prog=Gal['FirstProgGal']
                    currentgal=prog    
                    snap+=1
                    if prog==-1:                       
                        break
                
                if(Ngals>10):
                    if(ii%(Ngals/10)==0):
                        print(ii, Ngals)           
        
            #write to file
            bin = 0.2
            sel = cooling_heating>-99.
            (x_binned, median, mean, pc16, pc84, rms)= median_and_percentiles(bin, xlim[0], max(StellarMass)+bin/2., 
                                                                              StellarMass[sel], cooling_heating[sel]) 
           
            #(x_binned, median, mean, pc16, pc84, rms)= median_and_percentiles_fixed_xx(StellarMass, cooling_heating) 
            fa = open(Datadir+"cooling_heating_evo_"+model+'_snap_'+char_snap+".txt", "w")
            fa.write("%d\n" % len(x_binned))
            for kk in range (0,len(x_binned)):               
                fa.write("%0.2f " % x_binned[kk] + "%0.2f " % median[kk]  + "%0.2f " % pc16[kk]  + "%0.2f\n" % pc84[kk])   
            fa.close() 
            
        
        
           
        fa = Datadir+"cooling_heating_evo_"+model+'_snap_'+char_snap+".txt"
        (x_binned,median, pc16, pc84)=read_data_with_err(fa)  
        subplot.fill_between(x_binned, pc84, pc16, interpolate=True,alpha=0.8, facecolor='grey', edgecolor='grey') 
        subplot.plot(x_binned, median, linestyle='-', linewidth=2, color='black', alpha=0.8)
        #subplot.plot(ran_gal_z, ran_gal_Mvir, linestyle='-', linewidth=2, color='black')
        
        
    
        if(idx==0):
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.15, y_percentage=0.93, color='black', xlog=0, ylog=0,
                        label='$M_{200\mathrm{c}}$', fontsize=10)    
            plot_label (subplot, 'line', xlim, ylim, x_percentage=0.05, y_percentage=0.945, 
                        color='grey', x2_percentage=0.12, xlog=0, ylog=0, linestyle='-', linewidth=2) 
           
        label='$M_{\mathrm{200c,0}}\sim$%0.1f' % element
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.62, y_percentage=0.9, color='black', 
                    xlog=0, ylog=0, label=label, fontsize=10)    
        
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HWL18_cooling_heating_evo.pdf')
    plt.close()
    
    return 
#end  cooling_heating_evo  
    
       
    
    
        
def gradients_evo():
    
    #Get Redshift List
    file=Datadir+'Database_snapshots_table.fits'
    fits_table=fits.open(file)
    snap_table = fits_table[1]       
    FullRedshiftList=snap_table.data['z'][::-1]
    FullRedshiftList=FullRedshiftList[FullRedshiftList>0.]
       
    #Define Arrays to contain Evo of galaxy properties    
    DiskMass_Array = np.zeros([len(FullRedshiftList), RNUM],dtype=np.float32)
    Mass_Array = np.zeros([len(FullRedshiftList), RNUM],dtype=np.float32)
    ColdGas_Array  = np.zeros([len(FullRedshiftList), RNUM],dtype=np.float32)
    H2_Array  = np.zeros([len(FullRedshiftList), RNUM],dtype=np.float32)
    SFR_Array      = np.zeros([len(FullRedshiftList), RNUM],dtype=np.float32)
        
    #First time use these lines to select HaloID    
    #select subsample of galaxies
    '''G0_MR=G_MR[(G_MR['SnapNum']==58) & (G_MR['Type']==0) & 
               (np.log10(G_MR['Mvir']*1.e10/Hubble_h)>13)]  
    G0_MR=G0_MR[np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h))<-10.5]
    Ngals = len(G0_MR)
    print('Ngals={:d}'.format(Ngals))
    
    
    #select one galaxy randomly
    random.seed(a=1)    
    random_list = random.sample(range(0, Ngals), 1)     
    initial_gal=G0_MR['GalID'][random_list]  
    print(G0_MR['HaloID'][random_list], G0_MR['Mvir'][random_list]*1.e10)  
      
    G0_MR_tree=G_MR[(G_MR['GalID']>=G0_MR['GalID'][random_list]) & (G_MR['GalID']<=G0_MR['LastProgGal'][random_list])]'''  
        
    sel = (G_MR['HaloID']==5000069000005) & (G_MR['Type']==0)    
    initial_gal = G_MR['GalID'][sel]
    G0_MR_tree=G_MR[(G_MR['GalID']>=G_MR['GalID'][sel]) & (G_MR['GalID']<=G_MR['LastProgGal'][sel])]
    
    #print(np.log10(G_MR['StellarMass'][sel]*1.e10/Hubble_h))
    
    #Walk The Tree    
    i_snap=0  
    currentgal=initial_gal
    while True:
        sel=G0_MR_tree['GalID']==currentgal
        DiskMass_Array[i_snap,:] = G0_MR_tree['DiskMassRings'][sel,:][0]*1.e10/Hubble_h 
        Mass_Array[i_snap,:] = (G0_MR_tree['DiskMassRings'][sel,:][0]*1.e10/Hubble_h +
                                G0_MR_tree['BulgeMassRings'][sel,:][0]*1.e10/Hubble_h)
        ColdGas_Array[i_snap,:] = G0_MR_tree['ColdGasRings'][sel,:][0]*1.e10/Hubble_h  
        H2_Array[i_snap,:] = G0_MR_tree['ColdGasRings'][sel,:][0]*G0_MR_tree['H2fractionRings'][sel,:][0]*1.e10/Hubble_h  
        SFR_Array[i_snap,:] = G0_MR_tree['SfrRings'][sel,:][0]

        prog=G0_MR_tree['FirstProgGal'][G0_MR_tree['GalID']==currentgal]
        currentgal=prog    
        i_snap+=1
        if prog==-1:
            break
       
          
            
            
    #Area of different Rings        
    Area = np.zeros(RNUM,dtype=np.float32)
    Area[0] = (3.14*RingRadius[0]**2) 
    for jj in range(1, RNUM):
        Area[jj] = (3.14*(RingRadius[jj]**2-RingRadius[jj-1]**2))           
     
    
    #DO PLOTS
    fig = plt.figure(figsize=(two_one_size_small[0],two_one_size_small[1]))
    grid = gridspec.GridSpec(2, 1)
    grid.update(wspace=0.0, hspace=0.0)
    #grid.update(wspace=0.4, hspace=0.3)    
    #plt.subplots_adjust(left=0.1, right=0.97, top=0.97, bottom=0.09) 
    
    
    Nprops = 2
    for i_prop in range(0, Nprops):
        
        subplot=plt.subplot(grid[i_prop])    
              
        #if(i_prop==0):
        #    ylim=[1.5, 4.]
        #    ylab='$\log_{10}(\Sigma_{\mathrm{H_2}}[M_{\odot} \mathrm{pc^{-2}}])$' 
        #    subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
        #    subplot.yaxis.set_minor_locator(MultipleLocator(0.1))            
            
        if(i_prop==0):
            ylim=[-1.5, 1.]
            ylab='$\log_{10}(\Sigma_{\mathrm{SFR}}[M_{\odot} \mathrm{yr^{-1}} \mathrm{kpc^{-2}}])$' 
            subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
            subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
                
        if(i_prop==1):
            ylim=[7., 11.]
            ylab='$\log_{10}(\Sigma_*[M_{\odot} \mathrm{kpc^{-2}}])$' 
            subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
            subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
        
        #xlim=[0.0,30.] 
        xlim=[-0.5, 1.75]
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)   
        xlab='$\log_{10}(r[Kpc])$'  
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14) 
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25)) 
        
        if(i_prop==0):
            plt.tick_params(axis='x', which='both', bottom=True, labelbottom=False)
        
        #set color map
        n = len(FullRedshiftList)              
        subplot.set_prop_cycle('color',[plt.cm.Spectral(i) for i in np.linspace(0.95, 0.05, n-15)])    
        
        '''aa = FullRedshiftList[:n-15]
        bb = np.array([i for i in np.linspace(0, 1, n-15)])           
        labels = [aa[(bb<0.01)][0], aa[(bb>0.19) & (bb<0.21)][0], aa[(bb>0.39) & (bb<0.41)][0], 
                      aa[(bb>0.59) & (bb<0.61)][0], aa[(bb>0.79) & (bb<0.81)][0], aa[(bb>0.99)][0]]
       
        ticklabels = ["{:0.1f}".format(x) for x in labels] 
        
        if(i_prop==0):           
            plt.rcParams.update({'ytick.labelsize': 9, 'axes.linewidth': 0.5, 
                                 'ytick.major.size': 2, 'ytick.major.width': 0.5})             
            col_map = plt.get_cmap('Spectral')           
            c_map_ax = fig.add_axes([0.8, 0.65, 0.05, 0.25]) #left, bottom, width, height]
           
            # and create another colorbar with:
            colorbar.ColorbarBase(c_map_ax, cmap=col_map, orientation = 'vertical', ticks=[0,0.2,0.4,0.6,0.8,1])            
            c_map_ax.set_yticklabels(ticklabels)
            c_map_ax.set_ylabel('redshift', fontsize=9)
            plt.rcParams.update({'ytick.labelsize': 14, 'axes.linewidth': 2, 
                                 'ytick.major.size': 6, 'ytick.major.width': 1.5})'''  
            
        for i_snap in range(0,len(FullRedshiftList)):                           
            #if(i_prop==0):
            #    subplot.plot(RingRadius, np.log10(H2_Array[i_snap,:]/Area), linestyle='-', linewidth=1)    
            
            if(i_prop==0):
                subplot.plot(np.log10(RingRadius), np.log10(SFR_Array[i_snap,:]/Area), linestyle='-', linewidth=1)
                
            if(i_prop==1):
                subplot.plot(np.log10(RingRadius), np.log10(Mass_Array[i_snap,:]/Area), linestyle='-', linewidth=1)    
      
    plt.tight_layout()   
    current_function =  inspect.getframeinfo(inspect.currentframe()).function 
    plt.savefig('./fig/HYF19_'+current_function+'.pdf')
    plt.savefig('./fig/plots_'+current_function+'.pdf')   
    plt.close()
  
    return 
#gradients_evo


def gradients_mean_evo():
    
    plot_color=['orange', 'red'] 
    mass_low = 10.0
    mass_high = 11.0
    
    fig = plt.figure(figsize=(two_one_size_small[0],two_one_size_small[1]),dpi=200)
    grid = gridspec.GridSpec(2, 1)
    grid.update(wspace=0.0, hspace=0.0)    
    subplot_1=plt.subplot(grid[0])
    subplot_2=plt.subplot(grid[1])
  

    xlim=[-0.4, 1.5]
    ylim=[-2.5, 1.]    
    #ylim=[-14, -8.5]         
    subplot_1.set_ylim(ylim),subplot_1.set_xlim(xlim)
    
    xlab=''        
    ylab='$\log_{10}(\Sigma_{SFR}[M_{\odot} \mathrm{yr}^{-1}\mathrm{Kpc^{-2}}])$' 
    subplot_1.set_xlabel(xlab,fontsize=14), subplot_1.set_ylabel(ylab,fontsize=14)       
    subplot_1.tick_params(axis='x', which='both', bottom=True, labelbottom=False)

    subplot_1.yaxis.set_major_locator(MultipleLocator(1.0))    
    subplot_1.yaxis.set_minor_locator(MultipleLocator(0.1)) 
    subplot_1.xaxis.set_major_locator(MultipleLocator(0.5))    
    subplot_1.xaxis.set_minor_locator(MultipleLocator(0.1)) 
    
    xlim=[-0.4, 1.5]
    ylim=[7, 10.5]     
    subplot_2.set_ylim(ylim),subplot_2.set_xlim(xlim)

    xlab='$\log_{10}(r/\mathrm{kpc})$'        
    ylab='$\log_{10}(\Sigma_*[M_{\odot} \mathrm{Kpc^{-2}}])$'        
    subplot_2.set_xlabel(xlab,fontsize=14), subplot_2.set_ylabel(ylab,fontsize=14)       

    subplot_2.yaxis.set_major_locator(MultipleLocator(1.0))    
    subplot_2.yaxis.set_minor_locator(MultipleLocator(0.1)) 
    subplot_2.xaxis.set_major_locator(MultipleLocator(0.5))    
    subplot_2.xaxis.set_minor_locator(MultipleLocator(0.1)) 
       
    #G0_MR = G_MR       
    #G0_MR=np.random.choice(G_MR, size=16000000)     
    #G0_MR=np.random.choice(G_MR, size=2000000)    
    #G0_MR=np.random.choice(G_MR, size=500000)    
   
    #select subsample of galaxies
    G0_MR = G_MR
    G0_MR=G0_MR[(G0_MR['Sfr']>0.) & (G0_MR['SnapNum']==57) & (G0_MR['Type']==0) & 
               (np.log10(G0_MR['StellarMass']*1e10/Hubble_h)>mass_low) &
               (np.log10(G0_MR['StellarMass']*1e10/Hubble_h)<mass_high) ]  
    log_SSFR = np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h))
    sel = log_SSFR >np.log10((1+0.0)**2/(1.37e10/2.))-0.3  #bottom of main sequence
    #sel = log_SSFR <np.log10((1+0.0)**2/(1.37e10/2.))-1.0  #bottom of main sequence
    #sel = log_SSFR -10.5  #bottom of main sequence
 
    G0_MR=G0_MR[sel]   
    print('Ngals={:d}'.format(len(G0_MR)))
    Ngals = min(len(G0_MR),1000)    
    G0_MR=np.random.choice(G0_MR, size=Ngals, replace=False)    
   
                          
    
    #Loop through all galaxies and get a list of GalIDs and redshifts for the main progenitors    
    GalIDs = []
    GalSnaps = []
    GalRedshifts = []    
    for i_gal in range (0,Ngals):  
        #print("doing tree", i_gal)
        if(i_gal%10==0):
            print(i_gal)
        initial_gal = G0_MR['GalID'][i_gal] 
        LastProgGal = G0_MR['LastProgGal'][i_gal]
       
        #print(G0_MR['HaloID'][i_gal], np.log10(G0_MR['Mvir'][i_gal]*1.e10/Hubble_h))  
      
        G0_MR_tree=G_MR[(G_MR['GalID']>=initial_gal) & (G_MR['GalID']<=LastProgGal)]
        
        #print(f"Number of progs: {len(G0_MR_tree):d}")
        #Walk The Tree    
        i_snap=0  
        currentgal=initial_gal
        while True:
            sel=G0_MR_tree['GalID']==currentgal
            GalIDs.append(currentgal.item(0))
            GalSnaps.append(G0_MR_tree['SnapNum'][sel].item(0))            
            GalRedshifts.append(G0_MR_tree['Redshift'][sel].item(0))
            #print(np.log10(G0_MR_tree['StellarMass'][sel]*1.e10/Hubble_h))
            prog=G0_MR_tree['FirstProgGal'][G0_MR_tree['GalID']==currentgal]
            currentgal=prog    
            i_snap+=1
            if prog==-1:
                break
    
    sorted_Snaps = sorted(set(GalSnaps))[::3]  
    sorted_Redshifts = sorted(set(GalRedshifts), reverse=True)[::3]  
    #print(sorted(set(GalSnaps)))
    #print(sorted_Snaps)
    #print(sorted_Redshifts)
    #sorted_Snaps = [58,54,47,38,30,25,22,19,17,15,13,12]    
    #for element in sorted_Redshifts:
    #    print(f'{element:0.2f}')
   
    sorted_Snaps = [13, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57]
    sorted_Redshifts = [8.220924377441406, 6.969651699066162, 5.4567365646362305, 4.282843589782715, 
                        3.3652682304382324, 2.642801523208618, 2.0700461864471436, 1.612944483757019, 
                        1.245949387550354, 0.9497320055961609, 0.7090891003608704, 0.5132874250411987, 
                        0.35295113921165466, 0.2214486449956894, 0.11378046870231628, 0.025611571967601776]

    print('Ids selected')

    
        
    #DO PLOTS
    Area = np.zeros(RNUM,dtype=np.float32)
    Area[0] = (3.14*RingRadius[0]**2) 
    for jj in range(1, RNUM):
        Area[jj] = (3.14*(RingRadius[jj]**2-RingRadius[jj-1]**2))           
                
    n = len(sorted_Snaps)              
    plt_color = [plt.cm.Spectral(i) for i in np.linspace(0.05, 0.95, n)]   
  
    aa = np.array(sorted_Redshifts)  
    bb = np.array([i for i in np.linspace(0, 1, n)])   
    labels = [aa[(bb<0.01)][0], aa[(bb>0.19) & (bb<0.22)][0], aa[(bb>0.39) & (bb<0.43)][0], 
              aa[(bb>0.56) & (bb<0.61)][0], aa[(bb>0.78) & (bb<0.81)][0], aa[(bb>0.99)][0]] 
    labels = labels[::-1]
    ticklabels = ["{:0.1f}".format(x) for x in labels] 
        
    plt.rcParams.update({'ytick.labelsize': 9, 'axes.linewidth': 0.5, 
                         'ytick.major.size': 2, 'ytick.major.width': 0.5})             
    col_map = plt.get_cmap('Spectral_r')           
    c_map_ax = fig.add_axes([0.75, 0.65, 0.05, 0.25]) #left, bottom, width, height]

    # and create another colorbar with:
    colorbar.ColorbarBase(c_map_ax, cmap=col_map, orientation = 'vertical', ticks=[0,0.2,0.4,0.6,0.8,1])            
    c_map_ax.set_yticklabels(ticklabels)
    c_map_ax.set_ylabel('redshift', fontsize=9)
    plt.rcParams.update({'ytick.labelsize': 14, 'axes.linewidth': 2, 
                         'ytick.major.size': 6, 'ytick.major.width': 1.5})
            
    
    GalIDs = np.array(GalIDs)
    GalSnaps = np.array(GalSnaps)
    #Loop through all the snaps and plot one line at each
    i_snap=0
    for snap in sorted_Snaps:          
        sel = GalSnaps == snap         
        G = G_MR[np.isin(G_MR['GalID'], GalIDs[sel])]
        #print(snap, len(G))
        #print(np.log10(G['StellarMass']*1.e10/Hubble_h))
        
        for i_prop in range(0,2):
        
            NGals=len(G)           
            x_variable=np.zeros(int(RNUM*NGals),dtype=np.float32)             
            y_variable=np.zeros(int(RNUM*NGals),dtype=np.float32)   
            new_x_var = np.arange(xlim[0],xlim[1]+1.0,0.2)

            #loop on radial bins
            for jj in range(0,RNUM):  
                x_variable[NGals*jj:NGals*(jj+1)]=np.log10(RingRadius[jj])  
                
                #1e6 -> from kpc^2 to pc^2
                #sigma_SFR
                if(i_prop==0):
                    y_variable[NGals*jj:NGals*(jj+1)]=G['SfrRings'][:,jj]/Area[jj]                   
                #sigma_*
                else:                   
                    y_variable[NGals*jj:NGals*(jj+1)]=((G['DiskMassRings'][:,jj]*1.e10/Hubble_h +
                                                        G['BulgeMassRings'][:,jj]*1.e10/Hubble_h)/Area[jj])    
            #endfor RNUM      
            N_random_gals = 50000
            if(N_random_gals>NGals):
                N_random_gals=NGals

            random.seed(a=2)    
            random_list = random.sample(range(0, NGals), N_random_gals)
            #N_random_gals = 2
            #random_list = random_list[0:1]
            interpol_x_variable=np.zeros(int(len(new_x_var)*N_random_gals),dtype=np.float32)      
            interpol_y_variable=np.zeros(int(len(new_x_var)*N_random_gals),dtype=np.float32)   
            to_do = np.zeros(int(len(new_x_var)*N_random_gals),dtype=np.float32)   

            i_index = 0
            for ii in random_list:              
                slice_ii = [x*NGals+ii for x in range(0,12)]                
                xx = x_variable[slice_ii]
                yy = np.log10(y_variable[slice_ii])  

                #ignore galaxies without halflightradius or with nan on the y_variable
                sel = (~np.isnan(xx)) & (~np.isinf(xx)) & (~np.isnan(yy)) & (~np.isinf(yy))                 
                if(len(xx[sel])>3):  
                    to_do[i_index*len(new_x_var):(i_index+1)*len(new_x_var)]=1
                    f = interpolate.UnivariateSpline(xx[sel], yy[sel], s=0)                      
                    interpol_y_variable[i_index*len(new_x_var):(i_index+1)*len(new_x_var)] = f(new_x_var)
                    interpol_x_variable[i_index*len(new_x_var):(i_index+1)*len(new_x_var)] = new_x_var           
                i_index += 1

            sel = to_do ==1    
            (x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles_fixed_xx(interpol_x_variable[sel], 
                                                                                          interpol_y_variable[sel], non_zero=0)  
            if(i_prop==0):
                subplot_1.plot(x_binned, median, color=plt_color[i_snap], linewidth=1.5, linestyle='-')
                file = Datadir+file_to_write
                #file += '/GradientsMeanEvo/NoFed_NoInflowModel_SF_M10.0_11.0_SigmaSfr'+str(f'_Snap{snap}')+'.csv'  
                file += '/GradientsMeanEvo/FullModel_SF_M10.0_11.0_SigmaSfr'+str(f'_Snap{snap}')+'.csv'
                #file += '/GradientsMeanEvo/FullModel_Passive_M11.0_11.5_SigmaSfr'+str(f'_Snap{snap}')+'.csv'
            else:
                subplot_2.plot(x_binned, median, color=plt_color[i_snap], linewidth=1.5, linestyle='-')
                file = Datadir+file_to_write
                #file += '/GradientsMeanEvo/NoFed_NoInflowModel_SF_M10.0_11.0_SigmaStar'+str(f'_Snap{snap}')+'.csv' 
                file += '/GradientsMeanEvo/FullModel_SF_M10.0_11.0_SigmaStar'+str(f'_Snap{snap}')+'.csv' 
                #file += '/GradientsMeanEvo/FullModel_Passive_M11.0_11.5_SigmaStar'+str(f'_Snap{snap}')+'.csv'  
                
            #WRITE OUTPUT      
            if(write_to_file==1):
                df = pd.DataFrame({'log10_r':x_binned, 'Sigma':median})           
                df.to_csv(file,index=False)
            
        i_snap+=1
      
    
    #plot_label (subplot_2, 'label', xlim, ylim, x_percentage=0.35, y_percentage=0.87, color='black', 
    #            xlog=0, ylog=0, label='$11.0<\log_{10}(M_*/M_{\odot})<12.0$', fontsize=12)
    #plot_label (subplot_2, 'label', xlim, ylim, x_percentage=0.48, y_percentage=0.77, 
    #            color='black', xlog=0, ylog=0, label='Passive at $z=0$', fontsize=12) 
    
    plot_label (subplot_2, 'label', xlim, ylim, x_percentage=0.35, y_percentage=0.87, color='black', 
                xlog=0, ylog=0, label='$10.0<\log_{10}(M_*/M_{\odot})<11.0$', fontsize=12)
    plot_label (subplot_2, 'label', xlim, ylim, x_percentage=0.42, y_percentage=0.77, 
                color='black', xlog=0, ylog=0, label='Star-Forming at $z=0$', fontsize=12) 
      
        
      
        
    plt.tight_layout()   
    current_function =  inspect.getframeinfo(inspect.currentframe()).function 
    plt.savefig('./fig/HYF19_'+current_function+'.pdf')
    plt.savefig('./fig/plots_'+current_function+'.pdf')   
    plt.close()
  
    return 
#gradients_mean_evo
    
    

    
def gradients_mean_evo_combined():
    
    plot_color=['orange', 'red'] 
    
    
    fig = plt.figure(figsize=(two_three_size_large[0],two_three_size_large[1]))
    grid = gridspec.GridSpec(2, 3)
    grid.update(wspace=0.0, hspace=0.0)    
    
    Nprop = 2
    NSelect = 3
    
    label_mass_sel = ['$10.0<\log_{10}(M_*/M_{\odot})<11.0$', '$10.0<\log_{10}(M_*/M_{\odot})<11.0$', 
                      '$11.0<\log_{10}(M_*/M_{\odot})<11.5$']
    label_SF_sel = ['Star-Forming at $z=0$','Star-Forming at $z=0$','     Passive at $z=0$',]
    
    label_file = ['NoFed_NoInflowModel_SF_M10.0_11.0_SigmaSfr', 'FullModel_SF_M10.0_11.0_SigmaSfr',
                  'FullModel_Passive_M11.0_11.5_SigmaSfr', 'NoFed_NoInflowModel_SF_M10.0_11.0_SigmaStar',
                  'FullModel_SF_M10.0_11.0_SigmaStar','FullModel_Passive_M11.0_11.5_SigmaStar' ]
    
    label_model = ['No Feedback, No Inflow', 'Full Model', 'Full Model']
                
    i_grid = 0
    for i_prop in range(0, Nprop):
        
        for i_select in range(0, NSelect):
    
            subplot = plt.subplot(grid[i_grid])           
    
            xlim=[-0.4, 1.5]
            xlab='$\log_{10}(r/\mathrm{kpc})$'  
            
            if (i_prop==0):
                ylim=[-2.5, 1.]    
                ylab='$\log_{10}(\Sigma_{\mathrm{SFR}}/(M_{\odot} \mathrm{yr}^{-1}\mathrm{kpc^{-2}}))$' 
                subplot.tick_params(axis='x', which='both', bottom=True, labelbottom=False)
            else:
                ylim=[7, 10.5]     
                ylab='$\log_{10}(\Sigma_*/(M_{\odot} \mathrm{kpc^{-2}}))$'  
                
            if(i_select>0):
                subplot.tick_params(axis='y', which='both', left=True, labelleft=False)
                ylab=''
                
            subplot.set_ylim(ylim),subplot.set_xlim(xlim)
            subplot.set_xlabel(xlab,fontsize=14), subplot.set_ylabel(ylab,fontsize=14)      
            
            
            
            
            
            subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
            subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
            subplot.xaxis.set_major_locator(MultipleLocator(0.5))    
            subplot.xaxis.set_minor_locator(MultipleLocator(0.1)) 
            
            
            
            
            sorted_Snaps = [13, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57]
            sorted_Redshifts = [8.220924377441406, 6.969651699066162, 5.4567365646362305, 4.282843589782715, 
                                3.3652682304382324, 2.642801523208618, 2.0700461864471436, 1.612944483757019, 
                                1.245949387550354, 0.9497320055961609, 0.7090891003608704, 0.5132874250411987, 
                                0.35295113921165466, 0.2214486449956894, 0.11378046870231628, 0.025611571967601776]


            #COLOR BAR    
            n = len(sorted_Snaps)              
            plt_color = [plt.cm.Spectral(i) for i in np.linspace(0.05, 0.95, n)]   

            aa = np.array(sorted_Redshifts)  
            bb = np.array([i for i in np.linspace(0, 1, n)])   
            labels = [aa[(bb<0.01)][0], aa[(bb>0.19) & (bb<0.22)][0], aa[(bb>0.39) & (bb<0.43)][0], 
                      aa[(bb>0.56) & (bb<0.61)][0], aa[(bb>0.78) & (bb<0.81)][0], aa[(bb>0.99)][0]] 
            labels = labels[::-1]
            ticklabels = ["{:0.1f}".format(x) for x in labels] 

            plt.rcParams.update({'ytick.labelsize': 9, 'axes.linewidth': 0.5, 
                                 'ytick.major.size': 2, 'ytick.major.width': 0.5})             
            col_map = plt.get_cmap('Spectral_r')           
            c_map_ax = fig.add_axes([0.88, 0.68, 0.02, 0.25]) #left, bottom, width, height]

            # and create another colorbar with:
            colorbar.ColorbarBase(c_map_ax, cmap=col_map, orientation = 'vertical', ticks=[0,0.2,0.4,0.6,0.8,1])            
            c_map_ax.set_yticklabels(ticklabels)
            c_map_ax.set_ylabel('redshift', fontsize=9)
            plt.rcParams.update({'ytick.labelsize': 14, 'axes.linewidth': 2, 
                                 'ytick.major.size': 6, 'ytick.major.width': 1.5})
          
        
            #PLOT LINES    
            i_snap=0
            for snap in sorted_Snaps:  
                file = Datadir+file_to_write+'/GradientsMeanEvo/'+label_file[i_grid]+str(f'_Snap{snap}')+'.csv' 
                print(file)
                df = pd.read_csv(file)
                subplot.plot(df['log10_r'],df['Sigma'],color=plt_color[i_snap])
                i_snap+=1
               


            #plot_label (subplot_2, 'label', xlim, ylim, x_percentage=0.35, y_percentage=0.87, color='black', 
            #            xlog=0, ylog=0, label='$11.0<\log_{10}(M_*/M_{\odot})<12.0$', fontsize=12)
            #plot_label (subplot_2, 'label', xlim, ylim, x_percentage=0.48, y_percentage=0.77, 
            #            color='black', xlog=0, ylog=0, label='Passive at $z=0$', fontsize=12) 
            
            if (i_prop == 1):
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.13, y_percentage=0.87, color='black', 
                            xlog=0, ylog=0, label=label_mass_sel[i_select], fontsize=10)
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.21, y_percentage=0.77, 
                            color='black', xlog=0, ylog=0, label=label_SF_sel[i_select], fontsize=10) 
            if (i_prop == 0):    
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.07, y_percentage=0.85, 
                            color='black', xlog=0, ylog=0, label=label_model[i_select], fontsize=12) 
            
            
            i_grid+=1
    
     
    

    
    
     
   
    
      
        
      
        
    plt.tight_layout()   
    current_function =  inspect.getframeinfo(inspect.currentframe()).function 
    plt.savefig('./fig/HYF19_'+current_function+'.pdf')
    plt.savefig('./fig/plots_'+current_function+'.pdf')   
    plt.close()
  
    return 
#gradients_mean_evo_combined
    
    
    
def BT_evo():  
    
    plot_color=['orange', 'red']     
    fig, subplot = plt.subplots(figsize=(one_one_size_small[0],one_one_size_small[1]))
       
    xlim=[0., 5.0]
    ylim=[0.0,1.0]
    xlab='$redshift'  

   
    #ylim=[-2.5, 1.]    
    ylab='B/T' 
    #subplot.tick_params(axis='x', which='both', bottom=True, labelbottom=False)

    #subplot.tick_params(axis='y', which='both', left=True, labelleft=False)
      

    subplot.set_ylim(ylim),subplot.set_xlim(xlim)
    subplot.set_xlabel(xlab,fontsize=14), subplot.set_ylabel(ylab,fontsize=14)      





    subplot.yaxis.set_major_locator(MultipleLocator(0.5))    
    subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
    subplot.xaxis.set_major_locator(MultipleLocator(1.0))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.1)) 
            
    mass_low = 10.5
    mass_high = 11.0
    
    G0_MR = G_MR
    G0_MR=G0_MR[(G0_MR['Type']==0) & (G0_MR['SnapNum']==58) & 
               (np.log10(G0_MR['StellarMass']*1e10/Hubble_h)>mass_low) &
               (np.log10(G0_MR['StellarMass']*1e10/Hubble_h)<mass_high) ]  
    log_SSFR = np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h))
    #sel = log_SSFR >np.log10((1+0.0)**2/(1.37e10/2.))-0.3  #bottom of main sequence
    #sel = log_SSFR <np.log10((1+0.0)**2/(1.37e10/2.))-1.0  #bottom of main sequence
   
 
    #G0_MR=G0_MR[sel]   
    print('Ngals={:d}'.format(len(G0_MR)))
    Ngals = min(len(G0_MR),100)    
    G0_MR=np.random.choice(G0_MR, size=Ngals, replace=False)    
    
    sel = G0_MR['BulgeMass']/G0_MR['StellarMass']>0.5
    print(len(G0_MR['BulgeMass'][sel]))
    
    #Loop through all galaxies and get a list of GalIDs and redshifts for the main progenitors    
    
    for i_gal in range (0, Ngals):        
        Snaps = []
        Redshifts = []  
        log_StellarMass = []
        BT = []
        
        if(i_gal%10==0):
            print(i_gal)
        initial_gal = G0_MR['GalID'][i_gal] 
        LastProgGal = G0_MR['LastProgGal'][i_gal]
       
        G0_MR_tree=G_MR[(G_MR['GalID']>=initial_gal) & (G_MR['GalID']<=LastProgGal)]
                
        #Walk The Tree    
        i_snap=0  
        currentgal=initial_gal
        while True:
            sel=G0_MR_tree['GalID']==currentgal           
            Snaps.append(G0_MR_tree['SnapNum'][sel][0])            
            Redshifts.append(G0_MR_tree['Redshift'][sel][0])
            #print(G0_MR_tree['SnapNum'][sel])
            log_StellarMass.append(np.log10(G0_MR_tree['StellarMass'][sel][0]*1.e10/Hubble_h))
            BT.append(G0_MR_tree['BulgeMass'][sel][0]/G0_MR_tree['StellarMass'][sel][0])
        
            prog=G0_MR_tree['FirstProgGal'][G0_MR_tree['GalID']==currentgal]
            currentgal=prog    
            i_snap+=1
            if prog==-1:
                break
     
        subplot.plot(Redshifts,BT,linewidth=1)
        
        
        
        
        
        
        
    plt.tight_layout()   
    current_function =  inspect.getframeinfo(inspect.currentframe()).function 
    plt.savefig('./fig/HYF19_'+current_function+'.pdf')
    plt.savefig('./fig/plots_'+current_function+'.pdf')   
    plt.close()
  
    return         