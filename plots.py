
# coding: utf-8

# # Python Plots for LGalaxies

# ## Import Libraries and Read Catalogs

# <p>Use functions read_snap or read_tree to read catalogs. These are both defined in procedures.py. In case of read_snap, SnapshotList will be returned containing the list of snapshots read (usefull to later select galaxies in a given redshift).<p>

# In[1]:

import numpy as np
get_ipython().magic('matplotlib inline')

import pandas as pd

get_ipython().magic('pylab inline')
import seaborn as sns
sns.set_style('darkgrid')

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.table import Table
from importlib import reload
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm

import procedures
reload (procedures)
from procedures import *
import plots_input
reload (plots_input)
from plots_input import *

from LGalaxies_Henriques2015a_struct import LGalaxiesStruct
from LGalaxies_Henriques2015a_struct import PropertiesToRead


FirstFile = 40
LastFile =  41

Volume_MR = (BoxSize_MR**3.0) * (LastFile - FirstFile + 1) / MaxTreeFiles 
Volume_MRII = (BoxSize_MRII**3.0) * (LastFile - FirstFile + 1) / MaxTreeFiles 

(G_MR, SnapshotList) = read_snap(DirName_MR,FirstFile,LastFile,
                 PropertiesToRead,LGalaxiesStruct,RedshiftsToRead,RedshiftList)
 
#print (np.log10(gal['StellarMass'][1:5]*1.e10))
#help(gal)


# ## Plots

# In[33]:

plt.rcParams.update({'font.size': 18, 'xtick.labelsize': 18, 'ytick.labelsize': 18, 'axes.linewidth': 2, 
                     'xtick.major.size': 6, 'xtick.major.width': 1.5, 
                     'ytick.major.size': 6, 'ytick.major.width': 1.5, 
                     'xtick.minor.size': 3, 'xtick.minor.width': 1.,                   
                     'ytick.minor.size': 3, 'ytick.minor.width': 1.})


# In[ ]:

with PdfPages('./fig/plots.pdf') as pdf:

    
    
    if stellar_mass_function==1:
        
        xmin=7.0
        xmax=12.5
        ymin=-6.5
        ymax=0.5
        bin=0.1


        plt.rcParams.update({'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 
                             'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0})

        fig = plt.figure(figsize=(9,9))
        grid = gridspec.GridSpec(2, 2)
        grid.update(wspace=0.0, hspace=0.0)

        for ii in range(0,4):
    
            if ii == 0 :
                redshift=RedshiftList[ii]
            else :  #avoid using z=0.4
                redshift=RedshiftList[ii+1] 
    
            subplot=plt.subplot(grid[ii])

            subplot.set_ylim([ymin, ymax]), subplot.set_xlim([xmin, xmax])
            if ii==2 or ii == 3: 
                xlab='$log_{10}(M_*[h^{-2}M_{\odot}])$'
            else:
                xlab=''
            if ii==0 or ii == 2:
                ylab='$log_{10}(\phi [h^3 Mpc^{-3} log_{10}(M^{-1})])$'
            else:
                ylab=''      
            subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)
       
            majorFormatter = FormatStrFormatter('%d')
            subplot.xaxis.set_major_locator(MultipleLocator(1))    
            subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
            subplot.yaxis.set_minor_locator(MultipleLocator(0.25))
    
            if ii==1 or ii == 3:
                plt.tick_params(axis='y', which='both', left='on', labelleft='off')
           
    
            #OBSERVATIONS   
            char_redshift="%0.0f" % redshift       
            file = Datadir + '/ObsConstraints/StellarMassFunction_z'+char_redshift+'.00.txt'
            obs = Table.read(file, format='ascii')
    
            obs_xbin=obs['col1']+(obs['col2']-obs['col1'])/2.
            asy_yerror = [np.log10(obs['col3']/(obs['col3']-obs['col4'])), 
                          np.log10((obs['col3']+obs['col4'])/obs['col3'])]
            subplot.errorbar(obs_xbin, np.log10(obs['col3']),yerr=asy_yerror,
                     fmt='o', markersize=5, ecolor='blue', color='blue')
            #sub = plt.subplot(111)
    
    
            #MODEL
            if ii == 0 :
                sel= (G_MR['SnapNum']==SnapshotList[ii])
            else :  #avoid using z=0.4    
                sel=(G_MR['SnapNum']==SnapshotList[ii+1])
            G0_MR=G_MR[sel]   
            G0_MR=G0_MR[G0_MR['StellarMass']>0.]
            StellarMass=(np.log10(G0_MR['StellarMass']*1.e10*Hubble_h) +
                         np.random.randn(len(G0_MR['StellarMass']))*0.08*(1+redshift))

            bin_arr=np.arange(7.0,12.0+bin,bin)
            hist=np.histogram(StellarMass, bins=bin_arr, range=(7.0,12.0))   
            subplot.plot(hist[1][0:len(hist[1][:])-1]+bin/2.,np.log10(hist[0][:]/(Volume_MR*bin)),
                     color='red', linewidth=2)


            #LABELS
            if ii==0:
                subplot.text(7.4, 0.0, 'Observations used in MCMC', fontsize= 13)
                subplot.errorbar(7.3, 0.12, yerr=0.1, fmt='o', markersize=5, color='blue')
        
            if ii==3:
                subplot.text(7.7, -5.0, 'This Work', fontsize= 13)
                subplot.plot([7.3,7.6], [-4.9,-4.9], linestyle='-', linewidth=2, color='red')    
        #endfor


        plt.tight_layout()
        plt.savefig('./fig/plots_smf_evo.pdf')
        pdf.savefig()
        plt.close()
    #endif stellar_mass_function==1:

    
    
    
    if bluck_red_fractions==1:
        
        sel= (G_MR['SnapNum']==SnapshotList[0])       
        G0_MR=G_MR[sel]   
               
        xmin=8.0
        xmax=13.0
        ymin=-12.
        ymax=-8.0
        bin=0.1

        plt.rcParams.update({'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 
                             'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0})
        fig = plt.figure(figsize=(9,9))
        grid = gridspec.GridSpec(1, 2)
        #subplot=plt.subplot(grid[ii]) 
        subplot=plt.subplot()
               
        #FIND PASSIVE CUT
        subplot.set_ylim([ymin, ymax]), subplot.set_xlim([xmin, xmax])         
        xlab='$log_{10}(M_*[h^{-2}M_{\odot}])$'           
        ylab='$log_{10}(\phi [h^3 Mpc^{-3} log_{10}(M^{-1})])$'               
        subplot.set_xlabel(xlab, fontsize=16)
        subplot.set_ylabel(ylab, fontsize=16)     
        
        
        BHMass=np.log10(G0_MR['BlackHoleMass']*1.e10/Hubble_h) 
        StellarMass=np.log10(G0_MR['StellarMass']*1.e10/Hubble_h) 
        SFR=np.log10(G0_MR['Sfr'] ) 
        SSFR=np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h)) 
            
        sel=np.logical_and(SSFR>-100.,StellarMass>-10.)    
        StellarMass=StellarMass[sel]
        SSFR=SSFR[sel]
        #subplot.plot(StellarMass, SSFR, 'o', markersize=2, color='blue')
             
        #subplot.hist2d(StellarMass, SSFR, bins=30, norm=LogNorm())  
        #subplot.hexbin(StellarMass, SSFR, gridsize=200)
        #plt.colorbar()
        
        #plt.tight_layout()
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
        subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)   
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
        
        plot_colors=['blue','green','yellow','red']
        halo_min_mass=11.
        halo_bin=1.0
        for jj in range (0,4):
            print((halo_min_mass+jj*halo_bin),(halo_min_mass+(jj+1)*halo_bin))
            sel=np.logical_and(HaloMass > (halo_min_mass+jj*halo_bin), HaloMass < (halo_min_mass+(jj+1)*halo_bin))
            Nbins=int((xmax-xmin)/bin+1)  
            hist=np.array([],dtype=np.float64)
            x_array=np.array([],dtype=np.float64)
            for ii in range(0,Nbins):          
                passive  = SSFR[np.logical_and(np.logical_and(BHMass[sel] > (xmin+ii*bin), BHMass[sel] < (xmin+(ii+1)*bin)),
                                                SSFR[sel]<-10.5)]
                total  = SSFR[np.logical_and(BHMass[sel] > (xmin+ii*bin), BHMass[sel] < (xmin+(ii+1)*bin))]
                if len(total)>0. :
                    hist=np.append(hist,len(passive)/len(total))
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
        subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)   
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
        
        plot_colors=['blue','green','yellow','red']
        bh_min_mass=5.
        bh_bin=1.0
        
        for jj in range (0,4):
            
            print((bh_min_mass+jj*bh_bin),(bh_min_mass+(jj+1)*bh_bin))
            sel=np.logical_and(BHMass > (bh_min_mass+jj*bh_bin), BHMass < (bh_min_mass+(jj+1)*bh_bin))
            
            Nbins=int((xmax-xmin)/bin+1)  
            hist=np.array([],dtype=np.float64)
            x_array=np.array([],dtype=np.float64)
            for ii in range(0,Nbins):          
                passive  = SSFR[np.logical_and(np.logical_and(HaloMass[sel] > (xmin+ii*bin), HaloMass[sel] < (xmin+(ii+1)*bin)),
                                                SSFR[sel]<-10.5)]
                total  = SSFR[np.logical_and(HaloMass[sel] > (xmin+ii*bin), HaloMass[sel] < (xmin+(ii+1)*bin))]
                if len(total)>0. :
                    hist=np.append(hist,len(passive)/len(total))
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
        
        
    #endif bluck_red_fractions
    
    
if test_plots==1:
       
    fig = plt.figure(figsize=(10,10))
        
    sel= np.logical_and(np.logical_and(np.logical_and(G_MR['SnapNum']==SnapshotList[0], G_MR['StellarMass']>0.),
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
    sel= numpy.random.uniform(0.0,1.0,Ngals) < NN/Ngals 
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
#end test_plots


# In[ ]:




# In[10]:




# In[ ]:



