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
from astropy.io import fits
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
import sys
from scipy.ndimage import zoom
from importlib import reload
import inspect   
from scipy import interpolate
import random

import procedures
reload (procedures)
from procedures import *
import plots_input
reload (plots_input)
from plots_input import *

def milkyway_sfr_and_gas_profiles(ThisRedshiftList):
    
    labels_to_write=['sigma_star', 'sigma_HI', 'sigma_H2', 'sigma_Cold']
    
    plot_color=['red','purple']        
   
    fig = plt.figure(figsize=(two_two_size_small[0],two_two_size_small[1]))
    grid = gridspec.GridSpec(2, 2)
    grid.update(wspace=0.0, hspace=0.0)
     
    List_Names_Massive=['NGC0628','NGC3198','NGC3184','NGC4736','NGC3351','NGC6946',
                        'NGC3627','NGC5194','NGC3521','NGC2841','NGC5055','NGC7331'] 
    List_Vflat_Massive=np.array([217.,150.,210.,156.,196.,186.,
                                 192.,219.,227.,302.,192.,244.])  
    
    List_Names_Dwarfs=['DDO145','HoI','HoII','IC2574','NGC4214','NGC2976',
                       'NGC4449','NGC3077','NGC7793','NGC2403','NGC0925'] 
    List_Vflat_Dwarfs=np.array([50.,53.,36.,134.,37.,92.,
                               -99.,-99.,115.,134.,136.,])
    Ngals=12
    
    #READ Observations    
    
    #BLUEDISKS
    '''df_bd = pd.read_csv(Datadir+'bluedisks_props.csv', index_col=0)    
    #display(df_bd[:3])
    #df_bd.to_csv(Datadir+'bluedisks_props.csv',index=False)
    sel = (df_bd['log_Mstar']>10.8) & (df_bd['log_M_HI']>9.5)
    bd_indexes = df_bd[sel].index.values.tolist()   
    print(bd_indexes)
    df_bd_profs = pd.read_csv(Datadir+'bluedisks_profiles.csv')
    #display(df_bd_profs[:3])'''
    
    #LEROY
    index=0
    fa = open(Datadir+"/leroy2008_gradients.txt", "r") 
    Leroy_Names=[]
    for line in fa:           
        if(index==0):                
            fields = line.strip().split() 
            n_header=int(fields[0])                    
        if(index==n_header+1):
            fields = line.strip().split()            
            Aux_Leroy = np.zeros(int(fields[0]), dtype=[('radius',np.float32),('normradius',np.float32),
                                                    ('SigmaHI',np.float32),('SigmaHI_err',np.float32),
                                                    ('SigmaH2',np.float32),('SigmaH2_err',np.float32),
                                                    ('SigmaStar',np.float32),('SigmaStar_err',np.float32),
                                                    ('SigmaSFR',np.float32),('SigmaSFR_err',np.float32),
                                                    ('SigmaSFR_FUV',np.float32),('SigmaSFR_24',np.float32)])          
        if(index>n_header+1):                  
            fields = line.strip().split()         
            Leroy_Names.append(fields[0])
            #Leroy['name'][index-n_header-2]=str(fields[0])         
            Aux_Leroy['radius'][index-n_header-2]=float(fields[1])
            Aux_Leroy['normradius'][index-n_header-2]=float(fields[2])   
            Aux_Leroy['SigmaHI'][index-n_header-2]=float(fields[3])
            Aux_Leroy['SigmaHI_err'][index-n_header-2]=float(fields[4])                 
            Aux_Leroy['SigmaH2'][index-n_header-2]=float(fields[5])
            Aux_Leroy['SigmaH2_err'][index-n_header-2]=float(fields[6])
            Aux_Leroy['SigmaStar'][index-n_header-2]=float(fields[7])
            Aux_Leroy['SigmaStar_err'][index-n_header-2]=float(fields[8])
            Aux_Leroy['SigmaSFR'][index-n_header-2]=float(fields[9])
            Aux_Leroy['SigmaSFR_err'][index-n_header-2]=float(fields[10]) 
            Aux_Leroy['SigmaSFR_FUV'][index-n_header-2]=float(fields[11])                 
            Aux_Leroy['SigmaSFR_24'][index-n_header-2]=float(fields[12])                 
        index+=1  
    #endfor
    #print(Leroy['name'])    
    fa.close()    
    
    Obs_Gal=[('radius',np.float32,100),('normradius',np.float32,100),('SigmaHI',np.float32,100),
             ('SigmaHI_err',np.float32,100),('SigmaH2',np.float32,100),('SigmaH2_err',np.float32,100),
             ('SigmaStar',np.float32,100),('SigmaStar_err',np.float32,100),('SigmaSFR',np.float32,100),
             ('SigmaSFR_err',np.float32,100),('SigmaSFR_FUV',np.float32,100),('SigmaSFR_24',np.float32,100)] 
    
    Leroy = np.zeros(Ngals,dtype=Obs_Gal)
    #next lines are to avoid have a lot of zeros at zero radius
    Leroy['radius'][:,:]=-99.
    Leroy['normradius'][:,:]=-99.
    
    for ii in range(0,len(List_Names_Massive)): 
        index=0
        for jj in range(0,len(Leroy_Names)):                
            if(Leroy_Names[jj]==List_Names_Massive[ii]): 
                Leroy['radius'][ii,index]        = Aux_Leroy['radius'][jj]
                Leroy['normradius'][ii,index]    = Aux_Leroy['normradius'][jj]
                Leroy['SigmaHI'][ii,index]       = Aux_Leroy['SigmaHI'][jj]
                Leroy['SigmaHI_err'][ii,index]   = Aux_Leroy['SigmaHI_err'][jj]                 
                Leroy['SigmaH2'][ii,index]       = Aux_Leroy['SigmaH2'][jj]
                Leroy['SigmaH2_err'][ii,index]   = Aux_Leroy['SigmaH2_err'][jj]
                Leroy['SigmaStar'][ii,index]     = Aux_Leroy['SigmaStar'][jj]
                Leroy['SigmaStar_err'][ii,index] = Aux_Leroy['SigmaStar_err'][jj]
                Leroy['SigmaSFR'][ii,index]      = Aux_Leroy['SigmaSFR'][jj]
                Leroy['SigmaSFR_err'][ii,index]  = Aux_Leroy['SigmaSFR_err'][jj]
                Leroy['SigmaSFR_FUV'][ii,index]  = Aux_Leroy['SigmaSFR_FUV'][jj]               
                Leroy['SigmaSFR_24'][ii,index]   = Aux_Leroy['SigmaSFR_24'][jj]  
                
                index+=1
                #radius.append(Leroy['radius'][jj])
                #sigma.append(Leroy['SigmaStar'][jj])                   
                    
                    
    for ii in range(0,len(ThisRedshiftList)):        
        
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]       
        SSFR = np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h))   
        '''G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (G0_MR['DiskMass']>0.) &       
                    (G0_MR['Vvir']>150.) & (G0_MR['Vvir']<235.) & (G0_MR['Type']==0) & 
                    (SSFR>np.log10(2.*(1+ThisRedshiftList[ii])**2./(1.37e10))-1.5) &
                    (G0_MR['BulgeMass']/G0_MR['StellarMass']<0.3)]''' 
        G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (G0_MR['DiskMass']>0.) & (G0_MR['Type']==0) &                    
                    (SSFR>np.log10(2.*(1+ThisRedshiftList[ii])**2./(1.37e10))-1.0) &
                    #(G0_MR['Sfr']>0.01) &
                    (G0_MR['Vvir']>200.) & (G0_MR['Vvir']<235.) &
                    (G0_MR['BulgeMass']/G0_MR['StellarMass']<0.15) &
                    (np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)>10.3) &
                    (np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)<10.7) &
                    (G0_MR['H2fraction']>0.)]
        #180, 235
        #G0_MR=G0_MR[(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)>10.) & (G0_MR['DiskMass']>0.) &       
        #            (G0_MR['Vmax']>125.) & (G0_MR['MagDust'][:,1]<-20.) &
        #            ((G0_MR['BulgeMass']/G0_MR['StellarMass'])<0.15)]
     
        print("Ngals in Milkyway selection:",len(G0_MR))
        #print(G0_MR['Vvir'],G0_MR['Mvir']*1e10/Hubble_h)
        xlim=[0.0,19.0]
        ylim=[0.,3.8]         
        bin_obs=2.
        
        #0->Sigmastar, 1->SigmaHI, 2->SigmaH2, 3->SigmaGas
        for i_property in range(0,4):
                
            subplot=plt.subplot(grid[i_property])    
            subplot.set_ylim(ylim), subplot.set_xlim(xlim)   
         
            majorFormatter = FormatStrFormatter('%d')
            subplot.xaxis.set_major_locator(MultipleLocator(5))    
            subplot.xaxis.set_minor_locator(MultipleLocator(1.))
            subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
            subplot.yaxis.set_minor_locator(MultipleLocator(0.1))
            
            if((i_property==0) or (i_property==2)):
                ylab='$\log_{10}(\Sigma/(M_{\odot} \mathrm{pc^{-2}}))$'     
                subplot.set_ylabel(ylab, fontsize=14)                 
            if((i_property==2) or (i_property==3)):    
                xlab='$r/\mathrm{kpc}$'
                subplot.set_xlabel(xlab, fontsize=14)  
                
            if (i_property==0) or (i_property==1):    
                plt.tick_params(axis='x', which='both', bottom='on', labelbottom='off')
            if (i_property==1) or (i_property==3):
                plt.tick_params(axis='y', which='both', left='on', labelleft='off')
           
        
            #PLOT OBSERVATIONS           
            #MEDIAN 
            if(i_property==0):              
                Sigma=Leroy['SigmaStar']
                Sigma_err=Leroy['SigmaStar_err']
            if(i_property==1):
                Sigma=Leroy['SigmaHI']
                Sigma_err=Leroy['SigmaHI_err']
            if(i_property==2):
                Sigma=Leroy['SigmaH2']
                Sigma_err=Leroy['SigmaH2_err']
            if(i_property==3):
                Sigma=Leroy['SigmaHI']+Leroy['SigmaH2']
                Sigma_err=Leroy['SigmaHI_err']+Leroy['SigmaH2_err']
            sel=Sigma>0.
            '''(x_binned,median,mean,pc16,pc84,rms)=median_and_percentiles(bin_obs,xlim[0]+0.5,xlim[1],
                                                                        Leroy['radius'][sel],Sigma[sel])  
            subplot.plot(x_binned, np.log10(median),color='black', linewidth=2) 
            #subplot.plot(x_binned, pc16,color='black', linewidth=2) 
            #subplot.plot(x_binned, pc84,color='black', linewidth=2)     
            subplot.fill_between(x_binned,np.log10(pc16),np.log10(pc84),facecolor='lightgrey', 
                                         interpolate=True,alpha=0.8,edgecolor='black')'''
        
        
            #DATA POINTS
            #error=[np.log10((Leroy['SigmaStar']+Leroy['SigmaStar_err'])/Leroy['SigmaStar']),
            #       np.log10(Leroy['SigmaStar']/(Leroy['SigmaStar']-Leroy['SigmaStar_err']))]           
            #subplot.errorbar(Leroy['radius'], np.log10(Leroy['SigmaStar']), error, 
            #                 fmt='o', markersize=3, ecolor='blue', color='blue')  
        
            #ONE LINE PER GALAXY       
            for igal in range(0,len(List_Names_Massive)):
                if((List_Vflat_Massive[igal]>200.) & (List_Vflat_Massive[igal]<235.)):                   
                    subplot.plot(Leroy['radius'][igal,:], np.log10(Sigma[igal,:]),
                                 color='blue', linewidth=1, linestyle='-') 
                    y_err=[np.log10(Sigma[igal,:]/(Sigma[igal,:]-Sigma_err[igal,:])),
                           np.log10((Sigma[igal,:]+Sigma_err[igal,:])/Sigma[igal,:])]
                    subplot.errorbar(Leroy['radius'][igal,:], np.log10(Sigma[igal,:]),yerr=y_err,color='blue',zorder=-1) 
                           
            #BLUE DISKS
            if(i_property==1):
                for jj in range(0,5):
                    df = pd.read_csv(Datadir+'Bluediscs_HI_profiles_'+str(jj+1)+'.csv')      
                    subplot.plot(df['x'], np.log10(df['y']), color='seagreen', linewidth=2, linestyle='-') 
                    
            if(i_property==2):
                for jj in range(0,4):
                    df = pd.read_csv(Datadir+'Bluediscs_H2profiles_'+str(jj+1)+'.csv')      
                    subplot.plot(df['x'], np.log10(df['y']), color='seagreen', linewidth=2, linestyle='-')         
               
            
    
            #MODEL
            Sigma=np.zeros(RNUM,dtype=np.float32)
            Sigma_mean=np.zeros(RNUM,dtype=np.float32)
            Bulge_density=np.zeros(RNUM,dtype=np.float32)
            pc16=np.zeros(RNUM,dtype=np.float32) 
            pc84=np.zeros(RNUM,dtype=np.float32)             
            
            for kk in range(0,RNUM):
                if(kk==0):
                    r_in = 0.
                else:
                    r_in = (RingRadius[kk-1])                
                r_out = RingRadius[kk]
               
                        
                if(i_property==0):               
                    Mass=G0_MR['DiskMassRings'][:,kk]*1e10/Hubble_h + G0_MR['BulgeMassRings'][:,kk]*1e10/Hubble_h    
                if(i_property==1):
                    Mass=G0_MR['ColdGasRings'][:,kk]*1e10/Hubble_h*(1.-G0_MR['H2fractionRings'][:,kk])                   
                if(i_property==2):
                    Mass=G0_MR['ColdGasRings'][:,kk]*1e10/Hubble_h*G0_MR['H2fractionRings'][:,kk]                   
                if(i_property==3):
                    Mass=G0_MR['ColdGasRings'][:,kk]*1e10/Hubble_h
                    
                    
                y_variable=Mass/(3.14*(r_out**2-r_in**2)*1e6)                 
                Sigma[kk]=np.median(y_variable)           
                y_sorted = np.sort(y_variable)
                pc16[kk] = y_sorted[int(16*len(y_variable)/100)]      
                pc84[kk] = y_sorted[int(84*len(y_variable)/100)] 
                                           
                Sigma_mean[kk]=np.mean(y_variable)
            
            if(i_property==2):
                print(Sigma)
            subplot.plot(RingRadius, np.log10(Sigma_mean),color='red', linewidth=2, linestyle=':')
            subplot.plot(RingRadius, np.log10(Sigma),color='red', linewidth=2)     
            subplot.plot(RingRadius, np.log10(pc16),color='red', linestyle='--') 
            subplot.plot(RingRadius, np.log10(pc84),color='red', linestyle='--') 
            
            #WRITE OUTPUT      
            if(write_to_file==1):
                df = pd.DataFrame({'r_kpc':RingRadius, 'mean':np.log10(Sigma_mean), 'median':np.log10(Sigma), 
                                   'pc16':np.log10(pc16), 'pc84':np.log10(pc84)})         
                df.to_csv(Datadir + file_to_write + 'MilkyWay_Gradients_' + labels_to_write[i_property] + 
                          str(f'_z{ThisRedshiftList[ii]:0.2f}')+'.csv', index=False)
        
                #df = pd.read_csv(Datadir + file_to_write + 'MilkyWay_Gradients_'  + labels_to_write[i_property] + 
                #          str(f'_z{ThisRedshiftList[ii]:0.2f}')+'.csv')
                #subplot.plot(df['r_kpc'], df['mean'],color='black', linestyle='--') 
       
            if(i_property==0):              
                label='$\Sigma_{*}$'                
            if(i_property==1):
                label='$\Sigma_{\mathrm{HI}}$'
            if(i_property==2):
                label='$\Sigma_{\mathrm{H_2}}$'
            if(i_property==3):
                label='$\Sigma_{\mathrm{cold}}$'
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.75, y_percentage=0.84, 
                        color='black', xlog=0, ylog=0, label=label, 
                        fontsize=16, fontweight='normal') 
               
            
          
            if(i_property==2):                 
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.1, y_percentage=0.9, color='black', 
                            xlog=0, ylog=0, label=prefix_this_model, fontsize=13, fontweight='normal')             
                plot_label (subplot,'line',xlim,ylim,x_percentage=0.02,y_percentage=0.92,
                            color='red',x2_percentage=0.08,xlog=0,ylog=0,linestyle='-',linewidth=2)
          
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.1, y_percentage=0.81,
                            color='black', xlog=0, ylog=0,  label='Leroy 2008', fontsize=13, fontweight='normal')             
                plot_label (subplot,'line',xlim,ylim,x_percentage=0.02,y_percentage=0.83,
                            color='blue',x2_percentage=0.08,xlog=0,ylog=0,linestyle='-',linewidth=2)
                
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.1, y_percentage=0.72, 
                            color='black', xlog=0, ylog=0, label='BlueDiscs', fontsize=13, fontweight='normal')             
                plot_label (subplot,'line',xlim,ylim,x_percentage=0.02,y_percentage=0.74,
                            color='seagreen',x2_percentage=0.08,xlog=0,ylog=0,linestyle='-',linewidth=2)
            
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYF19_milkyway_sfr_and_gas_profiles.pdf')
    plt.close()

    return   
#end  milkyway_sfr_and_gas_profiles









 
     

def milkyway_gradients(ThisRedshiftList):
         
    #Model SELECTION
    ii=0   
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)                 
    G0_MR_unsel=G_MR[sel] 
   
    if(opt_detailed_enrichment==1): 
        G0_MR=G0_MR_unsel[(G0_MR_unsel['StellarMass']>0.) & (G0_MR_unsel['DiskMass']>0.) & 
                          ((G0_MR_unsel['MetalsDiskMass'][:,0] + G0_MR_unsel['MetalsDiskMass'][:,1] +
                            G0_MR_unsel['MetalsDiskMass'][:,2])>0.) & (G0_MR_unsel['Type']==0) & 
                          (np.log10(G0_MR_unsel['Sfr']/(G0_MR_unsel['StellarMass']*1e10/Hubble_h))>-10.5) & 
                          (np.log10(G0_MR_unsel['StellarMass']*1.e10/Hubble_h)>10.5) &
                          (np.log10(G0_MR_unsel['StellarMass']*1.e10/Hubble_h)<11.) &
                          (G0_MR_unsel['Vvir']>200.) & (G0_MR_unsel['Vvir']<235.) &
                          ((G0_MR_unsel['BulgeMass']/G0_MR_unsel['StellarMass'])<0.15)] 
    else:
        G0_MR=G0_MR_unsel[(G0_MR_unsel['StellarMass']>0.) & (G0_MR_unsel['DiskMass']>0.) & 
                          (G0_MR_unsel['MetalsDiskMass']>0.) & (G0_MR_unsel['Type']==0) & 
                          (np.log10(G0_MR_unsel['Sfr']/(G0_MR_unsel['StellarMass']*1e10/Hubble_h))>-10.5) & 
                          (np.log10(G0_MR_unsel['StellarMass']*1.e10/Hubble_h)>10.5) &
                          (np.log10(G0_MR_unsel['StellarMass']*1.e10/Hubble_h)<11.) &
                          (G0_MR_unsel['Vvir']>200.) & (G0_MR_unsel['Vvir']<235.) &
                          ((G0_MR_unsel['BulgeMass']/G0_MR_unsel['StellarMass'])<0.15)]
         
    fig = plt.figure(figsize=(two_one_size_small[0],two_one_size_small[1]))
    grid = gridspec.GridSpec(2,1)     
    grid.update(wspace=0.0, hspace=0.0)
    
    
    
    for i_property in range(0,2):
   
        subplot=plt.subplot(grid[i_property])  
        xlim=[0.0,15.0]
        ylim=[-0.8, 0.7]  
        subplot.set_ylim(ylim), subplot.set_xlim(xlim) 
        
        xlab='$r$[kpc]'                   
        subplot.set_xlabel(xlab, fontsize=14)
        if(i_property==0): 
            ylab='$\log_{10}$$(Z_{\mathrm{gas}}/Z_\odot)$'   
            plt.tick_params(axis='x', which='both', bottom='on', labelbottom='off')
        else:
            ylab='$\log_{10}$$(Z_*/Z_\odot)$'            
        subplot.set_ylabel(ylab, fontsize=14)
    
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(5))    
        subplot.xaxis.set_minor_locator(MultipleLocator(1))
        subplot.yaxis.set_major_locator(MultipleLocator(0.2))    
        subplot.yaxis.set_minor_locator(MultipleLocator(0.1))
    
        #OBSERVATIONS
        if(i_property==0):    
            file = Datadir + '/Moran12.txt'     
            Moran12 = Table.read(file, format='ascii')
            subplot.fill_between(Moran12['radius'],Moran12['log10_Zgas']-Moran12['err_down'],
                                 Moran12['log10_Zgas']+Moran12['err_up'], facecolor='lightgrey', 
                                 interpolate=True, alpha=0.8, edgecolor='black') 
        else:
            file = Datadir + 'Fu2013_Cepheids_Metals.txt'     
            Fu13_Ceph = Table.read(file, format='ascii')  
            subplot.errorbar(Fu13_Ceph['radius'], Fu13_Ceph['log10_Zstars'],
                             yerr=[Fu13_Ceph['err_down'],Fu13_Ceph['err_up']],
                             fmt='o', markersize=5, ecolor='blue', color='blue')    
        
        
        #MODEL
        mean_metallicity=np.zeros(RNUM,dtype=np.float32)  
        pc16=np.zeros(RNUM,dtype=np.float32)  
        pc84=np.zeros(RNUM,dtype=np.float32)  
        for jj in range(0,RNUM):
            #gas metallicity
            if(i_property==0):
                if(opt_detailed_enrichment==1):                  
                    MetalsRing=(G0_MR['MetalsColdGasRings'][:,jj,0] +
                                G0_MR['MetalsColdGasRings'][:,jj,1] +
                                G0_MR['MetalsColdGasRings'][:,jj,2])
                else:         
                    MetalsRing=G0_MR['MetalsColdGasRings'][:,jj]
                                               
                MassRing=G0_MR['ColdGasRings'][:,jj]
            #stellar metallicity    
            else:                
                if(opt_detailed_enrichment==1):                  
                    MetalsRing=(G0_MR['MetalsDiskMassRings'][:,jj,0] +
                                G0_MR['MetalsDiskMassRings'][:,jj,1] +
                                G0_MR['MetalsDiskMassRings'][:,jj,2])
                else:         
                    MetalsRing=G0_MR['MetalsDiskMassRings'][:,jj]                                               
                MassRing=G0_MR['DiskMassRings'][:,jj]
                
                if(opt_rings_in_bulges==1):
                    if(opt_detailed_enrichment==1):                  
                        MetalsRing+=(G0_MR['MetalsBulgeMassRings'][:,jj,0] +
                                    G0_MR['MetalsBulgeMassRings'][:,jj,1] +
                                    G0_MR['MetalsBulgeMassRings'][:,jj,2])
                    else:         
                        MetalsRing+=G0_MR['MetalsBulgeMassRings'][:,jj]                                               
                    MassRing+=G0_MR['BulgeMassRings'][:,jj]
                      
            metallicity=MetalsRing[MetalsRing>0.]/MassRing[MetalsRing>0.]/0.016         
            mean_metallicity[jj]=np.log10(np.median(metallicity))
            sel=metallicity>0.
            y_sorted = np.sort(metallicity[sel])
            if(len(metallicity[sel])>0.):
                pc16[jj] = y_sorted[int(16*len(metallicity[sel])/100)]      
                pc84[jj] = y_sorted[int(84*len(metallicity[sel])/100)] 
            else:
                pc16[jj] = 0.     
                pc84[jj] = 0.
            
        #endfor  
        subplot.plot(RingRadius,mean_metallicity,color='red',linewidth=2) 
        subplot.plot(RingRadius,np.log10(pc16),color='red',linewidth=2,linestyle='--') 
        subplot.plot(RingRadius,np.log10(pc84),color='red',linewidth=2,linestyle='--') 
        
        
        #LABELS       
        if(i_property==0):
            plot_label (subplot,'label',xlim,ylim,x_percentage=0.12,y_percentage=0.9, 
                        color='black',xlog=0,ylog=0,label=prefix_this_model,fontsize=13,fontweight='normal')             
            plot_label (subplot,'line',xlim,ylim,x_percentage=0.02,y_percentage=0.92,
                        color='red',x2_percentage=0.10,xlog=0,ylog=0,linestyle='-',linewidth=2)

            plot_label (subplot, 'label', xlim, ylim, 
                        x_percentage=0.12, y_percentage=0.82, color='black', xlog=0, ylog=0, 
                        label='Moran 2012', fontsize=13, fontweight='normal')             
            plot_label (subplot,'line',xlim,ylim,x_percentage=0.03,y_percentage=0.845,
                        color='lightgrey',x2_percentage=0.09,xlog=0,ylog=0,linestyle='-',linewidth=6)
        else:
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.08, y_percentage=0.9, 
                        color='blue', xlog=0, ylog=0, label='Fu 2013 (Cepheids)', 
                        fontsize=13, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.92, 
                        color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.04)
        
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYJ18_plots_milkyway_gradients.pdf')
    plt.close()
    
    return 
#end gas_milkyway_gradients


def gas_metallicity_gradients_mass_bins(ThisRedshiftList):
     
    ii=0   
    
    plot_color=['blue','green','red']
    
    fig = plt.figure(figsize=(two_one_size_small[0],two_one_size_small[1]))
    grid = gridspec.GridSpec(2, 1)
    grid.update(wspace=0.0, hspace=0.0)
    
    #OBSERVATIONS  
    Bresolin_Names=['NGC1512','NGC3621','M83','NGC4625']
    Bresolin_logMass=[11.4771,10.3010,]
    
    #haynes_MHI=haynes.data['HI']
    #haynes_Magr=haynes.data['mr']-5*(np.log10(haynes.data['distance']*1.e6)-1)
    '''Moran_Names=['G10019', 'G11845', 'G23315', 'G3505',   'G41323', 'G10218', 'G11956', 'G23408', 'G3509',  'G4137',
                 'G10358', 'G11989', 'G23419', 'G3519',   'G41783', 'G10404', 'G12002', 'G23443', 'G3524',  'G4195',
                 'G10447', 'G12025', 'G23450', 'G3537',   'G41969', 'G10813', 'G12069', 'G23453', 'G3593',  'G42013',
                 'G10817', 'G12318', 'G23496', 'G35981',  'G42025', 'G12460', 'G23685', 'G42140', 'G10827', 'G12970',
                 'G24064', 'G36307', 'G42141', 'G10831',  'G13037', 'G24094', 'G4216s', 'G10836', 'G13227', 'G24149',
                 'G3631',  'G4223',  'G10841', 'G13775',  'G24168', 'G3645',  'G4228',  'G10844', 'G13976', 'G24183',
                 'G3757',  'G4230',  'G10850', 'G14017',  'G24366', 'G3777',  'G4233',  'G10872', 'G14247', 'G24496',
                 'G3817',  'G4239',  'G10884',  'G14260', 'G25214', 'G3819',  'G42402', 'G10889', 'G14288', 'G4130',
                 'G25347', 'G3821',  'G47221', 'G10942',  'G14831', 'G25763', 'G38472', 'G48356', 'G10943', 'G14943',
                 'G26221', 'G3859',  'G51351', 'G10944',  'G15166', 'G26311', 'G38758', 'G51416', 'G10948', 'G15181',
                 'G26368', 'G3879',  'G51563', 'G10949',  'G15257', 'G26602', 'G3880',  'G51899', 'G10950', 'G16695',
                 'G26650', 'G38964', 'G52045', 'G11004',  'G17029', 'G26822', 'G39119', 'G52163', 'G11015', 'G17640',
                 'G27167', 'G39270', 'G52297', 'G11016b', 'G17659', 'G28168', 'G39548', 'G53795', 'G11016', 'G17684',
                 'G28461', 'G39567', 'G55745', 'G11019',  'G18202', 'G29487', 'G39595', 'G56307', 'G11071', 'G18335',
                 'G29510', 'G3971',  'G56312', 'G11087',  'G18421', 'G29555', 'G3976',  'G56320', 'G11109', 'G18900',
                 'G29596', 'G3977',  'G56375', 'G11112',  'G19132', 'G29699', 'G4017',  'G56612', 'G11120', 'G19672', 
                 'G29842', 'G4018',  'G56737', 'G11126',  'G1977',  'G29892', 'G40247', 'G57017', 'G11223', 'G19852',
                 'G30175', 'G40257', 'G6583',  'G11270',  'G19918', 'G30338', 'G4030',  'G7031',  'G11295', 'G19950',
                 'G30479', 'G40317', 'G7286',  'G11311',  'G20026', 'G30508', 'G4037',  'G8096',  'G11349', 'G20041',
                 'G30811', 'G4038',  'G8634',  'G11386',  'G20042', 'G3189',  'G4039',  'G9109',  'G11437', 'G20133',   
                 'G3261',  'G4040',  'G9463',  'G11488',  'G20144', 'G32937', 'G4048',  'G9551',  'G11513', 'G20183',   
                 'G3293',  'G40500', 'G9814',  'G11514',  'G20292', 'G3301',  'G40501', 'G9891',  'G11794', 'G21842',   
                 'G3435',  'G4054',  'G9948',  'G11808',  'G22999', 'G3439',  'G4057',  'G11817', 'G23120', 'G3465',    
                 'G4094',  'G11824', 'G23194', 'G3504'] 
    
    file=Datadir+'Moran/individual/'
    file=file+Moran_Names[10]+'spec_cat.fits'
    #file=file+'G10019spec_cat.fits'
    fits_table=fits.open(file)   
    cols = fits_table[1].columns
    cols.info()
    print("")
    print("")
    Moran = fits_table[1]
    print(Moran.data['METALLICITY'])
    print(Moran.data['R_ASEC_IN'])
     
        
    file=Datadir+'Moran/combined_cats/combined_cat_DR2DR3.fits'  
    fits_table=fits.open(file)   
    cols = fits_table[1].columns
    cols.info()    
    Moran = fits_table[1]
    print(Moran.data['GASS_ID'])'''
    
    
    for i_radius in range (0,1+1):
         
        subplot=plt.subplot(grid[i_radius])     
        if(i_radius==0):
            xlim=[0.0,15.0]
            xlab='$r$[kpc]'
            subplot.xaxis.set_label_position('top') 
            plt.tick_params(axis='x',which='both',top='on',labeltop='on',bottom='off',labelbottom='off')
            subplot.xaxis.set_major_locator(MultipleLocator(5.))    
            subplot.xaxis.set_minor_locator(MultipleLocator(1.))
        else:
            xlim=[-1.0,1.]
            xlab='$r/r_{d}$'
            plt.tick_params(axis='x',which='both',top='off',labeltop='off')
            subplot.xaxis.set_major_locator(MultipleLocator(0.5))    
            subplot.xaxis.set_minor_locator(MultipleLocator(0.1))
            
        ylim=[8.0, 9.1]     
        ylab='$12+\log_{10}(O/H)_{\mathrm{gas}}$'
        subplot.set_ylim(ylim),subplot.set_xlim(xlim)                 
        subplot.set_xlabel(xlab,fontsize=14), subplot.set_ylabel(ylab,fontsize=14)  
      
        subplot.yaxis.set_major_locator(MultipleLocator(0.2))    
        subplot.yaxis.set_minor_locator(MultipleLocator(0.1))
        
    
        median_metallicity=np.zeros(RNUM,dtype=np.float32)
        low_mass_limits=[9.5,10.0,10.5]       
        massbin=0.5   
    
        #***************       
        #*    MODEL    *
        #***************
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)                 
        G0_MR_unsel=G_MR[sel] 
    
        if(opt_detailed_enrichment==1):   
            G0_MR_unsel=G0_MR_unsel[((G0_MR_unsel['MetalsColdGas'][:,0] + 
                                      G0_MR_unsel['MetalsColdGas'][:,1] + 
                                      G0_MR_unsel['MetalsColdGas'][:,2])>.0) &
                                   (np.log10(G0_MR_unsel['Sfr']/(G0_MR_unsel['StellarMass']*1e10/Hubble_h))>-11.) & 
                                    (G0_MR_unsel['BulgeMass']/G0_MR_unsel['StellarMass']<0.15)]           
        else:
            G0_MR_unsel=G0_MR_unsel[(G0_MR_unsel['MetalsColdGas']>.0) & 
                                    (np.log10(G0_MR_unsel['Sfr']/(G0_MR_unsel['StellarMass']*1e10/Hubble_h))>-11.) & 
                                    (G0_MR_unsel['BulgeMass']/G0_MR_unsel['StellarMass']<0.15)] #& 
                                    #(G0_MR_unsel['BulgeMass']/G0_MR_unsel['StellarMass']<.3)]
       
        for kk in range(0,len(low_mass_limits)):
      
            G0_MR=G0_MR_unsel[(np.log10(G0_MR_unsel['StellarMass']*1.e10/Hubble_h) > low_mass_limits[kk]) & 
                              (np.log10(G0_MR_unsel['StellarMass']*1.e10/Hubble_h) < low_mass_limits[kk]+massbin)]
              
            NGals=len(G0_MR)
            x_variable=np.zeros(int(RNUM*NGals),dtype=np.float32) 
            y_variable=np.zeros(int(RNUM*NGals),dtype=np.float32)   
            
            for jj in range(0,RNUM):   
               
                ColdGasRing=G0_MR['ColdGasRings'][:,jj]
                if(opt_individual_elements==1):                  
                    '''MetalsColdGasRing=(G0_MR['MetalsColdGasRings'][:,jj,0] +
                                       G0_MR['MetalsColdGasRings'][:,jj,1] +                                       
                                       G0_MR['MetalsColdGasRings'][:,jj,2])
                    sel=MetalsColdGasRing>0.
                    Metallicity=MetalsColdGasRing[sel]/ColdGasRing[sel]/0.0134*10**8.69'''
                    N_H=G0_MR['ColdGasRings_elements'][:,jj,0]/1.
                    N_O=G0_MR['ColdGasRings_elements'][:,jj,4]/16.
                    sel=(N_O>0.) & (N_H>0.) 
                    Metallicity=1e12*(N_O[sel]/N_H[sel])                     
                else:            
                    MetalsColdGasRing=G0_MR['MetalsColdGasRings'][:,jj] 
                    sel=MetalsColdGasRing>0.
                    Metallicity=MetalsColdGasRing[sel]/ColdGasRing[sel]/0.0134*10**8.69
                
                if(i_radius==0):                    
                    median_metallicity[jj]=np.log10(np.median(Metallicity))                 
                else:
                    #aux needed to because len(Metallicity)!= NGals
                    aux_x=np.zeros(int(NGals),dtype=np.float32)
                    aux_y=np.zeros(int(NGals),dtype=np.float32)
                    aux_x[sel]=np.log10(RingRadius[jj]/(G0_MR['StellarHalfMassRadius'][sel]*1000./Hubble_h)) 
                    aux_y[sel]=Metallicity #no [sel] here because the selection was applied above
                    x_variable[NGals*jj:NGals*(jj+1)]=aux_x     
                    y_variable[NGals*jj:NGals*(jj+1)]=aux_y
                              
            #metallicity versus physical radius
            if(i_radius==0):               
                subplot.plot(RingRadius, median_metallicity,color=plot_color[kk], linewidth=2)  
            #metallicity versus physical radius/rd    
            else:  
                bin=0.1
                sel=y_variable>0. 
                if(len(y_variable[sel])>0.):
                    (x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles (bin, xlim[0], xlim[1], 
                                                                                      x_variable[sel], y_variable[sel])
                    #print(np.log10(median))
                    subplot.plot(x_binned, np.log10(median), color=plot_color[kk], linewidth=2)
            
                #median_radius=np.median(G0_MR['GasDiskRadius'][G0_MR['GasDiskRadius']>0.]*1000./Hubble_h)
                #Rings=np.log10(RingRadius/median_radius)
                #subplot.plot(Rings, median_metallicity, color=plot_color[kk], linewidth=2)   
        
        
            #labels
            if(i_radius==1):
                label="%0.1f" % low_mass_limits[kk] + "<$M_{\star}[M_{\odot}]$<" + "%0.1f" % (low_mass_limits[kk]+massbin) 
                plot_label (subplot, 'label',xlim,ylim,x_percentage=0.15,y_percentage=0.05+(kk*0.09), 
                            color='black',xlog=0,ylog=0,label=label,fontsize=12,fontweight='normal')             
                plot_label (subplot,'line',xlim,ylim,x_percentage=0.05,y_percentage=0.075+(kk*0.09),
                            color=plot_color[kk],x2_percentage=0.12,xlog=0,ylog=0,linestyle='-',linewidth=2)
           
        
        
            #****************       
            #* OBSERVATIONS *
            #****************
            if((low_mass_limits[kk]==10.) & (i_radius==0)):
                file = Datadir + '/bresolin12_gradients_NGC3621.txt'       
                obs = Table.read(file, format='ascii') 
                subplot.plot([obs['x_1'],obs['x_2']], [obs['y_1'],obs['y_2']], 
                             color=plot_color[kk],linestyle='--',linewidth=2)
                plot_label (subplot, 'label',xlim,ylim,x_percentage=0.05,y_percentage=0.32, 
                            color=plot_color[kk],xlog=0,ylog=0,label='B12-NGC3621',
                            fontsize=12,fontweight='normal')
              
                
            if((low_mass_limits[kk]==10.5) & (i_radius==0)):
                file = Datadir + '/bresolin12_gradients_NGC1512.txt'       
                obs = Table.read(file, format='ascii') 
                subplot.plot([obs['x_1'],obs['x_2']], [obs['y_1'],obs['y_2']], 
                             color=plot_color[kk],linestyle='--',linewidth=2)
                plot_label (subplot, 'label',xlim,ylim,x_percentage=0.6,y_percentage=0.32, 
                            color=plot_color[kk],xlog=0,ylog=0,label='B12-NGC1512',
                            fontsize=12,fontweight='normal')
                
            '''if((low_mass_limits[kk]==10.5) & (i_radius==0)):
                file = Datadir + '/bresolin12_gradients_M83.txt'       
                obs = Table.read(file, format='ascii') 
                subplot.plot([obs['x_1'],obs['x_2']], [obs['y_1'],obs['y_2']], color=plot_color[kk],linestyle=':')
            if((low_mass_limits[kk]==10.5) & (i_radius==0)):
                file = Datadir + '/bresolin12_gradients_NGC4625.txt'       
                obs = Table.read(file, format='ascii') 
                subplot.plot([obs['x_1'],obs['x_2']], [obs['y_1'],obs['y_2']], color=plot_color[kk],linestyle=':')'''
        
        
        #endfor -> MASS
    #endfor -> i_radius
              
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYW17_plots_gas_metallicity_gradients_mass_bins.pdf')
    plt.close()
    
    return 
    
#end gas_metallicity_gradients



def stellar_metallicity_gradients_mass_bins(ThisRedshiftList):
     
    ii=0   
    
    plot_color=['purple','blue','lightblue','green','orange','red','brown']    
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
     
   
    subplot=plt.subplot()     
       
    xlim=[0.0,3.]
    xlab='$r/r_{d}$'
    plt.tick_params(axis='x',which='both',top='off',labeltop='off')
    subplot.xaxis.set_major_locator(MultipleLocator(1.0))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.1))

    ylim=[-1.0, 1.0]     
    ylab='$\mathrm{log_{10}}(Z_*/Z_{\odot})$'
    subplot.set_ylim(ylim),subplot.set_xlim(xlim)                 
    subplot.set_xlabel(xlab,fontsize=14), subplot.set_ylabel(ylab,fontsize=14)  

    subplot.yaxis.set_major_locator(MultipleLocator(0.2))    
    subplot.yaxis.set_minor_locator(MultipleLocator(0.1))


    median_metallicity=np.zeros(RNUM,dtype=np.float32)
    mass_bins=[9.1,9.6,10.1,10.6,10.9,11.2,11.5,11.8] 


    
    #Model
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)                 
    G0_MR_unsel=G_MR[sel] 

    if(opt_detailed_enrichment==1):   
        G0_MR_unsel=G0_MR_unsel[((G0_MR_unsel['MetalsDiskMass'][:,0] + 
                                  G0_MR_unsel['MetalsDiskMass'][:,1] + 
                                  G0_MR_unsel['MetalsDiskMass'][:,2])>.0)] #&
                                #(np.log10(G0_MR_unsel['Sfr']/(G0_MR_unsel['StellarMass']*1e10/Hubble_h))>-10.5) & 
                               # (G0_MR_unsel['BulgeMass']/G0_MR_unsel['StellarMass']<0.15)]
    else:
        G0_MR_unsel=G0_MR_unsel[(G0_MR_unsel['MetalsDiskMass']>.0) & 
                                (np.log10(G0_MR_unsel['Sfr']/(G0_MR_unsel['StellarMass']*1e10/Hubble_h))>-10.5) & 
                                (G0_MR_unsel['BulgeMass']/G0_MR_unsel['StellarMass']<0.15)] #& 
        #(G0_MR_unsel['BulgeMass']/G0_MR_unsel['StellarMass']<.3)]

    for kk in range(0,len(mass_bins)-1):

        G0_MR=G0_MR_unsel[(np.log10(G0_MR_unsel['StellarMass']*1.e10/Hubble_h) > mass_bins[kk]) & 
                          (np.log10(G0_MR_unsel['StellarMass']*1.e10/Hubble_h) < mass_bins[kk+1])]

       
        for jj in range(0,RNUM):             
            if(opt_detailed_enrichment==1):                  
                MetalsColdGasRing=(G0_MR['MetalsDiskMassRings'][:,jj,0] +
                                   G0_MR['MetalsDiskMassRings'][:,jj,1] +
                                   G0_MR['MetalsDiskMassRings'][:,jj,2])
            else:            
                MetalsColdGasRing=G0_MR['MetalsDiskMassRings'][:,jj]                

            ColdGasRing=G0_MR['DiskMassRings'][:,jj]
            metallicity=MetalsColdGasRing[MetalsColdGasRing>0.]/ColdGasRing[MetalsColdGasRing>0.]/0.02  
            if(len(metallicity)>0.):
                median_metallicity[jj]=np.log10(np.median(metallicity))           

       
        median_radius=np.median(G0_MR['StellarHalfLightRadius'][G0_MR['StellarHalfLightRadius']>0.] *
                                1000./Hubble_h)
        Rings=np.log10(RingRadius/median_radius)
        Rings=RingRadius/median_radius
        subplot.plot(Rings, median_metallicity, color=plot_color[kk], linewidth=2)   
        
        #labels            
        label="%0.1f" % mass_bins[kk] + "<$M_{\star}[M_{\odot}]$<" + "%0.1f" % mass_bins[kk+1] 
        plot_label (subplot, 'label',xlim,ylim,x_percentage=0.15,y_percentage=0.05+(kk*0.09), 
                    color='black',xlog=0,ylog=0,label=label,fontsize=12,fontweight='normal')             
        plot_label (subplot,'line',xlim,ylim,x_percentage=0.05,y_percentage=0.075+(kk*0.09),
                    color=plot_color[kk],x2_percentage=0.12,xlog=0,ylog=0,linestyle='-',linewidth=2)

        
        
          
        
        #endfor -> MASS
    #endfor -> i_radius
              
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYW17_plots_stellar_metallicity_gradients_mass_bins.pdf')
    plt.close()
    
    return    
    
#end stellar_metallicity_gradients


'''
  for kk in range(0,RNUM):
                if(kk==0):
                    r_in=0.
                else:
                    r_in=(RingRadius[kk-1]+RingRadius[kk])/2.
                        
                if(kk==RNUM-1):
                    r_out=RingRadius[kk]+((RingRadius[kk]-RingRadius[kk-1])/2.)
                else:
                    r_out=(RingRadius[kk]+RingRadius[kk+1])/2.
                        
                if(i_property==0):            
                        
                    r_bulge=G0_MR['BulgeSize']*1000./Hubble_h #From Mpc/h to Kpc
                    if(opt_rings_in_bulges==1):
                        BulgeMass_this_ring=G0_MR['BulgeMassRings'][:,kk]*1e10/Hubble_h
                    else:
                        BulgeMass_this_ring=((G0_MR['BulgeMass']*1e10/Hubble_h)*
                                             ((1+(r_in/r_bulge))**(-1.) - (1.+(r_out/r_bulge))**(-1.)) )  
                   
                    Mass=G0_MR['DiskMassRings'][:,kk]*1e10/Hubble_h
                    sel=(r_bulge>0.)
                    Mass[sel]+=BulgeMass_this_ring[sel]                    
                       
                    
                y_variable=Mass/(3.14*(r_out**2-r_in**2)*1e6)    
                sel=y_variable>0.               
                if(len(y_variable[sel])>0.):
                    Sigma[kk]=np.median(y_variable[sel])
                    Sigma_mean[kk]=np.mean(y_variable[sel])
                    y_sorted = np.sort(y_variable[sel])
                    pc16[kk] = y_sorted[int(16*len(y_variable[sel])/100)]      
                    pc84[kk] = y_sorted[int(84*len(y_variable[sel])/100)]  
      '''  


def gradients_enci(ThisRedshiftList):
    
       
    '''fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()

    xlim=[9.0, 12.0]    
    ylim=[-3.5, 3.0]         
    subplot.set_ylim(ylim),subplot.set_xlim(xlim)

    xlab='Mass'
    ylab='SFR'  
    subplot.set_xlabel(xlab,fontsize=14), subplot.set_ylabel(ylab,fontsize=14)       

    subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
    subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
    subplot.xaxis.set_major_locator(MultipleLocator(0.5))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.1)) 

    i_z=0
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)                 
    G0_MR = G_MR[sel] 
    G0_MR = G0_MR[(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)>8.5) & (G0_MR['Sfr']>0) & (G0_MR['Type']==0)]
       
    Mass = np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)
    SFR = np.log10(G0_MR['Sfr'])
       
    bin=[0.2,0.2]
    Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]
    Ngals=len(Mass)    
    H, xedges, yedges = np.histogram2d(Mass, SFR, bins=Nbins)            
    extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
    plt.subplots_adjust(bottom=0.15, left=0.15)        
    mylevels = np.linspace(1., Nbins[0], Nbins[0])*Ngals/(Nbins[0]**2/3.)        
    H = zoom(H, 20)        
    cont=plt.contourf(H.transpose()[::], origin='lower', cmap='Greys_r', levels=mylevels, extent=extent)    
       
    a = 1.0
    b = -9.65
    
    xx = np.arange(xlim[0],xlim[1],0.1)
    yy = a*xx+b
    subplot.plot(xx,yy, color='blue')
    
    yy = a*xx+b+0.1
    subplot.plot(xx,yy, color='blue',linestyle='--')
    
    yy = a*xx+b+0.3
    subplot.plot(xx,yy, color='blue',linestyle='--')
    
    yy = a*xx+b+0.5
    subplot.plot(xx,yy, color='blue',linestyle='--')
    
            
    
    
    plt.savefig('./fig/HYJ18_main_sequence.pdf')
    plt.close()
        
        
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()

    xlim=[7.0, 10.0]    
    ylim=[-3.5, 0.0]         
    subplot.set_ylim(ylim),subplot.set_xlim(xlim)

    xlab='Mass'
    ylab='SFR'  
    subplot.set_xlabel(xlab,fontsize=14), subplot.set_ylabel(ylab,fontsize=14)       

    subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
    subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
    subplot.xaxis.set_major_locator(MultipleLocator(0.5))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.1)) 

    i_z=0
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)                 
    G0_MR = G_MR[sel] 
    G0_MR = G0_MR[(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)>8.5) & (G0_MR['Sfr']>0) & (G0_MR['Type']==0)]
   
    area = np.zeros(RNUM,dtype=np.float32)
    for jj in range(0,RNUM):
        if(jj==0):            
            area[jj]=(3.14*RingRadius[jj]*RingRadius[jj])
        else:
            area[jj]=(3.14*(RingRadius[jj]*RingRadius[jj]-RingRadius[jj-1]*RingRadius[jj-1]))
    
    #print(G0_MR['DiskMassRings'][:,:]*area.shape)
    Mass = np.log10(G0_MR['DiskMassRings'][:,10]/area[10]*1.e10/Hubble_h)
    SFR = np.log10(G0_MR['SfrRings'][:,10]/area[10])
    
    #for jj in range(0,RNUM):
    #     subplot.scatter(Mass[:1000,jj], SFR[:1000,jj], s=2)
    Mass_spaxels = Mass.flatten()
    SFR_spaxels = SFR.flatten()
    
    sel = (Mass_spaxels>7.0) & (SFR_spaxels>-4.0)
    
    bin=[0.2,0.2]
    Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]
    Ngals=len(Mass_spaxels)/100000.   
    H, xedges, yedges = np.histogram2d(Mass_spaxels[sel], SFR_spaxels[sel], bins=Nbins)            
    extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
    plt.subplots_adjust(bottom=0.15, left=0.15)        
    mylevels = np.linspace(1., Nbins[0], Nbins[0])*Ngals/(Nbins[0]**2/3.)        
    H = zoom(H, 20)        
    cont=plt.contourf(H.transpose()[::], origin='lower', cmap='Greys_r', levels=mylevels, extent=extent)    
       
    #subplot.scatter(Mass_spaxels[:5000], SFR_spaxels[:5000], s=2)
    
    a = 1.1
    b = -10.7
    
    xx = np.arange(xlim[0],xlim[1],0.1)
    yy = a*xx+b
    subplot.plot(xx,yy, color='black')
    
    a = [1.3, 1.17, 1.14, 1.14, 1.14, 1.1, 1.1, 0.9, 0.9, 0.85, 0.85, 0.85]
    b = [-12.45, -11.3, -11.05, -11.05, -11.05, -10.7, -10.6, -8.95, -8.95, -8.47, -8.47, -8.47]
    
    a = 0.85
    b = -8.47
    
    xx = np.arange(xlim[0],xlim[1],0.1)
    yy = a*xx+b
    subplot.plot(xx,yy, color='blue')
 
            
    
    
    plt.savefig('./fig/HYJ18_main_sequence_spaxels.pdf')
    plt.close()'''
            
        
        
        
        
    plot_color=['purple','red', 'orange', 'black', 'green', 'blue'] 
       
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    
    global_a = 1.0
    global_b = -9.65
    
    spaxel_a = [1.3, 1.17, 1.14, 1.14, 1.14, 1.1, 1.1, 0.9, 0.9, 0.85, 0.85, 0.85]
    spaxel_b = [-12.45, -11.3, -11.05, -11.05, -11.05, -10.7, -10.6, -8.95, -8.95, -8.47, -8.47, -8.47]
    
    area = np.zeros(RNUM,dtype=np.float32)
    for jj in range(0,RNUM):
        if(jj==0):            
            area[jj]=(3.14*RingRadius[jj]*RingRadius[jj])
        else:
            area[jj]=(3.14*(RingRadius[jj]*RingRadius[jj]-RingRadius[jj-1]*RingRadius[jj-1]))
    
    delta_low = [-5.0,-1.0,-0.3,-0.1,0.0,0.3]
    delta_high = [-1.0, -0.3,0.0,0.1, 0.3,1.0]
        
    for delta_i in range(0,len(delta_low)):           

        subplot=plt.subplot()

        xlim=[0., 2.]
        ylim=[-3, 0.5]          
        subplot.set_ylim(ylim),subplot.set_xlim(xlim)

        xlab='$r/r_e$'
        ylab='$\log_{10}(\Sigma_{SFR}[\mathrm{M}_{\odot}\mathrm{yr}^{-1}\mathrm{Kpc^{-2}}])$'  
        subplot.set_xlabel(xlab,fontsize=14), subplot.set_ylabel(ylab,fontsize=14)       

        subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
        subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
        subplot.xaxis.set_major_locator(MultipleLocator(0.5))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.1)) 

        i_z=0
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)                 
        G0_MR_unsel = G_MR[sel] 
        G0_MR_unsel = G0_MR_unsel[(np.log10(G0_MR_unsel['StellarMass']*1.e10/Hubble_h)>10.5) & 
                                  (np.log10(G0_MR_unsel['StellarMass']*1.e10/Hubble_h)<11.0) & 
                                  (G0_MR_unsel['Sfr']>0.) & (G0_MR_unsel['Type']==0)]  
      
        Mass = np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)
        SFR = np.log10(G0_MR_unsel['Sfr'])
        
        #global deviation from main squence
        G0_MR=G0_MR_unsel[(SFR>Mass*global_a+global_b+delta_low[delta_i]) & (SFR<Mass*global_a+global_b+delta_high[delta_i])]
                          
        
        NGals=len(G0_MR)            
        x_variable=np.zeros(int(RNUM*NGals),dtype=np.float32)             
        y_variable=np.zeros(int(RNUM*NGals),dtype=np.float32)   
        new_x_var = np.arange(xlim[0],xlim[1]+1.0,0.1)

        re = G0_MR['StellarDiskRadius']*1000./Hubble_h
        #loop on radial bins
        for jj in range(0,RNUM):  
            x_variable[NGals*jj:NGals*(jj+1)] = RingRadius[jj]/re
            #SFR_this_ring=np.log10(G0_MR['SfrRings'][:,jj]/area[jj])             
            #SFR_Main_Sequence = spaxel_a[delta_i] * np.log10(G0_MR['DiskMassRings'][:,jj]/area[jj]*1.e10/Hubble_h)+ spaxel_b[delta_i]            
            #y_variable[NGals*jj:NGals*(jj+1)]= np.log10(SFR_this_ring - SFR_Main_Sequence)
            
            #y_variable[NGals*jj:NGals*(jj+1)]= np.log10(G0_MR['SfrRings'][:,jj]/area[jj]) 
            y_variable[NGals*jj:NGals*(jj+1)]= np.log10(G0_MR['SfrRings'][:,jj] / area[jj] * ((0.2*re)**2)) 
                                     
        #endfor RNUM      
        N_random_gals = 1000
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
            yy = y_variable[slice_ii]

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
        subplot.plot(x_binned[1:], median[1:], color=plot_color[delta_i], linewidth=2, linestyle='-')
        #end loop on mass bins    

        
    #end loop on properties to plot
    
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYJ18_gradients_ellison.pdf')
    plt.close()
    
    return 
#end gradients_enci

def gradients_ellison(ThisRedshiftList):
    
       
    '''fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()

    xlim=[9.0, 12.0]    
    ylim=[-3.5, 3.0]         
    subplot.set_ylim(ylim),subplot.set_xlim(xlim)

    xlab='Mass'
    ylab='SFR'  
    subplot.set_xlabel(xlab,fontsize=14), subplot.set_ylabel(ylab,fontsize=14)       

    subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
    subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
    subplot.xaxis.set_major_locator(MultipleLocator(0.5))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.1)) 

    i_z=0
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)                 
    G0_MR = G_MR[sel] 
    G0_MR = G0_MR[(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)>8.5) & (G0_MR['Sfr']>0) & (G0_MR['Type']==0)]
       
    Mass = np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)
    SFR = np.log10(G0_MR['Sfr'])
       
    bin=[0.2,0.2]
    Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]
    Ngals=len(Mass)    
    H, xedges, yedges = np.histogram2d(Mass, SFR, bins=Nbins)            
    extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
    plt.subplots_adjust(bottom=0.15, left=0.15)        
    mylevels = np.linspace(1., Nbins[0], Nbins[0])*Ngals/(Nbins[0]**2/3.)        
    H = zoom(H, 20)        
    cont=plt.contourf(H.transpose()[::], origin='lower', cmap='Greys_r', levels=mylevels, extent=extent)    
       
    a = 1.0
    b = -9.65
    
    xx = np.arange(xlim[0],xlim[1],0.1)
    yy = a*xx+b
    subplot.plot(xx,yy, color='blue')
    
    yy = a*xx+b+0.1
    subplot.plot(xx,yy, color='blue',linestyle='--')
    
    yy = a*xx+b+0.3
    subplot.plot(xx,yy, color='blue',linestyle='--')
    
    yy = a*xx+b+0.5
    subplot.plot(xx,yy, color='blue',linestyle='--')
    
            
    
    
    plt.savefig('./fig/HYJ18_main_sequence.pdf')
    plt.close()
        
        
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()

    xlim=[7.0, 10.0]    
    ylim=[-3.5, 0.0]         
    subplot.set_ylim(ylim),subplot.set_xlim(xlim)

    xlab='Mass'
    ylab='SFR'  
    subplot.set_xlabel(xlab,fontsize=14), subplot.set_ylabel(ylab,fontsize=14)       

    subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
    subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
    subplot.xaxis.set_major_locator(MultipleLocator(0.5))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.1)) 

    i_z=0
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)                 
    G0_MR = G_MR[sel] 
    G0_MR = G0_MR[(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)>8.5) & (G0_MR['Sfr']>0) & (G0_MR['Type']==0)]
   
    area = np.zeros(RNUM,dtype=np.float32)
    for jj in range(0,RNUM):
        if(jj==0):            
            area[jj]=(3.14*RingRadius[jj]*RingRadius[jj])
        else:
            area[jj]=(3.14*(RingRadius[jj]*RingRadius[jj]-RingRadius[jj-1]*RingRadius[jj-1]))
    
    #print(G0_MR['DiskMassRings'][:,:]*area.shape)
    Mass = np.log10(G0_MR['DiskMassRings'][:,10]/area[10]*1.e10/Hubble_h)
    SFR = np.log10(G0_MR['SfrRings'][:,10]/area[10])
    
    #for jj in range(0,RNUM):
    #     subplot.scatter(Mass[:1000,jj], SFR[:1000,jj], s=2)
    Mass_spaxels = Mass.flatten()
    SFR_spaxels = SFR.flatten()
    
    sel = (Mass_spaxels>7.0) & (SFR_spaxels>-4.0)
    
    bin=[0.2,0.2]
    Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]
    Ngals=len(Mass_spaxels)/100000.   
    H, xedges, yedges = np.histogram2d(Mass_spaxels[sel], SFR_spaxels[sel], bins=Nbins)            
    extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
    plt.subplots_adjust(bottom=0.15, left=0.15)        
    mylevels = np.linspace(1., Nbins[0], Nbins[0])*Ngals/(Nbins[0]**2/3.)        
    H = zoom(H, 20)        
    cont=plt.contourf(H.transpose()[::], origin='lower', cmap='Greys_r', levels=mylevels, extent=extent)    
       
    #subplot.scatter(Mass_spaxels[:5000], SFR_spaxels[:5000], s=2)
    
    a = 1.1
    b = -10.7
    
    xx = np.arange(xlim[0],xlim[1],0.1)
    yy = a*xx+b
    subplot.plot(xx,yy, color='black')
    
    a = [1.3, 1.17, 1.14, 1.14, 1.14, 1.1, 1.1, 0.9, 0.9, 0.85, 0.85, 0.85]
    b = [-12.45, -11.3, -11.05, -11.05, -11.05, -10.7, -10.6, -8.95, -8.95, -8.47, -8.47, -8.47]
    
    a = 0.85
    b = -8.47
    
    xx = np.arange(xlim[0],xlim[1],0.1)
    yy = a*xx+b
    subplot.plot(xx,yy, color='blue')
 
            
    
    
    plt.savefig('./fig/HYJ18_main_sequence_spaxels.pdf')
    plt.close()'''
            
        
        
        
        
    plot_color=['red', 'orange', 'yellow', 'green', 'blue', 'darkblue','purple'] 
       
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    
    global_a = 1.0
    global_b = -9.65
    
    spaxel_a = [1.3, 1.17, 1.14, 1.14, 1.14, 1.1, 1.1, 0.9, 0.9, 0.85, 0.85, 0.85]
    spaxel_b = [-12.45, -11.3, -11.05, -11.05, -11.05, -10.7, -10.6, -8.95, -8.95, -8.47, -8.47, -8.47]
    
    area = np.zeros(RNUM,dtype=np.float32)
    for jj in range(0,RNUM):
        if(jj==0):            
            area[jj]=(3.14*RingRadius[jj]*RingRadius[jj])
        else:
            area[jj]=(3.14*(RingRadius[jj]*RingRadius[jj]-RingRadius[jj-1]*RingRadius[jj-1]))
    
    delta_low = [-3.0,-0.5,-0.3,-0.1,0.1,0.3,0.5]
    delta_high = [-0.5,-0.3,-0.1,0.1,0.3,0.5,3.0]
    
    #delta_low = [-0.1]
    #delta_high = [0.1]
            
    for delta_i in range(0,len(delta_low)):           

        subplot=plt.subplot()

        xlim=[0., 10.]
        ylim=[-4, 4.]          
        subplot.set_ylim(ylim),subplot.set_xlim(xlim)

        xlab='$\log_{10}(r \mathrm{[kpc]})$'
        ylab='$\log_{10}(\Sigma_{SFR}[\mathrm{yr}^{-1}\mathrm{Kpc^{-2}}])$'  
        subplot.set_xlabel(xlab,fontsize=14), subplot.set_ylabel(ylab,fontsize=14)       

        subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
        subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
        subplot.xaxis.set_major_locator(MultipleLocator(5.0))    
        subplot.xaxis.set_minor_locator(MultipleLocator(1.0)) 

        i_z=0
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)                 
        G0_MR_unsel = G_MR[sel] 
        G0_MR_unsel = G0_MR_unsel[(np.log10(G0_MR_unsel['StellarMass']*1.e10/Hubble_h)>8.5) & 
                                  (G0_MR_unsel['Sfr']>0.) & (G0_MR_unsel['Type']==0)]  
      
        
    
        Mass = np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)
        SFR = np.log10(G0_MR_unsel['Sfr'])
        
        #global deviation from main squence
        G0_MR=G0_MR_unsel[(SFR>Mass*global_a+global_b+delta_low[delta_i]) & (SFR<Mass*global_a+global_b+delta_high[delta_i])]
        
        NGals=len(G0_MR)            
        x_variable=np.zeros(int(RNUM*NGals),dtype=np.float32)             
        y_variable=np.zeros(int(RNUM*NGals),dtype=np.float32)   
        new_x_var = np.arange(xlim[0],xlim[1]+1.0,0.05)

        resolved_MS = np.array([-4.870e-01, -7.104e-01, -8.779e-01, -1.045e+00, -1.225e+00, -1.490e+00, 
                                -1.856e+00, -2.426e+00, -6.482e+00, -5.002e+01, -4.590e+02, -4.060e+03])
        re = G0_MR['StellarDiskRadius']*1000./Hubble_h
        #loop on radial bins
        for jj in range(0,RNUM):  
            x_variable[NGals*jj:NGals*(jj+1)] = RingRadius[jj]
            SFR_this_ring=G0_MR['SfrRings'][:,jj]/area[jj]            
            #SFR_Main_Sequence = spaxel_a[delta_i] * np.log10(G0_MR['DiskMassRings'][:,jj]/area[jj]*1.e10/Hubble_h)+ spaxel_b[delta_i]  
            
            #y_variable[NGals*jj:NGals*(jj+1)]= np.log10(10**resolved_MS[jj]-SFR_this_ring)
            y_variable[NGals*jj:NGals*(jj+1)]= - resolved_MS[jj] + np.log10(SFR_this_ring)
            #y_variable[NGals*jj:NGals*(jj+1)]= np.log10(SFR_this_ring)
            
        #endfor RNUM      
        N_random_gals = 1000
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
            yy = y_variable[slice_ii]

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
        
        '''for kk in range(0, RNUM):
            print(f'{RingRadius[kk]:0.3f}')
        for kk in range(0, len(x_binned)):
            print(f'{x_binned[kk]:0.3f} {median[kk]:0.3e}')'''
            
        subplot.plot(x_binned, median, color=plot_color[delta_i], linewidth=2, linestyle='-')
        #end loop on mass bins    

        
    #end loop on properties to plot
    
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYJ18_gradients_ellison.pdf')
    plt.close()
    
    return 
#end gradients_ellison

'''def gradients_mean_evo (ThisRedshiftList):

    plot_color=['orange', 'red'] 
    mass_low = 11.0
    mass_high = 11.5
    
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
        
    subplot=plt.subplot()

    xlim=[-0.4, 1.5]
    ylim=[-3, 2.]    
    #ylim=[-14, -8.5]         
    subplot.set_ylim(ylim),subplot.set_xlim(xlim)

    xlab='$\log_{10}(r \mathrm{[kpc]})$'        
    ylab='$\log_{10}(\Sigma_{SFR}[M_{\odot} \mathrm{yr}^{-1}\mathrm{Kpc^{-2}}])$'         
    #if(i_z==1):
    #    ylab=''
    #    plt.tick_params(axis='y', which='both', left=True, labelleft=False)

    subplot.set_xlabel(xlab,fontsize=14), subplot.set_ylabel(ylab,fontsize=14)       

    subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
    subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
    subplot.xaxis.set_major_locator(MultipleLocator(0.5))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.1)) 
    
    
    for i_z in range (0,len(ThisRedshiftList)):   
    
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)                 
        G0_MR_unsel = G_MR[sel] 
        G0_MR = G0_MR_unsel[(G0_MR_unsel['Sfr']>0.) & (G0_MR_unsel['Type']==0) &
                            (np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)>mass_low) &
                            (np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)<mass_high)]  
        print(ThisRedshiftList[i_z])
      
        #if(i_z==0):
        #    SFR_cut = -10.5
        #else:
        #    SFR_cut = -10.0

        #if(quenched_state==0):
        #    sel = np.log10(G0_MR_unsel['Sfr']/(G0_MR_unsel['StellarMass']*1.e10/Hubble_h))>SFR_cut
        #else:
        #    sel = np.log10(G0_MR_unsel['Sfr']/(G0_MR_unsel['StellarMass']*1.e10/Hubble_h))<SFR_cut

        #G0_MR_unsel = G0_MR_unsel[sel]    
        

        NGals=len(G0_MR)   
        print(NGals)
        x_variable=np.zeros(int(RNUM*NGals),dtype=np.float32)             
        y_variable=np.zeros(int(RNUM*NGals),dtype=np.float32)   
        new_x_var = np.arange(xlim[0],xlim[1]+1.0,0.2)


        #loop on radial bins
        for jj in range(0,RNUM):  
            x_variable[NGals*jj:NGals*(jj+1)]=np.log10(RingRadius[jj])  
            SFR_this_ring=G0_MR['SfrRings'][:,jj]
            #1e6 -> from kpc^2 to pc^2
            if(jj==0):
                y_variable[NGals*jj:NGals*(jj+1)]=SFR_this_ring/(3.14*RingRadius[0]**2) 
            else:
                y_variable[NGals*jj:NGals*(jj+1)]=SFR_this_ring/(3.14*(RingRadius[jj]**2-RingRadius[jj-1]**2))                            
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
        
        subplot.plot(x_binned, median, color='red', linewidth=2, linestyle='-')

        #WRITE OUTPUT      
        if(write_to_file==1):
            df = pd.DataFrame({'log10_r':x_binned, 'log10_SigmaSFR':median})
            file = Datadir+file_to_write+'gradients_mean_evo'+str(f'_M{mass_low:0.2f}') + \
                   str(f'_{mass_high:0.2f}') + str(f'_Qstate{quenched_state:d}') + \
                   str(f'_z{ThisRedshiftList[i_z]:0.2f}')+'.csv'
            df.to_csv(file,index=False)
            #df = pd.read_csv(file)
            #subplot.plot(df['log10_r'],df['log10_SigmaSFR'], color='black')    

        #end loop on mass bins    
        
        
        if(i_z==0):
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.5, y_percentage=0.5,  color='black', xlog=0, 
                        ylog=0,label='Star Forming', fontsize=11, fontweight='normal') 

            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.15, y_percentage=0.28,  color='black', xlog=0, 
                        ylog=0,label='Passive', fontsize=11, fontweight='normal') 
        else:
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.59, y_percentage=0.62,  color='black', xlog=0, 
                        ylog=0,label='Star Forming', fontsize=11, fontweight='normal') 

            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.28, y_percentage=0.3,  color='black', xlog=0, 
                        ylog=0,label='Passive', fontsize=11, fontweight='normal') 

        if(i_z==1):
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.2, y_percentage=0.9,  color='black', xlog=0, 
                        ylog=0,label='$10.0<\log_{10}(M_*[\mathrm{M}_{\odot}])<11.0$', fontsize=11, fontweight='normal')        
            plot_label (subplot,'line',xlim,ylim,x_percentage=0.1,y_percentage=0.93,
                        color='orange',x2_percentage=0.18,xlog=0,ylog=0,linestyle='-',linewidth=2)
            plot_label (subplot,'line',xlim,ylim,x_percentage=0.1,y_percentage=0.91,
                        color='orange',x2_percentage=0.18,xlog=0,ylog=0,linestyle='--',linewidth=2)

            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.2, y_percentage=0.82, color='black', xlog=0, 
                        ylog=0,label='$11.0<\log_{10}(M_*[\mathrm{M}_{\odot}])<12.0$', fontsize=11, fontweight='normal')        
            plot_label (subplot,'line',xlim,ylim,x_percentage=0.1,y_percentage=0.85,
                        color='red',x2_percentage=0.18,xlog=0,ylog=0,linestyle='-',linewidth=2)         
            plot_label (subplot,'line',xlim,ylim,x_percentage=0.1,y_percentage=0.83,
                        color='red',x2_percentage=0.18,xlog=0,ylog=0,linestyle='--',linewidth=2)

            
        if(i_z==0):
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.1, y_percentage=0.1,  color='black', xlog=0, 
                        ylog=0,label='z=0', fontsize=13, fontweight='normal')
        else:
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.1, y_percentage=0.1,  color='black', xlog=0, 
                        ylog=0,label='z=2', fontsize=13, fontweight='normal')
    #end loop on properties to plot
     
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYF19_gradients_mean_evo.pdf')
    plt.close()
    
    return 
#end gradients_mean_evo
'''

    
def gradients_insideout_quenching_combined(ThisRedshiftList):
    
    
    '''fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
    subplot=plt.subplot()

    xlim=[9.0, 12.0]    
    ylim=[-13.0, -8.0]         
    subplot.set_ylim(ylim),subplot.set_xlim(xlim)

    xlab='Mass'
    ylab='SFR'  
    subplot.set_xlabel(xlab,fontsize=14), subplot.set_ylabel(ylab,fontsize=14)       

    subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
    subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
    subplot.xaxis.set_major_locator(MultipleLocator(0.5))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.1)) 

    i_z=0
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)                 
    G0_MR = G_MR[sel] 
    G0_MR = G0_MR[(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)>8.5) & (G0_MR['Sfr']>0) & (G0_MR['Type']==0)]
       
    Mass = np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)
    SFR = np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h))
    
    subplot.scatter(Mass[:10000], SFR[:10000], s=2)
    
    bin=[0.2,0.1]
    Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]
    Ngals=len(Mass)/10.    
    H, xedges, yedges = np.histogram2d(Mass, SFR, bins=Nbins)            
    extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
    plt.subplots_adjust(bottom=0.15, left=0.15)        
    mylevels = np.linspace(1., Nbins[0], Nbins[0])*Ngals/(Nbins[0]**2/3.)        
    H = zoom(H, 20)        
    cont=plt.contourf(H.transpose()[::], origin='lower', cmap='Greys_r', levels=mylevels, extent=extent)    
                   
    
    
    plt.savefig('./fig/HYJ18_main_sequence.pdf')
    plt.close()'''
    
    plot_color=['orange', 'red'] 
       
    fig = plt.figure(figsize=(one_two_size_small[0],one_two_size_small[1]))
    grid = gridspec.GridSpec(1, 2)
    grid.update(wspace=0.0, hspace=0.0)
    
    for i_z in range (0,len(ThisRedshiftList)):   
    
        subplot=plt.subplot(grid[i_z])

        xlim=[-0.4, 1.5]
        ylim=[-4, 2.]    
        #ylim=[-14, -8.5]         
        subplot.set_ylim(ylim),subplot.set_xlim(xlim)

        xlab='$\log_{10}(r \mathrm{[kpc]})$'        
        ylab='$\log_{10}(\Sigma_{SFR}[M_{\odot} \mathrm{yr}^{-1}\mathrm{Kpc^{-2}}])$'         
        if(i_z==1):
            ylab=''
            plt.tick_params(axis='y', which='both', left=True, labelleft=False)
            
        subplot.set_xlabel(xlab,fontsize=14), subplot.set_ylabel(ylab,fontsize=14)       

        subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
        subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
        subplot.xaxis.set_major_locator(MultipleLocator(0.5))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.1)) 

    
        MassBins = 2
        mass_low = [10.,11.]
        mass_high = [11.,12.]

        for quenched_state in range(0,2):
            for k_mass in range(0,len(mass_low)):           

                (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)                 
                G0_MR_unsel = G_MR[sel] 
                G0_MR_unsel = G0_MR_unsel[(G0_MR_unsel['Sfr']>0.) & (G0_MR_unsel['Type']==0)]  
                
                if(i_z==0):
                    SFR_cut = -10.5
                else:
                    SFR_cut = -10.0
                    
                if(quenched_state==0):
                    sel = np.log10(G0_MR_unsel['Sfr']/(G0_MR_unsel['StellarMass']*1.e10/Hubble_h))>SFR_cut
                else:
                    sel = np.log10(G0_MR_unsel['Sfr']/(G0_MR_unsel['StellarMass']*1.e10/Hubble_h))<SFR_cut
                        
                G0_MR_unsel = G0_MR_unsel[sel]    
                   
                #mass bins
                G0_MR=G0_MR_unsel[(np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)>mass_low[k_mass]) &
                                  (np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)<mass_high[k_mass])]

                NGals=len(G0_MR)            
                x_variable=np.zeros(int(RNUM*NGals),dtype=np.float32)             
                y_variable=np.zeros(int(RNUM*NGals),dtype=np.float32)   
                new_x_var = np.arange(xlim[0],xlim[1]+1.0,0.2)


                #loop on radial bins
                for jj in range(0,RNUM):  
                    x_variable[NGals*jj:NGals*(jj+1)]=np.log10(RingRadius[jj])  
                    SFR_this_ring=G0_MR['SfrRings'][:,jj]
                    #1e6 -> from kpc^2 to pc^2
                    if(jj==0):
                        y_variable[NGals*jj:NGals*(jj+1)]=SFR_this_ring/(3.14*RingRadius[0]**2) 
                    else:
                        y_variable[NGals*jj:NGals*(jj+1)]=SFR_this_ring/(3.14*(RingRadius[jj]**2-RingRadius[jj-1]**2))                            
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
                if(quenched_state==0):
                    subplot.plot(x_binned, median, color=plot_color[k_mass], linewidth=2, linestyle='--')
                else:
                    subplot.plot(x_binned, median, color=plot_color[k_mass], linewidth=2, linestyle='-')
                    
                #WRITE OUTPUT      
                if(write_to_file==1):
                    df = pd.DataFrame({'log10_r':x_binned, 'log10_SigmaSFR':median})
                    file = Datadir+file_to_write+'Gradients_insideout_quenching'+str(f'_M{mass_low[k_mass]:0.2f}') + \
                           str(f'_{mass_high[k_mass]:0.2f}') + str(f'_Qstate{quenched_state:d}') + \
                           str(f'_z{ThisRedshiftList[i_z]:0.2f}')+'.csv'
                    df.to_csv(file,index=False)
                    #df = pd.read_csv(file)
                    #subplot.plot(df['log10_r'],df['log10_SigmaSFR'], color='black')    
                    
                #end loop on mass bins    
        if(i_z==0):
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.5, y_percentage=0.5,  color='black', xlog=0, 
                        ylog=0,label='Star Forming', fontsize=11, fontweight='normal') 

            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.15, y_percentage=0.28,  color='black', xlog=0, 
                        ylog=0,label='Passive', fontsize=11, fontweight='normal') 
        else:
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.59, y_percentage=0.62,  color='black', xlog=0, 
                        ylog=0,label='Star Forming', fontsize=11, fontweight='normal') 

            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.28, y_percentage=0.3,  color='black', xlog=0, 
                        ylog=0,label='Passive', fontsize=11, fontweight='normal') 

        '''plot_label (subplot, 'label', xlim, ylim, x_percentage=0.57, y_percentage=0.73,  color='black', xlog=0, 
                    ylog=0,label='$SSFR[\mathrm{yr}^{-1}]>10^{-10}$', fontsize=11, fontweight='normal') 
        plot_label (subplot,'line',xlim,ylim,x_percentage=0.48,y_percentage=0.76,
                    color='red',x2_percentage=0.55,xlog=0,ylog=0,linestyle='--',linewidth=2)
        plot_label (subplot,'line',xlim,ylim,x_percentage=0.48,y_percentage=0.74,
                    color='orange',x2_percentage=0.55,xlog=0,ylog=0,linestyle='--',linewidth=2)

        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.57, y_percentage=0.65,  color='black', xlog=0, 
                    ylog=0,label='$SSFR[\mathrm{yr}^{-1}]<10^{-11}$', fontsize=11, fontweight='normal')   
        plot_label (subplot,'line',xlim,ylim,x_percentage=0.48,y_percentage=0.68,
                    color='red',x2_percentage=0.55,xlog=0,ylog=0,linestyle='-',linewidth=2)
        plot_label (subplot,'line',xlim,ylim,x_percentage=0.48,y_percentage=0.66,
                    color='orange',x2_percentage=0.55,xlog=0,ylog=0,linestyle='-',linewidth=2)'''

        if(i_z==1):
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.2, y_percentage=0.9,  color='black', xlog=0, 
                        ylog=0,label='$10.0<\log_{10}(M_*[\mathrm{M}_{\odot}])<11.0$', fontsize=11, fontweight='normal')        
            plot_label (subplot,'line',xlim,ylim,x_percentage=0.1,y_percentage=0.93,
                        color='orange',x2_percentage=0.18,xlog=0,ylog=0,linestyle='-',linewidth=2)
            plot_label (subplot,'line',xlim,ylim,x_percentage=0.1,y_percentage=0.91,
                        color='orange',x2_percentage=0.18,xlog=0,ylog=0,linestyle='--',linewidth=2)

            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.2, y_percentage=0.82, color='black', xlog=0, 
                        ylog=0,label='$11.0<\log_{10}(M_*[\mathrm{M}_{\odot}])<12.0$', fontsize=11, fontweight='normal')        
            plot_label (subplot,'line',xlim,ylim,x_percentage=0.1,y_percentage=0.85,
                        color='red',x2_percentage=0.18,xlog=0,ylog=0,linestyle='-',linewidth=2)         
            plot_label (subplot,'line',xlim,ylim,x_percentage=0.1,y_percentage=0.83,
                        color='red',x2_percentage=0.18,xlog=0,ylog=0,linestyle='--',linewidth=2)

            
        if(i_z==0):
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.1, y_percentage=0.1,  color='black', xlog=0, 
                        ylog=0,label='z=0', fontsize=13, fontweight='normal')
        else:
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.1, y_percentage=0.1,  color='black', xlog=0, 
                        ylog=0,label='z=2', fontsize=13, fontweight='normal')
    #end loop on properties to plot
     
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYF19_gradients_insideout_quenching_combined.pdf')
    plt.close()
    
    return 
#end gradients_insideout_quenching_combined

def gradients_insideout_quenching_SSFR(ThisRedshiftList):
       
    plot_color=['blue', 'orange', 'red'] 
       
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
   
   
    MassBins = 3
    mass_low = [9., 10.,11.]
    mass_high = [10., 11.,12.]
    
    for quenched_state in range(0,2):
        for k_mass in range(0,len(mass_low)):           

            subplot=plt.subplot()

            xlim=[-0.5, 1.5]
            #ylim=[-3, 2.]    
            ylim=[-14, -8.5]         
            subplot.set_ylim(ylim),subplot.set_xlim(xlim)

            xlab='$\log_{10}(r \mathrm{[kpc]})$'
            ylab='$\log_{10}(\Sigma_{SSFR}[\mathrm{yr}^{-1}\mathrm{Kpc^{-2}}])$'  
            subplot.set_xlabel(xlab,fontsize=14), subplot.set_ylabel(ylab,fontsize=14)       

            subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
            subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
            subplot.xaxis.set_major_locator(MultipleLocator(0.5))    
            subplot.xaxis.set_minor_locator(MultipleLocator(0.1)) 

            i_z=0
            (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)                 
            G0_MR_unsel = G_MR[sel] 
            G0_MR_unsel = G0_MR_unsel[(G0_MR_unsel['Sfr']>0.) & (G0_MR_unsel['Type']==0)]  
            if(quenched_state==0):
                  G0_MR_unsel=G0_MR_unsel[np.log10(G0_MR_unsel['Sfr']/(G0_MR_unsel['StellarMass']*1.e10/Hubble_h))>-10.]
            else:
                G0_MR_unsel=G0_MR_unsel[np.log10(G0_MR_unsel['Sfr']/(G0_MR_unsel['StellarMass']*1.e10/Hubble_h))<-11.]

            #mass bins
            G0_MR=G0_MR_unsel[(np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)>mass_low[k_mass]) &
                              (np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)<mass_high[k_mass])]

            NGals=len(G0_MR)            
            x_variable=np.zeros(int(RNUM*NGals),dtype=np.float32)             
            y_variable=np.zeros(int(RNUM*NGals),dtype=np.float32)   
            new_x_var = np.arange(xlim[0],xlim[1]+1.0,0.2)


            #loop on radial bins
            for jj in range(0,RNUM):  
                x_variable[NGals*jj:NGals*(jj+1)]=np.log10(RingRadius[jj])  
                SFR_this_ring=(G0_MR['SfrRings'][:,jj]/(G0_MR['DiskMassRings'][:,jj]*1e10/Hubble_h+G0_MR['BulgeMassRings'][:,jj]*1e10/Hubble_h)) 
                #1e6 -> from kpc^2 to pc^2
                if(jj==0):
                    y_variable[NGals*jj:NGals*(jj+1)]=SFR_this_ring/(3.14*RingRadius[0]**2) 
                else:
                    y_variable[NGals*jj:NGals*(jj+1)]=SFR_this_ring/(3.14*(RingRadius[jj]**2-RingRadius[jj-1]**2))                            
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
            if(quenched_state==0):
                subplot.plot(x_binned, median, color=plot_color[k_mass], linewidth=2, linestyle='--')
            else:
                subplot.plot(x_binned, median, color=plot_color[k_mass], linewidth=2, linestyle='-')
            #end loop on mass bins    

        
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.35, y_percentage=0.88, color='black', xlog=0, 
                ylog=0,label='$9.0<\log_{10}(M_*[\mathrm{M}_{\odot}])<10.0$', fontsize=11, fontweight='normal')        
    plot_label (subplot,'line',xlim,ylim,x_percentage=0.27,y_percentage=0.9,
                color='blue',x2_percentage=0.33,xlog=0,ylog=0,linestyle='-',linewidth=2)
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.35, y_percentage=0.78,  color='black', xlog=0, 
                ylog=0,label='$10.0<\log_{10}(M_*[\mathrm{M}_{\odot}])<11.0$', fontsize=11, fontweight='normal')        
    plot_label (subplot,'line',xlim,ylim,x_percentage=0.27,y_percentage=0.8,
                color='orange',x2_percentage=0.33,xlog=0,ylog=0,linestyle='-',linewidth=2)
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.35, y_percentage=0.68, color='black', xlog=0, 
                ylog=0,label='$11.0<\log_{10}(M_*[\mathrm{M}_{\odot}])<12.0$', fontsize=11, fontweight='normal')        
    plot_label (subplot,'line',xlim,ylim,x_percentage=0.27,y_percentage=0.7,
                color='red',x2_percentage=0.33,xlog=0,ylog=0,linestyle='-',linewidth=2)
    
    #end loop on properties to plot
     
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYJ18_gradients_insideout_quenching_SSFR.pdf')
    plt.close()
    
    return 
#end gradients_insideout_quenching_SSFR

def gradients_insideout_quenching_SFR(ThisRedshiftList):
       
    plot_color=['orange', 'red'] 
       
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
   
   
    MassBins = 2
    mass_low = [10.,11.]
    mass_high = [11.,12.]
    
    for quenched_state in range(0,2):
        for k_mass in range(0,len(mass_low)):           

            subplot=plt.subplot()

            xlim=[-0.5, 1.5]
            ylim=[-4, 1.]    
            #ylim=[-14, -8.5]         
            subplot.set_ylim(ylim),subplot.set_xlim(xlim)

            xlab='$\log_{10}(r \mathrm{[kpc]})$'
            ylab='$\log_{10}(\Sigma_{SFR}[M_{\odot} \mathrm{yr}^{-1}\mathrm{Kpc^{-2}}])$'  
            subplot.set_xlabel(xlab,fontsize=14), subplot.set_ylabel(ylab,fontsize=14)       

            subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
            subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
            subplot.xaxis.set_major_locator(MultipleLocator(0.5))    
            subplot.xaxis.set_minor_locator(MultipleLocator(0.1)) 

            i_z=0
            (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)                 
            G0_MR_unsel = G_MR[sel] 
            G0_MR_unsel = G0_MR_unsel[(G0_MR_unsel['Sfr']>0.) & (G0_MR_unsel['Type']==0)]  
            if(quenched_state==0):
                  G0_MR_unsel=G0_MR_unsel[np.log10(G0_MR_unsel['Sfr']/(G0_MR_unsel['StellarMass']*1.e10/Hubble_h))>-10.]
            else:
                G0_MR_unsel=G0_MR_unsel[np.log10(G0_MR_unsel['Sfr']/(G0_MR_unsel['StellarMass']*1.e10/Hubble_h))<-11.]

            #mass bins
            G0_MR=G0_MR_unsel[(np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)>mass_low[k_mass]) &
                              (np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)<mass_high[k_mass])]

            NGals=len(G0_MR)            
            x_variable=np.zeros(int(RNUM*NGals),dtype=np.float32)             
            y_variable=np.zeros(int(RNUM*NGals),dtype=np.float32)   
            new_x_var = np.arange(xlim[0],xlim[1]+1.0,0.2)


            #loop on radial bins
            for jj in range(0,RNUM):  
                x_variable[NGals*jj:NGals*(jj+1)]=np.log10(RingRadius[jj])  
                SFR_this_ring=G0_MR['SfrRings'][:,jj]
                #1e6 -> from kpc^2 to pc^2
                if(jj==0):
                    y_variable[NGals*jj:NGals*(jj+1)]=SFR_this_ring/(3.14*RingRadius[0]**2) 
                else:
                    y_variable[NGals*jj:NGals*(jj+1)]=SFR_this_ring/(3.14*(RingRadius[jj]**2-RingRadius[jj-1]**2))                            
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
            if(quenched_state==0):
                subplot.plot(x_binned, median, color=plot_color[k_mass], linewidth=2, linestyle='--')
            else:
                subplot.plot(x_binned, median, color=plot_color[k_mass], linewidth=2, linestyle='-')
            #end loop on mass bins    
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.57, y_percentage=0.62,  color='black', xlog=0, 
                ylog=0,label='Star Forming', fontsize=11, fontweight='normal') 
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.23, y_percentage=0.3,  color='black', xlog=0, 
                ylog=0,label='Passive', fontsize=11, fontweight='normal') 
    
    '''plot_label (subplot, 'label', xlim, ylim, x_percentage=0.57, y_percentage=0.73,  color='black', xlog=0, 
                ylog=0,label='$SSFR[\mathrm{yr}^{-1}]>10^{-10}$', fontsize=11, fontweight='normal') 
    plot_label (subplot,'line',xlim,ylim,x_percentage=0.48,y_percentage=0.76,
                color='red',x2_percentage=0.55,xlog=0,ylog=0,linestyle='--',linewidth=2)
    plot_label (subplot,'line',xlim,ylim,x_percentage=0.48,y_percentage=0.74,
                color='orange',x2_percentage=0.55,xlog=0,ylog=0,linestyle='--',linewidth=2)
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.57, y_percentage=0.65,  color='black', xlog=0, 
                ylog=0,label='$SSFR[\mathrm{yr}^{-1}]<10^{-11}$', fontsize=11, fontweight='normal')   
    plot_label (subplot,'line',xlim,ylim,x_percentage=0.48,y_percentage=0.68,
                color='red',x2_percentage=0.55,xlog=0,ylog=0,linestyle='-',linewidth=2)
    plot_label (subplot,'line',xlim,ylim,x_percentage=0.48,y_percentage=0.66,
                color='orange',x2_percentage=0.55,xlog=0,ylog=0,linestyle='-',linewidth=2)'''
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.38, y_percentage=0.9,  color='black', xlog=0, 
                ylog=0,label='$10.0<\log_{10}(M_*[\mathrm{M}_{\odot}])<11.0$', fontsize=11, fontweight='normal')        
    plot_label (subplot,'line',xlim,ylim,x_percentage=0.3,y_percentage=0.93,
                color='orange',x2_percentage=0.36,xlog=0,ylog=0,linestyle='-',linewidth=2)
    plot_label (subplot,'line',xlim,ylim,x_percentage=0.3,y_percentage=0.91,
                color='orange',x2_percentage=0.36,xlog=0,ylog=0,linestyle='--',linewidth=2)
    
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.38, y_percentage=0.82, color='black', xlog=0, 
                ylog=0,label='$11.0<\log_{10}(M_*[\mathrm{M}_{\odot}])<12.0$', fontsize=11, fontweight='normal')        
    plot_label (subplot,'line',xlim,ylim,x_percentage=0.3,y_percentage=0.85,
                color='red',x2_percentage=0.36,xlog=0,ylog=0,linestyle='-',linewidth=2)         
    plot_label (subplot,'line',xlim,ylim,x_percentage=0.3,y_percentage=0.83,
                color='red',x2_percentage=0.36,xlog=0,ylog=0,linestyle='--',linewidth=2)
    
    #end loop on properties to plot
     
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYJ18_gradients_insideout_quenching_SFR.pdf')
    plt.close()
    
    return 
#end gradients_insideout_quenching_SFR

def gradients_insideout_quenching(ThisRedshiftList):
       
    plot_color=['red', 'blue'] 
    plot_color_obs=['darkorange', 'cornflowerblue'] 
    linestyles=['--','-']
    
    fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
   
   
    MassBins = 2
    mass_low = [9.5, 10.85]
    mass_high = [10.3, 11.7]
    
    for k_mass in range(0,len(mass_low)):           

        subplot=plt.subplot()

        xlim=[-0.5, 1.4]
        ylim=[7.0, 10.5]         
        subplot.set_ylim(ylim),subplot.set_xlim(xlim)

        xlab='$r \mathrm{(kpc)}$'
        ylab='$\log_{10}(\Sigma_*[M_{\odot} \mathrm{Kpc^{-2}}])$'  
        subplot.set_xlabel(xlab,fontsize=14), subplot.set_ylabel(ylab,fontsize=14)       

        subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
        subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
        subplot.xaxis.set_major_locator(MultipleLocator(0.5))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.1)) 

        for i_z in range(0, len(ThisRedshiftList)):
            (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)                 
            G0_MR_unsel = G_MR[sel] 
          
            #main sequence selection
            if(i_z==0):  
                G0_MR_unsel=G0_MR_unsel[np.log10(G0_MR_unsel['Sfr']/(G0_MR_unsel['StellarMass']*1.e10/Hubble_h))<-10.5]
            else:    
                G0_MR_unsel=G0_MR_unsel[np.log10(G0_MR_unsel['Sfr']/(G0_MR_unsel['StellarMass']*1.e10/Hubble_h))>-9.]
                
            #mass bins
            G0_MR=G0_MR_unsel[(np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)>mass_low[k_mass]) &
                              (np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)<mass_high[k_mass])]

            NGals=len(G0_MR)            
            x_variable=np.zeros(int(RNUM*NGals),dtype=np.float32)             
            y_variable=np.zeros(int(RNUM*NGals),dtype=np.float32)   
            new_x_var = np.arange(xlim[0],xlim[1]+1.0,0.1)
            

            #loop on radial bins
            for jj in range(0,RNUM):  
                x_variable[NGals*jj:NGals*(jj+1)]=np.log10(RingRadius[jj])  
                Mass_this_ring=(G0_MR['DiskMassRings'][:,jj]*1e10/Hubble_h)+(G0_MR['BulgeMassRings'][:,jj]*1e10/Hubble_h) 
                #1e6 -> from kpc^2 to pc^2
                if(jj==0):
                    y_variable[NGals*jj:NGals*(jj+1)]=Mass_this_ring/(3.14*RingRadius[0]**2) 
                else:
                    y_variable[NGals*jj:NGals*(jj+1)]=Mass_this_ring/(3.14*(RingRadius[jj]**2-RingRadius[jj-1]**2))                            
            #endfor RNUM      
            N_random_gals = 20000
            if(N_random_gals>NGals):
                N_random_gals=NGals
            random.seed(a=2)    
            random_list = random.sample(range(0, NGals), N_random_gals)
            
            interpol_x_variable=np.zeros(int(len(new_x_var)*N_random_gals),dtype=np.float32)      
            interpol_y_variable=np.zeros(int(len(new_x_var)*N_random_gals),dtype=np.float32)   
           
            i_index = 0
            for ii in random_list:                
                slice_ii = [x*NGals+ii for x in range(0,12)]                
                xx = x_variable[slice_ii]
                yy = np.log10(y_variable[slice_ii])                 
                #ignore galaxies without halflightradius or with nan on the y_variable
                sel = (~np.isnan(xx)) & (~np.isinf(xx)) & (~np.isnan(yy)) & (~np.isinf(yy))                 
                if(len(xx[sel])>3):                   
                    f = interpolate.UnivariateSpline(xx[sel], yy[sel], s=0)                      
                    interpol_y_variable[i_index*len(new_x_var):(i_index+1)*len(new_x_var)] = f(new_x_var)
                    interpol_x_variable[i_index*len(new_x_var):(i_index+1)*len(new_x_var)] = new_x_var
                i_index += 1
          
            (x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles_fixed_xx(interpol_x_variable, 
                                                                                      interpol_y_variable, non_zero=0)       
            subplot.plot(x_binned, median, color=plot_color[i_z], linewidth=2, linestyle=linestyles[1])
          
            #Observations
            if(k_mass==0):  
                df = pd.read_csv(Datadir+'tachella2015_highz_massbin1.csv')
                error_down = np.log10(df['y']/(df['y']+df['err_down']))
                error_up = np.log10((df['y']+df['err_down'])/df['y'])
                #subplot.errorbar(np.log10(df['x']), np.log10(df['y']), yerr=[error_down,error_up], fmt='o', markersize=5, 
                #                 color=plot_color[0], ecolor=plot_color[0])
                subplot.plot(np.log10(df['x']),np.log10(df['y']),color=plot_color_obs[0],linewidth=2, linestyle=linestyles[0])

                df = pd.read_csv(Datadir+'tachella2015_lowz_massbin1.csv')            
                #subplot.fill_between(np.log10(df['x']),np.log10(df['y']+df['err_down']),np.log10(df['y']+df['err_up']),
                #             facecolor='grey', interpolate=True, alpha=0.3)  
                subplot.plot(np.log10(df['x']),np.log10(df['y']),color=plot_color_obs[1],linewidth=2, linestyle=linestyles[0])
                
            if(k_mass==1):  
                df = pd.read_csv(Datadir+'tachella2015_highz_massbin2.csv')
                error_down = np.log10(df['y']/(df['y']+df['err_down']))
                error_up = np.log10((df['y']+df['err_down'])/df['y'])
                #subplot.errorbar(np.log10(df['x']), np.log10(df['y']), yerr=[error_down,error_up], fmt='o', markersize=5, 
                #                 color=plot_color[0], ecolor=plot_color[0])  
                subplot.plot(np.log10(df['x']),np.log10(df['y']),color=plot_color_obs[0],linewidth=2, linestyle=linestyles[0])

                df = pd.read_csv(Datadir+'tachella2015_lowz_massbin2.csv')            
                #subplot.fill_between(np.log10(df['x']),np.log10(df['y']+df['err_down']),np.log10(df['y']+df['err_up']),
                #             facecolor='grey', interpolate=True, alpha=0.3)  
                subplot.plot(np.log10(df['x']),np.log10(df['y']),color=plot_color_obs[1],linewidth=2, linestyle=linestyles[0])
                
        #end loop on mass bins    

    
        
    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.7, y_percentage=0.78, 
                color='black', xlog=0, ylog=0,label='Tachella 2015', fontsize=11, fontweight='normal')        
    plot_label (subplot,'line',xlim,ylim,x_percentage=0.62,y_percentage=0.81,
                color=plot_color_obs[0],x2_percentage=0.68,xlog=0,ylog=0,linestyle='--',linewidth=2)
    plot_label (subplot,'line',xlim,ylim,x_percentage=0.62,y_percentage=0.79,
                color=plot_color_obs[1],x2_percentage=0.68,xlog=0,ylog=0,linestyle='--',linewidth=2)

    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.7, y_percentage=0.88, 
                color='black', xlog=0, ylog=0,label='This Work', fontsize=11, fontweight='normal')        
    plot_label (subplot,'line',xlim,ylim,x_percentage=0.62,y_percentage=0.91,
                color='red',x2_percentage=0.68,xlog=0,ylog=0,linestyle='-',linewidth=2)
    plot_label (subplot,'line',xlim,ylim,x_percentage=0.62,y_percentage=0.89,
                color='blue',x2_percentage=0.68,xlog=0,ylog=0,linestyle='-',linewidth=2)
            
    #end loop on properties to plot
       
    #SSFR distribution at z=2 for selection    
    #subplot=plt.subplot(grid[2])         
    #xlim=[0.0, 12.0]
    #ylim=[-13.0, -8.0]         
    
    #G0_MR=G0_MR_unsel[np.log10(G0_MR_unsel['Sfr']/(G0_MR_unsel['StellarMass']*1.e10/Hubble_h))>-9.5]
    #subplot.scatter(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h), 
    #                    np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h)), marker='o',s=1)
        
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYJ18_gradients_insideout_quenching.pdf')
    plt.close()
    
    return 
#end gradients_insideout_quenching


def all_gradients(ThisRedshiftList):
    i_z=0   
    
    labels = ['SigmaStar','StellarMetallicity','SigmaGas','GasMetallicity']
    plot_color=['purple', 'blue', 'green', 'darkorange', 'red']    
    fig = plt.figure(figsize=(two_two_size_large[0],two_two_size_large[1]))
    grid = gridspec.GridSpec(2, 2)
    grid.update(wspace=0.4, hspace=0.3)    
    plt.subplots_adjust(left=0.12, right=0.97, top=0.97, bottom=0.09)             
    
    MassBins = 5
    mass_low = [9.0, 9.5, 10.0, 10.5, 11.0]
    mass_high = [9.5, 10.0, 10.5, 11.0, 11.5]
    
    #mass_low = [9.0]
    #mass_high = [9.5]
  
      
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, i_z, FullSnapshotList_MR)                 
    G0_MR_unsel = G_MR[sel] 
    G0_MR_unsel = G0_MR_unsel[np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)>9.0]
    
    #SFH structure needed for age gradients
    fa = open(DirName_MR+"SFH_Bins","rb")                
    nbins =  np.fromfile(fa,np.int32,1)
    template = np.dtype([('SnapNum',np.int32,1), ('Bin',np.int32,1), ('Lookbacktime',np.float64,1),                           
                         ('dt',np.float64,1),('nbins',np.int32,1)])
    SFH = np.fromfile(fa,template,int(nbins))    
    fa.close()                
    SFH=SFH[SFH['SnapNum']==G0_MR_unsel['SnapNum'][0]]             
   
    Nprops=4
    #loop on quantity to plot
    for plot_prop in range (0, Nprops):
                        
        subplot=plt.subplot(grid[plot_prop])
                 
        if(plot_prop==0):
            #***************************    
            #* stellar surface density *
            #*************************** 
            ylim=[1.0, 5.0]     
            ylab='$\log_{10}(\Sigma_*/(M_{\odot}\mathrm{pc^{-2}}))$'   
            subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
            subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
         
        if(plot_prop==1):                  
            #***************************    
            #*   stellar metallicity   *
            #*************************** 
            ylim=[-1., 0.4]     
            ylab='$\mathrm{log_{10}}(Z_*/Z_{\odot})$'
            subplot.yaxis.set_major_locator(MultipleLocator(0.5))    
            subplot.yaxis.set_minor_locator(MultipleLocator(0.1))
        
        if(plot_prop==2):
            #***************************    
            #*         ColdGas         *
            #***************************                  
            ylim=[-0.5, 2.5]       
            ylab='$\log_{10}(\Sigma_{\mathrm{cold}}/(M_{\odot} \mathrm{pc^{-2}}))$'      
            subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
            subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
                     
        if(plot_prop==3):                  
            #***************************    
            #*   Cold Gas metallicity   *
            #*************************** 
            ylim=[7.5, 9.4]     
            ylab='$12 + \log_{10}(\mathrm{O/H})_{\mathrm{cold}}$'
            subplot.yaxis.set_major_locator(MultipleLocator(0.5))    
            subplot.yaxis.set_minor_locator(MultipleLocator(0.1))
             
        xlim=[0.0, 3.0]
        subplot.set_ylim(ylim),subplot.set_xlim(xlim)
        xlab='$r/r_{e}$'
        subplot.set_xlabel(xlab,fontsize=14), subplot.set_ylabel(ylab,fontsize=14)                      
        subplot.xaxis.set_major_locator(MultipleLocator(0.5))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.1)) 
        
        #For comparison with MaNGA metallicity gradients select only galaxies with non-zero metallicity    
        if(plot_prop==2):  
            if(opt_detailed_enrichment==1):   
                G0_MR_unsel=G0_MR_unsel[((G0_MR_unsel['MetalsColdGas'][:,0] + 
                              G0_MR_unsel['MetalsColdGas'][:,1] + 
                              G0_MR_unsel['MetalsColdGas'][:,2])>.0)] 
            else:
                G0_MR_unsel=G0_MR_unsel[(G0_MR_unsel['MetalsColdGas']>.0)]             
        
        #if(plot_prop==0 or plot_prop==1):  
        #    G0_MR_unsel=G0_MR_unsel[(G0_MR_unsel['DiskMass']/G0_MR_unsel['StellarMass']>0.9)]
                
        '''if(plot_prop==0):  
            G0_MR=G0_MR_unsel[(np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)>mass_low[0]) &
                              (np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)<mass_high[0])]
            
            area=np.zeros(int(RNUM),dtype=np.float32)   
            for jj in range(0,12):
                #1e6 -> from kpc^2 to pc^2
                if(jj==0):
                    area[jj]=(3.14*(RingRadius[0]**2*1e6)) 
                else:
                    area[jj]=(3.14*(RingRadius[jj]**2-RingRadius[jj-1]**2)*1e6)      
                        
            #for i_gal in range(0,len(G0_MR)):
            for i_gal in range(0,500):
                x_variable=RingRadius/(G0_MR['StellarHalfLightRadius'][i_gal]*1000./Hubble_h)
                Mass = (G0_MR['DiskMassRings'][i_gal,:]*1e10/Hubble_h)+(G0_MR['BulgeMassRings'][i_gal,:]*1e10/Hubble_h)   
                y_variable=Mass/area
                subplot.plot(x_variable,np.log10(y_variable),color='blue')'''
                
        #loop on mass bins    
        for k_mass in range(0,len(mass_low)):             
            G0_MR=G0_MR_unsel[(np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)>mass_low[k_mass]) &
                              (np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)<mass_high[k_mass])]
                                
            NGals=len(G0_MR)            
            x_variable=np.zeros(int(RNUM*NGals),dtype=np.float32)             
            y_variable=np.zeros(int(RNUM*NGals),dtype=np.float32)   
            new_x_var = np.arange(xlim[0],xlim[1],0.1)
            interpol_x_variable=np.zeros(int(len(new_x_var)*NGals),dtype=np.float32)      
            interpol_y_variable=np.zeros(int(len(new_x_var)*NGals),dtype=np.float32)       
            
            #loop on radial bins
            for jj in range(0,RNUM):  
                
                if(plot_prop==0):  
                    x_variable[NGals*jj:NGals*(jj+1)]=RingRadius[jj]/(G0_MR['StellarHalfMassRadius']*1000./Hubble_h)                     
                elif(plot_prop==1):      
                    x_variable[NGals*jj:NGals*(jj+1)]=RingRadius[jj]/(G0_MR['StellarHalfMassRadius']*1000./Hubble_h)
                    #x_variable[NGals*jj:NGals*(jj+1)]=RingRadius[jj]/(G0_MR['StellarDiskRadius']*1000./Hubble_h)
                else:
                    x_variable[NGals*jj:NGals*(jj+1)]=RingRadius[jj]/(G0_MR['GasDiskRadius']*1000./Hubble_h)
                    
                x_variable[NGals*jj:NGals*(jj+1)]=RingRadius[jj]/(G0_MR['StellarHalfMassRadius']*1000./Hubble_h) 
               
                #***************************    
                #* Stellar surface density *
                #***************************
                if(plot_prop==0):                     
                    Mass_this_ring=(G0_MR['DiskMassRings'][:,jj]*1e10/Hubble_h)+(G0_MR['BulgeMassRings'][:,jj]*1e10/Hubble_h) 
                    #Mass_this_ring=(G0_MR['DiskMassRings'][:,jj]*1e10/Hubble_h)
                    #1e6 -> from kpc^2 to pc^2
                    if(jj==0):
                        y_variable[NGals*jj:NGals*(jj+1)]=Mass_this_ring/(3.14*(RingRadius[0]**2*1e6)) 
                    else:
                        y_variable[NGals*jj:NGals*(jj+1)]=Mass_this_ring/(3.14*(RingRadius[jj]**2-RingRadius[jj-1]**2)*1e6)          
                        
                   
                #***************************    
                #*   stellar metallicity   *
                #***************************   
                if(plot_prop==1): 
                    if(opt_detailed_enrichment==1):                  
                        Metals_this_ring = (G0_MR['MetalsDiskMassRings'][:,jj,0] + 
                                            G0_MR['MetalsDiskMassRings'][:,jj,1] +
                                            G0_MR['MetalsDiskMassRings'][:,jj,2])
                        Metals_this_ring += (G0_MR['MetalsBulgeMassRings'][:,jj,0] + 
                                             G0_MR['MetalsBulgeMassRings'][:,jj,1] +
                                             G0_MR['MetalsBulgeMassRings'][:,jj,2])
                    else:            
                        MetalsDiskMass_this_ring=G0_MR['MetalsDiskMassRings'][:,jj]  
                    Mass_this_Ring = G0_MR['DiskMassRings'][:,jj] + G0_MR['BulgeMassRings'][:,jj]   
                    y_variable[NGals*jj:NGals*(jj+1)]= Metals_this_ring/Mass_this_Ring/0.02    
                    
                #***************************    
                #* ColdGas surface density *
                #***************************
                if(plot_prop==2):                     
                    Mass_this_ring=(G0_MR['ColdGasRings'][:,jj]*1e10/Hubble_h)                      
                    #1e6 -> from kpc^2 to pc^2
                    if(jj==0):
                        y_variable[NGals*jj:NGals*(jj+1)]=Mass_this_ring/(3.14*RingRadius[0]**2*1e6) 
                    else:
                        y_variable[NGals*jj:NGals*(jj+1)]=Mass_this_ring/(3.14*(RingRadius[jj]**2-RingRadius[jj-1]**2)*1e6)                            
                
                   
                #***************************    
                #*     Gas metallicity     *
                #***************************   
                if(plot_prop==3): 
                    if(opt_detailed_enrichment==1):                  
                        MetalsColdGas_this_ring=(G0_MR['MetalsColdGasRings'][:,jj,0] + 
                                                  G0_MR['MetalsColdGasRings'][:,jj,1] +
                                                  G0_MR['MetalsColdGasRings'][:,jj,2])
                    else:            
                        MetalsColdGas_this_ring=G0_MR['MetalsColdGasRings'][:,jj]
                                    
                    y_variable[NGals*jj:NGals*(jj+1)] = 10**(np.log10(MetalsColdGas_this_ring / 
                                                                 G0_MR['ColdGasRings'][:,jj]/0.0134) + 8.69)                 
                    
            #endfor RNUM
            
             
            N_random_gals = 20000
            if(N_random_gals>NGals):
                N_random_gals=NGals
            random.seed(a=1)    
            random_list = random.sample(range(0, NGals), N_random_gals)
            
            for ii in random_list:
                slice_ii = [x*NGals+ii for x in range(0,RNUM)]                
                xx = x_variable[slice_ii]
                yy = np.log10(y_variable[slice_ii]) 
                             
                #ignore galaxies without halflightradius or with nan on the y_variable
                sel = (~np.isnan(xx)) & (~np.isinf(xx)) & (~np.isnan(yy)) & (~np.isinf(yy))                 
                if(len(xx[sel])>3):                   
                    f = interpolate.UnivariateSpline(xx[sel], yy[sel], s=0)          
                    interpol_y_variable[ii*len(new_x_var):(ii+1)*len(new_x_var)] = f(new_x_var)
                    interpol_x_variable[ii*len(new_x_var):(ii+1)*len(new_x_var)] = new_x_var
                        
            sel = interpol_x_variable>0.            
            (x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles_fixed_xx(interpol_x_variable[sel],
                                                                                      interpol_y_variable[sel])
        
            subplot.plot(x_binned, median, color=plot_color[k_mass], linewidth=2)
        
            #WRITE OUTPUT      
            if(write_to_file==1):
                df = pd.DataFrame({'log10_M':x_binned, 'HalfMassRadius':median, 'pc16':median-rms, 'pc84':median+rms})
                file = Datadir+file_to_write+'All_Gradients_'+labels[plot_prop] + str(f'_M{mass_low[k_mass]:0.2f}') + \
                       str(f'_{mass_high[k_mass]:0.2f}') + str(f'_z{ThisRedshiftList[i_z]:0.2f}') + '.csv'
                df.to_csv(file,index=False)
                #df = pd.read_csv(file)
                #subplot.plot(df['log10_M'],df['HalfMassRadius'], color='black')
                
        #end loop on mass bins    
                    
    
        #Observations
        if(k_mass == MassBins-1):
            if(plot_prop==0 or plot_prop==1):
                #CALIFA (stellar surface densities)
                CALIFA_mass_low =[9.1,9.6,10.1,10.6,11.2]
                CALIFA_mass_high=[9.6,10.1,10.6,10.9,11.5]
                CALIFA_mass_index=[0,1,2,3,5]  
                obs_rings=np.arange(0.05,2.8,0.1)   
           
                for k_mass in range(0,len(CALIFA_mass_low)):
                    char_mass_low="%0.1f" % CALIFA_mass_low[k_mass]
                    char_mass_high="%0.1f" % CALIFA_mass_high[k_mass]
                    char_k_mass="%d" %  CALIFA_mass_index[k_mass]
                    k_color = plot_color[k_mass]
                    
                    if(plot_prop==0):
                        file = Datadir+'/CALIFA_rosa_stellar_density_'+char_k_mass+'_'+char_mass_low+'_'+char_mass_high+'.txt'
                        obs = Table.read(file, format='ascii')      
                        sel=obs['stellar_dens']!=0.0
                        subplot.plot(obs_rings[sel],obs['stellar_dens'][sel], color=k_color, linestyle='--',linewidth=2) 
                    
                    if(plot_prop==1):                 
                        file = Datadir + '/CALIFA_rosa_metals_'+char_k_mass+'_'+char_mass_low+'_'+char_mass_high+'.txt'
                        obs = Table.read(file, format='ascii')      
                        sel=obs['metallicity']!=0.0      
                        subplot.plot(obs_rings[sel],obs['metallicity'][sel], color=k_color,linestyle='--',linewidth=2)
                        
                        '''file = Datadir + '/MANGA_gradients/goddard2016_LT_new_MAD.txt'
                        obs = Table.read(file, format='ascii')            
                        subplot.plot(obs['Radius'],obs['MW_Metallicity_M1'], linewidth=2, linestyle=':', color=plot_color[0])
                        subplot.plot(obs['Radius'],obs['MW_Metallicity_M2'], linewidth=2, linestyle=':', color=plot_color[2]) 
                        subplot.plot(obs['Radius'],obs['MW_Metallicity_M3'], linewidth=2, linestyle=':', color=plot_color[3])'''
            
            if(plot_prop==3):               
                df = pd.read_csv(Datadir+'ferrer2019/ferrer2019_dop_lowmass.csv')
                subplot.errorbar(df['radius']+0.1, df['metallicity'], xerr=df['err_x'], 
                                 yerr=[df['err_y_down'],df['err_y_up']], fmt='o', markersize=5, 
                                 color=plot_color[1], ecolor=plot_color[1],capsize=2)
                
                df = pd.read_csv(Datadir+'ferrer2019/ferrer2019_dop_intermass.csv')
                subplot.errorbar(df['radius']+0.1, df['metallicity'], xerr=df['err_x'], 
                                 yerr=[df['err_y_down'],df['err_y_up']], fmt='o', markersize=5, 
                                 color=plot_color[2], ecolor=plot_color[2],capsize=2)
                
                df = pd.read_csv(Datadir+'ferrer2019/ferrer2019_dop_highmass.csv')
                subplot.errorbar(df['radius']+0.1, df['metallicity'], xerr=df['err_x'], 
                                 yerr=[df['err_y_down'],df['err_y_up']], fmt='o', markersize=5, 
                                 color=plot_color[4], ecolor=plot_color[4],capsize=2)
            
        #labels 
        if(plot_prop==0):
            plot_label (subplot, 'label',xlim,ylim,x_percentage=0.26,y_percentage=0.9, 
                        color='black',xlog=0,ylog=0,label='$\log_{10}(M_*/M_{\odot})=$', fontsize=10,fontweight='normal') 
            
            for k_mass in range(0,MassBins):
                char_mass_low="%0.1f" % mass_low[k_mass]
                char_mass_high="%0.1f" % mass_high[k_mass]
                x_values=[0.72, 0.36, 0.655, 0.33, 0.655]
                y_values=[0.9,0.82,0.82,0.74,0.74]
        
                label='['+char_mass_low+','+char_mass_high+']'
                plot_label (subplot, 'label',xlim,ylim,x_percentage=x_values[k_mass],y_percentage=y_values[k_mass], 
                color=plot_color[k_mass],xlog=0,ylog=0,label=label,fontsize=10,fontweight='normal') 
        
        if(plot_prop==1):
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.7, y_percentage=0.85, 
                        color='black', xlog=0, ylog=0,label='CALIFA', fontsize=11, fontweight='normal')        
            plot_label (subplot,'line',xlim,ylim,x_percentage=0.57,y_percentage=0.88,
                        color='red',x2_percentage=0.67,xlog=0,ylog=0,linestyle='--',linewidth=2)
            
            '''plot_label (subplot, 'label', xlim, ylim, x_percentage=0.7, y_percentage=0.78, 
                        color='black', xlog=0, ylog=0,label='MaNGA', fontsize=11, fontweight='normal')        
            plot_label (subplot,'line',xlim,ylim,x_percentage=0.57,y_percentage=0.8,
                        color='darkorange',x2_percentage=0.67,xlog=0,ylog=0,linestyle=':',linewidth=2)'''
                 
        if(plot_prop==3):
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.55, y_percentage=0.87, 
                        color='black', xlog=0, ylog=0,label='MUSE - MAD', fontsize=11, fontweight='normal')        
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.5, y_percentage=0.9, 
                        color='red', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.1) 
            
    #end loop on properties to plot
       
        
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYF19_all_gradients.pdf')
    plt.close()
    
    return 
#end gas_gradients

def gas_gradients(ThisRedshiftList):
    ii=0   
    
    plot_color=['purple', 'blue', 'green', 'darkorange', 'red']    
    fig = plt.figure(figsize=(three_one_size_small[0],three_one_size_small[1]))
    grid = gridspec.GridSpec(3, 1)
    grid.update(wspace=0.0, hspace=0.0)
   
    MassBins = 5
    mass_low = [9.0, 9.5, 10.0, 10.5, 11.0]
    mass_high = [9.5, 10.0, 10.5, 11.0, 11.5]
  
      
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)                 
    G0_MR_unsel = G_MR[sel] 
    G0_MR_unsel = G0_MR_unsel[np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)>9.0]
    
    #SFH structure needed for age gradients
    fa = open(DirName_MR+"SFH_Bins","rb")                
    nbins =  np.fromfile(fa,np.int32,1)
    template = np.dtype([('SnapNum',np.int32,1), ('Bin',np.int32,1), ('Lookbacktime',np.float64,1),                           
                         ('dt',np.float64,1),('nbins',np.int32,1)])
    SFH = np.fromfile(fa,template,int(nbins))    
    fa.close()                
    SFH=SFH[SFH['SnapNum']==G0_MR_unsel['SnapNum'][0]]             
        
        
    i_grid=0
    Nprops=3
    #loop on quantity to plot
    for plot_prop in range (0, Nprops):
                        
        subplot=plt.subplot(grid[i_grid])
        i_grid+=1
                        
        if(plot_prop==0):
            #***************************    
            #* stellar surface density *
            #*************************** 
            ylim=[1.5, 5.0]     
            ylab='$\log_{10}(\Sigma_*[M_{\odot}/\mathrm{pc^2}])$'   
            subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
            subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
            
        if(plot_prop==1):
            #***************************    
            #*    ColdGas age    *
            #***************************                  
            ylim=[0.0, 2.5]       
            ylab='$\log_{10}(\Sigma_{\mathrm{cold}}[M_{\odot}/\mathrm{pc^2}])$'      
            subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
            subplot.yaxis.set_minor_locator(MultipleLocator(0.5)) 
                     
        if(plot_prop==2):                  
            #***************************    
            #*   Cold Gas metallicity   *
            #*************************** 
            ylim=[7.5, 9.4]     
            ylab='$12 + \log_{10}(\mathrm{O/H})_{\mathrm{cold}}$'
            subplot.yaxis.set_major_locator(MultipleLocator(0.5))    
            subplot.yaxis.set_minor_locator(MultipleLocator(0.1))
             
        xlim=[0.0, 2.0]
        subplot.set_ylim(ylim),subplot.set_xlim(xlim)
        xlab='$r/r_{e}$'
        subplot.set_xlabel(xlab,fontsize=14), subplot.set_ylabel(ylab,fontsize=14)                      
        if((plot_prop==0) or (plot_prop==1)):
                plt.tick_params(axis='x', which='both', bottom=True, labelbottom=False)   
          
        
        #For comparison with MaNGA metallicity gradients select only galaxies with non-zero metallicity    
        if(plot_prop==2):  
            if(opt_detailed_enrichment==1):   
                G0_MR_unsel=G0_MR_unsel[((G0_MR_unsel['MetalsColdGas'][:,0] + 
                              G0_MR_unsel['MetalsColdGas'][:,1] + 
                              G0_MR_unsel['MetalsColdGas'][:,2])>.0)] 
            else:
                G0_MR_unsel=G0_MR_unsel[(G0_MR_unsel['MetalsColdGas']>.0)]             
        
      
        #loop on mass bins    
        for k_mass in range(0,len(mass_low)):             
            G0_MR=G0_MR_unsel[(np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)>mass_low[k_mass]) &
                              (np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)<mass_high[k_mass])]
         
            NGals=len(G0_MR)            
            x_variable=np.zeros(int(RNUM*NGals),dtype=np.float32)             
            y_variable=np.zeros(int(RNUM*NGals),dtype=np.float32)   
            new_x_var = np.arange(xlim[0],xlim[1],0.1)
            interpol_x_variable=np.zeros(int(len(new_x_var)*NGals),dtype=np.float32)      
            interpol_y_variable=np.zeros(int(len(new_x_var)*NGals),dtype=np.float32)       
            
            #loop on radial bins
            for jj in range(0,RNUM):  
                
                if(plot_prop==0):  
                    x_variable[NGals*jj:NGals*(jj+1)]=RingRadius[jj]/(G0_MR['StellarHalfLightRadius']*1000./Hubble_h)
                else:
                    x_variable[NGals*jj:NGals*(jj+1)]=RingRadius[jj]/(G0_MR['GasDiskRadius']*1000./Hubble_h)
                #***************************    
                #* Stellar surface density *
                #***************************
                if(plot_prop==0):                     
                    Mass_this_ring=(G0_MR['DiskMassRings'][:,jj]*1e10/Hubble_h)+(G0_MR['BulgeMassRings'][:,jj]*1e10/Hubble_h) 
                    #1e6 -> from kpc^2 to pc^2
                    if(jj==0):
                        y_variable[NGals*jj:NGals*(jj+1)]=Mass_this_ring/(3.14*RingRadius[0]**2*1e6) 
                    else:
                        y_variable[NGals*jj:NGals*(jj+1)]=Mass_this_ring/(3.14*(RingRadius[jj]**2-RingRadius[jj-1]**2)*1e6)                            
                    
                #***************************    
                #* ColdGas surface density *
                #***************************
                if(plot_prop==1):                     
                    Mass_this_ring=(G0_MR['ColdGasRings'][:,jj]*1e10/Hubble_h)                      
                    #1e6 -> from kpc^2 to pc^2
                    if(jj==0):
                        y_variable[NGals*jj:NGals*(jj+1)]=Mass_this_ring/(3.14*RingRadius[0]**2*1e6) 
                    else:
                        y_variable[NGals*jj:NGals*(jj+1)]=Mass_this_ring/(3.14*(RingRadius[jj]**2-RingRadius[jj-1]**2)*1e6)                            
                
                   
                #***************************    
                #*   Gas metallicity   *
                #***************************   
                if(plot_prop==2): 
                    if(opt_detailed_enrichment==1):                  
                        MetalsColdGas_this_ring=(G0_MR['MetalsColdGasRings'][:,jj,0] + 
                                                  G0_MR['MetalsColdGasRings'][:,jj,1] +
                                                  G0_MR['MetalsColdGasRings'][:,jj,2])
                    else:            
                        MetalsColdGas_this_ring=G0_MR['MetalsColdGasRings'][:,jj]
                                    
                    y_variable[NGals*jj:NGals*(jj+1)] = 10**(np.log10(MetalsColdGas_this_ring / 
                                                                 G0_MR['ColdGasRings'][:,jj]/0.0134) + 8.69)                 
                    
            #endfor RNUM
            
             
            N_random_gals = 20000
            if(N_random_gals>NGals):
                N_random_gals=NGals
            random.seed(a=1)    
            random_list = random.sample(range(0, NGals), N_random_gals)
            
            for ii in random_list:
                slice_ii = [x*NGals+ii for x in range(0,12)]                
                xx = x_variable[slice_ii]
                yy = np.log10(y_variable[slice_ii]) 
              
                #ignore galaxies without halflightradius or with nan on the y_variable
                sel = (~np.isnan(xx)) & (~np.isinf(xx)) & (~np.isnan(yy)) & (~np.isinf(yy))                 
                if(len(xx[sel])>3):                   
                    f = interpolate.UnivariateSpline(xx[sel], yy[sel], s=0)          
                    interpol_y_variable[ii*len(new_x_var):(ii+1)*len(new_x_var)] = f(new_x_var)
                    interpol_x_variable[ii*len(new_x_var):(ii+1)*len(new_x_var)] = new_x_var
                        
            sel = interpol_x_variable>0.            
            (x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles_fixed_xx(interpol_x_variable[sel],
                                                                                      interpol_y_variable[sel])
        
            subplot.plot(x_binned, median, color=plot_color[k_mass], linewidth=2)
          
        #end loop on mass bins    
                    
    
        #Observations
        if(plot_prop==2 and k_mass ==MassBins-1):
            df = pd.read_csv(Datadir+'erroz_ferrer_dopita_bin1.csv')
            subplot.errorbar(df['radius'], df['metallicity'], xerr=0.1, yerr=[-df['error_down'],df['error_up']], 
                             fmt='o', markersize=5, color=plot_color[1], ecolor=plot_color[1])
            
            df = pd.read_csv(Datadir+'erroz_ferrer_dopita_bin3.csv')
            subplot.errorbar(df['radius'], df['metallicity'], xerr=0.1, yerr=[-df['error_down'],df['error_up']], 
                             fmt='o', markersize=5, color=plot_color[4], ecolor=plot_color[4])
      
        #labels 
        if(plot_prop==0):
            plot_label (subplot, 'label',xlim,ylim,x_percentage=0.17,y_percentage=0.9, 
                        color='black',xlog=0,ylog=0,label='$\log_{10}(M_*[M_{\odot}])=$', fontsize=12,fontweight='normal') 
            
            for k_mass in range(0,MassBins):
                char_mass_low="%0.1f" % mass_low[k_mass]
                char_mass_high="%0.1f" % mass_high[k_mass]
                x_values=[0.54, 0.745, 0.2, 0.46, 0.72]
                y_values=[0.9,0.9,0.8,0.8,0.8]
        
                label='['+char_mass_low+','+char_mass_high+']'
                plot_label (subplot, 'label',xlim,ylim,x_percentage=x_values[k_mass],y_percentage=y_values[k_mass], 
                color=plot_color[k_mass],xlog=0,ylog=0,label=label,fontsize=12,fontweight='normal') 
        
        
        
  
        
    #end loop on properties to plot
       
        
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYJ18_gas_gradients.pdf')
    plt.close()
    
    return 
#end gas_gradients

def MANGA_CALIFA_gradients(ThisRedshiftList):
    ii=0   
    
    plot_color=['blue', 'green', 'red']    
    fig = plt.figure(figsize=(three_one_size_small[0],three_one_size_small[1]))
    grid = gridspec.GridSpec(3, 1)
    grid.update(wspace=0.0, hspace=0.0)
   
    MassBins = 3
    mass_low = [9.0,10.0,10.5]
    mass_high = [10.0,10.5,11.0]
  
      
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)                 
    G0_MR_unsel = G_MR[sel] 
    G0_MR_unsel = G0_MR_unsel[np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)>9.0]
    
    #SFH structure needed for age gradients
    fa = open(DirName_MR+"SFH_Bins","rb")                
    nbins =  np.fromfile(fa,np.int32,1)
    template = np.dtype([('SnapNum',np.int32,1), ('Bin',np.int32,1), ('Lookbacktime',np.float64,1),                           
                         ('dt',np.float64,1),('nbins',np.int32,1)])
    SFH = np.fromfile(fa,template,int(nbins))    
    fa.close()                
    SFH=SFH[SFH['SnapNum']==G0_MR_unsel['SnapNum'][0]]             
        
        
    i_grid=0
    Nprops=3
    #loop on quantity to plot
    for plot_prop in range (0, Nprops):
                        
        subplot=plt.subplot(grid[i_grid])
        i_grid+=1
                        
        if(plot_prop==0):
            #***************************    
            #* stellar surface density *
            #*************************** 
            ylim=[0.5, 5.0]     
            ylab='$\log_{10}(\Sigma_*[M_{\odot}/\mathrm{pc^2}])$'   
            subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
            subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
            
        if(plot_prop==1):
            #***************************    
            #*    Mass-weighted age    *
            #***************************                  
            ylim=[3.0, 9.5]       
            ylab='$Age_{MW}(Gyr)$'   
            subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
            subplot.yaxis.set_minor_locator(MultipleLocator(0.5)) 
                     
        if(plot_prop==2):                  
            #***************************    
            #*   stellar metallicity   *
            #*************************** 
            ylim=[-1., 0.4]     
            ylab='$\mathrm{log_{10}}(Z_*/Z_{\odot})$'
            subplot.yaxis.set_major_locator(MultipleLocator(0.5))    
            subplot.yaxis.set_minor_locator(MultipleLocator(0.1))
             
        xlim=[0.0, 2.0]
        subplot.set_ylim(ylim),subplot.set_xlim(xlim)
        xlab='$r/r_{e}$'
        subplot.set_xlabel(xlab,fontsize=14), subplot.set_ylabel(ylab,fontsize=14)                      
        if((plot_prop==0) or (plot_prop==1)):
                plt.tick_params(axis='x', which='both', bottom=True, labelbottom=False)   
          
        
        #OBSERVATIONS 
        CALIFA_mass_low=[9.1,10.1,10.6]
        CALIFA_mass_high=[9.6,10.6,10.9]
        CALIFA_mass_index=[0,2,3]  
        obs_rings=np.arange(0.05,2.8,0.1)   
        
        #CALIFA (surface densities and light weighted ages)
        if(plot_prop==0):         
            for k_mass in range(0,len(CALIFA_mass_low)):
                char_mass_low="%0.1f" % CALIFA_mass_low[k_mass]
                char_mass_high="%0.1f" % CALIFA_mass_high[k_mass]
                char_k_mass="%d" %  CALIFA_mass_index[k_mass]
                  
                file = Datadir+'/CALIFA_rosa_stellar_density_'+char_k_mass+'_'+char_mass_low+'_'+char_mass_high+'.txt'
                obs = Table.read(file, format='ascii')      
                sel=obs['stellar_dens']!=0.0
                subplot.plot(obs_rings[sel],obs['stellar_dens'][sel], color=plot_color[k_mass],linestyle='--',linewidth=2)          
        #MANGA 
        if(plot_prop==1 or plot_prop==2):                
            file = Datadir + '/MANGA_gradients/goddard2016_LT_new_MAD.txt'
            obs = Table.read(file, format='ascii') 
           
            if(plot_prop==1):
                obs_radius = np.array([0.121, 0.407, 0.663, 0.836, 1.069, 1.240, 1.435])
                obs_MW_age  = np.array([0.765, 0.733, 0.717, 0.702, 0.693, 0.693, 0.691])
                subplot.plot(obs_radius, 10**obs_MW_age, linewidth=2, linestyle='--', color=plot_color[0])
                obs_radius = np.array([0.109, 0.245, 0.345, 0.523, 0.703, 0.894, 1.100, 1.281, 1.439])
                obs_MW_age  = np.array([0.736, 0.719, 0.720, 0.674, 0.676, 0.684, 0.689, 0.684, 0.681])
                subplot.plot(obs_radius, 10**obs_MW_age, linewidth=2, linestyle='--', color=plot_color[1])
                obs_radius = np.array([0.113, 0.225, 0.320, 0.472, 0.632, 0.765, 0.961, 1.109, 1.246, 1.389, 1.440])
                obs_MW_age  = np.array([0.684, 0.708, 0.732, 0.721, 0.699, 0.693, 0.692, 0.710, 0.724, 0.726, 0.731])
                subplot.plot(obs_radius, 10**obs_MW_age, linewidth=2, linestyle='--', color=plot_color[2])

            elif(plot_prop==2):
                subplot.plot(obs['Radius'],obs['MW_Metallicity_M1'], linewidth=2, linestyle='--', color=plot_color[0])
                subplot.plot(obs['Radius'],obs['MW_Metallicity_M2'], linewidth=2, linestyle='--', color=plot_color[1]) 
                subplot.plot(obs['Radius'],obs['MW_Metallicity_M3'], linewidth=2, linestyle='--', color=plot_color[2])
        
        #plot_color=['darkblue','blue','lightblue','green','orange','red','brown']    
        #fig = plt.figure(figsize=(three_one_size_small[0],three_one_size_small[1]))
        #grid = gridspec.GridSpec(3, 1)
        #grid.update(wspace=0.0, hspace=0.0)
        
        #For comparison with MaNGA metallicity gradients select only galaxies with non-zero metallicity    
        if(plot_prop==2):  
            if(opt_detailed_enrichment==1):   
                G0_MR_unsel=G0_MR_unsel[((G0_MR_unsel['MetalsDiskMass'][:,0] + 
                              G0_MR_unsel['MetalsDiskMass'][:,1] + 
                              G0_MR_unsel['MetalsDiskMass'][:,2])>.0)] 
            else:
                G0_MR_unsel=G0_MR_unsel[(G0_MR_unsel['MetalsDiskMass']>.0)]             
        
        #For ManGa comparison with metallicity and age gradients select only late-type galaxies
        if(plot_prop==1 or plot_prop==2):       
            G0_MR_unsel=G0_MR_unsel[(G0_MR_unsel['DiskMass']/G0_MR_unsel['StellarMass']>0.8)]
               
        #loop on mass bins    
        for k_mass in range(0,len(mass_low)):             
            G0_MR=G0_MR_unsel[(np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)>mass_low[k_mass]) &
                              (np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)<mass_high[k_mass])]
         
            NGals=len(G0_MR)            
            x_variable=np.zeros(int(RNUM*NGals),dtype=np.float32)             
            y_variable=np.zeros(int(RNUM*NGals),dtype=np.float32)   
            new_x_var = np.arange(xlim[0],xlim[1],0.1)
            interpol_x_variable=np.zeros(int(len(new_x_var)*NGals),dtype=np.float32)      
            interpol_y_variable=np.zeros(int(len(new_x_var)*NGals),dtype=np.float32)       
            
            #loop on radial bins
            for jj in range(0,RNUM):                              
                x_variable[NGals*jj:NGals*(jj+1)]=RingRadius[jj]/(G0_MR['StellarHalfLightRadius']*1000./Hubble_h)
                                        
                #***************************    
                #* stellar surface density *
                #***************************
                if(plot_prop==0):                     
                    StellarMass_this_ring=(G0_MR['DiskMassRings'][:,jj]*1e10/Hubble_h + 
                                           G0_MR['BulgeMassRings'][:,jj]*1e10/Hubble_h)                      
                    #1e6 -> from kpc^2 to pc^2
                    if(jj==0):
                        y_variable[NGals*jj:NGals*(jj+1)]=StellarMass_this_ring/(3.14*RingRadius[0]**2*1e6) 
                    else:
                        y_variable[NGals*jj:NGals*(jj+1)]=StellarMass_this_ring/(3.14*(RingRadius[jj]**2-RingRadius[jj-1]**2)*1e6)                            
                    
                #***************************    
                #*    Mass-weighted AGE    *
                #***************************
                if(plot_prop==1):                      
                    age=np.zeros(NGals)
                    for ii in range(0,len(SFH)):                        
                        age+=SFH['Lookbacktime'][ii]*(G0_MR['sfh_DiskMassRings'][:,jj,ii]*(1.-0.43))
                             
                    y_variable[NGals*jj:NGals*(jj+1)] = age / G0_MR['DiskMassRings'][:,jj] / 1.e9
                
                   
                #***************************    
                #*   stellar metallicity   *
                #***************************   
                if(plot_prop==2): 
                    if(opt_detailed_enrichment==1):                  
                        MetalsDiskMass_this_ring=(G0_MR['MetalsDiskMassRings'][:,jj,0] + 
                                                  G0_MR['MetalsDiskMassRings'][:,jj,1] +
                                                  G0_MR['MetalsDiskMassRings'][:,jj,2])
                    else:            
                        MetalsDiskMass_this_ring=G0_MR['MetalsDiskMassRings'][:,jj]
                
                    #x_variable[NGals*jj:NGals*(jj+1)]=RingRadius[jj]/(G0_MR['StellarDiskRadius']*1000./Hubble_h)
                    #MetalsStellarMass_this_ring=(MetalsDiskMass_this_ring+MetalsBulgeMass_this_ring)*1e10/Hubble_h
                    #y_variable[NGals*jj:NGals*(jj+1)]= MetalsStellarMass_this_ring/StellarMass_this_ring/0.02 
                    y_variable[NGals*jj:NGals*(jj+1)]= MetalsDiskMass_this_ring/G0_MR['DiskMassRings'][:,jj]/0.02    
                    
                    
            #endfor RNUM
            
            
            #do plots for different properties
            '''bin=0.25
            sel=y_variable>0. 
            if(len(y_variable[sel])>0.):
                (x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles (bin, xlim[0]+0.1, xlim[1], 
                                                                                      x_variable[sel], y_variable[sel])
              
            if(plot_prop==1):
                subplot.plot(x_binned, mean, color=plot_color[k_mass], linewidth=2)
                #subplot.plot(x_binned, pc16, color=plot_color[k_mass], linewidth=2, linestyle=':')
                #subplot.plot(x_binned, pc84, color=plot_color[k_mass], linewidth=2, linestyle=':')
            else:
                subplot.plot(x_binned, np.log10(mean), color=plot_color[k_mass], linewidth=2)'''
             
            N_random_gals = 5000
            if(N_random_gals>NGals):
                N_random_gals=NGals
            random_list = random.sample(range(0, NGals), N_random_gals)
            
            for ii in random_list:                  
                slice_ii = [x*NGals+ii for x in range(0,12)]
                
                xx = x_variable[slice_ii]
                if(plot_prop==1):
                    yy = y_variable[slice_ii] 
                else:
                    yy = np.log10(y_variable[slice_ii]) 
                    
                #ignore galaxies without halflightradius or with nan on the y_variable
                sel = (~np.isnan(xx)) & (~np.isinf(xx)) & (~np.isnan(yy)) & (~np.isinf(yy))               
                if(len(xx[sel])>0):                   
                    f = interpolate.UnivariateSpline(xx[sel], yy[sel], s=0)          
                    interpol_y_variable[ii*len(new_x_var):(ii+1)*len(new_x_var)] = f(new_x_var)
                    interpol_x_variable[ii*len(new_x_var):(ii+1)*len(new_x_var)] = new_x_var
            
            
            '''
            #without interpolation
            bin=0.25
            sel=y_variable>0. 
            if(len(y_variable[sel])>0.):
                (x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles (bin, xlim[0]+0.1, xlim[1], 
                                                                                  x_variable[sel], y_variable[sel])'''
          
            sel = interpol_x_variable>0.
            (x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles_fixed_xx(interpol_x_variable[sel],
                                                                                      interpol_y_variable[sel])
            subplot.plot(x_binned, median, color=plot_color[k_mass], linewidth=2)
            
        #end loop on mass bins    
                         
    
        #labels 
        if(plot_prop==0):
            plot_label (subplot, 'label',xlim,ylim,x_percentage=0.07,y_percentage=0.9, 
                        color='black',xlog=0,ylog=0,label='$M_{\odot}=$', fontsize=12,fontweight='normal') 
            
            for k_mass in range(0,MassBins):
                char_mass_low="%0.1f" % mass_low[k_mass]
                char_mass_high="%0.1f" % mass_high[k_mass]
                x_values=[0.2,0.44,0.71]
                y_values=[0.9,0.9,0.9]
        
                label='['+char_mass_low+','+char_mass_high+']'
                plot_label (subplot, 'label',xlim,ylim,x_percentage=x_values[k_mass],y_percentage=y_values[k_mass], 
                color=plot_color[k_mass],xlog=0,ylog=0,label=label,fontsize=12,fontweight='normal') 
        
        if(plot_prop==0): 
            yy = 0.8
            obs_label = 'CALIFA'
        else:
            yy = 0.9
            obs_label = 'MaNGA'
            
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.69, y_percentage=yy, 
                    color='black', xlog=0, ylog=0,label=prefix_this_model, fontsize=13, fontweight='normal')        
        plot_label (subplot,'line',xlim,ylim,x_percentage=0.6,y_percentage=yy+0.02,
                    color='red',x2_percentage=0.67,xlog=0,ylog=0,linestyle='-',linewidth=2)
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.69, y_percentage=yy-0.1, 
                    color='black', xlog=0, ylog=0,label=obs_label, fontsize=13, fontweight='normal')        
        plot_label (subplot,'line',xlim,ylim,x_percentage=0.6,y_percentage=yy-0.08,
                    color='red',x2_percentage=0.67,xlog=0,ylog=0,linestyle='--',linewidth=2)
  
        
    #end loop on properties to plot
       
        
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYJ18_MANGA_CALIFA_gradients.pdf')
    plt.close()
    
    return 
#end MANGA_CALIFA_gradients











def CALIFA_gradients_morph_types(ThisRedshiftList):
      
    ii=0   
    
    plot_color=['brown','red','orange','green','lightblue','blue','darkblue']    
    fig = plt.figure(figsize=(three_one_size_small[0],three_one_size_small[1]))
    grid = gridspec.GridSpec(3, 1)
    grid.update(wspace=0.0, hspace=0.0)
                
    disk_fraction=[0.0,0.05,0.1,0.2,0.4,0.6,0.8,1.0]    
    morph_types=['E','S0','Sa','Sb','Sbc','Sc','Sd'] 
    obs_morph_types=['0_E','1_S0','2_Sa','3_Sb','4_Sbc','5_Sc','6_Sd'] 
    mass_low=[11.2,10.9,10.9,10.9,10.6,10.1,9.1]
    mass_high=[11.5,11.2,11.2,11.2,10.9,10.6,9.6]
  
    Nprops=3
    for plot_prop in range (0,Nprops):
    
          
        subplot=plt.subplot(grid[plot_prop])     
       
        xlim=[0.0,3.]
        xlab='$r/r_{d}$'
        subplot.set_xlabel(xlab,fontsize=14) 
        subplot.xaxis.set_major_locator(MultipleLocator(1.0))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.1))

        if(plot_prop==0):
            #***************************    
            #*   stellar metallicity   *
            #*************************** 
            ylim=[-1., 0.5]     
            ylab='$\mathrm{log_{10}}(Z_*/Z_{\odot})$'
            subplot.yaxis.set_major_locator(MultipleLocator(0.5))    
            subplot.yaxis.set_minor_locator(MultipleLocator(0.1))
     
        if(plot_prop==1):
            #***************************    
            #* stellar surface density *
            #*************************** 
            ylim=[0.5, 4.5]     
            ylab='$\log_{10}(\Sigma_*[M_{\odot}/\mathrm{pc^2}])$'   
            subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
            subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
            
        if(plot_prop==2):
            #***************************    
            #*    mass-weighted age    *
            #*************************** 
            ylim=[8.0, 10.5] 
            #ylim=[9.0, 10.0] 
            ylab='$\log_{10}(Age_{LW}[yr])$'   
            subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
            subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
                
        subplot.set_ylim(ylim),subplot.set_xlim(xlim)                 
        subplot.set_ylabel(ylab,fontsize=14)   
        subplot.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        subplot.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        if((plot_prop==0) or (plot_prop==1)):
            plt.tick_params(axis='x', which='both', bottom='on', labelbottom='off')

            
        
        #OBSERVATIONS  
        for k_type in range(0,len(morph_types)):       
            obs_rings=np.arange(0.05,2.8,0.1)
            if(plot_prop==0):
                file = Datadir + '/CALIFA_rosa_metals_'+obs_morph_types[k_type]+'.txt'
                obs = Table.read(file, format='ascii')      
                sel=obs['metallicity']!=0.0      
                subplot.plot(obs_rings[sel],obs['metallicity'][sel], color=plot_color[k_type],linestyle='--',linewidth=2)
       
            if(plot_prop==1):              
                file = Datadir + '/CALIFA_rosa_stellar_density_'+obs_morph_types[k_type]+'.txt'
                obs = Table.read(file, format='ascii')      
                sel=obs['stellar_dens']!=0.0
                subplot.plot(obs_rings[sel],obs['stellar_dens'][sel], color=plot_color[k_type],linestyle='--',linewidth=2)  
        
            if(plot_prop==2):              
                file = Datadir + '/CALIFA_rosa_ageL_'+obs_morph_types[k_type]+'.txt'
                obs = Table.read(file, format='ascii')      
                sel=obs['age_L']!=0.0
                subplot.plot(obs_rings[sel],obs['age_L'][sel], color=plot_color[k_type],linestyle='--',linewidth=2)
        
            
        
        median_metallicity=np.zeros(RNUM,dtype=np.float32)
  
        #Model
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)                 
        G0_MR_unsel=G_MR[sel] 

        if(opt_detailed_enrichment==1):   
            G0_MR_unsel=G0_MR_unsel[((G0_MR_unsel['MetalsDiskMass'][:,0] + 
                                      G0_MR_unsel['MetalsDiskMass'][:,1] + 
                                      G0_MR_unsel['MetalsDiskMass'][:,2])>.0)] 
        else:
            G0_MR_unsel=G0_MR_unsel[(G0_MR_unsel['MetalsDiskMass']>.0)] #& 
           
        for k_type in range(0,len(morph_types)):       
            G0_MR=G0_MR_unsel[(G0_MR_unsel['DiskMass']/G0_MR_unsel['StellarMass']>disk_fraction[k_type]) &
                              (G0_MR_unsel['DiskMass']/G0_MR_unsel['StellarMass']<disk_fraction[k_type+1]) & 
                              (np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)>mass_low[k_type]) &
                              (np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)<mass_high[k_type])]
     
            NGals=len(G0_MR)
            if(NGals>0):
                x_variable=np.zeros(int(RNUM*NGals),dtype=np.float32) 
                y_variable=np.zeros(int(RNUM*NGals),dtype=np.float32)   
                
                r_bulge=G0_MR['BulgeSize']*1000./Hubble_h #From Mpc/h to kpc
       
                #SFH structure needed for age gradients
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
    
    
                for jj in range(0,RNUM):
            
                    x_variable[NGals*jj:NGals*(jj+1)]=RingRadius[jj]/(G0_MR['StellarHalfMassRadius']*1000./Hubble_h)
                
                    if(opt_rings_in_bulges==1):
                        BulgeMass_this_ring=G0_MR['BulgeMassRings'][:,jj]*1e10/Hubble_h
                    else:
                        if(jj==0):
                            r_bulge_m=1.-1./(1.+RingRadius[0]/r_bulge)
                        else:
                            r_bulge_m=(1/(1+RingRadius[jj-1]/r_bulge)-1/(1+RingRadius[jj]/r_bulge))
                        BulgeMass_this_ring=G0_MR['BulgeMass']*r_bulge_m*1e10/Hubble_h  
                        BulgeMass_this_ring[r_bulge==0.]=0.   
                    
                    StellarMass_this_ring=G0_MR['DiskMassRings'][:,jj]*1e10/Hubble_h+BulgeMass_this_ring    
            
                    #***************************    
                    #*   stellar metallicity   *
                    #*************************** 
                    if(plot_prop==0): 
                        #METALS    
                        if(opt_rings_in_bulges==1):
                            if(opt_detailed_enrichment==1):                  
                                MetalsBulgeMass_this_ring=(G0_MR['MetalsBulgeMassRings'][:,jj,0] + 
                                                          G0_MR['MetalsBulgeMassRings'][:,jj,1] +
                                                          G0_MR['MetalsBulgeMassRings'][:,jj,2])
                            else:            
                                MetalsBulgeMass_this_ring=G0_MR['MetalsBulgeMassRings'][:,jj]
                        else:
                            if(opt_detailed_enrichment==1): 
                                MetalsBulgeMass_this_ring=(G0_MR['MetalsBulgeMass'][:,0]+ 
                                                           G0_MR['MetalsBulgeMass'][:,1]+
                                                           G0_MR['MetalsBulgeMass'][:,2])*r_bulge_m
                            else:              
                                MetalsBulgeMass_this_ring=G0_MR['MetalsBulgeMass']*r_bulge_m            
                            MetalsBulgeMass_this_ring[r_bulge==0]=0. 
                
                        if(opt_detailed_enrichment==1):                  
                            MetalsDiskMass_this_ring=(G0_MR['MetalsDiskMassRings'][:,jj,0] + 
                                                      G0_MR['MetalsDiskMassRings'][:,jj,1] +
                                                      G0_MR['MetalsDiskMassRings'][:,jj,2])
                        else:            
                            MetalsDiskMass_this_ring=G0_MR['MetalsDiskMassRings'][:,jj]
                    
                        MetalsStellarMass_this_ring=(MetalsDiskMass_this_ring+MetalsBulgeMass_this_ring)*1e10/Hubble_h
                
                        y_variable[NGals*jj:NGals*(jj+1)]= MetalsStellarMass_this_ring/StellarMass_this_ring/0.02
              
    
            
                    #***************************    
                    #* stellar surface density *
                    #***************************
                    if(plot_prop==1): 
                        #1e6 -> from kpc^2 to pc^2
                        if(jj==0):
                            y_variable[NGals*jj:NGals*(jj+1)]=StellarMass_this_ring/(3.14*RingRadius[0]**2*1e6) 
                        else:
                            y_variable[NGals*jj:NGals*(jj+1)]=StellarMass_this_ring/(3.14*(RingRadius[jj]**2-RingRadius[jj-1]**2)*1e6)        
            
            
                    #***************************    
                    #*    mass-weighted AGE    *
                    #***************************
                    if(plot_prop==2):
                        #we only need the SFH strucutre from the current snap
                        SFH=SFH[SFH['SnapNum']==G0_MR['SnapNum'][0]]     
                        #if(jj==0):
                        #    print(np.log10(SFH['Lookbacktime']))
                        age=np.zeros(NGals)
                        for ii in range(0,len(SFH)):
                            sel=G0_MR['sfh_DiskMassRings'][:,jj,ii]>0.
                            age[sel]+=SFH['Lookbacktime'][ii]*(G0_MR['sfh_DiskMassRings'][sel,jj,ii]*
                                                               (1.-0.43))*1e10/Hubble_h
                            if(opt_rings_in_bulges==1):
                                sel=G0_MR['sfh_BulgeMassRings'][:,jj,ii]>0.
                                age[sel]+=SFH['Lookbacktime'][ii]*(G0_MR['sfh_BulgeMassRings'][sel,jj,ii]*
                                                                   (1.-0.43))*1e10/Hubble_h
                            else:     
                                age+=G0_MR['MassWeightAge']*1e9*BulgeMass_this_ring    
                            #sel=StellarMass_this_ring>0.
                        #age[sel]=age[sel]/StellarMass_this_ring[sel]               
                        #scale the massweighted ages by the global light weighted age
                        y_variable[NGals*jj:NGals*(jj+1)] = (age/StellarMass_this_ring *
                                                                 G0_MR['rBandWeightAge']/G0_MR['MassWeightAge'])                    
                        #y_variable[NGals*jj:NGals*(jj+1)]=age/(G0_MR['DiskMassRings'][:,jj]*1e10/Hubble_h)
                    
                        '''if(jj==0):
                            age=np.log10(np.mean(G0_MR['MassWeightAge']*1e9))        
                            subplot.scatter([(k_type+1)*0.1,(k_type+1)*0.1]
                                            ,[age,age],color=plot_color[k_type],marker='o',s=20)'''                    
                #endfor RNUM
            
                bin=0.1
                sel=y_variable>0. 
                if(len(y_variable[sel])>0.):
                    (x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles (bin, xlim[0], xlim[1], 
                                                                                      x_variable[sel], y_variable[sel])
                  
                subplot.plot(x_binned, np.log10(median), color=plot_color[k_type], linewidth=2)
              
        
                #labels    
                x_values=[0.1,0.16,0.25,0.34,0.43,0.56,0.65]
                if(plot_prop==0):
                    label=morph_types[k_type]
                    plot_label (subplot, 'label',xlim,ylim,x_percentage=x_values[k_type],y_percentage=0.85, 
                                color=plot_color[k_type],xlog=0,ylog=0,label=label,fontsize=15,fontweight='normal')            
                    
                if(plot_prop==1):
                    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.69, y_percentage=0.9, color='black', 
                                xlog=0, ylog=0,label=prefix_this_model, fontsize=13, fontweight='normal')        
                    plot_label (subplot,'line',xlim,ylim,x_percentage=0.6,y_percentage=0.92,
                                color='brown',x2_percentage=0.67,xlog=0,ylog=0,linestyle='-',linewidth=2)
                    plot_label (subplot, 'label', xlim, ylim, x_percentage=0.69, y_percentage=0.8, 
                                color='black', xlog=0, ylog=0,label='CALIFA', fontsize=13, fontweight='normal')        
                    plot_label (subplot,'line',xlim,ylim,x_percentage=0.6,y_percentage=0.82,
                                color='brown',x2_percentage=0.67,xlog=0,ylog=0,linestyle='--',linewidth=2)
                
                
        #endfor -> morph_types
    #endfor -> plot_prop
       
        
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYJ18_CALIFA_gradients_morph_types.pdf')
    plt.close()
    
    return 
#end gradients_morph_types





















def CALIFA_gradients_mass_bins(ThisRedshiftList):
    
    #plot_color=['darkblue','blue','lightblue','green','orange','red','brown']  
    plot_color=['darkblue','green','red'] 
    fig = plt.figure(figsize=(three_one_size_small[0],three_one_size_small[1]))
    grid = gridspec.GridSpec(3, 1)
    grid.update(wspace=0.0, hspace=0.0)
                   
    #mass_low=[9.1,9.6,10.1,10.6,10.9,11.2,11.5]
    #mass_high=[9.6,10.1,10.6,10.9,11.2,11.5,11.8]
    mass_low=[9.1,9.6,10.1]
    mass_high=[9.6,10.1,10.6]
  
    Nprops=3
    for plot_prop in range (0,Nprops):
    
          
        subplot=plt.subplot(grid[plot_prop])     
       
        xlim=[0.0,3.]
        xlab='$r/r_{e}$'
        subplot.set_xlabel(xlab,fontsize=14) 
        subplot.xaxis.set_major_locator(MultipleLocator(1.0))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.1))
    
        if(plot_prop==0):
            #***************************    
            #* stellar surface density *
            #*************************** 
            ylim=[0.5, 4.5]     
            ylab='$\log_{10}(\Sigma_*[M_{\odot}/\mathrm{pc^2}])$'   
            subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
            subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
            
        if(plot_prop==1):
            #***************************    
            #*    mass-weighted age    *
            #*************************** 
            ylim=[8.0, 10.5] 
            ylim=[0.0, 6.0] 
            ylab='$Age_{LW}(Gyr)$'   
            subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
            subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
        
        if(plot_prop==2):
            #***************************    
            #*   stellar metallicity   *
            #*************************** 
            ylim=[-1., 0.2]     
            ylab='$\mathrm{log_{10}}(Z_*/Z_{\odot})$'
            subplot.yaxis.set_major_locator(MultipleLocator(0.5))    
            subplot.yaxis.set_minor_locator(MultipleLocator(0.1))
        
        
        subplot.set_ylim(ylim),subplot.set_xlim(xlim)                 
        subplot.set_ylabel(ylab,fontsize=14)   
        subplot.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        subplot.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        if((plot_prop==0) or (plot_prop==1)):
            plt.tick_params(axis='x', which='both', bottom='on', labelbottom='off')

                           
        
        #OBSERVATIONS      
        obs_rings=np.arange(0.05,2.8,0.1) 
        k_mass_index=[0,1,2]
        for k_mass in range(0,len(mass_low)):      
            
            char_mass_low="%0.1f" % mass_low[k_mass]
            char_mass_high="%0.1f" % mass_high[k_mass]
            char_k_mass="%d" % k_mass_index[k_mass]
                   
            if(plot_prop==0):              
                file = Datadir+'/CALIFA_rosa_stellar_density_'+char_k_mass+'_'+char_mass_low+'_'+char_mass_high+'.txt'
                obs = Table.read(file, format='ascii')      
                sel=obs['stellar_dens']!=0.0
                subplot.plot(obs_rings[sel],obs['stellar_dens'][sel], color=plot_color[k_mass],linestyle='--',linewidth=2)  
        
            if(plot_prop==1):              
                file = Datadir + '/CALIFA_rosa_ageL_'+char_k_mass+'_'+char_mass_low+'_'+char_mass_high+'.txt'
                obs = Table.read(file, format='ascii')      
                sel=obs['age_L']!=0.0
                subplot.plot(obs_rings[sel],10**obs['age_L'][sel]/1.e9, color=plot_color[k_mass],linestyle='--',linewidth=2)
        
            if(plot_prop==2):
                file = Datadir + '/CALIFA_rosa_metals_'+char_k_mass+'_'+char_mass_low+'_'+char_mass_high+'.txt'
                obs = Table.read(file, format='ascii')      
                sel=obs['metallicity']!=0.0      
                subplot.plot(obs_rings[sel],obs['metallicity'][sel], color=plot_color[k_mass],linestyle='--',linewidth=2)
        
        
        median_metallicity=np.zeros(RNUM,dtype=np.float32)
  
        #Model       
        ii=0       
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)                 
        G0_MR_unsel = G_MR[sel] 
        G0_MR_unsel = G0_MR_unsel[np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)>9.0]
       
        #SFH structure needed for age gradients
        fa = open(DirName_MR+"SFH_Bins","rb")                
        nbins =  np.fromfile(fa,np.int32,1)
        template = np.dtype([('SnapNum',np.int32,1), ('Bin',np.int32,1), ('Lookbacktime',np.float64,1), 
                             ('dt',np.float64,1), ('nbins',np.int32,1)])
        SFH = np.fromfile(fa,template,int(nbins))    
        fa.close()          
                        
        for k_mass in range(0,len(mass_low)): 
            if(plot_prop==2):
                if(opt_detailed_enrichment==1):   
                    G0_MR_unsel=G0_MR_unsel[((G0_MR_unsel['MetalsDiskMass'][:,0] + 
                                              G0_MR_unsel['MetalsDiskMass'][:,1] + 
                                              G0_MR_unsel['MetalsDiskMass'][:,2])>.0)] 
                else:
                    G0_MR_unsel=G0_MR_unsel[(G0_MR_unsel['MetalsDiskMass']>.0)] 
            
            G0_MR=G0_MR_unsel[(np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)>mass_low[k_mass]) &
                              (np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)<mass_high[k_mass])]
            
            NGals=len(G0_MR)
            x_variable=np.zeros(int(RNUM*NGals),dtype=np.float32) 
            y_variable=np.zeros(int(RNUM*NGals),dtype=np.float32)    
            new_x_var = np.arange(xlim[0],xlim[1],0.1)
            interpol_x_variable=np.zeros(int(len(new_x_var)*NGals),dtype=np.float32)      
            interpol_y_variable=np.zeros(int(len(new_x_var)*NGals),dtype=np.float32)       
             
     
    
            for jj in range(0,RNUM):            
                x_variable[NGals*jj:NGals*(jj+1)]=RingRadius[jj]/(G0_MR['StellarHalfLightRadius']*1000./Hubble_h)            
                StellarMass_this_ring=G0_MR['DiskMassRings'][:,jj]*1e10/Hubble_h+G0_MR['BulgeMassRings'][:,jj]*1e10/Hubble_h
                              
                #***************************    
                #* stellar surface density *
                #***************************
                if(plot_prop==0): 
                    #1e6 -> from kpc^2 to pc^2
                    if(jj==0):
                        y_variable[NGals*jj:NGals*(jj+1)]=StellarMass_this_ring/(3.14*RingRadius[0]**2*1e6) 
                    else:
                        y_variable[NGals*jj:NGals*(jj+1)]=StellarMass_this_ring/(3.14*(RingRadius[jj]**2-RingRadius[jj-1]**2)*1e6)        
        
        
                #***************************    
                #*    mass-weighted AGE    *
                #***************************
                if(plot_prop==1):
                    #we only need the SFH strucutre from the current snap
                    SFH=SFH[SFH['SnapNum']==G0_MR['SnapNum'][0]]     
                    #if(jj==0):
                    #    print(np.log10(SFH['Lookbacktime']))
                    
                    age=np.zeros(NGals)
                    for ii in range(0,len(SFH)):
                        sel=G0_MR['sfh_DiskMassRings'][:,jj,ii]>0.
                        age[sel]+=SFH['Lookbacktime'][ii]*(G0_MR['sfh_DiskMassRings'][sel,jj,ii]*(1.-0.43))*1e10/Hubble_h
                       
                        sel=G0_MR['sfh_BulgeMassRings'][:,jj,ii]>0.
                        age[sel]+=SFH['Lookbacktime'][ii]*(G0_MR['sfh_BulgeMassRings'][sel,jj,ii]*(1.-0.43))*1e10/Hubble_h
                    
                    #to select only star forming rings
                    #sel_SFRrings = G0_MR['SfrRings'][:,jj]==0.
                    #age[sel_SFRrings] = 0.
                    
                    #scale the massweighted ages by the global light weighted age
                    y_variable[NGals*jj:NGals*(jj+1)] = (age/StellarMass_this_ring *
                                                         G0_MR['rBandWeightAge']/G0_MR['MassWeightAge'])/1.e9 
                    
                    #y_variable[NGals*jj:NGals*(jj+1)]=age/(G0_MR['DiskMassRings'][:,jj]*1e10/Hubble_h)/1.e9  
                
                
                #***************************    
                #*   stellar metallicity   *
                #*************************** 
                if(plot_prop==2): 
                    #METALS  
                    if(opt_detailed_enrichment==1):                  
                        MetalsBulgeMass_this_ring=(G0_MR['MetalsBulgeMassRings'][:,jj,0] + 
                                                  G0_MR['MetalsBulgeMassRings'][:,jj,1] +
                                                  G0_MR['MetalsBulgeMassRings'][:,jj,2])
                    else:            
                        MetalsBulgeMass_this_ring=G0_MR['MetalsBulgeMassRings'][:,jj]
                                
                    if(opt_detailed_enrichment==1):                  
                        MetalsDiskMass_this_ring=(G0_MR['MetalsDiskMassRings'][:,jj,0] + 
                                                  G0_MR['MetalsDiskMassRings'][:,jj,1] +
                                                  G0_MR['MetalsDiskMassRings'][:,jj,2])
                    else:            
                        MetalsDiskMass_this_ring=G0_MR['MetalsDiskMassRings'][:,jj]
                
                    MetalsStellarMass_this_ring=(MetalsDiskMass_this_ring+MetalsBulgeMass_this_ring)*1e10/Hubble_h             
                    y_variable[NGals*jj:NGals*(jj+1)]= MetalsStellarMass_this_ring/StellarMass_this_ring/0.02
          

                        
            #endfor RNUM
        
            #plot individual gradients
            '''if(plot_prop==2 and k_mass==0): 
                
                for ii in range(0,10):                     
                    slice_ii = [x*NGals+ii for x in range(0,12)]                 
                    xx = x_variable[slice_ii]
                    yy = np.log10(y_variable[slice_ii])                
                   
                    subplot.plot(xx, yy, color='b', linewidth=1, linestyle=':')
                 
                    sel = (~np.isnan(xx)) & (~np.isinf(xx)) & (~np.isnan(yy)) & (~np.isinf(yy)) 
                    print(xx[sel])
                    print(yy[sel])
                   
                    f = interpolate.UnivariateSpline(xx[sel], yy[sel], s=0)
                    subplot.plot(new_x_var, f(new_x_var), color='orange', linewidth=1, linestyle='-')'''
            
            N_random_gals = 5000
            if(N_random_gals>NGals):
                N_random_gals=NGals
            random_list = random.sample(range(0, NGals), N_random_gals)
            
            for ii in random_list:                  
                slice_ii = [x*NGals+ii for x in range(0,12)]
                
                xx = x_variable[slice_ii]
                if(plot_prop==1):
                    yy = y_variable[slice_ii] 
                else:
                    yy = np.log10(y_variable[slice_ii]) 
                    
                #ignore galaxies without halflightradius or with nan on the y_variable
                sel = (~np.isnan(xx)) & (~np.isinf(xx)) & (~np.isnan(yy)) & (~np.isinf(yy))               
                if(len(xx[sel])>0):                   
                    f = interpolate.UnivariateSpline(xx[sel], yy[sel], s=0)          
                    interpol_y_variable[ii*len(new_x_var):(ii+1)*len(new_x_var)] = f(new_x_var)
                    interpol_x_variable[ii*len(new_x_var):(ii+1)*len(new_x_var)] = new_x_var
            
            
            
            
            '''
            #without interpolation
            bin=0.25
            sel=y_variable>0. 
            if(len(y_variable[sel])>0.):
                (x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles (bin, xlim[0]+0.1, xlim[1], 
                                                                                  x_variable[sel], y_variable[sel])'''
            sel = interpol_x_variable>0.
            (x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles_fixed_xx(interpol_x_variable[sel],
                                                                                      interpol_y_variable[sel])
            subplot.plot(x_binned, median, color=plot_color[k_mass], linewidth=2)
               
          
        
            #labels   
            #char_mass_low="%0.1f" % mass_low[k_mass]
            #char_mass_high="%0.1f" % mass_high[k_mass]
            #x_values=[0.02,0.22,0.455,0.715,0.02,0.285,0.55]
            #y_values=[0.9,0.9,0.9,0.9,0.82,0.82,0.82]
            #if(plot_prop==0):
            #    label='['+char_mass_low+','+char_mass_high+']'
            #    plot_label (subplot, 'label',xlim,ylim,x_percentage=x_values[k_mass],y_percentage=y_values[k_mass], 
            #                color=plot_color[k_mass],xlog=0,ylog=0,label=label,fontsize=12,fontweight='normal') 
            
            #labels 
            if(plot_prop==0):
                plot_label (subplot, 'label',xlim,ylim,x_percentage=0.12,y_percentage=0.9, 
                            color='black',xlog=0,ylog=0,label='$M_{\odot}=$', fontsize=12,fontweight='normal') 
            
                for k_mass in range(0, len(mass_low)):
                    char_mass_low="%0.1f" % mass_low[k_mass]
                    char_mass_high="%0.1f" % mass_high[k_mass]
                    x_values=[0.26,0.47,0.71]
                    y_values=[0.9,0.9,0.9]
        
                    label='['+char_mass_low+','+char_mass_high+']'
                    plot_label (subplot, 'label',xlim,ylim,x_percentage=x_values[k_mass],y_percentage=y_values[k_mass], 
                    color=plot_color[k_mass],xlog=0,ylog=0,label=label,fontsize=12,fontweight='normal') 
            
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.69, y_percentage=0.8, 
                            color='black', xlog=0, ylog=0,label=prefix_this_model, fontsize=13, fontweight='normal')        
                plot_label (subplot,'line',xlim,ylim,x_percentage=0.6,y_percentage=0.82,
                            color='red',x2_percentage=0.67,xlog=0,ylog=0,linestyle='-',linewidth=2)
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.69, y_percentage=0.7, 
                            color='black', xlog=0, ylog=0,label='CALIFA', fontsize=13, fontweight='normal')        
                plot_label (subplot,'line',xlim,ylim,x_percentage=0.6,y_percentage=0.72,
                            color='red',x2_percentage=0.67,xlog=0,ylog=0,linestyle='--',linewidth=2)
        #endfor -> massbins
    #endfor -> plot_prop
       
        
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYJ18_CALIFA_gradients_mass_bins.pdf')
    plt.close()
    
    return 
#end gradients_mass_bins














def MANGA_gradients_late_types(ThisRedshiftList):
    
    ii=0   
    
    plot_color=['brown','red','orange','green','lightblue','blue','darkblue']    
    fig = plt.figure(figsize=(two_four_size_large[0],two_four_size_large[1]))
    grid = gridspec.GridSpec(2, 4)
    grid.update(wspace=0.0, hspace=0.0)
   
    MassBins=4
    mass_low=[9.0,9.935,10.552,11.054]
    mass_high=[9.935,10.552,11.054,11.5]
  
      
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)                 
    G0_MR_unsel=G_MR[sel] 

    if(opt_detailed_enrichment==1):   
        G0_MR_unsel=G0_MR_unsel[((G0_MR_unsel['MetalsDiskMass'][:,0] + 
                                  G0_MR_unsel['MetalsDiskMass'][:,1] + 
                                  G0_MR_unsel['MetalsDiskMass'][:,2])>.0)] 
    else:
        G0_MR_unsel=G0_MR_unsel[(G0_MR_unsel['MetalsDiskMass']>.0)]   
        
    #SFH structure needed for age gradients
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
    SFH=SFH[SFH['SnapNum']==G0_MR_unsel['SnapNum'][0]]             
        
    i_grid=0
    Nprops=2
    #loop on metallicity and age
    for plot_prop in range (0,Nprops):
        #loop on mass bins
        for i_mass in range(0,MassBins):       
            G0_MR=G0_MR_unsel[(G0_MR_unsel['DiskMass']/G0_MR_unsel['StellarMass']>0.8) &                             
                              (np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)>mass_low[i_mass]) &
                              (np.log10(G0_MR_unsel['StellarMass']*1e10/Hubble_h)<mass_high[i_mass])]
            
            
            subplot=plt.subplot(grid[i_grid])
            i_grid+=1
       
            xlim=[0.0,1.7]
            xlab='$r/r_{d}$'
            subplot.set_xlabel(xlab,fontsize=14) 
            subplot.xaxis.set_major_locator(MultipleLocator(1.0))    
            subplot.xaxis.set_minor_locator(MultipleLocator(0.1))

            if(plot_prop==0):
                #***************************    
                #*    mass-weighted age    *
                #*************************** 
                ylim=[0.5, 1.2]                
                ylab='$\log_{10}(Age_{MW}[Gyr])$'   
                subplot.yaxis.set_major_locator(MultipleLocator(1.0))    
                subplot.yaxis.set_minor_locator(MultipleLocator(0.1)) 
            
            if(plot_prop==1):
                #***************************    
                #*   stellar metallicity   *
                #*************************** 
                ylim=[-2., 0.2]     
                ylab='$\mathrm{log_{10}}(Z_*/Z_{\odot})$'
                subplot.yaxis.set_major_locator(MultipleLocator(0.5))    
                subplot.yaxis.set_minor_locator(MultipleLocator(0.1))
                            
            subplot.set_ylim(ylim),subplot.set_xlim(xlim)  
            if(i_mass==0):   
                subplot.set_ylabel(ylab,fontsize=14)   
            subplot.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            subplot.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            if(plot_prop==0):
                plt.tick_params(axis='x', which='both', bottom='on', labelbottom='off')
                
            if(i_mass>0):
                plt.tick_params(axis='y', which='both', left='on', labelleft='off')    

            NGals=len(G0_MR)           
            x_variable=np.zeros(int(RNUM*NGals),dtype=np.float32) 
            y_variable=np.zeros(int(RNUM*NGals),dtype=np.float32)   
                
            r_bulge=G0_MR['BulgeSize']*1000./Hubble_h #From Mpc/h to kpc
           
            for jj in range(0,RNUM):
            
                x_variable[NGals*jj:NGals*(jj+1)]=RingRadius[jj]/(G0_MR['StellarHalfMassRadius']*1000./Hubble_h)
                
                if(opt_rings_in_bulges==1):
                    BulgeMass_this_ring=G0_MR['BulgeMassRings'][:,jj]*1e10/Hubble_h
                else:
                    if(jj==0):
                        r_bulge_m=1.-1./(1.+RingRadius[0]/r_bulge)
                    else:
                        r_bulge_m=(1/(1+RingRadius[jj-1]/r_bulge)-1/(1+RingRadius[jj]/r_bulge))
                    BulgeMass_this_ring=G0_MR['BulgeMass']*r_bulge_m*1e10/Hubble_h  
                    BulgeMass_this_ring[r_bulge==0.]=0.      
                StellarMass_this_ring=G0_MR['DiskMassRings'][:,jj]*1e10/Hubble_h+BulgeMass_this_ring    
            
                #***************************    
                #*    mass-weighted AGE    *
                #***************************
                if(plot_prop==0):
                    #we only need the SFH strucutre from the current snap                               
                    age=np.zeros(NGals)
                    for ii in range(0,len(SFH)):
                        sel=G0_MR['sfh_DiskMassRings'][:,jj,ii]>0.
                        age[sel]+=SFH['Lookbacktime'][ii]/1e9*(G0_MR['sfh_DiskMassRings'][sel,jj,ii]
                                                               *(1.-0.43))*1e10/Hubble_h
                    age+=G0_MR['MassWeightAge']*BulgeMass_this_ring                   
                    y_variable[NGals*jj:NGals*(jj+1)] = age/StellarMass_this_ring              
                  
                #***************************    
                #*   stellar metallicity   *
                #*************************** 
                if(plot_prop==1): 
                    #METALS
                    if(opt_rings_in_bulges==1):
                        if(opt_detailed_enrichment==1):                  
                            MetalsBulgeMass_this_ring=(G0_MR['MetalsBulgeMassRings'][:,jj,0] + 
                                                      G0_MR['MetalsBulgeMassRings'][:,jj,1] +
                                                      G0_MR['MetalsBulgeMassRings'][:,jj,2])
                        else:            
                            MetalsBulgeMass_this_ring=G0_MR['MetalsBulgeMassRings'][:,jj]
                    else:
                        if(opt_detailed_enrichment==1): 
                            MetalsBulgeMass_this_ring=(G0_MR['MetalsBulgeMass'][:,0]+ 
                                                       G0_MR['MetalsBulgeMass'][:,1]+
                                                       G0_MR['MetalsBulgeMass'][:,2])*r_bulge_m
                        else:              
                            MetalsBulgeMass_this_ring=G0_MR['MetalsBulgeMass']*r_bulge_m            
                        MetalsBulgeMass_this_ring[r_bulge==0]=0.          
            
                    if(opt_detailed_enrichment==1):                  
                        MetalsDiskMass_this_ring=(G0_MR['MetalsDiskMassRings'][:,jj,0] + 
                                                  G0_MR['MetalsDiskMassRings'][:,jj,1] +
                                                  G0_MR['MetalsDiskMassRings'][:,jj,2])
                    else:            
                        MetalsDiskMass_this_ring=G0_MR['MetalsDiskMassRings'][:,jj]
                
                    MetalsStellarMass_this_ring=(MetalsDiskMass_this_ring+MetalsBulgeMass_this_ring)*1e10/Hubble_h
            
                    y_variable[NGals*jj:NGals*(jj+1)]= MetalsStellarMass_this_ring/StellarMass_this_ring/0.02
               
                              
            #endfor RNUM
        
            bin=0.1
            sel=y_variable>0. 
            if(len(y_variable[sel])>0.):
                (x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles (bin, xlim[0], xlim[1], 
                                                                                  x_variable[sel], y_variable[sel])
              
            subplot.plot(x_binned, np.log10(median), color='red', linewidth=2)
          
        
            #labels    
            '''x_values=[0.1,0.16,0.25,0.34,0.43,0.56,0.65]
            if(plot_prop==0):
                label=morph_types[k_type]
                plot_label (subplot, 'label',xlim,ylim,x_percentage=x_values[k_type],y_percentage=0.85, 
                            color=plot_color[k_type],xlog=0,ylog=0,label=label,fontsize=15,fontweight='normal')   '''         
                  
            #OBSERVATIONS
            if(plot_prop==0):
                if(i_mass==0):
                    file = Datadir + '/MANGA_gradients/goddard2016_LT_new_MAD.txt'
                    obs = Table.read(file, format='ascii')             
                    subplot.fill_between(obs['Radius'],obs['MW_Age_M1']+obs['MW_Age_Error_M1'],
                                         obs['MW_Age_M1']-obs['MW_Age_Error_M1'], facecolor='lightblue', 
                                         interpolate=True, alpha=0.4, edgecolor='steelblue') 
                    subplot.plot(obs['Radius'],obs['MW_Age_M1'],color='blue', linewidth=2)
                if(i_mass==1):    
                    subplot.fill_between(obs['Radius'],obs['MW_Age_M2']+obs['MW_Age_Error_M2'],
                                         obs['MW_Age_M2']-obs['MW_Age_Error_M2'], facecolor='lightblue', 
                                         interpolate=True, alpha=0.4, edgecolor='steelblue') 
                    subplot.plot(obs['Radius'],obs['MW_Age_M2'],color='blue', linewidth=2)
                if(i_mass==2):    
                    subplot.fill_between(obs['Radius'],obs['MW_Age_M3']+obs['MW_Age_Error_M3'],
                                         obs['MW_Age_M3']-obs['MW_Age_Error_M3'], facecolor='lightblue', 
                                         interpolate=True, alpha=0.4, edgecolor='steelblue') 
                    subplot.plot(obs['Radius'],obs['MW_Age_M3'],color='blue', linewidth=2)
                if(i_mass==3):    
                    subplot.fill_between(obs['Radius'],obs['MW_Age_M4']+obs['MW_Age_Error_M4'],
                                         obs['MW_Age_M4']-obs['MW_Age_Error_M4'], facecolor='lightblue', 
                                         interpolate=True, alpha=0.4, edgecolor='steelblue') 
                    subplot.plot(obs['Radius'],obs['MW_Age_M4'],color='blue', linewidth=2)
                    
              
            if(plot_prop==1):
                if(i_mass==0):
                    file = Datadir + '/MANGA_gradients/goddard2016_LT_new_MAD.txt'
                    obs = Table.read(file, format='ascii')             
                    subplot.fill_between(obs['Radius'],obs['MW_Metallicity_M1']+obs['MW_Metallicity_Error_M1'],
                                         obs['MW_Metallicity_M1']-obs['MW_Metallicity_Error_M1'], facecolor='lightblue', 
                                         interpolate=True, alpha=0.4, edgecolor='steelblue') 
                    subplot.plot(obs['Radius'],obs['MW_Metallicity_M1'],color='blue', linewidth=2)
                if(i_mass==1):    
                    subplot.fill_between(obs['Radius'],obs['MW_Metallicity_M2']+obs['MW_Metallicity_Error_M2'],
                                         obs['MW_Metallicity_M2']-obs['MW_Metallicity_Error_M2'], facecolor='lightblue', 
                                         interpolate=True, alpha=0.4, edgecolor='steelblue') 
                    subplot.plot(obs['Radius'],obs['MW_Metallicity_M2'],color='blue', linewidth=2)
                if(i_mass==2):    
                    subplot.fill_between(obs['Radius'],obs['MW_Metallicity_M3']+obs['MW_Metallicity_Error_M3'],
                                         obs['MW_Metallicity_M3']-obs['MW_Metallicity_Error_M3'], facecolor='lightblue', 
                                         interpolate=True, alpha=0.4, edgecolor='steelblue') 
                    subplot.plot(obs['Radius'],obs['MW_Metallicity_M3'],color='blue', linewidth=2)
                if(i_mass==3):    
                    subplot.fill_between(obs['Radius'],obs['MW_Metallicity_M4']+obs['MW_Metallicity_Error_M4'],
                                         obs['MW_Metallicity_M4']-obs['MW_Metallicity_Error_M4'], facecolor='lightblue', 
                                         interpolate=True, alpha=0.4, edgecolor='steelblue') 
                    subplot.plot(obs['Radius'],obs['MW_Metallicity_M4'],color='blue', linewidth=2)
                               
          
            #if(plot_prop==0):
            label="%0.2f" % mass_low[i_mass] + r'$<\mathrm{log_{10}}(M_*[M_{\odot}])<$' + "%0.2f" % mass_high[i_mass]
            plot_label (subplot, 'label',xlim,ylim,x_percentage=0.05,y_percentage=0.05, 
                        color='black',xlog=0,ylog=0,label=label,fontsize=11,fontweight='normal')
            
        #endfor -> morph_types
    #endfor -> plot_prop
       
        
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYJ18_MANGA_gradients_late_types.pdf')
    plt.close()
    
    return 
#end gradients_morph_types





def SFR_gradients(ThisRedshiftList):
  
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
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14) 
        subplot.plot(radius, mean_SFR,color='red', linewidth=2)      
                
        label="%0.1f" % low_mass_limits[kk] + "$<M_{\star}[M_{\odot}]<$" + "%0.1f"  % (low_mass_limits[kk]+massbin)
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.88, 
                    color='black', xlog=0, ylog=0, label=label, 
                    fontsize=13, fontweight='normal') 
        
        #endfor
    #endfor
   
        
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()
    
    return 
#end SFR_gradients









def gasfractions_Saintonge17(ThisRedshiftList):
          
    for ii in range(0,len(ThisRedshiftList)):      
         
        xlim=[9.5,11.5]
        ylim=[-2.0,0.5]
        bin=0.25
               
        fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))    
        subplot=plt.subplot()    
        subplot.set_ylim(ylim), subplot.set_xlim(xlim) 
        
        ylab='$\log_{10}(M_{\mathrm{gas}}/M_*)$' 
        subplot.set_ylabel(ylab, fontsize=14) 
            
        xlab='$\log_{10}(M_*/M_{\odot})$'           
        subplot.set_xlabel(xlab, fontsize=14)   
          
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(0.5))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.1))
        subplot.yaxis.set_minor_locator(MultipleLocator(0.1))
       
                      
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]    
        log_StellarMass=np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)
        log_SFR=np.log10(G0_MR['Sfr'])
        G0_MR=G0_MR[(log_StellarMass>7.) & (G0_MR['ColdGas']>-1e-30) & ~np.isnan(G0_MR['H2fraction']) &
                    #(np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h))>-11.)] 
                    ((log_SFR-log_StellarMass)>np.log10((1+ThisRedshiftList[ii])**2/(1.37e10/2.)) -1.0)]
        #G0_MR=G0_MR[(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)>7.) & (G0_MR['ColdGas']>0.)] 
       
      
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
        
        #HI
        #sel=(1.-G0_MR['H2fraction'])>0.               
        Fraction_HI=np.log10(G0_MR['ColdGas']*(1.-G0_MR['H2fraction'])*1.e10/Hubble_h)-StellarMass       
        (x_binned, median,mean,pc16,pc84,rms)=median_and_percentiles(bin,xlim[0],xlim[1],StellarMass,Fraction_HI)        
        subplot.plot(x_binned,median,color='red',linewidth=2)     
        subplot.plot(x_binned,pc16,color='red',linewidth=2,linestyle='--')
        subplot.plot(x_binned,pc84,color='red',linewidth=2,linestyle='--')
        
        #WRITE OUTPUT      
        if(write_to_file==1):
            df = pd.DataFrame({'Log10_M':x_binned, 'median':median, 'pc16':pc16, 'pc84':pc84})         
            df.to_csv(Datadir + file_to_write + 'GasFractions_Saintonge2017_HI' + 
                      str(f'_z{ThisRedshiftList[ii]:0.2f}')+'.csv', index=False)
        #df = pd.read_csv(Datadir + file_to_write + 'GasFractions_Saintonge2017_HI' + 
        #                 str(f'_z{ThisRedshiftList[ii]:0.2f}')+'.csv')
        #subplot.plot(df['Log10_M'], df['median'],color='black', linestyle='--') 
        
        #H2
        #sel=G0_MR['H2fraction']>0.
        Fraction_H2=np.log10(G0_MR['ColdGas']*(G0_MR['H2fraction'])*1.e10/Hubble_h)-StellarMass        
        (x_binned, median,mean,pc16,pc84,rms)=median_and_percentiles(bin,xlim[0],xlim[1],StellarMass,Fraction_H2)         
        subplot.plot(x_binned,median,color='blue',linewidth=2)     
        subplot.plot(x_binned,pc16,color='blue',linewidth=2,linestyle='--')
        subplot.plot(x_binned,pc84,color='blue',linewidth=2,linestyle='--')
            
        #WRITE OUTPUT      
        if(write_to_file==1):
            df = pd.DataFrame({'Log10_M':x_binned, 'median':median, 'pc16':pc16, 'pc84':pc84})         
            df.to_csv(Datadir + file_to_write + 'GasFractions_Saintonge2017_H2' + 
                      str(f'_z{ThisRedshiftList[ii]:0.2f}')+'.csv', index=False)
        #df = pd.read_csv(Datadir + file_to_write + 'GasFractions_Saintonge2017_H2' + 
        #                 str(f'_z{ThisRedshiftList[ii]:0.2f}')+'.csv')
        #subplot.plot(df['Log10_M'], df['median'],color='black', linestyle='--') 
            
     
       
       

        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.19, y_percentage=0.9, 
                    color='black', xlog=0, ylog=0, label='xGass/xCOLD GASS: ', fontsize=13, fontweight='normal')         
        plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.75, y_percentage=0.925, 
                    color='red', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.05) 
        plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.855, y_percentage=0.925, 
                    color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.05) 
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.77, y_percentage=0.9, 
                    color='black', xlog=0, ylog=0, label='HI   H$_2$', fontsize=13, fontweight='normal') 
        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.435, y_percentage=0.8, 
                    color='black', xlog=0, ylog=0, label=prefix_this_model+': ', fontsize=13, fontweight='normal')         
        plot_label (subplot, 'line', xlim, ylim,x_percentage=0.74, y_percentage=0.825, 
                    color='red', x2_percentage=0.76, xlog=0, ylog=0, linestyle='-', linewidth=2)
        plot_label (subplot, 'line', xlim, ylim,x_percentage=0.84, y_percentage=0.825, 
                    color='blue', x2_percentage=0.86, xlog=0, ylog=0, linestyle='-', linewidth=2)
        
        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.77, y_percentage=0.8, 
                    color='black', xlog=0, ylog=0, label='HI   H$_2$', fontsize=13, fontweight='normal') 

        #plot_label (subplot, 'label', xlim, ylim, x_percentage=0.6, y_percentage=0.85, 
        #            color='black', xlog=0, ylog=0, label='$M_{\mathrm{cold}}/M_*$', 
        #            fontsize=15, fontweight='normal') 
           
        
        #OBSERVATIONS 
        
        #HI      
        file = Datadir+"/Saintonge2017_HI.csv"   
        df = pd.read_csv(file)
        subplot.errorbar(df['x'], df['y'],yerr=[df['err_down'], df['err_up']],
                 fmt='o', markersize=5, ecolor='red', color='red',zorder=+3)  
        
        #HI2     
        file = Datadir+"/Saintonge2017_H2.csv"   
        df = pd.read_csv(file)
        subplot.errorbar(df['x'], df['y'],yerr=[df['err_down'], df['err_up']],
                 fmt='o', markersize=5, ecolor='blue', color='blue',zorder=+3)  

        
        '''file = Datadir+"/Saintonge2016_gasfraction.txt"   
        Saint16 = Table.read(file, format='ascii')  
        Saint16_mass=(Saint16['mass_bin_low']+Saint16['mass_bin_high'])/2.  
        #H2 
        #OBSERVATIONS PLOT        
        y_err=np.zeros(len(Saint16['fH2']),dtype=np.float32)
        y_err=[np.log10(Saint16['fH2']/(Saint16['fH2']-Saint16['fH2_err'])),
               np.log10((Saint16['fH2']+Saint16['fH2_err'])/Saint16['fH2'])]
        subplot.errorbar(Saint16_mass, np.log10(Saint16['fH2']),xerr=0.12,yerr=y_err,
                 fmt='o', markersize=5, ecolor='blue', color='blue',zorder=+3)

        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.6, y_percentage=0.85, 
                    color='black', xlog=0, ylog=0, label='$M_{\mathrm{H_2}}/M_*$', 
                    fontsize=15, fontweight='normal')'''
            
        '''xx = [9.407, 9.638, 9.848, 10.04, 10.24, 
              10.46, 10.67, 10.87, 11.07, 11.27]
        yy = [-1.01,-1.08,-0.97,-0.90,-1.08,
             -1.05,-1.34,-1.41,-1.66,-2.02]
        
        xx = [9.388,9.669, 9.915, 10.10, 10.25,
              10.41,10.58, 10.75,10.90,11.08]
        yy = [-1.11,-1.28,-1.23,-1.16,-1.39,
              -1.40,-1.62,-1.75,-1.83,-2.01]
        
        xx = [9.406, 9.643,9.854,10.07,10.25,
             10.46,10.68,10.82,11.10,11.31]
        yy = [-0.93,-0.89,-0.96,-0.89,-1.05,
             -0.99,-1.10,-1.26,-1.41,-1.68]
        
        xx = [9.366,9.465,9.467,9.882,10.01,
              10.14,10.37,10.56,10.77,11.12]
        yy=[-0.99,-1.05,-1.12,-1.10,-1.11,
           -1.20,-1.27,-1.39,-1.44,-1.70]
        subplot.scatter(xx,yy, color='black')'''
            
                
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYF19_gasfractions_Saintonge17.pdf')
    plt.close()

    return 
#end



def gasfractions_vs_stellarmass(ThisRedshiftList):
    
    
    labels_to_write=['ColdGas', 'HI', 'H2', 'H2_HI']
    
    plot_color=['red','purple']        
   
    fig = plt.figure(figsize=(one_four_size_large[0],one_four_size_large[1]))
    grid = gridspec.GridSpec(1, 4)
    grid.update(wspace=0.0, hspace=0.0)

    #OBSERVATIONS READ  
    file = Datadir+"/Saintonge2016_gasfraction.txt"   
    Saint16 = Table.read(file, format='ascii')  
    Saint16_mass=(Saint16['mass_bin_low']+Saint16['mass_bin_high'])/2.  
  
     
    
    for ii in range(0,len(ThisRedshiftList)):        
               
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]    
        #G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (G0_MR['ColdGas']>0.) & 
        #            (G0_MR['Vvir']>120.) & (G0_MR['BulgeMass']/G0_MR['StellarMass']<0.15)]
        G0_MR=G0_MR[(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)>7.) & (G0_MR['ColdGas']>0.) & 
                    (np.log10(G0_MR['Sfr']/(G0_MR['StellarMass']*1.e10/Hubble_h))>-11.)] 
        #G0_MR=G0_MR[(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h)>10.) & (G0_MR['ColdGas']>0.)]
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
        
        
        xlim=[9.5,11.5]
        ylim=[-2.0,0.5]
        bin=0.25
        
        for i_gas in range(0,4):
        
            subplot=plt.subplot(grid[i_gas])    
            subplot.set_ylim(ylim), subplot.set_xlim(xlim) 
        
            if i_gas>0:
                plt.tick_params(axis='y', which='both', left='on', labelleft='off')
                
            if i_gas==0:    
                ylab='$\log_{10}(M_{\mathrm{gas}}/M_*)$' 
                subplot.set_ylabel(ylab, fontsize=14) 
            
            xlab='$\log_{10}(M_*/M_{\odot})$'           
            subplot.set_xlabel(xlab, fontsize=14)   
          
            #format axis
            majorFormatter = FormatStrFormatter('%d')
            subplot.xaxis.set_major_locator(MultipleLocator(1))    
            subplot.xaxis.set_minor_locator(MultipleLocator(0.25))
        
            
            #MODEL
            if(i_gas==0):
                sel=G0_MR['ColdGas']>0.
                Fraction=np.log10(G0_MR['ColdGas']*1.e10/Hubble_h)-StellarMass 
            if(i_gas==1):
                sel=(1.-G0_MR['H2fraction'])>0.               
                Fraction=np.log10(G0_MR['ColdGas']*(1.-G0_MR['H2fraction'])*1.e10/Hubble_h)-StellarMass         
            if(i_gas==2):
                sel=G0_MR['H2fraction']>0.
                Fraction=np.log10(G0_MR['ColdGas']*(G0_MR['H2fraction'])*1.e10/Hubble_h)-StellarMass         
            if(i_gas==3):
                sel=G0_MR['H2fraction']>0.
                Fraction=np.log10(G0_MR['H2fraction']/(1.-G0_MR['H2fraction']))
        
    
            (x_binned, median,mean,pc16,pc84,rms)=median_and_percentiles(bin,xlim[0],xlim[1],StellarMass[sel],Fraction[sel])    
            sel=(median!=0)        
            subplot.plot(x_binned[sel],median[sel],color=plot_color[ii],linewidth=2)     
            subplot.plot(x_binned[sel],pc16[sel],color=plot_color[ii],linewidth=2,linestyle='--')
            subplot.plot(x_binned[sel],pc84[sel],color=plot_color[ii],linewidth=2,linestyle='--')
            
            #WRITE OUTPUT      
            if(write_to_file==1):
                df = pd.DataFrame({'Log10_M':x_binned[sel], 'median':median[sel], 'pc16':pc16[sel], 'pc84':pc84[sel]})         
                df.to_csv(Datadir + file_to_write + 'GasFractions_' + labels_to_write[i_gas] + 
                          str(f'_z{ThisRedshiftList[ii]:0.2f}')+'.csv', index=False)
        
                #df = pd.read_csv(Datadir + file_to_write + 'GasFractions_'  + labels_to_write[i_gas] + 
                #                 str(f'_z{ThisRedshiftList[ii]:0.2f}')+'.csv')
                #subplot.plot(df['Log10_M'], df['median'],color='black', linestyle='--') 
            
            #OBSERVATIONS PLOT  
            if(i_gas==0):
                #Cold
                y_err=np.zeros(len(Saint16['fHI']),dtype=np.float32)
                Saint16_H2plusHI=Saint16['fH2']+Saint16['fHI']
                Saint16_H2plusHI_err=Saint16['fH2_err']+Saint16['fHI_err']
       
                y_err=[np.log10(Saint16_H2plusHI/(Saint16_H2plusHI-Saint16_H2plusHI_err)),
                       np.log10((Saint16_H2plusHI+Saint16_H2plusHI_err)/Saint16_H2plusHI)]
                subplot.errorbar(Saint16_mass, np.log10(Saint16_H2plusHI),xerr=0.12,yerr=y_err,
                                 fmt='o', markersize=5, ecolor='blue', color='blue',zorder=+3,capsize=2)    
               
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.15, y_percentage=0.2, 
                            color='black', xlog=0, ylog=0, label=prefix_this_model, 
                            fontsize=13, fontweight='normal')             
                plot_label (subplot, 'line', xlim, ylim,x_percentage=0.05, y_percentage=0.2175, 
                            color='red', x2_percentage=0.12, xlog=0, ylog=0, linestyle='-', linewidth=2)
        
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.15, y_percentage=0.12, 
                            color='black', xlog=0, ylog=0, label='Saintonge2016', 
                            fontsize=13, fontweight='normal') 
                plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.12, y_percentage=0.14, 
                            color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.075)     
           
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.6, y_percentage=0.85, 
                            color='black', xlog=0, ylog=0, label='$M_{\mathrm{cold}}/M_*$', 
                            fontsize=15, fontweight='normal') 
           
            if(i_gas==1):   
                #HI                
                #OBSERVATIONS PLOT        
                y_err=np.zeros(len(Saint16['fHI']),dtype=np.float32)
                y_err=[np.log10(Saint16['fHI']/(Saint16['fHI']-Saint16['fHI_err'])),
                       np.log10((Saint16['fHI']+Saint16['fHI_err'])/Saint16['fHI'])]
                subplot.errorbar(Saint16_mass, np.log10(Saint16['fHI']),xerr=0.12,yerr=y_err,
                         fmt='o', markersize=5, ecolor='blue', color='blue',zorder=+3)  
                         
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.6, y_percentage=0.85, 
                            color='black', xlog=0, ylog=0, label='$M_{\mathrm{HI}}/M_*$', 
                            fontsize=15, fontweight='normal')
                
            if(i_gas==2):   
                #H2 
                #OBSERVATIONS PLOT        
                y_err=np.zeros(len(Saint16['fH2']),dtype=np.float32)
                y_err=[np.log10(Saint16['fH2']/(Saint16['fH2']-Saint16['fH2_err'])),
                       np.log10((Saint16['fH2']+Saint16['fH2_err'])/Saint16['fH2'])]
                subplot.errorbar(Saint16_mass, np.log10(Saint16['fH2']),xerr=0.12,yerr=y_err,
                         fmt='o', markersize=5, ecolor='blue', color='blue',zorder=+3)
        
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.6, y_percentage=0.85, 
                            color='black', xlog=0, ylog=0, label='$M_{\mathrm{H_2}}/M_*$', 
                            fontsize=15, fontweight='normal')
            
            if(i_gas==3):   
                #H2/HI                
                #OBSERVATIONS PLOT        
                y_err=np.zeros(len(Saint16['fHI']),dtype=np.float32)
                Saint16_H2overHI=Saint16['fH2']/Saint16['fHI']
                Saint16_H2overHI_err=Saint16['fH2_err']/Saint16['fH2']+Saint16['fHI_err']/Saint16['fHI']
                sel=Saint16_H2overHI_err>Saint16_H2overHI
                Saint16_H2overHI_err[sel]=Saint16_H2overHI[sel]-0.01

                y_err=[np.log10(Saint16_H2overHI/(Saint16_H2overHI-Saint16_H2overHI_err)),
                       np.log10((Saint16_H2overHI+Saint16_H2overHI_err)/Saint16_H2overHI)]
                subplot.errorbar(Saint16_mass, np.log10(Saint16_H2overHI),xerr=0.12,yerr=y_err,
                         fmt='o', markersize=5, ecolor='blue', color='blue',zorder=+3)
        
                plot_label (subplot, 'label', xlim, ylim, x_percentage=0.6, y_percentage=0.85, 
                            color='black', xlog=0, ylog=0, label='$\mathrm{H_2}/\mathrm{HI}$', 
                            fontsize=15, fontweight='normal')
            
            
                
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.savefig('./fig/HYF19_gasfractions_vs_stellarmass.pdf')
    plt.close()

    return 
#end


def H2fraction_vs_stellarmass(ThisRedshiftList):
  
    plot_color=['red','purple']        
  
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
            
        xlab='$log_{10}(M_*[M_{\odot}])$'
        ylab='$f_{\mathrm{H_2}}$'     
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   
            
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
        G0_MR=G_MR[sel]          
        G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (G0_MR['ColdGas']>0.) & (G0_MR['H2fraction']>0.)]
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
        #Fraction=np.log10(G0_MR['H2fraction'])
        Fraction=G0_MR['H2fraction']/(1.-G0_MR['H2fraction'])
 
        (x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles (bin, xlim[0], xlim[1], StellarMass, Fraction)          
        sel=(median!=0)        
        subplot.plot(x_binned[sel], np.log10(median[sel]),color=plot_color[ii], linewidth=2)     
        subplot.plot(x_binned[sel], np.log10(pc16[sel]),color=plot_color[ii], linewidth=2, linestyle='--')
        subplot.plot(x_binned[sel], np.log10(pc84[sel]),color=plot_color[ii], linewidth=2, linestyle='--')
      
        
       
            
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()

    return   
#end H2fraction_vs_stellarmass











def evo_milkyway_gas_profile(ThisRedshiftList):
  
    plot_color=['red','purple']        
   
    fig = plt.figure(figsize=(15,4))  
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
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)                       
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
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)           
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
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)        
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
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()

    return     
#end evo_milkyway_gas_profile




def evo_milkyway_stellar_profiles(ThisRedshiftList):
     
    #select galaxies at z=0 and then follow progenitors
    ii=0
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)        
    G0_MR=G_MR[sel]         
    
    
    #**********************
    #* SELECT DISC GALAXY *
    #**********************
    fig = plt.figure(figsize=(15,4))  
    grid = gridspec.GridSpec(1, 3)
    #grid.update(wspace=0.0, hspace=0.0)
    
    selected_Gal=G0_MR[(G0_MR['StellarMass']>0.) & (G0_MR['DiskMass']>0.) &  
                       (G0_MR['Vvir']/Hubble_h>200.) & (G0_MR['Vvir']/Hubble_h<235.) & (G0_MR['Type']==0) 
                       & (G0_MR['BulgeMass']/G0_MR['StellarMass']<0.15)] 
        
    print('Number of Galaxies Selected at z=0: ',len(selected_Gal))         
    selected_Gal=selected_Gal[0]              
    MainBranch=G_MR[(G_MR['GalID']>=selected_Gal['GalID']) & (G_MR['GalID']<=selected_Gal['MainLeafId'])]      
    print('First Galaxy selected ID:',selected_Gal['GalID'],'MainLeafID:',selected_Gal['MainLeafId'])  
    print('')
    print('diskmass=%0.4f stellarmass=%0.4f disk fraction=%0.4f' % 
          (np.log10(selected_Gal['DiskMass']*1.e10), np.log10(selected_Gal['StellarMass']*1.e10),
           selected_Gal['DiskMass']/selected_Gal['StellarMass'] ))  
         
    # Have a look at the colormaps here and decide which one you'd like:
    # http://matplotlib.org/1.2.1/examples/pylab_examples/show_colormaps.html
    num_plots=len(MainBranch)
    colormap = plt.cm.gist_ncar_r
    plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0.2, 1.0, num_plots)])
           
    xlim=[0.0,20.0]
    ylim=[0.,3.]
               
    #Stars Disk
    bin=2.           
    subplot=plt.subplot(grid[0])    
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)        
    xlab='$r[\mathrm{kpc}]$'
    ylab='$\Sigma_{\mathrm{DiskMass}}[M_{\odot}/\mathrm{pc^2}]$'     
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)        
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
          
    #Stars Bulge
    bin=2.           
    subplot=plt.subplot(grid[1])    
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)        
    xlab='$r[\mathrm{kpc}]$'
    ylab='$\Sigma_{\mathrm{BulgeMass}}[M_{\odot}/\mathrm{pc^2}]$'     
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)        
    colormap = plt.cm.gist_ncar_r
    plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0.2, 1.0, num_plots)])

    for jj in range (0,np.amax(MainBranch['SnapNum'])+1):        
        Gal=MainBranch[MainBranch['SnapNum']==jj]          
        if(len(Gal)>0):
            #print(Gal['SnapNum'])
            Sigma=np.zeros(RNUM,dtype=np.float32)          
            for ii in range(0,RNUM):
                Mass=Gal['BulgeMassRings'][0][ii]*1e10/Hubble_h   
                Sigma[ii]= Mass/(3.14*RingRadius[ii]*RingRadius[ii]*1e6)                     
            subplot.plot(RingRadius, np.log10(Sigma), linewidth=1)   
         
        
    #Stars Combined
    bin=2.           
    subplot=plt.subplot(grid[2])    
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)        
    xlab='$r[\mathrm{kpc}]$'
    ylab='$\Sigma_{\mathrm{stars}}[M_{\odot}/\mathrm{pc^2}]$'     
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)        
    colormap = plt.cm.gist_ncar_r
    plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0.2, 1.0, num_plots)])
              
    for jj in range (0,np.amax(MainBranch['SnapNum'])+1):        
        Gal=MainBranch[MainBranch['SnapNum']==jj]          
        if(len(Gal)>0):
            #print(Gal['SnapNum'])
            Sigma=np.zeros(RNUM,dtype=np.float32)          
            for ii in range(0,RNUM):
                Mass=(Gal['DiskMassRings'][0][ii]+Gal['BulgeMassRings'][0][ii])*1e10/Hubble_h   
                Sigma[ii]= Mass/(3.14*RingRadius[ii]*RingRadius[ii]*1e6)                     
            subplot.plot(RingRadius, np.log10(Sigma), linewidth=1)   
        
        
    plt.tight_layout()
    plt.savefig('./fig/plots_evo_milkyway_stellar_profiles_disc.pdf')
    #pdf.savefig()
    #plt.close()
    
    output.append(fig)
    
    
    
    
    #***********************
    #* SELECT BULGE GALAXY *    
    #***********************
    fig = plt.figure(figsize=(15,4))  
    grid = gridspec.GridSpec(1, 3)
    #grid.update(wspace=0.0, hspace=0.0)
    
    selected_Gal=G0_MR[(G0_MR['StellarMass']>0.) & (G0_MR['DiskMass']>0.) &  
                       (G0_MR['Vvir']/Hubble_h>200.) & (G0_MR['Vvir']/Hubble_h<235.) & (G0_MR['Type']==0) 
                       & (G0_MR['BulgeMass']/G0_MR['StellarMass']>0.95)] 
        
    print('Number of Galaxies Selected at z=0: ',len(selected_Gal))         
    selected_Gal=selected_Gal[0]              
    MainBranch=G_MR[(G_MR['GalID']>=selected_Gal['GalID']) & (G_MR['GalID']<=selected_Gal['MainLeafId'])]      
    print('First Galaxy selected ID:',selected_Gal['GalID'],'MainLeafID:',selected_Gal['MainLeafId'])  
    print('')
    print('diskmass=%0.4f stellarmass=%0.4f disk fraction=%0.4f' % 
          (np.log10(selected_Gal['DiskMass']*1.e10), np.log10(selected_Gal['StellarMass']*1.e10),
           selected_Gal['DiskMass']/selected_Gal['StellarMass'] ))  
         
    # Have a look at the colormaps here and decide which one you'd like:
    # http://matplotlib.org/1.2.1/examples/pylab_examples/show_colormaps.html
    num_plots=len(MainBranch)
    colormap = plt.cm.gist_ncar_r
    plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0.2, 1.0, num_plots)])
           
    xlim=[0.0,20.0]
    ylim=[0.,3.]
               
    #Stars Disk
    bin=2.           
    subplot=plt.subplot(grid[0])    
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)        
    xlab='$r[\mathrm{kpc}]$'
    ylab='$\Sigma_{\mathrm{DiskMass}}[M_{\odot}/\mathrm{pc^2}]$'     
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)        
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
          
    #Stars Bulge
    bin=2.           
    subplot=plt.subplot(grid[1])    
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)        
    xlab='$r[\mathrm{kpc}]$'
    ylab='$\Sigma_{\mathrm{BulgeMass}}[M_{\odot}/\mathrm{pc^2}]$'     
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)        
    colormap = plt.cm.gist_ncar_r
    plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0.2, 1.0, num_plots)])

    for jj in range (0,np.amax(MainBranch['SnapNum'])+1):        
        Gal=MainBranch[MainBranch['SnapNum']==jj]          
        if(len(Gal)>0):
            #print(Gal['SnapNum'])
            Sigma=np.zeros(RNUM,dtype=np.float32)          
            for ii in range(0,RNUM):
                Mass=Gal['BulgeMassRings'][0][ii]*1e10/Hubble_h   
                Sigma[ii]= Mass/(3.14*RingRadius[ii]*RingRadius[ii]*1e6)                     
            subplot.plot(RingRadius, np.log10(Sigma), linewidth=1)   
         
        
    #Stars Combined
    bin=2.           
    subplot=plt.subplot(grid[2])    
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)        
    xlab='$r[\mathrm{kpc}]$'
    ylab='$\Sigma_{\mathrm{stars}}[M_{\odot}/\mathrm{pc^2}]$'     
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)        
    colormap = plt.cm.gist_ncar_r
    plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0.2, 1.0, num_plots)])
              
    for jj in range (0,np.amax(MainBranch['SnapNum'])+1):        
        Gal=MainBranch[MainBranch['SnapNum']==jj]          
        if(len(Gal)>0):
            #print(Gal['SnapNum'])
            Sigma=np.zeros(RNUM,dtype=np.float32)          
            for ii in range(0,RNUM):
                Mass=(Gal['DiskMassRings'][0][ii]+Gal['BulgeMassRings'][0][ii])*1e10/Hubble_h   
                Sigma[ii]= Mass/(3.14*RingRadius[ii]*RingRadius[ii]*1e6)                     
            subplot.plot(RingRadius, np.log10(Sigma), linewidth=1)   
        
        
    
    plt.tight_layout()
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()

    return  
#end evo_milkyway_stellar_profiles

def test_H2_prescriptions(ThisRedshiftList):
  
    for ii in range(0,len(ThisRedshiftList)):        
                   
        #HII
        xlim=[0.0,4.0]        
        ylim=[-2.5,0.2]
                
        plot_color=['red','purple']        
     
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
            #Metallicity[ii*len(G0_MR):(ii+1)*len(G0_MR)]=np.log10(G0_MR['MetalsColdGasRings'][:,ii]/
            #                                                      G0_MR['ColdGasRings'][:,ii]/0.02)
            StellarDensity[ii*len(G0_MR):(ii+1)*len(G0_MR)]=np.log10((G0_MR['DiskMassRings'][:,ii]*1e10/Hubble_h)/area)
        #SigmaGas=np.log10((G0_MR['ColdGas']*1e10/Hubble_h)/(3.14*G0_MR['GasDiskRadius']*G0_MR['GasDiskRadius']*1e6*1e6))       
        #Fraction=np.log10(G0_MR['H2fraction'])
        
        #Bin by metallicity
        grid = gridspec.GridSpec(1, 2)
        subplot=plt.subplot(grid[0])
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)                    
        xlab='$log_{10}(\Sigma_{\mathrm{gas}}[M_{\odot}\mathrm{pc}^{-2}])$'
        ylab='$log_{10}(f_{\mathrm{H_2}})$'  
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   
        
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
                (x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles (bin, xlim[0], xlim[1], 
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
                (x_binned, median, mean, pc16, pc84, rms)=median_and_percentiles (bin, xlim[0], xlim[1], 
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
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()

    return   
#end test_H2_prescriptions









def test_rings(ThisRedshiftList):
    
    ii=0   
        
    plot_color=['blue','green','red']        
    
    (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii, FullSnapshotList_MR)                 
    G_MR=G_MR[sel]    
    N=min(len(G_MR),1000)
    G0_MR=np.random.choice(G_MR, size=N, replace=False)  
    
    fig = plt.figure(figsize=(15,10))
    grid = gridspec.GridSpec(2, 3)
   

    #********************
    #*        MASS   
    #********************
    
    xlim=[6.0,11.0]
    ylim=[6.0,11.0]        
    #Gas Mass (Total vs Rings)
    subplot=plt.subplot(grid[0])    
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)      
    xlab='$\mathrm{Total} - M_{\mathrm{Cold}}$'           
    ylab='$\mathrm{Rings Sum} - M_{\mathrm{Cold}}$'               
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)     
    MassRings=np.sum(G0_MR['ColdGasRings'],axis=1)    
    subplot.scatter(np.log10(G0_MR['ColdGas']*1e10),np.log10(MassRings*1e10),s=1, color='black')  
    
    #Disk Mass (Total vs Rings)   
    subplot=plt.subplot(grid[1])   
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)      
    xlab='$\mathrm{Total} - M_{\mathrm{DiskMass}}$'           
    ylab='$\mathrm{Rings Sum} - M_{\mathrm{DiskMass}}$'               
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   
    MassRings=np.sum(G0_MR['DiskMassRings'],axis=1)           
    subplot.scatter(np.log10(G0_MR['DiskMass']*1e10),np.log10(MassRings*1e10),s=5, color='black')  
          
    #Bulge Mass (Total vs Rings)   
    subplot=plt.subplot(grid[2])           
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)      
    xlab='$\mathrm{Total} - M_{\mathrm{BulgeMass}}$'           
    ylab='$\mathrm{Rings Sum} - M_{\mathrm{BulgeMass}}$'               
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)   
    MassRings=np.sum(G0_MR['BulgeMassRings'],axis=1)           
    subplot.scatter(np.log10(G0_MR['BulgeMass']*1e10),np.log10(MassRings*1e10),s=5, color='black')  
    
    
    
    #********************
    #*        METALS   
    #********************
    
    xlim=[5.0,9.0]
    ylim=[5.0,9.0]  
    
    #Mass of Metals in Gas (Total vs Rings)  
    subplot=plt.subplot(grid[3])             
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)      
    xlab='$\mathrm{Total} - M_{\mathrm{MetalsCold}}$'           
    ylab='$\mathrm{Rings Sum} - M_{\mathrm{MetalsCold}}$'               
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)  
  
    if(opt_detailed_enrichment==1):
        Metals=G0_MR['MetalsColdGas'][:,0]+G0_MR['MetalsColdGas'][:,1]+G0_MR['MetalsColdGas'][:,2]      
        MetalsRings=G0_MR['MetalsColdGasRings'][:,:,0] + G0_MR['MetalsColdGasRings'][:,:,1] + \
        G0_MR['MetalsColdGasRings'][:,:,2]        
        MetalsRings=np.sum(MetalsRings,axis=1)        
    else:    
        Metals=G0_MR['MetalsColdGas']
        MetalsRings=np.sum(G0_MR['MetalsColdGasRings'],axis=1)  
    subplot.scatter(np.log10(Metals*1e10), np.log10(MetalsRings*1e10),s=5, color='black')  

  
    #DiskMass Metallicity (Total vs Rings)    
    subplot=plt.subplot(grid[4])         
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)      
    xlab='$\mathrm{Total} - M_{\mathrm{MetalsDiskMass}}$'           
    ylab='$\mathrm{Rings Sum} - M_{\mathrm{MetalsDiskMass}}$'               
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)  
        
    if(opt_detailed_enrichment==1):
        Metals=G0_MR['MetalsDiskMass'][:,0]+G0_MR['MetalsDiskMass'][:,1]+G0_MR['MetalsDiskMass'][:,2]      
        MetalsRings=G0_MR['MetalsDiskMassRings'][:,:,0] + G0_MR['MetalsDiskMassRings'][:,:,1] + \
                           G0_MR['MetalsDiskMassRings'][:,:,2]       
        MetalsRings=np.sum(MetalsRings,axis=1)
    else:
        Metals=G0_MR['MetalsDiskMass']        
        MetalsRings=np.sum(G0_MR['MetalsDiskMassRings'],axis=1)        
    subplot.scatter(np.log10(Metals*1e10), np.log10(MetalsRings*1e10),s=5, color='black')  

    
    #BulgeMass Metallicity (Total vs Rings)  
    subplot=plt.subplot(grid[5])  
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)      
    xlab='$\mathrm{Total} - M_{\mathrm{MetalsBulgeMass}}$'           
    ylab='$\mathrm{Rings Sum} - M_{\mathrm{MetalsBulgeMass}}$'               
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)  

    if(opt_detailed_enrichment==1):
        Metals=G0_MR['MetalsBulgeMass'][:,0]+G0_MR['MetalsBulgeMass'][:,1]+G0_MR['MetalsBulgeMass'][:,2]      
        MetalsRings=G0_MR['MetalsBulgeMassRings'][:,:,0] + G0_MR['MetalsBulgeMassRings'][:,:,1] + \
                           G0_MR['MetalsBulgeMassRings'][:,:,2]      
        MetalsRings=np.sum(MetalsRings,axis=1)
    else:
        Metals=G0_MR['MetalsBulgeMass']        
        MetalsRings=np.sum(G0_MR['MetalsBulgeMassRings'],axis=1)        
    subplot.scatter(np.log10(Metals*1e10), np.log10(MetalsRings*1e10),s=5, color='black')  

       
    #pdf.savefig()
    #plt.close()       
  
    
    
    
    
    
    
    
    
    #**************************************************
    #*  Sum Mass in Elements vs Sum Mass in Rings  
    #**************************************************
    
    
    fig = plt.figure(figsize=(15,10))
    grid = gridspec.GridSpec(2, 3)
    
    xlim=[7.0,10.0]
    ylim=[7.0, 10.0]
        
    #Cold Gas elements vs cold gas(Total vs Rings)
    if(opt_detailed_enrichment==1):       
        subplot=plt.subplot(grid[0])               
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)      
        xlab='Total - Sum Elements in ColdGas'           
        ylab='Rings Sum - ColdGas'               
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)      
        
        ElementsSum=np.sum(G0_MR['ColdGas_elements'],axis=1)    
        MassRings=np.sum(G0_MR['ColdGasRings']*1e10/Hubble_h ,axis=1)        
        subplot.scatter(np.log10(ElementsSum), np.log10(MassRings),s=5, color='red')  
      
    #Disk elements vs DiskMass(Total vs Rings)
    if(opt_detailed_enrichment==1):       
        subplot=plt.subplot(grid[1])       
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)      
        xlab='Total - Sum Elements in DiskMass'           
        ylab='Rings Sum - DiskMass'               
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)      
        
        ElementsSum=np.sum(G0_MR['DiskMass_elements'],axis=1)    
        MassRings=np.sum(G0_MR['DiskMassRings']*1e10/Hubble_h ,axis=1)        
        subplot.scatter(np.log10(ElementsSum), np.log10(MassRings),s=5, color='red')  
        
    #Bulge elements vs BulgeMass(Total vs Rings)
    if(opt_detailed_enrichment==1):       
        subplot=plt.subplot(grid[2])       
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)      
        xlab='Total - Sum Elements in BulgeMass'           
        ylab='Rings Sum - BulgeMass'               
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)      
        
        ElementsSum=np.sum(G0_MR['BulgeMass_elements'],axis=1)    
        MassRings=np.sum(G0_MR['BulgeMassRings']*1e10/Hubble_h ,axis=1)        
        subplot.scatter(np.log10(ElementsSum), np.log10(MassRings),s=5, color='red')     
        
        
    #**************************************************
    #*  Mass in Element vs Element in Rings  
    #**************************************************    
    
    xlim=[4.0,8.0]
    ylim=[4.0, 8.0]        
    #Cold Gas Elements(Total vs Rings)
    if(opt_detailed_enrichment==1):
        subplot=plt.subplot(grid[3])        
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)      
        xlab='Total - ColdGas in Element Nr6'           
        ylab='Rings Sum - ColdGas in Element Nr6'               
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)    
      
        Elements=G0_MR['ColdGas_elements'][:,6]      
        ElementsRings=np.sum(G0_MR['ColdGasRings_elements'][:,:,6],axis=1)
        subplot.scatter(np.log10(Elements), np.log10(ElementsRings),s=5, color='red')  
              
    #DiskMass Elements(Total vs Rings)
    if(opt_detailed_enrichment==1):
        subplot=plt.subplot(grid[4])        
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)      
        xlab='Total - DiskMass in Element Nr6'           
        ylab='Rings Sum - DiskMass in Element Nr6'               
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)    
      
        Elements=G0_MR['DiskMass_elements'][:,6]      
        ElementsRings=np.sum(G0_MR['DiskMassRings_elements'][:,:,6],axis=1)
        subplot.scatter(np.log10(Elements), np.log10(ElementsRings),s=5, color='red')  
    
    #BulgeMass Elements(Total vs Rings)
    if(opt_detailed_enrichment==1):
        subplot=plt.subplot(grid[5])        
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)      
        xlab='Total - BulgeMass in Element Nr6'           
        ylab='Rings Sum - BulgeMass in Element Nr6'               
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)    
      
        Elements=G0_MR['BulgeMass_elements'][:,6]      
        ElementsRings=np.sum(G0_MR['BulgeMassRings_elements'][:,:,6],axis=1)
        subplot.scatter(np.log10(Elements), np.log10(ElementsRings),s=5, color='red') 
    
        plt.tight_layout()        
        #pdf.savefig()
        #plt.close()           
       
    #***************************
    #*     TEST SFH ARRAYS     *
    #***************************
    
    fig = plt.figure(figsize=(15,10))
    grid = gridspec.GridSpec(2, 3)
    
    xlim=[7.0,11.0]
    ylim=[7.0,11.0]   
    x=np.arange(xlim[0],xlim[1]+1.0,0.1)   
    
    
    #******************
    #*       DISC     *
    #******************
    #Disk Mass (Total vs SFH)
    subplot=plt.subplot(grid[0])         
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)      
    xlab='$\mathrm{Total} - M_{\mathrm{DiskMass}}$'           
    ylab='$\mathrm{SFH\;Sum} - M_{\mathrm{DiskMass}}$'               
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)  
    
    sfh_Mass=np.sum(G0_MR['sfh_DiskMass'],axis=1)           
    subplot.scatter(np.log10(G0_MR['DiskMass']*1e10),np.log10(sfh_Mass*1e10*(1-0.43)),s=5, color='black') 
    #subplot.scatter(np.log10(G0_MR['DiskMass']*1e10),np.log10(sfh_Mass*1e10),s=5, color='black') 
    subplot.plot(x,x) 
    
    
    #Disk Mass RING Nr 2 (Total vs SFH)
    subplot=plt.subplot(grid[1])        
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)      
    xlab='$\mathrm{Total} - M_{\mathrm{DiskMass}} \mathrm{Ring[2]}$'           
    ylab='$\mathrm{SFH\;Sum} - M_{\mathrm{DiskMass}} \mathrm{Ring[2]}$'               
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14) 
    
    sfh_MassRing=np.sum(G0_MR['sfh_DiskMassRings'][:,2,:],axis=1)           
    subplot.scatter(np.log10(G0_MR['DiskMassRings'][:,2]*1e10),np.log10(sfh_MassRing*1e10*(1-0.43)),s=5, color='black')  
    #subplot.scatter(np.log10(G0_MR['DiskMassRings'][:,2]*1e10),np.log10(sfh_MassRing*1e10),s=5, color='black')  
    subplot.plot(x,x) 
   
    
    #Disk Mass (Total vs SFHRings)
    subplot=plt.subplot(grid[2])       
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)      
    xlab='$\mathrm{Total} - M_{\mathrm{DiskMass}}$'           
    ylab='$\mathrm{SFH\;Sum\;All Rings} - M_{\mathrm{DiskMass}}$'               
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)  
   
    sfh_MassRings=np.zeros(len(G0_MR))
    for jj in range(0,RNUM):
        sfh_MassThisRing=G0_MR['sfh_DiskMassRings'][:,jj,:]    
        sfh_MassRings+=np.sum(sfh_MassThisRing,axis=1)            
    subplot.scatter(np.log10(np.sum(G0_MR['sfh_DiskMass'],axis=1)*1e10),np.log10(sfh_MassRings*1e10),s=5, color='black')
    subplot.plot(x,x)
    
    
    #******************
    #*       BULGE     *
    #******************
    #BulgeMass (Total vs SFH)
    subplot=plt.subplot(grid[3])         
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)      
    xlab='$\mathrm{Total} - M_{\mathrm{BulgeMass}}$'           
    ylab='$\mathrm{SFH\;Sum} - M_{\mathrm{BulgeMass}}$'               
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)  
    
    sfh_Mass=np.sum(G0_MR['sfh_BulgeMass'],axis=1)           
    subplot.scatter(np.log10(G0_MR['BulgeMass']*1e10),np.log10(sfh_Mass*1e10*(1-0.43)),s=5, color='black')  
    #subplot.scatter(np.log10(G0_MR['BulgeMass']*1e10),np.log10(sfh_Mass*1e10),s=5, color='black')  
    subplot.plot(x,x) 
    
    #Bulge Mass RING Nr 2 (Total vs SFH)
    subplot=plt.subplot(grid[4])        
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)      
    xlab='$\mathrm{Total} - M_{\mathrm{BulgeMass}} \mathrm{Ring[2]}$'           
    ylab='$\mathrm{SFH\;Sum} - M_{\mathrm{BulgeMass}} \mathrm{Ring[2]}$'               
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14) 
    
    sfh_MassRing=np.sum(G0_MR['sfh_BulgeMassRings'][:,0,:],axis=1)           
    subplot.scatter(np.log10(G0_MR['BulgeMassRings'][:,0]*1e10),np.log10(sfh_MassRing*1e10*(1-0.43)),s=5, color='black') 
    #subplot.scatter(np.log10(G0_MR['BulgeMassRings'][:,0]*1e10),np.log10(sfh_MassRing*1e10),s=5, color='black')  
    subplot.plot(x,x) 
   
    
    #Bulge Mass (Total vs SFHRings)
    subplot=plt.subplot(grid[5])       
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)      
    xlab='$\mathrm{Total SFH} - M_{\mathrm{BulgeMass}}$'           
    ylab='$\mathrm{SFH\;Sum\;All Rings} - M_{\mathrm{BulgeMass}}$'               
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)  
   
    sfh_MassRings=np.zeros(len(G0_MR))
    for jj in range(0,RNUM):
        sfh_MassThisRing=G0_MR['sfh_BulgeMassRings'][:,jj,:]    
        sfh_MassRings+=np.sum(sfh_MassThisRing,axis=1)            
    subplot.scatter(np.log10(np.sum(G0_MR['sfh_BulgeMass'],axis=1)*1e10),np.log10(sfh_MassRings*1e10),s=5, color='black')
    subplot.plot(x,x) 
    
    plt.tight_layout()       
    #pdf.savefig()
    #plt.close()    
   
    
    
    #***************************
    #*  TEST SFH METAL ARRAYS  *
    #***************************
    
    fig = plt.figure(figsize=(15,10))
    grid = gridspec.GridSpec(2, 3)
    
    xlim=[0.0,11.0]
    ylim=[0.0,11.0]   
    x=np.arange(xlim[0],xlim[1]+1.0,0.1)   
    
    
    #******************
    #*       DISC     *
    #******************
    #Metals Disk Mass (Total vs SFH)
    subplot=plt.subplot(grid[0])         
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)      
    xlab='$\mathrm{Total} - M_{\mathrm{MetalsDiskMass}}$'           
    ylab='$\mathrm{SFH\;Sum} - M_{\mathrm{MetalsDiskMass}}$'               
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)  
    
    if(opt_detailed_enrichment==1):
        sfh_metals=np.sum(G0_MR['sfh_MetalsDiskMass'],axis=2)
        metals=np.sum(G0_MR['MetalsDiskMass'],axis=1)
        #sfh_metals=np.sum(G0_MR['sfh_ElementsDiskMass'],axis=2)
        #metals=np.sum(G0_MR['DiskMass_elements'],axis=1)
    else:
        sfh_metals=G0_MR['sfh_MetalsDiskMass']
        metals=G0_MR['MetalsDiskMass']
    sfh_metals=np.sum(sfh_metals,axis=1)      
    subplot.scatter(np.log10(metals*1e10),np.log10(sfh_metals*1e10*(1-0.43)),s=5, color='black') 
    #subplot.scatter(np.log10(metals),np.log10(sfh_metals),s=5, color='black') 
    subplot.plot(x,x) 
    
    
    #******************
    #*       BULGE     *
    #******************
    #Metals BulgeMass (Total vs SFH)
    subplot=plt.subplot(grid[3])         
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)      
    xlab='$\mathrm{Total} - M_{\mathrm{MetalsBulgeMass}}$'           
    ylab='$\mathrm{SFH\;Sum} - M_{\mathrm{MetalsBulgeMass}}$'               
    subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)  
    
    if(opt_detailed_enrichment==1):
        sfh_metals=np.sum(G0_MR['sfh_MetalsBulgeMass'],axis=2)
        metals=np.sum(G0_MR['MetalsBulgeMass'],axis=1)
    else:
        sfh_metals=G0_MR['sfh_MetalsBulgeMass']
        metals=G0_MR['MetalsBulgeMass']
    sfh_metals=np.sum(sfh_metals,axis=1)       
    subplot.scatter(np.log10(metals*1e10),np.log10(sfh_metals*1e10*(1-0.43)),s=5, color='black')  
    subplot.plot(x,x) 
    
    
    
    if(opt_detailed_enrichment==1): 
        #******************
        #*       DISC     *
        #******************
        #Elements Disk Mass (Total vs SFH)
        subplot=plt.subplot(grid[1])         
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)      
        xlab='$\mathrm{Total} - M_{\mathrm{ElementsDiskMass}}$'           
        ylab='$\mathrm{SFH\;Sum} - M_{\mathrm{ElementsDiskMass}}$'               
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)  
    
        sfh_metals=np.sum(G0_MR['sfh_ElementsDiskMass'],axis=2)
        metals=np.sum(G0_MR['DiskMass_elements'],axis=1)  
        sfh_metals=np.sum(sfh_metals,axis=1)      
        subplot.scatter(np.log10(metals),np.log10(sfh_metals*(1-0.43)),s=5, color='black') 
        subplot.plot(x,x) 
    
    
        #******************
        #*       BULGE     *
        #******************
        #Elements BulgeMass (Total vs SFH)
        subplot=plt.subplot(grid[4])         
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)      
        xlab='$\mathrm{Total} - M_{\mathrm{ElementsBulgeMass}}$'           
        ylab='$\mathrm{SFH\;Sum} - M_{\mathrm{ElementsBulgeMass}}$'               
        subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)  
    
        sfh_metals=np.sum(G0_MR['sfh_ElementsBulgeMass'],axis=2)
        metals=np.sum(G0_MR['BulgeMass_elements'],axis=1)   
        sfh_metals=np.sum(sfh_metals,axis=1)       
        subplot.scatter(np.log10(metals),np.log10(sfh_metals*(1-0.43)),s=5, color='black')  
        subplot.plot(x,x) 
    
    #print(sfh_Mass)
    plt.tight_layout()       
    current_function =  inspect.getframeinfo(inspect.currentframe()).function   
    plt.savefig('./fig/plots_'+current_function+'.pdf')
    plt.close()    
       
    return       
#end test_rings        