# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

"""
read_snap
read_tree
redshift_to_time
select_current_redshift
stellar_mass_with_err
median_and_percentiles
median_and_percentiles_fixed_xx
grayify_cmap
plot_label_three_models
plot_label
plot_joint_MR_MRII
mag_to_lum
get_slope
smooth
convert_2d_array_into_region
abs_to_app_mag
comdist
plot_mass_function
"""
    
    
          
""" 
self.capacity *= 4
newdata = np.zeros((self.capacity,))
newdata[:self.size] = self.data
self.data = newdata
self.data[self.size] = x
     
data = self.data[:self.size]
return np.reshape(data, newshape=(len(data)/5, 5))"""


import numpy as np
from plots_input import *
import matplotlib.pyplot as plt

def read_snap(folder,FirstFile,LastFile,
              props,template,RedshiftsToRead,FullRedshiftList):    
    """ Reads L-Galaxy output files.
    Returns: (nTrees,nHalos,nTreeHalos,gals)
    Inputs: (folder,file_prefix,FirstFile,LastFile,props,template)
    props - list of properties to return
    template - structure dtype definition from database """     
    nTrees = 0
    nGals = 0    
    nTreeHalos = np.array([],dtype=np.int32)
    
    filter_list = []
    for prop in props:
        if props[prop]:
            filter_list.append((prop,template[prop]))
    filter_dtype = np.dtype(filter_list)  
            
    SnapshotList=np.array([],dtype=np.int32)
    
    #read only headers to figure out total nGals
    print ("\n\nReading Headers\n")
    for iredshift in range(0,len(FullRedshiftList)):
        if RedshiftsToRead[iredshift]:              
            for ifile in range(FirstFile,LastFile+1):
                char_redshift="%0.2f" % FullRedshiftList[iredshift]
                filename = folder+'/'+'SA_z'+char_redshift+"_"+"%d"%(ifile)               
                f = open(filename,"rb")
                
                this_nTrees =  np.fromfile(f,np.int32,1)
                nTrees += this_nTrees
                this_nGals = np.fromfile(f,np.int32,1)
                nGals += this_nGals
                f.close()
            print ("z=", char_redshift," nGals = ",nGals)  
               
    gals = np.zeros(nGals,dtype=filter_dtype)
  
    print ("\n")
    offset=0
    for iredshift in range(0,len(FullRedshiftList)):
        if RedshiftsToRead[iredshift]:  
            print ("\nReading redshift: ", FullRedshiftList[iredshift], "\n")
            for ifile in range(FirstFile,LastFile+1):
                char_redshift="%0.2f" % FullRedshiftList[iredshift]
                filename = folder+'/'+'SA_z'+char_redshift+"_"+"%d"%(ifile)
                #print(filename)
                f = open(filename,"rb")
                
                this_nTrees =  np.fromfile(f,np.int32,1)
                nTrees += this_nTrees
                this_nGals = np.fromfile(f,np.int32,1)
                nGals += this_nGals
                print ("File ", ifile," nGals = ",this_nGals)  
                
                addednTreeHalos = np.fromfile(f,np.int32,int(this_nTrees))
                nTreeHalos = np.append(nTreeHalos,addednTreeHalos)
                full_this_gals = np.fromfile(f,template,int(this_nGals)) # all properties
                this_gals = np.zeros(this_nGals,dtype=filter_dtype) # selected props
                
                for prop in template.names:
                    if props[prop]:
                        this_gals[prop] = full_this_gals[prop]
                              
                gals[offset:offset+int(this_nGals)] = this_gals[:]    
                offset+=int(this_nGals)
                f.close()       
                
            #endfor
        #endif    
        #assign snapshot of current redshift given by the last galaxy on the last file 
        SnapshotList=np.append(SnapshotList,gals['SnapNum'][offset-1])       
    #endfor
    
    return (gals, SnapshotList)




def read_tree(folder,FirstFile,LastFile,
              props,template):    
    """ Reads L-Galaxy output files.
    Returns: (nTrees,nHalos,nTreeHalos,gals)
    Inputs: (folder,file_prefix,FirstFile,LastFile,props,template)
    props - list of properties to return
    template - structure dtype definition from database """   
    nGals = 0    
    
    filter_list = []
    for prop in props:
        if props[prop]:
            filter_list.append((prop,template[prop]))
    filter_dtype = np.dtype(filter_list)  
            
    #SnapshotList=np.array([],dtype=np.int32)
    
    #read only headers to figure out total nGals
    print ("\n\nReading Headers\n")
    for ifile in range(FirstFile,LastFile+1):       
        filename = folder+'/'+'SA_galtree_'+"%d"%(ifile)               
        f = open(filename,"rb")
        one = np.fromfile(f,np.int32,1)       
        nbytes = np.fromfile(f,np.int32,1)
        this_nGals = np.fromfile(f,np.int32,1)
        #print("TotNgals=",this_nGals)
        #this_nGals = 35000000
        #this_nGals = 10000000
        nGals += this_nGals  
        f.close()
    gals = np.zeros(nGals,dtype=filter_dtype)
    
    print("TotNgals=",nGals)
    print ("\n")
    
    offset=0
    for ifile in range(FirstFile,LastFile+1):         
        filename = folder+'/'+'SA_galtree_'+"%d"%(ifile)             
        f = open(filename,"rb")
        one = np.fromfile(f,np.int32,1)
        nbytes = np.fromfile(f,np.int32,1)         
        nskip=nbytes/4-3
        this_nGals = int(np.fromfile(f,np.int32,1))  
        #print(one, nbytes,this_nGals)
        #this_nGals =35000000
        #this_nGals =3000
        #this_nGals = 10000000
        nGals += this_nGals       
        print ("File ", ifile," nGals = ",this_nGals)           
        ib=np.fromfile(f,np.float32,int(nskip))       
       
        full_this_gals = np.fromfile(f,template,this_nGals) # all properties 
        #full_this_gals = np.fromfile(f,template,this_nGals) # all properties        
        #full_this_gals = np.fromfile(f,template,this_nGals) # all properties 
             
        this_gals = np.zeros(this_nGals,dtype=filter_dtype) # selected props
               
       
        for prop in template.names:
            if props[prop]:
                this_gals[prop] = full_this_gals[prop]
                              
        gals[offset:offset+this_nGals] = this_gals[:]    
        offset+=this_nGals
        #f.close()          
    #endfor
   
    return (gals)

def redshift_to_time (z):
    Tyr = 977.8    ;# coefficent for converting 1/H into Gyr                               
    WM = 0.315
    WV = 0.685
    H0=67.3
    h = H0/100.
    WR = 4.165E-5/(h*h)   ;# includes 3 massless neutrino species, T0 = 2.72528            
    WK = 1-WM-WR-WV
    az = 1.0/(1+1.0*z)
    age = 0.
    n=1000        ; # number of points in integrals                                        
    a=0
    for i in range(0, n-1):
        a = az*(i+0.5)/n
        adot = np.sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
        age = age + 1./adot
    
    zage = az*age/n
    age_Gyr = (Tyr/H0)*zage
    
    return (age_Gyr)
#end redshift_to_time

def select_current_redshift(G_MR, ThisRedshiftList, ii, SnapshotList):
    
    found_redshift=0
                           
    for jj in range(0, len(FullRedshiftList)):           
        if(ThisRedshiftList[ii]<1.):
            if round(FullRedshiftList[jj],1)==round(ThisRedshiftList[ii],1):                
                sel= (G_MR['SnapNum']==SnapshotList[jj])
                found_redshift=1                  
        else:    
            if round(FullRedshiftList[jj],0)==round(ThisRedshiftList[ii],0):               
                sel= (G_MR['SnapNum']==SnapshotList[jj])
                found_redshift=1 
                    
    if found_redshift==0:
        sys.exit("redshift:",ThisRedshiftList[ii],"needed for stellar mass function not read.") 
        
    return (sel)

#end  select_current_redshift


def stellar_mass_with_err(G0_MR, Hubble_h, redshift):
    
    np.random.seed(seed=10)
    mass= np.log10(G0_MR['StellarMass']*1.e10/Hubble_h) + np.random.randn(len(G0_MR['StellarMass']))*0.08*(1+redshift)    
    #mass= np.log10(G0_MR['StellarMass']*1.e10*Hubble_h) #+ np.random.randn(len(G0_MR['StellarMass']))*0.08*(1+redshift)

    return mass

#end get_stellar_mass

def sfr_with_err(G0_MR, Hubble_h, redshift):
    
    np.random.seed(seed=10)
    sfr= np.log10(G0_MR['Sfr']) + np.random.randn(len(G0_MR['Sfr']))*0.08*(1+redshift)    
    sfr= np.log10(G0_MR['Sfr']) + np.random.randn(len(G0_MR['Sfr']))*0.4
    
    return sfr

#end get_stellar_mass


def median_and_percentiles (bin, xmin, xmax, x_variable, y_variable): 
  
    min_x=min(x_variable[x_variable > (-1e20)])
    
    if(min_x > xmin):
        xmin=min_x
    #print(xmin,xmax)
    Nbins=int((xmax-xmin)/bin+1)   
   
    median=np.zeros(Nbins, np.float32)
    mean=np.zeros(Nbins, np.float32)
    pc16=np.zeros(Nbins, np.float32) 
    pc84=np.zeros(Nbins, np.float32)  
    rms=np.zeros(Nbins, np.float32) 
    x_min=xmin-bin/2.
  
    for ii in range(0,Nbins): 

        x_variable_sel=x_variable[(x_variable > (x_min+(ii*bin))) & (x_variable < (x_min+(ii+1.0)*bin))] 
        y_variable_sel=y_variable[(x_variable > (x_min+(ii*bin))) & (x_variable < (x_min+(ii+1.0)*bin))] 
       
        if(len(x_variable_sel) > 0):           
            median[ii]=np.median(y_variable_sel) 
            mean[ii]=np.mean(y_variable_sel) 
            y_sorted = np.sort(y_variable_sel)
            pc16[ii] = y_sorted[int(16*len(y_variable_sel)/100)]      
            pc84[ii] = y_sorted[int(84*len(y_variable_sel)/100)]            
            #rms[ii]=np.sqrt(np.mean((np.log10(y_variable_sel))**2))
            rms[ii]=np.sqrt(1./len(y_variable_sel)*np.sum((mean[ii]-y_variable_sel)**2))
    #endfor

    x_binned=np.arange(Nbins)*((x_min+(Nbins*bin))-(x_min+(Nbins*0.0)))/(Nbins*1.)+x_min+bin/2.
 
    return (x_binned, median, mean, pc16, pc84, rms)

#end median_and_percentiles
    

def median_and_percentiles_fixed_xx (x_variable, y_variable, non_zero=1): 
  
    x_binned = np.unique(x_variable) 
    if(non_zero):
        x_binned=x_binned[x_binned>0.]
    Nbins = len(x_binned)
   
    median=np.zeros(Nbins, np.float32)
    mean=np.zeros(Nbins, np.float32)
    pc16=np.zeros(Nbins, np.float32) 
    pc84=np.zeros(Nbins, np.float32)  
    rms=np.zeros(Nbins, np.float32) 
 
    for idx, element in enumerate(x_binned):          
        #sel = (x_variable == element) & (y_variable>0.)
        sel = (x_variable == element) 
        x_variable_sel=x_variable[sel] 
        y_variable_sel=y_variable[sel] 
     
        if(len(x_variable_sel) > 0):           
            median[idx]=np.median(y_variable_sel) 
            mean[idx]=np.mean(y_variable_sel) 
            y_sorted = np.sort(y_variable_sel)
            pc16[idx] = y_sorted[int(16*len(y_variable_sel)/100)]      
            pc84[idx] = y_sorted[int(84*len(y_variable_sel)/100)]            
            #rms[ii]=np.sqrt(np.mean((np.log10(y_variable_sel))**2))
            rms[idx]=np.sqrt(1./len(y_variable_sel)*np.sum((mean[idx]-y_variable_sel)**2))
    #endfor

    return (x_binned, median, mean, pc16, pc84, rms)

#end median_and_percentiles

def plot_fraction_vmax_weighted (xmin,xmax,bin,property,mag):
        
    Nbins=(xmax-xmin)/bin+1.      
    frac=np.zeros(int(Nbins))
         
    max_d=10**((17.6-mag)/5.+1.)/1.e6 
    weight = 1./max_d**3
    
    nii=np.sum(weight)
    ind=0
      
    for jj in range (0,int(Nbins)):
        cid=xmin+bin*jj
        ww=((property > cid) & (property < cid+bin))
        if len(ww)==0:
            ind+=1
            continue
          
        frac[ind]=np.sum(weight[ww])/nii
        ind+=1
        
    x_arr=np.arange(xmin,xmax+bin,bin)        
    
    return (x_arr,frac)

#end plot_fraction_vmax_weighted

def grayify_cmap(cmap):
    """Return a grayscale version of the colormap"""
    cmap = plt.cm.get_cmap(cmap)
    colors = cmap(np.arange(cmap.N))
    
    # convert RGBA to perceived greyscale luminance
    # cf. http://alienryderflex.com/hsp.html
    RGB_weight = [0.299, 0.587, 0.114]
    luminance = np.sqrt(np.dot(colors[:, :3] ** 2, RGB_weight))
    colors[:, :3] = luminance[:, np.newaxis]
    
    return cmap.from_list(cmap.name + "_grayscale", colors, cmap.N)

#end grayify_cmap


def plot_label_three_models (subplot, xlim, ylim, position):
 
    if position=='top_left':
        x1=0.15        
        x21=0.04
        x22=0.13
        
        previous_model1_y1=0.9
        previous_model1_y2=0.92       
        previous_model2_y1=0.83
        previous_model2_y2=0.85        
        this_model_y1=0.76
        this_model_y2=0.78
        
    if position=='bottom_left':
        x1=0.15        
        x21=0.04
        x22=0.13
        
        previous_model1_y1=0.2
        previous_model1_y2=0.22       
        previous_model2_y1=0.13
        previous_model2_y2=0.15        
        this_model_y1=0.06
        this_model_y2=0.08
        
    plot_label (subplot, 'label', xlim, ylim, 
                x_percentage=x1, y_percentage=previous_model1_y1, color='black', xlog=0, ylog=0, 
                label=prefix_previous_model1, fontsize=10, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim,
                x_percentage=x21, y_percentage=previous_model1_y2, color='red', x2_percentage=x22, 
                xlog=0, ylog=0, linestyle=linestyle_previous_model1, linewidth=2)  
                    
    plot_label (subplot, 'label', xlim, ylim, 
                x_percentage=x1, y_percentage=previous_model2_y1, color='black', xlog=0, ylog=0, 
                label=prefix_previous_model2, fontsize=10, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim,
                x_percentage=x21, y_percentage=previous_model2_y2, color='red', x2_percentage=x22, 
                xlog=0, ylog=0, linestyle=linestyle_previous_model2, linewidth=2)
                    
    plot_label (subplot, 'label', xlim, ylim, 
                x_percentage=x1, y_percentage=this_model_y1, color='black', xlog=0, ylog=0, 
                label=prefix_this_model, fontsize=10, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim,
                x_percentage=x21, y_percentage=this_model_y2, color='red', x2_percentage=x22, 
                xlog=0, ylog=0, linestyle='-', linewidth=2)
                
#endf plot_label_three_models
            
                
def plot_label (subplot, label_type, xlim, ylim, x_percentage, y_percentage, color, 
                x2_percentage=0., xlog=0, ylog=0, label='', linestyle='-', linewidth=2, 
                fontsize=16, fontweight='normal', sym='o', sym_size=5, err_size=0.1, 
                rotation=0,backgroundcolor='none', alpha=1.,back_alpha=0.,mfc='black'):    
    
    if(mfc=='black'):
        mfc=color
    
    if xlog==0 & ylog==0:
      
        if label_type=='label':
            x=xlim[0]+(xlim[1]-xlim[0])*x_percentage
            y=ylim[0]+(ylim[1]-ylim[0])*y_percentage             
            t=subplot.text(x,y,label, fontsize=fontsize, fontweight=fontweight,rotation=rotation, 
                           color=color, backgroundcolor=backgroundcolor, alpha=alpha)
            t.set_bbox(dict(alpha=back_alpha, color=backgroundcolor))
        else:
            if label_type =='line':
                x1=xlim[0]+(xlim[1]-xlim[0])*x_percentage
                x2=xlim[0]+(xlim[1]-xlim[0])*x2_percentage
                y=ylim[0]+(ylim[1]-ylim[0])*y_percentage
                subplot.plot([x1,x2],[y,y],color=color,linestyle=linestyle, linewidth=linewidth)
                
            else:
                if label_type=='symbol':
                    x=xlim[0]+(xlim[1]-xlim[0])*x_percentage
                    y=ylim[0]+(ylim[1]-ylim[0])*y_percentage                     
                    subplot.errorbar(x,y, yerr=err_size, fmt=sym, markersize=sym_size, 
                                     color=color, markeredgecolor=color,mfc=mfc)
                  
 
#end plot_label


def plot_joint_MR_MRII(hist_MR, hist_MRII, cut_MR_MRII, Volume_MR, Volume_MRII, 
                       bin, subplot, color='red',linewidth=2, linestyle='-'):
        
    x_axis_MRII=hist_MRII[1][0:len(hist_MRII[1][:])-1]+bin/2.
    sel=x_axis_MRII<cut_MR_MRII
    hist_MRII=hist_MRII[0][sel]
    x_axis_MRII=x_axis_MRII[sel]
        
    x_axis_MR=hist_MR[1][0:len(hist_MR[1][:])-1]+bin/2.
    sel=x_axis_MR>cut_MR_MRII
    hist_MR=hist_MR[0][sel]
    x_axis_MR=x_axis_MR[sel]
        
    x_axis=np.concatenate((x_axis_MRII, x_axis_MR), axis=0)
    y_axis=np.concatenate((np.log10(hist_MRII/(Volume_MRII*bin)),np.log10(hist_MR/(Volume_MR*bin))), axis=0)
    subplot.plot(x_axis,y_axis, color=color, linewidth=linewidth, linestyle=linestyle) 
        
    return (x_axis, y_axis)    

#end join_MR_MRII
        
    
def mag_to_lum(mag):

    return 10.**(-0.4*(mag-4.68))

#end 

def get_slope(x1=1., y1=1., x2=1., y2=1.):

    slope=(y2-y1)/(x2-x1)
    b=y1-(y2-y1)/(x2-x1)*x1
    
    return (slope,b)
#end


def read_file(fa):
    fa = open(fa, "r") 
    index=0
    for line in fa:
        if(index==0):                
            fields = line.strip().split()  
            N_elements=int(fields[0])
            x_axis=np.zeros(N_elements,dtype=np.float32)
            y_axis=np.zeros(N_elements,dtype=np.float32)               
        else:
            fields = line.strip().split()               
            x_axis[index-1]=float(fields[0])
            y_axis[index-1]=float(fields[1])               
        index+=1    
        
    return(x_axis,y_axis)  

#end read_file

def read_data_with_one_err(fa):
    fa = open(fa, "r") 
    index=0
    for line in fa:
        if(index==0):                
            fields = line.strip().split()             
            N_elements=int(fields[0])          
            x_axis=np.zeros(N_elements,dtype=np.float32)
            y_axis=np.zeros(N_elements,dtype=np.float32)
            y_err=np.zeros(N_elements,dtype=np.float32)           
        else:
            fields = line.strip().split()               
            x_axis[index-1]=float(fields[0])
            y_axis[index-1]=float(fields[1])  
            y_err[index-1]=float(fields[2])          
        index+=1    
        
    return(x_axis,y_axis, y_err)  

def read_data_with_err(fa):
    fa = open(fa, "r") 
    index=0
    for line in fa:
        if(index==0):                
            fields = line.strip().split()             
            N_elements=int(fields[0])
            x_axis=np.zeros(N_elements,dtype=np.float32)
            y_axis=np.zeros(N_elements,dtype=np.float32)
            y_err_down=np.zeros(N_elements,dtype=np.float32)
            y_err_up=np.zeros(N_elements,dtype=np.float32)
        else:
            fields = line.strip().split()               
            x_axis[index-1]=float(fields[0])
            y_axis[index-1]=float(fields[1])  
            y_err_down[index-1]=float(fields[2])
            y_err_up[index-1]=float(fields[3])
        index+=1    
        
    return(x_axis,y_axis, y_err_down, y_err_up)  

#end read_file

# <codecell>


# <codecell>

def smooth(array, N):

    for jj in range(0, N):
        for ii in range (1, len(array)-2):
            array[ii]=(array[ii-1]+array[ii]+array[ii+1])/3.
        
    return(array)
#end

def convert_2d_array_into_region(xx, yy):

    Nbins=len(xx)    
    median=np.zeros(Nbins, np.float32)
    mean=np.zeros(Nbins, np.float32)
    pc16=np.zeros(Nbins, np.float32) 
    pc84=np.zeros(Nbins, np.float32)  
    rms=np.zeros(Nbins, np.float32) 
    
    for ii in range(0,Nbins):       
        aux_y=yy[ii,:]
        aux_y=aux_y[aux_y>0.]        
        if(len(aux_y)>0):
            median[ii]=np.median(aux_y) 
            mean[ii]=np.mean(aux_y) 
            y_sorted = np.sort(aux_y)
            pc16[ii] = y_sorted[int(16*len(aux_y)/100)]      
            pc84[ii] = y_sorted[int(84*len(aux_y)/100)]              
            rms[ii]=np.sqrt(1./len(aux_y)*np.sum((mean[ii]-aux_y)**2))
       
    sel=median>0.
    return(xx[sel],median[sel],mean[sel],pc16[sel],pc84[sel],rms[sel])
#end

    
      
           
    #endfor
    
def abs_to_app_mag(abs_M, redshift, Hubble_h, Omega_m):
    
    cdist=comdist(redshift, Hubble_h, Omega_m)
    lumdist=(1.+redshift)*cdist
    return(abs_M+5.*(np.log10(lumdist)+6-1))

#end


def comdist(redshift, Hubble_h, Omega_m):
 
    H0 = Hubble_h*100. 
    WM = Omega_m                       #Omega(matter)                
    WV = 1.0 - WM - 0.4165/(H0*H0)     #Omega(lambda)
    WR = 4.165E-5/(Hubble_h*Hubble_h)  #Omega(radiation)
    WK = 1-WM-WR-WV                    #Omega curvaturve = 1-Omega(total)
 
    c = 299792.458   #velocity of light in km/sec 
    DCMR = 0.0       #comoving radial distance in units of c/H0  
    a = 1.0          # 1/(1+z), the scale factor of the Universe
    az = 0.5         # 1/(1+z(object)) 
    az = 1.0/(1+1.0*redshift)
 
    n=1000   # number of points in integrals

  

    #do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
    for i in range(0, n):    
        a = az+(1-az)*(i+0.5)/n
        adot = np.sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))    
        DCMR = DCMR + 1./(a*adot)  
 
    DCMR = (1.-az)*DCMR/n 
    DCMR_Mpc = (c/H0)*DCMR

    return(DCMR_Mpc)

#end



def plot_mass_function(subplot, Mass_MR, Volume_MR, xlim, bin, MRII, Mass_MRII=0, Volume_MRII=0, linestyle='-', color='black'):
    
    bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
    hist_MR=np.histogram(Mass_MR, bins=bin_arr, range=(xlim[0],xlim[1]))   
       
    #MRII
    if(MRII==1):    
        bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
        hist_MRII=np.histogram(Mass_MRII, bins=bin_arr, range=(xlim[0],xlim[1]))   
           
        
    #join MR+MRII & plot     
    if(MRII==1):
        cut_MR_MRII=10.
        (x_axis,y_axis)=plot_joint_MR_MRII(hist_MR, hist_MRII, cut_MR_MRII, Volume_MR, Volume_MRII, 
                                           bin, subplot, color=color,linewidth=2, linestyle=linestyle)
    else:
        x_axis=hist_MR[1][0:len(hist_MR[1][:])-1]+bin/2.           
        hist_MR=hist_MR[0]       
        y_axis=np.log10(hist_MR/(Volume_MR*bin))
        subplot.plot(x_axis,y_axis, color=color, linewidth=2, linestyle=linestyle) 
        
    return(x_axis,y_axis)
            
#end plot_mass_function