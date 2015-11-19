# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

"""
read_snap
read_tree
redshift_to_time
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


def read_snap(folder,FirstFile,LastFile,
              props,template,RedshiftsToRead,RedshiftList):    
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
    for iredshift in range(0,len(RedshiftList)-1):
        if RedshiftsToRead[iredshift]:              
            for ifile in range(FirstFile,LastFile+1):
                char_redshift="%0.2f" % RedshiftList[iredshift]
                filename = folder+'/'+'SA_z'+char_redshift+"_"+"%d"%(ifile)               
                f = open(filename,"rb")
                
                this_nTrees =  np.fromfile(f,np.int32,1)
                nTrees += this_nTrees
                this_nGals = np.fromfile(f,np.int32,1)
                nGals += this_nGals
                
            print ("z=", char_redshift," nGals = ",nGals)  
               
    gals = np.zeros(nGals,dtype=filter_dtype)
  
    print ("\n")
    offset=0
    for iredshift in range(0,len(RedshiftList)-1):
        if RedshiftsToRead[iredshift]:  
            print ("\nReading redshift: ", RedshiftList[iredshift], "\n")
            for ifile in range(FirstFile,LastFile+1):
                char_redshift="%0.2f" % RedshiftList[iredshift]
                filename = folder+'/'+'SA_z'+char_redshift+"_"+"%d"%(ifile)
                #print(filename)
                f = open(filename,"rb")
                
                this_nTrees =  np.fromfile(f,np.int32,1)
                nTrees += this_nTrees
                this_nGals = np.fromfile(f,np.int32,1)
                nGals += this_nGals
                print ("File ", ifile," nGals = ",this_nGals)  
                
                addednTreeHalos = np.fromfile(f,np.int32,this_nTrees)
                nTreeHalos = np.append(nTreeHalos,addednTreeHalos)
                full_this_gals = np.fromfile(f,template,this_nGals) # all properties
                this_gals = np.zeros(this_nGals,dtype=filter_dtype) # selected props
                
                for prop in template.names:
                    if props[prop]:
                        this_gals[prop] = full_this_gals[prop]
                              
                gals[offset:offset+this_nGals] = this_gals[:]    
                offset+=this_nGals
                f.close()           
            #endfor
        #endif    
        #assign snapshot of current redshift given by the last galaxy on the last file 
        SnapshotList=np.append(SnapshotList,gals['SnapNum'][offset-1])       
    #endfor
    
    return (gals, SnapshotList)




def read_tree(folder,FirstFile,LastFile,
              props,template,RedshiftsToRead,RedshiftList):    
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
            
    SnapshotList=np.array([],dtype=np.int32)
    
    #read only headers to figure out total nGals
    print ("\n\nReading Headers\n")
    for ifile in range(FirstFile,LastFile+1):       
        filename = folder+'/'+'SA_galtree_'+"%d"%(ifile)               
        f = open(filename,"rb")
        one = np.fromfile(f,np.int32,1)
        nbytes = np.fromfile(f,np.int32,1)
        this_nGals = np.fromfile(f,np.int32,1)
        nGals += this_nGals               
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
        this_nGals = np.fromfile(f,np.int32,1)      
        nGals += this_nGals       
        print ("File ", ifile," nGals = ",this_nGals)      
        ib=np.fromfile(f,np.float32,int(nskip))          
       
        full_this_gals = np.fromfile(f,template,this_nGals) # all properties      
        this_gals = np.zeros(this_nGals,dtype=filter_dtype) # selected props
               
       
        for prop in template.names:
            if props[prop]:
                this_gals[prop] = full_this_gals[prop]
                              
        gals[offset:offset+this_nGals] = this_gals[:]    
        offset+=this_nGals
        f.close()           
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



# <codecell>


# <codecell>


