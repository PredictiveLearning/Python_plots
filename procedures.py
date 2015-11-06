# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import numpy as np

def read_snap(folder,file_prefix,firstfile,lastfile,props,template):    
    """ Reads L-Galaxy output files.
    Returns: (nTrees,nHalos,nTreeHalos,gals)
    Inputs: (folder,file_prefix,firstfile,lastfile,props,template)
    props - list of properties to return
    template - structure dtype definition from database """     
    nTrees = 0
    nHalos = 0    
    nTreeHalos = np.array([],dtype=np.int32)
    filter_list = []
    for prop in props:
        if props[prop]:
            filter_list.append((prop,template[prop]))
    filter_dtype = np.dtype(filter_list)
    gals = np.array([],dtype=filter_dtype)
    for ifile in range(firstfile,lastfile+1):
        filename = folder+'/'+file_prefix+"_"+"%d"%(ifile)
        f = open(filename,"rb")
        this_nTrees =  np.fromfile(f,np.int32,1)
        nTrees += this_nTrees
        this_nHalos = np.fromfile(f,np.int32,1)
        nHalos += this_nHalos
        #print (“File ", ifile," nGals = ",this_nHalos)
        print (1)
        addednTreeHalos = np.fromfile(f,np.int32,this_nTrees)
        nTreeHalos = np.append(nTreeHalos,addednTreeHalos)
        this_addedGalaxy = np.fromfile(f,template,this_nHalos) # all properties
        addedGalaxy = np.zeros(this_nHalos,dtype=filter_dtype) # selected props
        for prop in template.names:
            if props[prop]:
                addedGalaxy[prop] = this_addedGalaxy[prop]
        gals = np.append(gals,addedGalaxy)      
        f.close()
        
    return (nTrees,nHalos,nTreeHalos,gals)

# <codecell>


# <codecell>


