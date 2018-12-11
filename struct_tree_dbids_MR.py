import numpy as np

LGalaxiesStruct = np.dtype([
('HaloID',np.int64,1),
('FileTreeNrID',np.int64,1), 
('FirstProgenitorID',np.int64,1),
('LastProgenitorID',np.int64,1), 
('NextProgenitorID',np.int64,1),
('DescendantID',np.int64,1), 
('FirstHaloInFOFgroupID',np.int64,1),
('NextHaloInFOFgroupID',np.int64,1),  
#('MainLeafID',np.int64,1),      
('Redshift',np.float32,1),  
('PeanoKey',np.int32,1),
('dummy',np.int64,1)
])


PropertiesToRead = {}
for ii in LGalaxiesStruct.names:
	PropertiesToRead[ii] = False
          
PropertiesToRead['HaloID'] = True 
PropertiesToRead['FileTreeNrID'] = True 
PropertiesToRead['FirstProgenitorID'] = True 
PropertiesToRead['LastProgenitorID'] = True 
PropertiesToRead['NextProgenitorID'] = True 
PropertiesToRead['DescendantID'] = True 
PropertiesToRead['FirstHaloInFOFgroupID'] = True 
PropertiesToRead['NextHaloInFOFgroupID'] = True 
PropertiesToRead['Redshift'] = True 
#PropertiesToRead['PeanoKey'] = True        
        