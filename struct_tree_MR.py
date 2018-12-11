import numpy as np

LGalaxiesStruct = np.dtype([
('Descendant',np.int32,1),
('FirstProgenitor',np.int32,1),
('NextProgenitor',np.int32,1),
('FirstHaloInFOFgroup',np.int32,1),
('NextHaloInFOFgroup',np.int32,1),
('Len',np.int32,1),
('M_Mean200',np.float32,1),
('M_Crit200',np.float32,1),
('M_TopHat',np.float32,1),
('Pos',np.float32,3),
('Vel',np.float32,3), 
('VelDisp',np.float32,1),
('Vmax',np.float32,1),
('Spin',np.float32,3),
('MostBoundID',np.int64,1),
('SnapNum',np.int32,1),
('FileNr',np.int32,1),
('SubhaloIndex',np.int32,1),
('SubHalfMass',np.float32,1)
])

PropertiesToRead = {}
for ii in LGalaxiesStruct.names:
	PropertiesToRead[ii] = False
          
#PropertiesToRead['Descendant'] = True 
#PropertiesToRead['FirstProgenitor'] = True
#PropertiesToRead['NextProgenitor'] = True
#PropertiesToRead['FirstHaloInFOFgroup'] = True
#PropertiesToRead['NextHaloInFOFgroup'] = True
PropertiesToRead['Len'] = True
#PropertiesToRead['M_Mean200'] = True
PropertiesToRead['M_Crit200'] = True
PropertiesToRead['M_TopHat'] = True
PropertiesToRead['Pos'] = True
#PropertiesToRead['Vel'] = True
#PropertiesToRead['VelDisp',] = True
#PropertiesToRead['Vmax'] = True
#PropertiesToRead['Spin'] = True
#PropertiesToRead['MostBoundID'] = True
PropertiesToRead['SnapNum'] = True
#PropertiesToRead['FileNr'] = True
#PropertiesToRead['SubhaloIndex'] = True
#PropertiesToRead['SubHalfMass'] = True



        
