import numpy as np

LGalaxiesStruct = np.dtype([
('Type',np.int32,1),
('HaloIndex',np.int32,1),
('SnapNum',np.int32,1),
('CentralMvir',np.float32,1),
('Len',np.int32,1),
('Mvir',np.float32,1),
('Rvir',np.float32,1),
('Vvir',np.float32,1),
('Vmax',np.float32,1),
('InfallVmax',np.float32,1),
('ColdGas',np.float32,1),
('ColdGasRings',np.float32,12), 
('StellarMass',np.float32,1),        
('BulgeMass',np.float32,1),
('DiskMass',np.float32,1),
('DiskMassRings',np.float32,12),        
('HotGas',np.float32,1),
('EjectedMass',np.float32,1),
('BlackHoleMass',np.float32,1),
('BlackHoleGas',np.float32,1),
('H2fraction',np.float32,1), 
('H2fractionRings',np.float32,12),    
('MetalsColdGas',np.float32,1),
('MetalsColdGasRings',np.float32,12),  
('MetalsStellarMass',np.float32,1),
('MetalsBulgeMass',np.float32,1),
('MetalsDiskMass',np.float32,1),
('MetalsDiskMassRings',np.float32,12),        
('MetalsHotGas',np.float32,1),
('MetalsEjectedMass',np.float32,1),
('Sfr',np.float32,1),
('SfrRings',np.float32,12),        
('SfrBulge',np.float32,1),
('XrayLum',np.float32,1),
('BulgeSize',np.float32,1),
('StellarDiskRadius',np.float32,1),
('CosInclination',np.float32,1),
('DisruptOn',np.int32,1),
('MergeOn',np.int32,1)
])

PropertiesToRead = {}
for ii in LGalaxiesStruct.names:
	PropertiesToRead[ii] = False
           
PropertiesToRead['Type'] = True
#PropertiesToRead['HaloIndex'] = True
PropertiesToRead['SnapNum'] = True
#PropertiesToRead['LookBackTimeToSnap'] = True
#PropertiesToRead['CentralMvir'] = True
#PropertiesToRead['CentralRvir'] = True

#PropertiesToRead['Len'] = True
PropertiesToRead['Mvir'] = True
PropertiesToRead['Rvir'] = True
PropertiesToRead['Vvir'] = True
#PropertiesToRead['Vmax'] = True
#PropertiesToRead['GasSpin'] = True
#PropertiesToRead['StellarSpin'] = True
#PropertiesToRead['InfallVmax'] = True
#PropertiesToRead['InfallVmaxPeak'] = True
#PropertiesToRead['InfallSnap'] = True
#PropertiesToRead['InfallHotGas'] = True
#PropertiesToRead['HotRadius'] = True
#PropertiesToRead['OriMergTime'] = True
#PropertiesToRead['MergTime'] = True
PropertiesToRead['ColdGas'] = True
PropertiesToRead['ColdGasRings'] = True
PropertiesToRead['StellarMass'] = True
PropertiesToRead['H2fraction'] = True
PropertiesToRead['H2fractionRings'] = True
PropertiesToRead['BulgeMass'] = True
PropertiesToRead['DiskMass'] = True
PropertiesToRead['DiskMassRings'] = True
PropertiesToRead['HotGas'] = True
#PropertiesToRead['EjectedMass'] = True
PropertiesToRead['BlackHoleMass'] = True
#PropertiesToRead['ICM'] = True
PropertiesToRead['MetalsColdGas'] = True
PropertiesToRead['MetalsColdGasRings'] = True
PropertiesToRead['MetalsStellarMass'] = True
PropertiesToRead['MetalsBulgeMass'] = True
PropertiesToRead['MetalsDiskMass'] = True
PropertiesToRead['MetalsDiskMassRings'] = True
#PropertiesToRead['MetalsHotGas'] = True
#PropertiesToRead['MetalsEjectedMass'] = True
#PropertiesToRead['MetalsICM'] = True
#PropertiesToRead['PrimordialAccretionRate'] = True
#PropertiesToRead['CoolingRadius'] = True
#PropertiesToRead['CoolingRate'] = True
#PropertiesToRead['CoolingRate_beforeAGN'] = True
#PropertiesToRead['QuasarAccretionRate'] = True
#PropertiesToRead['RadioAccretionRate'] = True
PropertiesToRead['Sfr'] = True
PropertiesToRead['SfrRings'] = True
#PropertiesToRead['SfrBulge'] = True
#PropertiesToRead['XrayLum'] = True
PropertiesToRead['BulgeSize'] = True
PropertiesToRead['StellarDiskRadius'] = True
#PropertiesToRead['GasDiskRadius'] = True
#PropertiesToRead['CosInclination'] = True
#PropertiesToRead['DisruptOn'] = True
#PropertiesToRead['MergeOn'] = True

        

# <codecell>


