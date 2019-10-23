# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import numpy as np

LGalaxiesStruct = np.dtype([
('Type',np.int32,1),
('HaloIndex',np.int32,1),
('SnapNum',np.int32,1),
('LookBackTimeToSnap',np.float32,1),
('CentralMvir',np.float32,1),
('CentralRvir',np.float32,1),
('DistanceToCentralGal',np.float32,3),
('Pos',np.float32,3),
('Vel',np.float32,3),
('Len',np.int32,1),
('Mvir',np.float32,1),
('Rvir',np.float32,1),
('Vvir',np.float32,1),
('Vmax',np.float32,1),
('HaloSpin',np.float32,3),    
('InfallVmax',np.float32,1),
('InfallVmaxPeak',np.float32,1),
('InfallSnap',np.int32,1),
('InfallHotGas',np.float32,1),
('HotRadius',np.float32,1),
('OriMergTime',np.float32,1),
('MergTime',np.float32,1),
#('flagSplashBack',np.int32,1), 
#('TimeSinceSlashBack',np.float32,1),   
#('NMajorMergers',np.float32,1),  
#('NMinorMergers',np.float32,1),   
('ColdGas',np.float32,1),
('ColdGasRings',np.float32,12),
('H2fraction',np.float32,1), 
('H2fractionRings',np.float32,12),              
('StellarMass',np.float32,1),
('DiskMass',np.float32,1),
('BulgeMass',np.float32,1),
('DiskMassRings',np.float32,12),   
('BulgeMassRings',np.float32,12),   
('HotGas',np.float32,1),
#('ReheatedGas',np.float32,1),
('EjectedMass',np.float32,1),
('BlackHoleMass',np.float32,1),
('ICM',np.float32,1),
('MassFromInSitu',np.float32,1),
('MassFromMergers',np.float32,1),
('MassFromBursts',np.float32,1),
('MetalsColdGas',np.float32,3),
('MetalsColdGasRings',np.float32,[12,3]),
('MetalsStellarMass',np.float32,3),
('MetalsDiskMass',np.float32,3),
('MetalsBulgeMass',np.float32,3),
('MetalsDiskMassRings',np.float32,[12,3]),  
('MetalsBulgeMassRings',np.float32,[12,3]),  
('MetalsHotGas',np.float32,3),
#('MetalsReheatedGas',np.float32,3),
('MetalsEjectedMass',np.float32,3),
('MetalsICM',np.float32,3), 
('PrimordialAccretionRate',np.float32,1),
('CoolingRadius',np.float32,1),
('CoolingRate',np.float32,1),
('CoolingRate_beforeAGN',np.float32,1),
('QuasarAccretionRate',np.float32,1),
('RadioAccretionRate',np.float32,1),
('Sfr',np.float32,1),
('SfrRings',np.float32,12),        
('SfrBulge',np.float32,1),
('XrayLum',np.float32,1),
('BulgeSize',np.float32,1),
('StellarDiskRadius',np.float32,1),
('GasDiskRadius',np.float32,1),
('StellarHalfMassRadius',np.float32,1),  
('StellarHalfLightRadius',np.float32,1),           
('CosInclination',np.float32,1),
('DisruptOn',np.int32,1),
('MergeOn',np.int32,1), 
('MagDust',np.float32,40),
('Mag',np.float32,40),
('MagBulge',np.float32,40),
('MassWeightAge',np.float32,1),
('rBandWeightAge',np.float32,1),
('sfh_ibin',np.int32,1),      
('sfh_numbins',np.int32,1),        
('sfh_DiskMass',np.float32,20),
('sfh_BulgeMass',np.float32,20),
('sfh_DiskMassRings',np.float32,[12,20]),
('sfh_BulgeMassRings',np.float32,[12,20]),
('sfh_ICM',np.float32,20),
('sfh_MetalsDiskMass',np.float32,[20,3]),
('sfh_MetalsBulgeMass',np.float32,[20,3]),
('sfh_MetalsICM',np.float32,[20,3]),      
##('sfh_MassFromInSitu',np.float32,20),
##('sfh_MassFromMergers',np.float32,20),
##('sfh_MassFromBurst',np.float32,20),
('sfh_ElementsDiskMass',np.float32,[20,11]),
('sfh_ElementsBulgeMass',np.float32,[20,11]),
('sfh_ElementsICM',np.float32,[20,11]),
('DiskMass_elements',np.float32,11),
('BulgeMass_elements',np.float32,11),
('DiskMassRings_elements',np.float32,[12,11]),        
('BulgeMassRings_elements',np.float32,[12,11]),    
('ColdGas_elements',np.float32,11),
('ColdGasRings_elements',np.float32,[12,11]),        
('HotGas_elements',np.float32,11),
##('ReheatedGas_elements',np.float32,11),
('ICM_elements',np.float32,11),
('EjectedMass_elements',np.float32,11)
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
#PropertiesToRead['DistanceToCentralGal'] = True
#PropertiesToRead['Pos'] = True
#PropertiesToRead['Vel'] = True
#PropertiesToRead['Len'] = True
#PropertiesToRead['Mvir'] = True
#PropertiesToRead['Rvir'] = True
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
#PropertiesToRead['flagSplashBack'] = True
#PropertiesToRead['TimeSinceSlashBack'] = True
PropertiesToRead['ColdGas'] = True
PropertiesToRead['ColdGasRings'] = True
PropertiesToRead['H2fraction'] = True
PropertiesToRead['H2fractionRings'] = True
PropertiesToRead['StellarMass'] = True
PropertiesToRead['DiskMass'] = True
PropertiesToRead['BulgeMass'] = True
PropertiesToRead['DiskMassRings'] = True
PropertiesToRead['BulgeMassRings'] = True
PropertiesToRead['HotGas'] = True
PropertiesToRead['EjectedMass'] = True
PropertiesToRead['BlackHoleMass'] = True
#PropertiesToRead['ICM'] = True
#PropertiesToRead['MassFromInSitu'] = True
#PropertiesToRead['MassFromMergers'] = True
#PropertiesToRead['MassFromBursts'] = True
PropertiesToRead['MetalsColdGas'] = True
PropertiesToRead['MetalsColdGasRings'] = True
PropertiesToRead['MetalsStellarMass'] = True
PropertiesToRead['MetalsDiskMass'] = True
PropertiesToRead['MetalsBulgeMass'] = True
PropertiesToRead['MetalsDiskMassRings'] = True
PropertiesToRead['MetalsBulgeMassRings'] = True
PropertiesToRead['MetalsHotGas'] = True
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
PropertiesToRead['GasDiskRadius'] = True
PropertiesToRead['StellarHalfMassRadius'] = True
PropertiesToRead['StellarHalfLightRadius'] = True
#PropertiesToRead['CosInclination'] = True
#PropertiesToRead['DisruptOn'] = True
#PropertiesToRead['MergeOn'] = True
PropertiesToRead['MagDust'] = True
#PropertiesToRead['Mag'] = True
#PropertiesToRead['MagBulge'] = True
PropertiesToRead['MassWeightAge'] = True
PropertiesToRead['rBandWeightAge'] = True
#PropertiesToRead['sfh_ibin'] = True
#PropertiesToRead['sfh_numbins'] = True

#PropertiesToRead['sfh_DiskMass'] = True
#PropertiesToRead['sfh_BulgeMass'] = True
#PropertiesToRead['sfh_DiskMassRings'] = True
#PropertiesToRead['sfh_BulgeMassRings'] = True
#PropertiesToRead['sfh_ICM'] = True
#PropertiesToRead['sfh_MetalsDiskMass'] = True
#PropertiesToRead['sfh_MetalsBulgeMass'] = True
#PropertiesToRead['sfh_MetalsICM'] = True

#PropertiesToRead['sfh_ElementsDiskMass'] = True
#PropertiesToRead['sfh_ElementsBulgeMass'] = True

#PropertiesToRead['DiskMass_elements'] = True
#PropertiesToRead['BulgeMass_elements'] = True
#PropertiesToRead['DiskMassRings_elements'] = True
#PropertiesToRead['BulgeMassRings_elements'] = True
#PropertiesToRead['ColdGas_elements'] = True
PropertiesToRead['ColdGasRings_elements'] = True
#PropertiesToRead['HotGas_elements'] = True

# <codecell>

