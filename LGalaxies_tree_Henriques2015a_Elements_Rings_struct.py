import numpy as np

LGalaxiesStruct = np.dtype([
('GalID',np.int64,1),
('HaloID',np.int64,1),
('FirstProgGal',np.int64,1),
('NextProgGal',np.int64,1),
('LastProgGal',np.int64,1),
('FOFCentralGal',np.int64,1),
('FileTreeNr',np.int64,1),
('DescendantGal',np.int64,1),
('MainLeafId',np.int64,1),
('TreeRootId',np.int64,1),
('SubID',np.int64,1),
('MMSubID',np.int64,1),
('PeanoKey',np.int32,1),
('Redshift',np.float32,1),
('Type',np.int32,1),
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
#('GasSpin',np.float32,3),
#('StellarSpin',np.float32,3),
('DiskSpin',np.float32,3),        
('InfallVmax',np.float32,1),
('InfallVmaxPeak',np.float32,1),
('InfallSnap',np.int32,1),
('InfallHotGas',np.float32,1),
('HotRadius',np.float32,1),
('OriMergTime',np.float32,1),
('MergTime',np.float32,1),
('ColdGas',np.float32,1),
('ColdGasRings',np.float32,12),
('H2fraction',np.float32,1), 
('H2fractionRings',np.float32,12),              
('StellarMass',np.float32,1),
('BulgeMass',np.float32,1),
('DiskMass',np.float32,1),
('DiskMassRings',np.float32,12),        
('HotGas',np.float32,1),
('EjectedMass',np.float32,1),
('BlackHoleMass',np.float32,1),
('ICM',np.float32,1),
('MetalsColdGas',np.float32,3),
('MetalsColdGasRings',np.float32,[3,12]),        
('MetalsBulgeMass',np.float32,3),
('MetalsDiskMass',np.float32,3),
('MetalsDiskMassRings',np.float32,[3,12]),        
('MetalsHotGas',np.float32,3),
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
('sfh_ICM',np.float32,20),
('sfh_MetalsDiskMass',np.float32,[3,20]),
('sfh_MetalsBulgeMass',np.float32,[3,20]),
('sfh_MetalsICM',np.float32,[3,20]),
('sfh_ElementsDiskMass',np.float32,[11,20]),
('sfh_ElementsBulgeMass',np.float32,[11,20]),
('sfh_ElementsICM',np.float32,[11,20]),
('DiskMass_elements',np.float32,11),
('DiskMass_elementsRings',np.float32,[11,12]),        
('BulgeMass_elements',np.float32,11),
('ColdGas_elements',np.float32,11),
('ColdGas_elementsRings',np.float32,[11,12]),        
('HotGas_elements',np.float32,11),
('ICM_elements',np.float32,11),
('EjectedMass_elements',np.float32,11)        
])


PropertiesToRead_tree = {}
for ii in LGalaxiesStruct.names:
	PropertiesToRead_tree[ii] = False
           
PropertiesToRead_tree['GalID'] = True 
PropertiesToRead_tree['HaloID'] = True
PropertiesToRead_tree['FirstProgGal'] = True
PropertiesToRead_tree['NextProgGal'] = True
PropertiesToRead_tree['LastProgGal'] = True
PropertiesToRead_tree['FOFCentralGal'] = True
#PropertiesToRead_tree['FileTreeNr'] = True
PropertiesToRead_tree['DescendantGal'] = True
PropertiesToRead_tree['MainLeafId'] = True
#PropertiesToRead_tree['TreeRootId'] = True
#PropertiesToRead_tree['SubID'] = True
#PropertiesToRead_tree['MMSubID'] = True           
PropertiesToRead_tree['Redshift'] = True
           
PropertiesToRead_tree['Type'] = True
#PropertiesToRead_tree['HaloIndex'] = True
PropertiesToRead_tree['SnapNum'] = True
#PropertiesToRead_tree['LookBackTimeToSnap'] = True
#PropertiesToRead_tree['CentralMvir'] = True
#PropertiesToRead_tree['CentralRvir'] = True
PropertiesToRead_tree['DistanceToCentralGal'] = True
PropertiesToRead_tree['Pos'] = True
PropertiesToRead_tree['Vel'] = True
#PropertiesToRead_tree['Len'] = True
PropertiesToRead_tree['Mvir'] = True
PropertiesToRead_tree['Rvir'] = True
PropertiesToRead_tree['Vvir'] = True
#PropertiesToRead_tree['Vmax'] = True
#PropertiesToRead_tree['GasSpin'] = True
#PropertiesToRead_tree['StellarSpin'] = True
#PropertiesToRead_tree['InfallVmax'] = True
#PropertiesToRead_tree['InfallVmaxPeak'] = True
#PropertiesToRead_tree['InfallSnap'] = True
#PropertiesToRead_tree['InfallHotGas'] = True
#PropertiesToRead_tree['HotRadius'] = True
#PropertiesToRead_tree['OriMergTime'] = True
#PropertiesToRead_tree['MergTime'] = True
PropertiesToRead_tree['ColdGas'] = True
PropertiesToRead_tree['ColdGasRings'] = True
PropertiesToRead_tree['H2fraction'] = True
PropertiesToRead_tree['H2fractionRings'] = True
PropertiesToRead_tree['StellarMass'] = True
PropertiesToRead_tree['BulgeMass'] = True
PropertiesToRead_tree['DiskMass'] = True
PropertiesToRead_tree['DiskMassRings'] = True
PropertiesToRead_tree['HotGas'] = True
#PropertiesToRead_tree['EjectedMass'] = True
PropertiesToRead_tree['BlackHoleMass'] = True
#PropertiesToRead_tree['ICM'] = True
PropertiesToRead_tree['MetalsColdGas'] = True
PropertiesToRead_tree['MetalsColdGasRings'] = True
#PropertiesToRead_tree['MetalsStellarMass'] = True
PropertiesToRead_tree['MetalsBulgeMass'] = True
PropertiesToRead_tree['MetalsDiskMass'] = True
#PropertiesToRead_tree['MetalsHotGas'] = True
#PropertiesToRead_tree['MetalsEjectedMass'] = True
#PropertiesToRead_tree['MetalsICM'] = True
#PropertiesToRead_tree['PrimordialAccretionRate'] = True
#PropertiesToRead_tree['CoolingRadius'] = True
#PropertiesToRead_tree['CoolingRate'] = True
#PropertiesToRead_tree['CoolingRate_beforeAGN'] = True
#PropertiesToRead_tree['QuasarAccretionRate'] = True
#PropertiesToRead_tree['RadioAccretionRate'] = True
PropertiesToRead_tree['Sfr'] = True
PropertiesToRead_tree['SfrRings'] = True
#PropertiesToRead_tree['SfrBulge'] = True
#PropertiesToRead_tree['XrayLum'] = True
PropertiesToRead_tree['BulgeSize'] = True
PropertiesToRead_tree['StellarDiskRadius'] = True
PropertiesToRead_tree['GasDiskRadius'] = True
#PropertiesToRead_tree['CosInclination'] = True
#PropertiesToRead_tree['DisruptOn'] = True
#PropertiesToRead_tree['MergeOn'] = True
PropertiesToRead_tree['MagDust'] = True
#PropertiesToRead_tree['Mag'] = True
#PropertiesToRead_tree['MagBulge'] = True
PropertiesToRead_tree['MassWeightAge'] = True
PropertiesToRead_tree['rBandWeightAge'] = True
#PropertiesToRead_tree['sfh_ibin'] = True
#PropertiesToRead_tree['sfh_numbins'] = True
#PropertiesToRead_tree['sfh_DiskMass'] = True
#PropertiesToRead_tree['sfh_BulgeMass'] = True
#PropertiesToRead_tree['sfh_ICM'] = True
#PropertiesToRead_tree['sfh_MetalsDiskMass'] = True
#PropertiesToRead_tree['sfh_MetalsBulgeMass'] = True
#PropertiesToRead_tree['sfh_MetalsICM'] = True
#PropertiesToRead_tree['DiskMass_elements'] = True
#PropertiesToRead_tree['DiskMass_elementsRings'] = True
#PropertiesToRead_tree['BulgeMass_elements'] = True
#PropertiesToRead_tree['ColdGas_elements'] = True
#PropertiesToRead_tree['ColdGas_elementsRings'] = True
#PropertiesToRead_tree['HotGas_elements'] = True
#PropertiesToRead_tree['ICM_elements'] = True
#PropertiesToRead_tree['EjectedMass_elements'] = True

# <codecell>


