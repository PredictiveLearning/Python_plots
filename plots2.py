# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import numpy as np
%matplotlib inline
import matplotlib.pyplot as plt
from astropy.table import Table

reload(procedures)
reload(LGalaxies_Henriques2015a_struct)

from procedures import read_snap
from LGalaxies_Henriques2015a_struct import LGalaxiesStruct
from LGalaxies_Henriques2015a_struct import PropertiesToRead


Datadir = '/net/bootes/export/data1/data/'

folder = '/net/bootes/export/data1/SAM/test1/MR/'
file_prefix= 'SA_z0.00'
firstfile = 40
lastfile = 40

(nTrees,nHalos,nTreeHalos,gal) = read_snap(folder,file_prefix,firstfile,lastfile,PropertiesToRead,LGalaxiesStruct)
print np.log10(gal['StellarMass'][1:5]*1.e10)
#help(gal)

# <codecell>

plt.rcParams.update({'font.size': 18, 'xtick.labelsize': 18, 'ytick.labelsize': 18, 'axes.linewidth': 2, 
                     'xtick.major.size': 4, 'xtick.major.width': 2., 'xtick.minor.width': 1.,
                     'ytick.major.size': 4, 'ytick.major.width': 2., 'ytick.minor.width': 1.})

# <codecell>

obs = Table.read(Datadir + 'liwhite2009_SMF.txt', format='ascii')

fig = plt.figure(figsize=(9,7))
plt.ylim([-7.0,1.])
plt.xlim([7., 12.])
plt.ylabel('$log_{10}(\phi [h^3 Mpc^{-3} log_{10}(M^{-1})])$'), 
plt.xlabel('$log_{10}(M_*[h^{-2}M_{\odot}])$')

plt.errorbar(obs['col3'], np.log10(obs['col4']),yerr=obs['col5'], fmt='o', markersize=5)

sub = plt.subplot(111)

plt.text(7.3, 0.18, 'Observations used in MCMC')
plt.errorbar(7.2, 0.3, yerr=0.15, fmt='o', markersize=5, color='blue')



plt.tight_layout()
plt.savefig('./fig/plots.pdf')

# <codecell>


