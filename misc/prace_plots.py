
# coding: utf-8

# In[41]:

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.ticker as mtick

plt.rcParams.update({'font.size': 14, 'xtick.labelsize': 14, 'ytick.labelsize': 14, 'axes.linewidth': 2, 
                     'xtick.major.size': 6, 'xtick.major.width': 1.5, 
                     'ytick.major.size': 6, 'ytick.major.width': 1.5, 
                     'xtick.minor.size': 3, 'xtick.minor.width': 1.,                   
                     'ytick.minor.size': 3, 'ytick.minor.width': 1.})  


one_one_size_small=[5,4]

xlim=[9e4,6e5]
ylim=[9e4,6e5]
          
fig = plt.figure(figsize=(one_one_size_small[0],one_one_size_small[1]))
subplot=plt.subplot()
subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
            
#format axis
#subplot.xaxis.set_major_locator(style='sci')    
#subplot.xaxis.set_minor_locator(MultipleLocator(0.25)) 

#subplot.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
#subplot.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))

xlab='Number of MCMC steps'       
ylab='CPU time [hours]' 

subplot.set_xlabel(xlab, fontsize=14), subplot.set_ylabel(ylab, fontsize=14)  

x_array=np.array([100000.,200000.,300000.,400000.,500000.])
y_array=np.array([100000.,200000.,300000.,400000.,500000.])
subplot.scatter(x_array, y_array)
#subplot.set_yscale('log')
#subplot.set_xscale('log')

subplot.ticklabel_format(style='sci', scilimits=(0,0), axis='both')

plt.tight_layout()
plt.savefig('./prace_performance.pdf')   
plt.close()


# In[ ]:



