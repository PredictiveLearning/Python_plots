'''
basic_animation()
anime_mass_gr()
'''
AnimDir='./animations/'

import numpy as np
from pylab import *

import matplotlib.pyplot as plt
from importlib import reload
from scipy.ndimage import zoom 
from matplotlib import animation
from matplotlib.colors import LinearSegmentedColormap

import procedures
reload (procedures)
from procedures import *
import plots_input
reload (plots_input)
from plots_input import *


def basic_animation():
    
    fig = plt.figure()
    ax = plt.axes(xlim=(0, 2), ylim=(-2, 2))
    line, = ax.plot([], [], lw=2)    
   
    # initialization function: plot the background of each frame
    def init():
        line.set_data([], [])
        return line,

    # animation function.  This is called sequentially
    def animate(i):
        x = np.linspace(0, 2, 1000)
        y = np.sin(2 * np.pi * (x - 0.01 * i))
        line.set_data(x, y)
        return line,

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=200, interval=20, blit=True)
    
    anim.save(AnimDir+'basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
    plt.show()
    
#end basic_animation

def anime_mass_gr(G_MR):
    xlim=[9.7,11.5]
    ylim=[0.1, 1.0]
    bin=[0.1,0.05]
    Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]
    
    fig = plt.figure(figsize=(12,8))    
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)        
    #ax = plt.axes(xlim=xlim, ylim=ylim)      
    xlab='$log_{10}(M_*[M_{\odot}])$'           
    ylab='$g-r$'               
    subplot.set_xlabel(xlab, fontsize=22), subplot.set_ylabel(ylab, fontsize=22)  
    yy= np.linspace(0.,1.,1000.)    
    # animation function
    
    def animate(nframe): 
        print('doing frame:',nframe)  
        #clear previous image
        clf()
        subplot=plt.subplot()
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)         
        subplot.set_xlabel(xlab, fontsize=22), subplot.set_ylabel(ylab, fontsize=22)               
        x=[xlim[0],xlim[1]]    
        #fill the entire plot area to delete previous frame, plot axis
        subplot.fill_between(x, [ylim[1],ylim[1]], [ylim[0],ylim[0]], facecolor='lavender')                
        for jj in range (1,10):
            value=9.5++jj/4.
            shift=0.003
            subplot.fill_between([value-shift,value+shift], [ylim[1],ylim[1]], [ylim[0],ylim[0]], facecolor='white',edgecolor="w")
            value=jj/10.
            shift=0.002            
            subplot.fill_between(x, [value+shift,value+shift], [value-shift,value-shift], facecolor='white',edgecolor="w")
                  
        #label        
        G0=G_MR[G_MR['SnapNum']==nframe+20] 
        props = dict(boxstyle='round', facecolor='ghostwhite')
        subplot.text(10.83, 0.165, 'z='+'%0.2f\n' % G0['Redshift'][0] 
                     + 'Age='+'%0.2f' % redshift_to_time(G0['Redshift'][0])+' Gyr', 
                     fontsize= 24, bbox=props)
            
        #ignore first 20 snaps which have no galaxies 
        G0_MR=G_MR[(G_MR['StellarMass']>0.1) & (G_MR['SnapNum']==nframe+20)]  
        Ngals=len(G0_MR)
        StellarMass=(np.log10(G0_MR['StellarMass']*1.e10/Hubble_h))     
        color_gr=G0_MR['MagDust'][:,16]-G0_MR['MagDust'][:,17]                   
        H, xedges, yedges = np.histogram2d(StellarMass, color_gr, bins=Nbins)            
        extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
        subplots_adjust(bottom=0.15, left=0.15)        
        mylevels = np.linspace(1., Nbins[0], Nbins[0])*Ngals/(Nbins[0]**2/0.7)        
        H = zoom(H, 20)    
        cont=contourf(H.transpose()[::], origin='lower', cmap='rainbow', levels=mylevels, extent=extent)        
        plt.colorbar()        
        return cont  
    #end animate
    
    anim = animation.FuncAnimation(fig, animate, frames=39, blit=True)
    #plt.colorbar()
    anim.save(AnimDir+'mass_gr.mp4', fps=1)
    
print('finished')
