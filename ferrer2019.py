# !/usr/bin/env python
# -*- coding: utf-8 -*-

# import os.path
import numpy as np
from astropy.table import Table
# import astropy.units as u

import scipy.stats

# import multiprocessing as mp

# import pyezmad.cmap

# from pyezmad.utilities import (pixel_to_physical,
#                                pixel_to_arcsec,
#                                get_ned_velocity)

# import seaborn.apionly as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import ticker

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.font_manager
import matplotlib

matplotlib.rcParams['mathtext.fontset'] = 'stixsans'
matplotlib.rcParams['font.family'] = 'sans'
matplotlib.rcParams['font.serif'] = 'STIXGeneral'
matplotlib.rcParams['xtick.labelsize'] = 'xx-large'
matplotlib.rcParams['ytick.labelsize'] = 'xx-large'


def func_l68(x):
    if x.size > 50:
        return(np.nanpercentile(x, 16))
    else:
        return(np.nan)


def func_u68(x):
    if x.size > 50:
        return(np.nanpercentile(x, 84))
    else:
        return(np.nan)


def func_med(x):
    if x.size > 50:
        return(np.nanmedian(x))
    else:
        return(np.nan)


def func_std(x):
    if x.size > 50:
        return(np.nanstd(x))
    else:
        return(np.nan)


def get_stats(x, y, bins):

    med, bin_edges, binnum \
        = scipy.stats.binned_statistic(x, y,
                                       statistic=func_med,
                                       bins=bins,
                                       range=(bins[0], bins[-1]))
    l68, bin_edges, binnum \
        = scipy.stats.binned_statistic(x, y,
                                       statistic=func_l68,
                                       bins=bins,
                                       range=(bins[0], bins[-1]))
    u68, bin_edges, binnum \
        = scipy.stats.binned_statistic(x, y,
                                       statistic=func_u68,
                                       bins=bins,
                                       range=(bins[0], bins[-1]))

    return(med, l68, u68, bin_edges, binnum)

if __name__ == '__main__':


# ALL GALAXIES
    FerrerDir = '/scratch-ssd/data/ferrer2019/'
    objinfo_all = np.genfromtxt(FerrerDir + 'table_info_mass_all.txt', dtype=None, 
                                names=['name', '2a', 'b_a', 'pa','d', 're', 'xc', 'yc', 
                                       'mass', 'xplot', 'yplot','npix','deltaMS'])


    nobj = objinfo_all['name'].size    
    outfile='C_number_mass_both.pdf'

    fig = plt.figure(figsize=(6, 5))
    nh, nv = 1,1
    gs = gridspec.GridSpec(nv, nh)  
    gs.update(top=0.8, bottom=0.19, left=0.05, right=0.25,  wspace = 0)
   
   
    ax_each = []
    
    i_plot = 0

    med_each_obj = np.empty(nobj, dtype=np.object)
    l68_each_obj = np.empty(nobj, dtype=np.object)
    u68_each_obj = np.empty(nobj, dtype=np.object)

    x_all = np.array([])
    y_all = np.array([])
    z_all = np.array([])

    for i in range(nobj):

        objinfo = objinfo_all
        obj = objinfo['name'][i].decode('UTF-8').lower()       
        dir=FerrerDir
        tb = Table.read(dir+obj + '_allprop_halpha_sf.fits')
       
        z0=tb['lsmd']
        y0=tb['oh12_dopita']
        x0=tb['r_ell']*0.2/objinfo['re'][i]


        # mask = tb['npix'] > 100  # S/N=5 per pixel
        # mask = tb['npix'] > 277  # S/N=3 per pixel
        mask = np.isnan(x0)
        # mask = np.logical_or(mask, mask_param!=1)

        mask = np.logical_or(mask, np.isnan(y0))
        mask = np.logical_or(mask, np.isnan(z0))
        # mask = np.logical_or(mask, x0<0.5)
        # mask = np.logical_or(mask, x0>5.5)
        # mask = np.logical_or(mask, tb['nfev_hb'] >= 1500)
        # mask = np.logical_or(mask, ~tb['scs_ha'])
        # mask = np.logical_or(mask, ~tb['scs_hb'])
        # mask = np.logical_or(mask, tb['ppchisq1'] > 5)
        mask = np.logical_or(mask, tb['bpt_info'] !=1)
        # mask = np.logical_or(mask, tb['dig_info'] !=1)

        x = x0[~mask]
        y = y0[~mask]
        z = z0[~mask]
        x_all = np.concatenate((x_all, x))
        y_all = np.concatenate((y_all, y))
        z_all = np.concatenate((z_all, z))

        med_2d_gal, y_edges_gal, x_edges_gal, binnum_gal \
        = scipy.stats.binned_statistic_2d(
            y, x, z,
            statistic='count',
            # statistic='median',
            # statistic='sum',
            bins=(60, 60),
            # range=[[-6,2], [0, 8]])
                range=[[7.6,9.6], [0,5.5]])
        # med_2d_gal=med_2d_gal*0+1
        mask5=med_2d_gal==0
        med_2d_gal[~mask5]=1
        med_2d_gal[mask5]=0
        if i==0:
            med_2d_tot=med_2d_gal*0
        med_2d_tot=med_2d_tot+med_2d_gal



    d_bin = 0.25 / 2.
    bin_min = 0 + d_bin
    bin_max = 6

    nbin = int((bin_max - bin_min + 2 * d_bin) / (2 * d_bin)) + 1
    bins = np.linspace(bin_min - d_bin, bin_max + d_bin, nbin)


    med_each, l68_each, u68_each, bin_edges, binnum \
        = get_stats(x, y, bins)

    bin_width = (bin_edges[1] - bin_edges[0])
    bin_centers = bin_edges[1:] - bin_width / 2.


    #Number galaxy plot
    med_2d_super, y_edges_super, x_edges_super, binnum_super \
        = scipy.stats.binned_statistic_2d(
            y_all, x_all, z_all,
            statistic='count',
            # statistic='median',
            bins=(60, 60),
            # range=[[-6,2], [0, 8]])
            range=[[7.6,9.6], [0,5.5]])
    mask6=med_2d_super<5
    med_2d_tot[mask6]=0
    mask7=med_2d_tot==0
    med_2d_tot[mask7]='NaN'

    ax = fig.add_subplot(1,1,1)

    # cax1 = ax.imshow(med_2d_tot,
    #                   interpolation='nearest', origin='lower',
    #                   extent=[x_edges_gal[0], x_edges_gal[-1],
    #                           y_edges_gal[0], y_edges_gal[-1]],
    #                   vmin=1,
    #                   vmax=25,
    #                   cmap=cm.Greys,
    #                   aspect='auto')

    ax.set_xlim(0,2.95)
    ax.set_ylim(7.6,9.9)


    props = dict(facecolor='white', edgecolor='none', alpha=0.5)
    # ax.set_title("ALL GALAXIES",fontsize='xx-large')

    ax.set_ylabel('$12+log(O/H)$', fontsize=24)
    ax.set_xlabel('$R/R_{e}$', fontsize=24)

    # ax.set_yticklabels([])
    ax.set_yticks(ticks=[8.0,8.5,9.0,9.5], minor=False)
  
    # ax.set_xticklabels([])

    ax.set_xticks(ticks=[0,1,2,3,4,5], minor=False)
    ax.minorticks_on()
    # divider = make_axes_locatable(ax)
    # cbaxes = divider.append_axes("right", size = "5%", pad = 0.04)
    # cbar = plt.colorbar(cax1, cax = cbaxes, label="Number of galaxies")
    # # cbar.ax.set_xticklabels(['     SF', 'AGN/shock                    '])
    # # tick_locator = ticker.MaxNLocator(nbins=5)
    # # cbar.locator = tick_locator
    # cbar.ax.yaxis.label.set_font_properties(matplotlib.font_manager.FontProperties(size='xx-large'))
    # # cbar.update_ticks()
    # cbar.ax.tick_params(labelsize=16)

    # med_2d_tot=med_2d_tot-2
    med_2d_num, xedges_num, yedges_num = np.histogram2d(x_all, y_all, bins=(60, 60), range = [[0,5.5], [7.6,9.6]])
    x = .5 * (xedges_num[1:] + xedges_num[:-1])
    y = .5 * (yedges_num[1:] + yedges_num[:-1])

    # med_2d_num[mask6]=0
    # ax.contour(x, y, med_2d_num.T, levels = [5],  linestyles='dashed',colors = 'blue', linewidths=3)  #, zorder = 10)



####################GALAXIES WITH logM<10

    objinfo_all = np.genfromtxt(FerrerDir+'table_info_mass_low.txt', dtype=None,names=['name', '2a', 'b_a', 'pa','d', 're', 'xc', 'yc', 'mass','xplot','yplot','npix'])
    nobj = objinfo_all['name'].size

    # objinfo_all = np.genfromtxt('table_info_type'+type_num+bymass+'_mass.txt', dtype=None,
    #                             names=['name', '2a', 'b_a', 'pa',
    #                                    'd', 're', 'xc', 'yc','mass','sfr','rbreak'])
    # p.obj_ignore = 'PGC3853'
    med_each_obj = np.empty(nobj, dtype=np.object)
    l68_each_obj = np.empty(nobj, dtype=np.object)
    u68_each_obj = np.empty(nobj, dtype=np.object)

    x_all = np.array([])
    y_all = np.array([])
    z_all = np.array([])

    for i in range(nobj):

        objinfo = objinfo_all
        obj = objinfo['name'][i].decode('UTF-8').lower()  
        dir=FerrerDir

        tb = Table.read(dir+obj+ '_allprop_halpha_sf.fits')
      
        z0=tb['lsmd']
        y0=tb['oh12_dopita']
        x0=tb['r_ell']*0.2/objinfo['re'][i]

        mask = np.isnan(x0)        
        mask = np.logical_or(mask, np.isnan(y0))
        mask = np.logical_or(mask, np.isnan(z0))    
        mask = np.logical_or(mask, tb['bpt_info'] !=1)
        # mask = np.logical_or(mask, tb['dig_info'] !=1)

        x = x0[~mask]
        y = y0[~mask]
        z = z0[~mask]
        x_all = np.concatenate((x_all, x))
        y_all = np.concatenate((y_all, y))
        z_all = np.concatenate((z_all, z))

        med_2d_gal, y_edges_gal, x_edges_gal, binnum_gal \
        = scipy.stats.binned_statistic_2d(
            y, x, z,
            statistic='count',
            # statistic='median',
            # statistic='sum',
            bins=(60, 60),
            # range=[[-6,2], [0, 8]])
                range=[[7.6,9.6], [0,5.5]])
        # med_2d_gal=med_2d_gal*0+1
        mask5=med_2d_gal==0
        med_2d_gal[~mask5]=1
        med_2d_gal[mask5]=0
        if i==0:
            med_2d_tot=med_2d_gal*0
        med_2d_tot=med_2d_tot+med_2d_gal


    med_2d_super, y_edges_super, x_edges_super, binnum_super \
        = scipy.stats.binned_statistic_2d(
            y_all, x_all, z_all,
            statistic='count',
            # statistic='median',
            bins=(60, 60),
            # range=[[-6,2], [0, 8]])
            range=[[7.6,9.6], [0,5.5]])
    mask6=med_2d_super<10
    med_2d_tot[mask6]=0
    mask7=med_2d_tot==0
    med_2d_tot[mask7]='NaN'


    med, l68, u68, bin_edges, binnum \
        = get_stats(x_all, y_all, bin_edges)


    ax.errorbar(bin_centers, med, xerr=(bin_edges[1] - bin_edges[0]) / 2., yerr=[med - l68, u68 - med],
                fmt='s', ecolor='blue', elinewidth=2.5, capsize=0, ms=12, mew=1, mfc='blue',  mec='blue',
                zorder=3)
                # mec='blue',label=ur'log[M$_\mathregular{*}$/M$_\odot$] < 10')
    print('dop low mass')
    x_error = (bin_edges[1] - bin_edges[0]) / 2.
    print('radius,metallicity,err_x,err_y_down,err_y_up')
    for ii in range(0, len(bin_centers)):
        print('{:f},{:f},{:f},{:f},{:f}'.format(bin_centers[ii],med[ii],x_error, med[ii] - l68[ii], u68[ii] - med[ii]))
                   
    ax.set_xlim(0,2.95)
    ax.set_ylim(7.6,9.9)


####################GALAXIES WITH 10<logM<11

    objinfo_all = np.genfromtxt(FerrerDir+'table_info_mass_med.txt', dtype=None, 
                                names=['name', '2a', 'b_a', 'pa','d', 're', 'xc', 'yc', 'mass',
                                       'xplot','yplot', 'npix'])
    nobj = objinfo_all['name'].size

    nobj = objinfo_all['name'].size
    med_each_obj = np.empty(nobj, dtype=np.object)
    l68_each_obj = np.empty(nobj, dtype=np.object)
    u68_each_obj = np.empty(nobj, dtype=np.object)

    x_all = np.array([])
    y_all = np.array([])
    z_all = np.array([])

    for i in range(nobj):

        objinfo = objinfo_all
        obj =  obj = objinfo['name'][i].decode('UTF-8').lower()      
        dir=FerrerDir
        tb = Table.read(dir+obj+ '_allprop_halpha_sf.fits')
     
        z0=tb['lsmd']
        y0=tb['oh12_dopita']
        x0=tb['r_ell']*0.2/objinfo['re'][i]

        mask = np.isnan(x0)
        mask = np.logical_or(mask, np.isnan(y0))
        mask = np.logical_or(mask, np.isnan(z0))    
        mask = np.logical_or(mask, tb['bpt_info'] !=1)
        # mask = np.logical_or(mask, tb['dig_info'] !=1)

        x = x0[~mask]
        y = y0[~mask]
        z = z0[~mask]
        x_all = np.concatenate((x_all, x))
        y_all = np.concatenate((y_all, y))
        z_all = np.concatenate((z_all, z))

        med_2d_gal, y_edges_gal, x_edges_gal, binnum_gal \
        = scipy.stats.binned_statistic_2d(
            y, x, z,
            statistic='count',
            # statistic='median',
            # statistic='sum',
            bins=(60, 60),
            # range=[[-6,2], [0, 8]])
                range=[[7.6,9.6], [0,5.5]])
        # med_2d_gal=med_2d_gal*0+1
        mask5=med_2d_gal==0
        med_2d_gal[~mask5]=1
        med_2d_gal[mask5]=0
        if i==0:
            med_2d_tot=med_2d_gal*0
        med_2d_tot=med_2d_tot+med_2d_gal


    med, l68, u68, bin_edges, binnum \
        = get_stats(x_all, y_all, bin_edges)


    ax.errorbar(bin_centers, med, xerr=(bin_edges[1] - bin_edges[0]) / 2., yerr=[med - l68, u68 - med], 
                fmt='D', ecolor='blue', elinewidth=2.5,
                 capsize=0,
                 ms=12,
                 mew=1,
                 mfc='blue',                
                 mec='blue',zorder=3)
    print('dop inter mass')
    x_error = (bin_edges[1] - bin_edges[0]) / 2.
    print('radius,metallicity,err_x,err_y_down,err_y_up')
    for ii in range(0, len(bin_centers)):
        print('{:f},{:f},{:f},{:f},{:f}'.format(bin_centers[ii],med[ii],x_error, med[ii] - l68[ii], u68[ii] - med[ii]))
                     
    med_2d_num, xedges_num, yedges_num = np.histogram2d(x_all, y_all, bins=(60, 60), range = [[0,5.5], [7.6,9.6]])
    x = .5 * (xedges_num[1:] + xedges_num[:-1])
    y = .5 * (yedges_num[1:] + yedges_num[:-1])
    # ax.contour(x, y, med_2d_num.T, levels = [5],  linestyles='-',colors = 'blue', linewidths=3,label='Intermediate')  #, zorder = 10)


####################GALAXIES WITH logM>11

    objinfo_all = np.genfromtxt(FerrerDir+'table_info_mass_high.txt', dtype=None,names=['name', '2a', 'b_a', 'pa','d', 're', 'xc', 'yc', 'mass','xplot','yplot','npix'])
    nobj = objinfo_all['name'].size
    nobj = objinfo_all['name'].size

    med_each_obj = np.empty(nobj, dtype=np.object)
    l68_each_obj = np.empty(nobj, dtype=np.object)
    u68_each_obj = np.empty(nobj, dtype=np.object)

    x_all = np.array([])
    y_all = np.array([])
    z_all = np.array([])

    for i in range(nobj):

        objinfo = objinfo_all
        obj = objinfo['name'][i].decode('UTF-8').lower()     

        dir=FerrerDir
        tb = Table.read(dir+obj+ '_allprop_halpha_sf.fits')
   
        z0=tb['lsmd']
        y0=tb['oh12_dopita']
        x0=tb['r_ell']*0.2/objinfo['re'][i]

        mask = np.isnan(x0) 
        mask = np.logical_or(mask, np.isnan(y0))
        mask = np.logical_or(mask, np.isnan(z0))
     
        mask = np.logical_or(mask, tb['bpt_info'] !=1)
        # mask = np.logical_or(mask, tb['dig_info'] !=1)

        x = x0[~mask]
        y = y0[~mask]
        z = z0[~mask]
        x_all = np.concatenate((x_all, x))
        y_all = np.concatenate((y_all, y))
        z_all = np.concatenate((z_all, z))

        med_2d_gal, y_edges_gal, x_edges_gal, binnum_gal \
        = scipy.stats.binned_statistic_2d(
            y, x, z,
            statistic='count',
            # statistic='median',
            # statistic='sum',
            bins=(60, 60),
            # range=[[-6,2], [0, 8]])
                range=[[7.6,9.6], [0,5.5]])
        # med_2d_gal=med_2d_gal*0+1
        mask5=med_2d_gal==0
        med_2d_gal[~mask5]=1
        med_2d_gal[mask5]=0
        if i==0:
            med_2d_tot=med_2d_gal*0
        med_2d_tot=med_2d_tot+med_2d_gal

    med, l68, u68, bin_edges, binnum \
        = get_stats(x_all, y_all, bin_edges)


    ax.errorbar(bin_centers, med,
                 xerr=(bin_edges[1] - bin_edges[0]) / 2.,
                 yerr=[med - l68,
                       u68 - med],
                 fmt='o',
                 ecolor='blue',
                 elinewidth=2.5,
                 capsize=0,
                 ms=12,
                 mew=1,
                 mfc='blue',
                 mec='blue',zorder=3)
    print('dop high mass')
    x_error = (bin_edges[1] - bin_edges[0]) / 2.
    print('radius,metallicity,err_x,err_y_down,err_y_up')
    for ii in range(0, len(bin_centers)):
        print('{:f},{:f},{:f},{:f},{:f}'.format(bin_centers[ii],med[ii],x_error, med[ii] - l68[ii], u68[ii] - med[ii]))
              


    med_2d_num, xedges_num, yedges_num = np.histogram2d(x_all, y_all, bins=(60, 60), range = [[0,5.5], [7.6,9.6]])
    x = .5 * (xedges_num[1:] + xedges_num[:-1])
    y = .5 * (yedges_num[1:] + yedges_num[:-1])
    # ax.contour(x, y, med_2d_num.T, levels = [5],  linestyles='-',colors = 'blue', linewidths=3,label='High mass')  #, zorder = 10)

    ax.legend(loc=1, scatterpoints=1, frameon=False, fontsize=14)

    ax.set_xlim(0,2.95)
    ax.set_ylim(7.6,9.9)


# ALL GALAXIES

    objinfo_all = np.genfromtxt(FerrerDir + 'table_info_mass_all.txt', dtype=None,names=['name', '2a', 'b_a', 'pa','d', 're', 'xc', 'yc', 'mass','xplot','yplot','npix','deltaMS'])


    nobj = objinfo_all['name'].size
  
    
    # i_plot = 0

    med_each_obj = np.empty(nobj, dtype=np.object)
    l68_each_obj = np.empty(nobj, dtype=np.object)
    u68_each_obj = np.empty(nobj, dtype=np.object)

    x_all = np.array([])
    y_all = np.array([])
    z_all = np.array([])

    for i in range(nobj):

        objinfo = objinfo_all
        obj = objinfo['name'][i].decode('UTF-8').lower()      
        dir=FerrerDir

        tb = Table.read(dir+obj + '_allprop_halpha_sf.fits')
       
        z0=tb['lsmd']
        y0=tb['oh12_M13_O3N2']
        x0=tb['r_ell']*0.2/objinfo['re'][i]

        mask = np.isnan(x0)      
        mask = np.logical_or(mask, np.isnan(y0))
        mask = np.logical_or(mask, np.isnan(z0))     
        mask = np.logical_or(mask, tb['bpt_info'] !=1)
      
        x = x0[~mask]
        y = y0[~mask]
        z = z0[~mask]
        x_all = np.concatenate((x_all, x))
        y_all = np.concatenate((y_all, y))
        z_all = np.concatenate((z_all, z))

        med_2d_gal, y_edges_gal, x_edges_gal, binnum_gal \
        = scipy.stats.binned_statistic_2d(
            y, x, z,
            statistic='count',
            # statistic='median',
            # statistic='sum',
            bins=(60, 60),
            # range=[[-6,2], [0, 8]])
                range=[[7.6,9.6], [0,5.5]])
        # med_2d_gal=med_2d_gal*0+1
        mask5=med_2d_gal==0
        med_2d_gal[~mask5]=1
        med_2d_gal[mask5]=0
        if i==0:
            med_2d_tot=med_2d_gal*0
        med_2d_tot=med_2d_tot+med_2d_gal



    d_bin = 0.25 / 2.
    bin_min = 0 + d_bin
    bin_max = 6

    nbin = int((bin_max - bin_min + 2 * d_bin) / (2 * d_bin)) + 1
    bins = np.linspace(bin_min - d_bin, bin_max + d_bin, nbin)


    med_each, l68_each, u68_each, bin_edges, binnum \
        = get_stats(x, y, bins)

    bin_width = (bin_edges[1] - bin_edges[0])
    bin_centers = bin_edges[1:] - bin_width / 2.


    #Number galaxy plot
    med_2d_super, y_edges_super, x_edges_super, binnum_super \
        = scipy.stats.binned_statistic_2d(
            y_all, x_all, z_all,
            statistic='count',
            # statistic='median',
            bins=(60, 60),
            # range=[[-6,2], [0, 8]])
            range=[[7.6,9.6], [0,5.5]])
    mask6=med_2d_super<5
    med_2d_tot[mask6]=0
    mask7=med_2d_tot==0
    med_2d_tot[mask7]='NaN'


    props = dict(facecolor='white', edgecolor='none', alpha=0.5)
  

    # med_2d_tot=med_2d_tot-2
    med_2d_num, xedges_num, yedges_num = np.histogram2d(x_all, y_all, bins=(60, 60), range = [[0,5.5], [7.6,9.6]])
    x = .5 * (xedges_num[1:] + xedges_num[:-1])
    y = .5 * (yedges_num[1:] + yedges_num[:-1])


####################GALAXIES WITH logM<10

    objinfo_all = np.genfromtxt(FerrerDir+'table_info_mass_low.txt', dtype=None,names=['name', '2a', 'b_a', 'pa','d', 're', 'xc', 'yc', 'mass','xplot','yplot','npix'])
    nobj = objinfo_all['name'].size

    med_each_obj = np.empty(nobj, dtype=np.object)
    l68_each_obj = np.empty(nobj, dtype=np.object)
    u68_each_obj = np.empty(nobj, dtype=np.object)

    x_all = np.array([])
    y_all = np.array([])
    z_all = np.array([])

    for i in range(nobj):

        objinfo = objinfo_all
        obj = objinfo['name'][i].decode('UTF-8').lower()      
        dir=FerrerDir
        tb = Table.read(dir+obj+ '_allprop_halpha_sf.fits')
  
        z0=tb['lsmd']
        y0=tb['oh12_M13_O3N2']
        x0=tb['r_ell']*0.2/objinfo['re'][i]

        mask = np.isnan(x0)
        # mask = np.logical_or(mask, mask_param!=1)

        mask = np.logical_or(mask, np.isnan(y0))
        mask = np.logical_or(mask, np.isnan(z0))    
        mask = np.logical_or(mask, tb['bpt_info'] !=1)
        # mask = np.logical_or(mask, tb['dig_info'] !=1)

        x = x0[~mask]
        y = y0[~mask]
        z = z0[~mask]
        x_all = np.concatenate((x_all, x))
        y_all = np.concatenate((y_all, y))
        z_all = np.concatenate((z_all, z))

        med_2d_gal, y_edges_gal, x_edges_gal, binnum_gal \
        = scipy.stats.binned_statistic_2d(
            y, x, z,
            statistic='count',
            # statistic='median',
            # statistic='sum',
            bins=(60, 60),
            # range=[[-6,2], [0, 8]])
                range=[[7.6,9.6], [0,5.5]])
        # med_2d_gal=med_2d_gal*0+1
        mask5=med_2d_gal==0
        med_2d_gal[~mask5]=1
        med_2d_gal[mask5]=0
        if i==0:
            med_2d_tot=med_2d_gal*0
        med_2d_tot=med_2d_tot+med_2d_gal


    med_2d_super, y_edges_super, x_edges_super, binnum_super \
        = scipy.stats.binned_statistic_2d(
            y_all, x_all, z_all,
            statistic='count',
            # statistic='median',
            bins=(60, 60),
            # range=[[-6,2], [0, 8]])
            range=[[7.6,9.6], [0,5.5]])
    mask6=med_2d_super<10
    med_2d_tot[mask6]=0
    mask7=med_2d_tot==0
    med_2d_tot[mask7]='NaN'

    # cax1 = ax.imshow(med_2d_tot,
    #                   interpolation='nearest', origin='lower',
    #                   extent=[x_edges_gal[0], x_edges_gal[-1],
    #                           y_edges_gal[0], y_edges_gal[-1]],
    #                   vmin=0,
    #                   vmax=1,
    #                   cmap=cm.Spectral_r,
    #                   aspect='auto')


    med, l68, u68, bin_edges, binnum \
        = get_stats(x_all, y_all, bin_edges)


    ax.errorbar(bin_centers, med,
                 xerr=(bin_edges[1] - bin_edges[0]) / 2.,
                 yerr=[med - l68,
                       u68 - med],
                 fmt='s',
                 ecolor='red',
                 elinewidth=2.5,
                 capsize=0,
                 ms=12,
                 mew=1,
                 mfc='red',
                 mec='red',zorder=2)
                 # mec='red',label=ur'log[M$_\mathregular{*}$/M$_\odot$] < 10')
    print('M13 low mass')
    x_error = (bin_edges[1] - bin_edges[0]) / 2.
    print('radius,metallicity,err_x,err_y_down,err_y_up')
    for ii in range(0, len(bin_centers)):
        print('{:f},{:f},{:f},{:f},{:f}'.format(bin_centers[ii],med[ii],x_error, med[ii] - l68[ii], u68[ii] - med[ii]))
   
    
    med_2d_num, xedges_num, yedges_num = np.histogram2d(x_all, y_all, bins=(60, 60), range = [[0,5.5], [7.6,9.6]])
    x = .5 * (xedges_num[1:] + xedges_num[:-1])
    y = .5 * (yedges_num[1:] + yedges_num[:-1])

    # ax.contour(x, y, med_2d_num.T, levels = [5],  linestyles='-',colors = 'red', linewidths=3,label=ur'log[M$_\mathregular{*}$/M$_\odot$] < 10')  #, zorder = 10)

    ax.set_xlim(0,2.95)
    ax.set_ylim(7.6,9.9)


####################GALAXIES WITH 10<logM<11

    objinfo_all = np.genfromtxt(FerrerDir+'table_info_mass_med.txt', dtype=None,names=['name', '2a', 'b_a', 'pa','d', 're', 'xc', 'yc', 'mass','xplot','yplot','npix'])
    nobj = objinfo_all['name'].size

    # objinfo_all = np.genfromtxt('table_info_type'+type_num+bymass+'_mass.txt', dtype=None,
    #                             names=['name', '2a', 'b_a', 'pa',
    #                                    'd', 're', 'xc', 'yc','mass','sfr','rbreak'])
    # p.obj_ignore = 'PGC3853'

    nobj = objinfo_all['name'].size
    med_each_obj = np.empty(nobj, dtype=np.object)
    l68_each_obj = np.empty(nobj, dtype=np.object)
    u68_each_obj = np.empty(nobj, dtype=np.object)

    x_all = np.array([])
    y_all = np.array([])
    z_all = np.array([])

    for i in range(nobj):

        objinfo = objinfo_all
        obj = objinfo['name'][i].decode('UTF-8').lower()      
        dir=FerrerDir
        tb = Table.read(dir+obj + '_allprop_halpha_sf.fits')
        
        z0=tb['lsmd']
        y0=tb['oh12_M13_O3N2']
        x0=tb['r_ell']*0.2/objinfo['re'][i]

        mask = np.isnan(x0)
        # mask = np.logical_or(mask, mask_param!=1)

        mask = np.logical_or(mask, np.isnan(y0))
        mask = np.logical_or(mask, np.isnan(z0))
        mask = np.logical_or(mask, tb['bpt_info'] !=1)
        # mask = np.logical_or(mask, tb['dig_info'] !=1)

        x = x0[~mask]
        y = y0[~mask]
        z = z0[~mask]
        x_all = np.concatenate((x_all, x))
        y_all = np.concatenate((y_all, y))
        z_all = np.concatenate((z_all, z))

        med_2d_gal, y_edges_gal, x_edges_gal, binnum_gal \
        = scipy.stats.binned_statistic_2d(
            y, x, z,
            statistic='count',
            # statistic='median',
            # statistic='sum',
            bins=(60, 60),
            # range=[[-6,2], [0, 8]])
                range=[[7.6,9.6], [0,5.5]])
        # med_2d_gal=med_2d_gal*0+1
        mask5=med_2d_gal==0
        med_2d_gal[~mask5]=1
        med_2d_gal[mask5]=0
        if i==0:
            med_2d_tot=med_2d_gal*0
        med_2d_tot=med_2d_tot+med_2d_gal


    med, l68, u68, bin_edges, binnum \
        = get_stats(x_all, y_all, bin_edges)


    ax.errorbar(bin_centers, med,
                 xerr=(bin_edges[1] - bin_edges[0]) / 2.,
                 yerr=[med - l68,
                       u68 - med],
                 fmt='D',
                 ecolor='red',
                 elinewidth=2.5,
                 capsize=0,
                 ms=12,
                 mew=1,
                 mfc='red',
                 mec='red',zorder=2)
    print('M13 inter mass')
    x_error = (bin_edges[1] - bin_edges[0]) / 2.
    print('radius,metallicity,err_x,err_y_down,err_y_up')
    for ii in range(0, len(bin_centers)):
        print('{:f},{:f},{:f},{:f},{:f}'.format(bin_centers[ii],med[ii],x_error, med[ii] - l68[ii], u68[ii] - med[ii]))
       

    med_2d_num, xedges_num, yedges_num = np.histogram2d(x_all, y_all, bins=(60, 60), range = [[0,5.5], [7.6,9.6]])
    x = .5 * (xedges_num[1:] + xedges_num[:-1])
    y = .5 * (yedges_num[1:] + yedges_num[:-1])
    # ax.contour(x, y, med_2d_num.T, levels = [5],  linestyles='-',colors = 'red', linewidths=3,label='Intermediate')  #, zorder = 10)


####################GALAXIES WITH logM>11

    objinfo_all = np.genfromtxt(FerrerDir+'table_info_mass_high.txt', dtype=None,names=['name', '2a', 'b_a', 'pa','d', 're', 'xc', 'yc', 'mass','xplot','yplot','npix'])
    nobj = objinfo_all['name'].size

    # objinfo_all = np.genfromtxt('table_info_type'+type_num+bymass+'_mass.txt', dtype=None,
    #                             names=['name', '2a', 'b_a', 'pa',
    #                                    'd', 're', 'xc', 'yc','mass','sfr','rbreak'])
    # p.obj_ignore = 'PGC3853'

    nobj = objinfo_all['name'].size



    med_each_obj = np.empty(nobj, dtype=np.object)
    l68_each_obj = np.empty(nobj, dtype=np.object)
    u68_each_obj = np.empty(nobj, dtype=np.object)

    x_all = np.array([])
    y_all = np.array([])
    z_all = np.array([])

    for i in range(nobj):

        objinfo = objinfo_all
        obj = objinfo['name'][i].decode('UTF-8').lower() 
        dir=FerrerDir
        tb = Table.read(dir+obj+ '_allprop_halpha_sf.fits')

        z0=tb['lsmd']
        y0=tb['oh12_M13_O3N2']
        x0=tb['r_ell']*0.2/objinfo['re'][i]

        mask = np.isnan(x0)
        # mask = np.logical_or(mask, mask_param!=1)

        mask = np.logical_or(mask, np.isnan(y0))
        mask = np.logical_or(mask, np.isnan(z0))    
        mask = np.logical_or(mask, tb['bpt_info'] !=1)
        # mask = np.logical_or(mask, tb['dig_info'] !=1)

        x = x0[~mask]
        y = y0[~mask]
        z = z0[~mask]
        x_all = np.concatenate((x_all, x))
        y_all = np.concatenate((y_all, y))
        z_all = np.concatenate((z_all, z))

        med_2d_gal, y_edges_gal, x_edges_gal, binnum_gal \
        = scipy.stats.binned_statistic_2d(
            y, x, z,
            statistic='count',
            # statistic='median',
            # statistic='sum',
            bins=(60, 60),
            # range=[[-6,2], [0, 8]])
                range=[[7.6,9.6], [0,5.5]])
        # med_2d_gal=med_2d_gal*0+1
        mask5=med_2d_gal==0
        med_2d_gal[~mask5]=1
        med_2d_gal[mask5]=0
        if i==0:
            med_2d_tot=med_2d_gal*0
        med_2d_tot=med_2d_tot+med_2d_gal

    med, l68, u68, bin_edges, binnum \
        = get_stats(x_all, y_all, bin_edges)


    ax.errorbar(bin_centers, med, xerr=(bin_edges[1] - bin_edges[0]) / 2., yerr=[med - l68, u68 - med], 
                fmt='o', ecolor='red', elinewidth=2.5, capsize=0, ms=12, mew=1, mfc='red', mec='red',zorder=2)
    
    print('M13 high mass')
    x_error = (bin_edges[1] - bin_edges[0]) / 2.
    print('radius,metallicity,err_x,err_y_down,err_y_up')
    for ii in range(0, len(bin_centers)):
        print('{:f},{:f},{:f},{:f},{:f}'.format(bin_centers[ii],med[ii],x_error, med[ii] - l68[ii], u68[ii] - med[ii]))
   
    med_2d_num, xedges_num, yedges_num = np.histogram2d(x_all, y_all, bins=(60, 60), 
                                                        range = [[0,5.5], [7.6,9.6]])
    x = .5 * (xedges_num[1:] + xedges_num[:-1])
    y = .5 * (yedges_num[1:] + yedges_num[:-1])
    # ax.contour(x, y, med_2d_num.T, levels = [5],  linestyles='-',colors = 'red', linewidths=3,label='High mass')  #, zorder = 10)


    # FAKE LEGEND
    ax.errorbar(bin_centers, med+10,fmt='o', ecolor='black',elinewidth=2.5,capsize=0,ms=12,
                mew=1,mfc='black',mec='black',label='log[M$_\mathregular{*}$/M$_\odot$] > 10.8')
    ax.errorbar(bin_centers, med+10,fmt='D', ecolor='black',elinewidth=2.5,capsize=0,ms=12,
                mew=1,mfc='black',mec='black',label='10 < log[M$_\mathregular{*}$/M$_\odot$] < 10.8')
    ax.errorbar(bin_centers, med+10,fmt='s', ecolor='black',elinewidth=2.5,capsize=0,ms=12,
                mew=1,mfc='black',mec='black',label='log[M$_\mathregular{*}$/M$_\odot$] < 10')



    ax.legend(loc=1, scatterpoints=1, frameon=False, fontsize=14)

    ax.text(0.3,9.6,'DOP16',color='blue', fontsize=20)
    ax.text(0.3,9.4,'M13',color='red', fontsize=20)

    ax.set_xlim(0,2.95)
    ax.set_ylim(7.6,9.9)

    plt.savefig(outfile, bbox_inches='tight')



