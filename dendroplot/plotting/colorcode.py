#!/usr/bin/env python

# Produce color-coded maps of structure locations, and color-coded tree diagrams.

import sys
import csv
import numpy as np
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import pyplot as plt
from astrodendro import Dendrogram, analysis
from astropy import units as u
from astropy.io import fits
from astropy.table import Table, Column

def scale_values(cat=None, type=None, cubhd=None):
    if type == 'imean':  # Convert Jy*ch to K km/s
        if 'RESTFREQ' in cubhd.keys():
            freq = cubhd['RESTFREQ'] * u.Hz
        elif 'RESTFRQ' in cubhd.keys():
            freq = cubhd['RESTFRQ'] * u.Hz
        deltav = cubhd['cdelt3']/1000. * u.km / u.s
        as2 = 1 * u.arcsec**2
        kflux = deltav * as2 * cat['flux'].to(u.K,
            equivalencies=u.brightness_temperature(as2,freq))
        cat[type] = kflux/cat['area_exact']
        cat[type].unit = 'K km / s'
        label = 'mean intensity'
    elif type == 'v_rms':  # Convert m/s to km/s
        cat[type] = cat['v_rms'].to(u.km/u.s)
        cat[type].unit = 'km / s'
        label = 'velocity dispersion'
    elif type == 'v_cen':  # convert ch number to km/s
        cat[type] = 1.e-3*(cubhd['crval3']+cubhd['cdelt3']*
            (cat['v_cen']-cubhd['crpix3']))
        cat[type].unit = 'km / s'
        label = 'mean velocity'
    elif type == 'vrms_k':
        label = 'velocity dispersion'
    elif type == 'tmax':
        label = 'peak temperature'
    else:
        label = type
    return label
 
def get_limits(vmin=None, vmax=None, datavals=None, lognorm=False, i=0):
    if vmin is None:
        v0 = np.min(datavals[np.isfinite(datavals)])
    else:
        if not isinstance(vmin, list): vmin = [vmin]
        v0 = vmin[i]
    if vmax is None:
        v1 = np.max(datavals[np.isfinite(datavals)])
        if not v1 > v0:
            v0 = np.floor(v0)
            v1 = np.ceil(v1)
    else:
        if not isinstance(vmax, list): vmax = [vmax]
        v1 = vmax[i]
    # Choose the ticks
    print('v0, v1:', v0, v1)
    if lognorm:
        if (v0<=0): v0 = 1e-2
        if (v1<=0): v1 = 1.
        tick0 = np.floor(np.log10(v0))
        tick1 = np.ceil(np.log10(v1))
        pow = np.arange(tick0, tick1, 1)
        ticks=[]
        for j in pow:
            for tscale in [1, 2, 5]:
                if v0 <= tscale*10**j and v1 >= tscale*10**j:
                    ticks.extend([tscale*10**j])
        ticklab = ["%g" % val for val in ticks]
    else:
        dex = 10**np.floor(np.log10(abs(v1-v0)))
        if abs(v1-v0)/dex > 5:
            dt = 2*dex
        else:
            dt = dex
        tick0 = dt * np.ceil(v0/dt)
        tick1 = dt * np.ceil(v1/dt)
        ticks = np.arange(tick0, tick1, dt)
        if dt < 1:
            ticklab = ["%.1f" % val for val in ticks]
        else:
            ticklab = ["%d" % val for val in ticks]
    print('Ticks:', ticklab)
    return v0, v1, ticks, ticklab

#%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%&%&%&%&%&%&%&%
# Image structures color-coded by properties
#%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%&%&%&%&%&%&%&%

def props_colmap(dendrogram=None, subcat=None, img=None, cubhd=None, 
        props=['tmax'], vmin=None, vmax=None, lognorm=False,
        cmapname='jet', prefix='output', noshow=False, **kwargs):
    print("Image leaves and clusters colored by properties")
    srclist = subcat['_idx'].tolist()
    # Make a plot for each requested property
    for i, type in enumerate(props):
        print('\nWorking on', type)
        fig = plt.figure(figsize=(8, 8))
        # --- Plot the image as background
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        imax = np.nanmax(img)
        im = ax.matshow(img, origin='lower', cmap=plt.cm.Greys, vmax=imax)
        if 'xlims' in kwargs:
            ax.set_xlim(kwargs['xlims'])
        if 'ylims' in kwargs:
            ax.set_ylim(kwargs['ylims'])
        name = scale_values(cat=subcat, type=type, cubhd=cubhd)
        datavals=subcat[type]
        if lognorm:
            datavals[datavals<=0] = np.nan
        if (~np.isfinite(datavals)).all():
            print('All data for {} are NaN'.format(type))
            plt.close()
            continue
        plt.tick_params(axis='both', which='both', bottom=False, top=False, 
            left=False, labelleft=False, labeltop=False)
        v0, v1, ticks, tlbl = get_limits(vmin=vmin, vmax=vmax, datavals=datavals, 
                                        lognorm=lognorm, i=i)
        print('{} vmin and vmax: {} {}'.format(type,v0,v1))
        cmap = plt.cm.get_cmap(cmapname)
        for i, c in enumerate(srclist):
            s = analysis.PPVStatistic(dendrogram[c])
            if lognorm:
                cnorm = mpl.colors.LogNorm(vmin=v0, vmax=v1)
            else:
                cnorm = mpl.colors.Normalize(vmin=v0, vmax=v1)
            scaled_v = cnorm(subcat[type][i])
            ellipse = s.to_mpl_ellipse(color=cmap(scaled_v))
            ax.add_patch(ellipse)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.15)
        cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap,
             orientation='vertical', norm=cnorm)
        cbar.ax.tick_params(labelsize=9)
        cbar.set_ticks(ticks)
        cbar.ax.set_yticklabels(tlbl, rotation=0, va='center')
        if str(subcat[type].unit) == '':
            cbar.set_label(name,size=12,labelpad=10)
        else:
            cbar.set_label(name+' ['+str(subcat[type].unit)+']',size=12,labelpad=10)
        plt.savefig(prefix+'_'+type.replace("_","")+'.pdf', 
            bbox_inches='tight')
        if noshow:
            plt.close(fig)
    return

#%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%&%&%&%&%&%&%&%
# Draw tree diagram color-coded by properties
#%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%&%&%&%&%&%&%&%

def props_coltree(label=None, dendrogram=None, cat=None, cubhd=None, 
        props=['tmax'], cmapname='jet', vmin=None, vmax=None, lognorm=False,
        xmin=None, xmax=None, prefix='output', noshow=False):
    print("Draw tree diagram colored by properties")
    # Make a plot for each requested property
    for i, type in enumerate(props):
        fig = plt.figure(figsize=(8, 5))
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        ax.set_xlabel('Structure')
        ax.set_ylabel('Intensity [K]')
        if xmin is not None and xmax is not None:
            ax.set_xlim(xmin,xmax)
        p = dendrogram.plotter()
        name = scale_values(cat=cat, type=type, cubhd=cubhd)
        datavals = cat[type]
        if lognorm:
            datavals[datavals<=0] = np.nan
        if (~np.isfinite(datavals)).all():
            print('All data for {} are NaN'.format(type))
            plt.close()
            continue
        v0, v1, ticks, tlbl = get_limits(vmin=vmin, vmax=vmax, datavals=datavals,
                                        lognorm=lognorm, i=i)
        print('{} vmin and vmax: {} {}'.format(type,v0,v1))
        cmap = plt.cm.get_cmap(cmapname)
        for st in dendrogram.all_structures:
            if lognorm:
                cnorm = mpl.colors.LogNorm(vmin=v0, vmax=v1)
            else:
                cnorm = mpl.colors.Normalize(vmin=v0, vmax=v1)
            scaled_v = cnorm(cat[type][st.idx])
            p.plot_tree(ax, structure=[st], colors=cmap(scaled_v), subtree=False)
        ax.set_yscale('log')
        cax = fig.add_axes([0.93, 0.1, 0.02, 0.8])
        cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap,
             orientation='vertical', norm=cnorm)
        cbar.ax.tick_params(labelsize=9)
        cbar.set_ticks(ticks)
        cbar.ax.set_yticklabels(tlbl, rotation=0, va='center')
        if str(cat[type].unit) == '':
            cbar.set_label(name,size=12,labelpad=10)
        else:
            cbar.set_label(name+' ['+str(cat[type].unit)+']',size=12,labelpad=15)
        if label == 'pcc_12':
            ax.annotate('N', xy=(63, 2.5), xytext=(40, 60),
                arrowprops=dict(facecolor='gray',width=3,headwidth=10,headlength=10,
                alpha=0.7), xycoords='data', textcoords='offset points')
            ax.annotate('S', xy=(111, 9), xytext=(-60, 20),
                arrowprops=dict(facecolor='gray',width=3,headwidth=10,headlength=10,
                alpha=0.7), xycoords='data', textcoords='offset points')
        plt.savefig(prefix+'_dendrogram_'+type.replace("_","")+'.pdf', 
            bbox_inches='tight')
        if noshow:
            plt.close(fig)
    return

# -------------------------------------------------------------------------------

def colorcode(label='scimes', table='full_catalog', cubefile=None, mom0file=None, 
        types=['v_cen','v_rms','tmax','imean'], plotdir='plots', vmin=None, vmax=None,
        lognorm=False, noshow=False, **kwargs):
    # Header info
    hdu3 = fits.open(cubefile)[0]
    hd3 = hdu3.header
    # Moment image
    hdu2 = fits.open(mom0file)[0]
    img = hdu2.data
    # Load the dendrogram
    d = Dendrogram.load_from(label+'_dendrogram.hdf5')
    print('\nOpening '+label+'_'+table+'.txt')
    cat = Table.read(label+'_'+table+'.txt', format='ascii.ecsv')
    # Plot colored ellipses on maps
    for set in ['leaves', 'trunks', 'clusters']:
        print('Working on', set)
        try:
            idc = []
            with open(label+'_'+set+'.txt', 'r') as f:
                reader=csv.reader(f, delimiter=' ')
                for row in reader:
                    idc.append(int(row[0]))
            subcat = cat[idc]
            props_colmap(dendrogram=d, subcat=subcat, img=img, cubhd=hd3,
                props=types, prefix=plotdir+'/'+label+'_'+set, vmin=vmin, vmax=vmax, 
                lognorm=lognorm, noshow=noshow, **kwargs)
        except IOError:
            print(label,set,'not found')
    # Plot colored dendrogram
    props_coltree(label=label, dendrogram=d, cat=cat, cubhd=hd3, props=types, 
        prefix=plotdir+'/'+label, vmin=vmin, vmax=vmax, lognorm=lognorm, noshow=noshow)
    return
