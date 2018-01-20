#!/usr/bin/env python

import sys
import csv
import numpy as np
import matplotlib as mpl
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
    elif type == 'tmax':
        label = 'peak temperature'
    else:
        label = type
    return label
 
def get_limits(vmin=None, vmax=None, datavals=None, i=0):
    if vmin is None:
        v0 = np.floor(np.min(datavals[np.nonzero(datavals)]))
    else:
        if not isinstance(vmin, list): vmin = [vmin]
        v0 = vmin[i]
    if vmax is None:
        v1 = datavals.max()
    else:
        if not isinstance(vmax, list): vmax = [vmax]
        v1 = vmax[i]
    return v0, v1

#%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%&%&%&%&%&%&%
# Image structures color-coded by properties
#%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%&%&%&%&%&%&%

def props_colmap(dendrogram=None, subcat=None, img=None, cubhd=None, 
        props=['tmax'], vmin=None, vmax=None, cmapname='jet', prefix='output', **kwargs):
    print("Image leaves and clusters colored by properties")
    srclist = subcat['_idx'].tolist()
    # Make a plot for each requested property
    for i, type in enumerate(props):
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        imax = np.nanmax(img)
        im = ax.matshow(img, origin='lower', cmap=plt.cm.Blues, vmax=imax)
        if 'xlims' in kwargs:
            ax.set_xlim(kwargs['xlims'])
        if 'ylims' in kwargs:
            ax.set_ylim(kwargs['ylims'])

        plt.tick_params(axis='both', which='both', bottom='off', top='off', 
            labelbottom='off', right='off', left='off', labelleft='off',
            labeltop='off')
        name = scale_values(cat=subcat, type=type, cubhd=cubhd)
        v0, v1 = get_limits(vmin=vmin, vmax=vmax, datavals=subcat[type], i=i)
        print('{} vmin and vmax: {} {}'.format(type,v0,v1))
        cmap = plt.cm.get_cmap(cmapname)
        for i, c in enumerate(srclist):
            s = analysis.PPVStatistic(dendrogram[c])
            scaled_v = (subcat[type][i]-v0)/(v1-v0)
            col = cmap(scaled_v)
            ellipse = s.to_mpl_ellipse(color=col)
            ax.add_patch(ellipse)
        cax = fig.add_axes([0.92, 0.1, 0.03, 0.8])
        cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap,
             orientation='vertical', norm=mpl.colors.Normalize(vmin=v0, vmax=v1))
        cbar.ax.tick_params(labelsize=9)
        cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), rotation=90)
        cbar.set_label(name+' ['+str(subcat[type].unit)+']',size=12,labelpad=10)
        plt.savefig(prefix+'_'+type.replace("_","")+'.pdf', 
            bbox_inches='tight')
        plt.close()
    return

#%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%&%&%&%&%&%&%
# Draw tree diagram color-coded by properties
#%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%&%&%&%&%&%&%

def props_coltree(label=None, dendrogram=None, cat=None, cubhd=None, 
        props=['tmax'], cmapname='jet', vmin=None, vmax=None, prefix='output'):
    print("Draw tree diagram colored by properties")
    srclist = cat['_idx'].tolist()
    # Make a plot for each requested property
    for i, type in enumerate(props):
        fig = plt.figure(figsize=(8, 5))
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        ax.set_xlabel('Structure')
        ax.set_ylabel('Intensity [K]')
        ax.set_yscale('log')
        p = dendrogram.plotter()
        name = scale_values(cat=cat, type=type, cubhd=cubhd)
        v0, v1 = get_limits(vmin=vmin, vmax=vmax, datavals=cat[type], i=i)
        print('{} vmin and vmax: {} {}'.format(type,v0,v1))
        cmap = plt.cm.get_cmap(cmapname)
        for st in dendrogram.all_structures:
            scaled_v = (cat[type][st.idx]-v0)/(v1-v0)
            dcolr = cmap(scaled_v)
            p.plot_tree(ax, structure=[st], colors=dcolr, subtree=False)
        cax = fig.add_axes([0.93, 0.1, 0.02, 0.8])
        cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap,
             orientation='vertical', norm=mpl.colors.Normalize(vmin=v0, vmax=v1))
        cbar.ax.tick_params(labelsize=9)
        cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), rotation=90)
        cbar.set_label(name+' ['+str(cat[type].unit)+']',size=12,labelpad=10)
        if label == 'pcc_12':
            ax.annotate('N', xy=(63, 2.5), xytext=(40, 60),
                arrowprops=dict(facecolor='gray',width=3,headwidth=10,headlength=10,
                alpha=0.7), xycoords='data', textcoords='offset points')
            ax.annotate('S', xy=(111, 9), xytext=(-60, 20),
                arrowprops=dict(facecolor='gray',width=3,headwidth=10,headlength=10,
                alpha=0.7), xycoords='data', textcoords='offset points')
        plt.savefig(prefix+'_dendrogram_'+type.replace("_","")+'.pdf', 
            bbox_inches='tight')
        plt.close()
    return

# -------------------------------------------------------------------------------

def colorcode(label='scimes', table='full_catalog', cubefile=None, mom0file=None, 
        types=['v_cen','v_rms','tmax','imean'], outdir='plots', vmin=None, vmax=None,
        **kwargs):
    # Header info
    hdu3 = fits.open(cubefile)[0]
    hd3 = hdu3.header
    # Moment image
    hdu2 = fits.open(mom0file)[0]
    img = hdu2.data
    # Load the dendrogram
    d = Dendrogram.load_from(label+'_dendrogram.hdf5')
    print('\n')
    cat = Table.read(label+'_'+table+'.txt', format='ascii.ecsv')
    # Plot colored ellipses on maps
    for set in ['leaves', 'trunks', 'clusters']:
        idc = []
        with open(label+'_'+set+'.txt', 'r') as f:
            reader=csv.reader(f, delimiter=' ')
            for row in reader:
                idc.append(int(row[0]))
        subcat = cat[idc]
        props_colmap(dendrogram=d, subcat=subcat, img=img, cubhd=hd3,
            props=types, prefix=outdir+'/'+label+'_'+set, vmin=vmin, vmax=vmax, **kwargs)
    # Plot colored dendrogram
    props_coltree(label=label, dendrogram=d, cat=cat, cubhd=hd3, props=types, 
        prefix=outdir+'/'+label, vmin=vmin, vmax=vmax)
    return
