#!/usr/bin/env python

import sys
import csv
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from astrodendro import Dendrogram, ppv_catalog, analysis
from astropy import units as u
from astropy.io import fits
from astropy.table import Table, Column

#%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%&%&%&%&%&%&%
# Image structures color-coded by properties
#%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%&%&%&%&%&%&%

def props_colmap(dendrogram=None, subcat=None, img=None, cubhd=None, 
        props=['tmax'], vmin=None, vmax=None, cmapname='nipy_spectral_r',
        prefix='output'):

    print("Image leaves and clusters colored by properties")
    srclist = subcat['_idx'].tolist()
    # Make a plot for each requested property
    for type in props:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        imax = np.nanmax(img)
        im = ax.matshow(img, origin='lower', cmap=plt.cm.Blues, vmax=imax)
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.setp(ax.get_yticklabels(), visible=False)
        if type == 'flux':  # Convert Jy*ch to K km/s
            deltav = hd3['cdelt3']/1000. * u.km / u.s
            as2 = 1 * u.arcsec**2
            kflux = deltav * as2 * subcat[type].to(u.K,
                equivalencies=u.brightness_temperature(as2,freq))
            subcat[type] = kflux/subcat['area_exact']
            subcat[type].unit = 'K km / s'
        vmin = np.floor(subcat[type].min())
        vmax =  np.ceil(subcat[type].max())
        print(type,' vmin and vmax:',vmin,vmax)
        cmap = plt.cm.get_cmap(cmapname)
        for i, c in enumerate(srclist):
            s = analysis.PPVStatistic(dendrogram[c])
            scaled_v = (subcat[type][i]-vmin)/(vmax-vmin)
            col = cmap(scaled_v)
            ellipse = s.to_mpl_ellipse(color=col)
            ax.add_patch(ellipse)
        ax1 = fig.add_axes([0.3, 0.95, 0.45, 0.02])
        if type == 'v_cen':  # convert ch number to km/s
            v0 = 1.e-3*(cubhd['crval3']+cubhd['cdelt3']*(vmin-cubhd['crpix3']))
            v1 = 1.e-3*(cubhd['crval3']+cubhd['cdelt3']*(vmax-cubhd['crpix3']))
            subcat[type].unit = 'km / s'
        else:
            v0 = vmin
            v1 = vmax
        cbar = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
             orientation='horizontal',norm=mpl.colors.Normalize(vmin=v0, vmax=v1))
        cbar.ax.tick_params(labelsize=9)
        if type == 'flux':
            name = 'mean intensity'
        elif type == 'tmax':
            name = 'peak temperature'
        elif type == 'v_cen':
            name = 'mean velocity'
        elif type == 'v_rms':
            name = 'velocity dispersion'
        else:
            name = type
        cbar.set_label(name+' ['+str(subcat[type].unit)+']',size=9,labelpad=-35)
        plt.savefig(prefix+'_'+type.replace("_","")+'.pdf', 
            bbox_inches='tight')
        plt.close()
    return

#%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%&%&%&%&%&%&%
# Draw tree diagram color-coded by properties
#%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%&%&%&%&%&%&%

def props_coltree(dendrogram=None, cat=None, cubhd=None, props=['tmax'],
        cmapname='jet', vmin=None, vmax=None, prefix='output'):
    srclist = cat['_idx'].tolist()
    # Make a plot for each requested property
    for i, type in enumerate(props):
        fig = plt.figure(figsize=(8, 5))
        ax = fig.add_subplot(111)            
        ax.set_xlabel('Structure')
        ax.set_ylabel('Intensity')
        p = dendrogram.plotter()
        if vmin is None:
            v0 = np.floor(0.9*np.min(cat[type][np.nonzero(cat[type])]))
        else:
            if not isinstance(vmin, list): vmin = [vmin]
            v0 = vmin[i]
        if vmax is None:
            v1 = np.ceil(cat[type].max())
        else:
            if not isinstance(vmax, list): vmax = [vmax]
            v1 = vmax[i]
        print(type,' vmin and vmax:',v0,v1)
        cmap = plt.cm.get_cmap(cmapname)
        for st in dendrogram.all_structures:
            scaled_v = (cat[type][st.idx]-v0)/(v1-v0)
            dcolr = cmap(scaled_v)
            p.plot_tree(ax, structure=[st], colors=dcolr, subtree=False)
        ax1 = fig.add_axes([0.13, 0.95, 0.77, 0.03])
        if type == 'v_cen':  # convert ch number to km/s
            v0 = 1.e-3*(cubhd['crval3']+cubhd['cdelt3']*(v0-cubhd['crpix3']))
            v1 = 1.e-3*(cubhd['crval3']+cubhd['cdelt3']*(v1-cubhd['crpix3']))
            cat[type].unit = 'km / s'
        cbar = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
             orientation='horizontal',norm=mpl.colors.Normalize(vmin=v0, vmax=v1))
        cbar.ax.tick_params(labelsize=9)
        if type == 'tmax':
            name = 'peak temperature'
        elif type == 'v_cen':
            name = 'mean velocity'
        elif type == 'v_rms':
            name = 'velocity dispersion'
        else:
            name = type
        cbar.set_label(name+' ['+str(cat[type].unit)+']',size=9,labelpad=-35)
        plt.savefig(prefix+'_dendrogram_'+type.replace("_","")+'.pdf', 
            bbox_inches='tight')
        plt.close()
    return

# -------------------------------------------------------------------------------

def colorcode(label='scimes', table='full_catalog', cubefile=None, mom0file=None, 
        types=['v_cen','v_rms','tmax'], outdir='plots', vmin=None, vmax=None):
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
    for set in ['leaves', 'clusters']:
        with open(label+'_'+set+'.txt', 'r') as f:
            reader=csv.reader(f, delimiter=' ')
            a = zip(*reader)
        idc = map(int,a[0])
        subcat = cat[idc]
        props_colmap(dendrogram=d, subcat=subcat, img=img, cubhd=hd3,
            props=types, prefix=outdir+'/'+label+'_'+set, vmin=vmin, vmax=vmax)
    # Plot colored dendrogram
    props_coltree(dendrogram=d, cat=cat, cubhd=hd3, props=types, 
        prefix=outdir+'/'+label, vmin=vmin, vmax=vmax)
    return
