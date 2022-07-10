#!/usr/bin/env python

import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
from astropy.io import fits
import numpy as np
import matplotlib.colors as mcolors

# Plot the 12CO and 13CO cluster maps side by side.

fitsdir = 'images/'
analdir = 'struct/'
file0 = fitsdir+'30Dor_feather_mosaic_12CO_12meter.mom0.fits.gz'
hd0 = fits.getheader(file0)
wcs = WCS(hd0)
filer = fitsdir+'30Dor_feather_mosaic_12CO_12meter.rms.K.fits.gz'
rms = fits.getdata(filer)[0]
fov = rms*0 + 1
fov[np.isnan(fov)] = 0

#olays = np.array(wcs.world_to_pixel(SkyCoord(ra=[84.676650]*u.deg, dec=[-69.100933]*u.deg)))
olays = np.array(wcs.world_to_pixel(SkyCoord(ra=[84.67625, 84.5708]*u.deg, 
                                    dec=[-69.100917, -69.0667]*u.deg)))
skel = fits.getdata(analdir+'30Dor_feather_12CO_fp_80_st_80_at_10_gt_4.skel.fits.gz')

fig, axs = plt.subplots(ncols=2, figsize=(18,10), subplot_kw={"projection":wcs})

fig.subplots_adjust(right=0.8, wspace=0.02)

for i, line in enumerate(['12', '13']):
    clustasgn = fits.getdata(analdir+'30Dor_feather_mosaic_1p8_'
                             +line+'_clusters_asgn.fits.gz')
    clustlist = analdir+'30Dor_feather_mosaic_1p8_'+line+'_clusters.txt'
    with open(clustlist, 'r') as data:
        clust_id = []
        clust_co = []
        for line in data:
            p = line.split()
            clust_id.append(int(p[0]))
            clust_co.append((p[1]))
    print('Number of clusters:',len(clust_id))
    # Clusters
    for j, idnum in enumerate(clust_id):
        mask3d = (clustasgn == idnum).astype(int)
        mask2d = np.amax(mask3d, axis=0)
        if i == 0:
            axs[i].contourf(mask2d, colors=clust_co[j], levels=[0.5,1.5], alpha=0.6)
        elif i == 1:
            axs[i].contourf(mask2d, colors=clust_co[j], levels=[0.5,1.5], alpha=0.7)
            axs[i].contour(mask2d, colors='silver', levels=[1])
    # FOV contour
    axs[i].contour(fov, [1], colors='k', linestyles='dotted')
    # FilFinder skeleton
    if i == 0:
        axs[i].contour(skel, levels=[1], colors='k', linewidths=1)
    # Star cluster positions
    print('R136 is at',olays[0][0], olays[1][0])
    print('Hodge 301 is at',olays[0][1], olays[1][1])
    axs[i].plot(olays[0], olays[1], 'b*', ms=12, zorder=-1)
    axs[i].annotate('R136', (olays[0][0], olays[1][0]), textcoords="offset points",
                 xytext=(-5,-20), ha='right', c='b', size=13)
    axs[i].annotate('Hodge 301', (olays[0][1], olays[1][1]), textcoords="offset points",
                 xytext=(-5,-20), ha='center', c='b', size=13)
    if i == 1:
        axs[i].coords[1].set_ticklabel_visible(False)
    axs[i].coords[0].set_ticklabel(size=14)
    axs[i].coords[1].set_ticklabel(size=14)
    axs[i].set_xlim(-0.5,799.5)
    axs[i].set_ylim(-0.5,999.5)
    axs[i].set_xlabel('Right Ascension (J2000)', size=15)
    if i == 0:
        axs[i].text(0.03,0.95,'$^{12}$CO clumps & filaments',ha='left',va='center',
                 fontsize=18,transform=axs[i].transAxes, 
                    bbox=dict(facecolor='white', edgecolor='none'))
        axs[i].set_ylabel('Declination (J2000)', size=15)
    else:
        axs[i].text(0.03,0.95,'$^{13}$CO clumps',ha='left',va='center',
                 fontsize=18,transform=axs[i].transAxes)
        axs[i].set_ylabel(' ')
    axs[i].tick_params(direction='in')

plt.savefig('30Dor_12CO_13CO_clust.pdf', bbox_inches='tight')




