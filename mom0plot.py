#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u

def mom0plot(mom0file=None, fluxfile=None, cmap='hot_r', figsize=[6,6],
        labelax='none', xoff=[-60,60], yoff=[-60,60], v0=0., v1=None,
        cbar_tick=5., ra_tick=8., cmin=0, label=None, outfile=None):
    # --- Set up plot window
    hdu = fits.open(mom0file)[0]
    hdu.header.remove('PC03_01')
    hdu.header.remove('PC03_02')
    hdu.header.remove('PC01_03')
    hdu.header.remove('PC02_03')
    hdu.header.remove('PC03_03')
    wcs = WCS(hdu.header)
    fig = plt.figure(figsize=(figsize[0], figsize[1]))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=wcs)
    lon=ax.coords[0]
    lat=ax.coords[1]
    # --- Plot subimage centered on reference pixel
    xminpix = max(0., wcs.wcs.crpix[0] + xoff[0]/abs(3600*wcs.wcs.cdelt[0]) - 0.5)
    xmaxpix = min(hdu.data.shape[1], wcs.wcs.crpix[0] + xoff[1]/abs(3600*wcs.wcs.cdelt[0]) + 0.5)
    yminpix = max(0., wcs.wcs.crpix[1] + yoff[0]/abs(3600*wcs.wcs.cdelt[1]) - 0.5)
    ymaxpix = min(hdu.data.shape[0], wcs.wcs.crpix[1] + yoff[1]/abs(3600*wcs.wcs.cdelt[1]) + 0.5)
    print("xmin, xmax, ymin, ymax: {} {} {} {}".format(xminpix,xmaxpix,yminpix,ymaxpix))
    ax.set_xlim(xminpix, xmaxpix)
    ax.set_ylim(yminpix, ymaxpix)
    if v1 is None:
        v1 = np.nanmax(hdu.data)
    im = ax.imshow(hdu.data, origin='lower', vmin=v0, vmax=v1, cmap=cmap)
    clevs = np.logspace(cmin, cmin+6, 7, base=2)
    ax.contour(hdu.data, levels=clevs, colors='blue', alpha=0.7, linewidths=0.7)
    print("Contour levels: {}".format(clevs))
    # --- Plot gain
    hdug=fits.open(fluxfile)[0]
    hdug.header.remove('PC03_01')
    hdug.header.remove('PC03_02')
    hdug.header.remove('PC01_03')
    hdug.header.remove('PC02_03')
    hdug.header.remove('PC03_03')
    hdug.header.remove('CTYPE3')
    hdug.header.remove('CRVAL3')
    hdug.header.remove('CDELT3')
    hdug.header.remove('CRPIX3')
    hdug.header.remove('CUNIT3')
    ax.contour(hdug.data, transform=ax.get_transform(WCS(hdug.header)),
        levels=[0.6], colors='red', alpha=0.5, linewidths=1, linestyles='dashed')
    # --- Plot labels
    if labelax == 'empty':
        lon.set_ticks_visible(False)
        lon.set_ticklabel_visible(False)
        lat.set_ticks_visible(False)
        lat.set_ticklabel_visible(False)
    else:
        lon.set_major_formatter('hh:mm:ss')
        lat.set_major_formatter('dd:mm')
        lon.set_ticklabel(size=10)
        lat.set_ticklabel(size=10)
        lon.set_ticks(spacing=ra_tick * u.hourangle/3600.)
        if labelax == 'full':
            lon.set_axislabel('Right Ascension (J2000)', size=11)
            lat.set_axislabel('Declination (J2000)', size=11)
    if label is not None:
        ax.text(0.5,1.05,label,ha='center',va='top',fontsize=11,
            transform=ax.transAxes)
    # --- Plot colorbar
    cbar = plt.colorbar(im, aspect=30, pad=0.03)
    cbar.ax.tick_params(labelsize=10)
    if labelax != 'empty':
        cbar.set_label('Integrated Intensity [K km/s]',fontsize=11)
    # --- Save PDF file
    if outfile is None:
        outfile = 'mom0plot.pdf'
    plt.savefig(outfile, bbox_inches='tight')
    plt.close()
    return
