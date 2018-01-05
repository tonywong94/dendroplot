#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u

def mom0plot(mom0file=None, fluxfile=None, cmap='hot_r', figsize=[6,6],
        labelax='none', xoff=[-60,60], yoff=[-60,60], v0=0., v1=None,
        cbar_tick=5., ra_tick=8., cmin=0, fov=0.5, label=None, outfile=None):
    # --- Set up plot window
    hdu = fits.open(mom0file)[0]
    if hdu.header['NAXIS'] == 2:
        for key in ['PC03_01', 'PC03_02', 'PC01_03', 'PC02_03', 'PC03_03', 'CTYPE3', 
                'CRVAL3', 'CDELT3','CRPIX3', 'CUNIT3', 'NAXIS3']:
            if key in hdu.header.keys():
                hdu.header.remove(key)
    wcs = WCS(hdu.header)
    fig = plt.figure(figsize=(figsize[0], figsize[1]))
    ax = plt.subplot(projection=wcs)
    lon=ax.coords[0]
    lat=ax.coords[1]
    # Convert to K km/s if necessary
    if hdu.header['BUNIT'].lower() == 'jy/beam.km/s':
        bmaj = hdu.header['BMAJ']*3600.
        bmin = hdu.header['BMIN']*3600.
        if 'RESTFREQ' in hdu.header.keys():
            freq = hdu.header['RESTFREQ'] * u.Hz
        elif 'RESTFRQ' in hdu.header.keys():
            freq = hdu.header['RESTFRQ'] * u.Hz
        omega_B = np.pi/(4*np.log(2)) * bmaj * bmin * u.arcsec**2
        convfac = (u.Jy).to(u.K, equivalencies=u.brightness_temperature(omega_B,freq))
        print("Units of the input cube are {}".format(hdu.header['BUNIT']))
        print("The beam size is {0} x {1} arcsec".format(bmaj,bmin))
        print("Scaling data by {0} to convert to K".format(convfac))
        hdu.header['bunit'] = 'K'
        hdu.data = hdu.data*convfac
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
    if hdug.header['NAXIS'] == 2:
        gaindata = hdug.data
        for key in ['PC03_01', 'PC03_02', 'PC01_03', 'PC02_03', 'PC03_03', 'CTYPE3', 
                'CRVAL3', 'CDELT3','CRPIX3', 'CUNIT3', 'NAXIS3']:
            if key in hdug.header.keys():
                hdug.header.remove(key)
    elif hdug.header['NAXIS'] == 3:
        gaindata = hdug.data[0]
    ax.contour(gaindata, transform=ax.get_transform(WCS(hdu.header)),
        levels=[fov], colors='red', alpha=0.5, linewidths=1, linestyles='dashed')
    # --- Plot labels
    if labelax == 'empty':
        lon.set_ticks_visible(False)
        lon.set_ticklabel_visible(False)
        lat.set_ticks_visible(False)
        lat.set_ticklabel_visible(False)
    else:
        lon.set_major_formatter('hh:mm:ss')
        lon.set_ticks(spacing=ra_tick * u.hourangle/3600., exclude_overlapping=True)
        lon.set_ticks_position('b')
        lon.set_ticklabel(size=10)
        lat.set_major_formatter('dd:mm')
        lat.set_ticks_position('l')
        lat.set_ticklabel(size=10)
        if labelax == 'full':
            lon.set_axislabel('Right Ascension (J2000)', size=11)
            lat.set_axislabel('Declination (J2000)', size=11)
    if label is not None:
        ax.text(0.5,1.08,label,ha='center',va='top',fontsize=11,
            transform=ax.transAxes)
    # --- Plot colorbar
    #ax2 = plt.gca()
    #divider = make_axes_locatable(ax2)
    #cax = divider.append_axes("right", size="5%", pad=0.04)
    #cbar = plt.colorbar(im, cax=cax, orientation='vertical')
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
