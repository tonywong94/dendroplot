#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Ellipse
#from astropy.visualization import ImageNormalize


def mom0plot(mom0file=None, fluxfile=None, cmap='hot_r', distpc=5e4, 
        figsize=[6,6], labelax='full', xrange=[None,None], yrange=[None,None], 
        v0=0., v1=None, cmin=0, fov=0.5, cpad=0.1, label=None, outfile=None):
    # --- Set up plot window
    hdu = fits.open(mom0file)[0]
    if hdu.header['NAXIS'] == 2:
        for key in ['PC03_01', 'PC03_02', 'PC01_03', 'PC02_03', 'PC03_03',
                    'PC3_1', 'PC3_2', 'PC1_3', 'PC2_3', 'PC3_3', 'CTYPE3', 
                    'CRVAL3', 'CDELT3','CRPIX3', 'CUNIT3', 'NAXIS3']:
            if key in hdu.header.keys():
                hdu.header.remove(key)
    bmaj = hdu.header['BMAJ']*3600.
    bmin = hdu.header['BMIN']*3600.
    if 'BPA' in hdu.header.keys():
        bpa = hdu.header['BPA']
    else:
        bpa = 0.
    cdel = np.abs(hdu.header['CDELT2']*3600.)
    wcs = WCS(hdu.header)
    fig = plt.figure(figsize=(figsize[0], figsize[1]))
    print("Image dimensions in pixels:", hdu.data.shape)
    print("  Beam width in pixels: {:.2f} {:.2f}".format(bmaj/cdel, bmin/cdel))
    # --- Select the subregion
    datcut = hdu.data[slice(yrange[0],yrange[1],1),slice(xrange[0],xrange[1],1)]
    wcscut = wcs[slice(yrange[0],yrange[1],1),slice(xrange[0],xrange[1],1)]
    ax = plt.subplot(projection=wcscut)
    # Convert to K km/s if necessary
    if hdu.header['BUNIT'].lower() == 'jy/beam.km/s':
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
        datcut *= convfac
    # --- Plot image
    if v1 is None:
        v1 = np.nanmax(datcut)
        print("  Image range: {} to {}".format(v0,v1))
    im = ax.imshow(datcut, origin='lower', cmap=cmap, vmin=v0, vmax=v1)
    clevs = np.logspace(cmin, cmin+6, 7, base=2)
    ax.contour(datcut, levels=clevs, colors='blue', alpha=0.7, linewidths=0.7)
    print("  Contour levels: {}".format(clevs))
    # --- Plot gain
    hdug=fits.open(fluxfile)[0]
    if hdug.header['NAXIS'] == 2:
        gaincut = hdug.data[slice(yrange[0],yrange[1],1),slice(xrange[0],xrange[1],1)]
        for key in ['PC03_01', 'PC03_02', 'PC01_03', 'PC02_03', 'PC03_03',
                    'PC3_1', 'PC3_2', 'PC1_3', 'PC2_3', 'PC3_3', 'CTYPE3', 
                    'CRVAL3', 'CDELT3','CRPIX3', 'CUNIT3', 'NAXIS3']:
            if key in hdug.header.keys():
                hdug.header.remove(key)
    elif hdug.header['NAXIS'] == 3:
        gaincut = hdug.data[0][slice(xrange[0],xrange[1],1),slice(yrange[0],yrange[1],1)]
    ax.contour(gaincut, transform=ax.get_transform(wcscut),
        levels=[fov], colors='red', alpha=0.5, linewidths=2, linestyles='dashed')
    # --- Plot beam and parsec scale
    axis_to_data = ax.transAxes + ax.transData.inverted()
    beamx, beamy = axis_to_data.transform([0.04,0.04])
    beam = Ellipse(xy=(beamx,beamy), width=bmaj/cdel, height=bmin/cdel, angle=90+bpa, 
                edgecolor=None, facecolor='blue')
    beampc = bmaj*distpc/206265.
    #beamstr = f'{beampc:.1f} pc'
    ax.add_patch(beam)
    #ax.annotate(beamstr, (0.1,0.025), xycoords='axes fraction', ha='center')
    tenpc = (10/distpc)*206265
    barx, bary = axis_to_data.transform([0.05,0.97])
    plt.hlines(bary, barx, barx+tenpc/cdel, color='k', alpha=0.8, lw=3)
    midpt = barx+0.5*(tenpc/cdel)
    ax.annotate('10 pc', (midpt,bary-8), ha='center', va='top', fontsize=11)
    # --- Plot labels
    lon=ax.coords[0]
    lat=ax.coords[1]
    if labelax == 'empty':
        lon.set_ticks_visible(False)
        lon.set_ticklabel_visible(False)
        lat.set_ticks_visible(False)
        lat.set_ticklabel_visible(False)
    else:
        lon.set_major_formatter('hh:mm:ss')
        #lon.set_ticks(spacing=ra_tick * u.hourangle/3600.)
        lon.set_ticks_position('b')
        lon.set_ticklabel(size=10, exclude_overlapping=True)
        lat.set_major_formatter('dd:mm')
        #lat.set_major_formatter('dd:mm:ss')
        lat.set_ticks_position('l')
        #lat.set_ticks(spacing=dc_tick * u.degree/60.)
        lat.set_ticklabel(size=10)
        if labelax != 'tickonly':
            lon.set_axislabel('Right Ascension (J2000)', size=11)
            lat.set_axislabel('Declination (J2000)', size=11)
        else:
            lon.set_auto_axislabel(False)
            lat.set_auto_axislabel(False)
    if label is not None:
        ax.text(0.97,0.97,label,ha='right',va='top',fontsize=13,
            transform=ax.transAxes)
    # --- Plot colorbar
#     divider = make_axes_locatable(ax)
#     cax = divider.append_axes("right", size="4%", pad=cpad)
#     cbar = plt.colorbar(im, cax=cax, orientation='vertical')
#     cbar.ax.tick_params(labelsize=10)
#     cbar.ax.coords[0].set_ticks_visible(False)
#     cbar.ax.coords[0].set_ticklabel_visible(False)
#     cbar.ax.coords[0].grid(False)
#     cbar.ax.coords[1].grid(False)
#     cbar.ax.coords[1].set_ticklabel_position('r')
#     cbar.ax.coords[1].set_ticks_position('r')
#     cbar.ax.coords[1].set_axislabel_position('r')
#     if labelax == 'full':
#         cbar.set_label('Integrated Intensity [K km/s]',fontsize=12)
    # --- Save PDF file
    if outfile is None:
        outfile = 'mom0plot.pdf'
    plt.savefig(outfile, bbox_inches='tight')
    plt.close()
    return
