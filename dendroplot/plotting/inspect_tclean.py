#!/usr/bin/env python3

import re
import os.path
import numpy as np
import astropy.io.fits as pyfits
import astropy.units as u
from astropy.stats import mad_std
from spectral_cube import SpectralCube
import matplotlib.pyplot as plt

# Input parameters:
#      prename = root of image name (before .image.fits)
#      casalog = name of CASA log file from which to extract run parameters
#      sig_v0  = minimum velocity of signal range in km/s
#      sig_v1  = maximum velocity of signal range in km/s
#      xlims   = x range for plotting, in pixels, e.g. [100,500]
#      ylims   = y range for plotting, in pixels, e.g. [100,500]

def inspect_tclean(prename = 'GMC1_12CO_12m7m', casalog=None, 
    fileext='fits', sig_v0=None, sig_v1=None, xlims=None, ylims=None,
    outfile=None, merge=False):

    image   = prename + '.image.' + fileext                             # image cube
    resid   = prename + '.residual.' + fileext                          # residual cube
    model   = prename + '.convmodel.' + fileext                         # model cube
    mask    = prename + '.mask.' + fileext                              # mask cube

    # Integrate the cube over the chosen velocity range
    if os.path.isfile(image):
        data_im  = SpectralCube.read(image)
        print(data_im)
        if sig_v0 is not None:
            data_sub = data_im.spectral_slab(sig_v0 * u.km / u.s, sig_v1 * u.km / u.s)
        else:
            data_sub = data_im
        data_mom = np.asarray(data_sub.moment(order=0))
        imvals = np.asarray(data_sub).flatten()
        imvals = imvals[~np.isnan(imvals)]
        data_rms = mad_std(imvals)
    else:
        data_mom = None

    # Project the mask over the chosen velocity range
    if os.path.isfile(mask):
        mask_im  = SpectralCube.read(mask)
        if sig_v0 is not None:
            mask_sub = mask_im.spectral_slab(sig_v0 * u.km / u.s, sig_v1 * u.km / u.s)
        else:
            mask_sub = mask_im
        mask_max = np.asarray(mask_sub.max(axis=0))
    else:
        mask_max = None

    # Integrate the model over the chosen velocity range
    if os.path.isfile(model):
        model_im  = SpectralCube.read(model)
        if sig_v0 is not None:
            model_sub = model_im.spectral_slab(sig_v0 * u.km / u.s, sig_v1 * u.km / u.s)
        else:
            model_sub = model_im
        model_mom = np.asarray(model_sub.moment(order=0))
    else:
        model_mom = None

    # Maximum of residual over the chosen velocity range
    if os.path.isfile(resid):
        resid_im  = SpectralCube.read(resid)
        if sig_v0 is not None:
            resid_sub = resid_im.spectral_slab(sig_v0 * u.km / u.s, sig_v1 * u.km / u.s)
        else:
            resid_sub = resid_im
        resvals = np.asarray(resid_sub).flatten()
        resvals = resvals[~np.isnan(resvals)]
        resid_mom = np.asarray(resid_sub.moment(order=0))
    else:
        resid_mom = None

    fig = plt.figure(figsize=(12,8))

    # Plot cleaned image
    if data_mom is not None:
        ax1 = fig.add_subplot(231)
        im1 = plt.imshow(data_mom, cmap='gist_ncar')
        cbar = fig.colorbar(im1, ax=ax1)
        cbar.ax.tick_params(labelsize=8) 
        ax1.contour(mask_max, colors='w', linewidths=1, levels=[0.5])
        ax1.set_title('Moment 0')
        if xlims is not None:
            ax1.set_xlim(xlims)
        if ylims is not None:
            ax1.set_ylim(ylims)
        ax1.tick_params(labelsize=8) 
        #ax1.axes.get_xaxis().set_ticks([])
        #ax1.axes.get_yaxis().set_ticks([])
        # Histogram of values
        ax4 = fig.add_subplot(234)
        ax4.hist(imvals,100,histtype='bar',facecolor='k')
        ax4.set_yscale('log')
        ax4.set_ylabel('PIX number')
        ax4.set_xlabel('Intensity (Jy/beam)')
        ax4.annotate('image cube in vrange', xy=[0.45,0.95], xycoords='axes fraction')

    # Plot residual image
    if resid_mom is not None:
        ax2 = fig.add_subplot(232)
        im2 = plt.imshow(resid_mom, origin='lower', cmap='nipy_spectral')
        cbar = fig.colorbar(im2, ax=ax2)
        cbar.ax.tick_params(labelsize=8) 
        ax2.contour(mask_max, colors='w', linewidths=1, levels=[0.5])
        ax2.set_title('Residual')
        if xlims is not None:
            ax2.set_xlim(xlims)
        if ylims is not None:
            ax2.set_ylim(ylims)
        ax2.axes.get_xaxis().set_ticks([])
        ax2.axes.get_yaxis().set_ticks([])
        # Histogram of values
        ax5 = fig.add_subplot(235)
        ax5.hist(resvals,100,histtype='bar',facecolor='k')
        ax5.set_yscale('log')
        ax5.set_xlabel('Intensity (Jy/beam)')
        ax5.annotate('residual cube in vrange', xy=[0.4,0.95], xycoords='axes fraction')

    # Plot model image
    if model_mom is not None:
        ax3 = fig.add_subplot(233)
        im3 = plt.imshow(model_mom, origin='lower', cmap='gist_ncar')
        cbar = fig.colorbar(im3, ax=ax3)
        cbar.ax.tick_params(labelsize=8) 
        ax3.contour(mask_max, colors='w', linewidths=1, levels=[0.5])
        ax3.set_title('Model')
        if xlims is not None:
            ax3.set_xlim(xlims)
        if ylims is not None:
            ax3.set_ylim(ylims)
        ax3.axes.get_xaxis().set_ticks([])
        ax3.axes.get_yaxis().set_ticks([])


    # Show parameters
    ax6 = fig.add_subplot(236)
    ax6.axis('off')
    xleft = 0.00
    xind  = 0.05
    ytop  = 0.98
    ystep = 0.055
    ax6.text(xleft,ytop,'log = {0}'.format(casalog), transform=ax6.transAxes)
    ax6.text(xleft,ytop-ystep,'image = {0}'.format(image), transform=ax6.transAxes)
    if sig_v0 is not None:
        ax6.text(xleft,ytop-2*ystep,'vrange = {0} to {1} km/s'.format(sig_v0,sig_v1),
            transform=ax6.transAxes)
    else:
        ax6.text(xleft,ytop-2*ystep,'vrange = {0} to {1}'.format(
            data_im.spectral_extrema[0],data_im.spectral_extrema[1]),
            transform=ax6.transAxes)
    ax6.text(xleft,ytop-3*ystep,'bmaj, bmin, bpa = {0:.2f}, {1:.2f}, {2:.2f}'.format(
        3600*data_im.header['bmaj'],3600*data_im.header['bmin'],
        data_im.header['bpa']), transform=ax6.transAxes)
    ax6.text(xleft,ytop-4*ystep,'rms noise = {0:.3f} Jy/bm'.format(
        data_rms), transform=ax6.transAxes)

    # Read parameters from casa.log:
    if os.path.isfile(casalog):
        keep = ''
        with open(casalog) as f:
            lines=f.readlines()
            for line in lines:
                if re.search('INFO\ttclean', line):
                    keep += line
        # Capture values in brackets:
        p1=re.compile("([^\s,\(]+)=(\[.+\]),")
        matches = p1.findall(keep)
        # Capture values in quotes:
        p2=re.compile("([^\s,]+)=\"([^,]+)\",")
        matches += p2.findall(keep)
        # Capture numbers:
        p3=re.compile("([^\s,\(]+)=([\d.]+),")
        matches += p3.findall(keep)
        tclnpar = dict(matches)
        #print(tclnpar)
        # Output parameters
        ax6.text(xleft,ytop-5*ystep,'cell = {0}, imsize = {1}'.format(
            tclnpar['cell'],tclnpar['imsize']), transform=ax6.transAxes)
        ax6.text(xleft,ytop-6*ystep,'weighting = {0}'.format(
            tclnpar['weighting']), transform=ax6.transAxes)
        ax6.text(xleft,ytop-7*ystep,'niter = {0}'.format(
            tclnpar['niter']), transform=ax6.transAxes)
        ax6.text(xleft,ytop-8*ystep,'nsigma = {0}'.format(
            tclnpar['nsigma']), transform=ax6.transAxes)
        ax6.text(xleft,ytop-9*ystep,'pblimit = {0}'.format(
            tclnpar['pblimit']), transform=ax6.transAxes)
        if 'startmodel' in tclnpar:
            ax6.text(xleft,ytop-10*ystep,'startmodel = {0}'.format(
                tclnpar['startmodel']), transform=ax6.transAxes)
            ynext = ytop-11*ystep
        else:
            ynext = ytop-10*ystep
        ax6.text(xleft,ynext,'deconvolver = {0}'.format(
            tclnpar['deconvolver']), transform=ax6.transAxes)
        ynext = ynext-ystep
        if (tclnpar['deconvolver'] == "multiscale"):
            ax6.text(xleft+xind,ynext,'scales = {0}'.format(
                tclnpar['scales']), transform=ax6.transAxes)
            ax6.text(xleft+xind,ynext-ystep,'smallscalebias = {0}'.format(
                tclnpar['smallscalebias']), transform=ax6.transAxes)
            ynext = ynext-2*ystep
        ax6.text(xleft,ynext,'usemask = {0}'.format(
            tclnpar['usemask']), transform=ax6.transAxes)
        if (tclnpar['usemask'] != "user"):
            ax6.text(xleft+xind,ynext-ystep,'pbmask = {0}'.format(
                tclnpar['pbmask']), transform=ax6.transAxes)
        if tclnpar['usemask'].startswith('auto-thresh'):
            ax6.text(xleft+xind,ynext-2*ystep,'maskthreshold = {0}'.format(
                tclnpar['maskthreshold']), transform=ax6.transAxes)
            ax6.text(xleft+xind,ynext-3*ystep,'maskresolution = {0}'.format(
                tclnpar['maskresolution']), transform=ax6.transAxes)
        if tclnpar['usemask'].startswith('auto-multithresh'):
            ax6.text(xleft+xind,ynext-2*ystep,'sidelobethreshold = {0}'.format(
                tclnpar['sidelobethreshold']), transform=ax6.transAxes)
            ax6.text(xleft+xind,ynext-3*ystep,'noisethreshold = {0}'.format(
                tclnpar['noisethreshold']), transform=ax6.transAxes)
            ax6.text(xleft+xind,ynext-4*ystep,'lownoisethreshold = {0}'.format(
                tclnpar['lownoisethreshold']), transform=ax6.transAxes)
            ax6.text(xleft+xind,ynext-5*ystep,'smoothfactor = {0}'.format(
                tclnpar['smoothfactor']), transform=ax6.transAxes)
            ax6.text(xleft+xind,ynext-6*ystep,'minbeamfrac = {0}'.format(
                tclnpar['minbeamfrac']), transform=ax6.transAxes)
            ax6.text(xleft+xind,ynext-7*ystep,'cutthreshold = {0}'.format(
                tclnpar['cutthreshold']), transform=ax6.transAxes)

    plt.subplots_adjust(wspace=0.15, hspace=0.1)
    if outfile is None:
        outfile = prename
    print("output file is {0}".format(outfile + '.tclean_results.pdf'))
    plt.savefig(outfile + '.tclean_results.pdf', bbox_inches='tight')
    plt.close()

    return
