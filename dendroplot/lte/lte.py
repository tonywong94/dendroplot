#!/usr/bin/env python

# Original written in Python 3.5.2
# Based on various previous iterations of lte.py scripts written by Tony Wong and Evan Wojciechowski
# Approximates an LTE (local thermdynamic equilibrium) mass to determine the amount of CO12 and CO13 is present in a given region of the sky

from astropy.io import fits
from astropy import constants as const
from astropy import units as u
import numpy as np
import os
from os.path import join
import warnings
import sys

def lte(files = [], tfloor = 8., tx_method = 'peak', onlywrite = [], 
        indir = '', outdir = '', outname=None):
    """
    Calculate LTE column density from J=2-1 or J=1-0 transitions of 12CO and 13CO.

    Parameters
    ----------
    # files are in this order: [incube12, incube13, inrms12, inrms13, inmask12]
    
    tfloor : float, optional
        Floor (minimum value) to impose on excitation temperature, in K.
    tx_method : string, optional
        Method to use for 12CO excitation temperature.
        'cube': Use the 12CO temperature at each voxel (pixel and channel)
        'peak': Use the peak 12CO temperature at each pixel
        Default is 'peak'.
    indir : string, optional
        Directory where input files reside.
        Default: Read from the current directory.
    outdir : string, optional
        Directory to write the output files.
        Default: Write to the current directory.
    outname : string, optional
        Basename for output files.  For instance, outname='foo' produces files
        'foo_peak_tex12.fits.gz', etc.
        Default: Based on root name of files[0].
    """

    # Declarations of input and output files
    incube12 = join(indir, files[0])
    incube13 = join(indir, files[1])
    inrms12  = join(indir, files[2])
    inrms13  = join(indir, files[3])
    inmask12 = join(indir, files[4])

    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    if outname is not None:
        basename = outname
    else:
        basename = os.path.basename(files[0]).split('.')[0]
    outtex12      = join(outdir, basename + '_' + tx_method + '_tex12.fits.gz')
    outtau13      = join(outdir, basename + '_' + tx_method + '_tau13.fits.gz')
    outtau13err   = join(outdir, basename + '_' + tx_method + '_tau13err.fits.gz')
    outtau13pk    = join(outdir, basename + '_' + tx_method + '_tau13pk.fits.gz')
    outn13cube    = join(outdir, basename + '_' + tx_method + '_n13cube.fits.gz')
    outn13cubeerr = join(outdir, basename + '_' + tx_method + '_n13cubeerr.fits.gz')
    outn13col     = join(outdir, basename + '_' + tx_method + '_n13col.fits.gz')
    outn13colerr  = join(outdir, basename + '_' + tx_method + '_n13colerr.fits.gz')
    outsnr13      = join(outdir, basename + '_' + tx_method + '_n13snr.fits.gz')

    # Load 12CO cube [units K]
    print('\nReading {0}...'.format(incube12))
    t12cube, hd3d = fits.getdata(incube12, header = True)
    if 'RESTFREQ' in hd3d.keys():
        freq12 = hd3d['RESTFREQ'] * u.Hz
    elif 'RESTFRQ' in hd3d.keys():
        freq12 = hd3d['RESTFRQ'] * u.Hz
    print('The 12CO rest frequency is {0:.4f}'.format((freq12).to(u.GHz)))
    print('min/max values of 12CO [K] are {0:.2f} and {1:.2f}'.format(
        np.nanmin(t12cube), np.nanmax(t12cube)))

    # Load 12CO uncertainty [2D plane]
    print('\nReading {0}...'.format(inrms12))
#     t12err, hd2d = fits.getdata(inrms12, header = True)
#     print('min/max values of 12CO uncertainty are {0:.3f} and {1:.3f}'.format(
#         np.nanmin(t12err), np.nanmax(t12err)))
    hd2d = fits.getheader(inrms12)
    if 'datamin' in hd2d and 'datamax' in hd2d:
        print('min/max values of 12CO uncertainty are {0:.3f} and {1:.3f}'.format(
            hd2d['datamin'], hd2d['datamax']))
        if hd2d['naxis'] == 3:
            for k in list(hd2d['*3*'].keys()):
                hd2d.remove(k)
            for frq in list(hd2d['*frq*'].keys()):
                hd2d.remove(frq)

    # Load 12CO mask [3D cube or 2D plane]
    print('\nReading {0}...'.format(inmask12))
    mask = fits.getdata(inmask12)
    print('Number of mask == 1 values: {0}'.format(np.count_nonzero(
        mask[~np.isnan(mask)] > 0)))
    print('Number of mask == 0 values: {0}'.format(np.count_nonzero(
        mask[~np.isnan(mask)] < 1)))
    print('Number of mask == NaN values: {0}'.format(np.count_nonzero(
        np.isnan(mask))))
    mask3d = (mask == 0)
    mask2d = (np.nansum(mask, axis = 0) == 0)

    # Calculate Tex for Tex > Tfloor
    # Different methods are calculated slightly differently
    with np.errstate(invalid = 'ignore', divide = 'ignore'):
        if tx_method == 'peak':
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                t12 = np.nanmax(t12cube, axis = 0)
            hdtx = hd2d
            t12[mask2d] = np.nan
        elif tx_method == 'cube':
            t12 = t12cube
            hdtx = hd3d
            t12[mask3d] = np.nan

        print('\nCalculating Tex [excitation temperature]...')
        tcmb = 2.73 * u.K
        t0_12 = (const.h * freq12 / const.k_B).to(u.K)
        Jtcmb = t0_12/(np.exp(t0_12/tcmb)-1)
        #tex = 11.06 / (np.log(1 + 11.06/(t12 + 0.187)))
        tex = t0_12 / (np.log(1 + t0_12/((t12*u.K) + Jtcmb)))
        tex[tex < (tfloor * u.K)] = (tfloor * u.K)
        print('min/max values of Tex [K] are {0:.2f} and {1:.2f}'.format(
            np.nanmin(tex), np.nanmax(tex)))

    if (len(onlywrite) == 0) or ('outtex12' in onlywrite) == True:
        hdtx['datamin'] = np.nanmin(tex).value
        hdtx['datamax'] = np.nanmax(tex).value
        hdtx['tfloor'] = tfloor
        fits.writeto(outtex12, tex, hdtx, overwrite = True)
        print('File {0} successfully written'.format(outtex12))

    # Load 13CO cube [units K]
    print('\nReading {0}...'.format(incube13))
    t13, hd3d = fits.getdata(incube13, header = True)
    if 'RESTFREQ' in hd3d.keys():
        freq13 = hd3d['RESTFREQ'] * u.Hz
    elif 'RESTFRQ' in hd3d.keys():
        freq13 = hd3d['RESTFRQ'] * u.Hz
    print('The 13CO rest frequency is {0:8.4f}'.format((freq13).to(u.GHz)))
    print('min/max values of 13CO [K] are {0:.2f} and {1:.2f}'.format(
        np.nanmin(t13), np.nanmax(t13)))

    # Load 13CO uncertainty [2D plane]
    print('\nReading {0}...'.format(inrms13))
    t13err, hd2d = fits.getdata(inrms13, header = True)
    print('min/max values of 13CO uncertainty are {0:.3f} and {1:.3f}'.format(
        np.nanmin(t13err), np.nanmax(t13err)))

    # Calculate 13CO optical depth cube
    with np.errstate(invalid = 'ignore'):
        print('\nCalculating tau13 [13CO optical depth]...')
        t0_13 = (const.h * freq13 / const.k_B).to(u.K)
        #tau13 = -np.log(1-(t13/10.6)/(1/(np.exp(10.6/tex)-1)-1/(np.exp(10.6/2.73)-1)))
        tau13 = -np.log(1-((t13*u.K)/t0_13)/(1/(np.exp(t0_13/tex)-1)-1/
            (np.exp(t0_13/tcmb)-1)))
        print('min/max values of 13CO optical depth are {0:.2f} are {1:.2f}'.format(
            np.nanmin(tau13), np.nanmax(tau13)))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            tau13peak = np.nanmax(tau13, axis = 0)
        print('min/max values of peak 13CO optical depth are {0:.2f} are {1:.2f}'.format(
            np.nanmin(tau13peak), np.nanmax(tau13peak)))
    
    if (len(onlywrite) == 0) or ('outtau13' in onlywrite) == True:
        hd3d['datamin'] = np.nanmin(tau13).value
        hd3d['datamax'] = np.nanmax(tau13).value
        hd3d['bunit'] = ''
        hd3d['tfloor'] = tfloor
        fits.writeto(outtau13, tau13, hd3d, overwrite = True)
        print('File {0} successfully written'.format(outtau13))

    if (len(onlywrite) == 0) or ('outtau13pk' in onlywrite):
        hd2d['datamin'] = np.nanmin(tau13peak).value
        hd2d['datamax'] = np.nanmax(tau13peak).value
        hd2d['bunit'] = ''
        hd2d['tfloor'] = tfloor
        fits.writeto(outtau13pk, tau13peak, hd2d, overwrite = True)
        print('File {0} successfully written'.format(outtau13pk))

    # Uncertainty of 13CO optical depth [linear approximation]
    print('\nCalculating error in tau13...')
    #tau13err = (t13err/10.6)/(1/(np.exp(10.6/tex)-1)-1/(np.exp(10.6/2.73)-1))
    tau13err = (t13err*u.K/t0_13)/(1/(np.exp(t0_13/tex)-1)-1/(np.exp(t0_13/tcmb)-1))
    print('min/max values of 13CO tau uncertainty are {0:.3f} and {1:.3f}'.format(
        np.nanmin(tau13err), np.nanmax(tau13err)))

    if (len(onlywrite) == 0) or ('outtau13err' in onlywrite) == True:
        hdtx['datamin'] = np.nanmin(tau13err).value
        hdtx['datamax'] = np.nanmax(tau13err).value
        hdtx['bunit'] = ' '
        hdtx['tfloor'] = tfloor
        fits.writeto(outtau13err, tau13err, hdtx, overwrite = True)
        print('File {0} successfully written'.format(outtau13err))

    # Calculate 13CO column density and error cubes
    print('\nCalculating the 13CO column density cube and error...')
    # Equation from Bourke et al. (1997ApJ...476...781B) equation (A4)
    B = 55.101e9 * u.Hz             # Rotational constant of 13CO
    jbot = round((freq13/(2*B)).value - 1)      # J of bottom state
    print('J value for bottom state is {0}'.format(jbot))
    mu2 = (0.112 * 1e-18)**2 * u.cm**3 * u.erg  # 0.112 debye for dipole moment of 13CO
    hB_3k = (const.h * B/(3 * const.k_B)).to(u.K)
    cm2perKkms = u.cm**-2*u.s/(u.km*u.K)
    # Changed np.pi**2 to np.pi**3 on 3 May 2019
    prefac = (3*(const.h)/(8*np.pi**3*mu2)*const.k_B/((jbot+1)*const.h*B)).to(cm2perKkms)
    print('Pre-factor is {0}'.format(prefac))

    with np.errstate(invalid = 'ignore', divide = 'ignore'):
        ntotexp = np.exp(const.h * B * jbot * (jbot+1) / (const.k_B * tex))
        #n13 = prefac * (tex*u.K+hB_3k) * np.exp(5.29/tex) * tau13 / (1 - np.exp(-10.6/tex))
        n13 = prefac * (tex + hB_3k) * ntotexp * tau13 / (1 - np.exp(-t0_13/tex))
        if tx_method == 'peak' and t13err.ndim == 2:
            tau13ecube = np.tile(tau13err, (np.shape(n13)[0], 1, 1))
            n13ecube = (n13/tau13) * tau13ecube
        else:
            n13ecube = (n13/tau13) * tau13err
    n13[mask3d] = float('NaN')
    n13ecube[mask3d] = float('NaN')
    print('min/max values of N(13CO)/∆v are {0:.4E} and {1:.4E}'.format(
        np.nanmin(n13), np.nanmax(n13)))
    print('min/max values of N(13CO)/∆v uncertainty are {0:.4E} and {1:.4E}'.format(
        np.nanmin(n13ecube), np.nanmax(n13ecube)))

    # Write column density and error cubes
    if (len(onlywrite) == 0) or ('outn13cube' in onlywrite):
        hd3d['datamin'] = np.nanmin(n13).value
        hd3d['datamax'] = np.nanmax(n13).value
        hd3d['bunit'] = 'cm^-2 / (km/s)'
        hd3d['tfloor'] = tfloor
        fits.writeto(outn13cube, n13, hd3d, overwrite = True)
        print('File {0} successfully written'.format(outn13cube))

    if (len(onlywrite) == 0) or ('outn13cubeerr' in onlywrite):
        hd3d['datamin'] = np.nanmin(n13ecube).value
        hd3d['datamax'] = np.nanmax(n13ecube).value
        hd3d['tfloor'] = tfloor
        fits.writeto(outn13cubeerr, n13ecube, hd3d, overwrite = True)
        print('File {0} successfully written'.format(outn13cubeerr))

    # Calculate integrated column density and error maps
    print('\nCalculating integrated column density and error...')
    n13col = np.nansum(n13, axis = 0) * abs(hd3d['cdelt3']/1000.) * u.km/u.s
    n13col[np.all(np.isnan(n13), axis=0)] = np.nan
    with np.errstate(all = 'ignore'):
        n13colerr = np.sqrt(np.nansum(n13ecube**2, axis = 0)) * abs(
            hd3d['cdelt3']/1000.) * u.km/u.s
        n13colerr[np.all(np.isnan(n13ecube), axis=0)] = np.nan
    print('min/max values of N(13CO) are {0:.4E} and {1:.4E}'.format(
        np.nanmin(n13col), np.nanmax(n13col)))
    print('min/max values of N(13CO) uncertainty are {0:.4E} and {1:.4E}'.format(
        np.nanmin(n13colerr), np.nanmax(n13colerr)))

    # Write integrated column density and error maps
    if (len(onlywrite) == 0) or ('outn13col' in onlywrite):
        hd2d['datamin'] = np.nanmin(n13col).value
        hd2d['datamax'] = np.nanmax(n13col).value
        hd2d['bunit'] = 'cm^-2'
        hd2d['tfloor'] = tfloor
        fits.writeto(outn13col, n13col, hd2d, overwrite = True)
        print('File {0} successfully written'.format(outn13col))

    if (len(onlywrite) == 0) or ('outn13colerr' in onlywrite):
        hd2d['datamin'] = np.nanmin(n13colerr).value
        hd2d['datamax'] = np.nanmax(n13colerr).value
        hd2d['tfloor'] = tfloor
        fits.writeto(outn13colerr, n13colerr, hd2d, overwrite = True)
        print('File {0} successfully written'.format(outn13colerr))
    
    # Write SNR image
    print('\nCalculating signal-to-noise ratio...')
    with np.errstate(invalid = 'ignore'):
        n13snr = n13col / n13colerr
    print('min/max values of SNR are {0:.2f} and {1:.2f}'.format(
        np.nanmin(n13snr), np.nanmax(n13snr)))

    if (len(onlywrite) == 0) or ('outsnr13' in onlywrite):
        hd2d['datamin'] = np.nanmin(n13snr).value
        hd2d['datamax'] = np.nanmax(n13snr).value
        hd2d['bunit'] = ' '
        hd2d['tfloor'] = tfloor
        fits.writeto(outsnr13, n13snr, hd2d, overwrite = True)
        print('File {0} successfully written'.format(outsnr13))

### END OF FILE ###
