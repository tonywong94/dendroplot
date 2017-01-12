# Written for Python 3.5.2
# Based on various previous iterations of lte.py scripts written by Tony Wong and Evan Wojciechowski
# Approximates an LTE (local thermdynamic equilibrium) mass to determine the amount of CO12 and CO13 is present in a given region of the sky

from astropy.io import fits
from astropy import constants as const
from astropy import units as u
import numpy as np
import sys

def lte(files = [], tfloor = 8., datainfo = '', tx_method = '', onlywrite = []):

    # tx_methods accounted for are 'cube' and 'peak'
    # datainfo should provide info on what source data is from and possibly a number corresponding to some form of iteration

    # Declarations of input and output files

    incube12 = files[0]
    incube13 = files[1]
    inrms12  = files[2]
    inrms13  = files[3]
    inmask12 = files[4]

    outtex12      = datainfo + '_' + tx_method + '_tex12.fits.gz'
    outtau13      = datainfo + '_' + tx_method + '_tau13.fits.gz'
    outtau13err   = datainfo + '_' + tx_method + '_tau13err.fits.gz'
    outn13cube    = datainfo + '_' + tx_method + '_n13cube.fits.gz'
    outn13cubeerr = datainfo + '_' + tx_method + '_n13cubeerr.fits.gz'
    outn13col     = datainfo + '_' + tx_method + '_n13col.fits.gz'
    outn13colerr  = datainfo + '_' + tx_method + '_n13colerr.fits.gz'
    outsnr13      = datainfo + '_' + tx_method + '_n13snr.fits.gz'

    # Load 12CO cube [units K]
    print('Reading {0}...'.format(incube12))
    t12cube, hd3d = fits.getdata(incube12, header = True)
    print('min/max values of 12CO [K] are {0} and {1}\n'.format(np.nanmin(t12cube), np.nanmax(t12cube)))

    # Load 12CO uncertainty [2D plane]
    print('Reading {0}...'.format(inrms12))
    t12err, hd2d = fits.getdata(inrms12, header = True)
    print('min/max values of 12CO uncertainty are {0} and {1}\n'.format(np.nanmin(t12err), np.nanmax(t12err)))

    # Load 12CO mask [3D cube or 2D plane]
    print('Reading {0}...'.format(inmask12))
    mask = fits.getdata(inmask12)
    print('Number of mask == 1 values = {0}'.format(np.count_nonzero(mask[~np.isnan(mask)] > 0)))
    print('Number of mask == 0 values = {0}'.format(np.count_nonzero(mask[~np.isnan(mask)] < 1)))
    print('Number of mask == NaN values = {0}\n'.format(np.count_nonzero(np.isnan(mask))))
    mask3d = (mask == 0)
    mask2d = (np.nansum(mask, axis = 0) == 0)

    # Calculate Tex for Tex > Tfloor
    # Different methods are calculated slightly differently
    with np.errstate(invalid = 'ignore', divide = 'ignore'):
        if tx_method == 'peak':
            t12 = np.nanmax(t12cube, axis == 0)
            hdtx = hd2d
            t12[mask2d] = float('NaN')
        elif tx_method == 'cube':
            t12 = t12cube
            hdtx = hd3d
            t12[mask3d] = float('NaN')

        print('Calculating Tex [excitation temperature]...')
        tex = 11.06 / (np.log(1 + 11.06/(t12 + 0.187)))
        tex[tex < tfloor] = tfloor
        print('min/max values of Tex [K] are {0} and {1}\n'.format(np.nanmin(tex), np.nanmax(tex)))

    if (len(onlywrite) == 0) or ('outtex12' in onlywrite) == True:
        hdtx['datamin'] = np.nanmin(tex)
        hdtx['datamax'] = np.nanmax(tex)
        hdtx['tfloor'] = tfloor
        fits.writeto(outtex12, tex, hdtx, clobber = True)
        print('File {0} successfully written'.format(outtex12))

    # Load 13CO cube [units K]
    print('Reading {0}...'.format(incube13))
    t13, hd3d = fits.getdata(incube13, header = True)
    print('min/max values of 13CO [K] are {0} and {1}\n'.format(np.nanmin(t13), np.nanmax(t13)))

    # Load 13CO uncertainty [2D plane]
    print('Reading {0}...'.format(inrms13))
    t13err, hd2d = fits.getdata(inrms13, header = True)
    print('min/max values of 13CO undertainty are {0} and {1}\n'.format(np.nanmin(t13err), np.nanmax(t13err)))

    # Calculates 13CO optical depth cube
    with np.errstate(invalid = 'ignore'):
        print('Calculating tau13 [13CO optical depth]...')
        tau13 = -np.log(1-(t13/10.6)/(1/(np.exp(10.6/tex)-1)-1/(np.exp(10.6/2.73)-1)))
        print('min/max values of 13CO optical depth are {0} are {1}\n'.format(np.nanmin(tau13), np.nanmax(tau13)))
    
    if (len(onlywrite) == 0) or ('outtau13' in onlywrite) == True:
        hd3d['datamin'] = np.nanmin(tau13)
        hd3d['datamax'] = np.nanmax(tau13)
        hd3d['bunit'] = ''
        hd3d['tfloor'] = tfloor
        fits.writeto(outtau13, tau13, hd3d, clobber = True)
        print('File {0} successfully written'.format(outtau13))

    # Uncertainty of 13CO optical depth [linear approximation]
    print('Calculating error in tau13...')
    tau13err = (t13err/10.6)/(1/(np.exp(10.6/tex)-1)-1/(np.exp(10.6/2.73)-1))
    print('min/max values of 13CO optical depth uncertainty are {0} and {1}\n'.format(np.nanmin(tau13err), np.nanmax(tau13err)))

    if (len(onlywrite) == 0) or ('outtau13err' in onlywrite) == True:
        hdtx['datamin'] = np.nanmin(tau13err)
        hdtx['datamin'] = np.nanmax(tau13err)
        hdtx['bunit'] = ''
        hdtx['tfloor'] = tfloor
        fits.writeto(outtau13err, tau13err, hdtx, clobber = True)
        print('File {0} successfully written'.format(outtau13err))

    # Calculates 13CO column density and error cubes
    print('Calculating the 13CO column density and error...\n')
    # Equation from Bourke et al. (1997ApJ...476...781B) equation (A4)
    B = 55.101e9/u.s               # Rotational constant of 13CO
    jbot = 1                    # J of bottom state
    mu2 = (0.112 * 1e-18)**2 * u.cm**3 * u.erg  # Use 0.113 debye for dipole moment of 13CO
    hB_3k = const.h*B/(3*const.k_B)
    cm2perKkms = u.cm**-2*u.s/(u.km*u.K)
    prefac = (3*(const.h)/(8*np.pi**2*mu2)*const.k_B/((jbot+1)*const.h*B)).to(cm2perKkms)
    print('Pre-factor is {0}'.format(prefac))

    with np.errstate(invalid = 'ignore', divide = 'ignore'):
        n13 = prefac * (tex*u.K+hB_3k) * np.exp(5.29/tex) * tau13 / (1 - np.exp(-10.6/tex))
        if tx_method == 'peak':
            tau13ecube = np.tile(tau13err, (np.shape(n13)[0], 1, 1))
            n13ecube = (n13/tau13) * tau13ecube
        else:
            n13ecube = (n13/tau13) * tau13err
    n13[mask3d] = float('NaN')
    n13ecube[mask3d] = float('NaN')
    print('min/max values of 13CO column density are {0} and {1}\n'.format(np.nanmin(n13), np.nanmax(n13)))
    print('min/max values of 13CO column density uncertainty are {0} and {1}\n'.format(np.nanmin(n13ecube), np.nanmax(n13ecube)))

    # Write column density and error cubes
    if (len(onlywrite) == 0) or ('outn13cube' in onlywrite) == True:
        hd3d['datamin'] = np.nanmin(n13).value
        hd3d['datamax'] = np.nanmax(n13).value
        hd3d['bunit'] = 'cm^-2 / (km/s)'
        hd3d['tfloor'] = tfloor
        fits.writeto(outn13cube, n13, hd3d, clobber = True)
        print('File {0} successfully written'.format(outn13cube))

    if (len(onlywrite) == 0) or ('outn13cubeerr' in onlywrite) == True:
        hd3d['datamin'] = np.nanmin(n13ecube).value
        hd3d['datamax'] = np.nanmax(n13ecube).value
        hd3d['tfloor'] = tfloor
        fits.writeto(outn13cubeerr, n13ecube, hd3d, clobber = True)
        print('File {0} successfully written'.format(outn13cubeerr))

    # Make integrated column density and error maps
    print('Calculating integrated column density and error...')
    n13col = np.nansum(n13, axis = 0) * (hd3d['cdelt3']/1000.) * u.km/u.s
    masksum = np.nansum(mask, axis = 0)
    with np.errstate(all = 'ignore'):
#       n13colerr = np.nansum(n13ecube, axis = 0) * (hd3d['cdelt3']/1000.) * u.km/u.s / np.sqrt(masksum)
        n13colerr = np.sqrt(np.nansum(n13ecube**2, axis = 0)) * (hd3d['cdelt3']/1000.) * u.km/u.s # ?? / np.sqrt(masksum)

    # Apply 2D mask
    print('Applying mask to column density and error...')
    n13col[mask2d] = float('NaN')
    n13colerr[mask2d] = float('NaN')
    print('min/max values of 13CO column density are {0} and {1}\n'.format(np.nanmin(n13colerr), np.nanmax(n13colerr)))
    print('min/max values of 13CO column density uncertainty are {0} and {1}\n'.format(np.nanmin(n13colerr), np.nanmax(n13colerr)))

    # Write the FITS files
    if (len(onlywrite) == 0) or ('outn13col' in onlywrite) == True:
        hd2d['datamin'] = np.nanmin(n13col).value
        hd2d['datamax'] = np.nanmax(n13col).value
        hd2d['bunit'] = 'cm^-2'
        hd2d['tfloor'] = tfloor
        fits.writeto(outn13col, n13col, hd2d, clobber = True)
        print('File {0} successfully written'.format(outn13col))

    if (len(onlywrite) == 0) or ('outn13colerr' in onlywrite) == True:
        hd2d['datamin'] = np.nanmin(n13colerr).value
        hd2d['datamax'] = np.nanmax(n13colerr).value
        hd2d['tfloor'] = tfloor
        fits.writeto(outn13colerr, n13colerr, hd2d, clobber = True)
    
    # Write SNR image
    print('Calculating signal-to-noise ratio...')
    with np.errstate(invalid = 'ignore'):
        n13snr = n13col / n13colerr
    print('min/max values of SNR are {0} and {1}\n'.format(np.nanmin(n13snr), np.nanmax(n13snr)))
    if (len(onlywrite) == 0) or ('outsnr13' in onlywrite) == True:
        hd2d['datamin'] = np.nanmin(n13snr).value
        hd2d['datamax'] = np.nanmax(n13snr).value
        hd2d['bunit'] = ''
        hd2d['tfloor'] = tfloor
        fits.writeto(outsnr13, n13snr, hd2d, clobber = True)
        print('File {0} successfully written'.format(outsnr13))


### END OF FILE ###
