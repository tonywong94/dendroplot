#!/usr/bin/env python3

# Written by Tony Wong and Evan Wojciechowski
# Combined noise generation and rms calculation script
# Contains functinos noisegen and rms

from astropy.convolution import convolve_fft
from astropy.io import fits
from astropy.stats import mad_std
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import os
from radio_beam.beam import Beam

def noisegen(incube, gainname = '', outname='30Dor_13CO21.noiseadd.fits.gz', number = 1, 
    verbose = False, tomake = False, cd = ''):
    """
    From an input cube (with noise included) generate noise-added cubes

    PARAMETERS:
    incube   - Input 2D or 3D image, with constant noise across the field
    gainname - Input image with same shape as incube providing point source gain
    outname  - Output image(s) with noise added and gain correction applied
    number   - Number of output images to generate
    verbose  - True to provide additional output to terminal
    tomake   - True to generate additional FITS images for debugging
    
    """
    
    if cd != '':
        if os.path.isdir(cd) == 1:
            print('Found {}, changing directory...'.format(cd))
            os.chdir(cd)
        else:
            print('Directory {} doesn\'t exist, creating and changing...\n'.format(cd))
            os.mkdir(cd)
            os.chdir(cd)
    
    # file paths need to be absolute or defined properly in relation to working directory
    if os.path.exists(incube) == 1:
        print('Found {}...'.format(incube))
    else:
        print('File {} does not exist'.format(incube))
        return

    # Reads input data
    print('Reading {0}...'.format(incube))
    indata, inheader = fits.getdata(incube, 0, header = True)
    data = indata.flatten()
    #data = indata[0:16, :, :].flatten()
    rmask = mad_std(data[~np.isnan(data)])
    print('\nStandard deviation of input data is {0}'.format(rmask))
    if verbose == True:
        print('\nShape of data: {0}'.format(indata.shape))

    # Reads gain values
    if gainname != '':
        gain = fits.getdata(gainname)
        if verbose == True:
            print('\nShape of gain: {0}\n'.format(gain.shape))

    # Creates beam class
    bm = Beam.from_fits_header(inheader)
    if verbose == True:
        print('\nBeam: {0}'.format(bm))

    # Returns an elliptical Gaussian kernel of the beam
    gauss2d = bm.as_kernel(inheader['CDELT2']*u.deg)

    # Begin creating randomized data sets, writes fits images
    for n in range(number):
        # Generates random data
        random_data = np.random.normal(size = indata.shape)
        if verbose == True:
            print('\nStandard deviation before convolution is {0}'.format(
                np.nanstd(random_data)))

        if tomake == True:
            inheader['DATAMIN'] = random_data.min()
            inheader['DATAMAX'] = random_data.max()
            randomout = outname.replace('.fits', '.random_'+str(n + 1)+'.fits')
            fits.writeto(randomout, random_data, inheader, clobber = True)
            print('File {0} successfully written'.format(randomout))

        # Convolves data per plane in random data
        for plane in np.arange(indata.shape[0]):
            random_data[plane, :, :] = convolve_fft(random_data[plane, :, :], 
                gauss2d, normalize_kernel = True)

        # Calculates rms of smoothed random data
        noise_rms = np.nanstd(random_data)
        if verbose == True:
            print('\nStandard deviation after convolution is {0}'.format(
                noise_rms))

        if tomake == True:
            smoothout = outname.replace('.fits', '.ransmo_'+str(n + 1)+'.fits')
            fits.writeto(smoothout, random_data, inheader, clobber = True)
            print('File {0} successfully written'.format(smoothout))

        # Adds random noise to input data and writes a new fits image
        outdata = indata + (random_data * rmask) / noise_rms
        noiseout = mad_std(outdata[~np.isnan(outdata)])
        if gainname!= '':
            outdata = outdata / gain
        inheader['DATAMIN'] = np.nanmin(outdata)
        inheader['DATAMAX'] = np.nanmax(outdata)
        addout = outname.replace('.fits', '.'+str(n + 1)+'.fits')
        fits.writeto(addout, outdata, inheader, clobber = True)
        print('\nStandard deviation for {0} is {1}'.format(addout,noiseout))
    return

def rms(names, outname, cd = ''):
    """
    Generate a per-pixel RMS image from a set of input images.
    Includes proper handling of image masks.

    PARAMETERS:
    names   - List of input FITS files
    outname - Output image which is the RMS of the input images
    
    """
    # Computes rms from a set of fits images containing randomized data 

    if cd != '':
        if os.path.isdir(cd) == 1:
            print('Found {}, changing directory...'.format(cd))
            os.chdir(cd)
        else:
            print('Directory {} doesn\'t exist, creating and changing...\n'.format(cd))
            os.mkdir(cd)
            os.chdir(cd)
    
    # file paths need to be absolute or defined properly in relation to working directory
    for n in names:
        if os.path.exists(n) == 1:
            print('Found {}...'.format(n))
            continue
        else:
            print('File {} does not exist'.format(n))
            return   

    ### Extract data from first file ###
    print('Reading {0}...'.format(names[0]))
    data, header = fits.getdata(names[0], header = True)
    '''
    ### Need to create a hit-detection array; will record a 0 if a pixel has 
    a NaN value, will record a 1 if a pixel has an numerical value ###

    numpy.isnan(array)
        # Only returns booleans -> need to convert to integer 1 and 0
    numpy.isfinite(array)
        # Also only returns booleans, works without doing a flip of isnan
    array.astype(int)
        # Converts boolean true => 1, false => 0
    '''
    hitdet = (np.isfinite(data)).astype(int)

    data = np.nan_to_num(data)  # Eliminates NaN values for summation
    datasqr = data**2           # Running square of data
    ### End of first file data extraction ###

    ### Begin automated iteration ###
    for i in range(len(names) - 1):
        print('Reading {0}...'.format(names[i + 1]))
        newdata = fits.getdata(names[i + 1])            # Reads new cube
        hitdet += (np.isfinite(newdata)).astype(int)    # Adds new hits to hitdet array
        data += np.nan_to_num(newdata)                  # Adds new data to data array
        datasqr += (np.nan_to_num(newdata)**2)          # Adds squared new data to squared data array
    ### End of interating ###

    ## Begin variance calculation
    print('Calculating variance...')
    with np.errstate(divide = 'ignore', invalid = 'ignore'):
        variance = np.divide(datasqr, hitdet) - (np.divide(data, hitdet)**2)
    variance[hitdet == 0] = np.nan  # Should remask variance values
    rms = np.sqrt(variance)

    header['datamin'] = np.nanmin(rms)
    header['datamax'] = np.nanmax(rms)

    fits.writeto(outname, rms, header, clobber = True)
    print('File {0} written successfully.'.format(outname))
    return
