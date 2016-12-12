#!/usr/bin/env python

from lte import lte
from noisegen import noisegen, rms
import numpy as np

pre = '30Dor'
noiseiter = 25
img12   = '../12CO.combined.20130113.smooth.fits'
img13   = '../13CO.combined.20130113.hdfix.fits'
flat12  = '../12CO.combined.20130113.smooth.flat.fits'
flat13  = '../13CO.combined.20130113.flat.fits'
gain12  = '../20130113.flux.fits'      # gain12 != gain13
gain13  = '../20130113.flux.fits'      # in general
rms12   = '../mom/30Dor_combined_12CO_dil.rms.fits'
rms13   = '../mom/30Dor_combined_13CO_dil.rms.fits'
mask12  = '../mom/30Dor_combined_12CO_dil.mask.fits.gz'

lte_names   = [img12, img13, rms12, rms13, mask12]
lte(files = lte_names, tfloor = 8, datainfo = pre, tx_method = 'cube')

# Using noisegen function creates random data from 12CO and 13CO data sets
out12   = pre + '_12CO21.noiseadd.fits.gz'
out13   = pre + '_13CO21.noiseadd.fits.gz'
noisegen(incube = flat12, gainname = gain12, outname = out12, number = noiseiter)
noisegen(incube = flat13, gainname = gain13, outname = out13, number = noiseiter)

# Using lte function creates only 13CO column density images from noisegen-created randomized data sets
for n in range(noiseiter):
    cube12      = pre + '_12CO21.noiseadd.' + str(n + 1) + '.fits.gz'
    cube13      = pre + '_13CO21.noiseadd.' + str(n + 1) + '.fits.gz'
    temperature = (6 * np.random.rand()) + 6 # Returns a randomized temperature in range 6 - 12 K
    info        = pre + '_noise_' + str(n + 1)
    lte_names   = [cube12, cube13, rms12, rms13, mask12]
    lte(files = lte_names, tfloor = temperature, datainfo = info, tx_method = 'cube', onlywrite = ['outn13cube'])

# Using rms function creates a root-mean-square image based on lte-created 13CO column density images from random data
rms_names = [pre+'_noise_'+str(n+1)+'_cube_n13cube.fits.gz' for n in range(noiseiter)]
noiseout  = pre + '_noise_rms_cube_n13cube.fits.gz'
rms(names = rms_names, outname = noiseout)
