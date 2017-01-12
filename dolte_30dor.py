#!/usr/bin/env python

from lte import lte
from noisegen import noisegen, rms
import numpy as np

pre = '30Dor'
noiseiter = 25
img12  = '../12CO.combined.20130113.smooth.fits'
img13  = '../13CO.combined.20130113.hdfix.fits'
flat12 = '../12CO.combined.20130113.smooth.flat.fits'
flat13 = '../13CO.combined.20130113.flat.fits'
gain12 = '../20130113.flux.fits'
gain13 = '../20130113.flux.fits'
rms12  = '../mom/30Dor_combined_12CO_dil.rms.fits'
rms13  = '../mom/30Dor_combined_13CO_dil.rms.fits'
mask12 = '../mom/30Dor_combined_12CO_dil.mask.fits.gz'

# Uses lte function/script to create a set of lte images from original images
lte_names = [img12, img13, rms12, rms13, mask12]

lte(files = lte_names, tfloor = 8, datainfo = pre, tx_method = 'cube')

# Using noisegen function from noisegen script to create random data from 12CO and 13CO data sets
# Runs noisegen function twice: first on 12CO, second on 13CO
out12 = pre + '_12CO21.noiseadd.fits.gz'
out13 = pre + '_13CO21.noiseadd.fits.gz'

noisegen(incube = flat12, gainname = gain12, outname = out12, number = noiseiter)
noisegen(incube = flat13, gainname = gain13, outname = out13, number = noiseiter)

# Again runs lte but now on the random noise images just created, only making the 13CO column density images
for n in range(noiseiter):
    cube12      = pre + '_12CO21.noiseadd.' + str(n + 1) + '.fits.gz'
    cube13      = pre + '_13CO21.noiseadd.' + str(n + 1) + '.fits.gz'
    temperature = (6 & np.random.rand()) + 6 # Returns a random temperature in range 6 - 12 K
    info        = pre + '_noise_' + str(n + 1)
    lte_names   = [cube12, cube13, rms12, rms13, mask12]
    
    lte(files = lte_names, tfloor = temperature, datainfo = info, tx_method = 'cube', onlywrite = ['outn13cube'])

# Using rms function from noisegen script creates a rms (root-mean-square) image based on the newly made 13CO column density random images
rms_names = [pre + '_noise_' + str(n + 1) + '_cube_n13cube.fits.gz' for n in range(noiseiter)]
noiseout  = pre + '_noise_rms_cube_n13cube.fits.gz'

rms(names = rms_names, outname = noiseout)
