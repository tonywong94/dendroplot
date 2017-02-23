#!/usr/bin/env python3

from lte import lte
from noisegen import noisegen, rms
import numpy as np
import os

pre = 'PCC'
noiseiter = 25
img12   = '../PCC_12mTP_12CO21.pbcor.K.fits.gz'
img13   = '../PCC_12mTP_13CO21.pbcor.K.fits.gz'
flat12  = '../PCC_12mTP_12CO21.image.K.fits.gz'
flat13  = '../PCC_12mTP_13CO21.image.K.fits.gz'
gain12  = '../PCC_12mTP_12CO21.flux1.fits.gz'
gain13  = '../PCC_12mTP_13CO21.flux1.fits.gz'
rms12   = '../mom/PCC_12mTP_12CO_dil.rms.fits.gz'
rms13   = '../mom/PCC_12mTP_13CO_dil.rms.fits.gz'
mask12  = '../mom/PCC_12mTP_12CO_dil.mask.fits.gz'

# Uses lte function/script to create a set of lte images from original images
lte_names   = [img12, img13, rms12, rms13, mask12]

lte(files = lte_names, tfloor = 8, datainfo = pre, tx_method = 'cube', cd='../lte')

# Using noisegen function from noisegen script to create random data from 12CO and 13CO data sets
# Runs noisegen function twice: first on 12CO, second on 13CO
out12   = pre + '_12CO21.noiseadd.fits.gz'
out13   = pre + '_13CO21.noiseadd.fits.gz'

noisegen(incube = flat12, gainname = gain12, outname = out12, number = noiseiter)
noisegen(incube = flat13, gainname = gain13, outname = out13, number = noiseiter)

# Again runs lte but now on the random noise images just created, only making the 13CO column density images
for n in range(noiseiter):
    cube12      = pre + '_12CO21.noiseadd.' + str(n + 1) + '.fits.gz'
    cube13      = pre + '_13CO21.noiseadd.' + str(n + 1) + '.fits.gz'
    temperature = (6 * np.random.rand()) + 6 # Returns a random temperature in range 6 - 12 K
    info        = pre + '_noise_' + str(n + 1)
    lte_names   = [cube12, cube13, rms12, rms13, mask12]

    lte(files = lte_names, tfloor = temperature, datainfo = info, tx_method = 'cube', onlywrite = ['outn13cube'])

# Using rms function from noisegen script creates a rms (root-mean-square) image based on the newly made 13CO column density random images
rms_names = [pre + '_noise_' + str(n + 1) + '_cube_n13cube.fits.gz' for n in range(noiseiter)]
noiseout  = pre + '_noise_rms_cube_n13cube.fits.gz'

rms(names = rms_names, outname = noiseout)

# Clean up scratch files
input("Press enter to delete scratch files")
os.system('rm -f '+pre+'_12CO21.noiseadd.*.fits.gz')
os.system('rm -f '+pre+'_13CO21.noiseadd.*.fits.gz')
for f in rms_names:
    os.remove(f)
