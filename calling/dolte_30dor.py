#!/usr/bin/env python3

from lte import lte
from noisegen import noisegen, rms
from histplot import histplot
import numpy as np
import os

pre = '30Dor'
method = 'peak'
noiseiter = 25
tmin = 8
img12  = '../12CO.combined.20130113.smooth.fits.gz'
img13  = '../13CO.combined.20130113.hdfix.fits.gz'
flat12 = '../12CO.combined.20130113.smooth.flat.fits.gz'
flat13 = '../13CO.combined.20130113.flat.fits.gz'
gain12 = '../20130113.flux.fits.gz'
gain13 = '../20130113.flux.fits.gz'
rms12  = '../mom/30Dor_combined_12CO_dil.rms.fits.gz'
rms13  = '../mom/30Dor_combined_13CO_dil.rms.fits.gz'
mask12 = '../mom/30Dor_combined_12CO_dil.mask.fits.gz'

# Uses lte function/script to create a set of lte images from original images
lte_names   = [img12, img13, rms12, rms13, mask12]

old_dir = os.getcwd() # returns absolute path
lte(files = lte_names, tfloor = tmin, datainfo = pre, tx_method = method, cd='../lte')
#os.chdir('../lte')

# Using noisegen function from noisegen script to create random data
# Runs noisegen function twice: first on 12CO, second on 13CO
out12   = pre + '_12CO21.noiseadd.fits.gz'
out13   = pre + '_13CO21.noiseadd.fits.gz'

noisegen(incube = flat12, gainname = gain12, outname = out12, number = noiseiter)
noisegen(incube = flat13, gainname = gain13, outname = out13, number = noiseiter)

# Again runs lte but now on the random noise images just created, only making the 13CO column density images
for n in range(noiseiter):
    cube12      = pre + '_12CO21.noiseadd.' + str(n + 1) + '.fits.gz'
    cube13      = pre + '_13CO21.noiseadd.' + str(n + 1) + '.fits.gz'
    temperature = (4 * np.random.rand()) + 6 # Returns a random temperature in range 6 - 10 K
    info        = pre + '_noise_' + str(n + 1)
    lte_names   = [cube12, cube13, rms12, rms13, mask12]

    lte(files = lte_names, tfloor = temperature, datainfo = info, tx_method = method, 
        onlywrite = ['outn13cube', 'outn13col'])

# Using rms function from noisegen script creates a rms image based on the newly made 13CO column density random images
rms_names1 = [pre + '_noise_' + str(n + 1) + '_'+method+'_n13cube.fits.gz' for n in range(noiseiter)]
rms_names2 = [pre + '_noise_' + str(n + 1) + '_'+method+'_n13col.fits.gz' for n in range(noiseiter)]
noiseout1  = pre + '_' + method + '_noiserms_n13cube.fits.gz'
noiseout2  = pre + '_' + method + '_noiserms_n13col.fits.gz'
rms(names = rms_names1, outname = noiseout1)
rms(names = rms_names2, outname = noiseout2)

xname1 = pre + '_'+method+'_n13cubeerr.fits.gz'
xname2 = pre + '_'+method+'_n13colerr.fits.gz'
yname1 = pre + '_'+method+'_noiserms_n13cube.fits.gz'
yname2 = pre + '_'+method+'_noiserms_n13col.fits.gz'
oname1 = pre + '_'+method+'_noisecomp_cube'
oname2 = pre + '_'+method+'_noisecomp_col'

histplot(xname = xname1, yname = yname1, snrcut = 0, dolog2d = False, dolog1d = False, 
    nbins = 100, outname = oname1, extrema = [0, 0, 4e15, 4e15])
histplot(xname = xname2, yname = yname2, snrcut = 0, dolog2d = False, dolog1d = False, 
    nbins = 100, outname = oname2, extrema = [0, 0, 6e15, 6e15])

# Clean up scratch files
input("Press enter to delete scratch files")
os.system('rm -f '+pre+'_12CO21.noiseadd.*.fits.gz')
os.system('rm -f '+pre+'_13CO21.noiseadd.*.fits.gz')
for f in rms_names1+rms_names2:
    os.remove(f)

# Return to calling directory
os.chdir(old_dir)
