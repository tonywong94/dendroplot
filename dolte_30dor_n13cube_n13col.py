from histplot import histplot
from lte import lte
from noisegen import noisegen, rms
import numpy as np

pre = '30dor'
noiseiter = 25

img12  = '../12CO.combined.20130113.smooth.fits.gz'
img13  = '../13CO.combined.20130113.hdfix.fits.gz'
flat12 = '../12CO.combined.20130113.smooth.flat.fits.gz'
flat13 = '../13CO.combined.20130113.flat.fits.gz'
gain12 = '../20130113.flux.fits.gz'
gain13 = '../20130113.flux.fits.gz'
rms12  = '../mom/30Dor_combined_12CO_dil.rms.fits.gz'
rms13  = '../mom/30Dor_combined_13CO_dil.rms.fits.gz'
mask12 = '../mom/30Dor_combined_12CO_dil.mask.fits.gz'

lte_names = [img12, img13, rms12, rms13, mask12]
lte(files = lte_names, tfloor = 8., datainfo = pre, tx_method = 'cube')

out12 = pre + '_12CO21.noiseadd.fits.gz'
out13 = pre + '_13CO21.noiseadd.fits.gz'

noisegen(incube = flat12, gainname = gain12, outname = out12, number = noiseiter)
noisegen(incube = flat13, gainname = gain13, outname = out13, number = noiseiter)

for n in range(noiseiter):
    cube12      = pre + '_12CO21.noiseadd.' + str(n + 1) + '.fits.gz'
    cube13      = pre + '_13CO21.noiseadd.' + str(n + 1) + '.fits.gz'
    temperature  = (6 + np.random.rand()) + 6
    info        = pre + '_noise_' + str(n + 1)
    noise_names = [cube12, cube13, rms12, rms13, mask12]
    lte(files = noise_names, tfloor = temperature, datainfo = info, tx_method = 'cube', onlywrite = ['outn13cube', 'outn13col'])

rms_names1 = [pre + '_noise_' + str(n + 1) + '_cube_n13cube.fits.gz' for n in range(noiseiter)]
rms_names2 = [pre + '_noise_' + str(n + 1) + '_cube_n13col.fits.gz' for n in range(noiseiter)]
noiseout1  = pre + '_noise_rms_cube_n13cube.fits.gz'
noiseout2  = pre + '_noise_rms_cube_n13col.fits.gz'
rms(names = rms_names1, outname = noiseout1)
rms(names = rms_names2, outname = noiseout2)

xname1 = '../lte/30dor_cube_n13cubeerr.fits.gz'
xname2 = '../lte/30dor_cube_n13colerr.fits.gz'
yname1 = '30dor_noise_rms_cube_n13cube.fits.gz'
yname2 = '30dor_noise_rms_cube_n13col.fits.gz'
oname1 = '30dor_noise_25_trandom_n13cube.fits.gz'
oname2 = '30dor_noise_25_trandom_n13col.fits.gz'

histplot(xname = xname1, yname = yname1, snrcut = 0, dolog2d = False, dolog1d = False, nbins = 100, outname = oname1, extrema = [0, 0, 2.1e15, 2.1e15])
histplot(xname = xname2, yname = yname2, snrcut = 0, dolog2d = False, dolog1d = False, nbins = 100, outname = oname2)
