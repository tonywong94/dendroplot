#!/usr/bin/env python2

from run_scimes import *
from calc_phys_props import *
from add_ltemass import *

for line in ['12', '13']:
    label = '30dor_'+line
    if line == '12':
        cubefile = '../12CO.combined.20130113.smooth.flat.fits'
        mom0file = '../mom/30Dor_combined_12CO_dil.mom0.fits'
    else:
        cubefile = '../13CO.combined.20130113.flat.fits'
        mom0file = '../mom/30Dor_combined_13CO_dil.mom0.fits'

    criteria = ['volume']
    #criteria = ['volume', 'luminosity']
    #criteria = ['luminosity']

    #run_scimes(criteria=criteria, label=label, cubefile=cubefile, mom0file=mom0file)
    #explore_dendro(label=label, xaxis='radius', yaxis='v_rms')
    #calc_phys_props(label=label, cubefile=cubefile)

    # Run this after you have run lte.py
    add_ltemass(label=label, n13cub='../lte/30Dor_cube_n13cube.fits.gz',
        n13cub_uc=['../lte/30Dor_noise_25_itrs_rms_n13cube.fits.gz',
        '../lte/30Dor_cube_n13cubeerr.fits.gz'])
    


