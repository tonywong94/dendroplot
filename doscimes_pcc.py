#!/usr/bin/env python

from run_scimes import *
from calc_phys_props import *

for line in ['12', '13']:
    label = 'pcc_'+line
    if line == '12':
        cubefile = '../PCC_12mTP_12CO21.image.fits'
        mom0file = '../PCC_12mTP_12CO_dil.mom0.fits'
    else:
        cubefile = '../PCC_12mTP_13CO21.image.fits'
        mom0file = '../PCC_12mTP_13CO_dil.mom0.fits'

    criteria = ['volume']

    #run_scimes(criteria=criteria, label=label, cubefile=cubefile, mom0file=mom0file)
    calc_phys_props(label=label, cubefile=cubefile)
    #explore_dendro(label=label, xaxis='radius', yaxis='v_rms')
