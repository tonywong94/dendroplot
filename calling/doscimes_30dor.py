#!/usr/bin/env python2.7

import os
from run_scimes import *
from calc_phys_props import *
from add_ltemass import add_ltemass

workdir = '../props'
redo = 'n'   # whether to regenerate dendrogram.hdf file

for line in ['12', '13']:
    label = '30dor_'+line
    if line == '12':
        cubefile = '../12CO.combined.20130113.smooth.flat.fits.gz'
        mom0file = '../mom/30Dor_combined_12CO_dil.mom0.fits.gz'
    else:
        cubefile = '../13CO.combined.20130113.flat.fits.gz'
        mom0file = '../mom/30Dor_combined_13CO_mk12.mom0.fits.gz'
    n13cube = '30Dor_peak_n13cube.fits.gz'
    n13cube_uc1 = '30Dor_noise_rms_peak_n13cube.fits.gz'
    n13cube_uc2 = '30Dor_peak_n13cubeerr.fits.gz'

    criteria = ['volume']

    old_dir = os.getcwd() # returns absolute path
    if not os.path.isdir(workdir):
        os.makedirs(workdir)
    try:
        os.chdir(workdir)
        run_scimes(criteria=criteria, label=label, cubefile=cubefile, 
            mom0file=mom0file, redo=redo)
        calc_phys_props(label=label, cubefile=cubefile)
        # Must have already run lte.py
        if os.path.exists('../lte/'+n13cube):
            add_ltemass(label=label, n13cub='../lte/'+n13cube,
                n13cub_uc=['../lte/'+n13cube_uc1,'../lte/'+n13cube_uc2])
    finally:
        os.chdir(old_dir)

