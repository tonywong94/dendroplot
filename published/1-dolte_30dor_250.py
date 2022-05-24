#!/usr/bin/env python

import os
from dendroplot.lte import lte

# Directory where input images reside.
datadir = os.path.expanduser('~/Scratch3/30Dor/products_250')
# Directory to write output files.
workdir = os.path.expanduser('~/Scratch3/30Dor/analysis/lte')

tmin = 8
method = 'peak'
#method = 'cube'

doarr = ['7meter', '12meter']
domos = ['30Dor_mosaic', '30Dor_feather_mosaic']

for mos in domos:
    for arr in doarr:

        indir = datadir
        if arr == '7meter':
            outdir = workdir + '_7m'
        else:
            outdir = workdir

        incube12 = indir + '/' + mos + '_12CO_' + arr + '.pbcor.K.fits.gz'
        incube13 = indir + '/' + mos + '_13CO_' + arr + '.pbcor.K.fits.gz'
        inrms12  = indir + '/' + mos + '_12CO_' + arr + '.rms.K.fits.gz'
        inrms13  = indir + '/' + mos + '_13CO_' + arr + '.rms.K.fits.gz'
        inmask12 = indir + '/' + mos + '_12CO_' + arr + '.mask.fits.gz'
        lte_names   = [incube12, incube13, inrms12, inrms13, inmask12]

        old_dir = os.getcwd() # returns absolute path
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        try:
            os.chdir(outdir)
            lte(files=lte_names, tfloor=tmin, datainfo=mos, tx_method=method)
        finally:
            os.chdir(old_dir)

