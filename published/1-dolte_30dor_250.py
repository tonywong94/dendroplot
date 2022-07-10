#!/usr/bin/env python

import os
from dendroplot.lte import lte

## Compute the LTE column density.

# Directory where input images reside.
indir = 'images'
# Directory to write output files.
outdir = 'lte'

tmin = 8
method = 'peak'

doarr = ['12meter']
domos = ['30Dor_feather_mosaic']

for mos in domos:
    for arr in doarr:

        incube12 = indir + '/' + mos + '_12CO_' + arr + '.pbcor.K.fits.gz'
        incube13 = indir + '/' + mos + '_13CO_' + arr + '.pbcor.K.fits.gz'
        inrms12  = indir + '/' + mos + '_12CO_' + arr + '.rms.K.fits.gz'
        inrms13  = indir + '/' + mos + '_13CO_' + arr + '.rms.K.fits.gz'
        inmask12 = indir + '/' + mos + '_12CO_' + arr + '.mask.fits.gz'
        lte_names   = [incube12, incube13, inrms12, inrms13, inmask12]

        lte(files=lte_names, tfloor=tmin, outname=mos, outdir=outdir, tx_method=method)
