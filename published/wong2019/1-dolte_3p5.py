#!/usr/bin/env python

import os
from dendroplot.lte import lte

## Compute the LTE column density.

# Directory where input images reside.
indir = 'images'
# Directory to write output files.
outdir = 'lte'

tmin = 6
method = 'peak'

doclouds = ['A439', 'GMC1', 'GMC104', 'N59C', '30Dor', 'PCC']

for cloud in doclouds:
    if cloud in ['A439', 'GMC1', 'GMC104', 'N59C']:
        type = '_3p5as'
    else:
        type = '21_3p5as'

    incube12 = indir + '/' + cloud + '_12CO' + type + '.pbcor.fits.gz'
    incube13 = indir + '/' + cloud + '_13CO' + type + '.pbcor.fits.gz'
    inrms12  = indir + '/' + cloud + '_12CO' + type + '_dil.rms.fits.gz'
    inrms13  = indir + '/' + cloud + '_13CO' + type + '_dil.rms.fits.gz'
    inmask12 = indir + '/' + cloud + '_12CO' + type + '_dil.mask.fits.gz'
    lte_names   = [incube12, incube13, inrms12, inrms13, inmask12]

    lte(files=lte_names, tfloor=tmin, outname=cloud, outdir=outdir, tx_method=method)

