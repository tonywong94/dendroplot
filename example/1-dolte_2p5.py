#!/usr/bin/env python

import os
from dendroplot.lte import lte

# Directory where input images reside.
indir = os.getcwd()

# Directory to write output files.
outdir = os.getcwd()+'/lte'
if not os.path.isdir(outdir):
    os.makedirs(outdir)

tmin = 6
method = 'peak'

cloud = '30dor'
type  = '21_2p5as'

incube12 = indir + '/' + cloud + '_12CO' + type + '.pbcor.fits.gz'
incube13 = indir + '/' + cloud + '_13CO' + type + '.pbcor.fits.gz'
inrms12  = indir + '/' + cloud + '_12CO' + type + '.rms.fits.gz'
inrms13  = indir + '/' + cloud + '_13CO' + type + '.rms.fits.gz'
inmask12 = indir + '/mom/' + cloud + '_12CO' + type + '_dil.mask.fits.gz'
lte_names   = [incube12, incube13, inrms12, inrms13, inmask12]

lte(files=lte_names, tfloor=tmin, datainfo=outdir+cloud, tx_method=method)

