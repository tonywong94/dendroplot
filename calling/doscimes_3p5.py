#!/usr/bin/env python

import os
#import sys
from run_scimes import *
from calc_phys_props import *
from colorcode import *
from add_ltemass import add_ltemass

# Directory where input images reside.  CLOUD will be replaced by cloud name.
datadir = '/Volumes/Scratch3/tonywong/CLOUD/cubes/match_3p5/'
# Directory to write output files.
workdir = '/Volumes/Scratch3/tonywong/CLOUD/analysis/dendro/'

redo = 'y'   # whether to regenerate dendrogram.hdf file

#doclouds = ['A439']
doclouds = ['A439', 'GMC1', 'GMC104', 'N59C', '30Dor', 'PCC']
dolines = ['12', '13']
anclbl = '8um_avg'
ancimg = 'CLOUD_8um_3p5as.image.fits.gz'

for cloud in doclouds:
    if cloud in ['A439', 'GMC1', 'GMC104', 'N59C']:
        type = '_3p5as'
    else:
        type = '21_3p5as'
    if cloud in ['A439', 'GMC1', 'GMC104', 'N59C']:
        ascale = 2.5
    else:
        ascale = 4

    indir = datadir.replace('CLOUD', cloud)
    outdir = workdir.replace('CLOUD', cloud)
    ancfil = indir+ancimg.replace('CLOUD', cloud)

    for line in dolines:
        label = cloud+'_'+line
        cubefile = indir+cloud+'_'+line+'CO'+type+'.image.fits.gz'
        mom0file = indir+cloud+'_'+line+'CO'+type+'_dil.mom0.fits'

        criteria = ['volume']

        old_dir = os.getcwd() # returns absolute path
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        try:
            os.chdir(outdir)
            run_scimes(criteria=criteria, label=label, cubefile=cubefile, 
                mom0file=mom0file, redo=redo)
            calc_phys_props(label=label, cubefile=cubefile, efloor=0.1,
                alphascale=ascale, ancfile=ancfil, anclabel=anclbl)
            if os.path.isfile('../lte/'+cloud+'_peak_n13cube.fits.gz'):
                add_ltemass(label=label, n13cub='../lte/'+cloud+'_peak_n13cube.fits.gz', 
                    n13cub_uc='../lte/'+cloud+'_peak_n13cubeerr.fits.gz')
            colorcode(label=label, cubefile=cubefile, mom0file=mom0file, 
                outdir='plots', types=['v_cen','v_rms'])
        finally:
            os.chdir(old_dir)

