#!/usr/bin/env python2.7

import os
import sys
from run_scimes import *
from calc_phys_props import *
from colorcode import *
from add_ltemass import add_ltemass

workdir = '../props'
redo = 'y'   # whether to regenerate dendrogram.hdf file

doclouds = ['A439', 'GMC1', 'GMC104', 'N59C', '30Dor', 'PCC']
dolines = ['12', '13']

for cloud in doclouds:
    if cloud in ['A439', 'GMC1', 'GMC104', 'N59C']:
        type = '_12m7m'
    elif cloud == '30Dor':
        type = '21_12mAPEX'
    elif cloud == 'PCC':
        type = '21_12mTP'
    else:
        print('Unrecognized cloud!!')
    if cloud in ['A439', 'GMC1', 'GMC104', 'N59C']:
        type13 = '_12m'
        ascale = 2.5
    else:
        type13 = type
        ascale = 4

    for line in dolines:
        label = cloud+'_'+line
        if line == '12':
            cubefile = '../'+cloud+'_'+line+'CO'+type+'.image.fits.gz'
            mom0file = '../'+cloud+'_'+line+'CO'+type+'_dil.mom0.fits'
        else:
            cubefile = '../'+cloud+'_'+line+'CO'+type13+'.image.fits.gz'
            mom0file = '../'+cloud+'_'+line+'CO'+type13+'_dil.mom0.fits'

        criteria = ['volume']

        old_dir = os.getcwd() # returns absolute path
        if not os.path.isdir(workdir):
            os.makedirs(workdir)
        try:
            os.chdir(workdir)
            run_scimes(criteria=criteria, label=label, cubefile=cubefile, 
                mom0file=mom0file, redo=redo)
            calc_phys_props(label=label, cubefile=cubefile, efloor=0.1,
                alphascale=ascale)
            colorcode(label=label, cubefile=cubefile, mom0file=mom0file, 
                outdir='plots', types=['v_cen','v_rms'])
            #sys.exit("Stopping here")
        finally:
            os.chdir(old_dir)

