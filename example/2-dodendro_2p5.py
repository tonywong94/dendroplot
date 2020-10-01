#!/usr/bin/env python

import os
import sys
sys.path.append(os.path.expanduser('~/Work/bin/py-package/lmc_alma_analysis/'))
from run_dendro import run_dendro
from calc_phys_props import calc_phys_props
from colorcode import colorcode
from add_ltemass import add_ltemass

# Directory where input images reside.
indir = os.getcwd()
# Directory to write output files.
outdir = os.getcwd()+'/dendro'
if not os.path.isdir(outdir):
    os.makedirs(outdir)

redo = 'y'   # whether to regenerate dendrogram.hdf file

doclouds = ['30dor']
dolines  = ['12', '13']

for cloud in doclouds:
    type = '21_2p5as'
    alphascale = 3

    for line in dolines:
        label = cloud+'_'+line
        cubefile = indir + '/' + cloud+'_'+line+'CO'+type+'.image.fits.gz'
        cub12fil = indir + '/' + cloud+'_12CO'+type+'.pbcor.fits.gz'
        cub13fil = indir + '/' + cloud+'_13CO'+type+'.pbcor.fits.gz'
        mom0file = indir + '/mom/' + cloud+'_'+line+'CO'+type+'_dil.mom0.fits.gz'
        rmsfile  = indir + '/' + cloud+'_'+line+'CO'+type+'.rms.fits.gz'

        criteria = ['volume']

        os.chdir(outdir)
        run_dendro(criteria=criteria, label=label, cubefile=cubefile, 
            mom0file=mom0file, redo=redo)
        calc_phys_props(label=label, cubefile=cubefile, efloor=0.05,
            alphascale=alphascale, copbcor=cub12fil, conoise=rmsfile)
        if os.path.isfile('../lte/'+cloud+'_peak_n13cube.fits.gz'):
            add_ltemass(n13cub='../lte/'+cloud+'_peak_n13cube.fits.gz', 
                i12cub=cub12fil, i13cub=cub13fil, co13toh2 = 3.0e6,
                label=label,
                n13cub_uc='../lte/'+cloud+'_peak_n13cubeerr.fits.gz')
        colorcode(label=label, cubefile=cubefile, mom0file=mom0file, 
            outdir='plots', types=['v_cen','v_rms'])
