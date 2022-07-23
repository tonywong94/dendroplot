#!/usr/bin/env python

import os
from datetime import datetime
from dendroplot import run_dendro
from dendroplot.analysis import calc_phys_props, find_clusters, define_asgn
from dendroplot.plotting import colorcode
from dendroplot.lte import add_ltemass
import warnings
warnings.filterwarnings("ignore")

# 04may2019: Use co13toh2 = 3.0e6 (Mizuno:10, Fujii:14), alphascale=2.4 for
# CO(1-0) (Hughes:09) and alphascale=3 for CO(2-1) (2-1/1-0 = 0.8).

## Run the dendrogram analysis and calculate physical properties

redo = 'n'   # whether to regenerate dendrogram.hdf file

start = datetime.now()
old_dir = os.getcwd() # returns absolute path
datadir = old_dir + '/images/'
ltedir  = old_dir + '/lte/'
outdir  = 'dendro/'
if not os.path.isdir(outdir):
    os.makedirs(outdir)

doclouds = ['A439', 'GMC1', 'GMC104', 'N59C', '30Dor', 'PCC']
dolines  = ['12', '13']
anclbl   = '8um_avg'
ancimg   = 'CLOUD_8um_3p5as.image.fits.gz'

try:
    os.chdir(outdir)
    for cloud in doclouds:
        if cloud in ['A439', 'GMC1', 'GMC104', 'N59C']:
            ascale = 2.4
            type = '_3p5as'
        else:
            ascale = 3
            type = '21_3p5as'

        ancfil = datadir + ancimg.replace('CLOUD', cloud)

        for line in dolines:
            label = cloud+'_'+line
            cubefile = datadir + cloud+'_'+line+'CO'+type+'.addnse.fits.gz'
            cub12fil = datadir + cloud+'_12CO'+type+'.pbcor.fits.gz'
            cub13fil = datadir + cloud+'_13CO'+type+'.pbcor.fits.gz'
            mom0file = datadir + cloud+'_'+line+'CO'+type+'_dil.mom0.fits.gz'
            rmsfile  = datadir + cloud+'_'+line+'CO'+type+'_dil.rms.fits.gz'

            run_dendro(label=label, cubefile=cubefile, mom0file=mom0file, redo=redo)
            find_clusters(criteria=['volume'], label=label, cubefile=cubefile)
            calc_phys_props(label=label, cubefile=cubefile, efloor=0.05,
                alphascale=ascale, copbcor=cub12fil, conoise=rmsfile, 
                ancfile=ancfil, anclabel=anclbl)
            if os.path.isfile(ltedir+cloud+'_peak_n13cube.fits.gz'):
                add_ltemass(label=label, n13cub=ltedir+cloud+'_peak_n13cube.fits.gz', 
                            i12cub=cub12fil, i13cub=cub13fil, co13toh2=3e6, efloor=0.05,
                            n13cub_uc=ltedir+cloud+'_peak_n13cubeerr.fits.gz')
            colorcode(label=label, cubefile=cubefile, mom0file=mom0file, 
                      types=['v_cen'])
            colorcode(label=label, cubefile=cubefile, mom0file=mom0file,
                      types=['sigvir','alpha','vrms_k'], table='physprop',
                      lognorm=True)
            asgn = define_asgn(cubefile,label+'_dendrogram.hdf5', label_out=label)

finally:
    os.chdir(old_dir)

