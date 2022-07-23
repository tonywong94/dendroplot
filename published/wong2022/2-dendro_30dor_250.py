#!/usr/bin/env python

import os
from datetime import datetime
from dendroplot import run_dendro
from dendroplot.analysis import calc_phys_props, find_clusters, define_asgn
from dendroplot.plotting import colorcode
from dendroplot.lte import add_ltemass
from dendroplot.analysis.calc_phys_props import refdist_redo

## Run the dendrogram analysis and calculate physical properties

outdir  = 'dendro/'
redo = 'n'   # whether to regenerate dendrogram.hdf file

domos   = ['30Dor_feather_mosaic_']
res     = '1p8'
dolines = ['13', '12']

start = datetime.now()
old_dir = os.getcwd() # returns absolute path
datadir = old_dir + '/images/'

for mos in domos:
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    try:
        os.chdir(outdir)

        cubfile = ['', '']
        momfile = ['', '']
        for i, line in enumerate(dolines):
            cubfile[i] = datadir + mos + line + 'CO_12meter.image.fits.gz'
            cub12file  = datadir + mos + '12CO_12meter.pbcor.K.fits.gz'
            cub13file  = datadir + mos + '13CO_12meter.pbcor.K.fits.gz'
            conoise    = datadir + mos + '12CO_12meter.rms.K.fits.gz'
            momfile[i] = datadir + mos + line + 'CO_12meter.mom0.fits.gz'
            label = mos+res+'_'+line
            run_dendro(label=label, cubefile=cubfile[i], mom0file=momfile[i], redo=redo)
            find_clusters(criteria=['volume'], label=label, cubefile=cubfile[i])
            calc_phys_props(label=label, cubefile=cubfile[i],
                            copbcor=cub12file, conoise=conoise, 
                            alphascale=1.6, efloor=0.1,
                            refpos=[84.67625,-69.100917])
            if os.path.exists('lte'):
                n13cube    = mos + 'peak_n13cube.fits.gz'
                n13cube_uc = mos + 'peak_n13cubeerr.fits.gz'
                add_ltemass(label=label, n13cub='lte/'+n13cube, 
                            i12cub=cub12file, i13cub=cub13file, efloor=0.1,
                            n13cub_uc=['lte/'+n13cube_uc], co13toh2=3e6)
            colorcode(label=label, cubefile=cubfile[i], mom0file=momfile[i], 
                      types=['v_cen'])
            colorcode(label=label, cubefile=cubfile[i], mom0file=momfile[i],
                      types=['sigvir','alpha','vrms_k'], table='physprop',
                      lognorm=True)
            asgn = define_asgn(cubfile[i],label+'_dendrogram.hdf5', label_out=label)
    finally:
        os.chdir(old_dir)

end = datetime.now()
time_taken = end - start
print('Execution Time: ',time_taken)
