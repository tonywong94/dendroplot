#!/usr/bin/env python

import os
from pltprops import pltprops
from comp_props import comp_props

# Directories to write output files.
workdir = '/Volumes/Scratch3/tonywong/CLOUD/analysis/dendro/'
compdir = '/Volumes/Scratch3/tonywong/repository/comp_props/'

# How many clouds
doclouds = ['30Dor', 'N59C', 'A439', 'GMC104', 'GMC1', 'PCC']
markers = ['o', 'v', '^', 's', '<', 'D']

# Which lines
dolines = ['12','13']

# Which quantities to color by in comp_props
colorcodes = ['med_8u', 'med_co', 'avg_st', 'med_24', '8um_avg', 'siglum']

# ----------------------------------------------------------------------------------

old_dir = os.getcwd() # returns absolute path
for cloud in doclouds:
    outdir = workdir.replace('CLOUD', cloud)
    try:
        os.chdir(outdir)
    except OSError:
        print('Could not change directory to {}'.format(outdir))
    else:
        if cloud in ['A439', 'GMC1', 'GMC104', 'N59C']:
            f12 = 115.2712
            f13 = 110.2013
        else:
            f12 = 230.538
            f13 = 220.399
        for line in dolines:
            label = cloud+'_'+line
            if line == '12':
                pltprops(label, fghz=f12, dvkms=0.2, beam=3.5,
                    xplot=['rad_pc', 'vrms_k', 'area_pc2','area_pc2','siglum',   '8um_avg'],
                    yplot=['vrms_k', 'mlumco',   'mlumco',    'mvir','sigvir',   'sigvir'],
                    xlims=[[-1,1.5],   [-2,2],     [-1,3],    [-1,3],  [-1,4], [-1,1]],
                    ylims=[[-2,1.5],   [-1,5],     [-1,6],    [-1,6],  [-1,4], [-1,4]],
                    pltname=[ 'rdv', 'dvflux', 'areaflux',  'areamv',   'bnd',   'sfnorm'])
            else:
                pltprops(label, fghz=f13, dvkms=0.5, beam=3.5,
                    xplot=['rad_pc',   'vrms_k', 'area_pc2' ],
                    yplot=['vrms_k',   'mlumco',   'mlumco' ],
                    xlims=[[-1,1.5],     [-2,2],     [-1,3] ],
                    ylims=[[-2,1.5], [-1.5,4.5],     [-1,4] ],
                    pltname=[ 'rdv',   'dvflux', 'areaflux' ])
try:
    os.chdir(compdir)
except OSError:
    print('Could not change directory to {}'.format(compdir))
else:
        comp_props(dolines, dotypes=colorcodes, clouds=doclouds, 
            markers=markers, indir=workdir, leaves=False,
            xplot=['rad_pc', 'area_pc2', 'siglum', '8um_avg'],
            yplot=['vrms_k', 'mlumco'  , 'sigvir', 'sigvir'],
            xlims=[[-1,1.5], [-1,3], [-1,4], [-0.5,3.5]],
            ylims=[[-2,1.5], [-1,6], [-1,4], [-1,4]],
            pltname=['rdv', 'areaflux' , 'bnd', 'sfnorm'],
            slope=[0.5, 1, 1, 1],
            pad=[0.03, 0.01, 0.13, 0.13])
# Exit statement
os.chdir(old_dir)
