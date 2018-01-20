#!/usr/bin/env python

import os
from pltprops import pltprops

# Directory to write output files.
workdir = '/Volumes/Scratch3/tonywong/CLOUD/analysis/dendro/'

#doclouds = ['A439', 'GMC1', 'GMC104', 'N59C', '30Dor', 'PCC']
doclouds = ['PCC']

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
        for line in ['12','13']:
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

# Exit statement
os.chdir(old_dir)
