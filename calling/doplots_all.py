#!/usr/bin/env python

import os
from pltprops import pltprops

workdir = '../props'
doclouds = ['A439', 'GMC1', 'GMC104', 'N59C', '30Dor', 'PCC']

old_dir = os.getcwd() # returns absolute path
if not os.path.isdir(workdir):
    os.makedirs(workdir)
try:
    os.chdir(workdir)
    for cloud in doclouds:
        if cloud in ['A439', 'GMC1', 'GMC104', 'N59C']:
            f12 = 115.2712
            f13 = 110.2013
        else:
            f12 = 230.538
            f13 = 220.399
        for line in ['12','13']:
            label = cloud+'_'+line
            if line == '12':
                pltprops(label, fghz=f12, dvkms=0.2, beam=2,
                    xplot=['rad_pc', 'vrms_k', 'area_pc2','area_pc2','siglum',   'mlumco'],
                    yplot=['vrms_k', 'mlumco',   'mlumco',    'mvir','sigvir',     'mvir'],
                    xlims=[[-1,1.5],   [-2,2],     [-1,3],    [-1,3],  [-1,4], [-0.5,5.5]],
                    ylims=[[-2,1.5],   [-1,5],     [-1,6],    [-1,6],  [-1,4], [-0.5,5.5]],
                    pltname=[ 'rdv', 'dvflux', 'areaflux',  'areamv',   'bnd',   'virial'])
            else:
                pltprops(label, fghz=f13, dvkms=0.5, beam=2.5,
                    xplot=['rad_pc',   'vrms_k', 'area_pc2' ],
                    yplot=['vrms_k',   'mlumco',   'mlumco' ],
                    xlims=[[-1,1.5],     [-2,2],     [-1,3] ],
                    ylims=[[-2,1.5], [-1.5,4.5],     [-1,4] ],
                    pltname=[ 'rdv',   'dvflux', 'areaflux' ])
finally:
    os.chdir(old_dir)
