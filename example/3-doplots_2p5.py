#!/usr/bin/env python

import os
import sys
sys.path.append(os.path.expanduser('~/Work/bin/py-package/lmc_alma_analysis/'))
from pltprops import pltprops

# Which lines
cloud   = '30dor'
dolines = ['12', '13']
beam = 2.5

# Go to directory with dendrogram output
old_dir = os.getcwd()
os.chdir(os.getcwd()+'/dendro')

for line in dolines:
    label = cloud+'_'+line
    if line == '12':
        pltprops(label, dvkms=0.1, beam=beam, 
            xplot=['rad_pc','vrms_k','area_pc2','area_pc2','siglum'],
            yplot=['vrms_k','mlumco',  'mlumco',    'mvir','sigvir'],
            xlims=[[-1,1.5],  [-2,2],  [-1,3.5],  [-1,3.5],  [-1,4]],
            ylims=[[-2,1.5],  [-1,6],    [-1,6],    [-1,6],  [-1,4]],
            pltname=[ 'rdv','dvmlum','areamlum',  'areamv',   'bnd'])
    else:
        pltprops(label, dvkms=0.1, beam=beam, 
            xplot=['rad_pc','vrms_k','area_pc2','siglum','siglte'],
            yplot=['vrms_k','mlumco',  'mlumco','sigvir','sigvir'],
            xlims=[[-1,1.5],  [-2,2],  [-1,3.5],  [-1,4],  [-1,4]],
            ylims=[[-2,1.5],  [-1,6],    [-1,6],  [-1,4],  [-1,4]],
            pltname=[ 'rdv','dvmlum','areamlum',   'bnd','bndlte'])

os.chdir(old_dir)
