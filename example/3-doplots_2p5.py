#!/usr/bin/env python

import os
from dendroplot.plotting import pltprops

# Which lines
cloud   = '30dor'
dolines = ['12', '13']

# Go to directory with dendrogram output
old_dir = os.getcwd()
os.chdir(os.getcwd()+'/dendro')

for line in dolines:
    catalog = cloud+'_'+line+'_full_catalog.txt'
    if line == '12':
        pltprops(catalog, dvkms=0.1, beam=2.5, 
            xplot=['rad_pc','vrms_k','area_pc2','area_pc2','siglum'],
            yplot=['vrms_k','mlumco',  'mlumco',    'mvir','sigvir'],
            xlims=[[-1,1.5],  [-2,2],  [-1,3.5],  [-1,3.5],  [-1,4]],
            ylims=[[-2,1.5],  [-1,6],    [-1,6],    [-1,6],  [-1,4]],
            pltname=[ 'rdv','dvmlum','areamlum',  'areamv',   'bnd'],
            ccode=[    True,   False,     False,     False,   True], 
            colorcodes=['alpha'], nbin=10)
    else:
        pltprops(catalog, dvkms=0.1, beam=2.5, 
            xplot=['rad_pc','vrms_k','area_pc2','siglum','siglte'],
            yplot=['vrms_k','mlumco',  'mlumco','sigvir','sigvir'],
            xlims=[[-1,1.5],  [-2,2],  [-1,3.5],  [-1,4],  [-1,4]],
            ylims=[[-2,1.5],  [-1,6],    [-1,6],  [-1,4],  [-1,4]],
            pltname=[ 'rdv','dvmlum','areamlum',   'bnd','bndlte'],
            ccode=[    True,   False,     False,     False,   True],
            colorcodes=['alpha'], nbin=10)

os.chdir(old_dir)
