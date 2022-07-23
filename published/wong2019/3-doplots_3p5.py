#!/usr/bin/env python

from dendroplot.plotting import pltprops

## Plot the property correlations for individual GMCs

analdir = 'dendro/'

# Which clouds
doclouds = ['A439', 'GMC1', 'GMC104', 'N59C', '30Dor', 'PCC']

# Which lines
dolines = ['12','13']

# ----------------------------------------------------------------------------------

for cloud in doclouds:
    if cloud == '30Dor':
        delv = 0.5
    else:
        delv = 0.2
    for line in dolines:
        catalog = analdir+cloud+'_'+line+'_full_catalog.txt'
        if line == '12':
            pltprops(catalog, dvkms=delv, beam=3.5, plotdir=analdir+'plots/',
                xplot=['rad_pc','vrms_k','area_pc2','area_pc2','siglum'],
                yplot=['vrms_k','mlumco',  'mlumco',    'mvir','sigvir'],
                xlims=[[-1,1.5],  [-2,2],  [-1,3.5],  [-1,3.5],  [-1,4]],
                ylims=[[-2,1.5],  [-1,6],    [-1,6],    [-1,6],  [-1,4]],
                pltname=[ 'rdv','dvmlum','areamlum',  'areamv',   'bnd'],
                doline=[   True,    True,      True,      True,    True],
                ccode=[    True,   False,     False,     False,    True],
                colorcodes=['8um_avg'], nbin=10)
        else:
            pltprops(catalog, dvkms=delv, beam=3.5, plotdir=analdir+'plots/',
                xplot=['rad_pc','vrms_k','area_pc2','siglum','siglte','flux12' ],
                yplot=['vrms_k','mlumco',    'mlte','sigvir','sigvir','flux13' ],
                xlims=[[-1,1.5],  [-2,2],  [-1,3.5],  [-1,4],  [-1,4], [-1.5,5]],
                ylims=[[-2,1.5],  [-1,6],    [-1,6],  [-1,4],  [-1,4], [-2,4.5]],
                pltname=[ 'rdv','dvmlum','areamlte',   'bnd','bndlte','linerat'],
                doline=[   True,    True,      True,    True,    True,    False],
                ccode=[    True,   False,     False,    True,    True,     False],
                colorcodes=['8um_avg'], nbin=10)
