#!/usr/bin/env python

from os.path import join
from dendroplot.plotting import pltprops

## Plot the cloud property correlations

# Which lines
domos   = ['30Dor_feather_mosaic_']
dolines = ['13', '12']
res     = '1p8'
beam    = 1.75
dv      = 0.25

analdir = 'struct/'

for mos in domos:
    for line in dolines:
        catalog = join(analdir, mos+res+'_'+line+'_full_catalog.txt')
        if line == '12':
            pltprops(catalog, dvkms=dv, beam=beam, cmap='jet_r',
                xplot=['rad_pc_dcon','vrms_k','area_pc2','siglum','area_pc2'],
                yplot=['vrms_k','mlumco',  'mlumco','sigvir',    'mvir'],
                xlims=[[-1.5,1.5],  [-2,2],  [-1.5,3],   [0,4],  [-1.5,3]],
                ylims=[[-2,1.5],  [-1,6],    [-1,6],   [0,4],    [-1,6]],
                pltname=[ 'rdv','dvmlum','areamlum',   'bnd',  'areamv'],
                ccode=[    True,   False,     False,    True,     False],
                doline=[   True,    True,      True,   False,      True],
                panel=[   '(a)',    None,      None,    None,      None],
                colorcodes=['refdist'], lobin_col='salmon', hibin_col='cyan',
                nbin=10, plotdir=join(analdir,'plots'), resolve_output='none')
        else:
            pltprops(catalog, dvkms=dv, beam=beam, cmap='jet_r',
                xplot=['rad_pc_dcon','vrms_k','area_pc2','siglum','siglte'],
                yplot=['vrms_k','mlumco',  'mlumco','sigvir','sigvir'],
                xlims=[[-1.5,1.5],  [-2,2],  [-1.5,3],   [0,4],   [0,4]],
                ylims=[[-2,1.5],  [-1,6],    [-1,6],   [0,4],   [0,4]],
                pltname=[ 'rdv','dvmlum','areamlum',   'bnd','bndlte'],
                ccode=[    True,   False,     False,    True,    True],
                doline=[   True,    True,      True,   False,   False],
                panel=[   '(b)',    None,      None,    None,    None],
                colorcodes=['refdist'], lobin_col='salmon', hibin_col='cyan',
                nbin=10, plotdir=join(analdir,'plots'), resolve_output='none')

