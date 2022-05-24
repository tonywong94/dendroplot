#!/usr/bin/env python

from os.path import expanduser, join
from dendroplot.plotting import pltprops

# Which lines
# domos   = ['30Dor_mosaic100_', '30Dor_feather_mosaic100_']
domos   = ['30Dor_feather_mosaic100_']
dolines = ['13', '12']
res     = '1p8'
beam    = 1.75
dv      = 0.1

analdir = expanduser('~/Scratch3/30Dor/analysis/dendro/')

for mos in domos:
    if mos == '30Dor_mosaic100_':
        outdir = analdir + 'mosaic_100'
    else:
        outdir = analdir + 'feather_100'

    for line in dolines:
        catalog = join(outdir, mos+res+'_'+line+'_full_catalog.txt')
        if line == '12':
            pltprops(catalog, dvkms=dv, beam=beam, cmap='jet_r',
                xplot=['rad_pc_dcon','vrms_k','area_pc2','siglum','area_pc2'],
                yplot=['vrms_k','mlumco',  'mlumco', 'sigvir',   'mvir'],
                xlims=[[-1.5,1.5],  [-2,2],  [-1.5,3],  [-1,4],  [-1.5,3]],
                ylims=[[-2,1.5],  [-1,6],    [-1,6],  [-1,4],    [-1,6]],
                pltname=[ 'rdv','dvmlum','areamlum',   'bnd',  'areamv'],
                ccode=[    True,   False,     False,    True,     False],
                doline=[   True,    True,      True,   False,      True],
                panel=[    None,    None,      None,    None,      None],
                colorcodes=['refdist'], lobin_col='salmon', hibin_col='cyan',
                nbin=12, plotdir=join(outdir,'plots'), resolve_output='decon')
        else:
            pltprops(catalog, dvkms=dv, beam=beam, cmap='jet_r',
                xplot=['rad_pc_dcon','vrms_k','area_pc2','siglum','siglte'],
                yplot=['vrms_k','mlumco',  'mlumco','sigvir','sigvir'],
                xlims=[[-1.5,1.5],  [-2,2],  [-1.5,3],  [-1,4],  [-1,4]],
                ylims=[[-2,1.5],  [-1,6],    [-1,6],  [-1,4],  [-1,4]],
                pltname=[ 'rdv','dvmlum','areamlum',   'bnd','bndlte'],
                ccode=[    True,   False,     False,    True,    True],
                doline=[   True,    True,      True,   False,   False],
                panel=[   '(c)',    None,      None,    None,      None],
                colorcodes=['refdist'], lobin_col='salmon', hibin_col='cyan',
                nbin=12, plotdir=join(outdir,'plots'), resolve_output='decon')

