#!/usr/bin/env python

import os
import numpy as np
from dendroplot.plotting import color_code_bin
from astropy.table import Table
import matplotlib.pyplot as plt

# Figure 11 in ApJ paper
# Plot resolved dendro properties vs. distance from R136
# plotx: x-axis is refdist
# ploty: y-axis is sigvir or alpha(lte)
# plotc: colorcode is siglum or siglte

clustonly = False
analdir = 'struct/'
select = 'physprop_resolve'
mos = analdir+'30Dor_feather_mosaic_1p8_12_'+select+'.txt'

xlims   = [0, 320]

for line in ['12', '13']:
    if line == '13':
        mos = mos.replace('_12_','_13_')
        allploty = ['sigvir', 'alphalte']
    else:
        allploty = ['sigvir', 'alpha']
    cllist = mos.replace(select,'clusters')
    pcat = Table.read(mos, format='ascii.ecsv')
    pcat.add_index('_idx')
    mosname = os.path.basename(mos)
    print('\nWorking on', mosname)

    for ploty in allploty:
        for i, y_control in enumerate([False, True]):
            fig, axes = plt.subplots(figsize=(6,4))
            plotx = 'refdist'
            if line == '12':
                plotc = 'siglum'
                zlbl = 'CO-based'
            else:
                plotc = 'siglte'
                zlbl = 'LTE-based'
            if y_control:
                ploty = plotc
            print('\n{} vs {} color coded by {}'.format(ploty,plotx,plotc))
            if clustonly:
                clust = np.loadtxt(cllist, usecols=0, dtype=int)
                goodclust = np.intersect1d(pcat['_idx'], clust)
                pcat = pcat.loc[goodclust]
            valid = (pcat[ploty] > 0) & (pcat[plotc] > 0)
            x, y, z = [pcat[plotx], np.log10(pcat[ploty]), np.log10(pcat[plotc])]
            color_code_bin(x, y, z, colname='$\Sigma$, '+zlbl, cmap='gist_rainbow', 
                           axes=axes, lobin=25, lobin_col='salmon', hibin=75, 
                           hibin_col='cyan', xlims=xlims, zlog=True, median=True, nbin=8)
            choose_q1 = (z < np.nanpercentile(z, 25)) & valid
            print('Number in bottom quartile:',len(pcat[ploty][choose_q1]))
            if line == '13':
                print('Mass in bottom quartile:',np.sum(pcat['mlte'][choose_q1]))
            else:
                print('Mass in bottom quartile:',np.sum(pcat['mlumco'][choose_q1]))
            choose_q2 = (z > np.nanpercentile(z, 75)) & valid
            print('Number in top quartile:',len(pcat[ploty][choose_q2]))
            if line == '13':
                print('Mass in top quartile:',np.sum(pcat['mlte'][choose_q2]))
            else:
                print('Mass in top quartile:',np.sum(pcat['mlumco'][choose_q2]))
            axes.set_xlabel('Distance from R136 ['+str(x.unit)+']', fontsize=12)
            axes.set_xlim(xlims)
            if i == 0:
                if ploty == 'sigvir':
                    axes.set_ylabel(r'log $\Sigma$, virial'+' ['+str(y.unit)+']', size=12)
                else:
                    axes.set_ylabel(r'log $\alpha$, '+zlbl, fontsize=12)
            else:
                axes.set_ylabel(r'log $\Sigma$, '+zlbl+' ['+str(y.unit)+']', size=12)
            axes.text(0.97,0.92,'$^{'+line+'}$CO', size=14, ha='right', transform=axes.transAxes)
            if clustonly:
                axes.text(0.97,0.85,'clusters', size=14, ha='right', transform=axes.transAxes)
            # vertical line at 200" distance
            axes.axvline(200, c='b', ls='--', zorder=-1)
            # horizontal line at log alpha = 0
            if ploty.startswith('alpha'):
                axes.axhline(0, c='k', ls=':', zorder=-2)
            if clustonly:
                pdfname = plotx+'_'+ploty+'_'+mosname.split('_')[1]+'_'+line+'_clust.pdf'
            else:
                pdfname = plotx+'_'+ploty+'_'+mosname.split('_')[1]+'_'+line+'.pdf'
            print('Output to',pdfname)
            plt.savefig(pdfname, bbox_inches='tight')
            plt.close()

