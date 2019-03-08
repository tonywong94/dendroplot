#!/usr/bin/env python

import numpy as np
import os
import re
import csv
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import ticker
from astropy import units as u
from astropy import constants as const
from astropy.table import Table, Column
from matplotlib.colors import Normalize, LogNorm
#from mpl_toolkits.axes_grid1 import make_axes_locatable


# General Scatter Plot
def plot_ecsv(ecsvfile, xaxis, yaxis, zaxis=None, col='g', mark='o', mec='face', 
           zorder=-5, msize=6, linfit=None, label=None, leaves=False, **kwargs):
    cat = Table.read(ecsvfile, format='ascii.ecsv')
    xdata = cat[xaxis]
    ydata = cat[yaxis]
    if zaxis is not None:
        zdata = cat[zaxis]
    goodidx = (xdata>0) & (ydata>0)
    if leaves:
        leavelist=re.sub('physprop\w*', 'leaves', ecsvfile)
        with open(leavelist, 'r') as f:
            reader=csv.reader(f, delimiter=' ')
            a = zip(*reader)
        idcs = map(int,a[0])
        idsel = [val for val in idcs if goodidx[val]]
        xdata = xdata[idsel]
        ydata = ydata[idsel]
        if zaxis is not None:
            zdata = zdata[idsel]
    else:
        xdata = xdata[goodidx]
        ydata = ydata[goodidx]
        if zaxis is not None:
            zdata = zdata[goodidx]
    print('Plotting {} vs {} from file {}'.format(yaxis,xaxis,ecsvfile))
    # Specified colors
    if zaxis is None:
        axes.scatter(xdata, ydata, marker=mark, c=col, edgecolors=mec, 
            zorder=zorder, s=msize, linewidths=0.1, label=label, **kwargs)
    else:
        axes.scatter(xdata, ydata, marker=mark, c=zdata, edgecolors=mec, 
            zorder=zorder, s=msize, linewidths=0.1, label=label, **kwargs)
    axes.set_xscale('log')
    axes.set_yscale('log')
    axes.set_aspect('equal')
    mapping = { 'mvir':'virial mass', 
        'mlumco':'CO-based mass',
        'mlte':'LTE mass',
        'flux12':'Integrated $^{12}$CO flux',
        'flux13':'Integrated $^{13}$CO flux',
        'siglum':'$\Sigma$, CO-based',
        'siglte':'$\Sigma$, LTE-based',
        'sigvir':'$\Sigma$, virial',
        'rad_pc':'spherical radius',
        'vrms_k':'rms linewidth',
        'area_pc2':'projected area',
       }
    axlbl=['','']
    for i, axis in enumerate([xaxis, yaxis]):
        if axis in mapping.keys():
            axlbl[i] = mapping[axis]
        else:
            axlbl[i] = axis
    axes.set_xlabel('log '+axlbl[0]+' ['+str(xdata.unit)+']')
    axes.set_ylabel('log '+axlbl[1]+' ['+str(ydata.unit)+']')
    # Do a linear regression fit if requested
    if linfit is not None:
        sorted=np.argsort(xdata)
        m, b, rval, pval, std_err = stats.linregress(np.log10(xdata[sorted]),
            np.log10(ydata[sorted]))
        xmod = np.logspace(-3,6,100)
        ymod = b+m*np.log10(xmod)
        axes.plot(xmod, 10**(ymod), linestyle='--', color=linfit)
        axes.text(0.03,0.95,'slope = $%4.2f$' % m,size=9,transform=axes.transAxes)
        axes.text(0.03,0.90,'intercept = $%4.2f$' % b,size=9,transform=axes.transAxes)
    return

# -------------------------------------------------------------------------------

# Main program
def comp_props(dolines, dotypes=['sp8med'], clouds=None, markers=None,
            indir=None, leaves=False, cmap_name='gist_rainbow', mec='white',
            xplot=['rad_pc'],
            yplot=['vrms_k'],
            xlims=[[-1,1.5]],
            ylims=[[-2,1.5]],
            pltname=['rdv'], slope=[0.5], inter=[0],
            pad  = [0.03], magmacsv='islands.sageco.csv'):

    global axes
    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)
    cmap = plt.get_cmap(cmap_name)

    # Read table of cloud-averaged properties
    isle_id = {'N55':441, 'N11B':395, 'N159':94, 'GMC225':2, 'N166':50, 'N171':32,
          'N206':60, 'N206D':46, 'GMC2':259, 'GMC55':224, '30Dor':184,
          'PCC':39, 'N113':127, 'N103B':188, '30DorC':169, 'A439':29,
          'GMC104':57, 'GMC1':165, 'N59C':392}
    magmatab = Table.read(magmacsv, format='ascii.ecsv')
    keep = []
    for clname in clouds:
        keep.append(isle_id[clname]-1)
    cldtab = magmatab[keep]
#     cldtab = Table(names=('cld', 'sp8med', 'med_co', 'avg_st', 'sp24med'), 
#         dtype=('S6', 'f4', 'f4', 'f4', 'f4'))
#     cldtab.add_row(('30Dor', 36.8732, 2.0518, 0.8101, 386.8365))
#     cldtab.add_row(('N59C',   8.1936, 3.0875, 0.3671,  12.8014))
#     cldtab.add_row(('A439',   2.6004, 3.2009, 0.6587,   1.2010))
#     cldtab.add_row(('GMC104', 1.5529, 4.4368, 1.7162,   0.6574))
#     cldtab.add_row(('GMC1',   1.1792, 2.0274, 0.3059,   0.5231))
#     cldtab.add_row(('PCC',    0.7974, 1.4606, 0.5502,   0.1617))
#     # Convert stellar surf density to Msol/pc^2
#     cldtab['avg_st'] = cldtab['avg_st'] * 62.5

    # Generate plots
    for i in range(len(xplot)):
        for line in dolines:
            for type in dotypes:
                fig, axes = plt.subplots()
                for j, clname in enumerate(clouds):
                    dir = indir.replace('CLOUD', clname)
                    if os.path.isfile(dir+clname+'_'+line+'_physprop_add.txt'):
                        infile = dir+clname+'_'+line+'_physprop_add.txt'
                    else:
                        infile = dir+clname+'_'+line+'_physprop.txt'
                    if type in ['8um_avg', 'siglum', 'siglte']:
                        if 'sig' in type:
                            norm = LogNorm(vmin=10, vmax=1000)
                        else:
                            norm = LogNorm(vmin=1, vmax=100)
                        plot_ecsv(infile, xplot[i], yplot[i],
                            zaxis=type, cmap=cmap, mark=markers[j], msize=9,
                            zorder=i, label=clname, leaves=leaves, norm=norm)               
                    else:
                        if type == 'comean' or type == 'comax':
                            norm = Normalize(vmin=min(cldtab[type]), vmax=max(cldtab[type]))
                        else:
                            norm = LogNorm(vmin=min(cldtab[type]), vmax=max(cldtab[type]))
                        colr=np.array(cmap(norm(cldtab[type][j])))
                        plot_ecsv(infile, xplot[i], yplot[i], 
                            col=colr, mark=markers[j], msize=9, zorder=i, 
                            label=clname, leaves=leaves)
                axes.set_xlim(10**xlims[i][0], 10**xlims[i][1])
                axes.set_ylim(10**ylims[i][0], 10**ylims[i][1])
                logfmt = ticker.LogFormatterExponent(base=10.0, labelOnlyBase=True)
                axes.xaxis.set_major_formatter(logfmt)
                axes.yaxis.set_major_formatter(logfmt)
                axes.minorticks_off()
                # Reference slopes
                xmod = np.logspace(xlims[i][0], xlims[i][1], 100)
                ymod = 10**inter[i] * xmod**slope[i]
                if xplot[i] == 'rad_pc' and yplot[i] == 'vrms_k':
                    axes.plot(xmod, ymod, linestyle='-', color='r', lw=4, alpha=0.5, 
                        zorder=-1)
                    axes.text(10**(xlims[i][1]-0.05), 10**(ylims[i][1]-0.8), 'S87', 
                        horizontalalignment='right', color='r', rotation=30)
                else:
                    axes.plot(xmod, ymod, '--', marker=None, color='k')
                # Lines of constant external pressure
                if xplot[i].startswith('sig') and yplot[i] == 'sigvir':
                    ymod = xmod + (20/(3*np.pi*21.1))*1.e4/xmod
                    axes.plot(xmod, ymod, linestyle='-', color='g', lw=1)
                    axes.text(10**-0.97, 10**3.5, '$P_{ext}$ = $10^4$ cm$^{-3}$ K', 
                        color='g', rotation=-45)
                    ymod2 = xmod + (20/(3*np.pi*21.1))*1.e2/xmod
                    axes.plot(xmod, ymod2, linestyle='-', color='m', lw=1)
                    axes.text(10**-0.97, 10**2.1, '$P_{ext}$ = $10^2$ cm$^{-3}$ K', 
                        color='m', rotation=-45)
                    axes.text(10**-0.95, 10**-0.25, 'Virial Eq', 
                        color='k', rotation=45)
                # Legend and colorbar
                plt.legend(loc='lower right',fontsize='small',numpoints=1,markerscale=2)
                cax = fig.add_axes([pad[i]+0.7, 0.11, 0.02, 0.77])
                formatter = ticker.LogFormatter(10, labelOnlyBase=False, minor_thresholds=(4,3))
                if type == 'comean':
                    cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
                else:
                    cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, 
                        format=formatter, ticks=[1,2,5,10,20,50,100,200,500,1000])
                cbar.ax.tick_params(labelsize=9)
                if type == 'sp8med':
                    cbar.set_label('cloud median 8$\mu$m intensity [MJy/sr]', rotation=90)
                elif type == 'comean':
                    cbar.set_label('cloud mean CO intensity [K km/s]', rotation=90)
                elif type == 'comax':
                    cbar.set_label('cloud max CO intensity [K km/s]', rotation=90)
                elif type == 'stmean':
                    cbar.set_label('cloud mean $\Sigma_{*}$ [$M_\odot$ $pc^{-2}$]', 
                        rotation=90)
                elif type == 'sp24med':
                    cbar.set_label('cloud median 24$\mu$m intensity [MJy/sr]', rotation=90)
                elif type == '8um_avg':
                    cbar.set_label('local mean 8$\mu$m intensity [MJy/sr]', rotation=90)
                elif type == 'siglum':
                    cbar.set_label('local mean $\Sigma_{mol}$ [$M_\odot$ $pc^{-2}$]', rotation=90)
                elif type == 'siglte':
                    cbar.set_label('local mean $\Sigma_{LTE}$ [$M_\odot$ $pc^{-2}$]', rotation=90)
                plt.savefig('comp_'+line+'_'+pltname[i]+'_'+type+'.pdf', 
                    bbox_inches='tight')
                plt.close()
    return

