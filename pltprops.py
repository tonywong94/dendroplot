#!/usr/bin/env python2.7

import csv
import numpy as np
from scipy import stats
import os
import re
from matplotlib import pyplot as plt
from astropy import units as u
from astropy import constants as const
from astropy.table import Table, Column
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

# General Scatter Plot
def sctplot(xdata, ydata, zdata=None, col='g', mark='o', mec='k', 
           zorder=-5, msize=6, cmap=None, linfit=None, label=None):
    axes.set_xscale('log')
    axes.set_yscale('log')
    axes.set_aspect('equal')
    # Single color plot
    if cmap is None and np.size(col) == 1:
        axes.scatter(xdata, ydata, marker=mark, c=col, edgecolors=mec, 
            zorder=zorder, s=msize, linewidths=1, label=label)
    # Multi-color plot, array of colors
    elif cmap is None:
        for xp, yp, cp in zip(xdata, ydata, col):
            axes.scatter(xp, yp, marker=mark, c=cp, edgecolors=mec, 
                zorder=zorder, s=msize)
    # Continuous color map based on structure number
    elif zdata is None:
        xind = np.linspace(0., 1., len(xdata))
        axes.scatter(xdata, ydata, marker=mark, c=xind, zorder=zorder,
            cmap=cmap, edgecolors=mec, s=msize, label=None)
    # Continuous color map based on property value
    else:
        sc=axes.scatter(xdata, ydata, marker=mark, c=zdata, zorder=zorder,
            cmap=cmap, edgecolors=mec, s=msize, label=None)
        cbar = plt.colorbar(sc)
        cbar.ax.tick_params(labelsize=9) 
        cbar.set_label(label, rotation=90)
    # Do a linear regression fit if requested
    if linfit is not None:
        goodidx = (xdata>0) & (ydata>0)
        xdata = xdata[goodidx]
        ydata = ydata[goodidx]
        sorted=np.argsort(xdata)
        m, b, rval, pval, std_err = stats.linregress(np.log10(xdata[sorted]),
            np.log10(ydata[sorted]))
        xmod = np.logspace(-3,6,100)
        ymod = b+m*np.log10(xmod)
        axes.plot(xmod, 10**(ymod), linestyle='--', color=linfit)
        axes.text(0.03,0.95,'slope = $%4.2f$' % m,size=10,transform=axes.transAxes)
        axes.text(0.03,0.90,'intercept = $%4.2f$' % b,size=10,transform=axes.transAxes)
    return

# -------------------------------------------------------------------------------

def std_overlay(cat, xaxis, yaxis, xlims, ylims):
    # Plot gray shading indicating resolution limits
    rmstorad = 1.91
    radlim = ((avgbeam*rmstorad/np.sqrt(8*np.log(2))) * dist).to(
        u.pc, equivalencies=u.dimensionless_angles())
    dvlim = deltav.value/np.sqrt(8*np.log(2))
    if xaxis == 'rad_pc':
        axes.axvspan(1.e-3, radlim.value, fc='lightgray', alpha=0.3, lw=0)
        xname = 'spherical radius'
    elif xaxis == 'vrms_k':
        axes.axvspan(1.e-3, dvlim, fc='lightgray', alpha=0.3, lw=0)
        xname = 'rms linewidth'
    elif xaxis == 'mlumco':
        xname = 'luminous mass'
    elif xaxis == 'mvir':
        xname = 'virial mass'
    elif xaxis == 'mlte':
        xname = 'LTE mass'
    elif xaxis == 'siglum':
        xname = r'$\Sigma$, CO-based'
    elif xaxis == 'siglte':
        xname = r'$\Sigma$, LTE-based'
    else:
        xname = xaxis
    if yaxis == 'rad_pc':
        axes.axhspan(1.e-3, radlim.value, fc='lightgray', alpha=0.3, lw=0)
        yname = 'spherical radius'
    if yaxis == 'vrms_k':
        axes.axhspan(1.e-3, dvlim, fc='lightgray', alpha=0.3, lw=0)
        yname = 'rms linewidth'
    elif yaxis == 'mlumco':
        yname = 'luminous mass'
    elif yaxis == 'mvir':
        yname = 'virial mass'
    # Plot lines of constant external pressure when sigvir is on y-axis
    elif yaxis == 'sigvir':
        xmod = np.logspace(xlims[0],xlims[1],100)
        ymod = xmod + (20/(3*np.pi*21.1))*1.e4/xmod
        axes.plot(xmod, ymod, linestyle='-', color='g')
        axes.text(10**-0.6, 10**3.25, '$P_{ext}$ = $10^4$ cm$^{-3}$ K', color='g', rotation=-45)
        ymod2 = xmod + (20/(3*np.pi*21.1))*1.e2/xmod
        axes.plot(xmod, ymod2, linestyle='-', color='m')
        axes.text(10**-0.95, 10**1.6, '$P_{ext}$ = $10^2$ cm$^{-3}$ K', color='m', rotation=-45)
        yname = r'$\Sigma$, virial'
    else:
        yname = yaxis
    # If axes have identical units then plot y=x line
    if cat[xaxis].unit == cat[yaxis].unit:
        xmod = np.logspace(xlims[0],xlims[1],20)
        axes.plot(xmod, xmod, linestyle='-', color='k')
    # Set plot limits and labels
    axes.set_xlim(10**xlims[0], 10**xlims[1])
    axes.set_ylim(10**ylims[0], 10**ylims[1])
    axes.set_xlabel(xname+' ['+str(cat[xaxis].unit)+']')
    axes.set_ylabel(yname+' ['+str(cat[yaxis].unit)+']')
    return

# -------------------------------------------------------------------------------

def pltprops(label, fghz=230.538, distpc=5.e4, dvkms=0.2, beam=2,
            xplot=['radius', 'v_rms', 'area_exact'],
            xlims=[[-1.5,1], [-2,2], [-1,3]],
            yplot=['v_rms' , 'flux', 'flux'],
            ylims=[[-2,1.5], [-1.5,4.5], [-2,4]],
            pltname=['rdv', 'dvflux', 'areaflux']):

    global deltav, avgbeam, dist, freq, axes
    deltav  = dvkms * u.km / u.s
    avgbeam = beam * u.arcsec
    dist    = distpc * u.pc
    freq    = fghz * u.GHz

    # checks/creates directory to place plots
    if os.path.exists('plots') == 0:
        os.makedirs('plots')

    label = 'plots/' + label


    params = {'text.usetex': False, 'mathtext.fontset': 'stixsans'}
#    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)

    cat = Table.read(label+'_full_catalog.txt', format='ascii.ecsv')
    if os.path.isfile(label+'_physprop_add.txt'):
        pcat = Table.read(label+'_physprop_add.txt', format='ascii.ecsv')
    else:
        pcat = Table.read(label+'_physprop.txt', format='ascii.ecsv')
    newcol = Column(pcat['area_pc2']*0., name='e_area_pc2')
    newcol.unit = 'pc2'
    pcat.add_column(newcol)

    # Get the indices of trunks, branches, leaves, and clusters.
    # idc[0] is a list of trunk indices
    # idc[1] is a list of branch indices
    # idc[2] is a list of leaf indices
    # idc[3] is a list of cluster indices
    idc=[0,0,0,0]
    for i, typ in enumerate(['trunks', 'branches', 'leaves', 'clusters']):
        with open(label+'_'+typ+'.txt', 'r') as f:
            reader=csv.reader(f, delimiter=' ')
            a = zip(*reader)
        idc[i] = map(int,a[0])

    # Get the lists of trunk descendants
    f=open(label+'_trunks.txt','r')
    text=f.read()
    trd = []
    for line in text.splitlines():
        trd.append(map(int, line.split('|')[1].split(',')))    

    # Get the lists of cluster descendants and colors
    f=open(label+'_clusters.txt','r')
    text=f.read()
    cld = []
    clco = []
    for line in text.splitlines():
        cld.append(map(int, line.split('|')[1].split(',')))
        clco.append(line.split()[1]) 

    # Histogram of PAs
    val = 'position_angle'
    types = ['trunks', 'branches', 'leaves']
    bin_size = 15; min_edge = 0; max_edge = 180
    N = (max_edge-min_edge)/bin_size
    bin_list = np.linspace(min_edge, max_edge, N+1)
    pltdata = []
    for i in range(len(types)):
        xdata = cat[val][idc[i]]
        xdata[xdata<0] += 180.
        pltdata.append(xdata)
    fig, axes = plt.subplots()
    axes.hist(pltdata, bin_list, normed=0, histtype='bar', label=types)
    majorLocator = MultipleLocator(bin_size*2)
    minorLocator = MultipleLocator(bin_size)
    axes.xaxis.set_major_locator(majorLocator)
    axes.xaxis.set_minor_locator(minorLocator)
    #majorFormatter = FormatStrFormatter('%d')
    #ax.xaxis.set_major_formatter(majorFormatter)
    axes.set_xlabel(val+' ['+str(cat[val].unit)+']')
    axes.set_ylabel('Number')
    plt.legend(loc='best',fontsize='medium')
    plt.savefig(label+'_pahist.pdf', bbox_inches='tight')
    plt.close()

    # Size-linewidth relation, color coded
    plotx = 'rad_pc'
    ploty = 'vrms_k'
    eplotx = 'e_'+plotx
    eploty = 'e_'+ploty
    z    = ['x_cen', 'y_cen', 'v_cen', 'tpkav']
    cmap = plt.cm.get_cmap('nipy_spectral')
    for i in range(len(z)):
        fig, axes = plt.subplots()
        plt.errorbar(pcat[plotx], pcat[ploty], xerr=pcat[eplotx]*pcat[plotx], 
            yerr=pcat[eploty]*pcat[ploty], ecolor='dimgray', capsize=0, 
            zorder=1, marker=None, ls='None', lw=1, label=None)
        sctplot( pcat[plotx], pcat[ploty], cat[z[i]], 
            mec='none', msize=30, zorder=2, cmap=cmap, label=z[i] )
        std_overlay(pcat, plotx, ploty, xlims[0], ylims[0])
        shortname = re.sub('_', '', z[i])
        plt.savefig(label+'_rdv_'+shortname+'.pdf', bbox_inches='tight')
        plt.close()

    # Plot trunks, branches, leaves
    for i in range(len(xplot)):
        eplotx = 'e_'+xplot[i]
        eploty = 'e_'+yplot[i]
        fig, axes = plt.subplots()
        plt.errorbar(pcat[xplot[i]], pcat[yplot[i]], xerr=pcat[eplotx]*pcat[xplot[i]], 
            yerr=pcat[eploty]*pcat[yplot[i]], ecolor='dimgray', capsize=0, 
            zorder=1, marker=None, ls='None', lw=1, label=None)
        sctplot ( pcat[xplot[i]][idc[0]], pcat[yplot[i]][idc[0]], col='brown',
            mark='p', mec='k', msize=100, zorder=4, label='trunks' )
        sctplot ( pcat[xplot[i]][idc[1]], pcat[yplot[i]][idc[1]], col='w',
            mark='v', mec='k', msize=17, zorder=2, label='branches' )
        sctplot ( pcat[xplot[i]][idc[2]], pcat[yplot[i]][idc[2]], col='green',
            mark='o', mec='k', msize=20, zorder=3, label='leaves' )
        sctplot ( pcat[xplot[i]], pcat[yplot[i]], col='w',
            mark=None, mec='w', msize=1, zorder=0, label=None, linfit='b' )
        std_overlay(pcat, xplot[i], yplot[i], xlims[i], ylims[i])
        plt.legend(loc='lower right',fontsize='small',scatterpoints=1)
        plt.savefig(label+'_'+pltname[i]+'_full.pdf', bbox_inches='tight')
        plt.close()

    # Plot trunks and their descendants
    for i in range(len(xplot)):
        eplotx = 'e_'+xplot[i]
        eploty = 'e_'+yplot[i]
        fig, axes = plt.subplots()
        cmap = plt.cm.get_cmap('jet')
        plt.errorbar(pcat[xplot[i]][idc[0]], pcat[yplot[i]][idc[0]], 
            xerr=pcat[eplotx][idc[0]]*pcat[xplot[i]][idc[0]], 
            yerr=pcat[eploty][idc[0]]*pcat[yplot[i]][idc[0]], ecolor='dimgray', 
            capsize=0, zorder=2, marker=None, ls='None', lw=1, label=None)
        sctplot ( pcat[xplot[i]][idc[0]], pcat[yplot[i]][idc[0]], mark='s', 
            zorder=4, cmap=cmap, msize=30 )
        colors = plt.cm.jet(np.linspace(0, 1, len(idc[0])))
        for j, tno in enumerate(idc[0]):
            plt.errorbar(pcat[xplot[i]][trd[j]], pcat[yplot[i]][trd[j]], 
                xerr=pcat[eplotx][trd[j]]*pcat[xplot[i]][trd[j]], 
                yerr=pcat[eploty][trd[j]]*pcat[yplot[i]][trd[j]], ecolor='dimgray', 
                capsize=0, zorder=1, marker=None, ls='None', lw=1, label=None)
            sctplot ( pcat[xplot[i]][trd[j]], pcat[yplot[i]][trd[j]], col='w', 
                mec=colors[j], zorder=3, msize=10, label='trunk'+str(tno) )
        std_overlay(pcat, xplot[i], yplot[i], xlims[i], ylims[i])
        plt.legend(loc='lower right',fontsize='small',scatterpoints=1)
        plt.savefig(label+'_'+pltname[i]+'_trunks.pdf', bbox_inches='tight')
        plt.close()

    # Plot clusters and their descendants (get marker color from table)
    for i in range(len(xplot)):
        eplotx = 'e_'+xplot[i]
        eploty = 'e_'+yplot[i]
        fig, axes = plt.subplots()
        plt.errorbar(pcat[xplot[i]][idc[3]], pcat[yplot[i]][idc[3]], 
            xerr=pcat[eplotx][idc[3]]*pcat[xplot[i]][idc[3]], 
            yerr=pcat[eploty][idc[3]]*pcat[yplot[i]][idc[3]], ecolor='dimgray', 
            capsize=0, zorder=2, marker=None, ls='None', lw=1, label=None)
        sctplot ( pcat[xplot[i]][idc[3]], pcat[yplot[i]][idc[3]], mark='s', 
            zorder=4, col=clco, msize=25, linfit='b' )
        for j, tno in enumerate(idc[3]):
            plt.errorbar(pcat[xplot[i]][cld[j]], pcat[yplot[i]][cld[j]], 
                xerr=pcat[eplotx][cld[j]]*pcat[xplot[i]][cld[j]], 
                yerr=pcat[eploty][cld[j]]*pcat[yplot[i]][cld[j]], ecolor='dimgray', 
                capsize=0, zorder=1, marker=None, ls='None', lw=1, label=None)
            sctplot ( pcat[xplot[i]][cld[j]], pcat[yplot[i]][cld[j]], col='w', 
                mec=clco[j], zorder=3, msize=10, label='cluster'+str(tno) )
        std_overlay(pcat, xplot[i], yplot[i], xlims[i], ylims[i])
#        plt.legend(loc='best',fontsize='xx-small',numpoints=1)
        plt.savefig(label+'_'+pltname[i]+'_clusters.pdf', bbox_inches='tight')
        plt.close()

    return
