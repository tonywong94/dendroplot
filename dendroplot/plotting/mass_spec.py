#!/usr/bin/env python

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

# -------------------------------------------------------------------------------

def mspec_overlay(cat, xaxis, yaxis, xlims, ylims):
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

def mass_spec(label, val='mvir'):

#     global deltav, avgbeam, dist, freq, axes
#     deltav  = dvkms * u.km / u.s
#     avgbeam = beam * u.arcsec
#     dist    = distpc * u.pc
#     freq    = fghz * u.GHz

    params = {'text.usetex': False, 'mathtext.fontset': 'stixsans'}
#    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)

#    cat = Table.read(label+'_full_catalog.txt', format='ascii.ecsv')
    if os.path.isfile(label+'_physprop_add.txt'):
        pcat = Table.read(label+'_physprop_add.txt', format='ascii.ecsv')
    else:
        pcat = Table.read(label+'_physprop.txt', format='ascii.ecsv')
#    newcol = Column(pcat['area_pc2']*0., name='e_area_pc2')
#    newcol.unit = 'pc2'
#    pcat.add_column(newcol)

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

    # Histogram of Masses
#    val = 'mvir'
    types = ['trunks', 'branches', 'leaves', 'clusters']
    bin_size = 0.25; min_edge = 0; max_edge = 4
    N = (max_edge-min_edge)/bin_size
    bin_list = np.linspace(min_edge, max_edge, N+1)
    pltdata = []
    for i in range(len(types)):
        data  = pcat[val][idc[i]]
        xdata = np.log10(data[data>0])
        pltdata.append(xdata)
    fig, axes = plt.subplots()
    axes.hist(pltdata, bin_list, normed=0, log=True, histtype='bar', label=types)
    majorLocator = MultipleLocator(bin_size*2)
    minorLocator = MultipleLocator(bin_size)
    axes.xaxis.set_major_locator(majorLocator)
    axes.xaxis.set_minor_locator(minorLocator)
    #majorFormatter = FormatStrFormatter('%d')
    #ax.xaxis.set_major_formatter(majorFormatter)
    axes.set_xlabel('log '+val+' ['+str(pcat[val].unit)+']')
    axes.set_ylabel('Number')
    plt.legend(loc='best',fontsize='medium')
    plt.savefig(label+'_'+val+'hist.pdf', bbox_inches='tight')
    plt.close()

    return
