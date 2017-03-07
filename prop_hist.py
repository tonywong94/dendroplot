#!/usr/bin/env python2.7

# WIP script based on mass_spec.py by Tony Wong
# generalized histogram maker of any column from xxx_physprop_add.txt
# currently defined functions in mass_spec are mspec_overlay and mass_spec
# mspec_overlay is not used much, more of a variable translation table

import astropy.constants as const
import astropy.units as u
import csv
import matplotlib.pyplot as plt
import numpy as np
import re
import os
from astropy.table import Table, Column
from matplotlib.ticker import ScalarFormatter, FixedLocator, MultipleLocator
from scipy import stats


# ----------------------------------------------------------------------------------------
def prop_hist(label, dolog = True, val = 'mvir', lims = [], bin_size = 0, types = []):
    # new version of mass_spec function
    # label is the name of the data set (eg 30dor, pcc, etc.)
    # val is the input from the label_physprop_add.txt columns

    # global deltav, avgbeam, dist, freq, axes
    # deltav  = dvkms * u.km / u.s
    # avgbeam = beam * u.arcsec
    # dist    = distpc * u.pc
    # freq    = fghz * u.GHz

    # makes directory in which plots will be created if it doesn't already exist
    if os.path.isdir('plots') == 0:
        os.mkdir('plots')

    # default types is all of them, can be specified to only plot a few
    # loop below catches types not consistent with parameters, will not continue for illegal types
    if len(types) == 0:
        types = ['trunks', 'branches', 'leaves', 'clusters']
    else:
        for t in types:
            if t == 'trunks':
                continue
            elif t == 'branches':
                continue
            elif t == 'leaves':
                continue
            elif t == 'clusters':
                continue
            else:
                print('Type \'{}\' not recognized from types [trunks, branches, leaves, clusters], exiting...'.format(t))
                return

    # formatting parameters
    params = {'text.usetex' : False, 'mathtext.fontset' : 'stixsans'}
        # params = {'mathtext.default' : 'regular'}
    plt.rcParams.update(params)

    # read in data
    if os.path.isfile(label + '_physprop_add.txt') == True:
        pcat = Table.read(label + '_physprop_add.txt', format = 'ascii.ecsv')
    else:
        pcat = Table.read(label + '_physprop_add.txt', format = 'ascii.ecsv')


    # get indicies of trunks, branches, leaves, and clusters
    # defauls are as follows: idc[0] is a lis of trunk indices,
    # idc[1] of branches, idc[2] of leaves, idc[3] of clusters
    idc = [0, 0, 0, 0]
    for i, typ in enumerate(types):
        with open(label + '_' + typ + '.txt', 'r') as f:
            reader = csv.reader(f, delimiter = ' ')
            a = zip(*reader)
        idc[i] = map(int, a[0])

    # histogram of masses
    # val = 'mvir' as default
    # original default values: bin_size = 0.25, min_edge = 0, max_edge = 4

    pltdata = []
    for i in range(len(types)):
        data  = pcat[val][idc[i]]
        xdata = np.log10(data[data > 0])
        pltdata.append(xdata)

    # new automation loops
    # should auto generate limits, round to reasonable bin size
    if len(lims) == 0:
        limvals  = [item for sublist in pltdata for item in sublist]
        limmin   = np.nanmin(limvals)
        limmax   = np.nanmax(limvals)
        min_edge = float('%.3g' % limmin)
        max_edge = float('%.3g' % limmax)
        
        if bin_size == 0:
            limdif   = limmax - limmin
            bin_size = float('%.2g' % (limdif/20))

    elif len(lims) == 1:
        print('Only one limit ({}) specified, exiting...'.format(lims[0]))
        return

    else:
        min_edge = lims[0]
        max_edge = lims[1]

        if bin_size == 0:
            bin_size = (lims[1] - lims[0]) / 20

    N = (max_edge - min_edge) / bin_size
    bin_list = np.linspace(min_edge, max_edge, N + 1)   


    fig, axes = plt.subplots()
    n, bins, patches = axes.hist(pltdata, bin_list, normed = 0, log = dolog, histtype = 'bar', label = types, rwidth = .8)
    #majorLocator = MultipleLocator(bin_size * 2)
    #minorLocator = MultipleLocator(bin_size)
    #axes.xaxis.set_major_locator(majorLocator)
    #axes.xaxis.set_minor_locator(minorLocator)

    axes.xaxis.set_minor_locator(FixedLocator(bin_list[1::2] + bin_size/2))
    axes.xaxis.set_major_locator(FixedLocator(bin_list[::2] + bin_size/2))

    axes.tick_params(labelsize = 6)
    axes.set_xlabel('log ' + val + ' [' + str(pcat[val].unit) + ']')
    axes.set_ylabel('Number')
    
    axes.set_yscale('symlog')
    axes.yaxis.set_major_formatter(ScalarFormatter())
    
    #! some outputs have some blank space on lhs, still investigating
    plt.title('{0}_{1}'.format(label, val))
    plt.legend(loc = 'best', fontsize = 'medium')
    #plt.show()
    plt.savefig('plots/' + label + '_' + val + '_hist.pdf', bbox_inches = 'tight')
    plt.close()

    print('Plot created successfully for {0}_{1}'.format(label, val))

    return
