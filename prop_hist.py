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
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy import stats


# ----------------------------------------------------------------------------------------
def prop_hist(label, dolog = True, val = 'mvir'):
    # new version of mass_spec function
    # label is the name of the data set (eg 30dor, pcc, etc.)
    # val is the input from the label_physprop_add.txt columns

    # global deltav, avgbeam, dist, freq, axes
    # deltav  = dvkms * u.km / u.s
    # avgbeam = beam * u.arcsec
    # dist    = distpc * u.pc
    # freq    = fghz * u.GHz

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
    # idc[0] is a lis of trunk indices,
    # idc[1] of branches, idc[2] of leaves, idc[3] of clusters
    idc = [0, 0, 0, 0]
    for i, typ in enumerate(['trunks', 'branches', 'leaves', 'clusters']):
        with open(label + '_' + typ + '.txt', 'r') as f:
            reader = csv.reader(f, delimiter = ' ')
            a = zip(*reader)
        idc[i] = map(int, a[0])

    # histogram of masses
    # val = 'mvir' as default
    types = ['trunks', 'branches', 'leaves', 'clusters']
    bin_size = 0.25
    min_edge = 0
    max_edge = 4
    N = (max_edge - min_edge) / bin_size
    bin_list = np.linspace(min_edge, max_edge, N + 1)
    #print(bin_list)
    pltdata = []
    for i in range(len(types)):
        data  = pcat[val][idc[i]]
        xdata = np.log10(data[data > 0])
        pltdata.append(xdata)

    fig, axes = plt.subplots()
    axes.hist(pltdata, bin_list, normed = 0, log = dolog, histtype = 'bar', label = types, rwidth = .8)
    majorLocator = MultipleLocator(bin_size * 2)
    minorLocator = MultipleLocator(bin_size)
    #axes.xaxis.set_major_locator(majorLocator)
    #axes.xaxis.set_minor_locator(minorLocator)
    plt.xticks(bin_list + bin_size/2)
    axes.tick_params(labelsize = 6)
    axes.set_xlabel('log ' + val + ' [' + str(pcat[val].unit) + ']')
    axes.set_ylabel('Number')

    plt.legend(loc = 'best', fontsize = 'medium')
    #plt.show()
    plt.savefig(label + '_' + val + 'hist.pdf', bbox_inches = 'tight')
    plt.close()

    return
