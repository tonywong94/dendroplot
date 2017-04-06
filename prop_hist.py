#!/usr/bin/env python2.7

# WIP script based on mass_spec.py by Tony Wong
# generalized histogram maker of any column from xxx_physprop_add.txt
# currently defined functions in mass_spec are mspec_overlay and mass_spec
# mspec_overlay is not used much, more of a variable translation table

import astropy.constants as const
import astropy.units as u
    # not sure if above imports will be necessary
import csv
import matplotlib.pyplot as plt
import numpy as np
import re
import os
from astropy.table import Table, Column
from matplotlib.ticker import ScalarFormatter, FixedLocator, MultipleLocator
from scipy import stats


# ----------------------------------------------------------------------------------------
def prop_hist(label = '', dolog = True, val = 'mvir', lims = [], binnum = 20, binsize = 0, types = []):
    # new version of mass_spec function
    # label is the name of the data set (eg 30dor, pcc, etc.)
    # val is the input from the label_physprop_add.txt columns

    # global deltav, avgbeam, dist, freq, axes
    # deltav  = dvkms * u.km / u.s
    # avgbeam = beam * u.arcsec
    # dist    = distpc * u.pc
    # freq    = fghz * u.GHz

    # sanity check on label
    if label == '':
        label = str(input('No label [dataset_line] given, enter name to continue: '))

    # default types is all of them, can be specified to only plot a few
    # loop below catches types not consistent with parameters
    if len(types) == 0:
        types = ['trunks', 'branches', 'leaves', 'clusters']
    else:
        for t in types:
            if t not in ['trunks', 'branches', 'leaves', 'clusters']:
                print('Type \'{}\' not recognized from default types, exiting...'.format(t))
                return
            else:
                continue

    # makes directory in which plots will be created if it doesn't already exist
    if os.path.isdir('../prophists') == 0:
        os.mkdir('../prophists')

    # formatting parameters
    params = {'text.usetex' : False, 'mathtext.fontset' : 'stixsans'}
        # params = {'mathtext.default' : 'regular'}
    plt.rcParams.update(params)

    # read in data
    if os.path.isfile('../props/' + label + '_physprop_add.txt'):
        pcat = Table.read('../props/' + label + '_physprop_add.txt', format = 'ascii.ecsv')
    else:
        pcat = Table.read('../props/' + label + '_physprop.txt', format = 'ascii.ecsv')

    # get indicies of trunks, branches, leaves, and clusters
    # defauls are as follows: idc[0] is a lis of trunk indices,
    # idc[1] of branches, idc[2] of leaves, idc[3] of clusters
    idc = [0, 0, 0, 0]
    for i, typ in enumerate(types):
        with open('../props/' + label + '_' + typ + '.txt', 'r') as f:
            reader = csv.reader(f, delimiter = ' ')
            a = zip(*reader)
        idc[i] = map(int, a[0])

    # data flattening
    pltdata = []
    for i in range(len(types)):
        data  = pcat[val][idc[i]]
        xdata = np.log10(data[data > 0])
        pltdata.append(xdata)

    # limit and bin size determination
    if len(lims) == 0:
        limvals = [item for sublist in pltdata for item in sublist]
        limmin  = np.around(np.nanmin(limvals), 2)
        limmax  = np.around(np.nanmax(limvals), 2)
    elif len(lims) == 1:
        print('Only one limit ({}) specified, exiting...'.format(lims[0]))
        return
    else:
        if len(lims) > 2:
            print('Only first two limits will be used')
        limmin = lims[0]
        limmax = lims[1]

    if binsize == 0:
        limdif  = limmax - limmin
        optsize = np.around(limdif / binnum, 3)

        # choosing logic, preset sizes: .01, .02, .04, .08, .1, .2, .4, .8, 1
        # arbitrary but look good for the current datasets
        if (optsize < .01):
            binsize = .01
        elif (.01 <= optsize < .025):
            binsize = .02
        elif (.025 <= optsize < .06):
            binsize = .04
        elif (.06 <= optsize < .09):
            binsize = .08
        elif (.09 <= optsize < .15):
            binsize = .1
        elif (.15 <= optsize < .25):
            binsize = .2
        elif (.25 <= optsize < .6):
            binsize = .4
        elif (.6 <= optsize < .9):
            binsize = .8
        elif (.9 <= optsize):
            binsize = 1

    # bin spacing, gives enough room on either end of plot
    binlist = np.arange((limmin - 2*binsize), (limmax + 2*binsize), binsize)

    # plotting section
    fig, axes = plt.subplots()
    n, bins, patches = axes.hist(pltdata, binlist, normed = 0, log = dolog, histtype = 'bar', label = types, rwidth = .8)

    # changes y-axis to linear for small distributions (looks much cleaner)
    nummax = max([np.max(o) for o in n])
    if (nummax < 14):
        plt.cla()
        n, bins, patches = axes.hist(pltdata, binlist, normed = 0, log = 0, histtype = 'bar', label = types, rwidth = .8)

    axes.xaxis.set_minor_locator(FixedLocator(binlist[1::2] + binsize/2))
    axes.xaxis.set_major_locator(FixedLocator(binlist[0::2] + binsize/2))

    axes.tick_params(labelsize = 6)
    axes.set_xlabel('log ' + val + ' [' + str(pcat[val].unit) + ']')
    axes.set_ylabel('Number of objects binned')
    axes.yaxis.set_major_formatter(ScalarFormatter())

    plt.title('{0}_{1}'.format(label, val))
    plt.legend(loc = 'best', fontsize = 'medium')
    #plt.show()
    plt.savefig('../prophists/' + label + '_' + val + '_hist.pdf', bbox_inches = 'tight')
    plt.close()

    print('Plot created successfully for {0}_{1}'.format(label, val))
    
    return
