#!/usr/bin/env python2.7

# wip script to replace mass_spec and prop_hist scripts
# will include functions to plot individual or multiple datasets
# all 1-dimensional histograms

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

def hist1d_multiset(labels = ['30dor', 'pcc'], lines = ['12'], val = 'mvir', binnum = 20, binsize = 0, dolog = True, lims = [], types = []):
    # function will be able to run from any dataset's pyscript folder given that directory
    # structure is identical between datasets (i.e. parallel directories)

    # default types is all of them, can be specified to only plot a few
    # loop below catches types not consistent with parameters
    if len(types) == 0:
        types = ['trunks', 'branches', 'leaves', 'clusters']
    else:
        for t in types:
            if t not in ['trunks', 'branches', 'leaves', 'clusters']:
                print('Type \'{}\' not recognized from default types (trunks, branches, leaves, clusters), exiting...'.format(t))
                return
            else:
                continue

    # default output folder will be in parallel to running directory
    if os.path.isdir('../multiprop') == 0:
        os.mkdir('../multiprop')

    # formatting parameters
    params = {'text.usetex' : False, 'mathtext.fontset' : 'stixsans'}
    plt.rcParams.update(params)

    # dictionaries initialized for looping
    pcats    = {}   # dictionary of data from respepctive physprop_add files
    idcs     = {}   # dictionary of indicies of data for different types in physprop_add
    pltdatum = {}   # dictionary of data to plot

    # read in data
    for lb in range(len(labels)):
        for ln in range(len(lines)):
            fname = '../../' + labels[lb] + '/props/' + labels[lb] + '_' + lines[ln] + '_physprop_add.txt'
            lname = labels[lb] + '_' + lines[ln]

            if os.path.isfile(fname):
                pcats[lname] = Table.read(fname, format = 'ascii.ecsv')
            else:
                pcats[lname] = Table.read(fname.replace('prop_add', 'prop'), format = 'ascii.ecsv')

            # idc indicies are in the order of types list
            idc = [0,0,0,0]
            for i, typ in enumerate(types):
                with open(fname.replace('physprop_add', typ), 'r') as f:
                    reader = csv.reader(f, delimiter = ' ')
                    a = zip(*reader)
                idc[i] = map(int, a[0])
            idcs[lname] = idc

            # data flattening
            pltdata = []
            for i in range(len(types)):
                data  = pcats[lname][val][idcs[lname][i]]
                xdata = np.log10(data[data > 0])
                pltdata.append(xdata)
            pltdatum[lname] = pltdata

    # begin limit calculation
    # need to create list from all datasets that are being plotted
    if len(lims) == 0:
        limvals = []
        for p in pltdatum:
            #print(pltdatum[p])
            lim_swp = [item for sublist in pltdatum[p] for item in sublist]
            limvals.extend(lim_swp)
        limmin = np.around(np.nanmin(limvals), 2)
        limmax = np.around(np.nanmax(limvals), 2)
    elif len(lims) == 1:
        print('Only one limit ({}) specified, exiting...'.format(lims[0]))
        return
    else:
        if len(lims) > 2:
            print('Only first two limits will be used')
        limmin = lims[0]
        limmax = lims[1]

    if binsize == 0:
        limdif = limmax - limmin
        optsize = np.around(limdif / binnum, 3)

        # arbitrary choosing logic
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

    # bin spacing
    binlist = np.arange((limmin - 2*binsize), (limmax + 2*binsize), binsize)

    # plotting
    fig, axes = plt.subplots()
    for p in pltdatum:
        temp_types = [p + ' ' + t for t in types]
        axes.hist(pltdatum[p], binlist, normed = 0, log = dolog, histtype = 'bar', label = temp_types, rwidth = .8)
    axes.xaxis.set_minor_locator(FixedLocator(binlist[1::2] + binsize/2))
    axes.xaxis.set_major_locator(FixedLocator(binlist[::2] + binsize/2))

    axes.tick_params(labelsize = 6)
    axes.set_xlabel('log ' + val + ' [' + str(pcats[labels[0] + '_' + lines[0]][val].unit) + ']')
    axes.set_ylabel('Number of objects binned')

    axes.set_yscale('log')
    axes.yaxis.set_major_formatter(ScalarFormatter())
    plt.legend(loc = 'best', fontsize = 'medium')
    
    # automated looping for title, output name, output message
    title = '{} for '.format(val)
    for i in range(len(labels)):
        for j in range(len(lines)):
            title += '{}, '.format(labels[i] + '_' + lines[j])
    title = title[:len(title) - 2]
    plt.title(title)
    #plt.show()

    outname = '../multiprop/'
    msg = 'Plot created successfully for {} '.format(val)

    for i in range(len(labels)):
        outname += '{}_'.format(labels[i])
        msg += '{} '.format(labels[i])
    for j in range(len(lines)):
        outname += '{}CO_'.format(lines[j])
        msg += '{}CO '.format(lines[j])
    msg = msg[:len(msg)-1]

    outname += '{}.pdf'.format(val)
    plt.savefig(outname, bbox_inches = 'tight')
    plt.close()

    print(msg)

    return

def hist1d_oneset(label = '', dolog = True, val = 'mvir', lims = [], binnum = 20, binsize = 0, types = []):
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

    # default types is all of them, can be specified to only put a few
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
    plt.rcParams.update(params)

    # read in data
    if os.path.isfile('../props/' + label + '_physprop_add.txt'):
        pcat = Table.read('../props/' + label + '_physprop_add.txt', format = 'ascii.ecsv')
    else:
        pcat = Table.read('../props/' + label + '_physprop.txt', format = 'ascii.ecsv')

    # get indicies of types
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

        # choosing logic, presets are arbitrary but work well for 30dor and pcc
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

    # bin spacing, extra space for clarity
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

    plt.title('{0}_{1}'.format(label,val))
    plt.legend(loc = 'best', fontsize = 'medium')
    #plt.show()
    plt.savefig('../prophists/' + label + '_' + val + '_hist.pdf', bbox_inches = 'tight')
    plt.close()

    print('Plot created successfully for {0}_{1}'.format(label, val))

    return
