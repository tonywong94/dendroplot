#!/usr/bin/env python3

# Written for Python 3.5.2
# Originally authored by Evan Wojciechowski, based on scripts by Chaoyue Cui and Tony Wong
# Function file that is built to be called by script files; will take data from fits files and plot histograms of the data
# Derived from similar script histomom.py

from astropy.io import fits
import math
from matplotlib.colors import LogNorm
import matplotlib.figure
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import numpy as np
import os
from itertools import cycle


def histplot(xname=None, yname=None, snrcut=0, dolog2d=False, dolog1d=False, 
        xlbl=None, ylbl=None, shift=None, shiftlbl=None, 
        nbins=100, outname = '', extrema = [], cd = ''):
    # xname and yname are the two necessary inputs

    if cd != '':
        if os.path.exists(cd) == 1:
            print('Found {}, changing directory...'.format(cd))
            os.chdir(cd)
        else:
            print('Directory {} doesn\'t exist, creating and changing...\n'.format(cd))
            os.mkdir(cd)
            os.chdir(cd)

    ####### Data aquisition and error modification #######

    # Checks if file names are provided
    if xname == None:
        xname = input('No file specified, enter name of file for x-axis data: ')
    if yname == None:
        yname = input('No file specified, enter name of file for y-axis data: ')

    xfile = fits.open(xname)
    yfile = fits.open(yname)
    print('\nOpening {0} as xfile and {1} as yfile\n'.format(xname,yname))

    # Data extraction
    xread = xfile[0].data.flatten()
    yread = yfile[0].data.flatten()
    # Reject NaN values
    good = np.intersect1d(np.where(~np.isnan(xread))[0], np.where(~np.isnan(yread))[0])
    xdata = xread[good]
    ydata = yread[good]

    print('Raw xdata ranges from {0} to {1}'.format(min(xdata),max(xdata)))
    print('Raw ydata ranges from {0} to {1}\n'.format(min(ydata),max(ydata)))

    # When snrcut > 0, error data is required, and file names are inputted next
    if snrcut > 0:
        # Redefine snrcut as a float variable for ease later
        scut = float(snrcut)

        # Automatically assigns error data file names following a pattern
        xerrname = xname.replace('mom','emom')
        yerrname = yname.replace('mom','emom')

        # Checks and if necessary corrects error file names
        if os.path.exists(xerrname) == False:
            xerrname = input('Enter name of file for x-axis error data: ')
        if os.path.exists(yerrname) == False:
            yerrname = input('Enter name of file for y-axis error data: ')

        xerrfile = fits.open(xerrname)
        yerrfile = fits.open(yerrname)
        print('Opening {0} as xerrfile and {1} as yerrfile\n'.format(xerrname,yerrname))

        # Error data extraction
        xerrdata = xerrfile[0].data.flatten()[good]
        yerrdata = yerrfile[0].data.flatten()[good]
        print('Raw xerrdata ranges from {0} to {1}'.format(min(xerrdata),max(xerrdata)))
        print('Raw yerrdata ranges from {0} to {1}\n'.format(min(yerrdata),max(yerrdata)))

        # Sets up arrays for detections and non-detections in each axis and the four posibilites for each point
        xydet = np.intersect1d(np.where(xdata >= scut*xerrdata)[0],
            np.where(ydata >= scut*yerrdata)[0])
        xndet = np.intersect1d(np.where(xdata < scut*xerrdata)[0],
            np.where(ydata >= scut*yerrdata)[0])
        yndet = np.intersect1d(np.where(xdata >= scut*xerrdata)[0],
            np.where(ydata < scut*yerrdata)[0])
        xyndet = np.intersect1d(np.where(xdata < scut*xerrdata)[0],
            np.where(ydata < scut*yerrdata)[0])

        # Replace non-detections with upper limit values
        for i in np.concatenate([xndet,xyndet]):
            xdata[i] = scut*xerrdata[i]
        for i in np.concatenate([yndet,xyndet]):
            ydata[i] = scut*yerrdata[i]

        print('Number of detections: {0}'.format(len(xydet)))
        print('Non-detections only in xdata: {0}'.format(len(xndet)))
        print('Non-detections only in ydata: {0}'.format(len(yndet)))
        print('Non-detections in both xdata and ydata: {0}\n'.format(len(xyndet)))

        print('With upper limits xdata ranges from {0} to {1}'.format(min(xdata),max(xdata)))
        print('With upper limits ydata ranges from {0} to {1}'.format(min(ydata),max(ydata)))

    # Modify raw/error-corrected data into log scale if desired
    if dolog2d == True:
        pos = np.intersect1d(np.where(xdata>0)[0], np.where(ydata>0)[0])
        x = np.log10(xdata[pos])
        y = np.log10(ydata[pos])
    else:
        x = xdata
        y = ydata

    # Max and min boundaries for the graphs

    if extrema == []:
        xmin = math.floor(min(x))
        ymin = math.floor(min(y))
        xmax = math.ceil(max(x))
        ymax = math.ceil(max(y))
    else:
        xmin = extrema[0]
        ymin = extrema[1]
        xmax = extrema[2]
        ymax = extrema[3]

#     finitex = np.isfinite(x)
#     finitey = np.isfinite(y)
# 
#     fordeletion = []
#     for i in range(len(x)):
#         if (finitex[i] == 0 or finitey[i] == 0):
#             fordeletion.append(i)
# 
#     x = np.delete(x, fordeletion)
#     y = np.delete(y, fordeletion)

    ####### Set up geometry of three plots #######

    # Values pulled directly from histomom.py
    left, width = 0.12, 0.55
    bottom, height = 0.12, 0.55
    bottom_h = left_h = left+width+0.02

    img = [left, bottom, width, height] # Dimensions of temperature plot
    rect_histx = [left, bottom_h, width, 0.2] # Dimensions of x-histogram
    rect_histy = [left_h, bottom, 0.2, height] # Dimensions of y-histogram

    fig = plt.figure(1, figsize=(9.5,9))

    ax = plt.axes(img) # Temperature plot
    axHistx = plt.axes(rect_histx) # x-histogram
    axHisty = plt.axes(rect_histy) # y-histogram

    # Remove the inner axes numbers of this histograms to save space
    nullfmt = tck.NullFormatter()
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    ####### Plot 2-D histogram #######

    xbins = np.linspace(xmin, xmax, nbins)
    ybins = np.linspace(ymin, ymax, nbins)

    # Tuple[s] setup for data plotting
    ### Transposing counts after definition seems to be important to data processing, not yet understood
    if snrcut > 0:
        # 0 -> xydet, 1 -> xndet, 2 -> yndet, 3 -> xyndet
        counts0, xedge0, yedge0 = np.histogram2d(x[xydet], y[xydet], bins=(xbins,ybins))
        counts1, xedge1, yedge1 = np.histogram2d(x[xndet], y[xndet], bins=(xbins,ybins))
        counts2, xedge2, yedge2 = np.histogram2d(x[yndet], y[yndet], bins=(xbins,ybins))
        counts3, xedge3, yedge3 = np.histogram2d(x[xyndet], y[xyndet], bins=(xbins,ybins))

        counts0 = np.transpose(counts0)
        counts1 = np.transpose(counts1)
        counts2 = np.transpose(counts2)
        counts3 = np.transpose(counts3)

    else:
        # Also using 0 to streamline later plotting commands
        counts0, xedge0, yedge0 = np.histogram2d(x, y, bins=(xbins,ybins))

        counts0 = np.transpose(counts0)

    # Mesh edges to create coupled matricies
    xmesh, ymesh = np.meshgrid(xedge0, yedge0)

    ### Contour levels, need to figure out how to customize
        ### Change to np.logspace, add new input variable
        ### May need to compare the three ndet; whicever one has the largest amount of values, set that
    levels = [25, 100, 400, 1600, 6400, 25600]

    # Set limits on axes
    ax.set_xlim(xedge0[0], xedge0[-1])
    ax.set_ylim(yedge0[0], yedge0[-1])

    # Colorizes detections
    img = ax.pcolormesh(xmesh, ymesh, counts0, norm=LogNorm(), cmap=plt.cm.gist_ncar_r)
        ### If logscale on colorbar is not desired, change norm=None

    # Keeps contour boundaries in range of x and y edges
    extent = [xedge0[0], xedge0[-1], yedge0[0], yedge0[-1]]
    
    # Build contours
    #cset0 = ax.contour(counts0, levels, colors='green', extent=extent)
    if snrcut > 0:
        cset1 = ax.contour(counts1, levels, colors='black', extent=extent)
        cset2 = ax.contour(counts2, levels, colors='blue', extent=extent)
        cset3 = ax.contour(counts3, levels, colors='red', extent=extent)

        ### Section of histomom.py that is commented out
        ### Unknown reason why
        # cset1.collections[0].set_label('x non-det')
        # cset2.collections[0].set_label('y non-det')
        # cset3.collections[0].set_label('x and y non-det')
        # ax.clabel(cset1, inline=True, colors='k', fontsize=8)
        # legend = plt.legend(loc='upper left', fontsize='medium')
        ### End of commented section

    ####### Plot 2-D labels/ticks #######
    
    if xlbl is None:
        xlbl = xname
    if ylbl is None:
        ylbl = yname
    # Name labels
    if dolog2d == True:
        xlabel = 'log ({0})'.format(xlbl)
        ylabel = 'log ({0})'.format(ylbl)
    else:
        xlabel = xlbl
        ylabel = ylbl

    ax.set_xlabel(xlabel, fontsize=14, labelpad=10)
    ax.set_ylabel(ylabel, fontsize=14)

    # Format ticks
    ticklabels = ax.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(10)
        label.set_family('serif')
    
    ticklabels = ax.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(10)
        label.set_family('serif')

    ####### Plot 1-D histogram #######

    ### May need to add another log scale somewhere

    # Set limits
    axHistx.set_xlim(xmin, xmax)
    axHisty.set_ylim(ymin, ymax)
    
    # Set bins
    nxbins = np.arange(xmin, xmax, (xmax-xmin)/nbins)
    nybins = np.arange(ymin, ymax, (ymax-ymin)/nbins)

    # Plot histograms
    if dolog1d == True:
        if snrcut > 0:
            axHistx.hist(x[xydet], bins=nxbins, color='blue')
            axHisty.hist(y[xydet], bins=nybins, orientation='horizontal', color='red')
            axHistx.set_yscale('log')
            axHisty.set_xscale('log')
        else:
            axHistx.hist(x, bins=nxbins, color='blue')
            axHisty.hist(y, bins=nybins, orientation='horizontal', color='red')
            axHistx.set_yscale('log')
            axHisty.set_xscale('log')
    else:
        if snrcut > 0:
            axHistx.hist(x[xydet], bins=nxbins, color='blue')
            axHisty.hist(y[xydet], bins=nybins, orientation='horizontal', color='red')
        else:
            axHistx.hist(x, bins=nxbins, color='blue')
            axHisty.hist(y, bins=nybins, orientation='horizontal', color='red')

    ####### Plot 1-D ticks #######
    
    # Format ticks
    ticklabels = axHistx.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(10)
        label.set_family('serif')

    ticklabels = axHisty.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(10)
        label.set_family('serif')

    # Set location
    axHistx.yaxis.set_major_locator(tck.MaxNLocator(4))
    axHisty.xaxis.set_major_locator(tck.MaxNLocator(4))

    ####### Miscellaneous #######

    ### Arbitrary line plotted for limit I think?
    z = np.linspace(-10, 300)
    ax.plot(z, z)
    
    if shift is not None:
        ax.plot(z, z-shift, color='red', linewidth=4, alpha=0.5)
        ax.text((xmax-0.0*(xmax-xmin)), (ymax-shift-0.08*(ymax-ymin)), shiftlbl, 
            horizontalalignment='right', color='r', rotation=45, size=14)

    ####### Output file name #######

    if outname == '':
        figname = input('Enter output file name: ')
    else:
        figname = outname
    figname += '.pdf'
    plt.savefig(figname, bbox_inches='tight')
    plt.close()

### Output print statements as a file, .txt
### Output detection arrays?

def cumhist1d(xplot=[], nbins=100, xminlg=-2, xmaxlg=3, outname = 'cumhist1d'):
    lines = ["-","--","-.",":"]
    linecycler = cycle(lines)
    fig, axes = plt.subplots()
    for i in range(len(xplot)):
        xfile = fits.getdata(xplot[i])
        fname = os.path.basename(xplot[i]).replace('.mom0.fits.gz','')
        print('\nOpening {0}'.format(xplot[i]))
        xread = xfile.flatten()
        xdata = xread[np.where(~np.isnan(xread))[0]]
        print('Max value is {0}'.format(np.max(xdata)))
        values, base = np.histogram(xdata, bins=nbins, normed=True)
        dx = base[1] - base[0]
        print('Bin width is {0}\n'.format(dx))
        cumulative = np.cumsum(values)*dx
        plt.plot(base[:-1], cumulative, next(linecycler), label=fname)
    axes.set_xscale("log")
    axes.set_xlim(10**xminlg, 10**xmaxlg)
    axes.set_ylim(0, 1.)
    axes.set_xlabel('ICO [K km/s]')
    axes.set_ylabel('Fraction of Pixels < ICO')
    plt.legend(loc='best',fontsize='x-small', handlelength=3)
    plt.savefig(outname+'.pdf', bbox_inches='tight')
    plt.close()

def cumflux1d(xplot=[], nbins=100, xminlg=-2, xmaxlg=3, outname = 'cumflux1d'):
    lines = ["-","--","-.",":"]
    linecycler = cycle(lines)
    fig, axes = plt.subplots()
    for i in range(len(xplot)):
        xfile = fits.getdata(xplot[i])
        fname = os.path.basename(xplot[i]).replace('.mom0.fits.gz','')
        print('\nOpening {0}'.format(xplot[i]))
        xread = xfile.flatten()
        xdata = xread[np.where(~np.isnan(xread))[0]]
        print('Max value is {0}\n'.format(np.max(xdata)))
        bins = np.logspace(xminlg, xmaxlg, num=nbins)
        cumflux=np.empty(nbins)
        for i in range(len(bins)):
            cumflux[i] = sum(xdata[xdata<bins[i]])
        cumflux /= np.max(cumflux)
        plt.plot(bins, cumflux, next(linecycler), label=fname)
    axes.set_ylim(0, 1.)
    axes.set_xscale("log")
    axes.set_xlabel('ICO [K km/s]')
    axes.set_ylabel('Fraction of Flux < ICO')
    plt.legend(loc='best',fontsize='x-small', handlelength=3)
    plt.savefig(outname+'.pdf', bbox_inches='tight')
    plt.close()
