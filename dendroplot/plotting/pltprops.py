#!/usr/bin/env python

import numpy as np
from numpy.random import randint
from scipy import stats, odr
import os
from os.path import join
import re
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from astropy import units as u
from astropy import constants as const
from astropy.table import Table, Column
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from kapteyn import kmpfit

params = {'text.usetex': False, 'mathtext.fontset': 'stixsans'}
plt.rcParams.update(params)


def sctplot(xdata, ydata, xerr=None, yerr=None, zdata=None, cmap=None, 
            col='g', marker='o', mec='k', ms=6, linfitcolor=None, label=None,
            axes=None, **kwargs):
    '''
    General scatter plot between two variables

    === Parameters ===
    xdata : astropy.Column or numpy.array
        The values on the horizontal axis
    ydata : astropy.Column or numpy.array
        The values on the vertical axis
    xerr : astropy.Column or numpy.array
        Error bars for the horizontal axis
    yerr : astropy.Column or numpy.array
        Error bars for the vertical axis
    zdata : astropy.Column or numpy.array
        Optional values for color coding
    cmap : matplotlib.colors.Colormap
        Name of the color map for color coding
    col : matplotlib.colors
        Interior color(s) for plotted points.  Can be a single value or
        a list with the same length as xdata and ydata.
    marker : matplotlib.markers
        Marker type for plotted points
    mec : matplotlib.colors
        Edge color for plotted points.  Default is black.  Use 'none' to
        omit drawing the boundary.
    ms : float
        Marker size for plotted points
    linfitcolor : matplotlib.colors
        Color for regression fit.  Default is not to plot this.
    label : string
        Text label for colorbar or for legend.  Required.
    axes : matplotlib.axes
        Plot axes object
    **kwargs : dict
        Additional parameters for pyplot.scatter (e.g., zorder, alpha)

    === Returns ===
    cbar : matplotlib.pyplot.colorbar
        colorbar object if cmap and zdata are given
    '''
    if axes is None:
        axes = plt.gca()
    # Error bars at bottom layer if provided
    if yerr is not None:
        axes.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, 
                      ecolor='dimgray', capsize=0, zorder=0, marker=None, 
                      ls='None', elinewidth=1, label=None)
    # Single color plot
    if cmap is None and np.size(col) == 1:
        axes.scatter(xdata, ydata, marker=marker, c=col, edgecolors=mec, 
            s=ms, linewidths=1, label=label, **kwargs)
    # Multi-color plot, array of colors
    elif cmap is None:
        for xp, yp, cp in zip(xdata, ydata, col):
            axes.scatter(xp, yp, marker=marker, c=cp, edgecolors=mec, 
                s=ms, **kwargs)
    # Continuous color map based on structure number
    elif zdata is None:
        xind = np.linspace(0., 1., len(xdata))
        axes.scatter(xdata, ydata, marker=marker, c=xind,
            cmap=cmap, edgecolors=mec, s=ms, label=None, **kwargs)
    # Continuous color map based on property value
    else:
        sc=axes.scatter(xdata, ydata, marker=marker, c=zdata,
            cmap=cmap, edgecolors=mec, s=ms, label=None, **kwargs)
        cbar = plt.colorbar(sc)
        cbar.ax.tick_params(labelsize=9) 
        cbar.set_label(label, rotation=90)
    # Do a linear regression fit if requested
    if linfitcolor is not None:
        sorted=np.argsort(xdata)
        m, b, rval, pval, std_err = stats.linregress(xdata[sorted],
            ydata[sorted])
        xlims = axes.get_xlim()
        xmod = np.linspace(xlims[0],xlims[1],20)
        ymod = b + m * xmod
        axes.plot(xmod, ymod, linestyle='--', color=linfitcolor)
        axes.text(0.03,0.95,'slope = $%4.2f$' % m,size=9,transform=axes.transAxes)
        axes.text(0.03,0.90,'intercept = $%4.2f$' % b,size=9,transform=axes.transAxes)
    if cmap is not None and zdata is not None:
        return cbar
    else:
        return


def color_code_bin(x, y, z, xerr=None, yerr=None, lobin=25, lobin_col='cyan', 
            hibin=75, hibin_col='salmon', ms=6, xlims=None, nbin=8, binms=8,
            colname=None, alpha=None, cmap='jet', axes=None, zlog=False, median=False):
    '''
    Color coded scatter plot with binned values in upper and lower quartiles.

    === Parameters ===
    x : astropy.Column or numpy.array
        The values on the horizontal axis
    y : astropy.Column or numpy.array
        The values on the vertical axis
    z : astropy.Column
        The values to use for the plot symbol color
    xerr : astropy.Column or numpy.array
        Error bars for the horizontal axis.  Optional.
    yerr : astropy.Column or numpy.array
        Error bars for the vertical axis.  Optional.
    lobin : float
        Percentile for lower bin sequence.  Default plots bottom quartile.
    lobin_col : matplotlib.colors
        Symbol color for lower bin sequence.
    hibin : float
        Percentile for upper bin sequence.  Default plots top quartile.
    hibin_col : matplotlib.colors
        Symbol color for upper bin sequence.
    ms : float
        Marker size for plotted points; passed to scatter
    xlims : tuple of float
        The x limits of the plot.  Default is [x.min(), x.max()]
    nbin : int
        Number of bins to place uniformly across xlims
    binms : float
        Marker size for binned points; passed to errorbar
    colname : string
        Name of 'z' column, used for labeling colorbar only.
    alpha : float
        Transparency parameter for color-coded scatter plots, between 0 and 1
    cmap : matplotlib.colors.Colormap
        Name of the color map for color coding
    axes : matplotlib.axes
        Plot axes object
    zlog : boolean
        True if z-values are logarithmic, used for labeling colorbar only.
    median : boolean
        True to also plot binned medians as gray steps
    '''
    if axes is None:
        axes = plt.gca()
    if colname is None:
        colname = z.name
    if zlog:
        zlbl = 'log ' + colname +' ['+str(z.unit)+']'
    else:
        zlbl = colname +' ['+str(z.unit)+']'
    zlbl = zlbl.replace('[]','')
    cb = sctplot(x, y, xerr=xerr, yerr=yerr, zdata=z, mec='none', ms=ms, 
                 zorder=2, cmap=cmap, label=zlbl, alpha=alpha, axes=axes)
    q1 = np.nanpercentile(z, lobin)
    print('Percentile {} for {} has value {}'.format(lobin,colname,q1))
    q2 = np.nanpercentile(z, hibin)
    print('Percentile {} for {} has value {}'.format(hibin,colname,q2))
    if nbin > 0:
        if xlims is None:
            xlims = axes.get_xlim()
        ylocnt, xbinedge, _ = stats.binned_statistic(x[z < q1], y[z < q1],
            statistic='count', bins=nbin, range=xlims)
        ylomean, xbinedge, _ = stats.binned_statistic(x[z < q1], y[z < q1],
            statistic='mean', bins=nbin, range=xlims)
        ylostd, xbinedge, _  = stats.binned_statistic(x[z < q1], y[z < q1],
            statistic='std', bins=nbin, range=xlims)
        xbin = 0.5*(xbinedge[1:]+xbinedge[:-1])
        axes.errorbar(xbin[ylocnt>1], ylomean[ylocnt>1], yerr=ylostd[ylocnt>1], 
                    ecolor='k', marker='o', 
                    ms=binms, ls='', mfc=lobin_col, mec='k', zorder=3)
        if len(xbin[ylocnt>1]) > 2:
            b, a, rval, pval, std_err = stats.linregress(xbin[ylocnt>1], ylomean[ylocnt>1])
            print("Lower quartile fit: {:.2f} + x*({:.2f}+/-{:.2f})".format(
                   a, b, std_err))
            #plt.plot(xbin[ylocnt>1], a+b*xbin[ylocnt>1], 'r-', linewidth=3)
        yhicnt, xbinedge, _ = stats.binned_statistic(x[z > q2], y[z > q2],
            statistic='count', bins=nbin, range=xlims)
        yhimean, xbinedge, _ = stats.binned_statistic(x[z > q2], y[z > q2],
            statistic='mean', bins=nbin, range=xlims)
        yhistd, xbinedge, _  = stats.binned_statistic(x[z > q2], y[z > q2],
            statistic='std', bins=nbin, range=xlims)
        axes.errorbar(xbin[yhicnt>1], yhimean[yhicnt>1], yerr=yhistd[yhicnt>1], 
                    ecolor='k', marker='o', 
                    ms=binms, ls='', mfc=hibin_col, mec='k', zorder=3)
        if len(xbin[yhicnt>1]) > 2:
            b, a, rval, pval, std_err = stats.linregress(xbin[yhicnt>1], yhimean[yhicnt>1])
            print("Upper quartile fit: {:.2f} + x*({:.2f}+/-{:.2f})".format(
                   a, b, std_err))
            #plt.plot(xbin[yhicnt>1], a+b*xbin[yhicnt>1], 'b-', linewidth=3)
        if median:
            ymedian, xbinedge, _ = stats.binned_statistic(x, y,
                statistic='median', bins=nbin, range=xlims)
            axes.step(xbinedge , np.append(ymedian, ymedian[-1]), where='post', ls='-', 
                      color='grey', linewidth=4, alpha=0.5, zorder=3)
        # Colorbar annotations
        dnarrow = cb.ax.transLimits.transform([q1,q1])
        print('dnarrow:', dnarrow)
        cb.ax.errorbar( 0.5, dnarrow[1], yerr=0.04, color=lobin_col, uplims=True,
                         elinewidth=2, ecolor='grey', marker='o', mec='k', ms=binms,
                         transform=cb.ax.transAxes)
        uparrow = cb.ax.transLimits.transform([q2,q2])
        print('uparrow:', uparrow)
        cb.ax.errorbar( 0.5, uparrow[1], yerr=0.04, color=hibin_col, lolims=True,
                         elinewidth=2, ecolor='grey', marker='o', mec='k', ms=binms,
                         transform=cb.ax.transAxes)
    return


def std_overlay(cat, axvar, xlims=None, ylims=None, shade=None, axes=None, panel=None):
    '''
    Generate axis labels and custom overlays on log-log plots

    === Parameters ===
    cat : astropy.Table
        Table of cloud properties
    axvar : tuple of str
        Names of the plotted x and y columns
    xlims : tuple of float
        Range of x values to plot model lines, defined logarithmically.  
        For example, [1,3] = 10**1 to 10**3.
    ylims : tuple of float
        Range of y values to plot model lines, defined logarithmically.  
        For example, [1,3] = 10**1 to 10**3.
    shade : dict of name-value pairs
        Regions below these values are shaded in x and y respectively.
        Defined in linear space so the base-10 log is taken.
    axes : matplotlib.axes
        Plot axes object
    '''
    # Axis labels
    labels = { 
        'mvir':'virial mass', 
        'mlumco':'CO-based mass',
        'mlte':'LTE mass',
        'flux12':'Integrated $^{12}$CO flux',
        'flux13':'Integrated $^{13}$CO flux',
        'siglum':'$\Sigma$, CO-based',
        'siglte':'$\Sigma$, LTE-based',
        'sigvir':'$\Sigma$, virial',
        'rad_pc':'equivalent raw radius',
        'rad_pc_dcon':'equivalent radius',
        'vrms_k':'raw rms linewidth',
        'vrms_k_dcon':'rms linewidth',
        'area_pc2':'projected area',
       }
    axlbl=['','']
    for i in [0,1]:
        if axvar[i] in labels.keys():
            axlbl[i] = labels[axvar[i]]
        else:
            axlbl[i] = axvar[i]
    # Plot gray shading indicating resolution limits
    if axes is None:
        axes = plt.gca()
    if xlims is None:
        xlims = axes.get_xlim()
    if ylims is None:
        ylims = axes.get_ylim()
    if shade is not None:
        if axvar[0] in shade.keys():
            axes.axvspan(xlims[0], np.log10(shade[axvar[0]]), fc='lightgray', 
                         alpha=0.3, lw=0)
        if axvar[1] in shade.keys():
            axes.axhspan(ylims[0], np.log10(shade[axvar[1]]), fc='lightgray', 
                         alpha=0.3, lw=0)
    # Solomon et al. size-linewidth relation
    if axvar[0].startswith('rad_pc') and axvar[1].startswith('vrms_k'):
        xmod = np.linspace(xlims[0],xlims[1],20)
        ymod = np.log10(0.72) + 0.5*xmod
        axes.plot(xmod, ymod, linestyle='-', color='r', lw=4, alpha=0.5, 
            zorder=-1)
        axes.text((xlims[1]-0.05), (xlims[1]/2-0.23), 'S87', 
            horizontalalignment='right', color='r', rotation=25)
    # Lines of constant surface density and volume density
    if axvar[0].startswith('rad_pc') and axvar[1].startswith('m'):
        xmod = np.linspace(xlims[0],xlims[1],20)
        # 1, 10, and 100 Msol/pc^2
        axes.plot(xmod, np.log10(np.pi)+2*xmod+0, linestyle=':', color='k', 
            lw=1, zorder=-1)
        axes.plot(xmod, np.log10(np.pi)+2*xmod+1, linestyle=':', color='k', 
            lw=1, zorder=-1)
        axes.plot(xmod, np.log10(np.pi)+2*xmod+2, linestyle=':', color='k', 
            lw=1, zorder=-1)
        # Conversion from n to M/R^3
        ntom = (4*np.pi/3)*(1.36*2*const.m_p.cgs/u.cm**3).to(u.Msun/u.pc**3).value
        # 10**2, 10**3 and 10**4 mol cm^-3
        axes.plot(xmod, np.log10(ntom)+3*xmod+2, alpha=.5, color='r', ls=':', 
            lw=2, zorder=-1)
        axes.plot(xmod, np.log10(ntom)+3*xmod+3, alpha=.5, color='r', ls=':', 
            lw=2, zorder=-1)
        axes.plot(xmod, np.log10(ntom)+3*xmod+4, alpha=.5, color='r', ls=':', 
            lw=2, zorder=-1)
    # Lines of constant surface density
    if axvar[0] == 'area_pc2' and axvar[1].startswith('m'):
        xmod = np.linspace(xlims[0],xlims[1],20)
        axes.plot(xmod, xmod, linestyle=':', color='k', lw=1)
        axes.plot(xmod, xmod+1, linestyle=':', color='k', lw=1)
        axes.plot(xmod, xmod+2, linestyle=':', color='k', lw=1)
    # Lines of constant external pressure
    if axvar[0].startswith('sig') and axvar[1] == 'sigvir':
        xmod = np.linspace(xlims[0],xlims[1],100)
        ymod = np.log10(10**xmod + (20/(3*np.pi*21.1))*1.e4/10**xmod)
        ymod6 = np.log10(10**xmod + (20/(3*np.pi*21.1))*1.e6/10**xmod)
        axes.plot(xmod, ymod6, linestyle='-.', color='brown', lw=1)
        axes.plot(xmod, ymod, linestyle='-', color='g', lw=1)
        if xlims[0] == -1:
            ymod2 = np.log10(10**xmod + (20/(3*np.pi*21.1))*1.e2/10**xmod)
            axes.plot(xmod, ymod2, linestyle=':', color='m', lw=1)
            axes.text(-0.9, 2.40, '$P_{ext}$ = $10^4$ cm$^{-3}$ K', ha='left',
                color='g', rotation=-45)
            axes.text(-0.9, 0.50, '$P_{ext}$ = $10^2$ cm$^{-3}$ K', ha='left',
                color='m', rotation=-45)
    # If axes have identical units then plot y=x line
    if cat[axvar[0]].unit == cat[axvar[1]].unit:
        xmod = np.linspace(xlims[0],xlims[1],20)
        # Plot nominal R(12/13) of 8
        if axvar[0] == 'flux12' and axvar[1] == 'flux13':
            ymod = xmod - np.log10(8)
        else:
            ymod = xmod
        axes.plot(xmod, ymod, linestyle=':', color='k')
        if axvar[0].startswith('sig') and axvar[1] == 'sigvir':
            axes.text(0.05, 0.1, r'$\alpha$ > 1', ha='left',
                color='k', rotation=45, transform=axes.transAxes)            
            axes.text(0.1, 0.05, r'$\alpha$ < 1', ha='left',
                color='k', rotation=45, transform=axes.transAxes)            
    # Set plot limits and labels
    axes.set_xlim(xlims[0], xlims[1])
    axes.set_ylim(ylims[0], ylims[1])
    axes.set_xlabel('log '+axlbl[0]+' ['+str(cat[axvar[0]].unit)+']')
    axes.set_ylabel('log '+axlbl[1]+' ['+str(cat[axvar[1]].unit)+']')
    return

# -------------------------------------------------------------------------------

def model(p, x):
    a, b = p
    return a + b*x

def residuals(p, data):
    a, b = p
    x, y, errx, erry = data
    w = erry*erry + (b)**2*errx*errx
    wi = np.sqrt(np.where(w==0.0, 0.0, 1.0/(w)))
    d = wi*(y-model(p,x))
    return d

def linefitting(x, y, xerr=None, yerr=None, color='b', prob=.95,
                nbootstrap=0, doline=True, zorder=1, parprint=True, xlims=None):
    '''
    Perform linear fitting using kmpfit's effective variance method
    See https://kapteyn.readthedocs.io/en/latest/kmpfittutorial.html

    === Parameters ===
    x : astropy.Column or numpy.array
        The values on the horizontal axis
    y : astropy.Column or numpy.array
        The values on the vertical axis
    xerr : astropy.Column or numpy.array
        Error bars for the horizontal axis
    yerr : astropy.Column or numpy.array
        Error bars for the vertical axis
    color : matplotlib.colors
        Color for best-fit line
    prob : float
        Probability limit for drawing shaded confidence band.
    nbootstrap : int
        Number of bootstrap iterations for estimating errors of slope and
        intercept.  Default=0 (use standard errors).
    doline : boolean
        True to plot best-fit line and confidence band.  Otherwise the
        parameters are calculated but not plotted.
    zorder : int
        Layer for plotting the best-fit line.  Higher numbers show up easier.
    parprint : boolean
        True to show best-fit slope and intercept in upper left corner of plot
    xlims : list of tuple
        The x limits for drawing the fitted lines

    === Returns ===
    d : float
        fitted slope
    d_err : float
        uncertainty in slope.  Based on bootstrap if nbootstrap > 0.
    c : float
        fitted intercept
    c_err : float
        uncertainty in intercept.  Based on bootstrap if nbootstrap > 0.
    fitobj.rchi2_min : float
        reduced chi-squared of the fit.
    yscat : float
        rms of the residual (data - model)
    '''
    # --- Initial guess from simple linear regression
    sorted=np.argsort(x)
    xfit = x[sorted]
    yfit = y[sorted]
    if xerr is not None:
        xerrfit = xerr[sorted]
    else:
        xerrfit = None
    if yerr is not None:
        yerrfit = yerr[sorted]
    else:
        yerrfit = None
    b, a, rval, pval, std_err = stats.linregress(xfit, yfit)
    print("\nLineregress parameters: {:.2f} + x*({:.2f}+/-{:.2f})".format(
           a, b, std_err))
    # --- kmpfit approach
    fitobj = kmpfit.Fitter(residuals=residuals, data=(xfit, yfit, xerrfit, yerrfit))
    fitobj.fit(params0=[a, b])
    print("\n======== Results kmpfit with effective variance =========")
    print("Fitted parameters:      ", fitobj.params)
    print("Covariance errors:      ", fitobj.xerror)
    print("Standard errors:        ", fitobj.stderr)
    print("Chi^2 min:              ", fitobj.chi2_min)
    print("Reduced Chi^2:          ", fitobj.rchi2_min)
    print("Status Message:", fitobj.message)
    c, d = fitobj.params
    c_err, d_err = fitobj.stderr
    yscat = np.std(yfit-model([c, d], xfit))
    # --- Run the bootstrap if requested
    if nbootstrap > 0:
        slopes = []
        offsets = []
        for i in range(nbootstrap):
            indx = randint(0, len(x), len(x))
            # indx is an array of random indices. Use this array to create a new one.
            xfit[:] = x[indx]
            yfit[:] = y[indx]
            if xerr is not None:
                xerrfit[:] = xerr[indx]
            if yerr is not None:
                yerrfit[:] = yerr[indx]
            # Only do a regression if there are at least two different
            # data points in the pseudo sample
            ok = (xfit != xfit[0]).any()
            if (not ok):
                print("All elements are the same. Invalid sample.")
                print(xfit, yfit)
            else:
                fitobj.fit(params0=[a, b])
                offs, slope = fitobj.params
                slopes.append(slope)
                offsets.append(offs)
        slopes = np.array(slopes) - d
        offsets = np.array(offsets) - c
        d_err = slopes.std()
        c_err = offsets.std()
        print("Bootstrap errors:      ", c_err, d_err)
    # --- scipy ODR approach
    linear = odr.Model(model)
    mydata = odr.RealData(x, y, sx=xerr, sy=yerr)
    myodr  = odr.ODR(mydata, linear, beta0=[a,b])
    myoutput = myodr.run()
    print("\n======== Results from scipy.odr =========")
    myoutput.pprint()
    # Plot the results
    axes = plt.gca()
    if doline == True:
        if xlims is None:
            xlims = axes.get_xlim()
        xmod = np.linspace(xlims[0],xlims[1],50)
        ymod = model([c, d], xmod)
        axes.plot(xmod, ymod, linestyle='--', color=color, zorder=zorder)
        dfdp = [1, xmod]
        ydummy, upperband, lowerband = fitobj.confidence_band(xmod,dfdp,prob,model)
        verts = list(zip(xmod, lowerband)) + list(zip(xmod[::-1], upperband[::-1]))
        poly = Polygon(verts, closed=True, fc='c', ec='c', alpha=0.3) 
#            label="{:g}% conf.".format(prob*100))
        axes.add_patch(poly)
    # Overlay the results
    if parprint == True:
        axes.text(0.03,0.95,'$a_1$ = $%4.2f$ ' % d + u'\u00B1' + ' $%4.2f$' % 
            d_err,size=10,transform=axes.transAxes)
        axes.text(0.03,0.90,'$a_0$ = $%4.2f$ ' % c + u'\u00B1' + ' $%4.2f$' % 
            c_err,size=10,transform=axes.transAxes)
    return d, d_err, c, c_err, fitobj.rchi2_min, yscat

# -------------------------------------------------------------------------------

def pltprops(catalog, plotdir='plots', distpc=5e4, dvkms=0.2, beam=2, 
            alpha=1, cmap='jet', parprint=False,
            nbin=0, lobin_col='cyan', hibin_col='salmon',
            xplot=['rad_pc',   'vrms_k','area_pc2'],
            yplot=['vrms_k',   'mlumco',  'mlumco'],
            pltname=[ 'rdv',   'dvflux','areaflux'],
            xlims=[], ylims=[], doline=[], panel=[], ccode=[],
            physnam='physprop', resolve_output='none', 
            colorcodes=['alpha', 'sigvir'], noshow=False):
    '''
    Generate a set of summary plots for an astrodendro run

    === Parameters ===
    catalog : string
        Full path to the csv file with the full catalog.  The 'physprop'
        catalog and the lists of trunks, branches, and leaves are assumed
        to be in the same directory as this file.
    plotdir : string
        Name of folder to save plots.  Default is 'plots' in current directory.
    distpc : float
        Distance of cloud in pc.  Default is 50 kpc for LMC.
    dvkms : float
        Velocity resolution in km/s, used to exclude unresolved structures from fitting.
    beam : float
        Beam size in arseconds, used to exclude unresolved structures from fitting.
    alpha : float
        Transparency parameter for color-coded scatter plots, between 0 and 1.
    cmap : matplotlib.colors.Colormap
        Name of the color map for color coding
    parprint : boolean
        True to show best-fit slope and intercept in upper left corner of plot
    nbin : int
        Number of binned averages to generate across horizontal axis of each plot.
        Default is not to perform binning (nbin=0).
    lobin_col : matplotlib.colors
        Symbol color for lower bin sequence.
    hibin_col : matplotlib.colors
        Symbol color for upper bin sequence.
    xplot : list of str
        Columns to plot on x axis, should be in 'physprop' catalog
    yplot : list of str
        Columns to plot on y axis, should be in 'physprop' catalog
    pltname : list of str
        File name identifiers for the desired plots
    xlims : list of tuple
        The x limits of the desired plots.  List should match length of xplots.
        Default: autoscale
    ylims : list of tuple
        The y limits of the desired plots.  List should match length of xplots.
        Default: autoscale
    doline : list of boolean
        True to plot best-fit line and confidence band for full & cluster plots.
        Default: parameters are calculated but not plotted.
    ccode : list of boolean
        Whether to generate the color coded versions of each plot.  If given, this
        should be an array with the same length as xplot and yplot.
        Default: do not generate these
    colorcodes : list of str
        Columns to plot as color codes, these can be either in 'full' catalog
        or 'physprop' catalog.
    resolve_output: string
        Set to true to output the well-resolved structures to a new output file.
        'raw'  : choose based on all resolution limits
        'decon': choose based on valid deconvolved size only
        'none' : do not output
        Default is 'none'.
    '''

    # ------ Get the resolution limits
    deltav  = dvkms  * u.km / u.s
    avgbeam = beam   * u.arcsec
    dist    = distpc * u.pc
    # Minimum radius is FWHM beam converted to rms size then scaled by 1.91
    rmstorad = 1.91
    radlim = ((avgbeam*rmstorad/np.sqrt(8*np.log(2))) * dist).to(
               u.pc, equivalencies=u.dimensionless_angles())
    # Minimum area is 1 Gaussian beam
    arealim = np.pi/(4*np.log(2)) * ((avgbeam * dist).to(
               u.pc, equivalencies=u.dimensionless_angles()))**2
    # Minimum line width is channel width (~FWHM) divided by 2.35
    dvlim = deltav.value/np.sqrt(8*np.log(2))
    shade = {'rad_pc': radlim.value, 'rad_pc_dcon': radlim.value,
             'vrms_k': dvlim, 'vrms_k_dcon': dvlim,
             'area_pc2': arealim.value}

    # Create directory to place plots
    if not os.path.isdir(plotdir):
        try:
            os.makedirs(plotdir)
        except OSError:
            print('Unable to create output directory',plotdir)

    # ------ Read the input files
    cat = Table.read(catalog, format='ascii.ecsv')
    label = (os.path.basename(catalog)).replace('_full_catalog.txt','')
    indir = os.path.dirname(catalog)
    pcatalog = catalog.replace('_full_catalog.txt','_'+physnam+'_add.txt')
    if os.path.isfile(pcatalog):
        pcat = Table.read(pcatalog, format='ascii.ecsv')
    else:
        pcatalog = catalog.replace('_full_catalog.txt','_'+physnam+'.txt')
        pcat = Table.read(pcatalog, format='ascii.ecsv')

    # Get the indices of trunks, branches, leaves, and clusters.
    # idc[0] is a list of trunk indices
    # idc[1] is a list of branch indices
    # idc[2] is a list of leaf indices
    # idc[3] is a list of cluster indices
    idc=[[],[],[],[]]
    for i, typ in enumerate(['trunks', 'branches', 'leaves', 'clusters']):
        txtfile = join(indir,label+'_'+typ+'.txt')
        try:
            col1 = np.loadtxt(txtfile, usecols=0, dtype=int)
            idc[i] = list(np.atleast_1d(col1))
        except:
            print('{} not found'.format(txtfile))

    # Get the lists of trunk descendants
    try:
        with open(join(indir,label+'_trunks.txt'),'r') as f:
            text=f.read()
            trd = []
            for line in text.splitlines():
                trd.append(list(map(int, line.split('|')[1].split(','))))
    except:
        print('Failed reading',join(indir,label+'_trunks.txt'))
        pass
    
    # Get the lists of cluster descendants and colors
    try:
        with open(join(indir,label+'_clusters.txt'), 'r') as f:
            text=f.read()
            cld = []
            clco = []
            for line in text.splitlines():
                clco.append(line.split()[1])
                if line.split('|')[1].split(',') == [' ']:
                    cld.append([])
                else:
                    cld.append(list(map(int, line.split('|')[1].split(','))))
    except:
        print('Failed reading',join(indir,label+'_clusters.txt'))
        pass

    # ------ Plot histograms of selected properties.  No resolution filtering!
    # First plot is PA and is drawn linearly, subsequent plots are logarithmic.
    hist_struct = ['trunks', 'branches', 'leaves']
    hist_plots  = ['position_angle', 'flux', 'mvir']
    for i, histcol in enumerate(hist_plots):
        fig, axes = plt.subplots()
        histbins  = [12, 10, 10]
        histrange = [[0,180], None, None]
        if histcol in cat.keys():
            histval = cat[histcol]
        elif histcol in pcat.keys():
            histval = pcat[histcol]
        else:
            print('No column {} in {} or {}\n'.format(histcol,catalog,pcatalog))
        pltdata = []
        for j in range(len(hist_struct)):
            xdat = histval[idc[j]]
            if i == 0:
                xdat[xdat<0] += 180.
                pltdata.append(xdat)
            else:
                valid = (xdat > 0)
                pltdata.append(np.log10(xdat[valid]))
        if histrange[i] is None:
            histrange[i] = [np.floor(np.log10(min(histval[histval>0]))), 
                            np.ceil(np.log10(max(histval[histval>0]))) ]
        n, xbinedge, _ = axes.hist(pltdata, bins=histbins[i], range=histrange[i], 
                                   histtype='bar', label=hist_struct)
        axes.set_xlim(histrange[i])
        bin_size = xbinedge[1]-xbinedge[0]
        majorLocator = MultipleLocator(bin_size * 2)
        minorLocator = MultipleLocator(bin_size)
        axes.xaxis.set_major_locator(majorLocator)
        axes.xaxis.set_minor_locator(minorLocator)
        if i == 0:
            axes.set_xlabel(histcol+' ['+str(histval.unit)+']')
        else:
            axes.set_xlabel('log '+histcol+' ['+str(histval.unit)+']')
        axes.set_ylabel('Number')
        plt.legend(loc='best',fontsize='medium')
        plt.savefig(join(plotdir,label+'_'+histcol+'.pdf'), bbox_inches='tight')
        if noshow:
            plt.close()

    # ------ Get color code requests
    if isinstance(colorcodes, str):
        colorcodes = [colorcodes]

    # ------ Output table with line fitting parameters
    tab = Table(dtype=[('cloud', 'S10'), ('pltname', 'S10'), ('npts', 'i4'), 
                    ('a', 'f4'), ('a_err', 'f4'), ('b', 'f4'), ('b_err', 'f4'), 
                    ('chi2red', 'f4'), ('eps', 'f4')])
    for col in ['a', 'a_err', 'b', 'b_err', 'chi2red', 'eps']:
        tab[col].format = '.2f'

    # ------ Main set of scatter plots, as requested by user
    for i in range(len(xplot)):

        print('\n*** Plotting',xplot[i],'and',yplot[i])
        x, y = [pcat[xplot[i]], pcat[yplot[i]]]
        if 'e_'+xplot[i] in pcat.keys():
            xerr = pcat['e_'+xplot[i]]
        else:
            xerr = x*0 + 0.1
        if 'e_'+yplot[i] in pcat.keys():
            yerr = pcat['e_'+yplot[i]]
        else:
            yerr = y*0 + 0.1

        # --- Select points for display: default is all positive values
        dodisp = (x>0) & (y>0) & (xerr>0) & (yerr>0)
        print('{} points selected for display'.format(np.count_nonzero(dodisp)))
        # --- Display unresolved points (dodisp) but exclude from line fitting (dofit)
        dofit = dodisp
        if xplot[i] in shade.keys():
            print('Excluding points from {} below {}'.format(xplot[i],shade[xplot[i]]))
            dofit = dofit & (x > shade[xplot[i]])
        if yplot[i] in shade.keys():
            print('Excluding points from {} below {}'.format(yplot[i],shade[yplot[i]]))
            dofit = dofit & (y > shade[yplot[i]])
            
        # --- Restrict indices of subsets to positive or resolved values
        idsel = idc[:]
        idfit = idc[:]
        for j in range(4):
            idsel[j] = [val for val in idc[j] if val in np.where(dodisp)[0]]
            idfit[j] = [val for val in idc[j] if val in np.where(dofit)[0]]

        # --- Plot trunks, branches, and leaves together
        fig, axes = plt.subplots(figsize=(6.4,4.8))
        if xplot[i] == 'rad_pc' and yplot[i].startswith('m'):
            axes.set_aspect(0.6)
        else:
            axes.set_aspect('equal')
        # Get legend label
        reg  = label.split('_')[0]
        line = label.split('_')[-1]
        if line == '12':
            plt.plot([], [], ' ', label=reg+' CO')
        elif line == '13':
            plt.plot([], [], ' ', label=reg+' $^{13}$CO')
        # Plot the trunks as red pentagons
        sctplot( np.log10(x[idsel[0]]), np.log10(y[idsel[0]]), 
                 xerr=xerr[idsel[0]]/np.log(10), yerr=yerr[idsel[0]]/np.log(10), 
                 col='brown', marker='p', mec='k', ms=50, zorder=4, label='trunks' )
        # Plot the branches as white triangles
        sctplot( np.log10(x[idsel[1]]), np.log10(y[idsel[1]]), 
                 xerr=xerr[idsel[1]]/np.log(10), yerr=yerr[idsel[1]]/np.log(10), 
                 col='w', marker='v', mec='k', ms=17, zorder=2, label='branches' )
        # Plot the leaves as green circles
        sctplot( np.log10(x[idsel[2]]), np.log10(y[idsel[2]]), 
                 xerr=xerr[idsel[2]]/np.log(10), yerr=yerr[idsel[2]]/np.log(10), 
                 col='green', marker='o', mec='k', ms=15, zorder=3, label='leaves' )

        # Set the axis limits
        try:
            this_xlims = xlims[i]
            this_ylims = ylims[i]
        except IndexError:
            this_xlims = axes.get_xlim()
            this_ylims = axes.get_ylim()
        # Decide whether to plot a line
        try:
            drawline = doline[i]
        except IndexError:
            drawline = False
        # Decide whether to label as subpabel
        try:
            subpanel = panel[i]
        except IndexError:
            subpanel = None
        # Decide whether to generate the colorcode plots
        try:
            this_ccode = ccode[i]
        except IndexError:
            this_ccode = False
        
        # Plot the best-fitting line and confidence interval
        if len(x[dofit]) > 2:
            a1, a1_e, a0, a0_e, chi2, eps = linefitting( np.log10(x[dofit]), 
                np.log10(y[dofit]), xerr=xerr[dofit]/np.log(10), 
                yerr=yerr[dofit]/np.log(10), color='b',
                doline=drawline, parprint=parprint, prob=.997, xlims=this_xlims)
            tab.add_row([label, pltname[i], len(x[dofit]), a1, a1_e, a0, a0_e, chi2, eps])
# Old code that would do a fit including unresolved points
#         if pltname[i] == 'rdv':
#             a1, a1_e, a0, a0_e, chi2, eps = linefitting( np.log10(x[dodisp]), 
#                 np.log10(y[dodisp]), xerr=xerr[dodisp]/np.log(10), 
#                 yerr=yerr[dodisp]/np.log(10), color='b',
#                 doline=False, parprint=False)
#             tab.add_row([label, pltname[i]+'all', len(x[dodisp]), a1, a1_e, a0, a0_e, chi2, eps])
        # Plot the binned values if nbin > 0
        if nbin > 0:
            ycnt, xbinedge, _ = stats.binned_statistic(np.log10(x[dodisp]),
                np.log10(y[dodisp]), statistic='count', bins=nbin, range=this_xlims)
            ymean, xbinedge, _ = stats.binned_statistic(np.log10(x[dodisp]), 
                np.log10(y[dodisp]), statistic='mean', bins=nbin, range=this_xlims)
            ystd, xbinedge, _  = stats.binned_statistic(np.log10(x[dodisp]), 
                np.log10(y[dodisp]), statistic='std', bins=nbin, range=this_xlims)
            xbin = 0.5*(xbinedge[1:]+xbinedge[:-1])
            plt.errorbar(xbin[ycnt>1], ymean[ycnt>1], yerr=ystd[ycnt>1], ecolor='orange',
                         marker='o', ls='', mfc='yellow', mec='orange', zorder=10)
        # Make the labels and draw the gray shaded boxes
        std_overlay(pcat, [xplot[i], yplot[i]], this_xlims, this_ylims, shade)
        if subpanel is not None:
            axes.text(0.05,0.93, subpanel, size=12, transform=axes.transAxes)
        plt.legend(loc='lower right',fontsize='small',scatterpoints=1)
        plt.savefig(join(plotdir,label+'_'+pltname[i]+'_full.pdf'), bbox_inches='tight')
        if noshow:
            plt.close()

        # --- Plot trunks and their descendants
        if len(idsel[0]) > 0:
            print('Plotting {} {}(s)'.format(len(idsel[0]),'trunk'))
            fig, axes = plt.subplots(figsize=(6.4,4.8))
            if xplot[i] == 'rad_pc' and yplot[i].startswith('m'):
                axes.set_aspect(0.6)
            else:
                axes.set_aspect('equal')
            sctplot( np.log10(x[idsel[0]]), np.log10(y[idsel[0]]), 
                     xerr=xerr[idsel[0]]/np.log(10), yerr=yerr[idsel[0]]/np.log(10), 
                     marker='s', zorder=4, cmap=cmap, ms=30)
            colors = plt.cm.jet(np.linspace(0, 1, len(idsel[0])) )
            n_descend = 0
            for j, tno in enumerate(idsel[0]):
                descdnts = [val for val in trd[j] if val in np.where(dodisp)[0]]
                n_descend += len(descdnts)
                sctplot( np.log10(x[descdnts]), np.log10(y[descdnts]), 
                         xerr=xerr[descdnts]/np.log(10), yerr=yerr[descdnts]/np.log(10), 
                         col='w', mec=colors[j], zorder=3, ms=10, 
                         label='trunk'+str(tno) )
            print('Plotting {} {} descendants'.format(n_descend,'trunk'))
            std_overlay(pcat, [xplot[i], yplot[i]], this_xlims, this_ylims, shade)
            # Only show legend if there are 10 or fewer trunks
            if len(idsel[0]) <= 10:
                plt.legend(loc='lower right',fontsize='x-small',scatterpoints=1)
            if subpanel is not None:
                axes.text(0.05,0.93, subpanel, size=12, transform=axes.transAxes)
            plt.savefig(join(plotdir,label+'_'+pltname[i]+'_trunks.pdf'), 
                        bbox_inches='tight')
            if noshow:
                plt.close()

        # --- Plot clusters
        if len(idsel[3]) > 0:
            print('Plotting {} {}(s)'.format(len(idsel[3]),'cluster'))
            fig, axes = plt.subplots(figsize=(6.4,4.8))
            if xplot[i] == 'rad_pc' and yplot[i].startswith('m'):
                axes.set_aspect(0.6)
            else:
                axes.set_aspect('equal')
            sctplot( np.log10(x[idsel[3]]), np.log10(y[idsel[3]]), col=clco,
                     xerr=xerr[idsel[3]]/np.log(10), yerr=yerr[idsel[3]]/np.log(10), 
                     marker='s', zorder=4, ms=30)
            if len(x[idfit[3]]) > 2:
                a1, a1_e, a0, a0_e, chi2, eps = linefitting( np.log10(x[idfit[3]]), 
                    np.log10(y[idfit[3]]), xerr=xerr[idfit[3]]/np.log(10), 
                    yerr=yerr[idfit[3]]/np.log(10), color='b',
                    doline=drawline, parprint=parprint, prob=.997, xlims=this_xlims)
                tab.add_row([label, pltname[i]+'_clust', len(x[idfit[3]]), a1, a1_e, 
                            a0, a0_e, chi2, eps])
            std_overlay(pcat, [xplot[i], yplot[i]], this_xlims, this_ylims, shade)
            if subpanel is not None:
                axes.text(0.05,0.93, subpanel, size=12, transform=axes.transAxes)
            plt.savefig(join(plotdir,label+'_'+pltname[i]+'_clusters.pdf'), 
                        bbox_inches='tight')
            if noshow:
                plt.close()

        # ------ Plot all structures, color coded by 3rd variable and binned in quartiles
        if not this_ccode:
            continue
        for j in range(len(colorcodes)):
            fig, axes = plt.subplots(figsize=(6.4,4.8))
            axes.set_aspect('equal')
            print('Relation will be color coded by',colorcodes[j])
            if colorcodes[j] not in cat.keys() and colorcodes[j] not in pcat.keys():
                continue
            if colorcodes[j] in cat.keys():
                zcode = cat[colorcodes[j]][dodisp]
                zlog = False
            elif colorcodes[j] in ['axratio', 'refdist']:
                zcode = pcat[colorcodes[j]][dodisp]
                zlog = False
            else:
                zcode = np.log10(pcat[colorcodes[j]][dodisp])
                zlog = True
            color_code_bin(np.log10(x[dodisp]), np.log10(y[dodisp]), zcode, 
                           colname=colorcodes[j], alpha=alpha, cmap=cmap, axes=axes,
                           lobin_col=lobin_col, hibin_col=hibin_col, xlims=this_xlims, 
                           zlog=zlog, nbin=nbin)
            std_overlay(pcat, [xplot[i], yplot[i]], this_xlims, this_ylims, shade)
            shortname = re.sub('_', '', colorcodes[j])
            plt.savefig(join(plotdir,label+'_'+pltname[i]+'_'+shortname+'.pdf'), 
                        bbox_inches='tight')
            if noshow:
                plt.close()

    # ------ Write output files
    tab.write(join(indir, label+'_lfit.tex'), overwrite=True)    

    if resolve_output != 'none':
        if resolve_output == 'raw':
            mask = ((pcat['rad_pc'] > shade['rad_pc']) & 
                    (pcat['vrms_k'] > shade['vrms_k']) & 
                    (pcat['area_pc2'] > shade['area_pc2']))
        elif resolve_output == 'decon':
            mask = (pcat['rad_pc_dcon'] > 0) & (pcat['vrms_k'] > 0)
        print('Outputting {} out of {} dendros above resolution limits'.format(
              len(pcat[mask]), len(pcat)))
        pcat_out = catalog.replace('_full_catalog.txt','_physprop_resolve.txt')
        pcat[mask].write(pcat_out, format='ascii.ecsv', overwrite=True)
        cat_out = catalog.replace('_full_catalog.txt','_cat_resolve.txt')
        cat[mask].write(cat_out, format='ascii.ecsv', overwrite=True)

    return

