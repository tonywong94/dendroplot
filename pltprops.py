#!/usr/bin/env python2.7

import csv
import numpy as np
from scipy import stats
from scipy import odr
import os
import re
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon
from astropy import units as u
from astropy import constants as const
from astropy.table import Table, Column
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from kapteyn import kmpfit

# General Scatter Plot
def sctplot(xdata, ydata, zdata=None, col='g', mark='o', mec='k', 
           zorder=-5, msize=6, cmap=None, linfit=None, label=None, **kwargs):
    axes.set_aspect('equal')
    # Single color plot
    if cmap is None and np.size(col) == 1:
        axes.scatter(xdata, ydata, marker=mark, c=col, edgecolors=mec, 
            zorder=zorder, s=msize, linewidths=1, label=label, **kwargs)
    # Multi-color plot, array of colors
    elif cmap is None:
        for xp, yp, cp in zip(xdata, ydata, col):
            axes.scatter(xp, yp, marker=mark, c=cp, edgecolors=mec, 
                zorder=zorder, s=msize, **kwargs)
    # Continuous color map based on structure number
    elif zdata is None:
        xind = np.linspace(0., 1., len(xdata))
        axes.scatter(xdata, ydata, marker=mark, c=xind, zorder=zorder,
            cmap=cmap, edgecolors=mec, s=msize, label=None, **kwargs)
    # Continuous color map based on property value
    else:
        sc=axes.scatter(xdata, ydata, marker=mark, c=zdata, zorder=zorder,
            cmap=cmap, edgecolors=mec, s=msize, label=None, **kwargs)
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
        axes.text(0.03,0.95,'slope = $%4.2f$' % m,size=9,transform=axes.transAxes)
        axes.text(0.03,0.90,'intercept = $%4.2f$' % b,size=9,transform=axes.transAxes)
    return

# -------------------------------------------------------------------------------

def std_overlay(cat, axvar, xlims, ylims, shade=[0,0]):
    # Axis labels
    mapping = { 'mvir':'virial mass', 
        'mlumco':'luminous mass',
        'mlte':'LTE mass',
        'siglum':'$\Sigma$, CO-based',
        'siglte':'$\Sigma$, LTE-based',
        'sigvir':'$\Sigma$, virial',
        'rad_pc':'spherical radius',
        'vrms_k':'rms linewidth',
        'area_pc2':'projected area',
       }
    axlbl=['','']
    for i in [0,1]:
        if axvar[i] in mapping.keys():
            axlbl[i] = mapping[axvar[i]]
        else:
            axlbl[i] = axvar[i]
    # Plot gray shading indicating resolution limits
    if shade[0] > 0:
        axes.axvspan(-3, np.log10(shade[0]), fc='lightgray', alpha=0.3, lw=0)
    if shade[1] > 0:
        axes.axhspan(-3, np.log10(shade[1]), fc='lightgray', alpha=0.3, lw=0)
    # Solomon et al. size-linewidth relation
    if axvar[0] == 'rad_pc' and axvar[1] == 'vrms_k':
        xmod = np.linspace(xlims[0],xlims[1],20)
        ymod = np.log10(0.72) + 0.5*xmod
        axes.plot(xmod, ymod, linestyle='-', color='r', lw=4, alpha=0.5, 
            zorder=-1)
        axes.text((xlims[1]-0.05), (xlims[1]/2-0.1), 'S87', 
            horizontalalignment='right', color='r', rotation=30)
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
        axes.plot(xmod, ymod, linestyle='-', color='g', lw=1)
        axes.text(-0.6, 3.25, '$P_{ext}$ = $10^4$ cm$^{-3}$ K', 
            color='g', rotation=-45)
        ymod2 = np.log10(10**xmod + (20/(3*np.pi*21.1))*1.e2/10**xmod)
        axes.plot(xmod, ymod2, linestyle='-', color='m', lw=1)
        axes.text(-0.95, 1.6, '$P_{ext}$ = $10^2$ cm$^{-3}$ K', 
            color='m', rotation=-45)
    # If axes have identical units then plot y=x line
    if cat[axvar[0]].unit == cat[axvar[1]].unit:
        xmod = np.linspace(xlims[0],xlims[1],20)
        axes.plot(xmod, xmod, linestyle='-', color='k')
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

# -------------------------------------------------------------------------------

def residuals(p, data):
    a, b = p
    x, y, errx, erry = data
    w = erry*erry + (b)**2*errx*errx
    wi = np.sqrt(np.where(w==0.0, 0.0, 1.0/(w)))
    d = wi*(y-model(p,x))
    return d

# -------------------------------------------------------------------------------

# Perform linear fitting using kmpfit's effective variance method
def linefitting(x, y, xerr=None, yerr=None, xrange=[-5, 5], color='b', prob=.95):
    # Initial guess from simple linear regression
    sorted=np.argsort(x)
    b, a, rval, pval, std_err = stats.linregress(x[sorted], y[sorted])
    print "\nLineregress parameters: ", a, b
    # Run the fit
    fitobj = kmpfit.Fitter(residuals=residuals, data=(x, y, xerr, yerr))
    fitobj.fit(params0=[a, b])
    print "\n======== Results kmpfit with effective variance ========="
    print "Fitted parameters:      ", fitobj.params
    print "Covariance errors:      ", fitobj.xerror
    print "Standard errors:        ", fitobj.stderr
    print "Chi^2 min:              ", fitobj.chi2_min
    print "Reduced Chi^2:          ", fitobj.rchi2_min
    print "Status Message:", fitobj.message
    c, d = fitobj.params
    e, f = fitobj.stderr
    # Alternative method using Orthogonal Distance Regression
    linear = odr.Model(model)
    mydata = odr.RealData(x, y, sx=xerr, sy=yerr)
    myodr = odr.ODR(mydata, linear, beta0=[a,b])
    myoutput = myodr.run()
    myoutput.pprint()
    # Plot the results
    xmod = np.linspace(xrange[0],xrange[1],20)
    #ymod0 = model([a, b], xmod)
    ymod = model([c, d], xmod)
    axes.plot(xmod, ymod, linestyle='--', color=color, zorder=1)
    dfdp = [1, xmod]
    ydummy, upperband, lowerband = fitobj.confidence_band(xmod, dfdp, prob, model)
    verts = zip(xmod, lowerband) + zip(xmod[::-1], upperband[::-1])
    poly = Polygon(verts, closed=True, fc='c', ec='c', alpha=0.3, 
        label="{:g}% conf.".format(prob*100))
    axes.add_patch(poly)
    #axes.text(0.03,0.95,'slope = $%4.2f$' % b,size=10,transform=axes.transAxes)
    #axes.text(0.03,0.90,'intercept = $%4.2f$' % a,size=10,transform=axes.transAxes)
    axes.text(0.03,0.95,'slope = $%4.2f$ ' % d + u'\u00B1' + ' $%4.2f$' % 
        f,size=9,transform=axes.transAxes)
    axes.text(0.03,0.90,'intercept = $%4.2f$ ' % c + u'\u00B1' + ' $%4.2f$' % 
        e,size=9,transform=axes.transAxes)
    return

# -------------------------------------------------------------------------------

def pltprops(label, fghz=230.538, distpc=4.8e4, dvkms=0.2, beam=2,
            xplot=['rad_pc', 'vrms_k', 'area_pc2'],
            yplot=['vrms_k', 'mlumco',  'mlumco'],
            xlims=[[-1.5,1],   [-2,2],    [-1,3]],
            ylims=[[-2,1.5], [-1.5,4.5],  [-2,4]],
            pltname=['rdv', 'dvflux', 'areaflux']):

    global axes
    deltav  = dvkms * u.km / u.s
    avgbeam = beam * u.arcsec
    dist    = distpc * u.pc
    freq    = fghz * u.GHz
    # Min radius is FWHM beam converted to rms size then scaled by 1.91
    rmstorad = 1.91
    radlim = ((avgbeam*rmstorad/np.sqrt(8*np.log(2))) * dist).to(
        u.pc, equivalencies=u.dimensionless_angles())
    # Min area is 1 Gaussian beam
    arealim = np.pi/(4*np.log(2)) * ((avgbeam * dist).to(
        u.pc, equivalencies=u.dimensionless_angles()))**2
    # Min line width is channel width (~FWHM) divided by 2.35
    dvlim = deltav.value/np.sqrt(8*np.log(2))
    shade = {'rad_pc': radlim.value, 'vrms_k': dvlim, 'area_pc2': arealim.value}

    # checks/creates directory to place plots
    if os.path.isdir('plots') == 0:
        os.makedirs('plots')

    params = {'text.usetex': False, 'mathtext.fontset': 'stixsans'}
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
    plt.savefig('plots/'+label+'_pahist.pdf', bbox_inches='tight')
    plt.close()

    # Size-linewidth relation, color coded
    plotx = 'rad_pc'
    ploty = 'vrms_k'
    x, y, xerr, yerr = [pcat[plotx], pcat[ploty], pcat['e_'+plotx], pcat['e_'+ploty]]
    # Must be positive to take logarithm
    good = np.intersect1d(np.where(x>0)[0], np.where(y>0)[0])
    z    = ['x_cen', 'y_cen', 'v_cen', 'tpkav', 'siglum']
    cmap = plt.cm.get_cmap('nipy_spectral')
    for i in range(len(z)):
        fig, axes = plt.subplots()
        plt.errorbar( np.log10(x[good]), np.log10(y[good]), 
            xerr=xerr[good]/np.log(10), yerr=yerr[good]/np.log(10), 
            ecolor='dimgray', capsize=0, 
            zorder=1, marker=None, ls='None', lw=1, label=None)
        if z[i] in cat.keys():
            zlbl = z[i]+' ['+str(cat[z[i]].unit)+']'
            sctplot( np.log10(x[good]), np.log10(y[good]), cat[z[i]][good], 
                mec='none', msize=30, zorder=2, cmap=cmap, label=zlbl )
        elif z[i] in pcat.keys():
            zlbl = z[i]+' ['+str(pcat[z[i]].unit)+']'
            sctplot( np.log10(x[good]), np.log10(y[good]), pcat[z[i]][good], 
                mec='none', msize=30, zorder=2, cmap=cmap, label=zlbl )
        std_overlay(pcat, [plotx, ploty], xlims[0], ylims[0], 
            [shade['rad_pc'],shade['vrms_k']])
        shortname = re.sub('_', '', z[i])
        plt.savefig('plots/'+label+'_rdv_'+shortname+'.pdf', bbox_inches='tight')
        plt.close()

    # Main set of scatter plots, as requested by user
    for i in range(len(xplot)):
        x, y, xerr, yerr = [pcat[xplot[i]], pcat[yplot[i]], pcat['e_'+xplot[i]], 
            pcat['e_'+yplot[i]]]
        # Must be positive to take logarithm
        good = np.intersect1d(np.where(x>0)[0], np.where(y>0)[0])
        # Restrict indices of subsets to good values
        idsel = idc[:]
        for j in range(4):
            idsel[j] = [val for val in idc[j] if val in good.tolist()]
        # Exclude unresolved points from line fitting
        xmin = ymin = 0
        if xplot[i] in shade.keys():
            if shade[xplot[i]] > 0:
                xmin = shade[xplot[i]]
        if yplot[i] in shade.keys():
            if shade[yplot[i]] > 0:
                ymin = shade[yplot[i]]
        unshade = np.intersect1d(np.where(x>xmin)[0], np.where(y>ymin)[0])
        # --- Plot trunks, branches, leaves
        fig, axes = plt.subplots()
        # Plot the error bars of all points in gray
        plt.errorbar( np.log10(x[good]), np.log10(y[good]), 
            xerr=xerr[good]/np.log(10), yerr=yerr[good]/np.log(10), 
            ecolor='dimgray', capsize=0, 
            zorder=1, marker=None, ls='None', lw=1, label=None)
        # Plot the trunks as red pentagons
        sctplot ( np.log10(x[idsel[0]]), np.log10(y[idsel[0]]), col='brown',
            mark='p', mec='k', msize=80, zorder=4, label='trunks' )
        # Plot the branches as white triangles
        sctplot ( np.log10(x[idsel[1]]), np.log10(y[idsel[1]]), col='w',
            mark='v', mec='k', msize=17, zorder=2, label='branches' )
        # Plot the leaves as green circles
        sctplot ( np.log10(x[idsel[2]]), np.log10(y[idsel[2]]), col='green',
            mark='o', mec='k', msize=20, zorder=3, label='leaves' )
        # Plot the best-fitting line and confidence interval
        linefitting( np.log10(x[unshade]), np.log10(y[unshade]), 
            xerr=xerr[unshade]/np.log(10), yerr=yerr[unshade]/np.log(10), 
            xrange=xlims[i], color='b' )
        # Make the labels and draw the gray shaded boxes
        std_overlay(pcat, [xplot[i], yplot[i]], xlims[i], ylims[i], [xmin,ymin])
        plt.legend(loc='lower right',fontsize='small',scatterpoints=1)
        plt.savefig('plots/'+label+'_'+pltname[i]+'_full.pdf', bbox_inches='tight')
        plt.close()
        # --- Plot trunks and their descendants
        fig, axes = plt.subplots()
        cmap = plt.cm.get_cmap('jet')
        plt.errorbar( np.log10(x[idsel[0]]), np.log10(y[idsel[0]]), 
            xerr=xerr[idsel[0]]/np.log(10), yerr=yerr[idsel[0]]/np.log(10), 
            ecolor='dimgray', capsize=0, 
            zorder=1, marker=None, ls='None', lw=1, label=None)
        sctplot ( np.log10(x[idsel[0]]), np.log10(y[idsel[0]]), mark='s', 
            zorder=4, cmap=cmap, msize=30 )
        colors = plt.cm.jet(np.linspace(0, 1, len(idsel[0])))
        for j, tno in enumerate(idsel[0]):
            trsel = trd[j][:]
            trsel = [val for val in trd[j] if val in good.tolist()]
            plt.errorbar( np.log10(x[trsel]), np.log10(y[trsel]), 
                xerr=xerr[trsel]/np.log(10), yerr=yerr[trsel]/np.log(10),
                ecolor='dimgray', capsize=0, 
                zorder=1, marker=None, ls='None', lw=1, label=None)
            sctplot ( np.log10(x[trsel]), np.log10(y[trsel]), col='w', 
                mec=colors[j], zorder=3, msize=10, label='trunk'+str(tno) )
        std_overlay(pcat, [xplot[i], yplot[i]], xlims[i], ylims[i], [xmin,ymin])
        # Only show legend if there are 10 or fewer trunks
        if len(idsel[0]) <= 10:
            plt.legend(loc='lower right',fontsize='x-small',scatterpoints=1)
        plt.savefig('plots/'+label+'_'+pltname[i]+'_trunks.pdf', bbox_inches='tight')
        plt.close()
        # --- Plot clusters and their descendants (get marker color from table)
        fig, axes = plt.subplots()
        plt.errorbar( np.log10(x[idsel[3]]), np.log10(y[idsel[3]]), 
            xerr=xerr[idsel[3]]/np.log(10), yerr=yerr[idsel[3]]/np.log(10), 
            ecolor='dimgray', capsize=0, 
            zorder=2, marker=None, ls='None', lw=1, label=None)
        sctplot ( np.log10(x[idsel[3]]), np.log10(y[idsel[3]]), mark='s', 
            zorder=4, col=clco, msize=25 )
        unshade2 = idsel[3][:]
        unshade2 = [val for val in idsel[3] if val in unshade.tolist()]
        # Plot best-fitting line to clusters only
        linefitting( np.log10(x[unshade2]), np.log10(y[unshade2]), 
            xerr=xerr[unshade2]/np.log(10), yerr=yerr[unshade2]/np.log(10), 
            xrange=xlims[i], color='b' )
        for j, tno in enumerate(idsel[3]):
            clsel = cld[j][:]
            clsel = [val for val in cld[j] if val in good.tolist()]
            sctplot ( np.log10(x[clsel]), np.log10(y[clsel]), col='w', 
                mec=clco[j], zorder=3, msize=10, label='cluster'+str(tno), alpha=0.5 )
        std_overlay(pcat, [xplot[i], yplot[i]], xlims[i], ylims[i], [xmin,ymin])
        # Only show legend if there are 9 or fewer clusters
        if len(idsel[3]) <= 9:
            plt.legend(loc='lower right',fontsize='x-small',scatterpoints=1)
        plt.savefig('plots/'+label+'_'+pltname[i]+'_clusters.pdf', bbox_inches='tight')
        plt.close()
    
    return
