#!/usr/bin/env python

import numpy as np
from numpy.random import randint
from scipy import stats
from scipy import odr
import os
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

# General Scatter Plot
def sctplot(xdata, ydata, zdata=None, col='g', mark='o', mec='k', 
           zorder=-5, msize=6, cmap=None, linfit=None, label=None, axes=None,
           **kwargs):
    if axes is None:
        axes = plt.gca()
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

def std_overlay(cat, axvar, xlims, ylims, shade=[0,0], axes=None):
    # Axis labels
    mapping = { 
        'mvir':'virial mass', 
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
    for i in [0,1]:
        if axvar[i] in mapping.keys():
            axlbl[i] = mapping[axvar[i]]
        else:
            axlbl[i] = axvar[i]
    # Plot gray shading indicating resolution limits
    if axes is None:
        axes = plt.gca()
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
        axes.text((xlims[1]-0.05), (xlims[1]/2-0.15), 'S87', 
            horizontalalignment='right', color='r', rotation=25)
    # Lines of constant surface density and volume density
    if axvar[0] == 'rad_pc' and axvar[1].startswith('m'):
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
        axes.plot(xmod, ymod, linestyle='-', color='g', lw=1)
        axes.text(-0.9, 2.40, '$P_{ext}$ = $10^4$ cm$^{-3}$ K', ha='left',
            color='g', rotation=-45)
        ymod2 = np.log10(10**xmod + (20/(3*np.pi*21.1))*1.e2/10**xmod)
        axes.plot(xmod, ymod2, linestyle=':', color='m', lw=1)
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
# a: regression intercept
# b: regression slope
# c: fitted intercept
# d: fitted slope
def linefitting(x, y, xerr=None, yerr=None, xrange=[-5, 5], color='b', prob=.95,
                nbootstrap=0, zorder=1, doline=True, parprint=True):
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
    #print(np.amin(xfit),np.amax(xfit),np.amin(yfit),np.amax(yfit))
    #print(np.amin(xerrfit),np.amax(xerrfit),np.amin(yerrfit),np.amax(yerrfit))
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
    myodr = odr.ODR(mydata, linear, beta0=[a,b])
    myoutput = myodr.run()
    print("\n======== Results from scipy.odr =========")
    myoutput.pprint()
    # Plot the results
    axes = plt.gca()
    if doline == True:
        xmod = np.linspace(xrange[0],xrange[1],50)
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

def pltprops(label, distpc=5e4, dvkms=0.2, beam=2, alpha=1, nbin=16, colmap='jet',
            xplot=['rad_pc', 'vrms_k', 'area_pc2'],
            yplot=['vrms_k', 'mlumco',  'mlumco'],
            xlims=[[-1.5,1],   [-2,2],    [-1,3]],
            ylims=[[-2,1.5], [-1.5,4.5],  [-2,4]],
            pltname=['rdv', 'dvflux', 'areaflux']):

    deltav  = dvkms  * u.km / u.s
    avgbeam = beam   * u.arcsec
    dist    = distpc * u.pc
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

    # Get the indices of trunks, branches, leaves, and clusters.
    # idc[0] is a list of trunk indices
    # idc[1] is a list of branch indices
    # idc[2] is a list of leaf indices
    # idc[3] is a list of cluster indices
    idc=[[],[],[],[]]
    for i, typ in enumerate(['trunks', 'branches', 'leaves', 'clusters']):
        try:
            col1 = np.loadtxt(label+'_'+typ+'.txt', usecols=0, dtype=int)
            idc[i] = list(np.atleast_1d(col1))
        except:
            print('{} not found'.format(label+'_'+typ+'.txt'))

    # Get the lists of trunk descendants
    try:
        with open(label+'_trunks.txt','r') as f:
            text=f.read()
            trd = []
            for line in text.splitlines():
                trd.append(list(map(int, line.split('|')[1].split(','))))
    except:
        pass
    
    # Get the lists of cluster descendants and colors
    try:
        with open(label+'_clusters.txt', 'r') as f:
            text=f.read()
            cld = []
            clco = []
            for line in text.splitlines():
                cld.append(list(map(int, line.split('|')[1].split(','))))
                clco.append(line.split()[1]) 
    except:
        pass

    # Histogram of PAs
    val = 'position_angle'
    types = ['trunks', 'branches', 'leaves']
    bin_size = 15; min_edge = 0; max_edge = 180
    N = int((max_edge-min_edge)/bin_size)
    bin_list = np.linspace(min_edge, max_edge, N+1)
    pltdata = []
    for i in range(len(types)):
        xdata = cat[val][idc[i]]
        xdata[xdata<0] += 180.
        pltdata.append(xdata)
    fig, axes = plt.subplots()
    axes.hist(pltdata, bin_list, histtype='bar', label=types)
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
    postive = (x>0) & (y>0)
    z    = ['x_cen', 'y_cen', 'v_cen', 'tpkav', 'siglum', 'sigvir', '8um_avg', 'refdist']
    cmap = plt.cm.get_cmap(colmap)
    for i in range(len(z)):
        print('Size-linewidth relation color coded by',z[i])
        if z[i] not in cat.keys() and z[i] not in pcat.keys():
            continue
        fig, axes = plt.subplots()
        axes.set_aspect('equal')
#         plt.errorbar( np.log10(x[postive]), np.log10(y[postive]), 
#             xerr=xerr[postive]/np.log(10), yerr=yerr[postive]/np.log(10), 
#             ecolor='dimgray', capsize=0, 
#             zorder=1, marker=None, ls='None', lw=.5, label=None)
        if z[i] in cat.keys():
            zlbl = z[i]+' ['+str(cat[z[i]].unit)+']'
            sctplot( np.log10(x[postive]), np.log10(y[postive]), cat[z[i]][postive], 
                mec='none', msize=6, zorder=2, cmap=cmap, label=zlbl, alpha=alpha )
            q1 = np.nanpercentile(cat[z[i]], 25)
            q2 = np.nanpercentile(cat[z[i]], 75)
            print('Quartiles for {} go from {} to {}'.format(z[i],q1,q2))
            low = (cat[z[i]] < q1)
            high = (cat[z[i]] > q2)
        elif z[i] in pcat.keys():
            zlbl = z[i]+' ['+str(pcat[z[i]].unit)+']'
            sctplot( np.log10(x[postive]), np.log10(y[postive]), pcat[z[i]][postive], 
                mec='none', msize=6, zorder=2, cmap=cmap, label=zlbl, alpha=alpha )
            q1 = np.nanpercentile(pcat[z[i]], 25)
            q2 = np.nanpercentile(pcat[z[i]], 75)
            print('Quartiles for {} go from {} to {}'.format(z[i],q1,q2))
            low = (pcat[z[i]] < q1)
            high = (pcat[z[i]] > q2)
        if nbin > 0:
            lowsel = (postive) & (low)
            hisel = (postive) & (high)
            ymean, xbinedge, _ = stats.binned_statistic(np.log10(x[lowsel]), 
                np.log10(y[lowsel]), statistic='mean', bins=nbin, range=xlims[0])
            ystd, xbinedge, _  = stats.binned_statistic(np.log10(x[lowsel]), 
                np.log10(y[lowsel]), statistic='std', bins=nbin, range=xlims[0])
            xbin = 0.5*(xbinedge[1:]+xbinedge[:-1])
            plt.errorbar(xbin, ymean, yerr=ystd, ecolor='k', marker='o', 
                         ls='', mfc='salmon', mec='k', zorder=10)
            ymean, xbinedge, _ = stats.binned_statistic(np.log10(x[hisel]), 
                np.log10(y[hisel]), statistic='mean', bins=nbin, range=xlims[0])
            ystd, xbinedge, _  = stats.binned_statistic(np.log10(x[hisel]), 
                np.log10(y[hisel]), statistic='std', bins=nbin, range=xlims[0])
            plt.errorbar(xbin, ymean, yerr=ystd, ecolor='k', marker='o', 
                         ls='', mfc='cyan', mec='k', zorder=10)
        std_overlay(pcat, [plotx, ploty], xlims[0], ylims[0], 
            [shade['rad_pc'],shade['vrms_k']])
        shortname = re.sub('_', '', z[i])
        plt.savefig('plots/'+label+'_rdv_'+shortname+'.pdf', bbox_inches='tight')
        plt.close()

    # Main set of scatter plots, as requested by user
    tab = Table(dtype=[('cloud', 'S10'), ('pltname', 'S10'), ('a', 'f4'), 
                    ('a_err', 'f4'), ('b', 'f4'), ('b_err', 'f4'), 
                    ('chi2red', 'f4'), ('eps', 'f4')])
    for col in ['a', 'a_err', 'b', 'b_err', 'chi2red', 'eps']:
        tab[col].format = '.2f'
    for i in range(len(xplot)):
        print('\nPlotting',xplot[i],'and',yplot[i])
        x, y = [pcat[xplot[i]], pcat[yplot[i]]]
        if 'e_'+xplot[i] in pcat.keys():
            xerr = pcat['e_'+xplot[i]]
        else:
            xerr = x*0 + 0.1
        if 'e_'+yplot[i] in pcat.keys():
            yerr = pcat['e_'+yplot[i]]
        else:
            yerr = y*0 + 0.1
        # --- Must be positive to take logarithm
        postive = (x>0) & (y>0) & (xerr>0) & (yerr>0)
        # --- Restrict indices of subsets to positive values
        idsel = idc[:]
        for j in range(4):
            idsel[j] = [val for val in idc[j] if val in np.where(postive)[0].tolist()]
        # --- Exclude unresolved points from line fitting
        xmin = ymin = 0
        if xplot[i] in shade.keys():
            print('Excluding points from {0} below {1}'.format(xplot[i],shade[xplot[i]]))
            if shade[xplot[i]] > 0:
                xmin = shade[xplot[i]]
        if yplot[i] in shade.keys():
            print('Excluding points from {0} below {1}'.format(yplot[i],shade[yplot[i]]))
            if shade[yplot[i]] > 0:
                ymin = shade[yplot[i]]
        unshade = (postive) & (x > xmin) & (y > ymin)
        #
        # --- Plot trunks, branches, leaves
        fig, axes = plt.subplots()
        if xplot[i] == 'rad_pc' and yplot[i].startswith('m'):
            axes.set_aspect(0.6)
        else:
            axes.set_aspect('equal')
        # Get plot label
        reg = label.split('_')[0].upper()
        if reg == '30DOR':
            reg = '30Dor'
        line = label.split('_')[1]
        if line == '12':
            plt.plot([], [], ' ', label=reg+' CO')
        elif line == '13':
            plt.plot([], [], ' ', label=reg+' $^{13}$CO')
        # Plot the error bars of all points in gray
        plt.errorbar( np.log10(x[postive]), np.log10(y[postive]), 
            xerr=xerr[postive]/np.log(10), yerr=yerr[postive]/np.log(10), 
            ecolor='dimgray', capsize=0, 
            zorder=1, marker=None, ls='None', lw=1, label=None)
        # Plot the trunks as red pentagons
        sctplot ( np.log10(x[idsel[0]]), np.log10(y[idsel[0]]), col='brown',
            mark='p', mec='k', msize=50, zorder=4, label='trunks' )
        # Plot the branches as white triangles
        sctplot ( np.log10(x[idsel[1]]), np.log10(y[idsel[1]]), col='w',
            mark='v', mec='k', msize=17, zorder=2, label='branches' )
        # Plot the leaves as green circles
        sctplot ( np.log10(x[idsel[2]]), np.log10(y[idsel[2]]), col='green',
            mark='o', mec='k', msize=15, zorder=3, label='leaves' )
        # Plot the best-fitting line and confidence interval
        if pltname[i] not in ['bnd', 'bndlte']:
            if len(x[unshade]) > 2:
                a1, a1_e, a0, a0_e, chi2, eps = linefitting( np.log10(x[unshade]), 
                    np.log10(y[unshade]), xerr=xerr[unshade]/np.log(10), 
                    yerr=yerr[unshade]/np.log(10), xrange=xlims[i], color='b',
                    doline=True, parprint=False, prob=.997)
                tab.add_row([label, pltname[i], a1, a1_e, a0, a0_e, chi2, eps])
            if pltname[i] == 'rdv':
                a1, a1_e, a0, a0_e, chi2, eps = linefitting( np.log10(x[postive]), 
                    np.log10(y[postive]), xerr=xerr[postive]/np.log(10), 
                    yerr=yerr[postive]/np.log(10), xrange=xlims[i], color='b',
                    doline=False, parprint=False)
                tab.add_row([label, pltname[i]+'all', a1, a1_e, a0, a0_e, chi2, eps])
        # Plot the binned values if nbin > 0
        if nbin > 0:
            ymean, xbinedge, _ = stats.binned_statistic(np.log10(x[postive]), 
                np.log10(y[postive]), statistic='mean', bins=nbin, range=xlims[i])
            ystd, xbinedge, _  = stats.binned_statistic(np.log10(x[postive]), 
                np.log10(y[postive]), statistic='std', bins=nbin, range=xlims[i])
            xbin = 0.5*(xbinedge[1:]+xbinedge[:-1])
            plt.errorbar(xbin, ymean, yerr=ystd, ecolor='orange', marker='o', 
                         ls='', mfc='yellow', mec='orange', zorder=10)
        # Make the labels and draw the gray shaded boxes
        std_overlay(pcat, [xplot[i], yplot[i]], xlims[i], ylims[i], [xmin,ymin])
        plt.legend(loc='lower right',fontsize='small',scatterpoints=1)
        plt.savefig('plots/'+label+'_'+pltname[i]+'_full.pdf', bbox_inches='tight')
        plt.close()
        #
        # --- Plot trunks and their descendants
        fig, axes = plt.subplots()
        if xplot[i] == 'rad_pc' and yplot[i].startswith('m'):
            axes.set_aspect(0.6)
        else:
            axes.set_aspect('equal')
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
            trsel = [val for val in trd[j] if val in postive.tolist()]
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
        #
        # --- Plot clusters and their descendants (get marker color from table)
        if len(idsel[3]) > 0:
            fig, axes = plt.subplots()
            if xplot[i] == 'rad_pc' and yplot[i].startswith('m'):
                axes.set_aspect(0.6)
            else:
                axes.set_aspect('equal')
            plt.errorbar( np.log10(x[idsel[3]]), np.log10(y[idsel[3]]), 
                xerr=xerr[idsel[3]]/np.log(10), yerr=yerr[idsel[3]]/np.log(10), 
                ecolor='dimgray', capsize=0, 
                zorder=2, marker=None, ls='None', lw=1, label=None)
            sctplot ( np.log10(x[idsel[3]]), np.log10(y[idsel[3]]), mark='s', 
                zorder=4, col=clco, msize=25 )
            unshade2 = [val for val in idsel[3] if val in unshade.tolist()]
            # Plot best-fitting line to clusters only
            if pltname[i] not in ['bnd', 'bndlte']:
                if len(unshade2) > 2:
                    a1, a1_e, a0, a0_e, eps, chi2 = linefitting( np.log10(x[unshade2]), 
                        np.log10(y[unshade2]), xerr=xerr[unshade2]/np.log(10), 
                        yerr=yerr[unshade2]/np.log(10), xrange=xlims[i], color='b' )
            for j, tno in enumerate(idsel[3]):
                clsel = cld[j][:]
                clsel = [val for val in cld[j] if val in postive.tolist()]
                sctplot ( np.log10(x[clsel]), np.log10(y[clsel]), col='w', 
                    mec=clco[j], zorder=3, msize=10, label='cluster'+str(tno), alpha=0.5 )
            std_overlay(pcat, [xplot[i], yplot[i]], xlims[i], ylims[i], [xmin,ymin])
            # Only show legend if there are 9 or fewer clusters
            if len(idsel[3]) <= 9:
                plt.legend(loc='lower right',fontsize='x-small',scatterpoints=1)
            plt.savefig('plots/'+label+'_'+pltname[i]+'_clusters.pdf', bbox_inches='tight')
            plt.close()
    tab.write(label+'_lfit.tex', overwrite=True)
    
    return

# -------------------------------------------------------------------------------
# Hybrid plot with trunks and branches from 12CO and leaves from 13CO.
# -------------------------------------------------------------------------------

def plthybrid(cld, distpc=5e4, dvkms=0.2, beam=2,
            xplot=['rad_pc', 'vrms_k', 'area_pc2'],
            yplot=['vrms_k', 'mlumco',  'mlumco'],
            xlims=[[-1.5,1],   [-2,2],    [-1,3]],
            ylims=[[-2,1.5], [-1.5,4.5],  [-2,4]],
            pltname=['rdv', 'dvflux', 'areaflux']):

    deltav  = dvkms  * u.km / u.s
    avgbeam = beam   * u.arcsec
    dist    = distpc * u.pc
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

    if os.path.isfile(cld+'_12_physprop_add.txt'):
        pcat12 = Table.read(cld+'_12_physprop_add.txt', format='ascii.ecsv')
    else:
        pcat12 = Table.read(cld+'_12_physprop.txt', format='ascii.ecsv')
    newcol = Column(pcat12['area_pc2']*0., name='e_area_pc2')
    newcol.unit = 'pc2'
    pcat12.add_column(newcol)
    if os.path.isfile(cld+'_13_physprop_add.txt'):
        pcat13 = Table.read(cld+'_13_physprop_add.txt', format='ascii.ecsv')
    else:
        pcat13 = Table.read(cld+'_13_physprop.txt', format='ascii.ecsv')
    newcol = Column(pcat13['area_pc2']*0., name='e_area_pc2')
    newcol.unit = 'pc2'
    pcat13.add_column(newcol)

    # Get the indices of trunks, branches, leaves, and clusters.
    # idc[0] is a list of trunk indices
    # idc[1] is a list of branch indices
    # idc[2] is a list of leaf indices
    idc=[[],[],[]]
    for i, typ in enumerate(['trunks', 'branches']):
        col1 = np.loadtxt(cld+'_12_'+typ+'.txt', usecols=0, dtype=int)
        idc[i] = list(np.atleast_1d(col1))
    col1 = np.loadtxt(cld+'_13_leaves.txt', usecols=0, dtype=int)
    idc[2] = list(np.atleast_1d(col1))

    # Main set of scatter plots, as requested by user
    tab = Table(dtype=[('cloud', 'S10'), ('pltname', 'S10'), ('a', 'f4'), 
                    ('a_err', 'f4'), ('b', 'f4'), ('b_err', 'f4'), 
                    ('chi2red', 'f4'), ('eps', 'f4')])
    for col in ['a', 'a_err', 'b', 'b_err', 'chi2red', 'eps']:
        tab[col].format = '.2f'
    for i in range(len(xplot)):
        x = np.concatenate((pcat12[xplot[i]][idc[0]],
                            pcat12[xplot[i]][idc[1]],
                            pcat13[xplot[i]][idc[2]]))
        y = np.concatenate((pcat12[yplot[i]][idc[0]],
                            pcat12[yplot[i]][idc[1]],
                            pcat13[yplot[i]][idc[2]]))
        if 'e_'+xplot[i] in pcat12.keys() and 'e_'+xplot[i] in pcat13.keys():
            xerr = np.concatenate((pcat12['e_'+xplot[i]][idc[0]],
                                   pcat12['e_'+xplot[i]][idc[1]],
                                   pcat13['e_'+xplot[i]][idc[2]]))
        else:
            xerr = x*0 + 0.1
        if 'e_'+yplot[i] in pcat12.keys() and 'e_'+yplot[i] in pcat13.keys():
            yerr = np.concatenate((pcat12['e_'+yplot[i]][idc[0]],
                                   pcat12['e_'+yplot[i]][idc[1]],
                                   pcat13['e_'+yplot[i]][idc[2]]))
        else:
            yerr = y*0 + 0.1
        # --- Must be positive to take logarithm
        postive = (x > 0) & (y > 0)
        # --- Restrict indices of subsets to positive values
        idsel = idc[:]
        pos12 = (pcat12[xplot[i]] > 0) & (pcat12[yplot[i]] > 0)
        pos13 = (pcat13[xplot[i]] > 0) & (pcat13[yplot[i]] > 0)
        idsel[0] = [val for val in idc[0] if val in np.where(pos12)[0].tolist()]
        idsel[1] = [val for val in idc[1] if val in np.where(pos12)[0].tolist()]
        idsel[2] = [val for val in idc[2] if val in np.where(pos13)[0].tolist()]
        # --- Exclude unresolved points from line fitting
        xmin = ymin = 0
        if xplot[i] in shade.keys():
            print('Excluding points from {0} below {1}'.format(xplot[i],shade[xplot[i]]))
            if shade[xplot[i]] > 0:
                xmin = shade[xplot[i]]
        if yplot[i] in shade.keys():
            print('Excluding points from {0} below {1}'.format(yplot[i],shade[yplot[i]]))
            if shade[yplot[i]] > 0:
                ymin = shade[yplot[i]]
        unshade = (x > xmin) & (y > ymin)
        # --- Plot trunks, branches, leaves
        fig, axes = plt.subplots()
        if xplot[i] == 'rad_pc' and yplot[i].startswith('m'):
            axes.set_aspect(0.6)
        else:
            axes.set_aspect('equal')
        # Get plot label
        plt.plot([], [], ' ', label=cld)
        # Plot the error bars of all points in gray
        plt.errorbar( np.log10(x[postive]), np.log10(y[postive]), 
            xerr=xerr[postive]/np.log(10), yerr=yerr[postive]/np.log(10), 
            ecolor='dimgray', capsize=0, 
            zorder=1, marker=None, ls='None', lw=1, label=None)
        # Plot the trunks as red pentagons
        sctplot ( np.log10(pcat12[xplot[i]][idsel[0]]), 
                  np.log10(pcat12[yplot[i]][idsel[0]]), 
                  col='brown', mark='p', mec='k', msize=80, zorder=4, label='CO trunks' )
        # Plot the branches as white triangles
        sctplot ( np.log10(pcat12[xplot[i]][idsel[1]]), 
                  np.log10(pcat12[yplot[i]][idsel[1]]), 
                  col='w', mark='v', mec='k', msize=17, zorder=2, label='CO branches' )
        # Plot the leaves as green circles
        sctplot ( np.log10(pcat13[xplot[i]][idsel[2]]), 
                  np.log10(pcat13[yplot[i]][idsel[2]]), 
                  col='green', mark='o', mec='k', msize=20, zorder=3, label='$^{13}$CO leaves' )
        # Plot the best-fitting line and confidence interval
        if pltname[i] not in ['bnd', 'bndlte']:
            if len(x[unshade]) > 2:
                a1, a1_e, a0, a0_e, chi2, eps = linefitting( np.log10(x[unshade]), 
                    np.log10(y[unshade]), xerr=xerr[unshade]/np.log(10), 
                    yerr=yerr[unshade]/np.log(10), xrange=xlims[i], color='b',
                    doline=True, parprint=True, prob=.997)
                tab.add_row([cld, pltname[i], a1, a1_e, a0, a0_e, chi2, eps])
            if pltname[i] == 'rdv':
                a1, a1_e, a0, a0_e, chi2, eps = linefitting( np.log10(x[postive]), 
                    np.log10(y[postive]), xerr=xerr[postive]/np.log(10), 
                    yerr=yerr[postive]/np.log(10), xrange=xlims[i], color='b',
                    doline=False, parprint=False)
                tab.add_row([cld, pltname[i]+'all', a1, a1_e, a0, a0_e, chi2, eps])
        # Make the labels and draw the gray shaded boxes
        std_overlay(pcat12, [xplot[i], yplot[i]], xlims[i], ylims[i], [xmin,ymin])
        plt.legend(loc='lower right',fontsize='small',scatterpoints=1)
        plt.savefig('plots/'+cld+'_hyb_'+pltname[i]+'_full.pdf', bbox_inches='tight')
        plt.close()

    tab.write(cld+'_hyb_lfit.tex', overwrite=True)
    return
