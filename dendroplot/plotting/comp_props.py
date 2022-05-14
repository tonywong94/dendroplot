#!/usr/bin/env python

from .pltprops import linefitting
from scipy.stats import binned_statistic, spearmanr
import numpy as np
import os
import re
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import ticker
from astropy import units as u
from astropy import constants as const
from astropy.table import Table, Column
from matplotlib.colors import Normalize, LogNorm
import matplotlib.lines as mlines

# General Scatter Plot
def plot_ecsv(ecsvfile, xaxis, yaxis, zaxis=None, shade=None, col='g', 
              mark='o', mec='face', zorder=-5, msize=6, linfit=None, 
              label=None, include='all', **kwargs):
    cat = Table.read(ecsvfile, format='ascii.ecsv')
    goodidx = (cat[xaxis]>0) & (cat[yaxis]>0)
    if zaxis is not None:
        goodidx = goodidx & (~np.isnan(cat[zaxis]))
    # Uncomment these 2 lines to exclude unresolved structures from fitting
    if xaxis == 'rad_pc' and yaxis == 'vrms_k':
        goodidx = goodidx & (cat[xaxis]>shade[xaxis]) & (cat[yaxis]>shade[yaxis])
    if include != 'all':
        leavelist=re.sub('physprop\w*', include, ecsvfile)
        idcs = np.loadtxt(leavelist, usecols=0, dtype=int)
        goodidx = np.intersect1d(idcs, np.where(goodidx)[0])
    xdata = np.log10(cat[xaxis][goodidx])
    ydata = np.log10(cat[yaxis][goodidx])
    if 'e_'+xaxis in cat.keys() and xaxis != '8um_avg':
        x_err = cat['e_'+xaxis][goodidx]
    else:
        print('Using uniform error of 0.1 for x axis')
        x_err = np.zeros_like(xdata) + 0.1
    if 'e_'+yaxis in cat.keys() and not yaxis.startswith('sig'):
        y_err = cat['e_'+yaxis][goodidx]
    else:
        print('Using uniform error of 0.1 for y axis')
        y_err = np.zeros_like(ydata) + 0.1
    if zaxis is not None:
        zdata = cat[zaxis][goodidx]
    print('Plotting {} vs {} from file {}'.format(yaxis,xaxis,ecsvfile))
    # Specified colors
    if zaxis is None:
        axes.scatter(xdata, ydata, marker=mark, c=np.reshape(col,(1,-1)), edgecolors=mec, 
            zorder=zorder, s=msize, linewidths=0.1, label=label, **kwargs)
    else:
        axes.scatter(xdata, ydata, marker=mark, c=zdata, edgecolors=mec, 
            zorder=zorder, s=msize, linewidths=0.1, label=label, **kwargs)
    axes.set_aspect('equal')
    mapping = { 'mvir':'virial mass', 
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
    for i, axis in enumerate([xaxis, yaxis]):
        if axis in mapping.keys():
            axlbl[i] = mapping[axis]
        else:
            axlbl[i] = axis
    axes.set_xlabel('log '+axlbl[0]+' ['+str(xdata.unit)+']')
    axes.set_ylabel('log '+axlbl[1]+' ['+str(ydata.unit)+']')
    return np.column_stack((xdata,ydata,x_err,y_err))

# -------------------------------------------------------------------------------

# Main program
def comp_props(dolines, dotypes=['sp8med'], clouds=None, markers=None,
            analdir=None, include='all', binned=False, linefit=True, binfit=False,
            cmap_name='gist_rainbow', msize=10,
            xplot=['rad_pc'],
            yplot=['vrms_k'],
            xlims=[[-1,1.5]],
            ylims=[[-2,1.5]],
            pltname=['rdv'], slope=[0.5], inter=[0], beam=3.5,
            pad=[0.03], magmacsv='islands.sageco.csv'):

    global axes
    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)
    cmap = plt.get_cmap(cmap_name)

    # Read table of cloud-averaged properties
    isle_id = {'N55':441, 'N11B':395, 'N159':94, 'GMC225':2, 'N166':50, 'N171':32,
          'N206':60, 'N206D':46, 'GMC2':259, 'GMC55':224, '30Dor':184,
          'PCC':39, 'N113':127, 'N103B':188, '30DorC':169, 'A439':29,
          'GMC104':57, 'GMC1':165, 'N59C':392}
    magmatab = Table.read(magmacsv, format='ascii.ecsv')
    keep = []
    for clname in clouds:
        keep.append(isle_id[clname]-1)
    cldtab = magmatab[keep]

    # Sensitivity limits
    avgbeam = beam * u.arcsec
    dist    = 5e4 * u.pc
    # Min radius is FWHM beam converted to rms size then scaled by 1.91
    rmstorad = 1.91
    radlim = ((avgbeam*rmstorad/np.sqrt(8*np.log(2))) * dist).to(
        u.pc, equivalencies=u.dimensionless_angles())

    if linefit:
        fittab = Table(dtype=[('line', 'S2'), ('xplot', 'S10'), ('yplot', 'S10'), 
                        ('a1', 'f4'), ('a1_err', 'f4'), ('a0', 'f4'), 
                        ('a0_err', 'f4'), ('chi2red', 'f4'), ('eps', 'f4')])
        for col in ['a1', 'a1_err', 'a0', 'a0_err', 'chi2red', 'eps']:
            fittab[col].format = '.2f'

    # Generate plots
    for i in range(len(xplot)):
        for line in dolines:
            for type in dotypes:
                # Plot a single x-y pair with colorcode
                fig, axes = plt.subplots()
                merge_tbl = np.array([]).reshape(0,4)
                myhandles = []
                for j, clname in enumerate(clouds):
                    dir = analdir.replace('CLOUD', clname)
                    if os.path.isfile(dir+clname+'_'+line+'_physprop_add.txt'):
                        infile = dir+clname+'_'+line+'_physprop_add.txt'
                    else:
                        infile = dir+clname+'_'+line+'_physprop.txt'
                    if clname == '30Dor':
                        deltav = 0.5 * u.km / u.s
                    else:
                        deltav = 0.2 * u.km / u.s
                    # --- Min line width is channel width (~FWHM) divided by 2.35
                    dvlim = deltav.value/np.sqrt(8*np.log(2))
                    shade = {'rad_pc': radlim.value, 'vrms_k': dvlim}
                    if type in ['8um_avg', 'siglum', 'siglte', 'sigvir']:  # local color code
                        if 'sig' in type:
                            norm = LogNorm(vmin=10, vmax=500)
                        else:
                            norm = LogNorm(vmin=1, vmax=100)
                        new = plot_ecsv(infile, xplot[i], yplot[i], shade=shade,
                            zaxis=type, cmap=cmap, mark=markers[j], msize=msize,
                            zorder=i, label=clname, include=include, norm=norm)
                        ccode = 'local'
                        # Dummy handle for legend
                        hdl = mlines.Line2D([], [], color='C0', marker=markers[j], 
                                            ls='', label=clname)
                        myhandles.append(hdl)
                    else:   # cloud-averaged color code
                        if type == 'comean' or type == 'comax':
                            norm = Normalize(vmin=min(cldtab[type]), vmax=max(cldtab[type]))
                        else:
                            norm = LogNorm(vmin=min(cldtab[type]), vmax=max(cldtab[type]))
                        colr=np.array(cmap(norm(cldtab[type][j])))
                        new = plot_ecsv(infile, xplot[i], yplot[i], shade=shade,
                            col=colr, mark=markers[j], msize=msize, zorder=i, 
                            label=clname, include=include)
                        ccode = 'global'            
                    merge_tbl = np.concatenate((merge_tbl, new))
                axes.set_xlim(xlims[i][0], xlims[i][1])
                axes.set_ylim(ylims[i][0], ylims[i][1])
                axes.minorticks_off()
                # Reference slopes
                xmod = np.linspace(xlims[i][0], xlims[i][1], 100)
                ymod = inter[i] + xmod*slope[i]
                if xplot[i] == 'rad_pc' and yplot[i] == 'vrms_k':
                    parprint = True
                    axes.plot(xmod, ymod, linestyle='-', color='r', lw=4, alpha=0.5, 
                        zorder=-1)
                    axes.text((xlims[i][1]-0.05), (xlims[i][1]/2-0.15), 'S87', 
                        horizontalalignment='right', color='r', rotation=25)
                    if type in ['8um_avg', 'sp8med']:
                        lbltype = '8$\mu$m'
                    elif type in ['sp24med']:
                        lbltype = '24$\mu$m'
                    elif type in ['siglum', 'siglte', 'sigvir', 'comean']:
                        lbltype = '$\Sigma$'
                    axes.text(0.5, 1.03, r'$^{'+line+'}$CO colored by '+ccode+' '+
                        lbltype, size=13, ha='center', color='k', transform=axes.transAxes)
                else:
                    parprint = False
                    axes.plot(xmod, ymod, '--', marker=None, color='k')
                # Lines of constant external pressure
                if xplot[i].startswith('sig') and yplot[i] == 'sigvir':
                    ymod = np.log10(10**xmod + (20/(3*np.pi*21.1))*1.e4/10**xmod)
                    axes.plot(xmod, ymod, linestyle='-', color='g', lw=1)
                    axes.text(-0.9, 2.30, '$P_{ext}$ = $10^4$ cm$^{-3}$ K', 
                        color='g', rotation=-45)
                    ymod2 = np.log10(10**xmod + (20/(3*np.pi*21.1))*1.e2/10**xmod)
                    axes.plot(xmod, ymod2, linestyle='-', color='m', lw=1)
                    axes.text(-0.9, 0.90, '$P_{ext}$ = $10^2$ cm$^{-3}$ K', 
                        color='m', rotation=-45)
                    axes.text(-0.9, -0.6, 'Virial Eq', color='k', rotation=45)
                # Fit a line to all points
                sorted=np.argsort(merge_tbl[:,0])
                xdata = merge_tbl[:,0][sorted]
                ydata = merge_tbl[:,1][sorted]
                x_err = merge_tbl[:,2][sorted]/np.log(10)
                y_err = merge_tbl[:,3][sorted]/np.log(10)
                if binned:
                    # Bins are based on data unless plot limits are more constrained
                    makebin = np.linspace(np.amax([xdata[0],xlims[i][0]]), 
                                np.amin([xdata[-1],xlims[i][1]]), 10)
                    mu, edges, asgn = binned_statistic(xdata, ydata, statistic='mean', 
                                        bins=makebin)
                    sig, edges, asgn = binned_statistic(xdata, ydata, statistic='std', 
                                        bins=makebin)
                    x_bins = edges[:-1] + np.diff(edges)/2
                    axes.errorbar(x_bins, mu, yerr=sig, mfc='yellow', mec='k',
                                  ls='None', marker='o', markersize=5,
                                  ecolor='dimgray', capsize=0, zorder=11)
                    if binfit:
                        polyco, cov = np.polyfit(x_bins, mu, 1, w=1/sig, cov=True)
                        a1 = polyco[0]
                        a0 = polyco[1]
                        a1_e, a0_e = np.sqrt(np.diag(cov))
                        xmod = np.linspace(xlims[i][0], xlims[i][1], 10)
                        ymod = a1 * xmod + a0
                        axes.plot(xmod, ymod, linestyle='--', color='g', zorder=10,
                                  label='bin slope=$%4.2f$' % a1)
                        print("Binned slope, intercept:", a1, a0)
                # Fit a line to all points
                if linefit:
                    a1, a1_e, a0, a0_e, chi2, eps = linefitting(xdata, ydata,
                        xerr=x_err, yerr=y_err, xlims=xlims[i], color='b', 
                        parprint=parprint, prob=0.997, zorder=10)
                    fittab.add_row([line, xplot[i], yplot[i], a1, a1_e, a0, a0_e, chi2, eps])
                # Spearman correlation coefficient
                if yplot[i].startswith('sig') and (xplot[i] == '8um_avg' 
                                                   or xplot[i].startswith('sig')):
                    spear, prob = spearmanr(xdata, ydata)
                    if xplot[i].startswith('sig'):
                        axes.text(0.98,0.31,'$^{'+line+'}$CO structures',
                            size=12, color='k', ha='right', transform=axes.transAxes)
                        xpos = 0.7
                    else:
                        xpos = 0.6
                    axes.text(xpos,0.02,'$r_s$={:.2f}'.format(spear),
                          size=10, color='k', ha='right', transform=axes.transAxes)
                    if linefit:
                        axes.text(xpos,0.07,'$a_1$=$%4.2f$' % a1, size=10, color='k',
                            ha='right', transform=axes.transAxes)
                # Legend and colorbar
                if ccode == 'local':
                    plt.legend(loc='lower right',fontsize=8,handles=myhandles)
                else:
                    plt.legend(loc='lower right',fontsize=8,numpoints=1,markerscale=2)
                cax = fig.add_axes([pad[i]+0.7, 0.11, 0.02, 0.77])
                formatter = ticker.LogFormatter(10, labelOnlyBase=False, minor_thresholds=(4,3))
                if type == 'comean':
                    cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
                else:
                    cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, 
                        format=formatter, ticks=[1,2,5,10,20,50,100,200,500,1000])
                cbar.ax.tick_params(labelsize=9)
                if type == 'sp8med':
                    cbartext = 'cloud median 8$\mu$m intensity [MJy/sr]'
                elif type == 'comean':
                    cbartext = 'cloud mean CO intensity [K km/s]'
                elif type == 'comax':
                    cbartext = 'cloud max CO intensity [K km/s]'
                elif type == 'stmean':
                    cbartext = 'cloud mean $\Sigma_{*}$ [$M_\odot$ $pc^{-2}$]'
                elif type == 'sp24med':
                    cbartext = 'cloud median 24$\mu$m intensity [MJy/sr]'
                elif type == '8um_avg':
                    cbartext = 'local mean 8$\mu$m intensity [MJy/sr]'
                elif type == 'siglum':
                    cbartext = 'local mean $\Sigma_{lum}$ [$M_\odot$ $pc^{-2}$]'
                elif type == 'siglte':
                    cbartext = 'local mean $\Sigma_{LTE}$ [$M_\odot$ $pc^{-2}$]'
                elif type == 'sigvir':
                    cbartext = 'local mean $\Sigma_{vir}$ [$M_\odot$ $pc^{-2}$]'
                else:
                    cbartext = ''
                if cbartext != '':
                    cbar.set_label(cbartext, rotation=90)
                plt.savefig('comp_'+line+'_'+pltname[i]+'_'+type+'.pdf', 
                    bbox_inches='tight')
                plt.close()
    if linefit:
        fittab.write('comp_props_lfit.tex', overwrite=True)
    return

