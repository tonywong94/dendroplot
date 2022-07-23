#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
import csv
from astropy import units as u
from astropy.table import Table
from scipy import stats
from itertools import cycle
from dendroplot.plotting import sctplot
from matplotlib.colors import Normalize, LogNorm
from scipy import odr

# Plot the size-linewidth intercept vs. global cloud properties

def model(p, x):
    a, b = p
    return a + b*x

def calc_intercept(doclouds, dolines):
    xlims = [-1,1.5]
    ylims = [-2,1.5]
    plotx = 'rad_pc'
    ploty = 'vrms_k'
    tablist = []

    for j, run in enumerate(['all', 'leaves']):
        output = Table(data=None,names=('Label','Cloud','Line','Intercept','Uncertainty'),
                   dtype=('S11','S6','S6','f4','f4'))
        tablist.append(output)

    # Sensitivity limits
    avgbeam = 3.5 * u.arcsec
    dist    = 5e4 * u.pc
    # Min radius is FWHM beam converted to rms size then scaled by 1.91
    rmstorad = 1.91
    radlim = ((avgbeam*rmstorad/np.sqrt(8*np.log(2))) * dist).to(
        u.pc, equivalencies=u.dimensionless_angles())
    # Min area is 1 Gaussian beam
    arealim = np.pi/(4*np.log(2)) * ((avgbeam * dist).to(
        u.pc, equivalencies=u.dimensionless_angles()))**2

    # Make the individual cloud plots
    fig, axes = plt.subplots(4, len(doclouds),figsize=(15,10))
    cols = ['g', 'b']

    for i, cldname in enumerate(doclouds):
        if cldname in ['A439', 'GMC1', 'GMC104', 'N59C']:
            f12 = 115.2712
            f13 = 110.2013
        else:
            f12 = 230.538
            f13 = 220.399
        if cldname == '30Dor':
            deltav = 0.5 * u.km / u.s
        else:
            deltav = 0.2 * u.km / u.s
        for k, line in enumerate(dolines):
            if line == '12':
                if f12 > 200:
                    linename = '12CO21'
                else:
                    linename = '12CO10'
            else:
                if f13 > 200:
                    linename = '13CO21'
                else:
                    linename = '13CO10'

            # --- Min line width is channel width (~FWHM) divided by 2.35
            dvlim = deltav.value/np.sqrt(8*np.log(2))
            shade = {'rad_pc': radlim.value, 'vrms_k': dvlim, 'area_pc2': arealim.value}
            # --- Read the data table
            label = cldname+'_'+line
            indir = 'dendro/'
            print('Working on {}'.format(label))
        
            cat = Table.read(indir+label+'_full_catalog.txt', format='ascii.ecsv')
            if os.path.isfile(indir+label+'_physprop_add.txt'):
                pcat = Table.read(indir+label+'_physprop_add.txt', format='ascii.ecsv')
            else:
                pcat = Table.read(indir+label+'_physprop.txt', format='ascii.ecsv')
            leafidx = []
            with open(indir+label+'_leaves.txt', 'r') as f:
                reader=csv.reader(f, delimiter=' ')
                for row in reader:
                    leafidx.append(int(row[0]))

            # --- Set up the plot
            x, y, xerr, yerr = [pcat[plotx], pcat[ploty], pcat['e_'+plotx], pcat['e_'+ploty]]
            for j, run in enumerate(['all', 'leaves']):           
                if j == 0:
                    choose = (x>shade[plotx]) & (y>shade[ploty])   # Must be positive to take logarithm
                else:
                    choose = (x>shade[plotx]) & (y>shade[ploty]) & np.isin(pcat['_idx'], np.array(leafidx))
                axes[2*j+k,i].set_xlim(xlims[0], xlims[1])
                axes[2*j+k,i].set_ylim(ylims[0], ylims[1])
                axes[2*j+k,i].set_aspect('equal')
                axes[2*j+k,i].errorbar( np.log10(x[choose]), np.log10(y[choose]), 
                    xerr=xerr[choose]/np.log(10), yerr=yerr[choose]/np.log(10), 
                    ecolor='dimgray', capsize=0, 
                    zorder=1, marker=None, ls='None', lw=1, label=None)
                sctplot( np.log10(x[choose]), np.log10(y[choose]), 
                    mec='none', ms=30, zorder=2, col=cols[k], axes=axes[2*j+k,i] )
                bmean = np.mean(np.log10(y[choose]) - 0.5*np.log10(x[choose]))
                bstd = np.std(np.log10(y[choose]) - 0.5*np.log10(x[choose]))
                # bstd /= np.sqrt(len(x[choose]))
                xmod = np.linspace(xlims[0],xlims[1],num=10)
                ymod = 0.5*xmod + bmean
                axes[2*j+k,i].plot(xmod,ymod,color='b',linestyle='--' )
                smod = np.log10(0.72) + 0.5*xmod
                axes[2*j+k,i].plot(xmod, smod, linestyle='-', color='r', lw=4, alpha=0.5, zorder=-1)
                axes[2*j+k,i].text((xlims[1]-0.05), (xlims[1]/2-0.1), 'S87', 
                    horizontalalignment='right', color='r', rotation=30)
                axes[2*j+k,i].text((xlims[1]-0.1), (ylims[0]+0.1), run+line+'/'+cldname, 
                    horizontalalignment='right', color='k', fontsize=8)
                print('The mean intercept is {}'.format(bmean))
                tablist[j].add_row([label,cldname,linename,bmean,bstd])

    #fig.subplots_adjust(hspace=0.1)
    fig.savefig('rdv_allclouds.pdf')
    plt.close()
    for j, run in enumerate(['all', 'leaves']):
        tablist[j].write('rdv_intercept_'+run+'.csv',format='ascii.csv',overwrite=True)
    return

# Make the cloud-averaged plots

doclouds = ['30Dor', 'N59C', 'A439', 'GMC104', 'GMC1', 'PCC']
dolines  = ['12','13']

calc_intercept(doclouds, dolines)

params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)
# colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
# colcycle = cycle(colors)

magmafile = 'islands.sageco.csv'
magmatab = Table.read(magmafile, format='ascii.ecsv')
isle_id = {'N55':441, 'N11B':395, 'N159':94, 'GMC225':2, 'N166':50, 'N171':32,
          'N206':60, 'N206D':46, 'GMC2':259, 'GMC55':224, '30Dor':184,
          'PCC':39, 'N113':127, 'N103B':188, '30DorC':169, 'A439':29,
          'GMC104':57, 'GMC1':165, 'N59C':392}
keep = []
for clname in doclouds:
    keep.append(isle_id[clname]-1)
cldtab = magmatab[keep]
norm = LogNorm(vmin=min(cldtab['sp8med']), vmax=max(cldtab['sp8med']))
cmap = plt.get_cmap('gist_rainbow')

for run in ['all', 'leaves']:
    rdvtab = Table.read('rdv_intercept_'+run+'.csv')
    ydata = rdvtab['Intercept']
    xdata = ydata * 0.
    xerr  = ydata * 0.
    yerr = rdvtab['Uncertainty']
    for type in ['sp8', 'sp24', 'coav', 'comax', 'st', 'tdust', 'ndust']:
        colr = []
        for i in range(len(rdvtab)):
            cldname  = rdvtab['Cloud'][i]
            if type == 'sp8':
                xdata[i] = magmatab['sp8med'][isle_id[cldname]-1]
                xerr[i]  = magmatab['sp8std'][isle_id[cldname]-1]
                xlabel = 'median 8$\mu$m intensity'
                xunit = magmatab['sp8med'].unit
            if type == 'sp24':
                xdata[i] = magmatab['sp24med'][isle_id[cldname]-1]
                xerr[i]  = magmatab['sp24std'][isle_id[cldname]-1]
                xlabel = 'median 24$\mu$m intensity'
                xunit = magmatab['sp24med'].unit
            if type == 'coav':
                xdata[i] = magmatab['comean'][isle_id[cldname]-1]
                xerr[i]  = magmatab['costd'][isle_id[cldname]-1]
                xlabel = 'mean MAGMA CO intensity'
                xunit = magmatab['comean'].unit
            if type == 'comax':
                xdata[i] = magmatab['comax'][isle_id[cldname]-1]
                xerr[i]  = magmatab['costd'][isle_id[cldname]-1]
                xlabel = 'maximum MAGMA CO intensity'
                xunit = magmatab['comax'].unit
            if type == 'st':
                xdata[i] = magmatab['stmean'][isle_id[cldname]-1]
                xerr[i]  = magmatab['ststd'][isle_id[cldname]-1]
                xlabel = 'mean stellar surface density'
                xunit = magmatab['stmean'].unit
            if type == 'tdust':
                xdata[i] = magmatab['tdustmean'][isle_id[cldname]-1]
                xerr[i]  = magmatab['tduststd'][isle_id[cldname]-1]
                xlabel = 'mean dust temperature'
                xunit = magmatab['tdustmean'].unit
            if type == 'ndust':
                xdata[i] = magmatab['ndustmean'][isle_id[cldname]-1]
                xerr[i]  = magmatab['nduststd'][isle_id[cldname]-1]
                xlabel = 'mean dust column density'
                xunit = magmatab['ndustmean'].unit
            colr.append(np.array(cmap(norm(magmatab['sp8med'][isle_id[cldname]-1]))))
        # Analyze the correlations
        sorted=np.argsort(xdata)
        m, b, rval, pval, std_err = stats.linregress(np.log10(xdata[sorted]),
            ydata[sorted])
        linear = odr.Model(model)
        mydata = odr.RealData(np.log10(xdata), ydata, sx=xerr/(xdata*np.log(10)), sy=yerr)
        myodr  = odr.ODR(mydata, linear, beta0=[b,m])
        myoutput = myodr.run()
        print("\n======== Results from scipy.odr for",xlabel,run)
        myoutput.pprint()
        spear, prob = stats.spearmanr(np.log10(xdata[sorted]), ydata[sorted])
        xmod = np.logspace(np.log10(xdata[sorted])[0],np.log10(xdata[sorted])[-1],10)
        #ymod = b+m*np.log10(xmod)
        ymod = myoutput.beta[0] + myoutput.beta[1]*np.log10(xmod)

        fig=plt.figure(figsize=(7,5))
        for i in range(0,len(rdvtab),2):

            markcolor = colr[i]
            plt.errorbar(xdata[i],ydata[i],xerr=xerr[i],yerr=yerr[i],
                         ecolor='dimgray', capsize=0, ls='None', lw=1, 
                         label=rdvtab['Label'][i],marker='^',mec=markcolor,mfc=markcolor)
            plt.errorbar(xdata[i+1],ydata[i+1],xerr=xerr[i+1],yerr=yerr[i+1],
                         ecolor='dimgray', capsize=0, ls='None', lw=1, 
                         label=rdvtab['Label'][i+1],marker='s',mec=markcolor,mfc=markcolor)

        axes = plt.gca()
        handles, labels = axes.get_legend_handles_labels()
        # remove the errorbars from legend
        handles = [h[0] for h in handles]
        axes.legend(handles, labels, loc='upper left',ncol=2,fontsize=8)
        axes.set_xscale('log')
        axes.plot(xmod, ymod, linestyle='--', color='k')
        if run == 'leaves':
            axes.text(0.97,0.13,'leaves only',size=12,ha='right',
                    transform=axes.transAxes)
        axes.text(0.97,0.08,'$a_1$ = $%4.2f$' % myoutput.beta[1] + u'\u00B1' + '$%4.2f$'
            % myoutput.sd_beta[1], size=12, ha='right', transform=axes.transAxes)
        axes.text(0.97,0.03,'$r_s$ = $%4.2f$' % spear,size=12,ha='right',
                    transform=axes.transAxes)
        plt.xlabel('cloud '+xlabel+' ['+str(xunit)+']',fontsize=12)
        plt.ylabel('mean log[$\sigma_v/\sqrt{R}$]',fontsize=12)
        plt.savefig('rdv_int_'+type+'_'+run+'.pdf', bbox_inches='tight')

