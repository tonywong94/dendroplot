import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.table import Table, Column
from os.path import expanduser, join
from dendroplot.plotting import pltprops, sctplot, std_overlay
from astropy.io import fits

# Plot associations between 12CO and 13CO dendros.
# Must first run cluster12_overlap.py.
# Sensitivity limits calculated in "Tex floor.ipynb" in analysis/lte.

# Which lines
domos   = ['30Dor_feather_mosaic_']
dolines = ['12', '13']
res     = '1p8'

# Plotting parameters
pltcol = ['orange', 'green', 'blue', 'black']
pltmrk = ['s', '^', 'o', 'o']
z      = [3, 4, 5, 6]
pltms  = [20, 30, 30, 30]
pltlbl = ['no $^{13}$CO match', '$^{13}$CO dendro match', '$^{13}$CO clump match', 'filament match']

dendir = 'struct'
fildir = 'struct'

for mos in domos:

    for line in dolines:

        # ------ Read the text files
        physprop_ext = '_physprop_resolve'
        pcatalog = join(dendir, mos+res+'_'+line+physprop_ext+'.txt')
        label = (os.path.basename(pcatalog)).replace(physprop_ext+'.txt','')
        pcat = Table.read(pcatalog, format='ascii.ecsv')
        pcat.add_index('_idx')

        # ------ Read the assignment images
        if line == '12':
            xplot = 'siglum'
            # 12CO clusters with 13CO cluster counterparts
            yesover = fits.getdata(dendir+'/'+mos+res+'_12_clusters_asgn_13clustr_y.fits.gz')
            id_yes = np.intersect1d(pcat['_idx'], np.unique(yesover[yesover>-1]))
            # 12CO clusters with 13CO dendro counterparts
            dendover = fits.getdata(dendir+'/'+mos+res+'_12_clusters_asgn_13dendro_y.fits.gz')
            id_yes2 = np.intersect1d(pcat['_idx'], np.unique(dendover[dendover>-1]))
            # 12CO clusters with 13CO dendro (but not cluster) counterparts
            id_den = np.setdiff1d(id_yes2,id_yes)
            print('Number of co clusters with 13co counterparts',len(id_yes2))
            # 12CO clusters with no 13CO dendro counterparts
            notover = fits.getdata(dendir+'/'+mos+res+'_12_clusters_asgn_13dendro_n.fits.gz')
            id_no = np.intersect1d(pcat['_idx'], np.unique(notover[notover>-1]))
            print('Number of co clusters with no 13co counterparts',len(id_no))
        else:
            xplot = 'siglte'
            allclust = fits.getdata(dendir+'/'+mos+res+'_13_clusters_asgn.fits.gz')
            id_13 = np.intersect1d(pcat['_idx'], np.unique(allclust[allclust>-1]))
        yplot = 'sigvir'

        # ------ Make the plots
        fig, axes = plt.subplots(figsize=(6.4,4.8))
        axes.set_aspect('equal')
        if line == '12':
            plt.plot([], [], ' ', label='CO clumps')
            allarr = [id_no, id_den, id_yes]
            
            for i, arr in enumerate(allarr):
                xp = pcat.loc[arr][xplot]
                yp = pcat.loc[arr][yplot]
                xerr = pcat.loc[arr]['e_'+xplot]
                yerr = pcat.loc[arr]['e_'+yplot]
                sctplot( np.log10(xp), np.log10(yp), col=pltcol[i],
                         xerr=xerr/np.log(10), yerr=yerr/np.log(10), 
                         marker=pltmrk[i], zorder=z[i], ms=pltms[i], label=pltlbl[i])

            axes.axvline(0.55, ls='--', color='dimgrey')
            axes.text(0.4,3.3,'CO sens.',rotation=90)
            outfile = label+'_bnd_clusters_on13.pdf'
            print('Mean siglum for 13co counterparts:',np.mean(np.log10(pcat.loc[id_yes2][xplot])))
            print('Mean siglum for no 13co counterparts:',np.mean(np.log10(pcat.loc[id_no][xplot])))
            print('Mean alpha for 13co counterparts:',np.mean(np.log10(pcat.loc[id_yes2][yplot])-np.log10(pcat.loc[id_yes2][xplot])))
            print('Mean alpha for no 13co counterparts:',np.mean(np.log10(pcat.loc[id_no][yplot])-np.log10(pcat.loc[id_no][xplot])))
        else:
            sctplot( np.log10(pcat.loc[id_13][xplot]), 
                     np.log10(pcat.loc[id_13][yplot]), col='red',
                     xerr=pcat.loc[id_13]['e_'+xplot]/np.log(10), 
                     yerr=pcat.loc[id_13]['e_'+yplot]/np.log(10), 
                     marker='o', zorder=4, ms=30, label='$^{13}$CO clumps')
            axes.axvline(1.5, ls='--', color='dimgrey')
            axes.text(1.55,0.1,'13CO sens.',rotation=90)
            outfile = label+'_bndlte_clusters.pdf'
        std_overlay(pcat, [xplot, yplot], [0,4], [0,4])
        plt.legend(handletextpad=0.1)
        axes.text(0.1, 2.5, '$10^4$ cm$^{-3}$ K', ha='left',
            color='g', rotation=-45)
        axes.text(0.8, 3.3, '$10^6$ cm$^{-3}$ K', ha='left',
            color='brown', rotation=-45)
        plt.savefig(outfile, bbox_inches='tight')
        plt.close()
