#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
import csv
from astropy import units as u
from astropy.table import Table
from dendroplot.plotting import sctplot, std_overlay

## Plot the specific kinetic energy as a function of mass or radius

doclouds = ['30Dor', 'N59C', 'A439', 'GMC104', 'GMC1', 'PCC']
dolines  = ['12', '13']
indir    = 'dendro/'

for xplot in ['mass', 'radius']:

    fig, ax = plt.subplots(2, len(doclouds), figsize=(15,8))

    for i, cldname in enumerate(doclouds):
        if cldname in ['A439', 'GMC1', 'GMC104', 'N59C']:
            f12 = 115.2712
            f13 = 110.2013
        else:
            f12 = 230.538
            f13 = 220.399
        for k, line in enumerate(dolines):
            if line == '12':
                deltav = 0.2 * u.km / u.s
                if f12 > 200:
                    linename = '12CO21'
                else:
                    linename = '12CO10'
            else:
                deltav = 0.5 * u.km / u.s
                if f13 > 200:
                    linename = '13CO21'
                else:
                    linename = '13CO10'

            label = cldname+'_'+line
            print('Working on {}'.format(label))
        
            if os.path.isfile(indir+label+'_physprop_add.txt'):
                pcat = Table.read(indir+label+'_physprop_add.txt', format='ascii.ecsv')
            else:
                pcat = Table.read(indir+label+'_physprop.txt', format='ascii.ecsv')
            types=['trunks','branches','leaves']
            idc=[0,0,0]
            rad=[[],[],[]]
            mass=[[],[],[]]
            energy=[[],[],[]]
            for ii, typ in enumerate(types):
                idc[ii] = []
                with open(indir+label+'_'+typ+'.txt', 'r') as f:
                    reader=csv.reader(f, delimiter=' ')
                    for row in reader:
                        idc[ii].append(int(row[0]))
                rad[ii] = pcat['rad_pc'][idc[ii]]
                mass[ii] = pcat['mlumco'][idc[ii]]
                energy[ii] = pcat['vrms_k'][idc[ii]]**2
            if xplot == 'mass':
                sctplot ( np.log10(mass[0]), np.log10(energy[0]), col='brown',
                    marker='p', mec='k', ms=80, zorder=4, label='trunks', axes=ax[k,i] )
                sctplot ( np.log10(mass[1]), np.log10(energy[1]), col='w',
                    marker='v', mec='k', ms=17, zorder=2, label='branches', axes=ax[k,i] )
                sctplot ( np.log10(mass[2]), np.log10(energy[2]), col='green',
                    marker='o', mec='k', ms=20, zorder=3, label='leaves', axes=ax[k,i] )
                ax[k,i].set_xlim(0, 6)
                ax[k,i].set_ylim(-2.2, 2)
            else:
                sctplot ( np.log10(rad[0]), np.log10(energy[0]), col='brown',
                    marker='p', mec='k', ms=80, zorder=4, label='trunks', axes=ax[k,i] )
                sctplot ( np.log10(rad[1]), np.log10(energy[1]), col='w',
                    marker='v', mec='k', ms=17, zorder=2, label='branches', axes=ax[k,i] )
                sctplot ( np.log10(rad[2]), np.log10(energy[2]), col='green',
                    marker='o', mec='k', ms=20, zorder=3, label='leaves', axes=ax[k,i] )
                ax[k,i].set_xlim(-0.8, 1.5)
                ax[k,i].set_ylim(-2.2, 2)
            if i>0:
                ax[k,i].set_yticklabels([])
            ax[k,i].set_title(label)
    if xplot == 'mass':
        ax[1,0].set_xlabel('log mass [solMass]')
    else:
        ax[1,0].set_xlabel('log radius [pc]')
    ax[1,0].set_ylabel('log specific K.E. [km$^2$/s$^2$]')
    ax[1,0].legend(loc='lower right',fontsize='small',scatterpoints=1)
    plt.savefig('energy_'+xplot+'.pdf', bbox_inches='tight')
    plt.close()

