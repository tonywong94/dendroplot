#!/usr/bin/env python

import os
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy import units as u

import warnings
warnings.filterwarnings('ignore')

# Calculate the fraction of cluster mass at high alpha (alpha > 10)

clustonly = True
analdir = 'struct/'
select = 'physprop_resolve'
mos = analdir+'30Dor_feather_mosaic_1p8_12_'+select+'.txt'

xlims   = [0, 320]

for line in ['12', '13']:
    if line == '13':
        mos = mos.replace('_12_','_13_')
        mass = 'mlte'
        allploty = ['sigvir', 'alphalte']
    else:
        mass = 'mlumco'
        allploty = ['sigvir', 'alpha']
    pcat = Table.read(mos, format='ascii.ecsv')
    pcat.add_index('_idx')
    mosname = os.path.basename(mos)
    print('\nWorking on', mosname)

    if clustonly:
        cllist = mos.replace(select,'clusters')
        clust = np.loadtxt(cllist, usecols=0, dtype=int)
        goodclust = np.intersect1d(pcat['_idx'], clust)
        subcat = pcat.loc[goodclust]
        pcat = subcat
    valid = (pcat[mass] > 0)
    if line == '13':
        high_alpha = valid & (np.log10(pcat['alphalte']) > 1)
        low_alpha = valid & (np.log10(pcat['alphalte']) < 0.5)
    else:
        high_alpha = valid & (np.log10(pcat['alpha']) > 1)
        low_alpha = valid & (np.log10(pcat['alphalte']) < 0.5)
    print('Mass at log alpha > 1:',np.sum(pcat[mass][high_alpha])*u.solMass)
    print('Mass at log alpha < 0.5:',np.sum(pcat[mass][low_alpha])*u.solMass)
    print('Total Mass:',np.sum(pcat[mass][valid])*u.solMass)
    print('High Fraction:',np.sum(pcat[mass][high_alpha])/np.sum(pcat[mass][valid]))
    print('Low Fraction:',np.sum(pcat[mass][low_alpha])/np.sum(pcat[mass][valid]))
    ke_tot = np.sum(3*pcat[mass][valid]*pcat['vrms_k'][valid]**2)*u.solMass*u.km**2/u.s**2
    ke_hi = np.sum(3*pcat[mass][high_alpha]*pcat['vrms_k'][high_alpha]**2)*u.solMass*u.km**2/u.s**2
    ke_lo = np.sum(3*pcat[mass][low_alpha]*pcat['vrms_k'][low_alpha]**2)*u.solMass*u.km**2/u.s**2
    print('K.E. at log alpha > 1:',ke_hi.cgs)
    print('K.E. at log alpha < 0.5:',ke_lo.cgs)
    print('Total K.E.:',ke_tot.cgs)
    power = 1.2e39 * u.erg/u.s  # Bestenlehner+20
    print('Time needed for log alpha > 1:',(ke_hi/power).to(u.yr))
    print('Time needed for log alpha < 0.5:',(ke_lo/power).to(u.yr))
    print('Time needed for all:',(ke_tot/power).to(u.yr))
    
