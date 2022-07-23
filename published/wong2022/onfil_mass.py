#!/usr/bin/env python

import os
import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy import units as u

import warnings
warnings.filterwarnings('ignore')

# Calculate the fraction of mass on and off filaments.
# This should be run after Dendros_on_fils.ipynb

analdir = 'struct/'

for line in ['12', '13']:
    pcat = Table.read(analdir+'30Dor_feather_mosaic_1p8_'+line+
                      '_physprop_resolve.txt', format='ascii.ecsv')
    pcat.add_index('_idx')
    for dendrotyp in ['clusters', 'leaves']:
        asgn_all = analdir+'30Dor_feather_mosaic_1p8_'+line+'_'+dendrotyp+'_asgn.fits.gz'
        alldat = fits.getdata(asgn_all)
        asgn_onfil = analdir+'30Dor_feather_mosaic_1p8_'+line+'_'+dendrotyp+'_asgn_onfil.fits.gz'
        yesover = fits.getdata(asgn_onfil)
        id_all = np.intersect1d(pcat['_idx'], np.unique(alldat[alldat>-1]))
        id_yes = np.intersect1d(pcat['_idx'], np.unique(yesover[yesover>-1]))
        id_no  = np.setdiff1d(id_all,id_yes)

        if line == '13':
            mass = 'mlte'
        else:
            mass = 'mlumco'
        print('\nWorking on', os.path.basename(asgn_all))

        allcat = pcat.loc[id_all]
        print('Number of',dendrotyp,'is',len(allcat[mass]))
        oncat  = pcat.loc[id_yes]
        print('Number on filaments is',len(oncat[mass]))
        if len(id_no)>0:
            offcat = Table(pcat.loc[id_no])
            print('Number off filaments is',len(offcat[mass]))
        print('Total Mass:',np.sum(allcat[mass])*u.solMass)
        print('Mass on filaments:',np.sum(oncat[mass])*u.solMass)
        print('Frac on filaments:',np.sum(oncat[mass])/np.sum(allcat[mass]))
        if len(id_no)>0:
            print('Mass off filaments:',np.sum(offcat[mass])*u.solMass)
            print('Frac off filaments:',np.sum(offcat[mass])/np.sum(allcat[mass]))
