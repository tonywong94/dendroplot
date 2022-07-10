#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# Loop over 12CO clusters and identify overlaps with (1) 13CO clusters; (2) 13CO dendros.
# Results are plotted in clusterbnd.py

analdir = 'struct/'

clustr12, hd12 = fits.getdata(analdir+"30Dor_feather_mosaic_1p8_12_clusters_asgn.fits.gz", 
                             header=True)
id_c12 = np.unique(clustr12[clustr12>-1])
print('There are',len(id_c12),'12CO clusters:',id_c12)

for type13 in ['clustr', 'dendro']:
    yesover  = np.copy(clustr12)
    notover  = np.copy(clustr12)
    notoverlapping = 0
    yesoverlapping = 0
    if type13 == 'clustr':
        array13 = fits.getdata(analdir+"30Dor_feather_mosaic_1p8_13_clusters_asgn.fits.gz")
    else:
        array13 = fits.getdata(analdir+"30Dor_feather_mosaic_1p8_13.asgn.fits.gz")
    id_13 = np.unique(array13[array13>-1])
    print('There are',len(id_13),'13CO '+type13+'s:',id_13)
    mask3d13 = (array13 > -1).astype(int)
    for idnum in id_c12:
        mask3d12 = (clustr12 == idnum).astype(int)
        overlap  = np.sum(mask3d12*mask3d13)
        if overlap == 0.:
            print('Cluster',idnum,'in 12CO does not overlap a',type13,'in 13CO')
            yesover[yesover == idnum] = -1
            notoverlapping += 1
        else:
            print('Cluster',idnum,'overlaps on',overlap,'pixels')
            notover[notover == idnum] = -1
            yesoverlapping += 1
    fits.writeto(analdir+"30Dor_feather_mosaic_1p8_12_clusters_asgn_13"+type13+"_y.fits.gz", 
                 yesover, hd12, overwrite=True)
    fits.writeto(analdir+"30Dor_feather_mosaic_1p8_12_clusters_asgn_13"+type13+"_n.fits.gz", 
                 notover, hd12, overwrite=True)
