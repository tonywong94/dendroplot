#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# Loop over 13CO clusters and identify overlaps with (1) 12CO clusters; (2) 12CO dendros.
# Results are not used, since overlap is complete.

analdir = 'struct/'

clustr13, hd13 = fits.getdata(analdir+"30Dor_feather_mosaic_1p8_13_clusters_asgn.fits.gz", 
                             header=True)
id_c13 = np.unique(clustr13[clustr13>-1])
print('There are',len(id_c13),'13CO clusters:',id_c13)

for type12 in ['clustr', 'dendro']:
    yesover  = np.copy(clustr13)
    notover  = np.copy(clustr13)
    notoverlapping = 0
    yesoverlapping = 0
    if type12 == 'clustr':
        array12 = fits.getdata(analdir+"30Dor_feather_mosaic_1p8_12_clusters_asgn.fits.gz")
    else:
        array12 = fits.getdata(analdir+"30Dor_feather_mosaic_1p8_12.asgn.fits.gz")
    id_12 = np.unique(array12[array12>-1])
    print('There are',len(id_12),'12CO '+type12+'s:',id_12)
    mask3d12 = (array12 > -1).astype(int)
    for idnum in id_c13:
        mask3d13 = (clustr13 == idnum).astype(int)
        overlap  = np.sum(mask3d12*mask3d13)
        if overlap == 0.:
            print('Cluster',idnum,'in 13CO does not overlap a',type12,'in 12CO')
            yesover[yesover == idnum] = -1
            notoverlapping += 1
        else:
            print('Cluster',idnum,'overlaps on',overlap,'pixels')
            notover[notover == idnum] = -1
            yesoverlapping += 1
#     fits.writeto("30Dor_feather_mosaic_1p8_12_clusters_asgn_13"+type13+"_y.fits.gz", 
#                  yesover, hd12, overwrite=True)
#     fits.writeto("30Dor_feather_mosaic_1p8_12_clusters_asgn_13"+type13+"_n.fits.gz", 
#                  notover, hd12, overwrite=True)
