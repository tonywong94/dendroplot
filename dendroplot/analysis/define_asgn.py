import numpy as np
from astropy import units as u
from astropy import constants as const
from astrodendro import Dendrogram
from astropy.io import fits
from astropy.io.fits import getheader, getdata
import os

def define_asgn(image,dendrogram,label_out='',write=True):
    '''
    PURPOSE: Writes and returns a 3D fits file with the index from the dendrogram as the BVALUE
        Required keywords:
            image: location of fits cube
            dendrogram: location of dendrogram hdf5 file
        Optional keywords:
            label_out: string to save assignment cube. automatically adds .asgn.fits.gz to end of string, 
                        defaults to 'current working directory + assignment_cube.asgn.fits.gz' 
            write: set to False to prevent writing fits file
    '''
    if label_out == '':
        label_out = os.getcwd()+'/assignment_cube'
    d = Dendrogram.load_from(dendrogram)
    cube, hd3 = getdata(image, header=True)
    #alternatively, we only need to load the header:
    #hd3 = getheader(image)
    #asgn = np.ones((hd3['NAXIS3'], hd3['NAXIS2'], hd3['NAXIS1'])).astype(np.float32)
    asgn = np.ones(cube.shape).astype(np.float32)
    asgn[:] = np.NaN
    for j in range(len(d)):
        if j % 10 == 1:
            print("Adding indices for structure {} to assignment cube".format(j))
        asgn[d[j].get_mask(shape = asgn.shape,subtree=False)] = d[j].idx
    if write:
        new_header = hd3
        new_header.set('BUNIT','Index')
        hdu = fits.PrimaryHDU(asgn)
        hdu.header = new_header
        hdu.writeto(label_out+'.asgn.fits.gz')
        print('Wrote assignment cube to ',label_out)
    return asgn
