import numpy as np
from astropy import units as u
from astropy import constants as const
from astrodendro import Dendrogram
from astropy.io import fits
from astropy.io.fits import getheader, getdata
import os

def define_asgn(image,dendrogram, label_out='', write=True, dtype=np.short,
                check_structures=False):
    '''
    PURPOSE: Writes and returns a 3D fits file with the index from the dendrogram as the BVALUE
        Required parameters:
            image: location of fits cube
            dendrogram: string/location of dendrogram hdf5 file
        Optional keywords:
            label_out: string to save assignment cube. automatically adds .asgn.fits.gz to end of string, 
                       defaults to {current working directory}+'assignment_cube.asgn.fits.gz' 
            write: set to False to prevent writing fits file
            dtype: data type for assignment cube (default is np.short)
            check_structures: set to True to check if all of the structures in the dendrogram have
                              been written to the assignment cube (will take a long time for large dendrograms)
    '''
    if label_out == '':
        label_out = os.getcwd()+'/assignment_cube'
    d = Dendrogram.load_from(dendrogram)
    cube, hd3 = getdata(image, header=True)
    asgn = np.ones(cube.shape).astype(dtype)
    asgn *= -1
    def recursive_structures(trunks):
        for trunk in trunks:
            j = trunk.idx
            asgn[d[j].get_mask(shape=asgn.shape,subtree=False)] = j
            if len(trunk.children) > 0:
                recursive_structures(trunk.children)
    trunks = d.trunk
    recursive_structures(trunks)
    if check_structures:
        all_structures_present(asgn, d)
    if write:
        new_header = hd3
        new_header.set('BUNIT','Index')
        del new_header['history']
        hdu = fits.PrimaryHDU(asgn)
        hdu.header = new_header
        hdu.header['DATAMIN'] = -1
        hdu.header['DATAMAX'] = len(d) - 1
        hdu.writeto(label_out+'.asgn.fits.gz', overwrite=True)
        print('Wrote assignment cube to ', label_out)
    return asgn

def all_structures_present(asgn, dendrogram):
    '''
    PURPOSE: Check if all structures in the dendrogram are present in the assignment array
        Required parameters:
            asgn: assignment array
            dendrogram: Dendrogram object
    '''
    num = len(np.unique(asgn[asgn > -1]))
    if num != len(dendrogram):
        print('not all structures present in asgn cube')
        return False
    else:
        print('all structures present in asgn cube')
        return True
