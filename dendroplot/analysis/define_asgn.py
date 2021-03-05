import numpy as np
from astropy import units as u
from astropy import constants as const
from astrodendro import Dendrogram
from astropy.io import fits
from astropy.io.fits import getheader, getdata
import matplotlib.pyplot as plt

def define_asgn(indir,target,line,res='2p5as',outdir='',write=True):
    '''
    PURPOSE: Writes and returns a 3D fits file with the index from the dendrogram as the BVALUE
        Required keywords:
            indir: location of dendrogram and image cube
            target: name of object (e.g. '30dor')
            line: '12' for 12CO or '13' for 13CO
        Optional keywords:
            res: angular resolution of cube, defaults to 2.5 arcseconds
            outdir: location to save assignment cube, defaults to input directory
            write: set to False to prevent writing
    '''
    if outdir == '':
        outdir = indir
    label = indir+target+'_'+line
    d = Dendrogram.load_from(label+'_dendrogram.hdf5')
    cube, hd3 = getdata(label+'CO21_2p5as.image.fits.gz', header=True)
    asgn = np.ones(cube.shape)
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
        hdu.writeto(outdir+target+'_'+line+'_'+res+'.asgn.fits.gz')
        print('Wrote assignment cube to ',label+'asgn.fits.gz')
    return asgn