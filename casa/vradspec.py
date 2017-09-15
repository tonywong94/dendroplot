#python 2.7

import numpy as np
from math import *
import matplotlib.pyplot as plt
from astropy.io import fits

def findflux(fitsfile, pb=0):
    #spectrum of body, using radio velocity
    data, hd = fits.getdata(fitsfile, header=True)
    data[np.isnan(data)] = 0
    nvel, ncol, nrow= len(data), len(data[0]), len(data[0][0])

    flist, vlist = [], []
    fpblist = [0]*nvel
    if 'image' in fitsfile:
        vstart, vint = hd['CRVAL3']/1000, hd['CDELT3']/1000.
        a = 0
        BMAJ, BMIN = hd['BMAJ'], hd['BMIN']
        CDELT1, CDELT2 = hd['CDELT1'], hd['CDELT2']
        pb = (pi*BMAJ*BMIN)/(4*log(2)*abs(CDELT1*CDELT2))
        
        for i in range(nvel):
            vlist.append(vstart+a)
            a += vint
            for j in range(ncol):
                for k in range(nrow):
                    fpblist[i] += data[i][j][k]                        #(Jy*pixels)/beam
    
        flist, vlist = np.array(fpblist)/pb, np.array(vlist)
        return flist, vlist, pb

    else:
        for i in range(nvel):
            for j in range(ncol):
                for k in range(nrow):
                    fpblist[i] += data[i][j][k]
    
        flist = np.array(fpblist)/pb
        return flist

def vradplot(prename = 'GMC1_12CO_12m7m', fileext = 'fits', sig_v0=235, sig_v1=248):
    image = prename + '.image.' + fileext
    residual = prename + '.residual.' + fileext
    convmodel = prename + '.convmodel.' + fileext
    
    flist_i, vlist, pb = findflux(image)
    flist_r = findflux(residual, pb=pb)
    flist_c = findflux(convmodel, pb=pb)
    
    plt.figure(1)
    lin1 = plt.plot(vlist, flist_i, 'b', label = 'Image')
    lin2 = plt.plot(vlist, flist_r, 'g', label = 'Residual')
    lin3 = plt.plot(vlist, flist_c, 'r', label = 'Convmodel')
    plt.legend(loc='upper right')
    plt.xlim([sig_v0, sig_v1])
    plt.xlabel("Radio velocity [km/s]")
    plt.ylabel("Flux [Jy]")
    plt.savefig(prename+'.vrad_results.pdf')
    return
