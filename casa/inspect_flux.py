#python 2.7

import numpy as np
from math import *
import matplotlib.pyplot as plt
from astropy.io import fits

def findflux(fitsfile, beamadj=None):
    #flux at each slice of a FITs cube, as well as the total flux of said cube
    data, hd = fits.getdata(fitsfile, header=True)
    data[np.isnan(data)] = 0
    
    flist = []
    for slice in data:
        flux = np.sum(slice)
        flist.append(flux)
    flist = np.array(flist)
    totflux = np.sum(flist)

    if beamadj is None:
        BMAJ, BMIN = hd['BMAJ'], hd['BMIN']
        CDELT1, CDELT2 = hd['CDELT1'], hd['CDELT2']
        beamadj = (pi*BMAJ*BMIN)/(4*log(2)*abs(CDELT1*CDELT2))
        flist = flist/beamadj
        totflux = totflux/beamadj
        
        vstart, vint, vlist = hd['CRVAL3']/1000., hd['CDELT3']/1000., []
        for int in range(len(data)):
            vel = vstart+(int*vint)
            vlist.append(vel)
        vlist = np.array(vlist)
        return flist, totflux, vlist, vint, beamadj

    else:
        return flist/beamadj, totflux/beamadj

def fluxplot(prename = 'GMC1_12CO_12m7m', fileext = 'fits', outfile=None):
    image = prename + '.image.' + fileext
    residual = prename + '.residual.' + fileext
    convmodel = prename + '.convmodel.' + fileext
    
    flist_i, flux_i, vlist, vint, beamadj = findflux(image)
    print "Flux of the Image =                      %4.2f [Jy km/s]" %flux_i
    
    flist_r, flux_r = findflux(residual, beamadj=beamadj)
    print "Flux of the Residual =                   %4.2f [Jy km/s]" %flux_r
    
    flist_c, flux_c = findflux(convmodel, beamadj=beamadj)
    print "Flux of the Convolved Model =            %4.2f [Jy km/s]" %flux_c

    plt.step(vlist, flist_i, 'b', label = 'Image')
    plt.step(vlist, flist_r, 'g', label = 'Residual')
    plt.step(vlist, flist_c, 'r', label = 'Convmodel')
    plt.axhline(y=0, color='k', linestyle='--')
    plt.legend(loc='upper right')
    plt.xlabel("Radio velocity [km/s]")
    plt.ylabel("Flux [Jy]")

    if outfile==None:
        outfile=prename
    plt.savefig(outfile+'.flux_results.pdf')
    return

