#!/usr/bin/env python

import numpy as np
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
        beamadj = (np.pi*BMAJ*BMIN)/(4*np.log(2)*abs(CDELT1*CDELT2))
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
    txtflux_i = "%.1f" % flux_i
    print("Flux of the Image = ", txtflux_i, "[Jy km/s]")
    flist_r, flux_r = findflux(residual, beamadj=beamadj)
    txtflux_r = "%.1f" % flux_r
    print("Flux of the Residual = ", txtflux_r, "[Jy km/s]")
    flist_c, flux_c = findflux(convmodel, beamadj=beamadj)
    txtflux_c = "%.1f" % flux_c
    print("Flux of the Convolved Model = ", txtflux_c, "[Jy km/s]")

    plt.step(vlist, flist_i, 'b', label = 'Image '+txtflux_i)
    plt.step(vlist, flist_r, 'g', label = 'Residual '+txtflux_r)
    plt.step(vlist, flist_c, 'r', label = 'Convmodel '+txtflux_c)
    plt.axhline(y=0, color='k', linestyle='--')
    plt.legend(loc='upper right')
    plt.xlabel("Radio velocity [km/s]")
    plt.ylabel("Flux [Jy]")

    if outfile==None:
        outfile=prename
    plt.savefig(outfile+'.flux_results.pdf')
    return

