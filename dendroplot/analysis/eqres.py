#!/usr/bin/env python2.7

import numpy as np
import math as m
import matplotlib.pyplot as plt
import radio_beam
from spectral_cube import SpectralCube
from astropy.io import fits
from astropy import units as u
import os
#from reproject import reproject_interp
from astropy.wcs import WCS

"""
    using data cubes and putting them at the same resolution
    
    """

def flatten(unflat, flux, name):
    udata, uhd=fits.getdata(unflat, header=True)
    fdata, fhd=fits.getdata(flux, header=True)
    udata[np.isnan(udata)]=0
    fdata[np.isnan(fdata)]=0
    flat=udata*fdata
    #zerotonan
    flat[flat==0]=np.nan
    name=name+'.image.flat.fits.gz'
    
    return fits.writeto(name, flat, uhd, clobber=True)

"""
def aligner(image_1, image_2, name_1='12CO', name_2='13CO'):
    #align first image to second
    data_1,hd_1=fits.getdata(image_1, header=True)
    data_2,hd_2=fits.getdata(image_2, header=True)
    array, footprint = reproject_interp(image_1, hd_2)
    aligned='../GMC1/'+name_1+'.aligned.fits.gz'
    fits.writeto(aligned,array,hd_1)
    return aligned
"""

"""
    From here, we plug the flattened images into myriad's regrid option
"""

def convolver(image_1, image_2, name_1='12CO', name_2='13CO', size=None):
    resimages=[]
    images=[image_1,image_2]
    names=[name_1,name_2]
    if size!=None:
        # put both images at specified resolution
        newbeam=radio_beam.Beam(major=size*u.arcsec, minor=size*u.arcsec, pa=0*u.deg)
        print newbeam
        for i in range(2):
            cube=SpectralCube.read(images[i])
            print cube
            old_unit = cube._header['bunit']
            print cube.beam
            new_cube=cube.convolve_to(newbeam)
            # Keep the cube in single precision if it was originally
            if cube._header['bitpix'] == -32:
                y = new_cube._data.astype(np.float32)
                new_cube._data = y
            str1 = "{:.1f}".format(size)
            str2 = str1.replace(".0","")
            str3 = str2.replace(".","p")
            convolved=names[i] + '_' + str3 + 'as.image.fits.gz'
            new_cube.write(convolved,format='fits',overwrite=True)
            resimages.append(convolved)
            # fix brightness unit because spectralcube drops 'per beam'
            if old_unit == 'Jy/beam':
                hdulist=fits.open(convolved,mode='update')
                hdulist[0].header['BUNIT']='Jy/beam'
                hdulist.close()
        return resimages
    elif size==None:
        # put first image at resolution of second image
        cube=SpectralCube.read(image_1)
        print cube
        old_unit = cube._header['bunit']
        print cube.beam
        newbeam=radio_beam.Beam.from_fits_header(image_2)
        print newbeam
        new_cube=cube.convolve_to(newbeam)
        # Keep the cube in single precision if it was originally
        if cube._header['bitpix'] == -32:
            y = new_cube._data.astype(np.float32)
            new_cube._data = y
        convolved=name_1+'.convolved.fits.gz'
        new_cube.write(convolved,format='fits',overwrite=True)
        # fix brightness unit because spectralcube drops 'per beam'
        if old_unit == 'Jy/beam':
            hdulist=fits.open(convolved,mode='update')
            hdulist[0].header['BUNIT']='Jy/beam'
            hdulist.close()
        return [convolved,image_2]

def vradspec(image, name):
    #spectrum of body, using radio velocity
    data,hd=fits.getdata(image,header=True)
    data[np.isnan(data)]=0
    #replace nan with 0
    nlayer=len(data)
    #amount of vrad intervals
    ncolumn,nrow=len(data[0]),len(data[0][0])
    fpb_list=[0]*nlayer
    #flux/beam on respective vrad
    
    vrad_list,flux_list=[],[]
    vrad_start=hd['CRVAL3']/1000
    #starting radio velocity in [km/s]
    
    vrad_int=hd['CDELT3']/1000
    #radio velocity interval size in [km/s]
    
    BMAJ,BMIN=hd['BMAJ'],hd['BMIN']
    CDELT1,CDELT2=hd['CDELT1'],hd['CDELT2']
    
    pb=((m.pi)/(4*m.log(2)))*((BMAJ*BMIN)/abs(CDELT1*CDELT2))
    #pixels/beam
    a=0
    for i in range(nlayer):
        vrad_list.append(vrad_start+a)
        a+=vrad_interval
        for j in range(ncolumn):
            for k in range(nrow):
                fpb_list[i]+=data[i][j][k]
    #(Jy*pixels)/beam

    for i in range(len(fpb_list)):
        #Jy/beam -> Jy
        flux_list.append((fpb_list[i]/pb))
    plt.plot(vrad_list,flux_list,'b')
    #Flux v Vrad
    plt.xlabel("Radio velocity [km/s]")
    plt.ylabel("Flux [Jy]")
    plt.savefig(name+'.vradplot.pdf')
    plt.close()
    return


