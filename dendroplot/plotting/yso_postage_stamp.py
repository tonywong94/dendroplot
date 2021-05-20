import numpy as np
from spectral_cube import SpectralCube
from astropy import units as u
from matplotlib import pyplot as plt
from astropy.io.fits import getheader, getdata
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord as SkyCoord
import astropy

def yso_postage_stamp(cube_location, YSO_RA, YSO_DEC, RA_width=1*u.arcmin, DEC_width=1*u.arcmin,slice_=None):
    cube,header = (getdata(cube_location,header=True))
    wcs = WCS(header)
    xl = min([YSO_RA-RA_width, YSO_RA+RA_width])
    xh = max([YSO_RA-RA_width, YSO_RA+RA_width])
    yl = min([YSO_DEC-DEC_width, YSO_DEC+DEC_width])
    yh = max([YSO_DEC-DEC_width, YSO_DEC+DEC_width])
    if slice_ is not None:
        spec_cube = SpectralCube(data=cube, wcs=wcs)
        zl = spec_cube.spectral_axis[slice_]
        zh = spec_cube.spectral_axis[slice_]
        cutout = spec_cube.subcube(xlo=xl,xhi=xh,ylo=yl,yhi=yh,zlo = zl, zhi = zh)
        ax = plt.subplot(projection=wcs, slices=('x', 'y', slice_))
        im = plt.imshow(cutout[0].value,extent=[xl.value,xh.value, yl.value,yh.value])
        plt.scatter(YSO_RA, YSO_DEC, s=10, c='red', marker='o')
        plt.colorbar(im)
        ax.set_aspect('equal')
        plt.title(spec_cube.spectral_axis[slice_])
        plt.show()
    else:
        cutout = astropy.nddata.utils.Cutout2D(cube, SkyCoord(YSO_RA, YSO_DEC), (RA_width,DEC_width), wcs=wcs)
        plt.subplot(projection=wcs)
        im = plt.imshow(cutout.data,extent=[xl.value,xh.value, yl.value,yh.value])
        plt.scatter(YSO_RA, YSO_DEC, s=5, c='red', marker='o')
        plt.colorbar(im)
        ax.set_aspect('equal')
        plt.show()
    return