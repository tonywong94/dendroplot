import ipyvolume
from astropy.io.fits import getdata
from astropy import units as u
from astropy.wcs import WCS
from spectral_cube import SpectralCube
import numpy as np
def visualize_volume(cube_location,RA,DEC,RA_width=2.0*u.arcmin,DEC_width=0.5*u.arcmin,level_=None):
    cube,header=(getdata(cube_location,header=True))
    w = WCS(header)
    spec_cube = SpectralCube(data=cube,wcs=w)
    xl = min([RA-RA_width, RA+RA_width])
    xh = max([RA-RA_width, RA+RA_width])
    yl = min([DEC-DEC_width, DEC+DEC_width])
    yh = max([DEC-DEC_width, DEC+DEC_width])
    zl = spec_cube.spectral_axis[0]
    zh = spec_cube.spectral_axis[-1]
    cutout = spec_cube.subcube(xlo=xl,xhi=xh,ylo=yl,yhi=yh,zlo = zl, zhi = zh)
    cutout_cube = cutout.hdu.data
    if level_ is None:
        level_ = np.median(cutout_cube[~np.isnan(cutout_cube)])
    ipyvolume.pylab.plot_isosurface(cutout_cube, level=level_)
    ipyvolume.pylab.xlabel('VRAD')
    ipyvolume.pylab.ylabel('DEC')
    ipyvolume.pylab.zlabel('RA')
    ipyvolume.pylab.gcf()
    ipyvolume.pylab.show()