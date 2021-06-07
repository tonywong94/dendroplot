from matplotlib import pyplot as plt
from astropy.io.fits import getheader, getdata
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord as SkyCoord
from spectral_cube import SpectralCube
import astropy
import pickle
import astropy.units as u
from astropy.wcs.utils import skycoord_to_pixel
import numpy as np

def yso_asgn_structures(cube_location, YSO_RA, YSO_DEC,name):
    cube,header = (getdata(cube_location,header=True))
    wcs = WCS(header)
    pixel_ra,pixel_dec = skycoord_to_pixel(SkyCoord(YSO_RA,YSO_DEC), wcs)
    pixel_ra = pixel_ra.data.tolist()
    pixel_dec = pixel_dec.data.tolist()
    p1 = [cube[i][int(np.floor(pixel_dec))][int(np.floor(pixel_ra))] for i in range(len(cube))]
    p2 = [cube[i][int(np.ceil(pixel_dec))][int(np.ceil(pixel_ra))] for i in range(len(cube))]
    p3 = [cube[i][int(np.ceil(pixel_dec))][int(np.floor(pixel_ra))] for i in range(len(cube))]
    p4 = [cube[i][int(np.floor(pixel_dec))][int(np.ceil(pixel_ra))] for i in range(len(cube))]
    p1_d = dict(zip(np.arange(len(p1)),p1))
    p2_d = dict(zip(np.arange(len(p2)),p2))
    p3_d = dict(zip(np.arange(len(p3)),p3))
    p4_d = dict(zip(np.arange(len(p4)),p4))
    dicts = p1_d,p2_d,p3_d,p4_d 
    with open(name+".txt", "wb") as f:
        pickle.dump(dicts, f)
    return
def get_yso_asgn(name):
    '''returns tuple of dictionaries of the non-nan values of asgn cube
    
    '''
    with open(name+".txt","rb") as f:
        out = pickle.load(f)
    p1 = {i:int(out[0][i]) for i in range(len(out[0])) if ~np.isnan(out[0][i])}
    p2 = {i:int(out[1][i]) for i in range(len(out[1])) if ~np.isnan(out[1][i])}
    p3 = {i:int(out[2][i]) for i in range(len(out[2])) if ~np.isnan(out[2][i])}
    p4 = {i:int(out[3][i]) for i in range(len(out[3])) if ~np.isnan(out[3][i])}
    return p1,p2,p3,p4
    
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