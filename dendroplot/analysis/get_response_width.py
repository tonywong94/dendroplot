import numpy as np
from astropy.io import fits
from scipy.stats import pearsonr
from spectral_cube import SpectralCube

def calc_channel_corr(cube, mask=None, channel_lag=1):
    """
    Calculate the channel-to-channel correlation coefficient (Pearson's r)
    Courtesy Jiayi Sun
    https://github.com/astrojysun/SpectralCubeTools/blob/master/spectral_cube_tools/characterize.py

    Parameters
    ----------
    cube : SpectralCube object
        Input spectral cube
    mask : `np.ndarray` object, optional
        User-specified mask, within which to calculate r
    channel_lag : int
        Number of channel lag at which the correlation is calculated.
        Default is 1, which means to estimate correlation at one
        channel lag (i.e., between immediately adjacent channels)

    Returns
    -------
    r : float
        Pearson's correlation coefficient
    p-value : float
        Two-tailed p-value
    """
    if mask is None:
        mask = cube.mask.include()
    mask[-1, :] = False
    for i in np.arange(channel_lag):
        mask &= np.roll(mask, -1, axis=0)
    return pearsonr(
        cube.filled_data[mask],
        cube.filled_data[np.roll(mask, channel_lag, axis=0)])
        

def get_response_width(imcubename, maskname=None, edgech=None, clip=0.5, 
                       logic=True, channel_lag=1):
    """
    Calculate the effective spectral response width following the estimate by
    Leroy et al. (2016ApJ...831...16L).
    
    Parameters
    ----------
    imcubename : file name or SpectralCube object, required
        Input spectral cube (name of FITS file or SpectralCube object)
    maskname : file name or SpectralCube object or `np.ndarray` object, optional
        User-specified mask to exclude from noise calculations.
        Default is to use all non-NaN values in the image cube.
    edgech : int
        If maskname is not given can use this parameter to choose a certain
        number of edge channels on both ends to represent the noise.
    clip : float
        Value taken to separate noise and signal regions in mask image
        Default: 0.5 (suitable for an signal mask cube)
    logic: boolean
        When True, mask image values < clip are considered noise regions.
        When False, mask image values > clip are considered noise regions.
        Default: True (suitable for an assignment cube)
    channel_lag : int
        Number of channel lag at which the correlation is calculated.
        Default is 1, which means to estimate correlation at one
        channel lag (i.e., between immediately adjacent channels)

    Returns
    -------
    resp : float
        Effective spectral response width, in channels
        A value of 1 corresponds to the approximation in Rosolowsky & Leroy (2006).
        This is not divided by sqrt(2*pi), which is done in calc_phys_props.py.

    """
    if isinstance(imcubename, str):
        imcube = SpectralCube.read(imcubename)
    elif isinstance(imcubename, SpectralCube):
        imcube = imcubename
    else:
        raise TypeError('Not a FITS file or SpectralCube:',imcubename)
    if maskname is None:  # use the whole cube or edge channels as noise
        if edgech is not None:
            mask = np.ones_like(imcube._data)
            nchan = imcube._data.shape[0]
            lfree = np.r_[0:edgech,nchan-edgech:nchan]
            mask[lfree,:,:] = 0
        else:
            mask = np.zeros_like(imcube._data)
    elif isinstance(maskname, str):
        mask = fits.getdata(maskname)
    elif isinstance(maskname, SpectralCube):
        mask = maskname.unitless_filled_data
    elif isinstance(maskname, np.ndarray):
        mask = maskname
    else:  
        raise TypeError('Not a FITS file, SpectralCube, or ndarray:',maskname)

    # Apply logic
    if logic:
        noise = (mask < clip) & np.isfinite(imcube)
    else: 
        noise = (mask > clip) & np.isfinite(imcube)
    print('Identified {} of {} pixels as noise'.format(np.count_nonzero(noise),
          np.count_nonzero(np.isfinite(imcube))))

    # Do the calculation
    r = calc_channel_corr(imcube, mask=noise, channel_lag=channel_lag)[0]
    if r > 0.65:
        print('WARNING: Pearson r is very large:',r)
    else:
        print('Calculated Pearson r:',r)
    k = 0.47*r - 0.23*r**2 - 0.16*r**3 + 0.43*r**4
    print('Calculated k value:',k)
    resp = 1 + 1.18*k + 10.4*k**2
    print('Effective response width [channels]:', resp)
    return resp
