import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
from astropy.io import fits
from astropy.visualization.wcsaxes import SphericalCircle
import numpy as np
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle

# Plot the 12CO and 13CO moment maps side by side.

skel = False

# Get FOV from rms image
fitsdir = 'images/'
fildir  = 'struct/'
filer = fitsdir+'30Dor_feather_mosaic_12CO_12meter.rms.K.fits.gz'
rms = fits.getdata(filer)[0]
fov = rms*0 + 1
fov[np.isnan(fov)] = 0

for plttyp in ['snrpk', 'mom01', 'mom2']:
    if plttyp == 'snrpk':
        file0 = fitsdir+'30Dor_feather_mosaic_12CO_12meter.snrpk.fits.gz'
        file1 = fitsdir+'30Dor_feather_mosaic_13CO_12meter.snrpk.fits.gz'
        label = ['$^{12}$CO (2-1) peak SNR', '$^{13}$CO (2-1) peak SNR']
        cmap  = ['gist_ncar_r', 'gist_ncar_r']
        outfile = '30Dor_12CO_13CO_snrpk.pdf'
    elif plttyp == 'mom01':
        file0 = fitsdir+'30Dor_feather_mosaic_12CO_12meter.mom0.fits.gz'
        file1 = fitsdir+'30Dor_feather_mosaic_12CO_12meter.mom1.fits.gz'
        label = ['$^{12}$CO intensity', '$^{12}$CO velocity']
        cmap  = ['gist_ncar_r', 'RdBu_r']
        outfile = '30Dor_12CO_mom0_mom1.pdf'
    elif plttyp == 'mom2':
        file0 = fitsdir+'30Dor_feather_mosaic_12CO_12meter.mom2.fits.gz'
        file1 = fitsdir+'30Dor_feather_mosaic_13CO_12meter.mom2.fits.gz'
        label = ['$^{12}$CO mom-2', '$^{13}$CO mom-2']
        cmap  = ['jet_r', 'jet_r']
        outfile = '30Dor_12CO_13CO_mom2.pdf'

    im0, hd0 = fits.getdata(file0, header=True)
    im1, hd1 = fits.getdata(file1, header=True)
    wcs = WCS(hd0)
    olays = np.array(wcs.world_to_pixel(SkyCoord(ra=[84.67625, 84.5708]*u.deg, 
                                        dec=[-69.100917, -69.0667]*u.deg)))

    if plttyp == 'mom01':
        fig = plt.figure(figsize=(18,10))
    else:
        fig = plt.figure(figsize=(18,10))

    for i, image in enumerate([im0, im1]):
        ax = plt.subplot(1, 2, i+1, projection=wcs)
        if plttyp == 'snrpk':
            if i == 0:
                norm = mcolors.Normalize(vmin=3, vmax=90)
            else:
                norm = mcolors.Normalize(vmin=3, vmax=30)
        elif plttyp == 'mom01':
            if i == 0:
                norm = mcolors.Normalize(vmin=3, vmax=120)
            else:
                norm = mcolors.TwoSlopeNorm(vmin=235, vcenter=255, vmax=275)
        elif plttyp == 'mom2' or plttyp == 'lwidth':
            norm = mcolors.PowerNorm(gamma=0.5, vmin=0, vmax=9)
        im = ax.imshow(image, origin='lower', cmap=cmap[i], norm=norm)
        print('R136 is at',olays[0][0], olays[1][0])
        print('Hodge 301 is at',olays[0][1], olays[1][1])
        ax.plot(olays[0], olays[1], 'b*', ms=12)
        ax.annotate('R136', (olays[0][0], olays[1][0]), textcoords="offset points",
                     xytext=(-5,-20), ha='right', c='b', size=13)
        ax.annotate('Hodge 301', (olays[0][1], olays[1][1]), textcoords="offset points",
                     xytext=(-5,-20), ha='center', c='b', size=13)
        if plttyp == 'snrpk':
            r136circ = SphericalCircle((84.67625 * u.deg, -69.100917 * u.deg), 
                             200*u.arcsec, edgecolor='b', facecolor='none', 
                             ls='--', transform=ax.get_transform('fk5'))
            ax.add_patch(r136circ)
            from astropy.coordinates import SkyCoord
            coord = SkyCoord('05h38m47.8s', '-69d04m50.3s', frame='fk5')
            pixels = wcs.world_to_pixel(coord)  
            print(pixels)
            cycle0 = Rectangle((378.5-70,473.298), 120, 108, angle=-40,
                               edgecolor='k', facecolor='none', linestyle='--')
            ax.add_patch(cycle0)
        else:
            ax.contour(fov, [1], colors='k', linestyles='dotted')
        if (plttyp == 'mom01') and skel:
            skelfile = fits.getdata(fildir+'30Dor_feather_12CO_fp_80_st_80_at_10_gt_4.skel.fits.gz')
            ax.contour(skelfile, levels=[1], colors='k', linewidths=1)
        im_ratio = image.shape[0]/image.shape[1]
        print('image ratio is',im_ratio)
        cb = fig.colorbar(im, ax=ax, orientation='vertical', pad=0.04, 
                          fraction=0.045*im_ratio)
        cb.ax.tick_params(labelsize=11)
        ax.coords[0].set_ticklabel(size=14)
        ax.coords[1].set_ticklabel(size=14)
        if plttyp == 'mom01':
            plt.subplots_adjust(wspace=0.1)
            ax.coords[1].set_ticklabel_visible(False)

        ax.set_xlim(-0.5,799.5)
        ax.set_ylim(-0.5,999.5)
        ax.set_xlabel('Right Ascension (J2000)', size=15)
        ax.text(0.03,0.95,label[i],ha='left',va='center',size=18,transform=ax.transAxes)
        if i == 0:
            ax.set_ylabel('Declination (J2000)', size=15)
        else:
            ax.set_ylabel(' ')
        ax.tick_params(direction='in')

    plt.savefig(outfile, bbox_inches='tight')
