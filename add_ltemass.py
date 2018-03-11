#!/usr/bin/env python

import numpy as np
from astropy import units as u
from astropy import constants as const
from astrodendro import Dendrogram, ppv_catalog, analysis
from astropy.io.fits import getheader, getdata
from astropy.table import Table, Column

def add_ltemass(label = 'pcc_12', n13cub = None, n13cub_uc = None):

    # Make the uncertainty input a list
    if not isinstance(n13cub_uc, list): n13cub_uc = [n13cub_uc]
# Input parameters
# lines = ['12', '13']
# ltedir = '../lte/'
# n13cube    = ltedir+'PCC_cube_n13cub.fits.gz'
# n13cub_uc1 = ltedir+'PCC_cube_n13cuberr.fits.gz'
# n13cub_uc2 = ltedir+'PCC_noise_25_itrs_rms_n13cub.fits.gz'

    # Adopted parameters
    distpc = 4.8e4 * u.pc
    co13toh2 = 5.0e6  # Indebetouw+ 13

    # Get basic info from header
    hd = getheader(n13cub)
    deltav = hd['cdelt3']/1000.
    pixdeg = hd['cdelt2']
    pix2cm = (np.radians(pixdeg) * distpc).to(u.cm)
    ppbeam = np.abs((hd['bmaj']*hd['bmin'])/(hd['cdelt1']*hd['cdelt2'])
        *2*np.pi/(8*np.log(2)))
    osamp  = np.sqrt(ppbeam)

    # Total the LTE masses
    d = Dendrogram.load_from(label+'_dendrogram.hdf5')
    cat = Table.read(label+'_physprop.txt', format='ascii.ecsv')
    srclist = cat['_idx'].tolist()
    for col in ['mlte', 'e_mlte', 'siglte', 'e_siglte', 'e_mlte_alt']:
        newcol = Column(name=col, data=np.zeros(np.size(srclist)))
        
        if col == 'mlte':
            data = getdata(n13cub)
            newcol.description = 'LTE mass using H2/13CO='+str(co13toh2)
        elif col == 'e_mlte':
            data = getdata(n13cub_uc[0])
            newcol.description = 'fractional unc in mlte'
        elif col == 'siglte':
            data == getdata(n13cub)
            newcol.description = 'LTE mass divided by area in pc2'
        elif col == 'e_siglte':
            data = getdata(n13cub_uc[0])
            newcol.description = 'fractional unc in siglte [same as e_lte]'
        elif col == 'e_mlte_alt':
            if len(n13cub_uc) > 1:
                data = getdata(n13cub_uc[1])
                newcol.description = 'fractional unc in mlte from alt approach'
            else:
                break
        
        for i, c in enumerate(srclist):
            mask = d[c].get_mask()
            if (col == 'mlte' or col == 'siglte'):
                newcol[i] = np.nansum(data[np.where(mask)])
            else:
                newcol[i] = np.sqrt(np.nansum(data[np.where(mask)]**2)) * osamp 

        # Multiply by channel width in km/s and area in cm^2 to get molecule number 
        newcol *= deltav * pix2cm.value**2
        # Convert from molecule number to solar masses including He
        at_to_solmass = 8.411579469870226e-58
        newcol *= co13toh2 * 2 * 1.36 * at_to_solmass

        if col == 'mlte':
            newcol.unit = 'solMass'
        elif col == 'siglte':
            newcol /= cat['area_pc2']
            newcol.unit = 'solMass/pc2'
        else:
            newcol /= cat['mlte']
            newcol.unit = ''

        cat.add_column(newcol)
      
    #cat.pprint(show_unit=True)
    cat.write(label+'_physprop_add.txt', format='ascii.ecsv', overwrite=True)

    return
