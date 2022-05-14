#!/usr/bin/env python

import numpy as np
from astropy import units as u
from astropy import constants as const
from astrodendro import Dendrogram
from astropy.io.fits import getheader, getdata
from astropy.table import Table, Column

'''
PURPOSE: Add columns to physprop.txt table based on LTE analysis.
    Required keywords:
        label: prefix for table, e.g. 'pcc_12' for pcc_12_physprop.txt
        n13cub: name of 13CO column density cube FITS image
        n13cub_uc: 13CO column density uncertainty cube, can be a list of up to 2 cubes
    Optional keywords:
        i12cub: name of 12CO intensity cube FITS image (pb corrected)
        i13cub: name of 13CO intensity cube FITS image (pb corrected)
        distpc: distance in pc for mass calculation, defaults to 50000 (LMC)
        co13toh2: H2/13CO abundance ratio, defaults to 5.0e6 (Indebeouw+ 13)
'''

def add_ltemass(label = 'pcc_12', n13cub = None, i12cub = None, i13cub = None,
                efloor=0, n13cub_uc = None, distpc = 5e4, co13toh2 = 5.0e6):

    # Make the uncertainty input a list
    if not isinstance(n13cub_uc, list): n13cub_uc = [n13cub_uc]

    # Adopted parameters
    dist = distpc * u.pc

    # Get basic info from header
    hd = getheader(n13cub)
    deltav = np.abs(hd['cdelt3']/1000.)
    pixdeg = np.abs(hd['cdelt2'])
    pix2cm = (np.radians(pixdeg) * dist).to(u.cm)
    ppbeam = np.abs((hd['bmaj']*hd['bmin'])/(hd['cdelt1']*hd['cdelt2'])
        *2*np.pi/(8*np.log(2)))
    osamp  = np.sqrt(ppbeam)

    # Total the LTE masses (and optionally, 12CO and 13CO fluxes)
    d = Dendrogram.load_from(label+'_dendrogram.hdf5')
    print('\n')
    pcat = Table.read(label+'_physprop.txt', format='ascii.ecsv')
    srclist = pcat['_idx'].tolist()
    datcol = np.zeros(np.size(srclist))

    # Note that the order in which the columns are processed is important!
    for col in ['flux12', 'flux13', 'mlte', 'siglte', 'e_mlte', 
                'e_siglte', 'e_mlte_alt']:
        print('Extracting data for column', col)
        newcol = Column(name=col, data=np.zeros(np.size(srclist)))
        
        if col == 'flux12':
            if i12cub is not None:
                data, ihd = getdata(i12cub, header=True)
                if 'RESTFREQ' in ihd.keys():
                    rfreq = ihd['RESTFREQ'] * u.Hz
                elif 'RESTFRQ' in ihd.keys():
                    rfreq = ihd['RESTFRQ'] * u.Hz
                newcol.description = '12CO flux within the structure'
            else:
                continue
        elif col == 'flux13':
            if i13cub is not None:
                data, ihd = getdata(i13cub, header=True)
                if 'RESTFREQ' in ihd.keys():
                    rfreq = ihd['RESTFREQ'] * u.Hz
                elif 'RESTFRQ' in ihd.keys():
                    rfreq = ihd['RESTFRQ'] * u.Hz
                newcol.description = '13CO flux within the structure'
            else:
                continue
        elif col == 'mlte':
            data = getdata(n13cub)
            newcol.description = 'LTE mass using H2/13CO='+str(co13toh2)
            newcol.unit = 'solMass'
        elif col == 'siglte':
            newcol.description = 'LTE mass divided by area in pc2'
        elif col == 'e_mlte':
            data = getdata(n13cub_uc[0])
            newcol.description = 'fractional unc in mlte'
        elif col == 'e_siglte':
            newcol.description = 'fractional unc in siglte [same as e_mlte]'
        elif col == 'e_mlte_alt':
            if len(n13cub_uc) > 1:
                data = getdata(n13cub_uc[1])
                newcol.description = 'fractional unc in mlte from alt approach'
            else:
                continue
        
        # We use previously read data from mlte and e_mlte for siglte and e_siglte
        if 'flux' in col or 'mlte' in col:
            for i, c in enumerate(srclist):
                mask = d[c].get_mask()
                if not col.startswith('e_'):
                    datcol[i] = np.nansum(data[np.where(mask)])
                    # nansum returns zero if all are NaN, want NaN
                    chknan = np.asarray(np.isnan(data[np.where(mask)]))
                    if chknan.all():
                        datcol[i] = np.nan
                else:
                    datcol[i] = np.sqrt(np.nansum(data[np.where(mask)]**2)) * osamp 

        if col in ['flux12', 'flux13']:
            # Convert from K*pix*ch to Jy*km/s
            convfac = (1*u.K).to(u.Jy/u.deg**2, equivalencies=u.brightness_temperature(rfreq))
            newcol.data[:] = datcol * deltav * convfac.value * (pixdeg)**2
            newcol.unit = 'Jy km / s'
        else:
            # Multiply by channel width in km/s and area in cm^2 to get molecule number 
            newcol.data[:] = datcol * deltav * pix2cm.value**2
            # Convert from molecule number to solar masses including He
            newcol *= co13toh2 * 2 * 1.36 * const.m_p.value / const.M_sun.value
            if col == 'siglte':
                newcol /= pcat['area_pc2']
                newcol.unit = 'solMass/pc2'
            elif col == 'e_siglte' or 'e_mlte' in col:
                newcol /= np.abs(pcat['mlte'])
                #newcol.unit = ''

        # ---- apply a floor to the fractional uncertainty if requested
        if col.startswith('e_') and col.endswith('lte') and efloor > 0:
            print( "Applying a minimum fractional error of {:2.3f}".format(efloor) )
            newcol[newcol < efloor] = efloor

        pcat.add_column(newcol)
      
    #pcat.pprint(show_unit=True)
    pcat['alphalte'] = pcat['mvir'] / pcat['mlte']
    pcat['alphalte'].description = 'virial parameter from mvir and mlte'
    pcat['alphalte'].unit = ''
    pcat['e_alphalte'] = np.sqrt(pcat['e_mlte']**2 + pcat['e_mvir']**2)
    pcat['e_alphalte'].description = 'fractional error in alphalte'
    pcat.write(label+'_physprop_add.txt', format='ascii.ecsv', overwrite=True)

    return
