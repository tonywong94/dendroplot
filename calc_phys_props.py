#!/usr/bin/env python

import os
import numpy as np
from astropy import units as u
from astropy import constants as const
from astrodendro import Dendrogram, ppv_catalog
from astrodendro.analysis import ScalarStatistic, PPVStatistic
from astropy.table import Table, Column
from astropy.io.fits import getdata
from astropy.stats import mad_std
from astropy.wcs import WCS
from matplotlib import pyplot as plt

'''
PURPOSE: Create the physprop.txt table from the dendrogram catalog.
    Required keywords:
        label: prefix for dendrogram & catalog, e.g. 'pcc_12' for pcc_12_dendrogram.hdf5
        cubefile: FITS cube from which dendrogram was derived
    Optional keywords:
        boot_iter: number of bootstrap iterations, defaults to 400
        efloor: minimum fractional error to assume, defaults to 0
        alphascale: scaling of \alpha_CO vs. Galactic, defaults to 1
        distpc: distance in pc for mass calculation, defaults to 50000 
            (LMC)
        copbcor: primary-beam corrected 12CO cube to measure XCO-based mass
        conoise: RMS noise cube or map corresponding to copbcor, required if
            copbcor is given
        ancfile: another wavelength image (e.g. 8um) in which to calculate
            mean brightness of each structure.  Can be 2D image or cube, but
            either way should be already regridded to match 'cubefile'
        anclabel: label for column corresponding to ancfile (e.g. '8um_avg')
'''

def clustbootstrap(sindices, svalues, meta, bootstrap):
    bootiters = range(bootstrap)
    emmajs = []; emmins = []; emomvs = []; eflux = []; pa = []
    npts = len(svalues)
    for bootit in bootiters:
        # Shuffle the indices and values as in CPROPS...
        randidx = np.random.random(npts)*npts
        randidx = randidx.astype('int')
        bvalues = svalues[randidx]
        bindices = sindices[0][randidx],sindices[1][randidx],sindices[2][randidx]
        # ... and generate new statistics
        bscalstats = ScalarStatistic(bvalues,bindices)
        bstats = PPVStatistic(bscalstats, meta)
        emmajs.append(bstats.major_sigma.value)
        emmins.append(bstats.minor_sigma.value)
        emomvs.append(bstats.v_rms.value)
        eflux.append(bstats.flux.value)
        pa.append(bstats.position_angle.value)
    return emmajs, emmins, emomvs, eflux, pa


def calc_phys_props(label='pcc_12', cubefile=None, boot_iter=400, efloor=0,
        alphascale=1, distpc=5e4, copbcor=None, conoise=None, ancfile=None, 
        anclabel=None, refpos=None):

    rmstorad= 1.91
    alphaco = 4.3 * u.solMass * u.s / (u.K * u.km * u.pc**2) # Bolatto+ 13
    dist    = distpc * u.pc
    as2     = 1 * u.arcsec**2
    asarea  = (as2*dist**2).to(u.pc**2,equivalencies=u.dimensionless_angles())

    # ---- load the dendrogram and catalog
    d = Dendrogram.load_from(label+'_dendrogram.hdf5')
    cat = Table.read(label+'_full_catalog.txt', format='ascii.ecsv')
    srclist = cat['_idx'].tolist()

    # ---- load the cube and extract the metadata
    cube, hd3 = getdata(cubefile, header=True)
    metadata = {}
    if hd3['BUNIT'].upper()=='JY/BEAM':
        metadata['data_unit'] = u.Jy / u.beam
    elif hd3['BUNIT'].upper()=='K':
        metadata['data_unit'] = u.K
    else:
        print("\nWarning: Unrecognized brightness unit")
    metadata['vaxis'] = 0
    if 'RESTFREQ' in hd3.keys():
        freq = hd3['RESTFREQ'] * u.Hz
    elif 'RESTFRQ' in hd3.keys():
        freq = hd3['RESTFRQ'] * u.Hz
    cdelt1 = abs(hd3['cdelt1']) * 3600. * u.arcsec
    cdelt2 = abs(hd3['cdelt2']) * 3600. * u.arcsec
    deltav = abs(hd3['cdelt3'])/1000. * u.km / u.s
    metadata['wavelength'] = freq.to(u.m,equivalencies=u.spectral())
    metadata['spatial_scale']  =  cdelt2
    metadata['velocity_scale'] =  deltav
    bmaj = hd3['bmaj']*3600. * u.arcsec # FWHM
    bmin = hd3['bmin']*3600. * u.arcsec # FWHM
    metadata['beam_major'] = bmaj
    metadata['beam_minor'] = bmin
    ppbeam = np.abs((bmaj*bmin)/(cdelt1*cdelt2)*2*np.pi/(8*np.log(2)))
    print("\nPixels per beam: {:.2f}".format(ppbeam))
    # Assume every 2 channels have correlated noise
    indfac = np.sqrt(ppbeam*2)
    if refpos is not None:
        pixcrd = np.array([[hd3['CRPIX1'], hd3['CRPIX2'], 1]])
        w = WCS(hd3)
        # Debug: only works with hd3 having 3 axes
        wcscrd = w.wcs_pix2world(pixcrd, 1)
        refpix = w.wcs_world2pix([[refpos[0],refpos[1],wcscrd[0][2]]],0)
        xref = refpix[0][0]
        yref = refpix[0][1]
        print('Ref pixel is',xref,yref)

    # ---- read in ancillary files
    if copbcor is not None and conoise is not None:
        cube12, hd12 = getdata(copbcor, header=True)
        ecube12, ehd12 = getdata(conoise, header=True)
        if 'RESTFREQ' in hd12.keys():
            freq12 = hd12['RESTFREQ'] * u.Hz
        elif 'RESTFRQ' in hd12.keys():
            freq12 = hd12['RESTFRQ'] * u.Hz
        else:
            print('Rest frequency missing from file '+copbcor)
            raise
        if hd12['BUNIT'].upper() != 'K':
            print('Non-Kelvin units not yet supported for copbcor')
            raise
        if ehd12['NAXIS'] == 2:
            tmpcube = np.broadcast_to(ecube12, np.shape(cube12))
            ecube12 = tmpcube
    if ancfile is not None:
        ancdata,anchd = getdata(ancfile, header=True)

    # ---- call the bootstrapping routine
    emaj, emin, epa, evrms, errms, eaxra, eflux, emvir, ealpha, tb12, ancmean, ancrms = [
        np.zeros(len(srclist)) for _ in range(12)]
    print("Calculating property errors...")
    for j, clust in enumerate(srclist):
        if j % 10 == 1:
            print("Calculating property errors for structure",j)
        asgn = np.zeros(cube.shape)
        asgn[d[clust].get_mask(shape = asgn.shape)] = 1
        sindices = np.where(asgn == 1)
        svalues = cube[sindices]

        emmajs, emmins, emomvs, emom0s, pa = clustbootstrap(
            sindices, svalues, metadata, boot_iter)
        # bin_list=np.linspace(np.floor(min(emmajs)),np.ceil(max(emmajs)),50)
        # fig, axes = plt.subplots()
        # axes.hist(emmajs, bin_list, normed=0, histtype='bar')
        # plt.savefig('emmajs_histo.pdf', bbox_inches='tight')

        bootmaj   = np.asarray(emmajs) * u.arcsec
        bootmin   = np.asarray(emmins) * u.arcsec
        bootpa    = np.asarray(pa) * u.deg
        bootvrms  = (np.asarray(emomvs)*u.m/u.s).to(u.km/u.s)
        bootrrms  = (np.sqrt(np.asarray(emmajs)*np.asarray(emmins))*u.arcsec*dist).to(
            u.pc, equivalencies=u.dimensionless_angles())
        bootaxrat = bootmin/bootmaj
        bootflux  = np.asarray(emom0s) * u.Jy
        bootmvir  = (5*rmstorad*bootvrms**2*bootrrms/const.G).to(u.solMass)
        bootmlum  = alphaco*alphascale*deltav*asarea*(bootflux).to(
            u.K, equivalencies=u.brightness_temperature(freq, beam_area=as2))
        bootalpha = bootmvir/bootmlum

        emaj[j]   = indfac * mad_std(bootmaj)  / np.median(bootmaj)
        emin[j]   = indfac * mad_std(bootmin)  / np.median(bootmin)
        epa[j]    = indfac * mad_std(bootpa)   / np.median(bootpa)
        evrms[j]  = indfac * mad_std(bootvrms) / np.median(bootvrms)
        errms[j]  = indfac * mad_std(bootrrms) / np.median(bootrrms)    
        eaxra[j]  = indfac * mad_std(bootaxrat)/ np.median(bootaxrat)
        emvir[j]  = indfac * mad_std(bootmvir) / np.median(bootmvir)
        ealpha[j] = indfac * mad_std(bootalpha)/ np.median(bootalpha)

        if copbcor is not None and conoise is not None:
            tb12[j]  = np.nansum(asgn * cube12)
            eflux[j] = indfac * np.sqrt(np.nansum(asgn * ecube12**2)) / tb12[j]
            #print(j, len(svalues), tb12[j], etb12[j])
        else:
            eflux[j]  = indfac * mad_std(bootflux) / np.median(bootflux)
        if ancfile is not None:
            if anchd['NAXIS'] == 2:
                collapsedmask = np.amax(asgn, axis = 0)
                collapsedmask[collapsedmask==0] = np.nan
                ngood = np.nansum(collapsedmask)
                #print('Pixels in region:',ngood)
                ancmean[j] = np.nanmean(ancdata*collapsedmask)
                ancrms[j] = np.sqrt(np.nanmean((ancdata*collapsedmask)**2)-ancmean[j]**2)
                ancrms[j] /= np.sqrt(ngood/ppbeam)
                #print('Mean and uncertaintity:',ancmean[j],ancrms[j])
            else:
                ancmean[j] = np.nanmean(ancdata*asgn)
                ancrms[j] = np.sqrt(np.nanmean((ancdata*asgn)**2)-ancmean[j]**2)

    # ---- report the median uncertainties
    print( "The median fractional error in rad_pc is {:2.4f}"
        .format(np.nanmedian(errms)) )
    print( "The median fractional error in vrms_k is {:2.4f}"
        .format(np.nanmedian(evrms)) )
    print( "The median fractional error in mlumco is {:2.4f}"
        .format(np.nanmedian(eflux)) )
    
    # ---- apply a floor if requested
    if efloor > 0:
        print( "Applying a minimum fractional error of {:2.3f}".format(efloor) )
        errms[errms<efloor] = efloor
        evrms[evrms<efloor] = efloor
        eflux[eflux<efloor] = efloor

    # ---- calculate the physical properties
    rms_pc  = (cat['radius'] * dist).to(
        u.pc, equivalencies=u.dimensionless_angles())
    rad_pc  = rmstorad * rms_pc
    v_rms   = cat['v_rms'].to(u.km/u.s)
    axrat   = cat['minor_sigma'] / cat['major_sigma']
    ellarea = (cat['area_ellipse']*dist**2).to(
        u.pc**2,equivalencies=u.dimensionless_angles())
    xctarea = (cat['area_exact']*dist**2).to(
        u.pc**2,equivalencies=u.dimensionless_angles())
    if copbcor is not None and conoise is not None:
        # Convert from K*pix*ch to Jy*km/s
        #convfac = (1*u.K).to(u.Jy/u.arcsec**2, equivalencies=u.brightness_temperature(freq12))
        #flux12 = tb12 * deltav.value * convfac.value * cdelt2.value**2
        mlumco = alphaco * alphascale * tb12*u.K * deltav * asarea * cdelt2.value**2
    else:
        # lumco = Luminosity in K km/s pc^2
        lumco  = deltav * asarea * (cat['flux']).to(
            u.K,equivalencies=u.brightness_temperature(freq, beam_area=as2))
        mlumco = alphaco * alphascale * lumco
    siglum = mlumco/xctarea
    mvir   = (5*rmstorad*v_rms**2*rms_pc/const.G).to(u.solMass)  # Rosolowsky+ 08
    sigvir = mvir / xctarea
    alpha  = mvir / mlumco

    # ---- calculate the projected distance from the reference position in arcsec
    if refpos is not None:
        refdist = abs(cdelt2/u.pix) * np.sqrt(
                  (xref-cat['x_cen'])**2 + (yref-cat['y_cen'])**2 )
    
    # ---- make the physical properties table
    ptab = Table()
    ptab['_idx']      = Column(srclist)
    ptab['area_pc2']  = Column(xctarea, description='projected area of structure')
    ptab['rad_pc']    = Column(rad_pc, description='equivalent radius in pc')
    ptab['e_rad_pc']  = Column(errms, description='frac error in radius')
    ptab['vrms_k']    = Column(v_rms, description='rms linewidth in km/s')
    ptab['e_vrms_k']  = Column(evrms, description='frac error in linewidth')
    ptab['axratio']   = Column(axrat, unit='', description='minor to major axis ratio')
    ptab['e_axratio'] = Column(eaxra, description='frac error in axis ratio')
    if copbcor is not None and conoise is not None:
        #ptab['flux12']    = Column(flux12, unit='Jy km / s', description='CO flux in structure')
        #ptab['e_flux12']  = Column(eflux, description='frac error in CO flux')
        ptab['mlumco']    = Column(mlumco, description='CO-based mass with alphascale='+str(alphascale)+' from '+os.path.basename(copbcor))
    else:
        ptab['mlumco']    = Column(mlumco, description='Mass from scaling luminosity with alphascale='+str(alphascale))
    ptab['e_mlumco']  = Column(eflux, description='frac error in luminous mass')
    ptab['siglum']    = Column(siglum, description='average surface density from mlumco')
    ptab['e_siglum']  = Column(eflux, description='same as e_mlumco')
    ptab['mvir']      = Column(mvir, description='virial mass')
    ptab['e_mvir']    = Column(emvir, description='frac error in virial mass')
    ptab['sigvir']    = Column(sigvir, description='virial surface density')
    ptab['e_sigvir']  = Column(emvir, description='same as e_mvir')
    ptab['alpha']     = Column(alpha, unit='', description='virial parameter')
    ptab['e_alpha']   = Column(ealpha, description='frac error in virial parameter')
    if ancfile is not None:
        if anclabel is None:
            anclabel = ancimg.replace('.','_').split('_')[1]
        ptab[anclabel] = Column(ancmean, unit=anchd['BUNIT'])
        ancferr = ancrms / ancmean
        ptab['e_'+anclabel] = Column(ancferr)
    if refpos is not None:
        ptab['refdist'] = Column(refdist, description='distance from '+str(refpos))
    ptab.write(label+'_physprop.txt', format='ascii.ecsv', overwrite=True)


def refdist_redo(label='pcc_12', cubefile=None, refpos=None, 
                 outfile='_physprop_newref.txt'):
    cube, hd3 = getdata(cubefile, header=True)
    cat = Table.read(label+'_full_catalog.txt', format='ascii.ecsv')
    ptab = Table.read(label+'_physprop.txt', format='ascii.ecsv')
    if refpos is not None:
        cdelt2 = abs(hd3['CDELT2']) * 3600. * u.arcsec
        pixcrd = np.array([[hd3['CRPIX1'], hd3['CRPIX2'], 1]])
        w = WCS(hd3)
        # Debug: only works with hd3 having 3 axes
        wcscrd = w.wcs_pix2world(pixcrd, 1)
        refpix = w.wcs_world2pix([[refpos[0],refpos[1],wcscrd[0][2]]],0)
        xref = refpix[0][0]
        yref = refpix[0][1]
        print('Ref pixel is',xref,yref)
        refdist = abs(cdelt2/u.pix) * np.sqrt(
                  (xref-cat['x_cen'])**2 + (yref-cat['y_cen'])**2 )
        ptab['refdist'] = Column(refdist, description='distance from '+str(refpos))
        ptab.write(label+outfile, format='ascii.ecsv', overwrite=True)
