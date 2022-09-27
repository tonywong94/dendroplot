#!/usr/bin/env python

import os
import numpy as np
from astropy import units as u
from astropy import constants as const
from astrodendro import Dendrogram, ppv_catalog
from astrodendro.analysis import ScalarStatistic, PPVStatistic
from astropy.table import QTable, Table, Column
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

def deconvolve_gauss(meas_maj, meas_min, meas_pa, beam_maj, beam_min, beam_pa):
    # output result of deconvolving one Gaussian beam from another
    # ADAPTED FROM gaupar.for in MIRIAD via K. Sandstrom
    theta_in = meas_pa.to_value(u.rad)
    theta_bm = beam_pa.to_value(u.rad)
    alpha = ((meas_maj*np.cos(theta_in))**2 + (meas_min*np.sin(theta_in))**2 - 
             (beam_maj*np.cos(theta_bm))**2 - (beam_min*np.sin(theta_bm))**2)
    beta  = ((meas_maj*np.sin(theta_in))**2 + (meas_min*np.cos(theta_in))**2 - 
             (beam_maj*np.sin(theta_bm))**2 - (beam_min*np.cos(theta_bm))**2)
    gamma = 2.*((meas_min**2-meas_maj**2)*np.sin(theta_in)*np.cos(theta_in) 
              - (beam_min**2-beam_maj**2)*np.sin(theta_bm)*np.cos(theta_bm))
    s = alpha + beta
    t = np.sqrt((alpha-beta)**2 + gamma**2)
    src_maj = np.full_like(meas_maj, np.nan)
    src_min = np.full_like(meas_min, np.nan)
    src_pa  = np.full_like(theta_in, np.nan)
    good = (alpha > 0) & (beta > 0) & (s > t)
    src_maj[good] = np.sqrt(0.5*(s[good]+t[good]))
    src_min[good] = np.sqrt(0.5*(s[good]-t[good]))
    src_pa[good]  = 0.5 * np.arctan2(-1*gamma[good],alpha[good]-beta[good])
    return src_maj, src_min, src_pa


def clustbootstrap(sindices, svalues, meta, bootiter):
    # recalculate structure properties using bootstrap resampling
    # Adapted from code by D. Colombo
    emmajs, emmins, emomvs, emom0s, pas = np.zeros((5, bootiter))
    npts = len(svalues)
    for i in range(bootiter):
        # Shuffle the indices and values as in CPROPS...
        randidx = np.random.random(npts)*npts
        randidx = randidx.astype('int')
        bvalues = svalues[randidx]
        bindices = sindices[0][randidx],sindices[1][randidx],sindices[2][randidx]
        # ... and generate new statistics
        bscalstats = ScalarStatistic(bvalues, bindices)
        bstats = PPVStatistic(bscalstats, meta)
        emmajs[i] = bstats.major_sigma.value
        emmins[i] = bstats.minor_sigma.value
        emomvs[i] = bstats.v_rms.value
        emom0s[i] = bstats.flux.value
        pas[i]    = bstats.position_angle.value
    return emmajs, emmins, emomvs, emom0s, pas


def calc_phys_props(label='pcc_12', cubefile=None, boot_iter=100, efloor=0,
        alphascale=1, distpc=5e4, rmstorad=1.91, copbcor=None, conoise=None, 
        ancfile=None, anclabel=None, refpos=None, dv_resp=1):

    fwhmfac = np.sqrt(8*np.log(2))
    alphaco = 4.3 * u.solMass * u.s / (u.K * u.km * u.pc**2) # Bolatto+ 13
    alphaco *= alphascale
    dist    = distpc * u.pc
    as2     = 1 * u.arcsec**2
    asarea  = (as2*dist**2).to(u.pc**2,equivalencies=u.dimensionless_angles())

    # ---- load the dendrogram and catalog
    d = Dendrogram.load_from(label+'_dendrogram.hdf5')
    cat = QTable.read(label+'_full_catalog.txt', format='ascii.ecsv')
    cat.add_index('_idx')
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
    deltav = abs(hd3['cdelt3']) / 1000. * u.km / u.s
    sigvchan = deltav/np.sqrt(2*np.pi)
    metadata['wavelength'] = freq.to(u.m,equivalencies=u.spectral())
    metadata['spatial_scale']  =  cdelt2
    metadata['velocity_scale'] =  deltav
    bmaj = hd3['bmaj']*3600. * u.arcsec # FWHM
    bmin = hd3['bmin']*3600. * u.arcsec # FWHM
    bpa  = hd3['bpa'] * u.deg
    metadata['beam_major'] = bmaj
    metadata['beam_minor'] = bmin
    metadata['beam_pa'] = bpa
    beam_maj = bmaj/fwhmfac  # sigma
    beam_min = bmin/fwhmfac  # sigma
    ppbeam = 2*np.pi*np.abs((beam_maj*beam_min)/(cdelt1*cdelt2))
    print("\nPixels per beam: {:.2f}".format(ppbeam))

    # Assume every 2 channels have correlated noise
    indfac = np.sqrt(ppbeam * 2)
    if refpos is not None:
        pixcrd = np.array([[hd3['CRPIX1'], hd3['CRPIX2'], 1]])
        w = WCS(hd3)
        # Debug: only works with hd3 having 3 axes
        wcscrd = w.wcs_pix2world(pixcrd, 1)
        refpix = w.wcs_world2pix([[refpos[0],refpos[1],wcscrd[0][2]]],0)
        xref = refpix[0][0]
        yref = refpix[0][1]
        print('Ref position is at pixel',xref,yref)

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

    # ---- estimate property errors and flux in ancillary image
    emaj, emin, epa, evrms, evrms_d, errms, errms_d, eaxra, eflux, \
          tb12, ancmean, ancrms = [np.zeros(len(srclist)) for _ in range(12)]
    print("Calculating property errors...")
    for j, clust in enumerate(srclist):
        if j % 10 == 1:
            print("Calculating property errors for structure",j)
        asgn = np.zeros(cube.shape)
        asgn[d[clust].get_mask(shape = asgn.shape)] = 1
        sindices = np.where(asgn == 1)
        svalues  = cube[sindices]
        # call the bootstrapping routine        
        emmajs, emmins, emomvs, emom0s, pas = clustbootstrap(
            sindices, svalues, metadata, boot_iter)
        # bin_list=np.linspace(np.floor(min(emmajs)),np.ceil(max(emmajs)),50)
        # fig, axes = plt.subplots()
        # axes.hist(emmajs, bin_list, normed=0, histtype='bar')
        # plt.savefig('emmajs_histo.pdf', bbox_inches='tight')
        bootmaj    = emmajs * u.arcsec
        bootmin    = emmins * u.arcsec
        bootflux   = emom0s * u.Jy
        bootpa     = pas * u.deg
        bootvrms   = (emomvs * u.km/u.s)
        # Deconvolve channel width
        valid = (bootvrms > dv_resp*sigvchan)
        bootvrms_d = np.full_like(bootvrms, np.nan)
        bootvrms_d[valid] = np.sqrt(bootvrms[valid]**2-(dv_resp*sigvchan)**2)
        # Deconvolve beam from bootstrap results
        bootmaj_d, bootmin_d, bootpa_d = deconvolve_gauss(bootmaj, bootmin, bootpa, 
                                                          beam_maj, beam_min, bpa)
        bootaxrat  = bootmin / bootmaj
        bootrrms   = (np.sqrt(bootmaj*bootmin)*dist).to(u.pc, 
                      equivalencies=u.dimensionless_angles())
        bootrrms_d = (np.sqrt(bootmaj_d*bootmin_d)*dist).to(u.pc, 
                      equivalencies=u.dimensionless_angles())
        bootmlum   = alphaco * deltav * asarea * (bootflux).to(
                     u.K, equivalencies=u.brightness_temperature(freq, beam_area=as2))
#         bootmvir   = (5*rmstorad*bootvrms**2*bootrrms/const.G).to(u.solMass)
#         bootmvir_d = (5*rmstorad*bootvrms**2*bootrrms_d/const.G).to(u.solMass)
#         bootalpha  = bootmvir_d / bootmlum

        emaj[j]    = indfac * mad_std(bootmaj)    / np.median(bootmaj)
        emin[j]    = indfac * mad_std(bootmin)    / np.median(bootmin)
        epa[j]     = indfac * mad_std(bootpa)     / np.median(bootpa)
        evrms[j]   = indfac * mad_std(bootvrms)   / np.median(bootvrms)
        evrms_d[j] = indfac * mad_std(bootvrms_d, ignore_nan=True) / np.nanmedian(bootvrms_d)
        errms[j]   = indfac * mad_std(bootrrms)   / np.median(bootrrms)    
        errms_d[j] = indfac * mad_std(bootrrms_d, ignore_nan=True) / np.nanmedian(bootrrms_d)    
        eaxra[j]   = indfac * mad_std(bootaxrat)  / np.median(bootaxrat)
#         emvir[j]   = indfac * mad_std(bootmvir)   / np.median(bootmvir)
#         emvir_d[j] = indfac * mad_std(bootmvir_d, ignore_nan=True) / np.nanmedian(bootmvir_d)
#         ealpha[j]  = indfac * mad_std(bootalpha)  / np.median(bootalpha)

        # Use CO noise cube for flux error estimation if available
        if copbcor is not None and conoise is not None:
            tb12[j]  = np.nansum(asgn * cube12)
            eflux[j] = indfac * np.sqrt(np.nansum(asgn * ecube12**2)) / tb12[j]
        else:
            eflux[j] = indfac * mad_std(bootflux) / np.median(bootflux)

        # Measure flux in ancillary file (e.g. 8um or 13CO)
        # Assumed to have matched resolution for purposes of error estimation
        if ancfile is not None:
            if anchd['NAXIS'] == 2:
                collapsedmask = np.amax(asgn, axis = 0)
                collapsedmask[collapsedmask==0] = np.nan
                ngood = np.nansum(collapsedmask)
                #print('Pixels in region:',ngood)
                ancmean[j] = np.nanmean(ancdata*collapsedmask)
                ancrms[j]  = np.sqrt(np.nanmean((ancdata*collapsedmask)**2)-ancmean[j]**2)
                ancrms[j] /= np.sqrt(ngood/ppbeam)
                #print('Mean and uncertaintity:',ancmean[j],ancrms[j])
            else:
                # This probably still needs fixing to scale the rms properly
                ancmean[j] = np.nanmean(ancdata*asgn)
                ancrms[j]  = np.sqrt(np.nanmean((ancdata*asgn)**2)-ancmean[j]**2)

    # ---- report the median uncertainties
    print( "The median fractional error in rad_pc is {:2.4f}"
        .format(np.nanmedian(errms)) )
    print( "The median fractional error in rad_pc_dcon is {:2.4f}"
        .format(np.nanmedian(errms_d)) )
    print( "The median fractional error in vrms_k is {:2.4f}"
        .format(np.nanmedian(evrms)) )
    print( "The median fractional error in vrms_k_dcon is {:2.4f}"
        .format(np.nanmedian(evrms_d)) )
    print( "The median fractional error in mlumco is {:2.4f}"
        .format(np.nanmedian(eflux)) )
    
    # ---- apply a floor to the fractional uncertainty if requested
    if efloor > 0:
        print( "Applying a minimum fractional error of {:2.3f}".format(efloor) )
        errms[errms < efloor]     = efloor
        errms_d[errms_d < efloor] = efloor
        evrms[evrms < efloor]     = efloor
        evrms_d[evrms_d < efloor] = efloor
        eflux[eflux < efloor]     = efloor

    # ---- calculate the physical properties
    rms_pc  = (cat['radius'] * dist).to(u.pc, equivalencies=u.dimensionless_angles())
    rad_pc  = rmstorad * rms_pc
    deconmaj, deconmin, deconpa = deconvolve_gauss(
            cat['major_sigma'], cat['minor_sigma'], cat['position_angle'], 
            beam_maj, beam_min, bpa)
    rms_pc_dcon = (np.sqrt(deconmaj * deconmin) * dist).to(
                   u.pc, equivalencies=u.dimensionless_angles())
    rad_pc_dcon = rmstorad * rms_pc_dcon
    v_rms   = cat['v_rms'].to(u.km/u.s)
    # Deconvolve the linewidth (new 25-sep-2022)
    valid = (v_rms > dv_resp*sigvchan)
    v_rms_dcon = np.full_like(v_rms, np.nan)
    v_rms_dcon[valid] = np.sqrt(v_rms[valid]**2 - (dv_resp*sigvchan)**2)
    axrat   = cat['minor_sigma'] / cat['major_sigma']
    ellarea = (cat['area_ellipse']*dist**2).to(
               u.pc**2,equivalencies=u.dimensionless_angles())
    xctarea = (cat['area_exact']*dist**2).to(
               u.pc**2,equivalencies=u.dimensionless_angles())
    # Get CO luminosity from primary beam corrected data if available
    if copbcor is not None and conoise is not None:
        mlumco = alphaco * tb12*u.K * deltav * asarea * cdelt2.value**2
    else:
        # lumco = Luminosity in K km/s pc^2
        lumco  = deltav * asarea * (cat['flux']).to(
                 u.K,equivalencies=u.brightness_temperature(freq, beam_area=as2))
        mlumco = alphaco * lumco
    siglum = mlumco / xctarea
    # Virial mass using Rosolowsky+ 08
    mvir   = (5*rmstorad*v_rms**2*rms_pc/const.G).to(u.solMass)
    mvir_d = (5*rmstorad*v_rms_dcon**2*rms_pc_dcon/const.G).to(u.solMass)
    # Prefer deconvolved values for virial calculations
    sigvir = mvir_d / xctarea
    alpha  = mvir_d / mlumco
    # Match the error NaNs to the value NaNs
    errms[np.isnan(rms_pc)] = np.nan
    errms_d[np.isnan(rad_pc_dcon)] = np.nan
    evrms[np.isnan(v_rms)] = np.nan
    evrms_d[np.isnan(v_rms_dcon)] = np.nan
    # Standard error propagation for products and ratios
    emvir    = np.sqrt(2*evrms**2 + errms**2)
    emvir_d  = np.sqrt(2*evrms_d**2 + errms_d**2)
    ealpha   = np.sqrt(eflux**2 + emvir_d**2)

    # ---- calculate the projected distance from the reference position in arcsec
    if refpos is not None:
        refdist = abs(cdelt2) * np.sqrt(
                  (xref-cat['x_cen'].value)**2 + (yref-cat['y_cen'].value)**2 )
    
    # ---- make the physical properties table
    ptab = Table()
    ptab['_idx']      = Column(srclist)
    ptab['area_pc2']  = Column(xctarea, description='projected area of structure')
    ptab['rad_pc']    = Column(rad_pc, description='equivalent radius in pc')
    ptab['e_rad_pc']  = Column(errms, description='frac error in radius')
    ptab['rad_pc_dcon']    = Column(rad_pc_dcon, description='deconvolved radius in pc')
    ptab['e_rad_pc_dcon']  = Column(errms_d, description='frac error in deconvolved radius')
    ptab['vrms_k']    = Column(v_rms, description='rms linewidth in km/s')
    ptab['e_vrms_k']  = Column(evrms, description='frac error in linewidth')
    ptab['vrms_k_dcon'] = Column(v_rms_dcon, description='deconvolved rms linewidth in km/s')
    ptab['e_vrms_k_dcon']  = Column(evrms_d, description='frac error in deconvolved linewidth')
    ptab['axratio']   = Column(axrat, unit='', description='minor to major axis ratio')
    ptab['e_axratio'] = Column(eaxra, description='frac error in axis ratio')
    if copbcor is not None and conoise is not None:
        ptab['mlumco']= Column(mlumco, description='CO-based mass with alphascale='+str(alphascale)+' from '+os.path.basename(copbcor))
    else:
        ptab['mlumco']= Column(mlumco, description='Mass from scaling luminosity with alphascale='+str(alphascale))
    ptab['e_mlumco']  = Column(eflux, description='frac error in luminous mass')
    ptab['siglum']    = Column(siglum, description='average surface density from mlumco')
    ptab['e_siglum']  = Column(eflux, description='same as e_mlumco')
    ptab['mvir_raw']  = Column(mvir, description='virial mass using measured radius and lwidth')
    ptab['e_mvir_raw']= Column(emvir, description='frac error in mvir_raw')
    ptab['mvir']      = Column(mvir_d, description='virial mass using deconvolved radius and lwidth')
    ptab['e_mvir']    = Column(emvir_d, description='frac error in virial mass')
    ptab['sigvir']    = Column(sigvir, description='virial surface density from mvir')
    ptab['e_sigvir']  = Column(emvir_d, description='same as e_mvir')
    ptab['alpha']     = Column(alpha, unit='', description='virial parameter from mvir')
    ptab['e_alpha']   = Column(ealpha, description='frac error in virial parameter')
    if ancfile is not None:
        if anclabel is None:
            anclabel = ancimg.replace('.','_').split('_')[1]
        ptab[anclabel] = Column(ancmean, unit=anchd['BUNIT'])
        ancferr = ancrms / ancmean
        ptab['e_'+anclabel] = Column(ancferr)
    if refpos is not None:
        ptab['refdist'] = Column(refdist, description='angular offset from '+str(refpos))
    ptab.write(label+'_physprop.txt', format='ascii.ecsv', overwrite=True)
    return

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
        refdist = abs(cdelt2) * np.sqrt(
                  (xref-cat['x_cen'])**2 + (yref-cat['y_cen'])**2 )
        ptab['refdist'] = Column(refdist, description='angular offset from '+str(refpos))
        ptab.write(label+outfile, format='ascii.ecsv', overwrite=True)
    return
