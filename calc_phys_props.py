#!/usr/bin/env python

import numpy as np
from astropy import units as u
from astropy import constants as const
from astrodendro import Dendrogram, ppv_catalog
from astrodendro.analysis import ScalarStatistic, PPVStatistic
from astropy.table import Table, Column
from astropy.io.fits import getdata
from astropy.stats import mad_std

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
        #print bstats.major_sigma
        emmajs.append(bstats.major_sigma.value)
        emmins.append(bstats.minor_sigma.value)
        emomvs.append(bstats.v_rms.value)
        eflux.append(bstats.flux.value)
        pa.append(bstats.position_angle.value)
    return emmajs, emmins, emomvs, eflux, pa


def calc_phys_props(label='pcc_12', cubefile=None, boot_iter=100):

    rmstorad= 1.91
    alphaco = 4.3 * u.solMass * u.s / (u.K * u.km * u.pc**2) # Bolatto+ 13
    dist    = 5.e4 * u.pc
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
        print "\nWarning: Unrecognized brightness unit"
    metadata['vaxis'] = 0
    freq = hd3['restfrq'] * u.Hz
    metadata['wavelength'] = freq.to(u.m,equivalencies=u.spectral())
    metadata['spatial_scale']  =  hd3['cdelt2'] * 3600. * u.arcsec
    metadata['velocity_scale'] =  hd3['cdelt3'] * u.meter / u.second
    cdelt1 = hd3['cdelt1']*3600. * u.arcsec
    cdelt2 = hd3['cdelt1']*3600. * u.arcsec
    deltav = hd3['cdelt3']/1000. * u.km / u.s
    bmaj = hd3['bmaj']*3600. * u.arcsec # FWHM
    bmin = hd3['bmin']*3600. * u.arcsec # FWHM
    metadata['beam_major'] = bmaj
    metadata['beam_minor'] = bmin
    ppbeam = np.abs((bmaj*bmin)/(cdelt1*cdelt2)*2*np.pi/(8*np.log(2)))
    print "\nPixels per beam:", ppbeam
    indfac = np.sqrt(ppbeam)

    # ---- call the bootstrapping routine
    emaj, emin, epa, evrms, errms, eaxra, eflux, emlum, emvir, ealpha = [
        np.zeros(len(srclist)) for _ in xrange(10)]
    print "Calculating property errors..."
    for j, clust in enumerate(srclist):
        asgn = np.zeros(cube.shape)
        asgn[d[clust].get_mask(shape = asgn.shape)] = 1
        sindices = np.where(asgn == 1)
        svalues = cube[sindices]

        emmajs, emmins, emomvs, emom0s, pa = clustbootstrap(
            sindices, svalues, metadata, boot_iter)

        bootmaj   = np.asarray(emmajs) * u.arcsec
        bootmin   = np.asarray(emmins) * u.arcsec
        bootpa    = np.asarray(pa) * u.deg
        bootvrms  = (np.asarray(emomvs)*u.m/u.s).to(u.km/u.s)
        bootrrms  = (np.sqrt(np.asarray(emmajs)*np.asarray(emmins))*u.arcsec*dist).to(
            u.pc, equivalencies=u.dimensionless_angles())
        bootaxrat = bootmin/bootmaj
        bootflux  = np.asarray(emom0s) * u.Jy
        bootmvir  = (5*rmstorad*bootvrms**2*bootrrms/const.G).to(u.solMass)
        bootmlum  = alphaco*deltav*asarea*(bootflux).to(
            u.K, equivalencies=u.brightness_temperature(as2,freq))
        bootalpha = bootmvir/bootmlum

        emaj[j]   = indfac * mad_std(bootmaj)  / np.median(bootmaj)
        emin[j]   = indfac * mad_std(bootmin)  / np.median(bootmin)
        epa[j]    = indfac * mad_std(bootpa)   / np.median(bootpa)
        evrms[j]  = indfac * mad_std(bootvrms) / np.median(bootvrms)
        errms[j]  = indfac * mad_std(bootrrms) / np.median(bootrrms)    
        eaxra[j]  = indfac * mad_std(bootaxrat)/ np.median(bootaxrat)
        eflux[j]  = indfac * mad_std(bootflux) / np.median(bootflux)
        emlum[j]  = eflux[j]
        emvir[j]  = indfac * mad_std(bootmvir) / np.median(bootmvir)
        ealpha[j] = indfac * mad_std(bootalpha)/ np.median(bootalpha)

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
    # lumco = Luminosity in K km/s pc^2
    lumco  = deltav * asarea * (cat['flux']).to(
        u.K,equivalencies=u.brightness_temperature(as2,freq))
    mlumco = alphaco * lumco
    siglum = mlumco/xctarea
    mvir   = (5*rmstorad*v_rms**2*rms_pc/const.G).to(u.solMass)  # Rosolowsky+ 08
    sigvir = mvir / xctarea
    alpha  = mvir / mlumco

    # ---- make the physical properties table
    ptab = Table()
    ptab['_idx']      = Column(srclist)
    ptab['area_pc2']  = Column(xctarea)
    ptab['rad_pc']    = Column(rad_pc)
    ptab['e_rad_pc']  = Column(errms)
    ptab['vrms_k']    = Column(v_rms)
    ptab['e_vrms_k']  = Column(evrms)
    ptab['axratio']   = Column(axrat, unit='')
    ptab['e_axratio'] = Column(eaxra)
    ptab['mlumco']    = Column(mlumco)
    ptab['e_mlumco']  = Column(eflux)
    ptab['siglum']    = Column(siglum, description='average surface density')
    ptab['e_siglum']  = Column(eflux)
    ptab['mvir']      = Column(mvir)
    ptab['e_mvir']    = Column(emvir)
    ptab['sigvir']    = Column(sigvir, description='virial surface density')
    ptab['e_sigvir']  = Column(emvir)
    ptab['alpha']     = Column(alpha, unit='', description='virial parameter')
    ptab['e_alpha']   = Column(ealpha)
    ptab.write(label+'_physprop.txt', format='ascii.ecsv')

