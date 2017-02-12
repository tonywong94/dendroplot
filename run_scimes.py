#!/usr/bin/env python2

import os
import csv
import sys
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Ellipse
from astrodendro.scatter import Scatter
from astrodendro import Dendrogram, ppv_catalog, analysis
from astropy import stats
from astropy import units as u
from astropy.io import fits
from astropy.table import Table, Column
from scimes import SpectralCloudstering

def run_scimes(criteria=['volume'], label='scimes', cubefile=None, mom0file=None):

    #%&%&%&%&%&%&%&%&%&%&%&%
    #    Make dendrogram
    #%&%&%&%&%&%&%&%&%&%&%&%
    print 'Make dendrogram from the full cube'
    hdu3 = fits.open(cubefile)[0]
    hd3 = hdu3.header

    # Deal with oddities in 30 Dor cube
    if hd3['NAXIS'] == 3:
        for key in ['CTYPE4', 'CRVAL4', 'CDELT4', 'CRPIX4', 'CUNIT4', 'NAXIS4']:
            if key in hd3.keys():
                hd3.remove(key)

    # Get cube parameters
    sigma = stats.mad_std(hdu3.data[~np.isnan(hdu3.data)])
    print 'Robustly estimated RMS: ',sigma
    ppb = 1.133*hd3['bmaj']*hd3['bmin']/(abs(hd3['cdelt1']*hd3['cdelt2']))
    print 'Pixels per beam: ',ppb

    # Make the dendrogram
    d = Dendrogram.compute(hdu3.data, min_value=3*sigma, \
                    min_delta=2.5*sigma, min_npix=2*ppb, verbose = 1)
    d.save_to(label+'_dendrogram.hdf5')

    # Plot the tree
    fig = plt.figure(figsize=(14, 8))
    ax = fig.add_subplot(111)            
    ax.set_yscale('log')
    ax.set_xlabel('Structure')
    ax.set_ylabel('Intensity ['+hd3['BUNIT']+']')
    p = d.plotter()
    p.plot_tree(ax, color='black')
    plt.savefig(label+'_dendrogram.pdf')

    #%&%&%&%&%&%&%&%&%&%&%&%&%&%
    #   Generate the catalog
    #%&%&%&%&%&%&%&%&%&%&%&%&%&%
    print "Generate a catalog of dendrogram structures"
    metadata = {}
    if hd3['BUNIT'].upper()=='JY/BEAM':
        metadata['data_unit'] = u.Jy / u.beam
    elif hd3['BUNIT'].upper()=='K':
        metadata['data_unit'] = u.K
    else:
        print "Warning: Unrecognized brightness unit"
    metadata['vaxis'] = 0
    freq = hd3['restfrq'] * u.Hz
    metadata['wavelength'] = freq.to(u.m,equivalencies=u.spectral())
    metadata['spatial_scale']  =  hd3['cdelt2'] * 3600. * u.arcsec
    metadata['velocity_scale'] =  hd3['cdelt3'] * u.meter / u.second
    bmaj = hd3['bmaj']*3600. * u.arcsec # FWHM
    bmin = hd3['bmin']*3600. * u.arcsec # FWHM
    metadata['beam_major'] = bmaj
    metadata['beam_minor'] = bmin

    cat = ppv_catalog(d, metadata)
    print cat.info()

    # Add additional properties: Average Peak Tb and Maximum Tb
    srclist = cat['_idx'].tolist()
    tmax  = np.zeros(len(srclist), dtype=np.float64)
    tpkav = np.zeros(len(srclist), dtype=np.float64)
    for i, c in enumerate(srclist):
        peakim = np.nanmax(hdu3.data*d[c].get_mask(), axis=0)
        peakim[peakim==0] = np.nan
        tmax[i]  = np.nanmax(peakim)
        tpkav[i] = np.nanmean(peakim)
    if hd3['BUNIT'].upper()=='JY/BEAM':
        omega_B = np.pi/(4*np.log(2)) * bmaj * bmin
        convfac = (u.Jy).to(u.K, equivalencies=u.brightness_temperature(omega_B,freq))
        tmax *= convfac
        tpkav *= convfac
    newcol = Column(tmax, name='tmax')
    newcol.unit = 'K'
    cat.add_column(newcol)
    newcol = Column(tpkav, name='tpkav')
    newcol.unit = 'K'
    cat.add_column(newcol)

    cat.write(label+'_full_catalog.txt', format='ascii.ecsv')

    #%&%&%&%&%&%&%&%&%&%&%&%&%&%
    #     Running SCIMES
    #%&%&%&%&%&%&%&%&%&%&%&%&%&%
    print "Running SCIMES"
    dclust = SpectralCloudstering(d, cat, criteria = criteria)
    print dclust.clusters

    print "Visualize the clustered dendrogram"
    dclust.showdendro()
    plt.savefig(label+'_clusters_tree.pdf')

    print "Produce the assignment cube"
    dclust.asgncube(hd3)
    try:
        os.remove(label+'_asgncube.fits.gz')
    except OSError:
        pass
    dclust.asgn.writeto(label+'_asgncube.fits.gz')

    #sys.exit("Stopping here")

    #%&%&%&%&%&%&%&%&%&%&%&%&%&%
    #     Image the trunks
    #%&%&%&%&%&%&%&%&%&%&%&%&%&%
    print "Image the trunks"

    hdu2 = fits.open(mom0file)[0]

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    vmax = np.nanmax(hdu2.data)/2.
    im = ax.matshow(hdu2.data, origin='lower', cmap=plt.cm.Blues, vmax=vmax)

    # Make a trunk list
    tronly = [s for s in d.trunk if s not in d.leaves]
    f = open(label+'_trunks.txt', 'w')

    for c in tronly:
        f.write('{:<4d} | '.format(c.idx))
        # Plot the actual structure boundaries
        mask = d[c.idx].get_mask()
        mask_coll = np.amax(mask, axis = 0)
        plt.contour(mask_coll, colors='red', linewidths=1, levels = [0])
        # Plot the ellipse fits
        s = analysis.PPVStatistic(d[c.idx])
        ellipse = s.to_mpl_ellipse(edgecolor='black', facecolor='none')
        ax.add_patch(ellipse)
        # Make sub-lists of descendants
        print 'Finding descendants of trunk ',c.idx
        desclist = []
        if len(d[c.idx].descendants) > 0:
            for s in d[c.idx].descendants:
                desclist.append(s.idx)
            desclist.sort()
            liststr=','.join(map(str, desclist))
            f.write(liststr)
        f.write("\n")
    f.close()

    fig.colorbar(im, ax=ax)
    plt.savefig(label+'_trunks.pdf')
    plt.close()

    # Make a branch list
    branch = [s for s in d.all_structures if s not in d.leaves and s not in d.trunk]
    slist = []
    for c in branch:
        slist.append(c.idx)
    slist.sort()
    with open(label+'_branches.txt', 'w') as output:
        writer = csv.writer(output)
        for val in slist:
            writer.writerow([val])    

    #%&%&%&%&%&%&%&%&%&%&%&%&%&%
    #     Image the leaves
    #%&%&%&%&%&%&%&%&%&%&%&%&%&%
    print "Image the leaves"

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    vmax = np.nanmax(hdu2.data)/2.
    im = ax.matshow(hdu2.data, origin='lower', cmap=plt.cm.Blues, vmax=vmax)

    # Make a leaf list
    slist = []
    for c in d.leaves:
        slist.append(c.idx)
        # Plot the actual structure boundaries
        mask = d[c.idx].get_mask()
        mask_coll = np.amax(mask, axis = 0)
        plt.contour(mask_coll, colors='green', linewidths=1, levels = [0])
        # Plot the ellipse fits
        s = analysis.PPVStatistic(d[c.idx])
        ellipse = s.to_mpl_ellipse(edgecolor='black', facecolor='none')
        ax.add_patch(ellipse)
    slist.sort()
    with open(label+'_leaves.txt', "w") as output:
        writer = csv.writer(output)
        for val in slist:
            writer.writerow([val])    

    fig.colorbar(im, ax=ax)
    plt.savefig(label+'_leaves.pdf')
    plt.close()

    #%&%&%&%&%&%&%&%&%&%&%&%&%&%
    #     Image the clusters
    #%&%&%&%&%&%&%&%&%&%&%&%&%&%
    print "Image the resulting clusters"

    clusts = np.array(dclust.clusters)
    colors = np.array(dclust.colors)
    inds = np.argsort(clusts)
    clusts = clusts[inds]
    colors = colors[inds]

    print clusts
    print colors

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    vmax = np.nanmax(hdu2.data)/2.
    im = ax.matshow(hdu2.data, origin='lower', cmap=plt.cm.Blues, vmax=vmax)

    # Make a cluster list
    f = open(label+'_clusters.txt', 'w')
    for i, c in enumerate(clusts):
        f.write('{:<4d} {:7} | '.format(c,colors[i]))
        # Plot the actual structure boundaries
        mask = d[c].get_mask()
        mask_coll = np.amax(mask, axis = 0)
        plt.contour(mask_coll, colors=colors[i], linewidths=1, levels = [0])
        # Plot the ellipse fits
        s = analysis.PPVStatistic(d[c])
        ellipse = s.to_mpl_ellipse(edgecolor='black', facecolor='none')
        ax.add_patch(ellipse)
        # Make sub-lists of descendants
        print 'Finding descendants of cluster ',c
        desclist = []
        if len(d[c].descendants) > 0:
            for s in d[c].descendants:
                desclist.append(s.idx)
            desclist.sort()
            liststr=','.join(map(str, desclist))
            f.write(liststr)
        f.write("\n")
    f.close()

    fig.colorbar(im, ax=ax)
    plt.savefig(label+'_clusters_map.pdf')
    plt.close()

    #%&%&%&%&%&%&%&%&%&%&%&%&%&%
    #     Image structures color-coded by properties
    #%&%&%&%&%&%&%&%&%&%&%&%&%&%
    print "Image leaves and clusters colored by properties"

    for set in ['leaves', 'clusters']:
        with open(label+'_'+set+'.txt', 'r') as f:
            reader=csv.reader(f, delimiter=' ')
            a = zip(*reader)
        idc = map(int,a[0])
        subcat = cat[idc]
        #subcat = Table(rows=cat[idc])
        #cat = Table.read(label+'_'+set+'.txt', format='ascii.ecsv')
        srclist = subcat['_idx'].tolist()
        for type in ['v_cen', 'v_rms', 'flux']:
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            imax = np.nanmax(hdu2.data)
            im = ax.matshow(hdu2.data, origin='lower', cmap=plt.cm.Blues, vmax=imax)
            plt.setp(ax.get_xticklabels(), visible=False)
            plt.setp(ax.get_yticklabels(), visible=False)
            if type == 'flux':  # Convert Jy*ch to K km/s
                deltav = hd3['cdelt3']/1000. * u.km / u.s
                as2 = 1 * u.arcsec**2
                kflux = deltav * as2 * subcat[type].to(u.K,
                    equivalencies=u.brightness_temperature(as2,freq))
                subcat[type] = kflux/subcat['area_exact']
                subcat[type].unit = 'K km / s'
            vmin=np.floor(subcat[type].min())
            vmax=np.ceil(subcat[type].max())
            print label,set,type,' vmin and vmax:',vmin,vmax
            for i, c in enumerate(srclist):
                s = analysis.PPVStatistic(d[c])
                scaled_v = (subcat[type][i]-vmin)/(vmax-vmin)
                col=plt.cm.nipy_spectral_r(scaled_v)
                ellipse = s.to_mpl_ellipse(color=col)
                ax.add_patch(ellipse)
            ax1 = fig.add_axes([0.3, 0.95, 0.45, 0.02])
            cmap = plt.cm.get_cmap('nipy_spectral_r')
            if type == 'v_cen':  # convert ch number to km/s
                v0 = 1.e-3*(hd3['crval3'] + hd3['cdelt3']*(vmin-hd3['crpix3']))
                v1 = 1.e-3*(hd3['crval3'] + hd3['cdelt3']*(vmax-hd3['crpix3']))
                subcat[type].unit = 'km / s'
            else:
                v0 = vmin
                v1 = vmax
            cbar = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                 orientation='horizontal',norm=mpl.colors.Normalize(vmin=v0, vmax=v1))
            cbar.ax.tick_params(labelsize=9)
            if type == 'flux':
                name = 'mean brightness'
            elif type == 'v_cen':
                name = 'mean velocity'
            elif type == 'v_rms':
                name = 'velocity dispersion'
            cbar.set_label(name+' ['+str(subcat[type].unit)+']',size=9,labelpad=-35)
            plt.savefig(label+'_'+set+'_'+type.replace("_","")+'.pdf', 
                bbox_inches='tight')
            plt.close()

    return

# -------------------------------------------------------------------------------

def explore_dendro(label='scimes', xaxis='radius', yaxis='v_rms'):
    d = Dendrogram.load_from(label+'_dendrogram.hdf5')
    cat = Table.read(label+'_full_catalog.txt', format='ascii.ecsv')
    dv = d.viewer()
    ds = Scatter(d, dv.hub, cat, xaxis, yaxis)
    ds.set_loglog()
    dv.show()
    return

# -------------------------------------------------------------------------------

if __name__ == "__main__":
    run_scimes(label=sys.argv[1], cubefile=sys.argv[2], mom0file=sys.argv[3])
