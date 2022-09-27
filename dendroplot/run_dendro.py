import os
from os.path import join
import csv
import sys
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Ellipse
from astrodendro.scatter import Scatter
from astrodendro import Dendrogram, ppv_catalog, analysis
from astropy import stats
from astropy import units as u
from astropy.io import fits
from astropy.table import Table, Column

def run_dendro(label='mycloud', cubefile=None, mom0file=None, 
               redo='n', rfreq=None, verbose=True, noshow=False, 
               plotdir='plots', **kwargs):

    #%&%&%&%&%&%&%&%&%&%&%&%
    #    Make dendrogram
    #%&%&%&%&%&%&%&%&%&%&%&%
    hdu3 = fits.open(cubefile)[0]
    hd3 = hdu3.header

    # Deal with oddities in 30 Dor cube
    if hd3['NAXIS'] == 3:
        for key in ['CTYPE4', 'CRVAL4', 'CDELT4', 'CRPIX4', 'CUNIT4', 'NAXIS4']:
            if key in hd3.keys():
                hd3.remove(key)

    # Get cube parameters
    sigma = stats.mad_std(hdu3.data[~np.isnan(hdu3.data)])
    print('\nRobustly estimated RMS: {:.3f}'.format(sigma))
    ppb = 1.133*hd3['bmaj']*hd3['bmin']/(abs(hd3['cdelt1']*hd3['cdelt2']))
    print('Pixels per beam: {:.2f}'.format(ppb))

    # Make the dendrogram if not present or redo=y
    if redo == 'n' and os.path.isfile(label+'_dendrogram.hdf5'):
        print('Loading pre-existing dendrogram')
        d = Dendrogram.load_from(label+'_dendrogram.hdf5')
    else:
        print('Make dendrogram from the full cube')
        d = Dendrogram.compute(hdu3.data, min_value=3*sigma,
            min_delta=2.5*sigma, min_npix=2*ppb, verbose=verbose)
        d.save_to(label+'_dendrogram.hdf5')

    # checks/creates directory to place plots
    if not os.path.isdir(plotdir):
        os.makedirs(plotdir)

    # Plot the tree
    fig = plt.figure(figsize=(14, 8))
    ax = fig.add_subplot(111)            
    ax.set_xlabel('Structure')
    ax.set_ylabel('Intensity ['+hd3['BUNIT']+']')
    p = d.plotter()
    branch = [s for s in d.all_structures if s not in d.leaves and s not in d.trunk]
    tronly = [s for s in d.trunk if s not in d.leaves]
    for st in tronly:
        p.plot_tree(ax, structure=[st], color='brown', subtree=False)
    for st in branch:
        p.plot_tree(ax, structure=[st], color='black', subtree=False)
    for st in d.leaves:
        p.plot_tree(ax, structure=[st], color='green')
    ax.set_yscale('log')
    plt.savefig(join(plotdir, label+'_dendrogram.pdf'), bbox_inches='tight')
    if noshow:
        plt.close()

    #%&%&%&%&%&%&%&%&%&%&%&%&%&%
    #   Generate the catalog
    #%&%&%&%&%&%&%&%&%&%&%&%&%&%
    print('Generate a catalog of dendrogram structures')
    metadata = {}
    if hd3['BUNIT'].upper()=='JY/BEAM':
        metadata['data_unit'] = u.Jy / u.beam
    elif hd3['BUNIT'].upper()=='K':
        metadata['data_unit'] = u.K
    else:
        print('Warning: Unrecognized brightness unit')
    metadata['vaxis'] = 0
    if rfreq is None:
        if 'RESTFREQ' in hd3.keys():
            freq = hd3['RESTFREQ'] * u.Hz
        elif 'RESTFRQ' in hd3.keys():
            freq = hd3['RESTFRQ'] * u.Hz
    else:
        freq = rfreq * u.GHz
    metadata['wavelength'] = freq.to(u.m,equivalencies=u.spectral())
    metadata['spatial_scale']  =  abs(hd3['cdelt2']) * 3600. * u.arcsec
    metadata['velocity_scale'] =  abs(hd3['cdelt3']) * u.meter / u.second
    bmaj = hd3['bmaj']*3600. * u.arcsec # FWHM
    bmin = hd3['bmin']*3600. * u.arcsec # FWHM
    metadata['beam_major'] = bmaj
    metadata['beam_minor'] = bmin

    cat = ppv_catalog(d, metadata, verbose=verbose)
    print(cat.info)

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
        convfac = (u.Jy).to(u.K, equivalencies=u.brightness_temperature(freq,
                            beam_area=omega_B))
        tmax *= convfac
        tpkav *= convfac
    newcol = Column(tmax, name='tmax')
    newcol.unit = 'K'
    cat.add_column(newcol)
    newcol = Column(tpkav, name='tpkav')
    newcol.unit = 'K'
    cat.add_column(newcol)

    cat.write(label+'_full_catalog.txt', format='ascii.ecsv', overwrite=True)

    #%&%&%&%&%&%&%&%&%&%&%&%&%&%
    #     Image the trunks
    #%&%&%&%&%&%&%&%&%&%&%&%&%&%
    print("Image the trunks")

    hdu2 = fits.open(mom0file)[0]

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    vmax = np.nanmax(hdu2.data)/2.
    im = ax.matshow(hdu2.data, origin='lower', cmap=plt.cm.Blues, vmax=vmax)
    ax.axes.get_xaxis().set_ticks([])
    ax.axes.get_yaxis().set_ticks([])
    if 'xlims' in kwargs:
        ax.set_xlim(kwargs['xlims'])
    if 'ylims' in kwargs:
        ax.set_ylim(kwargs['ylims'])

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
        print('Finding descendants of trunk {}'.format(c.idx))
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
    plt.savefig(join(plotdir, label+'_trunks_map.pdf'), bbox_inches='tight')
    if noshow:
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
    print("Image the leaves")

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    vmax = np.nanmax(hdu2.data)/2.
    im = ax.matshow(hdu2.data, origin='lower', cmap=plt.cm.Blues, vmax=vmax)
    ax.axes.get_xaxis().set_ticks([])
    ax.axes.get_yaxis().set_ticks([])
    if 'xlims' in kwargs:
        ax.set_xlim(kwargs['xlims'])
    if 'ylims' in kwargs:
        ax.set_ylim(kwargs['ylims'])

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
    plt.savefig(join(plotdir, label+'_leaves_map.pdf'), bbox_inches='tight')
    if noshow:
        plt.close()


# -------------------------------------------------------------------------------

def explore_dendro(label='mycloud', xaxis='radius', yaxis='v_rms'):
    d = Dendrogram.load_from(label+'_dendrogram.hdf5')
    cat = Table.read(label+'_full_catalog.txt', format='ascii.ecsv')
    dv = d.viewer()
    ds = Scatter(d, dv.hub, cat, xaxis, yaxis)
    ds.set_loglog()
    dv.show()
    return

# -------------------------------------------------------------------------------

if __name__ == "__main__":
    run_dendro(label=sys.argv[1], cubefile=sys.argv[2], mom0file=sys.argv[3])
