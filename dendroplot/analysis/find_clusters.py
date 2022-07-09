import os
import numpy as np
from astrodendro import Dendrogram, analysis
from astropy.io import fits
from astropy.table import Table
from scimes import SpectralCloudstering
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm, PowerNorm

## Take an existing dendrogram analysis and identify clusters in the data
## using the SCIMES package.
## The clusters are a subset of the dendrogram structures and are listed in
## an output text file.
## Also generates assignment cubes, a tree plot, and a cluster map.

def find_clusters(label='pcc_12', criteria=['volume'], doellipse=False, 
                  plotdir='plots', cubefile=None, rms=np.nan, noshow=False):

    if not os.path.isdir(plotdir):
        try:
            os.makedirs(plotdir)
        except OSError:
            print('Unable to create output directory',plotdir)

    d = Dendrogram.load_from(label+'_dendrogram.hdf5')
    cat = Table.read(label+'_full_catalog.txt', format='ascii.ecsv')
    srclist = cat['_idx'].tolist()

    # ---- load the cube and extract the metadata
    cubedata, hd3 = fits.getdata(cubefile, header=True)

    dclust = SpectralCloudstering(d, cat, criteria=criteria, save_branches=True, 
                                  rms=rms, header=hd3)
    dclust.showdendro(savefile=plotdir+'/'+label+'_dendrogram_clust.pdf')

    clusts = np.array(dclust.clusters)
    colors = np.array(dclust.colors)
    inds = np.argsort(clusts)
    clusts = clusts[inds]
    colors = colors[inds]

    # Prepare the cluster map
    fig, ax = plt.subplots()
    tpeak = np.amax(cubedata, axis=0)
    vmin = np.nanmin(tpeak)*1.5
    vmax = np.nanmax(tpeak)/1.5
    im = ax.imshow(tpeak, origin='lower', cmap=plt.cm.Greys, 
                   norm=PowerNorm(gamma=0.5, vmin=vmin, vmax=vmax))
    ax.axes.get_xaxis().set_ticks([])
    ax.axes.get_yaxis().set_ticks([])

    # Make a cluster list
    f = open(label+'_clusters.txt', 'w')
    for i, c in enumerate(clusts):
        f.write('{:<4d} {:7} | '.format(c,colors[i]))
        # Plot the actual structure boundaries
        mask = d[c].get_mask()
        mask_coll = np.amax(mask, axis=0)
        ax.contour(mask_coll, colors=colors[i], linewidths=1, levels=[0])
        # Plot the ellipse fits
        if doellipse:
            stat = analysis.PPVStatistic(d[c])
            ellipse = stat.to_mpl_ellipse(edgecolor='black', facecolor='none')
            ax.add_patch(ellipse)
        # Make sub-lists of descendants
        print('Finding descendants of cluster {}'.format(c))
        desclist = []
        if len(d[c].descendants) > 0:
            for s in d[c].descendants:
                desclist.append(s.idx)
            desclist.sort()
            liststr=','.join(map(str, desclist))
            f.write(liststr)
        f.write("\n")
    f.close()

    # Finalize the cluster map
    fig.colorbar(im, ax=ax)
    plt.savefig(plotdir+'/'+label+'_clusters_map.pdf', bbox_inches='tight')
    if noshow:
        plt.close()

    for struct in ['clusters', 'trunks', 'leaves']:
        hdu = getattr(dclust, struct+'_asgn')
        hdu.header['datamin'] = hdu.data.min()
        hdu.header['datamax'] = hdu.data.max()
        hdu.writeto(label+'_'+struct+'_asgn.fits.gz', overwrite=True)

    return
