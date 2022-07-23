#!/usr/bin/env python

from dendroplot.plotting.mom0plot import mom0plot

## Plot the integrated intensity maps.
## CO(2-1) listed first here

imdir = 'images/'

dolines  = ['12', '13']
doclouds = [ '30Dor',   'PCC',   'A439',   'GMC1','GMC104',   'N59C']
cmin     = [  [2, 0], [0, -2],  [1, -1],  [1, -1], [1, -1],  [1, -1]]
xrange   = [[70,400],[30,380],[150,650],[150,850],[50,530],[120,680]]
yrange   = [[50,380],[50,710],[150,650], [70,730],[60,540], [50,750]]

for i, cld in enumerate(doclouds):
    if cld == 'GMC1':
        fov = 0.6
    else:
        fov = 0.5
    for j, line in enumerate(dolines):
        if i < 2:
            linestr = line+'CO21'
            linelbl = '$^{'+line+'}$CO(2-1)'
        else:
            linestr = line+'CO'
            linelbl = '$^{'+line+'}$CO(1-0)'
        mom0plot(mom0file=imdir+cld+'_'+linestr+'_3p5as_dil.mom0.fits.gz',
            fluxfile=imdir+cld+'_'+linestr+'_3p5as.pb1.fits.gz',
            xrange=xrange[i], yrange=yrange[i], fov=fov,
            label=cld+' '+linelbl, labelax='tickonly', cmin=cmin[i][j],
            outfile=cld+'_'+linestr+'_3p5as_mom0.pdf')
