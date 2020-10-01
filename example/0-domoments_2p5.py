#!/usr/bin/env python

import os
import sys
sys.path.append(os.path.expanduser('~/Work/bin/py-package/maskmoment/'))
from maskmoment import maskmoment

# Moment maps by dilated masking

# Directory to write output files.
outdir = os.getcwd()+'/mom'
if not os.path.isdir(outdir):
    os.makedirs(outdir)

lines = ['12CO', '13CO']

for line in lines:
    maskmoment(img_fits='30dor_'+line+'21_2p5as.pbcor.fits.gz', 
        rms_fits='30dor_'+line+'21_2p5as.rms3d.fits.gz', outdir=outdir,
        snr_hi=4, snr_lo=2, minbeam=2, snr_hi_minch=2, nguard=[1,1], 
        outname='30dor_'+line+'21_2p5as_dil', output_snr_peak=True)

