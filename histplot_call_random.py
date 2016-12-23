from histplot import histplot
import numpy as np
np.set_printoptions(threshold=np.inf)

x1 = '30Dor_noise_25_itrs_rms_n13col.fits.gz'
x2 = '30Dor_noise_25_itrs_rms_n13cube.fits.gz'
x3 = '30Dor_noise_25_itrs_rms_tau13.fits.gz'
y1 = '30Dor_noise_25_itrs_oldrms_n13col.fits.gz'
y2 = '30Dor_noise_25_itrs_oldrms_n13cube.fits.gz'
y3 = '30Dor_noise_25_itrs_oldrms_tau13.fits.gz'

out1 = 'old_vs_new_rms_n13col_nolog'
out2 = 'old_vs_new_rms_n13cube_nolog'
out3 = 'old_vs_new_rms_tau13_nolog'

histplot(xname = x1, yname = y1, snrcut = 0, dolog2d = False, dolog1d = False, nbins = 100, outname = out1)
histplot(xname = x2, yname = y2, snrcut = 0, dolog2d = False, dolog1d = False, nbins = 100, outname = out2)
histplot(xname = x3, yname = y3, snrcut = 0, dolog2d = False, dolog1d = False, nbins = 100, outname = out3)
