#split(vis="../12m/calibrated_final.ms",
#        outputvis="../12m/calibrated_src.ms",
#        datacolumn="data",field="GMC1",spw="1,5,9")

### Define clean parameters
vis='../12m/calibrated_src.ms'
prename='GMC1_12CO_12m'
imsize=[800, 800]
cell='0.5arcsec'
minpb=0.2
restfreq='115.2712GHz'
outframe='LSRK'
spw=''
width='0.2km/s'
start='220km/s'
nchan=200
phasecenter='J2000 04h47m30.8s -69d10m32s'
scales=[0,4,12]

### Make initial dirty image
os.system('rm -rf '+prename+'.* ' +prename+'_*')
tclean(vis=vis,
       imagename=prename,
       gridder='mosaic',
       pblimit=minpb,
       imsize=imsize,
       cell=cell,spw=spw,
       weighting='natural',
       phasecenter=phasecenter,
       specmode='cube',
       width=width,
       start=start,
       nchan=nchan,
       restoringbeam='common',
       restfreq=restfreq,
       outframe=outframe,
       veltype='radio',
       niter=5000,
       scales = scales,
       smallscalebias=0.6,
       deconvolver = 'multiscale',
       maskresolution = 2.0,
       usemask = 'auto-multithresh',
       sidelobethreshold=2.0,
       noisethreshold=2.0,
       lownoisethreshold=1.5,
       smoothfactor = 0.7,
       minbeamfrac=0.0,
       interactive=False)
