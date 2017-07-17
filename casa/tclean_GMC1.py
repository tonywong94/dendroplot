#split(vis="../12m/calibrated_final.ms",
#        outputvis="../12m/calibrated_src.ms",
#        datacolumn="data",field="GMC1",spw="1,5,9")

### Generic parameters for all runs
vis=['../12m/calibrated_src.ms','../7m/calibrated_src.ms']
prename=['GMC1_12CO_12m','GMC1_12CO_7m','GMC1_12CO_12m7m']
weighting=['briggs', 'natural']
deconvolver=['clark', 'multiscale']
usemask=['pb','auto-thresh','auto-thresh2','auto-multithresh']
imsize=[[1000, 800], [250, 200]]
cell=['0.5arcsec', '2arcsec']
minpb=0.2
restfreq='115.2712GHz'
outframe='LSRK'
spw=''
width='0.2km/s'
start='230km/s'
nchan=100
phasecenter='J2000 04h47m30.8s -69d10m32s'

### Specific parameters for this run
niter=10000
threshold=0.05
pbmask=0.5
thisweighting=weighting[0]
thisname=prename[2]
thisdecon=deconvolver[0]
thismask=usemask[0]

### Parameters for deconvolver='multiscale'
scales=[0,4,12]
smallscalebias=0.6

### Parameters for usemask='auto-thresh' and 'auto-thresh2'
maskresolution=2.
maskthreshold=4.
       
### Parameters for usemask='auto-multithresh'
sidelobethreshold=3.
noisethreshold=4.
lownoisethreshold=2.
smoothfactor=2.
cutthreshold=0.1
minbeamfrac=0.5
       
if thisname == prename[0]:
    thisvis=vis[0]
    thissize=imsize[0]
    thiscell=cell[0]
if thisname == prename[1]:
    thisvis=vis[1]
    thissize=imsize[1]
    thiscell=cell[1]
if thisname == prename[2]:
    thisvis=vis[0:2]
    thissize=imsize[0]
    thiscell=cell[0]

### Make the image
os.system('rm -rf '+prename+'.* ' +prename+'_*')
tclean(vis=thisvis,
       imagename=thisname,
       gridder='mosaic',
       pblimit=minpb,
       imsize=thissize,
       cell=thiscell,
       spw=spw,
       weighting=thisweighting,
       phasecenter=phasecenter,
       specmode='cube',
       width=width,
       start=start,
       nchan=nchan,
       restoringbeam='common',
       restfreq=restfreq,
       outframe=outframe,
       veltype='radio',
       niter=niter,
       scales=scales,
       smallscalebias=smallscalebias,
       deconvolver=thisdecon,
       usemask=thismask,
       maskresolution=maskresolution,
       maskthreshold=maskthreshold,
       sidelobethreshold=sidelobethreshold,
       noisethreshold=noisethreshold,
       lownoisethreshold=lownoisethreshold,
       smoothfactor=smoothfactor,
       cutthreshold=cutthreshold,
       threshold=threshold,
       minbeamfrac=minbeamfrac,
       pbmask=pbmask,
       pbcor=True,
       interactive=False)

# Convolve the model by the beam
imhdr = imhead(prename+'.image',mode='list')
imsmooth(imagename=prename+'.model',targetres=True,
    major=str(imhdr['beammajor']['value'])+imhdr['beammajor']['unit'],
    minor=str(imhdr['beamminor']['value'])+imhdr['beamminor']['unit'],
    pa=str(imhdr['beampa']['value'])+imhdr['beampa']['unit'],
    outfile=prename+'.convmodel', overwrite=True)

# Export to FITS
exportfits(imagename=prename+'.image',fitsimage=prename+'.image.fits',
    dropdeg=True,velocity=True,overwrite=True)
exportfits(imagename=prename+'.pb',fitsimage=prename+'.pb.fits',
    dropdeg=True,velocity=True,overwrite=True)
exportfits(imagename=prename+'.image.pbcor',fitsimage=prename+'.pbcor.fits',
    dropdeg=True,velocity=True,overwrite=True)
exportfits(imagename=prename+'.residual',fitsimage=prename+'.residual.fits',
    dropdeg=True,velocity=True,overwrite=True)
exportfits(imagename=prename+'.convmodel',fitsimage=prename+'.convmodel.fits',
    dropdeg=True,velocity=True)
exportfits(imagename=prename+'.mask',fitsimage=prename+'.mask.fits',
    dropdeg=True,velocity=True,overwrite=True)

