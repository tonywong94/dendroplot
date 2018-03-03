import os
import casac
from tasks import *

'''
Required parameters and possible options:
name = ['GMC1' | 'GMC104' | 'A439' | 'N59C' | 'N113']
line = ['12CO' | '13CO' | 'CS' | 'C18O']
vis12m = '../12m/calibrated_final.ms'
vis7m  = '../7m/calibrated_final.ms'
weighting = ['briggs' | 'natural']
deconvolver = ['clark' | 'multiscale']
usemask = ['pb' | 'auto-thresh' | 'auto-thresh2' | 'auto-multithresh']
'''

def run_tclean(name=None, line=None, level=None, vis12m=None, vis7m=None, 
        startmodel=None, weighting=None, deconvolver=None, usemask=None, mask=None,
        pblimit=0.2, outframe='LSRK', spw='', width='0.2km/s', niter=10000, 
        threshold=0.05, pbmask=0.5, scales=[0,4,12], smallscalebias=0.6,
        maskresolution=2., maskthreshold=4., sidelobethreshold=3.,
        noisethreshold=4., lownoisethreshold=2., smoothfactor=2.,
        cutthreshold=0.1, minbeamfrac=0.5):

    # Determine output prefix
    if name is None:
        print("Assuming region name is GMC1")
        name = 'GMC1'
    if line is None:
        print("Assuming line is 12CO")
        line = '12CO'
    if level is None:
        print("Assuming level is 10")
        line = '10'	
    restfreq = {'12CO(10)': '115.2712GHz',
				'12CO(21)': '230.538GHz',
                '13CO(10)': '110.201354GHz',
				'13CO(21)': '220.398684GHz',
				'CS(21)': '97.980953GHz',
				'C18O(10)': '109.782176GHz',
				'C18O(21)': '219.560358GHz'}
    if vis12m is None:
        if vis7m is None:
            raise ValueError("Must specify vis12m or vis7m")
        else:
            print ("Imaging 7m data")
            thisvis = vis7m
            arrcode = '7m'
    else:
        if vis7m is None:
            print ("Imaging 12m data")
            thisvis = vis12m
            arrcode = '12m'
        else:
            print ("Imaging both 12m and 7m data together")
            thisvis = vis12m + vis7m
            arrcode = '12m7m'
    if startmodel is not None:
        print ("Imaging with an initial TP model %s" % startmodel)
        arrcode = '12m7mTPM'
    else:
        startmodel = ''
	
	# print (spw)
	# print (thisvis)
	thisname = name+'_'+line+'_'+arrcode

    # Check parameter choices
    weight_types=['briggs', 'natural']
    if weighting not in weight_types:
        raise ValueError("Invalid weighting. Expected one of: %s" % weight_types)
    decon_types=['clark', 'multiscale']
    if deconvolver not in decon_types:
        raise ValueError("Invalid deconvolver. Expected one of: %s" % decon_types)
    mask_types=['pb','auto-thresh','auto-thresh2','auto-multithresh','user']
    if usemask not in mask_types:
        raise ValueError("Invalid usemask. Expected one of: %s" % mask_types)

    # Determine cloud-specific imaging parameters
    imsize1 = { 'GMC1': [1000,800], 'GMC104': [800, 800], 
                'A439': [800, 800], 'N59C': [800, 800],
				'N113': [500, 500]  }		#was [500, 500]
    imsize2 = { 'GMC1': [250, 250], 'GMC104': [250, 250], 
                'A439': [250, 250], 'N59C': [250, 250], 
				'N113': [128, 128] }

    if arrcode == '7m':
        thissize = imsize2[name]
        thiscell = '2arcsec'
    else:
        thissize = imsize1[name]
        thiscell = '0.5arcsec'
	if name == 'N113' and level == '21':
		if arrcode == '7m':
			thissize = [128, 128]
			thiscell = '2arcsec'
		else:
			thissize = [1000, 1000]
			thiscell = '0.2arcsec'

    nchan = { 'GMC1': 100, 'GMC104': 100, 
              'A439': 150, 'N59C': 250, 
			  'N113': 150 }                #was 'N59C': 200
    vstart = { 'GMC1': '230km/s', 'GMC104': '216km/s', 
               'A439': '210km/s', 'N59C': '263km/s',
			   'N113': '220km/s' }   #was 'N59C': '266km/s'
    phasecenter = { 'GMC1': 'J2000 04h47m30.8s -69d10m32s',
                    'GMC104': 'J2000 05h21m05.5s -70d13m36s',
                    'A439': 'J2000 05h47m26.1s -69d52m46s',
                    'N59C': 'J2000 05h35m18.8s -67d36m12s',
					'N113': 'J2000 05h13m21.0s -69d22m21s' }
    ### Make the image
    #os.system('rm -rf '+thisname+'.* ' +thisname+'_*')
    tclean(vis=thisvis,
           imagename=thisname,
           startmodel=startmodel,
           gridder='mosaic',
           pblimit=pblimit,
           imsize=thissize,
           cell=thiscell,
           spw=spw,
           weighting=weighting,
           phasecenter=phasecenter[name],
           specmode='cube',
           width=width,
           start=vstart[name],
           nchan=nchan[name],
           restoringbeam='common',
           restfreq=restfreq[line+'('+level+')'],
           outframe=outframe,
           veltype='radio',
           niter=niter,
           scales=scales,
           smallscalebias=smallscalebias,
           deconvolver=deconvolver,
           usemask=usemask,
		   mask=mask,
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
    imhdr = imhead(thisname+'.image',mode='list')
    imsmooth(imagename=thisname+'.model',targetres=True,
        major=str(imhdr['beammajor']['value'])+imhdr['beammajor']['unit'],
        minor=str(imhdr['beamminor']['value'])+imhdr['beamminor']['unit'],
        pa=str(imhdr['beampa']['value'])+imhdr['beampa']['unit'],
        outfile=thisname+'.convmodel', overwrite=True)
    # Export to FITS
    exportfits(imagename=thisname+'.image',fitsimage=thisname+'.image.fits',
        dropdeg=True,velocity=True,overwrite=True)
    exportfits(imagename=thisname+'.pb',fitsimage=thisname+'.pb.fits',
        dropdeg=True,velocity=True,overwrite=True)
    exportfits(imagename=thisname+'.image.pbcor',fitsimage=thisname+'.pbcor.fits',
        dropdeg=True,velocity=True,overwrite=True)
    exportfits(imagename=thisname+'.residual',fitsimage=thisname+'.residual.fits',
        dropdeg=True,velocity=True,overwrite=True)
    exportfits(imagename=thisname+'.convmodel',fitsimage=thisname+'.convmodel.fits',
        dropdeg=True,velocity=True,overwrite=True)
    exportfits(imagename=thisname+'.mask',fitsimage=thisname+'.mask.fits',
        dropdeg=True,velocity=True,overwrite=True)

    return
