import os
import casac
from tasks import *

def mergemask(name='GMC104', line='12CO', file7m=None, file12m=None, cropping=0.5):
	#cropping goes out to two decimals
    mask7m = file7m+'/'+name+'_'+line+'_7m.mask'
    mask12m = file12m+'/'+name+'_'+line+'_12m.mask'
    regrid = name+'_'+line+'_12m.maskrgd'
    merge = name+'_'+line+'_12m7m.mergemask'

    imregrid(imagename=mask7m,template=mask12m,output=regrid,interpolation='nearest',overwrite=True)
    os.system('rm -rf '+merge)
    immath(imagename=[mask12m,regrid],expr='max(IM0,IM1)',outfile=merge)

    if cropping is not None:
        maskpb = file12m+'/'+name+'_'+line+'_12m.pb'
        crop = name+'_'+line+'_12m7m.croppedmask'
        os.system('rm -rf '+ crop)
        immath(imagename=[merge,maskpb],expr='iif(IM1>%4.2f,IM0,0.)'%cropping,outfile=crop)
        #exportfits(imagename=crop,fitsimage=crop+'.fits',dropdeg=True,velocity=True,overwrite=True)

    #exportfits(imagename=merge,fitsimage=merge+'.fits',dropdeg=True,velocity=True,overwrite=True)
    os.system('rm -rf '+ regrid)

