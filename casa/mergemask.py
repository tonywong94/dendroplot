import os
import casac
from tasks import *

def mergemask(mask1=None, mask2=None,mask3=None,name='GMC104'):
	regrid=name+'_12CO_12m.maskrgd'
	merge=name+'_12CO_12m7m.mergemask'
	imregrid(imagename=mask1,template=mask2,output=regrid,interpolation='nearest',overwrite=True)

	os.system('rm -rf '+merge)
	immath(imagename=[mask2,regrid],expr='max(IM0,IM1)',outfile=merge)
	if mask3 is not None:
		imregrid(imagename=mask3,template=mask2,output=regrid+'2',interpolation='nearest',overwrite=True)

		os.system('rm -rf '+merge+'2')
		immath(imagename=[merge,regrid+'2'],expr='IM0*IM1',outfile=merge+'2')

		os.system('rm -rf '+regrid+'2')
	#exportfits(imagename=merge,fitsimage=merge+'.fits',dropdeg=True,velocity=True,overwrite=True)
	os.system('rm -rf '+regrid)
