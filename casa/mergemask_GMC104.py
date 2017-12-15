from mergemasks import mergemask

mergemask(mask1='../7m_5/GMC104_12CO_7m.mask',mask2='../12m_4/GMC104_12CO_12m.mask',name='GMC104')

#gzip GMC104_12CO_12m7m.mergemask.fits

"""
names=['GMC1','A439','GMC104','N59C']
name=names[2]
mask1='../7m_5/'+name+'_12CO_7m.mask'
mask2='../12m_4/'+name+'_12CO_12m.mask'
regrid=name+'_12CO_12m.maskrgd'
merge=name+'_12CO_12m7m.mergemask'
imregrid(imagename=mask1,template=mask2,output=regrid,interpolation='nearest',overwrite=True)

os.system('rm -rf '+merge)
immath(imagename=[mask2,regrid],expr='max(IM0,IM1)',outfile=merge)
exportfits(imagename=merge,fitsimage=merge+'.fits',dropdeg=True,velocity=True,overwrite=True)
"""
