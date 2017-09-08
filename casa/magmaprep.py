
from astropy import units as u

# Get region name
for file in os.listdir('.'):
    if file.endswith('_12CO_12m7m.image'):
        regname=file.split('_')[0]
        print("The region name is %s" % regname)

# Get scaling from K to Jy/pix
cdelt2=imhead(imagename=regname+'_12CO_12m7m.image',mode='get',hdkey='cdelt2')
aspp = (cdelt2['value'] * u.Unit(cdelt2['unit'])).to(u.arcsec)
freq=115.2712*u.GHz
convfac = (u.K).to(u.Jy, equivalencies=u.brightness_temperature(aspp**2,freq))
print("The conversion factor to Jy/pix is %s" % convfac)

# Prepare single dish image for use as startmodel in tclean
importfits(fitsimage='../magma/'+regname+'.cocub.fits', imagename=regname+'.magma.im',
    beam=['45arcsec','45arcsec','0deg'], defaultaxes=True,
    defaultaxesvalues=['','','','I',], overwrite=True)
os.system('rm -rf '+regname+'magma.lsrk')
imreframe(regname+'.magma.im',output=regname+'.magma.lsrk',outframe='lsrk') 
imregrid(regname+'.magma.lsrk',template=regname+'_12CO_12m7m.image',axes=[0,1,2,3], 
    output=regname+'.magma.rgd',interpolation='cubic',overwrite=True)
os.system('rm -rf '+regname+'magma.rgd.jypix')
immath(regname+'.magma.rgd',outfile=regname+'.magma.rgd.jypix',expr='IM0*'+str(convfac))
imhead(imagename=regname+'.magma.rgd.jypix',mode='put',hdkey='BUNIT',hdvalue='Jy/pixel')
imhead(imagename=regname+'.magma.rgd.jypix')
