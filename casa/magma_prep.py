import os
import math
import casac
from taskinit import *
from tasks import *

'''
    Takes an ALMA template image and a MAGMA cube in Jy/beam and applies the
    ALMA pb gain, converts to Jy/pixel for tclean, and optionally does a feather.

Required files:
    <prefix>.image: The ALMA cube to be used as the template (and in feathering)
        <prefix> should begin with <regname>_
    <prefix>.pb: The pb pattern (peaking at 1) corresponding to <prefix>.image
    ../magma/<regname>.cocub.jyb.fits: The MAGMA cube in Jy/beam units
'''

def magma_prep(prefix=None, dofeather=True):

    regname=prefix.split('_')[0]
    print("The region name is %s" % regname)

    # Regrid MAGMA image and apply ALMA pb gain
    importfits(fitsimage='../magma/'+regname+'.cocub.jyb.fits', 
        imagename=regname+'.magma.im',
        beam=['45arcsec','45arcsec','0deg'], defaultaxes=True,
        defaultaxesvalues=['','','','I',], overwrite=True)
    imhead(imagename=regname+'.magma.im')
    os.system('rm -rf '+regname+'magma.lsrk')
    imreframe(regname+'.magma.im',output=regname+'.magma.lsrk',outframe='lsrk') 
    imregrid(regname+'.magma.lsrk',template=prefix+'.image',axes=[0,1,2,3], 
        output=regname+'.magma.rgd',interpolation='cubic',overwrite=True)
    os.system('rm -rf '+regname+'.magma.rgd.pb')
    immath(imagename=[regname+'.magma.rgd', prefix+'.pb'], expr='IM0*IM1',
        outfile=regname+'.magma.rgd.pb')

    # Get conversion factor to Jy/pix
    bmaj=imhead(imagename=regname+'.magma.im',mode='get',hdkey='bmaj')
    bmin=imhead(imagename=regname+'.magma.im',mode='get',hdkey='bmin')
    cd1=qa.convert(imhead(imagename=prefix+'.image',mode='get',hdkey='cdelt1'),'arcsec')
    cd2=qa.convert(imhead(imagename=prefix+'.image',mode='get',hdkey='cdelt2'),'arcsec')
    convfac = (4*math.log(2)*abs(cd1['value']*cd2['value']))/(math.pi*bmaj['value']*bmin['value'])
    print 'The conversion factor from Jy/bm to Jy/pix is ',convfac
    os.system('rm -rf '+regname+'.magma.rgd.jypix')
    immath(regname+'.magma.rgd.pb',outfile=regname+'.magma.rgd.jypix',expr='IM0*'+str(convfac))
    imhead(imagename=regname+'.magma.rgd.jypix',mode='put',hdkey='BUNIT',hdvalue='Jy/pixel')
    imhead(imagename=regname+'.magma.rgd.jypix')

    # Do the feathering (optional)
    if dofeather == True:
        os.system('rm -rf '+prefix+'TPF.image')
        feather(imagename=prefix+'TPF.image',
                highres=prefix+'.image', lowres=regname+'.magma.rgd.pb')
        # Export to FITS:
        exportfits(imagename=prefix+'TPF.image', fitsimage=prefix+'TPF.image.fits',
            dropdeg=True, velocity=True, overwrite=True)

    return
