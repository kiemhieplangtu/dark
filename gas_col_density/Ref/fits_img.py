import sys, os
sys.path.insert(0, r'/home/vnguyen/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import matplotlib        as mpl
import healpy            as hp
import pylab             as pl
import operator

from   numpy             import array
from   restore           import restore
from   astropy.io        import fits

#================= MAIN ========================#
pth      = os.getenv("HOME")+'/hdata/hi/galfa/'
imgfile  = pth + 'GALFA_HI_RA+DEC_324.00+10.35_N.fits'
imgfile  = pth + 'GALFA_HI_RA+DEC_324.00+18.35_N.fits'
hdu_list = fits.open(imgfile)
hdu_list.info()
image_data = hdu_list[0].data
print type(image_data)
print hdu_list[0].header
print(image_data.shape)
plt.imshow(image_data[100,:,:], origin='lower')
plt.colorbar()
plt.show()

print image_data[:,0,0]

plt.plot(image_data[:,0,0])
plt.show()