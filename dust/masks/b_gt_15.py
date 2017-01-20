import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class
import matplotlib.pyplot as plt
import matplotlib        as mpl
import numpy             as np
import healpy            as hp

import operator
from restore             import restore

## Define constants ##
deg2rad = np.pi/180.
beam    = 3.5
dbeam   = beam/120.0 # Beam = 3.5' -> dbeam = beam/60/2
dbeam   = 0.5
edge    = dbeam/40. # To plot a rectangle about the source

# lowest  = hp.read_map(pth+'mask_Lowest1perc_ns256.fits', field = 0, h=False)
# low_hi  = hp.read_map(pth+'mask_LowNHI_ns256.fits',      field = 0, h=False)
# sth_cap = hp.read_map(pth+'mask_SouthCap_ns256.fits',    field = 0, h=False)
# gt15    = hp.read_map(pth+'mask_GLATgt15_ns256.fits',    field = 0, h=False)
# whole   = hp.read_map(pth+'mask_WholeSky_ns256.fits',    field = 0, h=False)

pth     = os.getenv("HOME")+'/hdata/dust/'
gt15    = hp.read_map(pth+'mask_Lowest1perc_ns256.fits', field = 0, h=False)
msk     = hp.ud_grade(gt15, nside_out=2048)  # weight = 2

## Color map ##
nside = 2048
cmap  = plt.cm.get_cmap('spring')
cmap  = mpl.colors.ListedColormap(['lightgray', 'firebrick'])
hp.mollview(msk, title='Lowest 1%', coord='G', unit='', rot=[0,0,0], norm=None, xsize=800, cmap=cmap)

hp.graticule()
# hp.write_map("planck_mask.fits", msk)
plt.show()