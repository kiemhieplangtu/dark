import numpy as np
from matplotlib import pyplot as plt

# warning: due to a bug in healpy, importing it before pylab can cause
#  a segmentation fault in some circumstances.
import healpy as hp

from astroML.datasets import fetch_wmap_temperatures

#------------------------------------------------------------
# Fetch the wmap data
wmap_unmasked = fetch_wmap_temperatures(masked=False)

#------------------------------------------------------------
# plot the unmasked map
fig = plt.figure(1)
hp.mollview(wmap_unmasked, min=-1, max=1, title='Raw WMAP data',
            fig=1, cmap=plt.cm.jet, unit=r'$\Delta$T (mK)')
plt.show()
