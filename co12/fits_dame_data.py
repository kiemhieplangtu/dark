import numpy             as np
import matplotlib.pyplot as plt
from   astropy.io        import fits
############## MAIN ###########

hdulist = fits.open('/home/vnguyen/Downloads/COGAL_all_interp_231.fits')
print hdulist[0].header
data = hdulist[0].data
print data.shape

T = data[:, 199, 0] # z y x
dv = 1.3

wco = sum(T[245:257])*dv
wco = round(wco, 2)

print wco

# plt.imshow(data[3,:,:], origin='lower')
plt.plot(T)
plt.show()


