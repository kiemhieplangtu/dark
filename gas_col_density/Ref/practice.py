import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
nside = 2048
size = 12*nside**2

# test_map = hp.read_map('HFI_CompMap_ThermalDustModel_2048_R1.20.fits')
# print type(test_map)
# print len(test_map)
aa = np.random.rand(size)*np.float(1000.0)
hp.mollview(aa, title='sdfh', min=20, max=500)

#hp.mollview(aa)


plt.show()
#hp.nest2ring()