import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import pyfits
m = hp.read_map('HFI_CompMap_ThermalDustModel_2048_R1.20.fits')
#m.mask = np.isnan(m)
hp.mollview(m, min=-1e-5, max=1e-5, xsize=2000)
plt.title("gll_iem_v02_p6_V11_DIFFUSE")
plt.loglog(hp.anafast(m))
plt.show()
