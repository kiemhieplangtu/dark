import matplotlib.pyplot as plt
import numpy as np
import healpy as hp

## MAP - FILES ##
## LFI_SkyMap_070_1024_R2.01_survey-8.fits
## HFI_CompMap_ThermalDustModel_2048_R1.20.fits
## HFI_CompMap_CO-Type1_2048_R1.10.fits
## LFI_SkyMap_070_1024_R2.01_survey-8.fits
## HFI_CompMap_CO-Type2_2048_R1.10.fits

# NSIDE = 32
fact     = 1.42144524614e-05 
filename = 'data/HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
filename = 'data/HFI_SkyMap_353_2048_R2.02_full.fits'

test_map = hp.read_map(filename)
#hp.mollview(test_map, title=filename, coord='G', unit='K', norm='hist', min=1e-7,max=1e-3, xsize=800)
hp.mollview(test_map, title=filename, coord='G', unit='K', norm='hist', xsize=800)

#hp.mollview(test_map)
#hp.cartview(test_map, title=filename, coord='G', rot=[0,0], unit='K', norm='hist', min=-1375,max=2687, xsize=800, lonra=[-1,1], latra=[-1,1])
#hp.orthview(test_map)
#hp.gnomview(test_map)

print hp.get_nside(test_map)
print hp.maptype(test_map)
print hp.get_map_size(test_map)
print len(test_map)
print test_map[0:10]*np.float(fact)

equateur_lon = [10.,0.]
equateur_lat = [10.,0.]
hp.projplot(equateur_lon, equateur_lat, lonlat=True, coord='G')

# plt.loglog(hp.anafast(test_map))
# plt.grid()
# plt.xlabel("$\ell$")
# plt.ylabel("$C_\ell$")

plt.grid()
plt.show()
