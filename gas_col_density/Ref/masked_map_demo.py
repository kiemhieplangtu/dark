import matplotlib.pyplot as plt
import numpy             as np
import healpy            as hp

m = np.arange(12.)
m[3] = hp.UNSEEN
m[11] = hp.UNSEEN
m[4] = hp.UNSEEN

hp.ma(m)
hp.mollview(m, title='Masked map demo', 
	coord='G', unit='K', norm='hist', 
	min=0,max=12, xsize=800)
plt.grid()
plt.show()