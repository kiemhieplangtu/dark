import matplotlib.pyplot as plt
import numpy             as np
import healpy            as hp
from   astropy.io        import fits
import operator

# Read info of 79 sources #
#
# params string fname Filename
#
# return void
# 
# Author Van Hiep
##
def read_info(fname = ''):
	ret = {}

	ret['src'] = []
	ret['yn']  = []
	ret['l']  = []
	ret['b']  = []
	ret['ra_icrs']  = []
	ret['de_icrs']  = []
	ret['ra_j']  = []
	ret['de_j']  = []

	file = open (fname,'r')
	file.readline() # comment line
	file.readline() # comment line

	for line in file:
	    line    = line.strip()
	    columns = line.split()

	    ret['src'].append(columns[1])
	    ret['yn'].append(int(columns[2]))
	    ret['l'].append(float(columns[3]))
	    ret['b'].append(float(columns[4]))
	    ret['ra_icrs'].append(float(columns[5]))
	    ret['de_icrs'].append(float(columns[6]))
	    ret['ra_j'].append(str(columns[7]))
	    ret['de_j'].append(str(columns[8]))

	file.close()

	return ret

#================= MAIN ========================#
image_file = 'data/HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
hdu_list = fits.open(image_file)
hdu_list.info()
image_data = hdu_list[1].data
print type(image_data)
print(image_data.shape)
plt.imshow(image_data, cmap='gray')
plt.colorbar()