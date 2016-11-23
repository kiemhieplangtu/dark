import matplotlib.pyplot as plt
import numpy             as np
import healpy            as hp
import operator

# Read info of 79 sources #
#
# params string fname Filename
#
# return void
# 
# Author Van Hiep
##
def read_info(fname = '79src_info.txt'):
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

# Find 26 sources with no CO #
#
# params string fname Filename
#
# return void
# 
# Author Van Hiep
##
def read_info_no_co(fname = '79src_info.txt'):
	info = read_info()
	#print map(operator.sub, info['b'], info['bc'])

	j=0
	for i in range(0,79):
		if (info['yn'][i] == 0) :
			print('{0}\t{1}\t\t{2}\t\t{3}\t\t{4}\t\t{5}\t\t{6}\t\t{6}'.format(j, info['src'][i], info['l'][i], info['b'][i], info['ra_icrs'][i], info['de_icrs'][i], info['ra_j'][i], info['de_j'][i]))
			j = j + 1


# Read info of 79 sources #
#
# params string fname Filename
#
# return void
# 
# Author Van Hiep
##
def read_info_no_co(fname = '26src_no_co.dat'):
	ret = {}

	ret['src'] = []
	ret['l']  = []
	ret['b']  = []
	ret['ra_icrs']  = []
	ret['de_icrs']  = []
	ret['ra_j']  = []
	ret['de_j']  = []

	file = open (fname,'r')
	for line in file:
	    line    = line.strip()
	    columns = line.split()

	    ret['src'].append(columns[1])
	    ret['l'].append(float(columns[2]))
	    ret['b'].append(float(columns[3]))
	    ret['ra_icrs'].append(float(columns[4]))
	    ret['de_icrs'].append(float(columns[5]))
	    ret['ra_j'].append(str(columns[6]))
	    ret['de_j'].append(str(columns[7]))

	file.close()

	return ret	

#================= MAIN ========================#	

info = read_info_no_co()

# NSIDE = 2048
fact     = 1.42144524614e-05 
filename = 'HFI_CompMap_ThermalDustModel_2048_R1.20.fits'

test_map = hp.read_map(filename, field = 0)
#hp.mollview(test_map, title=filename, coord='G', unit='K', norm='hist', min=1e-7,max=1e-3, xsize=800)

#hp.mollview(test_map)
tmap = hp.cartview(test_map, title=filename, coord='G', rot=[0,0], unit='K', 
	norm='hist', min=1e-7,max=1e-3, xsize=800, lonra=[-180,-160], latra=[63,90],
	return_projected_map=True, cbar=True)
#hp.orthview(test_map)
#hp.gnomview(test_map)

# print hp.get_nside(test_map)
# print hp.maptype(test_map)
# print hp.get_map_size(test_map)
# print len(test_map)
# print test_map[0:10]*np.float(fact)

# print test_map
n     = 16
rang  = range(0,12*n**2)
npix  = hp.nside2npix(n)
angs  = hp.pix2ang(n,rang)
phi   = angs[0]
theta = angs[1]

nb    = hp.get_all_neighbours(n, 3.1, 888.8, nest=True)

# print npix
# print type(angs)
# print len(angs)
# print nb

print tmap
print len(tmap)
print type(tmap[0])
print tmap[0]
print len(tmap[0])
print np.shape(tmap)

print np.average(tmap)

#plt.loglog(hp.anafast(test_map))

# plt.grid()
# plt.xlabel("$\ell$")
# plt.ylabel("$C_\ell$")

#plt.hist(theta)
plt.grid()
plt.show()