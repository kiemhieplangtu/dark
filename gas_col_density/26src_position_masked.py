import matplotlib.pyplot as plt
import numpy             as np
import healpy            as hp
import pylab             as pl
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
def read_info_no_co_79(fname = '79src_info.txt'):
	info = read_info()
	#print map(operator.sub, info['b'], info['bc'])

	j=0
	for i in range(0,79):
		if (info['yn'][i] == 0) :
			print('{0}\t{1}\t\t{2}\t\t{3}\t\t{4}\t\t{5}\t\t{6}\t\t{6}'.format(j, info['src'][i], info['l'][i], info['b'][i], info['ra_icrs'][i], info['de_icrs'][i], info['ra_j'][i], info['de_j'][i]))
			j = j + 1


# Read info of 26 sources #
#
# params string fname Filename
#
# return void
# 
# Author Van Hiep
##
def read_info_no_co(fname = '26src_no_co_info.dat'):
	ret = {}

	ret['src']      = []
	ret['l']        = []
	ret['b']        = []
	ret['ra_icrs']  = []
	ret['de_icrs']  = []
	ret['ra_j']     = []
	ret['de_j']     = []

	ret['nhi_heiles'] = []

	file = open (fname,'r')
	file.readline() # comment line
	file.readline() # comment line
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
	    ret['nhi_heiles'].append(float(columns[8]))

	file.close()

	return ret	

#================= MAIN ========================#	
deg2rad  = np.pi/180.
map_file = 'data/HFI_CompMap_ThermalDustModel_2048_R1.20.fits'

info  = read_info_no_co()
dbeam = 3.5/120.0 # Beam = 3.5' -> dbeam = beam/60/2
#dbeam = 0.1

fukui_cf  = 2.10 #2.10e26
planck_cf = 1.18 #1.18e26

test_map = hp.read_map(map_file, field = 0)
err_map  = hp.read_map(map_file, field = 1)
nside    = hp.get_nside(test_map)
res      = hp.nside2resol(nside, arcmin = False)
dd       = res/deg2rad/2.0

#============= For Mask ========= #
offset = 2.0 #degree

for i in range(0,26):

	l     = info['l'][i]
	b     = info['b'][i]

	theta = (90.0 - b)*deg2rad
	phi   = l*deg2rad
	pix   = hp.ang2pix(nside, theta, phi, nest=False)
	test_map[pix] = hp.UNSEEN

	for x in pl.frange(l-offset, l+offset, dd):
		for y in pl.frange(b-offset, b+offset, dd):
			theta = (90.0 - y)*deg2rad
			phi   = x*deg2rad
			pix   = hp.ang2pix(nside, theta, phi, nest=False)
			if (test_map[pix] > -1.0e30) : # Some pixels not defined
				test_map[pix] = hp.UNSEEN


#====== For Plotting ======#
hp.ma(test_map)
hp.mollview(test_map, title='Masked map demo', 
	coord='G', unit='K', norm='hist', 
	min=7e-10,max=0.025, xsize=800)

# hp.cartview(test_map, title=map_file, coord='G', unit='K', rot=[0,0],
# 	norm='hist', min=7e-10,max=0.025, xsize=800)
# 	, lonra=[118.,119.], latra=[-53,-51])

plt.grid()
plt.show()










#==================== An Example for Mask ==========#
# l = 118.6
# b = -52.7

# theta = (90.0 - b)*deg2rad
# phi   = l*deg2rad
# pix   = hp.ang2pix(nside, theta, phi, nest=False)
# test_map[pix] = hp.UNSEEN

# for x in pl.frange(l-offset, l+offset, dd):
# 	for y in pl.frange(b-offset, b+offset, dd):
# 		theta = (90.0 - y)*deg2rad
# 		phi   = x*deg2rad
# 		pix   = hp.ang2pix(nside, theta, phi, nest=False)
# 		#print x, y, test_map[pix]
# 		if (test_map[pix] > -1.0e30) :
# 			test_map[pix] = hp.UNSEEN


#============== Use cartview to average the Tau353 ========#
	# if (sl > 180.) :
	# 	sl = sl - 360.

	# long1 = sl - dbeam
	# long2 = sl + dbeam

	# lat1  = sb - dbeam
	# lat2  = sb + dbeam

	
	# tmap  = hp.cartview(test_map, title=map_file, coord='G', unit='K', rot=[sl,sb],
	# 			norm='hist', min=1e-7,max=1e-3, xsize=800, lonra=[-dbeam,dbeam], latra=[-dbeam,dbeam],
	# 			return_projected_map=True)

	# #print np.shape(tmap)
	# t353 = np.average(tmap)
	# nh   = cf*t353
	
	# print("{:2d}\t{:s}\t\t{:.2e}".format(i, info['src'][i], t353))