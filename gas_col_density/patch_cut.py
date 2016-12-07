import sys, os
sys.path.insert(0, r'/home/vnguyen/dark/common') # add folder of Class
import matplotlib.pyplot as plt
import matplotlib        as mpl
import numpy             as np
import healpy            as hp
import pylab             as pl

import operator
from restore             import restore

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


## Read info of 26 no-CO sources #
 # l,b, nhi, and nhi_error
 #
 # params string fname Filename
 # return dict info
 # 
 # version 11/2016
 # Author Van Hiep ##
def read_info_no_co(fname = '../26src_no_co_info.dat'):
	cols = ['idx','src','l','b','ra_icrs','de_icrs','ra_j','de_j','nhi','wnm','cnm','e_nhi','nhi_er','oh']
	fmt  = ['i',  's',  'f','f', 'f',    'f',       's',    's',   'f',  'f',  'f',  'f',    'f',    'f']
	data = restore(fname, 2, cols, fmt)
	dat  = data.read()
	return dat

##================= MAIN ========================##
pth      = os.getenv("HOME")+'/hdata/dust/'
map_file = pth + 'HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
quantity = 'tau353'
map_unit = ''

# map_file = 'data/HFI_SkyMap_353_2048_R2.02_full.fits'
# quantity = 'w353'
# map_unit = 'K_CMB'


info  = read_info_no_co('26src_no_co_info.dat')

# Define the width of area #
beam   = 3.5            # Beam = 3.5'
dbeam  = beam/120.0     # Beam = 3.5' -> dbeam = beam/60/2 in degree
offset = dbeam          # degree

# Define constants #
nside     = 2048
deg2rad   = np.pi/180.

tau_map   = hp.read_map(map_file, field = 0)
nside     = hp.get_nside(tau_map)
res       = hp.nside2resol(nside, arcmin=False)
dd        = res/deg2rad/15.0

for i in range(0,26):

	l  = info['l'][i]
	b  = info['b'][i]

	# Plot cartview a/o mollview #
	ll = l
	if (l>180):
		ll = ll-360.

	hp.cartview(tau_map, title=info['src'][i], coord='G', unit='',
			norm=None, xsize=800, lonra=[ll-offset-0.1*offset,ll+offset+0.1*offset], latra=[b-offset-0.1*offset,b+offset+0.1*offset],
			return_projected_map=True)

	hp.projplot(ll, b, 'bo', lonlat=True, coord='G')
	hp.projtext(ll, b, ' (' + str(round(ll,2)) + ',' + str(round(b,2)) + ')', lonlat=True, coord='G', fontsize=18, weight='bold')

	# hp.mollview(tau_map, title=info['src'][i]+'('+str(info['l'][i])+','+str(info['b'][i])+') - '+map_file,
	# 	coord='G', unit='', rot=[l,b,0], norm='hist', min=1e-7,max=1e-3, xsize=800)

	# Cal. #
	theta = (90.0-b)*deg2rad
	phi   = l*deg2rad
	pix   = hp.ang2pix(nside, theta, phi, nest=False)

	for x in pl.frange(l-offset, l+offset, dd):
		for y in pl.frange(b-offset, b+offset, dd):
			cosb = np.cos(b*deg2rad)
			cosy = np.cos(y*deg2rad)
			if ( (((x-l)**2 + (y-b)**2) <= offset**2) ):
				# hp.projtext(x, y, '.', lonlat=True, coord='G')
				hp.projplot(x, y, 'kx', lonlat=True, coord='G')

	mpl.rcParams.update({'font.size':30})
	plt.show()