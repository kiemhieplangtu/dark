import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

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
def read_info_no_co_79(fname = '79src_info.txt'):
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
 # return dict infocd 
 # 
 # version 1/2017
 # Author Van Hiep ##
def read_info_no_co(fname = '../../co12/result/26src_no_co_with_sponge.dat'):
	cols = ['idx','src','l','b','ra_icrs','de_icrs','ra_j','de_j', 'oh', 'nhi','nhi_er','thin','thin_er', 'cnm','cnm_er','wnm','wnm_er']
	fmt  = ['i',  's',  'f','f', 'f',    'f',       's',    's',    'i', 'f',   'f',     'f',    'f'    , 'f',   'f',     'f',  'f'    ]
	data = restore(fname, 2, cols, fmt)
	dat  = data.read()
	return dat

#================= MAIN ========================#	
deg2rad  = np.pi/180.
pth      = os.getenv("HOME")+'/hdata/dust/'
map_file = pth + 'HFI_CompMap_ThermalDustModel_2048_R1.20.fits'

info     = read_info_no_co('../../co12/result/26src_no_co_with_sponge.dat')
dbeam    = 3.5/120.0 # Beam = 3.5' -> dbeam = beam/60/2
#dbeam = 0.1

fukui_cf  = 2.10 #2.10e26
planck_cf = 1.18 #1.18e26

tau_map   = hp.read_map(map_file, field = 0)
err_map   = hp.read_map(map_file, field = 1)
nside     = hp.get_nside(tau_map)
res       = hp.nside2resol(nside, arcmin = False)
dd        = res/deg2rad/2.0

#============= For Mask ========= #
offset = 2.0 #degree

#====== For Plotting ======#
hp.mollview(tau_map, title=r'$\tau_{353}$', coord='G', unit='', norm='hist', min=7e-10,max=0.025)

for i in range(0,26):
	src   = info['src'][i]
	l     = info['l'][i]
	b     = info['b'][i]

	theta = (90.0 - b)*deg2rad
	phi   = l*deg2rad
	pix   = hp.ang2pix(nside, theta, phi, nest=False)
	val   = tau_map[pix]
	print src, l,b,pix, val
	hp.projplot(l, b, 'bo', lonlat=True, coord='G')
	# hp.projtext(l, b, src+','+str(l)+','+str(b), lonlat=True, coord='G')
	if (b<60):
		hp.projtext(l, b, src, lonlat=True, coord='G', fontsize=13, weight='bold')
	else:
		hp.projtext(l, b, src, lonlat=True, coord='G', fontsize=13, color='r', weight='bold')

	if(src == '3C109'):
		hp.projtext(l, b, src, lonlat=True, coord='G', fontsize=13, color='b', weight='bold')

mpl.rcParams.update({'font.size':30})
plt.grid()
plt.show()