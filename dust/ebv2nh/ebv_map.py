import sys, os
sys.path.insert(0, r'/home/vnguyen/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import matplotlib        as mpl
import numpy             as np
import healpy            as hp
import pylab             as pl
import operator
from   restore           import restore

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
 # return dict info
 # 
 # version 11/2016
 # Author Van Hiep ##
def read_info_no_co(fname = '../sub_data/26src_no_co_info.dat'):
	cols = ['idx','src','l','b','ra_icrs','de_icrs','ra_j','de_j','nhi','wnm','cnm','e_nhi','nhi_er','oh']
	fmt  = ['i',  's',  'f','f', 'f',    'f',       's',    's',   'f',  'f',  'f',  'f',    'f',    'f']
	data = restore(fname, 2, cols, fmt)
	dat  = data.read()
	return dat

#================= MAIN ========================#	
deg2rad  = np.pi/180.
pth      = os.getenv("HOME")+'/hdata/dust/'
map_file = pth + 'HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
# map_file = pth + 'lambda_sfd_ebv.fits'
# map_file = pth + 'lambda_green_dust_map_2d.fits'

info     = read_info_no_co('../sub_data/26src_no_co_info.dat')
dbeam    = 3.5/120.0 # Beam = 3.5' -> dbeam = beam/60/2
#dbeam = 0.1

fukui_cf  = 2.10 #2.10e26
planck_cf = 1.18 #1.18e26

## E(B-V) map ##
ci_map = hp.read_map(map_file, verbose=False, field = 2)
nside  = hp.get_nside(ci_map)
res    = hp.nside2resol(nside, arcmin = False)
dd     = res/deg2rad/2.0

#============= For Mask ========= #
offset = 2.0 #degree

#====== For Plotting ======#
hp.mollview(ci_map, title=r'$E(B-V)$', coord='G', unit='mag', norm=None) #, min=0.,max=452)

for i in range(0,26):
	src   = info['src'][i]
	l     = info['l'][i]
	b     = info['b'][i]

	theta = (90.0 - b)*deg2rad
	phi   = l*deg2rad
	pix   = hp.ang2pix(nside, theta, phi, nest=False)
	val   = ci_map[pix]
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
hp.graticule()
plt.grid()
plt.show()