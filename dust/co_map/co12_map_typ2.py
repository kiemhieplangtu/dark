import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

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

## Read info of 16 CO sources #
 # l,b, nhi, and nhi_error
 #
 # params string fname Filename
 # return dict info
 # 
 # version 11/2016
 # Author Van Hiep ##
def read_info_co(fname = '../result/nh_16src_from_dust.dat'):
	cols = ['idx','src','l','b','nhi','nhi_er','nh2','nh2_er','nhco','nhco_er','nhd','nhd_er']
	fmt  = ['i',  's',  'f','f', 'f',    'f',       's',    's',   'f',  'f',  'f',  'f']
	data = restore(fname, 3, cols, fmt)
	dat  = data.read()
	return dat

#================= MAIN ========================#	
deg2rad  = np.pi/180.
pth      = os.getenv("HOME")+'/hdata/co/'
map_file = pth + 'HFI_CompMap_CO-Type2_2048_R1.10.fits'

info     = read_info_no_co('../../co12/result/26src_no_co_with_sponge.dat')
inco     = read_info_co('../result/nh_16src_from_dust.txt')
dbeam    = 3.5/120.0 # Beam = 3.5' -> dbeam = beam/60/2
#dbeam = 0.1

fukui_cf  = 2.10 #2.10e26
planck_cf = 1.18 #1.18e26

## CO map ##
ir_map = hp.read_map(map_file, field = 0)
nside  = hp.get_nside(ir_map)
res    = hp.nside2resol(nside, arcmin = False)
dd     = res/deg2rad/2.0

#============= For Mask ========= #
offset = 2.0 #degree
# print any(x <0 for x in ir_map) ## No pixel having Negative value 

#====== For Plotting ======#
# hp.mollview(np.log10(ir_map), title=r'$I_{100}$', coord='G', unit='MJy/sr', norm=None) #, min=0.,max=452)
hp.mollview(ir_map, title='CO - High-res type 3', coord='G', unit='Krj km/s', norm='hist', min=-830.,max=4048)

print '============ 26 src No CO ============='
for i in range(0,26):
	src   = info['src'][i]
	l     = info['l'][i]
	b     = info['b'][i]

	theta = (90.0 - b)*deg2rad
	phi   = l*deg2rad
	pix   = hp.ang2pix(nside, theta, phi, nest=False)
	val   = ir_map[pix]
	print src, l,b,pix, val
	hp.projplot(l, b, 'bo', lonlat=True, coord='G')
	# hp.projtext(l, b, src+','+str(l)+','+str(b), lonlat=True, coord='G')
	if (b<60):
		hp.projtext(l, b, src, lonlat=True, coord='G', fontsize=13, weight='bold')
	else:
		hp.projtext(l, b, src, lonlat=True, coord='G', fontsize=13, color='r', weight='bold')

	if(src == '3C109'):
		hp.projtext(l, b, src, lonlat=True, coord='G', fontsize=13, color='b', weight='bold')

print '=========16 src with CO ================'
for i in range(0,16):
	src   = inco['src'][i]
	l     = inco['l'][i]
	b     = inco['b'][i]

	theta = (90.0 - b)*deg2rad
	phi   = l*deg2rad
	pix   = hp.ang2pix(nside, theta, phi, nest=False)
	val   = ir_map[pix]
	print src, l,b,pix, val
	hp.projplot(l, b, 'bo', lonlat=True, coord='G')
	# hp.projtext(l, b, src+','+str(l)+','+str(b), lonlat=True, coord='G')
	hp.projtext(l, b, src, lonlat=True, coord='G', fontsize=13, weight='bold')

mpl.rcParams.update({'font.size':30})
hp.graticule()
plt.grid()
plt.show()

l  = 51.641648
b  = -9.6750019

# Plot cartview a/o mollview #
ll = l
if (l>180):
	ll = ll-360.

offset = 1.
lonr = [51.25, 52.0]
latr = [-10., -9.35]
m    = hp.cartview(ir_map, title=r'$CO$', coord='G', unit='Krj km/s', min=-2.,max=3.,
		norm=None, xsize=800, lonra=lonr, latra=latr, #lonra=[ll-offset,ll+offset], latra=[b-offset,b+offset],
		return_projected_map=True)

# hp.mollview(tau_map, title=info['src'][i]+'('+str(info['l'][i])+','+str(info['b'][i])+') - '+map_file,
# 	coord='G', unit='', rot=[l,b,0], norm='hist', min=1e-7,max=1e-3, xsize=800)

# Cal. #
hp.projplot(ll, b, 'bo', lonlat=True, coord='G')
hp.projtext(ll, b, ' (' + str(round(ll,2)) + ',' + str(round(b,2)) + ')', lonlat=True, coord='G', fontsize=18, weight='bold')

# theta = (90.0-b)*deg2rad
# phi   = l*deg2rad
# pix   = hp.ang2pix(nside, theta, phi, nest=False)

# for x in pl.frange(l-offset, l+offset, dd):
# 	for y in pl.frange(b-offset, b+offset, dd):
# 		cosb = np.cos(b*deg2rad)
# 		cosy = np.cos(y*deg2rad)
# 		if ( (((x-l)**2 + (y-b)**2) <= offset**2) ):
# 			# hp.projtext(x, y, '.', lonlat=True, coord='G')
# 			hp.projplot(x, y, 'kx', lonlat=True, coord='G')

mpl.rcParams.update({'font.size':30})
hp.graticule()
plt.grid()
# plt.imshow(m, origin='lower',extent=(lonr[1],lonr[0],latr[0],latr[1]), interpolation = 'none')
plt.show()

## LOW CO map ##
lowco = ir_map
for i in range(len(lowco)):
	if (lowco[i] > 0.):
		lowco[i] = None

#====== For Plotting ======#
# hp.mollview(np.log10(ir_map), title=r'$I_{100}$', coord='G', unit='MJy/sr', norm=None) #, min=0.,max=452)
hp.mollview(lowco, title='CO - High-res type 3', coord='G', unit='Krj km/s', norm='hist')

print '============ 26 src No CO ============='
for i in range(0,26):
	src   = info['src'][i]
	l     = info['l'][i]
	b     = info['b'][i]

	theta = (90.0 - b)*deg2rad
	phi   = l*deg2rad
	pix   = hp.ang2pix(nside, theta, phi, nest=False)
	# val   = ir_map[pix]
	# print src, l,b,pix, val
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
# plt.imshow(m, origin='lower',extent=(lonr[1],lonr[0],latr[0],latr[1]), interpolation = 'none')
plt.show()