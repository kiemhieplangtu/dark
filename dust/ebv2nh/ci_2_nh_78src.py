import sys, os
sys.path.insert(0, r'/home/vnguyen/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import matplotlib        as mpl
import healpy            as hp
import pylab             as pl
import operator

from numpy    import array
from restore  import restore
from plotting import cplot

## Read info of 78 sources #
 # l,b, nhi, and nhi_error
 #
 # params string fname Filename
 # return dict info
 # 
 # version 04/10/2016
 # Author Van Hiep ##
def read_info_78src(fname = '../../hi/rearrange/nhi_lb_78src.txt'):
	cols = ['indx','src','l','b','nhi','nhi_er']
	fmt  = ['i',   's',  'f','f', 'f',  'f']
	data = restore(fname, 2, cols, fmt)
	return data.read()	

## Cal. uncertainty in the mean #
 #
 # params list lst List of numbers
 # return float ret Uncertainty in the mean
 # 
 # Author Van Hiep ##
def cal_uncertainty_in_mean(lst):
	n    = len(lst)
	mean = sum(lst)/float(n)

	s    = 0
	for i in range(n):
		s = s + (lst[i] - mean)**2

	s = np.sqrt(s)
	return s/n

## Read info of 26 sources #
 #
 # params string fname Filename
 #
 # return void
 # 
 # Author Van Hiep ##
def read_info_no_co(fname = '../sub_data/26src_no_co_info.dat'):
	cols = ['idx','src','l','b','ra_icrs','de_icrs','ra_j','de_j','nhi','wnm','cnm','er_nhi','err_nhi','oh', 'thin']
	fmt  = ['i','s','f','f','f','f','s','s','f','f','f','f','f','i','f']
	data = restore(fname, 2, cols, fmt)
	return data.read()

## Read info of 23 LOW NHI sources #
 # l,b, nhi, and nhi_error
 #
 # params string fname Filename
 # return dict info
 # 
 # version 12/2016
 # Author Van Hiep ##
def read_lownhi_23src(fname = '../../hi/result/lownhi_thin_cnm_wnm.txt'):
	cols = ['idx','src','l', 'b', 'nhi','nhi_er','thin','thin_er', 'cnm','cnm_er','wnm','wnm_er']
	fmt  = ['i',  's',  'f', 'f', 'f',   'f',     'f',    'f'    , 'f',   'f',     'f',  'f'    ]
	data = restore(fname, 2, cols, fmt)
	dat  = data.read()
	return dat

## Get tau353 values and err_tau353 values #
 #
 # params str map_file File of maps
 # params dict info Information of sources
 #
 # return void
 # 
 # Author Van Hiep ##
def plot_patches(map_file, info):
	src = info['src']

	# Define constants #
	deg2rad   = np.pi/180.

	# Define the width of area #
	beam   = 30.            # Beam = 30'
	dbeam  = beam/120.0     # Beam = 30' -> dbeam = beam/60/2 in degree
	offset = dbeam          # degree

	# TTYPE1  = 'TAU353  '           / Optical depth at 353GHz                        
	# TTYPE2  = 'ERR_TAU '           / Error on optical depth                         
	# TTYPE3  = 'EBV     '           / E(B-V) color excess                            
	# TTYPE4  = 'RADIANCE'           / Integrated emission                            
	# TTYPE5  = 'TEMP    '           / Dust equilibrium temperature                   
	# TTYPE6  = 'ERR_TEMP'           / Error on T                                     
	# TTYPE7  = 'BETA    '           / Dust emission spectral index                   
	# TTYPE8  = 'ERR_BETA'           / error on Beta  
	ci_map = hp.read_map(map_file, field = 2)
	nside  = hp.get_nside(ci_map)
	res    = hp.nside2resol(nside, arcmin=False)
	dd     = res/deg2rad/2.0

	# OK - Go #
	ebv    = []
	nh     = []
	nhi    = []

	for i in range(0,len(src)):

		# if( i != 15):
		# 	continue

		l  = info['l'][i]
		b  = info['b'][i]

		# Plot cartview a/o mollview #
		ll = l
		if (l>180):
			ll = ll-360.

		# offset = 1.
		hp.cartview(ci_map, title=info['src'][i], coord='G', unit='',
				norm='hist', xsize=800, lonra=[ll-offset-0.1*offset,ll+offset+0.1*offset], latra=[b-offset-0.1*offset,b+offset+0.1*offset],
				return_projected_map=True)

		# hp.mollview(tau_map, title=info['src'][i]+'('+str(info['l'][i])+','+str(info['b'][i])+') - '+map_file,
		# 	coord='G', unit='', rot=[l,b,0], norm='hist', min=1e-7,max=1e-3, xsize=800)

		# Cal. #
		hp.projplot(ll, b, 'bo', lonlat=True, coord='G')
		hp.projtext(ll, b, ' (' + str(round(ll,2)) + ',' + str(round(b,2)) + ')', lonlat=True, coord='G', fontsize=18, weight='bold')

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


## Get tau353 values and err_tau353 values from 16 sources with CO, to compare N(H) derived from CO #
 # Than calculate N(H) from Dust
 #
 # params str map_file File of maps
 # params dict info Information of 26 no CO sources
 # params dict info Information of 23 Low NHI sources
 #
 # return void
 # 
 # Author Van Hiep ##	
def cal_nh_from_dust(map_file, info, lowhi):
	src = info['src']  ## 78 src
	nhi = info['nhi']

	# Define constants #
	deg2rad   = np.pi/180.

	# Define the width of area #
	beam   = 5.0            # Beam = 30'
	dbeam  = beam/120.0     # Beam = 30' -> dbeam = beam/60/2 in degree
	offset = dbeam          # degree

	# TTYPE1  = 'TAU353  '           / Optical depth at 353GHz                        
	# TTYPE2  = 'ERR_TAU '           / Error on optical depth                         
	# TTYPE3  = 'EBV     '           / E(B-V) color excess                            
	# TTYPE4  = 'RADIANCE'           / Integrated emission                            
	# TTYPE5  = 'TEMP    '           / Dust equilibrium temperature                   
	# TTYPE6  = 'ERR_TEMP'           / Error on T                                     
	# TTYPE7  = 'BETA    '           / Dust emission spectral index                   
	# TTYPE8  = 'ERR_BETA'           / error on Beta  
	ci_map = hp.read_map(map_file, field = 2)
	# ci_map = hp.read_map(map_file, field = 0)
	nside  = hp.get_nside(ci_map)
	res    = hp.nside2resol(nside, arcmin=False)
	dd     = res/deg2rad/10.0

	# OK - Go 1: 26 src without CO #
	ebv    = []
	nh     = []

	ci     = {}
	ci_err = {}

	for i in range(0, len(src)):
		# Find the values of Tau353 and Err_tau353 in small area #
		ci[i]   = []

		l = info['l'][i]
		b = info['b'][i]

		# Plot cartview a/o mollview #
		# ll = l
		# if (l>180):
		# 	ll = ll-360.

		# hp.cartview(tau_map, title=info['src'][i]+'('+str(info['l'][i])+','+str(info['b'][i])+') - '+map_file, coord='G', unit='',
		# 		norm='hist', xsize=800, lonra=[ll-offset-0.1*offset,ll+offset+0.1*offset], latra=[b-offset-0.1*offset,b+offset+0.1*offset],
		# 		return_projected_map=True)
		## End Plot cartview a/o mollview ##

		# Cal. #
		theta = (90.0-b)*deg2rad
		phi   = l*deg2rad
		pix   = hp.ang2pix(nside, theta, phi, nest=False)

		if (ci_map[pix] > -0.000001) : # Some pixels not defined
			ci[i].append(ci_map[pix])

		for x in pl.frange(l-offset, l+offset, dd):
			for y in pl.frange(b-offset, b+offset, dd):
				cosb = np.cos(b*deg2rad)
				cosy = np.cos(y*deg2rad)

				if ( ((x-l)**2 + (y-b)**2) <= offset**2 ):
					# hp.projtext(x, y, '.', lonlat=True, coord='G')
					theta = (90.0 - y)*deg2rad
					phi   = x*deg2rad
					pix   = hp.ang2pix(nside, theta, phi, nest=False)

					if (ci_map[pix] > -0.000001) :
						ci[i].append(ci_map[pix])

		# plt.show()
		# continue

		vci = list(set(ci[i]))
		cnt = len(vci)

		# Calculate mean values of tau353 #
		# ebv.append(sum(ci[i])/float(cnt))
		val = sum(vci)/float(cnt)
		err = cal_uncertainty_in_mean(vci)
	   
		# Calculate the NH from E(B-V) #
		# n_h = val/1.44 # 1e22; (NH = 1.e22*EBV/1.44)
		n_h = 0.58*val # 1e22
		err = 0.58*err

		print("{}  {}\t{:08.4f}  {:08.4f}   {}   {}   {}   {}"
			.format(i, src[i],l,b, info['nhi'][i], info['nhi_er'][i],  n_h*100., err*100  ))

		# print src[i], val, nhi[i], n_h*100.
		nh.append(n_h*100.)
		# cols = ['idx','src', 'l', 'b', 'nhi','nhi_er', 'nh','nh_er']

	

	nhi = np.asarray(nhi)
	nh  = np.asarray(nh)
	x   = np.log10(nhi)
	y   = nh/nhi 
	plt.plot(x,y, 'r.', label='Ratio $f = N_{H}$/$N_{HI}$')
	plt.title('Correlation between $f = N_{H}$/$N_{HI}$ and HI column density $N_{HI}$ \nalong 26 lines-of-sight without the presence of CO line', fontsize=30)
	plt.ylabel('$Ratio f = N_{H}$/$N_{HI}$', fontsize=35)
	plt.xlabel('log$_{10}$($N_{HI}/10^{20}$ cm$^{-2}$)', fontsize=35)
	plt.xlim(0, 1.6)
	plt.ylim(0, 3)
	plt.grid(True)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=18)

	plt.text(0.2, 0.31, '(Available sources with the presence of OH are shown)', color='k', fontsize=17)
	# plt.text(0.2, 0.4, '$f = ['+str(m)+'\pm'+str(ea) +']\cdot log_{10}(N_{HI}/10^{20})+['+str(b)+'\pm'+str(eb)+']$', fontsize=20 )

	plt.legend(loc='upper left', fontsize=18)
	# plt.savefig("test.png",bbox_inches='tight')
	for i in range(len(src)):
		# if (oh[i] > 0) :
		plt.annotate('('+str(src[i])+')', xy=(x[i], y[i]), xycoords='data',
               xytext=(-50.,30.), textcoords='offset points',
               arrowprops=dict(arrowstyle="->"),fontsize=18,
               )
	plt.show()

	plt.plot(nhi,nh, 'rd', label='data no CO', ms=10)
	plt.plot([0,130],[0,130], 'k--', label='$N_{H} = N_{HI}$')
	plt.title('Correlation between $N_{H}$ and $N_{HI}$ \nalong 78 lines-of-sight', fontsize=30)
	plt.ylabel('$N_{H}[10^{20}$ cm$^{-2}]$', fontsize=35)
	plt.xlabel('$N_{HI} [10^{20}$ cm$^{-2}]$', fontsize=35)
	# plt.xlim(0, 1.6)
	# plt.ylim(0, 3)
	plt.grid(True)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=18)

	plt.text(15., 2., '(Available sources with the presence of OH are shown)', color='k', fontsize=17)
	plt.text(15., 3., r'$N_{H} = 5.8\cdot10^{21}[cm^{-2}mag^{-1}]\cdot E(B-V)$', color='k', fontsize=17)
	plt.text(15., 4., r'E(B-V) from Planck data R1.2', color='k', fontsize=17)

	plt.legend(loc='upper left', fontsize=18)
	# plt.savefig("test.png",bbox_inches='tight')
	for i in range(len(src)):
		# if (oh[i] > 0) :
		plt.annotate('('+str(src[i])+')', xy=(nhi[i], nh[i]), xycoords='data',
               xytext=(-50.,30.), textcoords='offset points',
               arrowprops=dict(arrowstyle="->"),fontsize=18,
               )
	plt.show()

#================= MAIN ========================#
pth      = os.getenv("HOME")+'/hdata/dust/'
map_file = pth + 'HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
map_file = pth + 'HFI_CompMap_DustOpacity_2048_R1.10.fits'
# map_file = pth + 'lambda_sfd_ebv.fits'  ## E(B-V) from SFD et al. 1998

## Infor of 78 src without CO && 23 Low NHI sources ##
info     = read_info_78src(fname = '../../hi/rearrange/nhi_lb_78src.txt')
lowhi    = read_lownhi_23src(fname = '../../hi/result/lownhi_thin_cnm_wnm.txt')

cal_nh_from_dust(map_file, info, lowhi)
# plot_patches(map_file, info)