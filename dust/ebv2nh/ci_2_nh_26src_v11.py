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

## Read info of 26 no-CO sources #
 # l,b, nhi, and nhi_error
 #
 # params string fname Filename
 # return dict info
 # 
 # version 1/2017
 # Author Van Hiep ##
def read_info_no_co(fname = '../../co12/result/26src_no_co_with_sponge.dat'):
	cols = ['idx','src','l','b','ra_icrs','de_icrs','ra_j','de_j', 'oh', 'nhi','nhi_er','thin','thin_er', 'cnm','cnm_er','wnm','wnm_er']
	fmt  = ['i',  's',  'f','f', 'f',    'f',       's',    's',    'i', 'f',   'f',     'f',    'f'    , 'f',   'f',     'f',  'f'    ]
	data = restore(fname, 2, cols, fmt)
	dat  = data.read()
	return dat

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

## Cal NH from E(B-V), Planck R1.2 #
 #
 # params array ebv_map  Map of E(B-V)
 # params array err_map  Error Map of E(B-V)
 # params dict  info     Infor of the sources
 #
 # return list info
 # 
 # version 12/2016
 # Author Van Hiep ##
def nh_from_ebv(ebv_map, err_map, info):
	# Define constants #
	deg2rad = np.pi/180.

	## 26 sources without CO
	src = info['src']  ## 26 src without CO	

	# Define the width of area #
	beam   = 5.             # Beam = 5'
	dbeam  = beam/120.0     # Beam = 5' -> dbeam = beam/60/2 in degree
	offset = dbeam          # degree

	nside  = hp.get_nside(ebv_map)
	res    = hp.nside2resol(nside, arcmin=False)
	dd     = res/deg2rad/10.0

	# OK - Go #
	ebv    = []
	nh     = []
	nher   = []

	ci     = {}
	ci_err = {}

	for i in range(0, len(src)):
		# Find the values of ebv353 and Err_ebv353 in small area #
		ci[i]     = []
		ci_err[i] = []

		l = info['l'][i]
		b = info['b'][i]

		# Cal. #
		theta = (90.0-b)*deg2rad
		phi   = l*deg2rad
		pix   = hp.ang2pix(nside, theta, phi, nest=False)

		if (ebv_map[pix] > -0.000001) : # Some pixels not defined
			ci[i].append(ebv_map[pix])
			ci_err[i].append(err_map[pix])

		for x in pl.frange(l-offset, l+offset, dd):
			for y in pl.frange(b-offset, b+offset, dd):
				cosb = np.cos(b*deg2rad)
				cosy = np.cos(y*deg2rad)

				if ( ((x-l)**2 + (y-b)**2) <= offset**2 ):
					theta = (90.0 - y)*deg2rad
					phi   = x*deg2rad
					pix   = hp.ang2pix(nside, theta, phi, nest=False)

					if (ebv_map[pix] > -0.000001) :
						ci[i].append(ebv_map[pix])
						ci_err[i].append(err_map[pix])

		temp_ebv  = list(set(ci[i]))
		ci_err[i] = list(set(ci_err[i]))

		cnt_ebv = len(temp_ebv)
		npix    = len(ci_err[i])

		if (cnt_ebv != npix) :
			tmp_ebv = []
			for k in range(npix):
				tmp_ebv.append(ci[i][k])

			ci[i] = tmp_ebv
		else:
			ci[i] = list(set(ci[i]))

		# Calculate mean values of EBV #
		val = sum(ci[i])/float(npix)
		# Uncertainties of mean values of EBV #
		sd = 0.
		for j in range(npix):
			sd = sd + (ci_err[i][j])**2

		sd = (sd**0.5)/npix # Uncertainty of a Sum

		# Calculate the NH from E(B-V) #
		# n_h = val/1.44 # 1e22; (NH = 1.e22*EBV/1.44)
		n_h = 0.58*val # 1e22
		err = 0.58*sd

		nh.append(n_h*100.)
		nher.append(err*100.)

		# print("{}  {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}"
		# 	.format(i, src[i], round((nhfk[i]), 4), round((fukui_sd), 4), round((nh[i]), 4), round((planck_sd), 4), \
		# 		nhi[i], nhier[i], wnm, cnm, rat1, rat2, oh[i], thin[i] ) )

	return nh, nher

## Get NH from E(B-V) #
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
	beam   = 5.            # Beam = 30'
	dbeam  = beam/120.0     # Beam = 30' -> dbeam = beam/60/2 in degree
	offset = dbeam          # degree

	# TTYPE1  = 'ci353  '           / Optical depth at 353GHz                        
	# TTYPE2  = 'ERR_ci '           / Error on optical depth                         
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

		if( src[i] != '3C345'):
			continue

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

		# hp.mollview(ci_map, title=info['src'][i]+'('+str(info['l'][i])+','+str(info['b'][i])+') - '+map_file,
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


## Get ebv values and err_ebv values from 26 sources without CO#
 # Than calculate N(H) from Dust
 #
 # params str map_file File of maps
 # params dict info   Information of 26 no CO sources
 # params dict lownhi Information of 23 Low NHI sources
 #
 # return void
 # 
 # version 12/2016
 # Author Van Hiep ##	
def cal_nh_from_dust(map_file, info, lownhi):
	## 26 sources without CO
	src   = info['src']  ## 26 src without CO
	nhi   = info['nhi']
	thin  = info['thin']
	thinr = info['thin_er']
	oh    = info['oh']
	nhier = info['nhi_er']

	## 23 lownhi sources
	hi     = lownhi['nhi']
	hier   = lownhi['nhi_er']
	lthin  = lownhi['thin']
	lthinr = lownhi['thin_er']

	# TTYPE1  = 'ci353  '           / opacity 353GHz                                 
	# TTYPE2  = 'ci353ERR'          / Error on opacity                               
	# TTYPE3  = 'EBV     '           / E(B-V)                                         
	# TTYPE4  = 'EBV_ERR '           / Error on E(B-V)                                
	# TTYPE5  = 'T_HF    '           / T for high freq correction                     
	# TTYPE6  = 'T_HF_ERR'           / Error on T                                     
	# TTYPE7  = 'BETAHF  '           / Beta for high freq correction                  
	# TTYPE8  = 'BETAHFERR'          / Error on Beta  
	ci_map     = hp.read_map(map_file, field = 2)
	er_map     = hp.read_map(map_file, field = 3)
	nh, nher   = nh_from_ebv(ci_map, er_map, info)
	xnh, lnher = nh_from_ebv(ci_map, er_map, lowhi)	

	nhi = np.asarray(nhi)
	nh  = np.asarray(nh)
	x   = np.log10(nhi)
	y   = nh/nhi 
	plt.plot(x,y, 'r.', label='Ratio $f = N_{H}$/$N_{HI}$')
	plt.title('$f = N_{H}$/$N_{HI}$ and HI column density $N_{HI}$ \nalong 26 lines-of-sight without the presence of CO line', fontsize=30)
	plt.ylabel('$Ratio f = N_{H}$/$N_{HI}$', fontsize=35)
	plt.xlabel('log$_{10}$($N_{HI}/10^{20}$ cm$^{-2}$)', fontsize=35)
	plt.xlim(0, 1.6)
	plt.ylim(0, 3)
	plt.grid(True)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=18)

	plt.text(0.2, 0.31, '(Available sources with the presence of OH are shown)', color='k', fontsize=17)

	plt.legend(loc='upper left', fontsize=18)
	for i in range(len(src)):
		if (oh[i] > 0) :
			plt.annotate('('+str(src[i])+')', xy=(x[i], y[i]), xycoords='data',
	               xytext=(-50.,30.), textcoords='offset points',
	               arrowprops=dict(arrowstyle="->"),fontsize=18,
	               )
	plt.show()

	plt.errorbar(nhi, nh,xerr=nhier, yerr=nher, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data without CO')
	plt.errorbar(hi,xnh,xerr=hier, yerr=lnher, color='b', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data low $N_{HI}, N_{HI}<3.0\cdot10^{20} cm^{-2}$')
	plt.plot([0,40],[0,40], 'k--', label=r'$N_{H} = N_{HI}, N_{H} = 5.8\cdot10^{21}[cm^{-2}mag^{-1}]\cdot E(B-V), R1.1$')
	plt.title('$N_{H}$ vs $N_{HI}$ along 26 lines-of-sight without CO & 23 LOS with low $N_{HI}$', fontsize=30)
	plt.ylabel('$N_{H}[10^{20}$ cm$^{-2}]$', fontsize=35)
	plt.xlabel('$N_{HI} [10^{20}$ cm$^{-2}]$', fontsize=35)
	# plt.xlim(0, 1.6)
	# plt.ylim(0, 3)
	plt.grid(True)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=18)

	plt.text(15., 2., '(Available sources with the presence of OH are shown)', color='k', fontsize=17)
	plt.text(15., 3., r'$N_{H} = 5.8\cdot10^{21}[cm^{-2}mag^{-1}]\cdot E(B-V)$', color='k', fontsize=17)
	plt.text(15., 4., r'E(B-V) from Planck data R1.1', color='k', fontsize=17)

	plt.legend(loc='upper left', fontsize=18)
	# plt.savefig("test.png",bbox_inches='tight')
	for i in range(len(src)):
		if (oh[i] > 0) :
			plt.annotate('('+str(src[i])+')', xy=(nhi[i], nh[i]), xycoords='data',
	               xytext=(-50.,30.), textcoords='offset points',
	               arrowprops=dict(arrowstyle="->"),fontsize=18,
	               )

	# for i in range(len(lsc)):
	# 	plt.annotate('('+str(lsc[i])+')', xy=(lhi[i], xnh[i]), xycoords='data',
 #               xytext=(-50.,30.), textcoords='offset points',
 #               arrowprops=dict(arrowstyle="->"),fontsize=18,
 #               )
	plt.show()

	## NH vs Thin, LAB data => Thin ##
	plt.errorbar(thin, nh,xerr=thinr, yerr=nher, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data without CO')
	plt.errorbar(lthin, xnh,xerr=lthinr, yerr=lnher, color='b', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data with low $N_{HI},  N_{HI}<3.0\cdot10^{20} cm^{-2}$')
	plt.plot([0,40],[0,40], 'k--', label=r'$N_{H} = N_{HI}, N_{H} = 5.8\cdot10^{21}[cm^{-2}mag^{-1}]\cdot E(B-V), R1.1$')
	plt.title('$N_{H}$ vs $N^{thin}_{HI}$ along 26 lines-of-sight without CO & 23 LOS with low $N_{HI}$', fontsize=30)
	plt.ylabel('$N_{H}[10^{20}$ cm$^{-2}]$', fontsize=35)
	plt.xlabel('$N^{thin}_{HI} [10^{20}$ cm$^{-2}]$', fontsize=35)
	plt.xlim(-1.0, 30.0)
	# plt.ylim(0, 3)
	plt.grid(True)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=18)

	plt.text(15., 2., '(Available sources with the presence of OH are shown)', color='k', fontsize=17)
	plt.text(15., 3., r'$N_{H} = 5.8\cdot10^{21}[cm^{-2}mag^{-1}]\cdot E(B-V)$', color='k', fontsize=17)
	plt.text(15., 4., r'E(B-V) from Planck data R1.1', color='k', fontsize=17)

	plt.legend(loc='upper left', fontsize=18)
	# plt.savefig("test.png",bbox_inches='tight')
	for i in range(len(src)):
		if (oh[i] > 0) :
			plt.annotate('('+str(src[i])+')', xy=(thin[i], nh[i]), xycoords='data',
	               xytext=(-50.,30.), textcoords='offset points',
	               arrowprops=dict(arrowstyle="->"),fontsize=15,
	               )
	plt.show()

#================= MAIN ========================#
pth      = os.getenv("HOME")+'/hdata/dust/'
# map_file = pth + 'HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
map_file = pth + 'HFI_CompMap_DustOpacity_2048_R1.10.fits'
# map_file = pth + 'lambda_sfd_ebv.fits'  ## E(B-V) from SFD et al. 1998, IRAS ~5'

## Infor of 26 src without CO && 23 Low NHI sources ##
info     = read_info_no_co('../../co12/result/26src_no_co_with_sponge.dat')
lowhi    = read_lownhi_23src(fname = '../../hi/result/lownhi_thin_cnm_wnm.txt')

cal_nh_from_dust(map_file, info, lowhi)
# plot_patches(map_file, lowhi)