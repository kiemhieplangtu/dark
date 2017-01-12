import os, sys
sys.path.insert(0, r'/home/vnguyen/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import healpy            as hp
import pylab             as pl
import operator
import matplotlib.cm     as cm

from restore             import restore
from masks               import masks

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

# Plot the relation between CNM and WNM #
#
# params dict info Infor of 26 sources read from file 26_src_no_co.txt
#
# return void
# 
# Author Van Hiep
##
def warm_cold_relation(info):
	nhi  = info['nhi_heiles']
	warm = info['nhi_warm']
	cold = info['nhi_cold']

	x =[]
	for i in range(0,26):
		x.append((warm[i] - cold[i])/warm[i])

	a = np.array(x)

	bins=np.histogram(np.hstack((a)), bins=144)[1] #get the bin edges

	plt.hist(a, bins, label='(NHI_Warm - NHI_Cold)/NHI_Warm')

	plt.xlabel('N(HI) diff ratio - [Warm-Cold]/Warm')
	plt.ylabel('Counts')
	plt.title('Histogram of (NHI_Warm - NHI_Cold)/NHI_Warm')
	plt.grid(True)

	plt.legend(loc='upper right')
	plt.xlim(-7, 2)
	plt.show()

## Cal NH from tau353 #
 #
 # params array tau_map  Map of tau353
 # params array err_map  Error map of tau353
 # params array msk      Mask to find sigma353
 # params dict  info     Infor of the sources
 #
 # return list info
 # 
 # version 12/2016
 # Author Van Hiep ##
def nh_from_tau353(tau_map, err_map, msk, info):
	## Classes ##
	msks = masks()

	## 26 sources without CO
	src   = info['src']  ## 26 src without CO
	nhi   = info['nhi']
	thin  = info['thin']
	nhier = info['nhi_er']
	oh    = info['oh']

	# Define constants #
	deg2rad     = np.pi/180.
	fukui_cf    = 2.10e26
	fk_fact_err = 0.0 #unknown

	# Define the width of area #
	beam   = 5.             # Beam = 5'
	dbeam  = beam/120.0     # Beam = 5' -> dbeam = beam/60/2 in degree
	offset = dbeam          # degree

	nside  = hp.get_nside(tau_map)
	res    = hp.nside2resol(nside, arcmin=False)
	dd     = res/deg2rad/10.0

	# OK - Go #
	tau353 = []
	nh     = []
	nher   = []
	nhfk   = []
	nhfker = []

	tau    = {}
	t_err  = {}

	rat_fk = 0.
	rat_pl = 0.

	for i in range(len(src)):
		# Find the values of Tau353 and Err_tau353 in small area #
		tau[i]   = []
		t_err[i] = []

		l = info['l'][i]
		b = info['b'][i]

		## find Planck Conversion Factor (Dust opacity and Its error) ## 
		planck_cf,pl_fact_err = msks.get_dust_opacity(msk,l,b) # 0.84 #8.40e-27, 0.84e-26

		# Cal. #
		theta = (90.0-b)*deg2rad
		phi   = l*deg2rad
		pix   = hp.ang2pix(nside, theta, phi, nest=False)

		if ( (err_map[pix] >= 6.9e-11) and (err_map[pix] <= 0.00081) and (tau_map[pix] > -1.0e30) ): # Checked Invalid error & Some pixels not defined
			tau[i].append(tau_map[pix])
			t_err[i].append(err_map[pix])

		for x in pl.frange(l-offset, l+offset, dd):
			for y in pl.frange(b-offset, b+offset, dd):
				cosb = np.cos(b*deg2rad)
				cosy = np.cos(y*deg2rad)

				if ( ((x-l)**2 + (y-b)**2) <= offset**2 ):
					# hp.projplot(x, y, 'kx', lonlat=True, coord='G')
					theta = (90.0 - y)*deg2rad
					phi   = x*deg2rad
					pix   = hp.ang2pix(nside, theta, phi, nest=False)

					if ( (err_map[pix] >= 6.9e-11) and (err_map[pix] <= 0.00081) and (tau_map[pix] > -1.0e30) ): # Checked Invalid error & Some pixels not defined
						tau[i].append(tau_map[pix])
						t_err[i].append(err_map[pix])

		temp_tau = list(set(tau[i]))
		t_err[i] = list(set(t_err[i]))

		cnt_tau = len(temp_tau)
		npix    = len(t_err[i])

		if (cnt_tau != npix) :
			tmp_tau = []
			for k in range(npix):
				tmp_tau.append(tau[i][k])

			tau[i] = tmp_tau
		else:
			tau[i] = list(set(tau[i]))

		# Calculate mean values of tau353 #
		tau353.append(sum(tau[i])/float(npix))

		# Calculate the N(HI) from Fukui factor #
		nhfk_i = fukui_cf*tau353[i]
		nhfk.append(nhfk_i)
	   
		# Calculate the NH from Planck factor #
		nh_i = tau353[i]/planck_cf
		nh.append(nh_i)

		# Uncertainties of mean values of tau353 #
		sd1_tau353 = 0.
		sd2_tau353 = 0.
		for j in range(npix):
			sd1_tau353 = sd1_tau353 + (tau[i][j]-tau353[i])**2
			sd2_tau353 = sd2_tau353 + (t_err[i][j])**2

		sd1_tau353 = (sd1_tau353/npix)**0.5 # Just to test, nothing to do with this
		sd2_tau353 = (sd2_tau353**0.5)/npix # Uncertainty of a Sum

		# Uncertainties of mean values of N(HI) and N(H) #
		fukui_sd  = (sd2_tau353*fukui_cf)
		planck_sd = (sd2_tau353/tau353[i])**2 + (pl_fact_err/planck_cf)**2
		planck_sd = (nh_i*planck_sd**0.5)

		nhfk[i]   = nhfk[i]*1.0e-20
		fukui_sd  = fukui_sd*1.0e-20
		nh[i]     = nh[i]*1.0e-20
		planck_sd = planck_sd*1.0e-20
		nher.append(planck_sd)

		rat1   = nhfk[i]/nhi[i]  ## For Fukui
		rat2   = nh[i]/nhi[i]    ## for Planck

		rat_fk = rat_fk + rat1 ## to cal. Mean of ratio
		rat_pl = rat_pl + rat2 ## to cal. Mean of ratio

		rat1   = round(rat1, 2)
		rat2   = round(rat2, 2)		

		wnm    = info['wnm'][i]
		cnm    = info['cnm'][i]

		print("{}  {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}"
			.format(i, src[i], round((nhfk[i]), 4), round((fukui_sd), 4), round((nh[i]), 4), round((planck_sd), 4), \
				nhi[i], nhier[i], wnm, cnm, rat1, rat2, oh[i], thin[i] ) )

	return nh, nher

## Get tau353 values and err_tau353 values #
 #
 # params str map_file File of maps
 # params dict info   Information of sources
 # params dict lownhi Information of 23 low-nhi sources
 #
 # return void
 #
 # version 11/2016 
 # Author Van Hiep
 ##	
def get_gas_column_density(map_file, info, lownhi):
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

	## Mask, to find Planck Conversion Factor (Dust opacity and Its error) ##
	msk    = hp.read_map(os.getenv("HOME")+'/hdata/dust/planck_mask.fits', field = 0, h=False)

	# tau353 map, err_tau353 map and resolution #
	tau_map  = hp.read_map(map_file, field = 0)
	err_map  = hp.read_map(map_file, field = 1)

	nh, nher   = nh_from_tau353(tau_map, err_map, msk, info)
	lnh, lnher = nh_from_tau353(tau_map, err_map, msk, lownhi)

	plt.errorbar(nhi, nh,xerr=nhier, yerr=nher, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data without CO')
	plt.errorbar(hi, lnh,xerr=hier, yerr=lnher, color='b', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data with low $N_{HI},  N_{HI}<3.0\cdot10^{20} cm^{-2}$')
	plt.plot([0,40],[0,40], 'k--', label='$N_{H} = N_{HI}$')
	plt.title('$N_{H}$ vs $N_{HI}$ along 26 lines-of-sight without CO & 23 LOS with low $N_{HI}$', fontsize=30)
	plt.ylabel('$N_{H}[10^{20}$ cm$^{-2}]$', fontsize=35)
	plt.xlabel('$N_{HI} [10^{20}$ cm$^{-2}]$', fontsize=35)
	plt.xlim(-1.0, 40.0)
	# plt.ylim(0, 3)
	plt.grid(True)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=18)

	plt.text(25., 2., '(Available sources with the presence of OH are shown)', color='k', fontsize=17)
	plt.text(25., 8., r'$N_{H} = \frac{\tau_{353}}{\sigma_{353}}$', color='k', fontsize=17)
	plt.text(25., 5., r'$\tau_{353}$ from Planck data R1.2', color='k', fontsize=17)

	plt.legend(loc='upper left', fontsize=18)
	# plt.savefig("test.png",bbox_inches='tight')
	for i in range(len(src)):
		if (oh[i] > 0) :
			plt.annotate('('+str(src[i])+')', xy=(nhi[i], nh[i]), xycoords='data',
	               xytext=(-50.,30.), textcoords='offset points',
	               arrowprops=dict(arrowstyle="->"),fontsize=15,
	               )
	plt.show()

	## NH vs Thin, LAB data => Thin ##
	plt.errorbar(thin, nh,xerr=thinr, yerr=nher, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data without CO')
	plt.errorbar(lthin, lnh,xerr=lthinr, yerr=lnher, color='b', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data with low $N_{HI},  N_{HI}<3.0\cdot10^{20} cm^{-2}$')
	plt.plot([0,40],[0,40], 'k--', label='$N_{H} = N^{thin}_{HI}$')
	plt.title('$N_{H}$ vs $N^{thin}_{HI}$ along 26 lines-of-sight without CO & 23 LOS with low $N_{HI}$', fontsize=30)
	plt.ylabel('$N_{H}[10^{20}$ cm$^{-2}]$', fontsize=35)
	plt.xlabel('$N^{thin}_{HI} [10^{20}$ cm$^{-2}]$', fontsize=35)
	plt.xlim(-1.0, 30.0)
	# plt.ylim(0, 3)
	plt.grid(True)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=18)

	plt.text(15., 2., '(Available sources with the presence of OH are shown)', color='k', fontsize=17)
	plt.text(15., 8., r'$N_{H} = \frac{\tau_{353}}{\sigma_{353}}$', color='k', fontsize=17)
	plt.text(15., 5., r'$\tau_{353}$ from Planck data R1.2', color='k', fontsize=17)

	plt.legend(loc='upper left', fontsize=18)
	# plt.savefig("test.png",bbox_inches='tight')
	for i in range(len(src)):
		if (oh[i] > 0) :
			plt.annotate('('+str(src[i])+')', xy=(thin[i], nh[i]), xycoords='data',
	               xytext=(-50.,30.), textcoords='offset points',
	               arrowprops=dict(arrowstyle="->"),fontsize=15,
	               )
	plt.show()


##================= MAIN ========================##
## Filename of the map
pth      = os.getenv("HOME")+'/hdata/dust/'
map_file = pth + 'HFI_CompMap_ThermalDustModel_2048_R1.20.fits'

# Info of 26 sources with no CO - l/b/name && 23 src low NHI #
info   = read_info_no_co('../../co12/result/26src_no_co_with_sponge.dat')
lownhi = read_lownhi_23src(fname = '../../hi/result/lownhi_thin_cnm_wnm.txt')

## cal N(H)
get_gas_column_density(map_file, info, lownhi)