import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import healpy            as hp
import pylab             as pl
import operator

from numpy               import array
from restore             import restore
from plotting            import cplot
from masks               import masks

## Read info of 78 sources #
 # l,b, nhi, and nhi_error
 #
 # params string fname Filename
 # return dict info
 # 
 # version 04/10/2016
 # Author Van Hiep ##
def read_info_78src(fname = ''):
	cols = ['indx','src','l','b','nhi','nhi_er']
	fmt  = ['i',   's',  'f','f', 'f',  'f']
	data = restore(fname, 2, cols, fmt)
	return data.read()	

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
	fukui_cf  = 2.10 #2.10e26
	planck_cf = 0.84 #8.40e-27, 0.84e-26

	pl_fact_err = 0.3 #3.0e-27
	fk_fact_err = 0.0 #unknown

	# Define the width of area #
	beam   = 3.5            # Beam = 3.5'
	dbeam  = beam/120.0     # Beam = 3.5' -> dbeam = beam/60/2 in degree
	offset = dbeam          # degree

	# tau353 map, err_tau353 map and resolution #
	tau_map  = hp.read_map(map_file, field = 0)
	err_map  = hp.read_map(map_file, field = 1)
	nside    = hp.get_nside(tau_map)
	res      = hp.nside2resol(nside, arcmin=False)
	dd       = res/deg2rad/10.0

	# OK - Go #
	tau353 = []
	nh     = []
	nhi    = []

	tau    = {}
	t_err  = {}

	rat_fk = 0.
	rat_pl = 0.

	for i in range(0,len(src)):

		if( i != 15):
			continue

		# Find the values of Tau353 and Err_tau353 in small area #
		tau[i]   = []
		t_err[i] = []

		l  = info['l'][i]
		b  = info['b'][i]

		# Plot cartview a/o mollview #
		ll = l
		if (l>180):
			ll = ll-360.

		hp.cartview(tau_map, title=info['src'][i]+'('+str(info['l'][i])+','+str(info['b'][i])+') - '+map_file, coord='G', unit='',
				norm='hist', xsize=800, lonra=[ll-offset-0.1*offset,ll+offset+0.1*offset], latra=[b-offset-0.1*offset,b+offset+0.1*offset],
				return_projected_map=True)

		# hp.mollview(tau_map, title=info['src'][i]+'('+str(info['l'][i])+','+str(info['b'][i])+') - '+map_file,
		# 	coord='G', unit='', rot=[l,b,0], norm='hist', min=1e-7,max=1e-3, xsize=800)

		# Cal. #
		theta = (90.0-b)*deg2rad
		phi   = l*deg2rad
		pix   = hp.ang2pix(nside, theta, phi, nest=False)

		if (tau_map[pix] > -1.0e30) : # Some pixels not defined
			tau[i].append(tau_map[pix])

		if (err_map[pix] >= 6.9e-11 and err_map[pix] <= 0.00081): # Checked Invalid error
			t_err[i].append(err_map[pix])

		for x in pl.frange(l-offset, l+offset, dd):
			for y in pl.frange(b-offset, b+offset, dd):
				cosb = np.cos(b*deg2rad)
				cosy = np.cos(y*deg2rad)
				if ( (((x-l)**2 + (y-b)**2) <= offset**2) ):
					hp.projtext(x, y, '.', lonlat=True, coord='G')

		plt.show()


## Get tau353 values and err_tau353 values from 16 sources with CO, to compare N(H) derived from CO #
 # Than calculate N(H) from Dust
 #
 # params str map_file File of maps
 # params dict info Information of sources
 #
 # return void
 # 
 # Author Van Hiep ##	
def cal_nh_from_dust(map_file, src_num, info):
	## Classes ##
	msks = masks()

	# Define constants #
	deg2rad     = np.pi/180.
	fukui_cf    = 2.10e26
	fk_fact_err = 0.0 #unknown

	# Define the width of area #
	beam   = 3.5            # Beam = 3.5'
	dbeam  = beam/120.0     # Beam = 3.5' -> dbeam = beam/60/2 in degree
	offset = dbeam          # degree

	## g35+g45+g56 Mask, to find Planck Conversion Factor (Dust opacity and Its error) ##
	msk    = hp.read_map('../gas_col_density/data/planck_mask.fits', field = 0, h=False)

	# tau353 map, err_tau353 map and resolution #
	tau_map  = hp.read_map(map_file, field = 0)
	err_map  = hp.read_map(map_file, field = 1)
	nside    = hp.get_nside(tau_map)
	res      = hp.nside2resol(nside, arcmin=False)
	dd       = res/deg2rad/10.0

	# OK - Go #
	tau353 = []
	nh     = []
	nhi    = []

	tau    = {}
	t_err  = {}

	rat_fk = 0.
	rat_pl = 0.

	for i in range(0,src_num):
		# Find the values of Tau353 and Err_tau353 in small area #
		tau[i]   = []
		t_err[i] = []

		l = info['l'][i]
		b = info['b'][i]

		# # Plot cartview a/o mollview #
		# ll = l
		# if (l>180):
		# 	ll = ll-360.

		# hp.cartview(tau_map, title=info['src'][i]+'('+str(info['l'][i])+','+str(info['b'][i])+') - '+map_file, coord='G', unit='',
		# 		norm='hist', xsize=800, lonra=[ll-offset-0.1*offset,ll+offset+0.1*offset], latra=[b-offset-0.1*offset,b+offset+0.1*offset],
		# 		return_projected_map=True)
		# hp.cartview(err_map, title='Err', coord='G', unit='',
		# 		norm='hist', xsize=800, lonra=[ll-offset-0.1*offset,ll+offset+0.1*offset], latra=[b-offset-0.1*offset,b+offset+0.1*offset],
		# 		return_projected_map=True)
		# # End Plot cartview a/o mollview ##

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

		# plt.show()
		# continue
		# Make tau353 and err_tau353 correspondent #
		# print info['src'][i], len(tau[i]), len(t_err[i])

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
		nhi_i = fukui_cf*tau353[i]
		nhi.append(nhi_i)
	   
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

		nhi[i]    = nhi[i]*1.0e-20
		fukui_sd  = fukui_sd*1.0e-20
		nh[i]     = nh[i]*1.0e-20
		planck_sd = planck_sd*1.0e-20

		nhi_hi = info['nhi'][i] ## From Carl
		rat1   = nhi[i]/nhi_hi  ## For Fukui
		rat2   = nh[i]/nhi_hi   ## for Planck

		rat_fk = rat_fk + rat1 ## to cal. Mean of ratio
		rat_pl = rat_pl + rat2 ## to cal. Mean of ratio

		rat1   = round(rat1, 2)
		rat2   = round(rat2, 2)

		print("{}  {}\t{:08.4f}  {:08.4f}   {}\t{}\t{}\t\t{}"
			.format(i, info['src'][i],l,b, info['nhi'][i], info['nhi_er'][i], round(nh[i],4), round(planck_sd,4)    ))

##================= MAIN ========================##
## Define constants ##
map_file = '../gas_col_density/data/HFI_CompMap_ThermalDustModel_2048_R1.20.fits'

## Infor of 78 src ##
info = read_info_78src(fname = '../hi/rearrange/nhi_lb_78src.txt')
num_of_src = len(info['src'])

cal_nh_from_dust(map_file, num_of_src, info)
# plot_patches(map_file, info)