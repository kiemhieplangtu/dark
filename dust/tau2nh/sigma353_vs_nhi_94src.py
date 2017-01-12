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

## Read NHI from 94src #
 #
 # params string fname Filename
 # return dict info of N(HI)
 # 
 # version 1/2017
 # Author Van Hiep ##
def read_nhi_94src(fname = '../../hi/result/nhi_thin_cnm_wnm_94src_sponge_prior.txt'):
	cols = ['indx', 'src', 'l', 'b', 'nhi', 'nhi_er', 'thin', 'thin_er', 'cnm', 'cnm_er', 'wnm', 'wnm_er']
	fmt  = ['i',     's',   'f','f',  'f',   'f',     'f',     'f',        'f',   'f',      'f',   'f'   ]
	data = restore(fname, 2, cols, fmt)
	return data.read()

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

## Plot the correlation between Sigma353 and NHI #
 #
 # params 
 # return Void
 # 
 # version 9/2016
 # Author Van Hiep ##
def plt_sigma353_vs_nhi():
	## Classes ##
	msks = masks()
	## Mask, to find Planck Conversion Factor (Dust opacity and Its error) ##
	msk  = hp.read_map(os.getenv("HOME")+'/hdata/dust/planck_mask.fits', field = 0, h=False)

	## Read info of 94 src
	info   = read_nhi_94src(fname = '../../hi/result/nhi_thin_cnm_wnm_94src_sponge_prior.txt')
	src    = info['src']
	nhi    = info['nhi']
	nhi_er = info['nhi_er']

	## Read info of 94 src
	info26 = read_info_no_co(fname = '../../co12/result/26src_no_co_with_sponge.dat')
	sc26   = info26['src']

	# Define constants #
	deg2rad  = np.pi/180.
	pth      = os.getenv("HOME")+'/hdata/dust/'
	map_file = pth + 'HFI_CompMap_ThermalDustModel_2048_R1.20.fits'

	# Define the width of area #
	beam   = 5.             # Beam = 5'
	dbeam  = beam/120.0     # Beam = 5' -> dbeam = beam/60/2 in degree
	offset = dbeam          # degree

	# tau353 map, err_tau353 map and resolution #
	tau_map  = hp.read_map(map_file, field = 0)
	err_map  = hp.read_map(map_file, field = 1)
	nside    = hp.get_nside(tau_map)
	res      = hp.nside2resol(nside, arcmin=False)
	dd       = res/deg2rad/10.0

	# OK - Go #
	tau353  = []
	tau_er  = []
	nh      = []
	nher    = []
	sg353   = []
	sg353er = []

	tau    = {}
	t_err  = {}
	for i in range(94):
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
		tau_val = sum(tau[i])/float(npix)
		tau353.append(tau_val)

		# Calculate the NH from Planck factor #
		nh_i = tau353[i]/planck_cf
		nh.append(nh_i)

		# Uncertainties of mean values of tau353 #
		sd_tau353 = 0.
		for j in range(npix):
			sd_tau353 = sd_tau353 + (t_err[i][j])**2

		sd_tau353 = (sd_tau353**0.5)/npix # Uncertainty of a Sum
		tau_er.append(sd_tau353)

		# Uncertainties of mean values of N(HI) and N(H) #
		planck_sd = (sd_tau353/tau353[i])**2 + (pl_fact_err/planck_cf)**2
		planck_sd = (nh_i*planck_sd**0.5)

		nh[i]     = nh[i]*1.0e-20
		planck_sd = planck_sd*1.0e-20
		sig353    = tau_val/(nhi[i]*1.0e20)
		sig353er  = uncertainty_of_ratio(tau_val, nhi[i]*1.0e20, sd_tau353, nhi_er[i]*1.0e20)
		sg353.append(sig353/1.0e-27)
		sg353er.append(sig353er/1.0e-27)
		nher.append(planck_sd)

		# print("{}  {}\t{}\t{}\t{}\t{}"
		# 	.format(i, src[i], nhi[i], nhi_er[i], round((nh[i]), 4), round((planck_sd), 4) ) )

	# plt.plot(nhi,nh, 'rd', label='data no CO', ms=10)
	# plt.errorbar(nhi,nh,xerr=nhi_er, yerr=nher, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data')
	# plt.plot([0,130],[0,130], 'k--', label='$N_{H} = N_{HI}$')
	# plt.title('Correlation between $N_{H}$ and $N_{HI}$ \nalong 94 lines-of-sight', fontsize=30)
	# plt.ylabel('$N_{H}[10^{20}$ cm$^{-2}]$', fontsize=35)
	# plt.xlabel('$N_{HI} [10^{20}$ cm$^{-2}]$', fontsize=35)
	# plt.xlim(-10.0, 165.0)
	# # plt.ylim(0, 3)
	# plt.grid(True)
	# plt.tick_params(axis='x', labelsize=18)
	# plt.tick_params(axis='y', labelsize=18)

	# plt.text(100., 30., r'$N_{H} = \tau_{353}/\sigma_{353}$', color='k', fontsize=17)
	# plt.text(100., 10., r'$\tau_{353}$ from Planck data R1.2', color='k', fontsize=17)

	# plt.legend(loc='upper left', fontsize=18)
	# # plt.savefig("test.png",bbox_inches='tight')
	# for i in range(len(src)):
	# 	# if (oh[i] > 0) :
	# 	plt.annotate('('+str(src[i])+')', xy=(nhi[i], nh[i]), xycoords='data',
 #               xytext=(-50.,30.), textcoords='offset points',
 #               arrowprops=dict(arrowstyle="->"),fontsize=18,
 #               )
	# plt.show()


	##===== color for 26src without CO =====##
	nhi26     = []
	nhi_er26  = []
	sg26      = []
	sg353er26 = []
	lsig      = []
	for i in range(len(nhi)):
		if(src[i] in sc26):
			nhi26.append(nhi[i])
			nhi_er26.append(nhi_er[i])
			sg26.append(sg353[i])
			sg353er26.append(sg353er[i])

		if(nhi[i]<3.0):
			lsig.append(sg353[i])

	lsig_mean = np.mean(np.array(lsig))
	## Plot Sigma353 vs NHI
	plt.xscale("log", nonposx='clip')
	plt.errorbar(nhi,sg353,xerr=nhi_er, yerr=sg353er, color='r', marker='o', ls='None', markersize=12, markeredgecolor='b', markeredgewidth=1, label='data')
	plt.errorbar(nhi26,sg26,xerr=nhi_er26, yerr=sg353er26, color='b', marker='o', ls='None', markersize=12, markeredgecolor='b', markeredgewidth=1, label='data, LOS without CO')
	plt.plot([0,1000],[6.2,6.2], 'k-.', label='From Planck Fig.20: \nAverage values found in the low $N_{HI}$ mask', lw=3)
	plt.plot([0,1000],[lsig_mean,lsig_mean], 'b--', label='From Hiep: \nAverage values found in the low $N_{HI}$ mask, ~' + str(round(lsig_mean,2)), lw=3 )
	plt.title('$\sigma_{353}$ vs $N_{HI}$', fontsize=35)
	plt.ylabel('$\sigma_{353}[10^{-27}$ cm$^{2}H^{-1}]$', fontsize=35)
	plt.xlabel('$N_{HI} [10^{20}$ cm$^{-2}]$', fontsize=35)
	# plt.xlim(0, 1.6)
	plt.ylim(0., 20.)
	plt.grid(True)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=18)

	# plt.text(15., 2., '(Available sources with the presence of OH are shown)', color='k', fontsize=17)

	plt.legend(loc='upper left', fontsize=18)
	# plt.savefig("test.png",bbox_inches='tight')
	# for i in range(len(src)):
	# 	if (oh[i] > 0) :
	# 		plt.annotate('('+str(src[i])+')', xy=(nhi[i], nh_pl[i]), xycoords='data',
	#                xytext=(-50.,30.), textcoords='offset points',
	#                arrowprops=dict(arrowstyle="->"),fontsize=18,
	#                )

	plt.show()

	## Plot Tau353 vs NHI
	tau353 = np.array(tau353)*1e6
	plt.errorbar(nhi,tau353,xerr=nhi_er, yerr=tau_er, color='r', marker='o', ls='None', markersize=12, markeredgecolor='b', markeredgewidth=1, label='data')
	plt.title(r'$\tau_{353}$ vs $N_{HI}$', fontsize=35)
	plt.plot([0,3],[-0.02,1.89], 'k--', label='Linear fit from Planck paper', lw=3)
	plt.ylabel(r'$\tau_{353}[10^{-6}]$', fontsize=35)
	plt.xlabel('$N_{HI} [10^{20}$ cm$^{-2}]$', fontsize=35)
	# plt.ylim(0., 2.5)
	# plt.xlim(0., 3.0)
	plt.grid(True)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=18)

	plt.text(0.5, 0.2, r'$\tau_{353} = (6.3\pm0.1) \cdot 10^-{27} - 0.02 \cdot 10^{-6}$', color='k', fontsize=17)

	plt.legend(loc='upper left', fontsize=20)
	# plt.savefig("test.png",bbox_inches='tight')
	# for i in range(len(src)):
	# 	if (oh[i] > 0) :
	# 		plt.annotate('('+str(src[i])+')', xy=(nhi[i], nh_pl[i]), xycoords='data',
	#                xytext=(-50.,30.), textcoords='offset points',
	#                arrowprops=dict(arrowstyle="->"),fontsize=18,
	#                )

	plt.show()

## Calculate the uncertainties of ratio a/b #
 #
 # params float a
 # params float b
 # params float aer
 # params float ber
 #
 # return float ret Uncertainty of a/b
 # 
 # Author Van Hiep ##
def uncertainty_of_ratio(a, b, aer, ber):
	r  = np.abs(a/b)
	d1 = aer/a
	d1 = d1*d1

	d2 = ber/b
	d2 = d2*d2

	d  = r * np.sqrt(d1+d2)

	return d	

## Calculate the uncertainties of factors #
 #
 # params 1D-array factr y-axis
 # params 1D-array nhi N(HI)
 # params 1D-array nhi_er Uncertainties of N(HI)
 # params 1D-array nh N(H)
 # params 1D-array nh_er Uncertainties of N(H)
 #
 # return factor_uncertainties
 # 
 # Author Van Hiep ##
def uncertainty_of_factors(factr, nhi, nhi_er, nh, nh_er):
	d1 = nhi_er/nhi
	d1 = d1*d1

	d2 = nh_er/nh
	d2 = d2*d2

	d  = np.sqrt(d1+d2)*factr

	return d

## linear fit #
 #
 # params x list x-data
 # params y list y-data
 #
 # return fit parameters and errors
 # 
 # Author Van Hiep ##
def linear_fit(x,y):
	sxy = 0.
	sx  = 0.
	sy  = 0.
	sx2 = 0.
	n   = len(x)
	for i in range(0,n) :
		sxy = sxy + x[i]*y[i]
		sx  = sx + x[i]
		sy  = sy + y[i]
		sx2 = sx2 + x[i]**2

	denom = (n*sx2 - sx**2)
	a = (n*sxy - sx*sy)/denom
	b = (sx2*sy - sx*sxy)/denom

	t    = n*sx2 - sx**2
	er_a = np.sqrt(n/t) 
	er_b = np.sqrt(sx2/t) 

	chi2 = 0.
	for i in range(0,n) :
		chi2 = chi2 + (y[i]-a*x[i]-b)**2

	chi2 = np.sqrt(chi2/(n-2))
	er_a = chi2*er_a
	er_b = chi2*er_b

	return a,b,er_a,er_b

#================= MAIN ========================#
plt_sigma353_vs_nhi()