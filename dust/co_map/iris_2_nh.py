import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import matplotlib        as mpl
import healpy            as hp
import pylab             as pl
import operator

from numpy    import array
from restore  import restore
from plotting import cplot

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
	beam   = 4.3            # Beam = 4.3'
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
	ir_map = hp.read_map(map_file, field = 0)
	nside  = hp.get_nside(ir_map)
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
		hp.cartview(ir_map, title=info['src'][i], coord='G', unit='',
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


## calculate N(H) from IRIS i100
 #
 # params str map_file File of maps
 # params dict info Information of sources
 #
 # return void
 # 
 # Author Van Hiep ##	
def cal_nh_from_i100(map_file, info):
	src = info['src']
	nhi = info['nhi']
	oh  = info['oh']

	# Define constants #
	deg2rad   = np.pi/180.

	# Define the width of area #
	beam   = 4.3            # Beam = 4.3'
	dbeam  = beam/120.0     # Beam = 30' -> dbeam = beam/60/2 in degree
	offset = dbeam          # degree

	ir_map = hp.read_map(map_file, field = 0)
	nside  = hp.get_nside(ir_map)
	res    = hp.nside2resol(nside, arcmin=False)
	dd     = res/deg2rad/10.0

	# OK - Go #
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

		if (ir_map[pix] > -0.000001) : # Some pixels not defined
			ci[i].append(ir_map[pix])

		for x in pl.frange(l-offset, l+offset, dd):
			for y in pl.frange(b-offset, b+offset, dd):
				cosb = np.cos(b*deg2rad)
				cosy = np.cos(y*deg2rad)

				if ( ((x-l)**2 + (y-b)**2) <= offset**2 ):
					# hp.projtext(x, y, '.', lonlat=True, coord='G')
					theta = (90.0 - y)*deg2rad
					phi   = x*deg2rad
					pix   = hp.ang2pix(nside, theta, phi, nest=False)

					if (ir_map[pix] > -0.000001) :
						ci[i].append(ir_map[pix])

		# plt.show()
		# continue

		vci = list(set(ci[i]))
		cnt = len(vci)

		# Calculate mean values of tau353 #
		# ebv.append(sum(ci[i])/float(cnt))
		val = sum(vci)/float(cnt)

		# Calculate the N(HI) from Fukui factor #
		# nhi_i = fukui_cf*tau353[i]
		# nhi.append(nhi_i)
	   
		# Calculate the NH from E(B-V) #
		# nh_i = tau353[i]/planck_cf
		# nh.append(nh_i)
		n_h = 1.13*val # 1e20; (NH = (1.0 +/- 0.3) e20*I100)

		# Uncertainties of mean values of tau353 #
		# sd1_tau353 = 0.
		# sd2_tau353 = 0.
		# for j in range(cnt_err):
		# 	sd1_tau353 = sd1_tau353 + (tau[i][j]-tau353[i])**2
		# 	sd2_tau353 = sd2_tau353 + (t_err[i][j])**2

		# sd2_tau353 = (sd2_tau353**0.5)/cnt_err # Uncertainty of a Sum

		# Uncertainties of mean values of N(HI) and N(H) #
		# fukui_sd  = (sd2_tau353*fukui_cf)*1e6

		# planck_sd = (sd2_tau353/tau353[i])**2 + (pl_fact_err/planck_cf)**2
		# planck_sd = (nh_i*planck_sd**0.5)*1e6

		# print("{}  {}\t{:08.4f}  {:08.4f}   {}   {}   {}   {}   {}   {}   {}   {}"
		# 	.format(i, src[i],l,b, info['nhi'][i], info['nhi_er'][i], info['nh2'][i],info['nh2_er'][i], info['nh'][i],info['nh_er'][i], nh_i*1e6, planck_sd   ))

		print src[i], val, nhi[i], n_h
		nh.append(n_h)

	nhi = np.asarray(nhi)
	nh  = np.asarray(nh)
	# x   = np.log10(nhi)
	# y   = nh/nhi 
	# plt.plot(x,y, 'r.', label='Ratio $f = N_{H}$/$N_{HI}$')
	# plt.title('Correlation between $f = N_{H}$/$N_{HI}$ and HI column density $N_{HI}$ \nalong 26 lines-of-sight without the presence of CO line', fontsize=30)
	# plt.ylabel('$Ratio f = N_{H}$/$N_{HI}$', fontsize=35)
	# plt.xlabel('log$_{10}$($N_{HI}/10^{20}$ cm$^{-2}$)', fontsize=35)
	# plt.xlim(0, 1.6)
	# plt.ylim(0, 3)
	# plt.grid(True)
	# plt.tick_params(axis='x', labelsize=18)
	# plt.tick_params(axis='y', labelsize=18)

	# plt.text(0.2, 0.31, '(Available sources with the presence of OH are shown)', color='k', fontsize=17)
	# # plt.text(0.2, 0.4, '$f = ['+str(m)+'\pm'+str(ea) +']\cdot log_{10}(N_{HI}/10^{20})+['+str(b)+'\pm'+str(eb)+']$', fontsize=20 )

	# plt.legend(loc='upper left', fontsize=18)
	# # plt.savefig("test.png",bbox_inches='tight')
	# for i in range(len(src)):
	# 	if (oh[i] > 0) :
	# 		plt.annotate('('+str(src[i])+')', xy=(x[i], y[i]), xycoords='data',
	#                xytext=(-50.,30.), textcoords='offset points',
	#                arrowprops=dict(arrowstyle="->"),fontsize=18,
	#                )
	# plt.show()

	plt.plot(nhi,nh, 'rd', label='data', ms=10)
	plt.plot([0,30],[0,30], 'k--', label='$N_{H} = N_{HI}$')
	plt.title('Correlation between $N_{H}$ and $N_{HI}$ \nalong 26 lines-of-sight without the presence of CO line', fontsize=30)
	plt.ylabel('$N_{H}[10^{20}$ cm$^{-2}]$', fontsize=35)
	plt.xlabel('$N_{HI} [10^{20}$ cm$^{-2}]$', fontsize=35)
	# plt.xlim(0, 1.6)
	# plt.ylim(0, 3)
	plt.grid(True)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=18)

	plt.text(15., 2., '(Available sources with the presence of OH are shown)', color='k', fontsize=17)

	plt.legend(loc='upper left', fontsize=18)
	# plt.savefig("test.png",bbox_inches='tight')
	for i in range(len(src)):
		if (oh[i] > 0) :
			plt.annotate('('+str(src[i])+')', xy=(nhi[i], nh[i]), xycoords='data',
	               xytext=(-50.,30.), textcoords='offset points',
	               arrowprops=dict(arrowstyle="->"),fontsize=18,
	               )
	plt.show()

#================= MAIN ========================#
pth      = os.getenv("HOME")+'/hdata/iris/'
map_file = pth + 'IRIS_4_2048.fits'

## Infor of 26 src without CO ##
info     = read_info_no_co('../../co12/result/26src_no_co_with_sponge.dat')

cal_nh_from_i100(map_file, info)
# plot_patches(map_file, info)