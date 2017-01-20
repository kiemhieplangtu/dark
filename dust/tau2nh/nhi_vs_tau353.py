import os, sys
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import healpy            as hp
import pylab             as pl
import pymc3             as pm
import operator
import matplotlib.cm     as cm
import copy

from restore             import restore
from masks               import masks
from mpfit               import mpfit

## Linear ##
 #
 # params 
 # params 
 #
 # return 
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
# def tb_exp(x,y,bg,tau,v0,wid,tex,cont,err=1.):
def myfunc(p, fjac=None, x=None, y=None, err=None):
	model  = p[0] * x + p[1]
	status = 0
	return [status, (y - model) / err]

## Read info of 26 no-CO sources #
 # l,b, nhi, and nhi_error
 #
 # params string fname Filename
 # return dict infocd 
 # 
 # version 11/2016
 # Author Van Hiep ##
def read_info_78src(fname = '../hi/rearrange/nhi_lb_78src.txt'):
	cols = ['idx','src','l','b','nhi','nhi_er']
	fmt  = ['i',  's',  'f','f', 'f',    'f']
	data = restore(fname, 2, cols, fmt)
	dat  = data.read()
	return dat

## Get tau353 values and err_tau353 values #
 #
 # params str map_file File of maps
 # params int src_num Number of sources
 # params dict info Information of sources
 #
 # return void
 # 
 # Author Van Hiep ##	
def plot_patches(map_file, src_num, info):
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

	for i in range(0,src_num):

		# Find the values of Tau353 and Err_tau353 in small area #
		tau[i]   = []
		t_err[i] = []

		l  = info['l'][i]
		b  = info['b'][i]

		sr = 15

		# Plot cartview a/o mollview #
		ll = l
		if (l>180):
			ll = ll-360.

		if (i == sr):
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
				if ( (((x-l)**2 + (y-b)**2) <= offset**2) and (i == sr) ):
					hp.projtext(x, y, '.', lonlat=True, coord='G')

	plt.show()

## Get tau353 values and err_tau353 values #
 #
 # params str map_file File of maps
 # params int src_num Number of sources
 # params dict info Information of sources
 #
 # return void
 #
 # version 11/2016 
 # Author Van Hiep
 ##	
def get_gas_column_density(map_file, src_num, info):
	## Classes ##
	msks = masks()

	# Define constants #
	deg2rad     = np.pi/180.
	fukui_cf    = 2.10e26
	fk_fact_err = 0.0 #unknown

	## Read info
	arcb_nhi = np.asarray(info['nhi'])
	arcb_er  = np.asarray(info['nhi_er'])

	# Define the width of area #
	beam   = 3.5            # Beam = 3.5'
	dbeam  = beam/120.0     # Beam = 3.5' -> dbeam = beam/60/2 in degree
	offset = dbeam          # degree

	## g35+g45+g56 Mask, to find Planck Conversion Factor (Dust opacity and Its error) ##
	msk    = hp.read_map('../gas_col_density/data/mask_g35g45g56.fits', field = 0, h=False)

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

		# print i, info['src'][i], nhi_hi, info['nhi_er'][i], tau353[i], sd2_tau353

	x_filt = []
	y_filt = []
	x_er   = []
	y_er   = []
	for i in range(len(arcb_nhi)):
		if(arcb_nhi[i] < 57.):
			x_filt.append(arcb_nhi[i]*1e20)
			y_filt.append(tau353[i])
			x_er.append(arcb_er[i]*1e20)
			y_er.append(sd2_tau353)

	xdata = np.asarray(x_filt)
	ydata = np.asarray(y_filt)
	a,b   = np.polyfit(xdata, ydata, 1)
	print a,b
	# plt.plot(np.asarray(arcb_nhi)*1e20,np.asarray(tau353), 'b.')
	# plt.plot(xdata, ydata, 'b.')
	# plt.xlim(0.,57.e20)
	# plt.ylim(0.,1.5e-4)
	# plt.grid()
	# plt.show()

	# sys.exit()

	# Fit and Plot #
	# Error bar for x-axis and y-axis
	xerr = np.asarray(x_er)
	yerr = np.asarray(y_er)

	## Fit  Tau, V0, Width ##
	lguess  = [ 1.5e-26, -1.4e-6]

	npar    =  len(lguess)
	guessp  = np.array(lguess, dtype='float64')
	plimd   = [[False,False]]*npar
	plims   = [[0.,0.]]*npar
	parbase = {'value': 0., 'fixed': 0, 'parname': '', 'limited': [0, 0], 'limits': [0., 0.]}
	pname   = ['slope','offset']
	pfix    = [False]*npar

	parinfo = []
	for i in range(len(guessp)):
		parinfo.append(copy.deepcopy(parbase))

	for i in range(len(guessp)):
		parinfo[i]['value']   = guessp[i]
		parinfo[i]['fixed']   = pfix[i]
		parinfo[i]['parname'] = pname[i]
		parinfo[i]['limited'] = plimd[i]

	##  1665 ###
	x    = xdata.astype(np.float64)
	y    = ydata.astype(np.float64)
	er   = yerr.astype(np.float64)

	fa = {'x':x, 'y':y, 'err':er}
	mp = mpfit(myfunc, guessp, parinfo=parinfo, functkw=fa, quiet=True)

	## ********* Results ********* ##
	print '********* Results *********'
	abp   = mp.params
	abper = mp.perror
	for i in range(len(parinfo)):
		print "%s = %f +/- %f" % (parinfo[i]['parname'],abp[i],abper[i])
	## Plot ##
	a     = np.array([ abp[0]-abper[0], abp[0]+abper[0] ])
	b     = np.array([ abp[1]-abper[1], abp[1]+abper[1] ])
	xfit  = np.linspace(xdata.min(), xdata.max(), 20)
	yfit  = a[:, None] * xfit + b[:, None]
	mu    = yfit.mean(0)
	sig   = 1.0*yfit.std(0)
	fit   = abp[0]*x+abp[1]

	m  = round(abp[0],2)
	b  = round(abp[1],2)
	ea = round(abper[0],2)
	eb = round(abper[1],2)

	plt.plot(xdata,ydata, 'ok', ls='None', marker='.', lw=1, label='Factor $f = N_{HI}$/$N^*_{HI}$')
	plt.errorbar(xdata,ydata,xerr=xerr, yerr=yerr, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='Ratio $f = N_{H}$/$N_{HI}$')
	plt.plot(xfit, mu, '-b', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='MPFIT linear fit')
	plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)
	# plt.plot(xfit, alpha*xfit[:, None] + beta, c='blue', alpha=0.01)
	# plt.plot(xfit, alpha_stats['mean']*xfit + beta_stats['mean'], linewidth=2, c='red')

	plt.title('Correlation between $f = N_{H}$/$N_{HI}$ and HI column density $N_{HI}$ \nalong 26 lines-of-sight without the presence of CO line', fontsize=30)
	plt.ylabel('$Ratio f = N_{H}$/$N_{HI}$', fontsize=35)
	plt.xlabel('log$_{10}$($N_{HI}/10^{20}$ cm$^{-2}$)', fontsize=35)
	plt.xlim(0.,57.e20)
	plt.ylim(0.,1.5e-4)
	plt.grid(True)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=18)

	# plt.text(0.2, 0.31, '(Available sources with the presence of OH are shown)', color='k', fontsize=17)
	# plt.text(0.2, 0.4, '$f = ['+str(m)+'\pm'+str(ea) +']\cdot log_{10}(N_{HI}/10^{20})+['+str(b)+'\pm'+str(eb)+']$', fontsize=20 )

	plt.legend(loc='upper left', fontsize=18)
	# plt.savefig("test.png",bbox_inches='tight')
	# for i in range(len(src)):
	# 	if (oh[i] > 0) :
	# 		plt.annotate('('+str(src[i])+')', xy=(xdata[i], ydata[i]), xycoords='data',
	#                xytext=(-50.,30.), textcoords='offset points',
	#                arrowprops=dict(arrowstyle="->"),fontsize=18,
	#                )
	plt.show()


#================= MAIN ========================#
# Define constants #
map_file = '../gas_col_density/data/HFI_CompMap_ThermalDustModel_2048_R1.20.fits'

# Info of 26 sources with no CO - l/b/name #
info       = read_info_78src('../hi/rearrange/nhi_lb_78src.txt')
num_of_src = len(info['src'])

get_gas_column_density(map_file, num_of_src, info)
# plot_patches(map_file, num_of_src, info)