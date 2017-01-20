import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import matplotlib        as mpl
import healpy            as hp
import pylab             as pl
import operator
import copy

from numpy               import array
from restore             import restore
from plotting            import cplot
from   mpfit             import mpfit

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
	d1 = d1**2

	d2 = nh_er/nh
	d2 = d2**2

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
	ret ={}

	sxy = 0.
	sx  = 0.
	sy  = 0.
	sx2 = 0.
	n   = len(x)
	for i in range(0,n) :
		sxy = sxy + x[i]*y[i]
		sx  = sx  + x[i]
		sy  = sy  + y[i]
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

	ret['a']  = a
	ret['b']  = b
	ret['ea'] = er_a
	ret['eb'] = er_b

	return ret

## Read info of 30 SPONGE sources #
 # l,b, nhi, and nhi_error
 #
 # params string fname Filename
 # return dict info
 # 
 # version 12/2016
 # Author Van Hiep ##
def read_info_ms_78src(fname = '../rearrange/nhi_lb_thin_78src.txt'):
	cols = ['idx','src','l', 'b', 'nhi','nhi_er','thin','thin_er', 'cnm','cnm_er','wnm','wnm_er']
	fmt  = ['i',  's',  'f', 'f', 'f',   'f',     'f',    'f'    , 'f',   'f',     'f',  'f'    ]
	data = restore(fname, 2, cols, fmt)
	dat  = data.read()
	return dat

## Read NHI from 94src #
 #
 # params string fname Filename
 # return dict info of N(HI)
 # 
 # version 1/2017
 # Author Van Hiep ##
def read_nhi_94src(fname = '../result/nhi_thin_cnm_wnm_94src_sponge_prior.txt'):
	cols = ['indx', 'src', 'l', 'b', 'nhi', 'nhi_er', 'thin', 'thin_er', 'cnm', 'cnm_er', 'wnm', 'wnm_er']
	fmt  = ['i',     's',   'f','f',  'f',   'f',     'f',     'f',        'f',   'f',      'f',   'f'   ]
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
	whi_map = hp.read_map(map_file, field = 2)
	nside  = hp.get_nside(whi_map)
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
		hp.cartview(whi_map, title=info['src'][i], coord='G', unit='',
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


## Get NH from E(B-V) #
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
	src   = info['src']  ## 78 src
	nhi   = info['nhi']
	nhier = info['nhi_er']

	# Define constants #
	deg2rad  = np.pi/180.
	dv       = 1.030571969
	cst      = 1.8224e18

	# Define the width of area #
	beam   = 30.0           # Beam = 30'
	dbeam  = beam/120.0     # Beam = 30' -> dbeam = beam/60/2 in degree
	offset = dbeam          # degree

	# TTYPE1  = 'WHI  '           / https://arxiv.org/PS_cache/arxiv/pdf/0706/0706.1703v2.pdf, page 1                
	# TTYPE2  = 'SOMETHING '      / Error on optical depth
	st_map  = hp.read_map(map_file, field = 0) ## Map of sum(T)/Number_of_channels
	ch_map  = hp.read_map(map_file, field = 1) ## Map of Number of channels
	whi_map = cst * st_map*ch_map*dv/1e20      ## NHI map in [1e20]
	nside   = hp.get_nside(whi_map)
	res     = hp.nside2resol(nside, arcmin=False)
	dd      = res/deg2rad/10.0

	# OK - Go 1: 26 src without CO #
	whi    = []
	lhi    = []
	lhier  = []

	wi     = {}
	wi_err = {}
	for i in range(len(src)):
		# Find the values of WHI and Err_WHI in small area #
		wi[i] = []

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

		if (whi_map[pix] > -0.000001) : # Some pixels not defined
			wi[i].append(whi_map[pix])

		for x in pl.frange(l-offset, l+offset, dd):
			for y in pl.frange(b-offset, b+offset, dd):
				cosb = np.cos(b*deg2rad)
				cosy = np.cos(y*deg2rad)

				if ( ((x-l)**2 + (y-b)**2) <= offset**2 ):
					# hp.projtext(x, y, '.', lonlat=True, coord='G')
					theta = (90.0 - y)*deg2rad
					phi   = x*deg2rad
					pix   = hp.ang2pix(nside, theta, phi, nest=False)

					if (whi_map[pix] > -0.000001) :
						wi[i].append(whi_map[pix])

		# plt.show()
		# continue

		vwi = list(set(wi[i]))
		cnt = len(vwi)

		# Calculate mean values of tau353 #
		val = sum(vwi)/float(cnt)
		err = cal_uncertainty_in_mean(vwi)
	   
		# Calculate the NHI #
		n_hi = val
		err  = err

		print("{}  {}\t{:08.4f}  {:08.4f}   {}   {}   {}   {}"
			.format(i, src[i],l,b, info['nhi'][i], info['nhi_er'][i],  n_hi/1e20, err/1e20  ))

		# print src[i], val, nhi[i], n_h*100.
		lhi.append(n_hi)
		lhier.append(err)

	# Fit and Plot #
	params = linear_fit(lhi,nhi)
	a      = params['a']
	b      = params['b']
	ea     = params['ea']
	eb     = params['eb']

	plt.errorbar(lhi,nhi,xerr=lhier, yerr=nhier, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data')
	plt.plot([0,200],[0,200], 'k--', label='$N_{H} = N_{HI}$')
	plt.plot(lhi, a*np.array(lhi) + b, 'r-', linewidth=4, label='Best linear fit')

	a  = round(a, 2)
	b  = round(b, 2)
	ea = round(ea, 2)
	eb = round(eb, 2)

	plt.title('$N_{HI}$ and $N^{thin}_{HI}$ along 94 lines-of-sight', fontsize=30)
	plt.ylabel('$N_{HI}[10^{20}$ cm$^{-2}]$, from 21-SPONGE and MS data', fontsize=35)
	plt.xlabel('$N^{thin}_{HI} [10^{20}$ cm$^{-2}]$, LAB data', fontsize=35)
	plt.xlim(-1.0, 60.0)
	# plt.ylim(0, 3)
	plt.grid(True)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=18)

	plt.text(0.0, 50.0, '$f = ['+str(a)+'\pm'+str(ea) +']\cdot log_{10}(N^*_{HI}/10^{20}) + ['+str(b)+'\pm'+str(eb)+']$', color='blue', fontsize=20)

	plt.legend(loc='upper left', fontsize=18)
	# plt.savefig("test.png",bbox_inches='tight')
	plt.show()

	## 94 MS+SPONGE Sources, 21SPONGE* prior, R vs log10(N*HI)
	fact   = []
	lognhi = []
	for i in range(len(lhi)):
		temp = round(nhi[i]/lhi[i], 3)
		lognhi.append(np.log10(lhi[i]))
		fact.append(temp)

	# Error bar for x-axis and y-axis
	thin   = np.asarray(lhi)
	hi     = np.asarray(nhi)
	xdata  = np.asarray(lognhi)
	ydata  = np.asarray(fact)
	hier   = np.array(nhier)
	thiner = np.array(lhier)
	xerr   = thiner/thin/np.log(10.0)
	yerr   = uncertainty_of_factors(ydata, thin, thiner, hi, hier)

	# Fit and Plot #
	params = linear_fit(lognhi,fact)
	a      = params['a']
	b      = params['b']
	ea     = params['ea']
	eb     = params['eb']

	# plt.plot(lognhi, fact, 'b^', label='Ratio $f = N_{HI}$/$N^*_{HI}$', markersize=10)
	plt.errorbar(xdata,ydata,xerr=xerr, yerr=yerr, color='b', marker='o', ls='None', markersize=4, markeredgecolor='b', markeredgewidth=1, label='Ratio $f = N_{HI}$/$N^*_{HI}$')
	plt.plot(lognhi, a*np.array(lognhi) + b, 'r-', linewidth=4, label='Best linear fit')

	a  = round(a, 2)
	b  = round(b, 2)
	ea = round(ea, 2)
	eb = round(eb, 2)

	plt.ylabel('$Ratio f = N_{HI}$/$N^*_{HI}$', fontsize=35)
	plt.xlabel('$log_{10}(N^*_{HI}/10^{20} cm^{-2}$)', fontsize=35)
	plt.title('Correlation between Total HI column densities $N_{HI}$ and \n HI optically thin column densities $N^*_{HI}$ along 94 (Millennium Survey & 21-SPONGE) lines-of-sight', fontsize=30)
	plt.grid(True)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=15)

	plt.text(0.0, 3.0, '$f = ['+str(a)+'\pm'+str(ea) +']\cdot log_{10}(N^*_{HI}/10^{20}) + ['+str(b)+'\pm'+str(eb)+']$', color='blue', fontsize=20)
	plt.text(0.0, 3.2, r'$f = [0.32\pm0.06]\cdot log_{10}(N^*_{HI}/10^{20}) + [0.81\pm0.05]$, Lee et al.', color='blue', fontsize=20)

	plt.legend(loc='upper left', fontsize=18)
	# plt.savefig("test.png",bbox_inches='tight')
	plt.show()

	########### MPFIT ############
	xdata = xdata
	ydata = ydata

	# Error bar for x-axis and y-axis
	xerr = xerr
	yerr = yerr

	## Fit ##
	lguess  = [ 0.33, 0.85]

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

	# plt.plot(xdata,ydata, 'ok', ls='None', marker='.', lw=1, label='Factor $f = N_{HI}$/$N^*_{HI}$')
	plt.errorbar(xdata,ydata,xerr=xerr, yerr=yerr, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='Ratio $f = N_{HI}$/$N^*_{HI}$')
	plt.plot(xfit, mu, '-b', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='MPFIT linear fit')
	plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)
	# plt.plot(xfit, alpha*xfit[:, None] + beta, c='blue', alpha=0.01)
	# plt.plot(xfit, alpha_stats['mean']*xfit + beta_stats['mean'], linewidth=2, c='red')

	plt.title('Correlation between Total HI column densities $N_{HI}$ and \n HI optically thin column densities $N^*_{HI}$ along 94 (Millennium Survey & 21-SPONGE) lines-of-sight', fontsize=30)
	plt.ylabel('$Ratio f = N_{HI}$/$N^*_{HI}$', fontsize=35)
	plt.xlabel('$log_{10}(N^*_{HI}/10^{20} cm^{-2}$)', fontsize=35)
	plt.xlim(0.0, 2.0)
	plt.ylim(-1.0, 6.0)
	plt.grid(True)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=15)

	plt.text(0.0, 3.0, '$f = ['+str(m)+'\pm'+str(ea) +']\cdot log_{10}(N^*_{HI}/10^{20}) + ['+str(b)+'\pm'+str(eb)+']$', color='blue', fontsize=20)
	plt.text(0.0, 3.2, r'$f = [0.32\pm0.06]\cdot log_{10}(N^*_{HI}/10^{20}) + [0.81\pm0.05]$, Lee et al.', color='blue', fontsize=20)
	plt.legend(loc='lower right', fontsize=18)
	# plt.savefig("test.png",bbox_inches='tight')
	# for i in range(len(sc)):
	# 	plt.annotate('('+str(sc[i])+')', xy=(xdata[i], ydata[i]), xycoords='data',
	#             xytext=(-50.,30.), textcoords='offset points',
	#             arrowprops=dict(arrowstyle="->"),fontsize=18,
	#             )
	plt.show()
	########### END - MPFIT ############

#================= MAIN ========================#
pth      = os.getenv("HOME")+'/hdata/dust/LAB_Healpix/'
map_file = pth + 'LAB_fullvel.fits'

## Infor of 94 src without CO && 23 Low NHI sources ##
info     = read_nhi_94src(fname = '../../hi/result/nhi_thin_cnm_wnm_94src_sponge_prior.txt')
info     = read_info_ms_78src(fname = '../result/nhi_lb_thin_cnm_wnm_78src.txt')
lowhi    = read_lownhi_23src(fname = '../../hi/result/lownhi_thin_cnm_wnm.txt')

cal_nh_from_dust(map_file, info, lowhi)
# plot_patches(map_file, info)