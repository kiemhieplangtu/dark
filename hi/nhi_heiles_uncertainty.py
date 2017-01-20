import os, sys, shutil
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt
import pylab             as pl
import copy

from numpy   import array
from restore import restore
import operator

## Read N(HI)_Heiles for each component #
 # Then calculate the uncertainties for each component
 #
 # params string fname Filename
 #
 # return dict Gaussian-fit parameters of each component
 # including the uncertainties of each component
 # 
 # version 10/2016
 # Author Van Hiep ##
def read_fit_params(fname = 'result/component_fit_params.txt'):
	## Define constants ##
	const = 0.018224 #1.8224/100. -> N(HI) [1e20/cm2]
	fct   = 2.0*((np.log(2))**0.5)

	## Read fit params of HI components from Carl papers ##
	cols  = ['t_peak','err_t_peak','tau','err_tau','v0','err_v0','del_v','err_del_v', 'tspin', 'err_tspin', 'tkmax', 'nhi', 'cnm', 'frac', 'err_frac', 'src']
	fmt   = ['f',     'f',          'f',  'f',      'f', 'f',     'f',    'f',         'f',      'f',         'f',     'f',   'i',    'f',   'f',        's']
	dat   = restore(fname, 3, cols, fmt)
	inf   = dat.read()

	src       = inf['src']
	frac      = inf['frac']
	nhi       = inf['nhi']
	err_tau   = inf['err_tau']
	err_tspin = inf['err_tspin']

	del_v     = inf['del_v']
	err_del_v = inf['err_del_v']
	tau       = inf['tau']
	tspin     = inf['tspin']
	tb        = inf['t_peak']
	err_tb    = inf['err_t_peak']

	## Cal. uncertainties for each component ##
	ret = {}
	ret['nhi_uncertainty'] = []
	ret['cnm_uncertainty'] = []
	ret['wnm_uncertainty'] = []
	ret['src']             = src
	ret['nhi']             = nhi
	ret['cnm']             = []
	ret['wnm']             = []
	for i in range(len(src)):
		del_v[i]      = del_v[i]/fct
		err_del_v [i] = err_del_v[i]/fct
		if(frac[i] < 0.) : # uncertainties for CNM component
			# if((nhi[i] == 0.) and (err_tau[i] == 0.) and (err_tspin[i] == 0.)):
			# 	nhi_uncertainty = 0.
			# 	cnm_uncertainty = 0.
			# else:
			d1 = del_v[i]*tau[i]*err_tspin[i]
			d1 = d1**2

			d2 = del_v[i]*tspin[i]*err_tau[i]
			d2 = d2**2

			d3 = tspin[i]*tau[i]*err_del_v[i]
			d3 = d3**2

			df2 = (d1+d2+d3)*np.pi
			nhi_uncertainty = const*df2**0.5
			cnm_uncertainty = nhi_uncertainty
			wnm_uncertainty = 0.
			ret['cnm'].append(nhi[i])
			ret['wnm'].append(0.)
		else: # uncertainties for WNM component
			t1 = del_v[i]*err_tb[i]
			t1 = t1**2

			t2 = tb[i]*err_del_v[i]
			t2 = t2**2

			nhi_uncertainty = const*(np.pi*(t1+t2))**0.5
			cnm_uncertainty = 0.
			wnm_uncertainty = nhi_uncertainty
			ret['cnm'].append(0.)
			ret['wnm'].append(nhi[i])

		ret['nhi_uncertainty'].append(nhi_uncertainty)
		ret['cnm_uncertainty'].append(cnm_uncertainty)
		ret['wnm_uncertainty'].append(wnm_uncertainty)

	return ret

## Calculate the uncertainties of N(HI) for each source 
 # from Carl's paper
 #
 # params void
 #
 # return dict Gaussian-fit parameters of each component
 # including the uncertainties of each component
 # 
 # version 10/2016
 # Author Van Hiep ##
def cal_nhi_error():
	## Note: Cal. uncertainties of N(HI), sum from each component ##
	## N(HI) = sum[N(HI_i)] ##
	## sigma_NHI = sqrt[ sum[sigma_i^2] ] ##
	params = read_fit_params('result/component_fit_params.txt')
	er_nhi = {}
	er_cnm = {}
	nhi_er = {}

	for i in range(len(params['src'])):
		if params['src'][i] not in nhi_er.keys():
			nhi_er[params['src'][i]]             = {}
			nhi_er[params['src'][i]]['nhi']      = params['nhi'][i]
			nhi_er[params['src'][i]]['cnm']      = params['cnm'][i]
			nhi_er[params['src'][i]]['wnm']      = params['wnm'][i]

			nhi_er[params['src'][i]]['nhi_er_i'] = [params['nhi_uncertainty'][i]]
			nhi_er[params['src'][i]]['nhi_er']   = (params['nhi_uncertainty'][i])**2
			nhi_er[params['src'][i]]['cnm_er']   = (params['cnm_uncertainty'][i])**2
			nhi_er[params['src'][i]]['wnm_er']   = (params['wnm_uncertainty'][i])**2
		else:
			nhi_er[params['src'][i]]['nhi']      = nhi_er[params['src'][i]]['nhi'] + params['nhi'][i]
			nhi_er[params['src'][i]]['cnm']      = nhi_er[params['src'][i]]['cnm'] + params['cnm'][i]
			nhi_er[params['src'][i]]['wnm']      = nhi_er[params['src'][i]]['wnm'] + params['wnm'][i]

			nhi_er[params['src'][i]]['nhi_er_i'] += [params['nhi_uncertainty'][i]]
			nhi_er[params['src'][i]]['nhi_er']   = nhi_er[params['src'][i]]['nhi_er'] + (params['nhi_uncertainty'][i])**2
			nhi_er[params['src'][i]]['cnm_er']   = nhi_er[params['src'][i]]['cnm_er'] + (params['cnm_uncertainty'][i])**2
			nhi_er[params['src'][i]]['wnm_er']   = nhi_er[params['src'][i]]['wnm_er'] + (params['wnm_uncertainty'][i])**2

	j = 0
	for sc in nhi_er:
		nhi_er[sc]['nhi_er'] = nhi_er[sc]['nhi_er']**0.5
		nhi_er[sc]['cnm_er'] = nhi_er[sc]['cnm_er']**0.5
		nhi_er[sc]['wnm_er'] = nhi_er[sc]['wnm_er']**0.5
		# print '{}\t{}\t\t{}\t{}\t{}\t{}\t{}\t{}'\
		# .format(j, sc, nhi_er[sc]['nhi'], round(nhi_er[sc]['nhi_er'],2), nhi_er[sc]['cnm'], round(nhi_er[sc]['cnm_er'],4), nhi_er[sc]['wnm'], round(nhi_er[sc]['wnm_er'],4) )
		j = j +1

	return nhi_er

##================= MAIN ========================##
## Calculate the uncertainties of N(HI) for each source from Carl's paper ##
## 04 Oct 2016 ##
nhi_info = cal_nhi_error()

## Get infor from 26 sources without CO ##
cols  = ['indx','src','nhi_fk','err_fk','nh_pl','err_pl','nhi_hl','nhi_warm', 'nhi_cold', 'fk_hl', 'pl_hl']
fmt   = ['i','s','f','f','f','f','f','f','f','f','f']
fname = 'result/26nhi_nh_halfbeam_2errbar.txt'
dat   = restore(fname, 4, cols, fmt)
inf   = dat.read()
src   = inf['src'] # 26 no-CO-sources

## Write in order of 26-no-CO-sources ##
for i in range(len(src)):
	print('{}\t{}\t\t{}\t{}\t{}\t{}\t{}\t{}'.\
	# format( i, src[i], inf['nhi_hl'][i], inf['nhi_warm'][i], inf['nhi_cold'][i], er_nhi[src[i]], er_cnm[src[i]]) )
	# format( i, src[i], inf['nhi_hl'][i], round(nhi_info[src[i]]['nhi_er'],4), inf['nhi_cold'][i], round(nhi_info[src[i]]['cnm_er'],4), inf['nhi_warm'][i], round(nhi_info[src[i]]['wnm_er'],4)       ))
	format( i, src[i], inf['nhi_hl'][i], round(nhi_info[src[i]]['nhi_er'],4), inf['nhi_cold'][i], round(nhi_info[src[i]]['cnm_er'],4), inf['nhi_warm'][i], round(nhi_info[src[i]]['wnm_er'],4)       ))
