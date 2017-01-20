import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import operator

from numpy    import array
from restore  import restore
from plotting import cplot

# Calculate the uncertainties of factors #
#
# params 
# params 
# params 
#
# return list of factor_uncertainties
# 
# Author Van Hiep
##
def uncertainty_of_factors(factr, nhi_hl, err_hl, nhi, err):
	d1 = np.array(err_hl)/np.array(nhi_hl)
	d1 = np.square(d1)

	d2 = np.array(err)/np.array(nhi)
	d2 = np.square(d2)

	d  = np.sqrt(d1+d2)*np.array(factr)

	return d.tolist()

# linear fit #
#
# params x list x-data
# params y list y-data
#
# return fit parameters and errors
# 
# Author Van Hiep
##
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

## Divide two lists #
 #
 # params x list 
 # params y list 
 #
 # return x/y
 # 
 # Author Van Hiep
 ##
def list_div(x,y):
	ret = []
	for i in range(0, len(x)):
		ret.append(x[i]/y[i])

	return ret

## Compute N(OH) ##
 #
 # params dict dat Data on N(OH)
 # return dict ret N(OH) of each source
 #
 # version 09/2016 
 # author Nguyen Van Hiep ##
def get_noh_for_each_src(fname):
	cols = ['idx','src','tau','v0','wid','tex','tex_er','noh','noh_er']
	fmt  = ['i','s','f','f','f','f','f','f','f']
	data = restore(fname, 2, cols, fmt)
	dat  = data.read()
	noh  = dat['noh']
	er2  = dat['noh_er']
	src  = dat['src']

	ret  = {}
	for i in range(0,len(dat['src'])):
		if dat['src'][i] not in ret.keys():
			ret[src[i]] = {}
			ret[src[i]]['noh'] = noh[i]
			ret[src[i]]['er2'] = er2[i]**2
		else:
			ret[src[i]]['noh'] = ret[src[i]]['noh'] + noh[i]
			ret[src[i]]['er2'] = ret[src[i]]['er2'] + er2[i]**2

	return ret

## plot N(OH) vs N(HI) #
 #
 # params dict dat Data on N(HI)
 # params dict noh1 Data on N(OH1)
 # params dict noh2 Data on N(OH2)
 # return void
 # 
 # version 09/2016 
 # author Nguyen Van Hiep ##
def plot_oh_vs_hi(dat,noh1,noh2):
	nhi    = dat['nhi']
	nhi_er = dat['nhi_er']
	src    = dat['src']
	idx    = dat['idx']

	nhi1  = []
	nhi2  = []
	n_hi1 = []
	n_hi2 = []
	nohlog1 = []
	nohlog2 = []
	noh65   = []
	noh67   = []
	yerr1 = [] # Sigma**2
	yerr2 = [] # Sigma**2
	xer1  = []
	xer2  = []
	for sc in noh1:
		if((sc in src) and (noh1[sc] > 0.)):
			nohlog1.append(np.log10(noh1[sc]['noh']))
			noh65.append(noh1[sc]['noh'])
			yerr1.append(noh1[sc]['er2']) # Sigma**2
			n_hi1.append(np.log10(nhi[src.index(sc)]))
			nhi1.append(nhi[src.index(sc)])
			xer1.append(nhi_er[src.index(sc)])

	for sc in noh2:
		if((sc in src) and (noh2[sc] > 0.)):
			nohlog2.append(np.log10(noh2[sc]['noh']))
			noh67.append(noh2[sc]['noh'])
			yerr2.append(noh2[sc]['er2']) # Sigma**2
			n_hi2.append(np.log10(nhi[src.index(sc)]))
			nhi2.append(nhi[src.index(sc)])
			xer2.append(nhi_er[src.index(sc)])

	# Error bar for x-axis and y-axis
	xer1 = np.array(xer1)/np.array(nhi1)/np.log(10.0)
	xer2 = np.array(xer2)/np.array(nhi2)/np.log(10.0)
	yer1 = np.array(np.sqrt(yerr1))/np.array(noh65)/np.log(10.0)
	yer2 = np.array(np.sqrt(yerr2))/np.array(noh67)/np.log(10.0)

	fig    = cplot()	
	trace2 = fig.error_bar(n_hi1,nohlog1,label='log(N(OH1665)/1e14)',
			err_x=dict(val=xer1),
			err_y=dict(val=yer1),			
			prop=dict(fmt='bo',
				markersize=8,
				markeredgecolor='b',
				markeredgewidth=1)
			)
	trace1 = fig.error_bar(n_hi2,nohlog2,label='log(N(OH1667)/1e14)',
			err_x=dict(val=xer2),
			err_y=dict(val=yer2),			
			prop=dict(fmt='ro',
				markersize=8,
				markeredgecolor='r',
				markeredgewidth=1)
			)

	data   = [trace2,trace1]
	layout = dict(title  = 'Correlation between N(OH) and N(HI)',
				  title_fontsize=30,
	              grid   = True,
	              legend = dict(loc='upper right', fontsize=18),
	              xaxis  = dict(label='log(N(HI)/1e20)',tick_size=18,fontsize=35,xlim=[0.5,2.5]),
	              yaxis  = dict(label='log(N(OH)/1e14)',tick_size=18,fontsize=35,ylim=[-5.,2.]),
	              # text   = [dict(loc=[0.2,0.4],text='a = '+str(m)+'$\pm$'+str(ea) +',  b = '+str(b)+'$\pm$'+str(eb),color='blue',fontsize=17),
	              # 			dict(loc=[0.2,0.3],text='(Available sources with the presence of OH line are shown)',color='red',fontsize=19)
	              # 		   ],
	             )

	for sc in noh1:
		if((sc in src) and (noh1[sc] > 0.)):
			fig.annotation(dict(label='('+sc+')',
				x=np.log10(nhi[src.index(sc)]),y=np.log10(noh1[sc]['noh']),
				fontsize=12,
				xycoords='data',
				xoff=-50.,yoff=20.,
				textcoords='offset points',
				arrowprops=dict(arrowstyle="->") )
				)

	fig.iplot(data,layout)

#================= MAIN ========================#

## N(HI) ##
cols = ['idx','src','nhi','nhi_er']
fmt  = ['i','s','f','f']
data = restore('sub_data/nhi_heiles_uncert_78_src.txt', 2, cols, fmt)
dat  = data.read()

## N(OH) ##
noh1 = get_noh_for_each_src(fname = 'result/noh1_src96_er.txt')
noh2 = get_noh_for_each_src(fname = 'result/noh2_src96_er.txt')

plot_oh_vs_hi(dat,noh1,noh2)

sys.exit()