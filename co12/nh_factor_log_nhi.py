import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import operator

from numpy    import array
from restore  import restore
from plotting import cplot

## Read nh (derived from CO), nhi and errors #
 #
 # params string fname Filename
 #
 # return dict info of N(H)
 # 
 # version 9/2016
 # Author Van Hiep ##
def read_nh_nhi(fname = ''):
	cols = ['indx', 'src','file', 'v1', 'v2', 'wco', 'wco_er', 'nhi', 'nhi_er', 'nh2', 'nh2_er', 'nh', 'nh_er']
	fmt  = ['i','s','s','f','f','f','f','f','f','f','f','f','f']
	data = restore(fname, 2, cols, fmt)
	return data.read()

##  Read N(HI)_fukui N(H)_Planck N(HI)_Heiles #
 #
 # params string fname Filename
 #
 # return dict info of N(H)
 # 
 # version 9/2016
 # Author Van Hiep ##
def read_planck_nh_nhi(fname = ''): 
	cols = ['idx','src','nhi_fk','err_fk','nh_pl','err_pl','nhi_hl','er_hl','err_hl','wnm','cnm','fk_hl','pl_hl','oh','nhi_thin']
	fmt  = ['i','s','f','f','f','f','f','f','f','f','f','f','f','f','f']
	data = restore(fname, 4, cols, fmt)
	return data.read()	

## Plot the correlation between Factors and log(N(HI)) only from CO #
 # f = n(H)/n(HI)
 #
 # params 
 # return Void
 # 
 # version 9/2016
 # Author Van Hiep ##
def nh_factor_vs_nhi():
	data    = read_nh_nhi(fname = 'result/nh_with_er_16src.txt')
	src1    = data['src']
	idx1    = np.asarray(data['indx'])
	nhi1    = np.asarray(data['nhi'])
	nhi_er1 = np.asarray(data['nhi_er'])
	nh1     = np.asarray(data['nh'])
	nh_er1  = np.asarray(data['nh_er'])

	x1 = np.log10(nhi1)
	y1 = nh1/nhi1

	# Error bar for x-axis and y-axis
	xerr1 = nhi_er1/nhi1/np.log(10.0)
	yerr1 = uncertainty_of_factors(y1, nhi1, nhi_er1, nh1, nh_er1)

	dat     = read_planck_nh_nhi(fname = 'sub_data/nhi_and_uncertainties_full.txt')
	src2    = dat['src']
	nhi2    = np.asarray(dat['nhi_hl'])
	nhi_er2 = np.asarray(dat['err_hl'])
	nh2     = np.asarray(dat['nh_pl'])
	nh_er2  = np.asarray(dat['err_pl'])

	## Append 2 set of data ##
	src    = np.append(src1, src2)
	nhi    = np.append(nhi1, nhi2 )
	nhi_er = np.append(nhi_er1, nhi_er2)
	nh     = np.append(nh1, nh2)
	nh_er  = np.append(nh_er1, nh_er2)

	x      = np.log10(nhi)
	y      = nh/nhi
	## End ##


	x2 = np.log10(nhi2)
	y2 = nh2/nhi2

	xerr2 = nhi_er2/nhi2/np.log(10.0)
	yerr2 = uncertainty_of_factors(y2, nhi2, nhi_er2, nh2, nh_er2)

	# Fit and Plot #
	m,b,ea,eb = linear_fit(x,y)

	fig    = cplot(filename='mmm')	
	trace2 = fig.error_bar(x2,y2,label='Factor $f = N_{H}$/$N_{HI}$, $N_{H}$ from thermal dust',
			err_x=dict(val=xerr2),
			err_y=dict(val=yerr2),			
			prop=dict(fmt='ro',
				markersize=8,
				markeredgecolor='b',
				markeredgewidth=1)
			)
	trace1 = fig.lines(x,m*x+b,label='Best linear fit',
			prop=dict(color='purple',
				      linewidth=3,
				      linestyle='solid',
				      marker='o',
				      markerfacecolor='b',
				      markersize=0
			))

	trace3 = fig.error_bar(x1,y1,label='Factor $f = N_{H}$/$N_{HI}$, $N_{H}$ from CO(1-0)',
			err_x=dict(val=xerr1),
			err_y=dict(val=yerr1),			
			prop=dict(fmt='b^',
				markersize=8,
				markeredgecolor='b',
				markeredgewidth=1)
			)

	m  = round(m,2)
	b  = round(b,2)
	ea = round(ea,2)
	eb = round(eb,2)

	data   = [trace2,trace1, trace3]
	layout = dict(title  = 'Correlation between Factor $f = N_{H}$/$N_{HI}$ and HI column density $N_{HI}$ \nalong 42 lines-of-sight',
				  title_fontsize=30,
	              grid   = True,
	              legend = dict(loc='upper left', fontsize=18),
	              xaxis  = dict(label='log$_{10}$($N_{HI}/10^{20}$ cm$^{-2}$)',tick_size=18,fontsize=35,xlim=[0.,2.5]),
	              yaxis  = dict(label='$Factor f = N_{H}$/$N_{HI}$',tick_size=18,fontsize=35,ylim=[0.,6.0]),
	              text   = [dict(loc=[0.2,0.4],text='a = '+str(m)+'$\pm$'+str(ea) +',  b = '+str(b)+'$\pm$'+str(eb),color='blue',fontsize=17),
	              			dict(loc=[1.8,5.4],text='$X = [2.54\pm0.13]10^{20} [H_{2}cm^{-2}/(K km/s)]$',color='blue',fontsize=17),
	              			dict(loc=[0.2,0.2],text='(Available sources with the presence of OH line are shown)',color='red',fontsize=19)
	              		   ],
	             )

	for i in range(len(src)):
		# if (oh[i] > 0) :
		fig.annotation(dict(label='('+str(src[i])+')',
			x=x[i],y=y[i],
			fontsize=18,
			xycoords='data',
			xoff=-50.,yoff=30.,
			textcoords='offset points',
			arrowprops=dict(arrowstyle="->") )
			)

	fig.iplot(data,layout)

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

# plot N(HI)/N(HI)_thin vs N(HI)_thin #
#
# params dict data data to plot
#
# return void
# 
# Author Van Hiep
##
def plot_factor_vs_nhi(data):

	fact       = []
	lognhi     = []
	nhi        = data['nhi_thin'] # Optically-thin assumption
	nhi_heiles = data['nhi_hl'] # N(HI) from Carl Heiles Paper

	for i in range(0, len(nhi)):
		temp   = round(nhi_heiles[i]/nhi[i], 3)
		lognhi.append(np.log10(nhi[i]))
		fact.append(temp)

	# Fit and Plot #
	a, b = np.polyfit(lognhi,fact,1)

	plt.plot(lognhi, fact, 'r.', label='Factor = N(HI)/N(HI)_thin')
	plt.plot(lognhi, a*np.array(lognhi) + b, 'k-', label='')

	a = round(a, 2)
	b = round(b, 2)

	plt.xlabel('Factor')
	plt.ylabel('log(N(HI)_opt.thin).1e20')
	plt.title('Correlation between N(HI)_Heiles and Optically thin N(HI)')
	plt.grid(True)

	plt.text(0.21, 1.31, 'Arecibo beam at 1.4GHz = 3.5\'', color='blue', fontsize=12)
	plt.text(0.21, 0.92, 'a = '+str(a)+'  b = '+str(b), color='blue', fontsize=12)

	plt.legend(loc='upper right')
	plt.show()

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
nh_factor_vs_nhi()