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
def read_dust_nh_nhi(fname = ''):
	cols = ['indx', 'src', 'l', 'b', 'nhi', 'nhi_er', 'nh_d', 'nhd_er']
	fmt  = ['i',     's',   'f','f',  'f',   'f',     'f',     'f']
	data = restore(fname, 3, cols, fmt)
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
	data    = read_dust_nh_nhi(fname = '../result/nh_nhi_uncert_78src.txt')
	src1    = data['src']
	idx1    = np.asarray(data['indx'])
	nhi1    = np.asarray(data['nhi'])
	nhi_er1 = np.asarray(data['nhi_er'])
	
	nh1     = np.asarray(data['nh_d'])
	nh_er1  = np.asarray(data['nhd_er'])

	x1 = np.log10(nhi1)
	y1 = nh1/nhi1

	# Error bar for x-axis and y-axis
	xerr1 = nhi_er1/nhi1/np.log(10.0)
	yerr1 = uncertainty_of_factors(y1, nhi1, nhi_er1, nh1, nh_er1)

	# Linear Fit#
	m,b,ea,eb = linear_fit(x1,y1)

	# Plot #
	fig    = cplot()	
	trace2 = fig.error_bar(x1,y1,label='Factor $f = N_{H}$/$N_{HI}$, $N_{H}$ from thermal dust',
			err_x=dict(val=xerr1),
			err_y=dict(val=yerr1),			
			prop=dict(fmt='ro',
				markersize=8,
				markeredgecolor='b',
				markeredgewidth=1)
			)
	trace1 = fig.lines(x1,m*x1+b,label='Best linear fit',
			prop=dict(color='purple',
				      linewidth=3,
				      linestyle='solid',
				      marker='o',
				      markerfacecolor='b',
				      markersize=0
			))

	m  = round(m,2)
	b  = round(b,2)
	ea = round(ea,2)
	eb = round(eb,2)

	data   = [trace2,trace1]
	layout = dict(title  = 'Correlation between Factor $f = N_{H}$/$N_{HI}$ and HI column density $N_{HI}$ \nalong 78 lines-of-sight',
				  title_fontsize=30,
	              grid   = True,
	              legend = dict(loc='upper left', fontsize=18), #nh_i = tau353[i]/planck_cf
	              xaxis  = dict(label='log$_{10}$($N_{HI}/10^{20}$ cm$^{-2}$)',tick_size=18,fontsize=35,xlim=[0.,2.5]),
	              yaxis  = dict(label='$Factor f = N_{H}$/$N_{HI}$',tick_size=18,fontsize=35,ylim=[0.,11.0]),
	              text   = [dict(loc=[0.01,4.0],text='a = '+str(m)+'$\pm$'+str(ea) +',  b = '+str(b)+'$\pm$'+str(eb),color='blue',fontsize=17),
	              			dict(loc=[0.01,9.0],text='$tau_{353} = [8.4\pm3.0]10^{-27}N_{H} $',color='blue',fontsize=17),
	              			# dict(loc=[0.2,0.2],text='(Available sources with the presence of OH line are shown)',color='red',fontsize=19)
	              		   ],
	             )

	# for i in range(len(src1)):
	# 	# if (oh[i] > 0) :
	# 	fig.annotation(dict(label='('+str(src1[i])+')',
	# 		x=x1[i],y=y1[i],
	# 		fontsize=18,
	# 		xycoords='data',
	# 		xoff=-50.,yoff=30.,
	# 		textcoords='offset points',
	# 		arrowprops=dict(arrowstyle="->") )
	# 		)

	fig.iplot(data,layout)

	plt.plot(nhi1,nh1, 'rd', label='data no CO', ms=10)
	plt.plot([0,130],[0,130], 'k--', label='$N_{H} = N_{HI}$')
	plt.title('Correlation between $N_{H}$ and $N_{HI}$ \nalong 78 lines-of-sight', fontsize=30)
	plt.ylabel('$N_{H}[10^{20}$ cm$^{-2}]$', fontsize=35)
	plt.xlabel('$N_{HI} [10^{20}$ cm$^{-2}]$', fontsize=35)
	# plt.xlim(0, 1.6)
	# plt.ylim(0, 3)
	plt.grid(True)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=18)

	plt.text(15., 2., '(Available sources with the presence of OH are shown)', color='k', fontsize=17)
	plt.text(15., 3., r'$N_{H} = 5.8\cdot10^{21}[cm^{-2}mag^{-1}]\cdot E(B-V)$', color='k', fontsize=17)
	plt.text(15., 4., r'E(B-V) from Planck data R1.2', color='k', fontsize=17)

	plt.legend(loc='upper left', fontsize=18)
	# plt.savefig("test.png",bbox_inches='tight')
	for i in range(len(src1)):
		# if (oh[i] > 0) :
		plt.annotate('('+str(src1[i])+')', xy=(nhi1[i], nh1[i]), xycoords='data',
               xytext=(-50.,30.), textcoords='offset points',
               arrowprops=dict(arrowstyle="->"),fontsize=18,
               )
	plt.show()

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