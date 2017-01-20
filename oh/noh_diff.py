import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import operator

from numpy    import array
from restore  import restore
from plotting import cplot

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
	for i in range(0,len(src)):
		if dat['src'][i] not in ret.keys():
			ret[src[i]] = {}
			ret[src[i]]['noh'] = noh[i]
			ret[src[i]]['er2'] = er2[i]**2
		else:
			ret[src[i]]['noh'] = ret[src[i]]['noh'] + noh[i]
			ret[src[i]]['er2'] = ret[src[i]]['er2'] + er2[i]**2

		print i, src[i], ret[src[i]]['noh']

	return ret

## plot diff of N(OH1) and N(OH2) #
 #
 # params dict noh1 Data on N(OH1)
 # params dict noh2 Data on N(OH2)
 # return void
 # 
 # version 09/2016 
 # author Nguyen Van Hiep ##
def hist_noh_diff(noh1,noh2):
	n_oh1 = []
	n_oh2 = []
	ind   = []
	src   = []
	count = 0
	for sc in noh1:
		print count, sc
		ind.append(count)
		src.append(sc)
		n_oh1.append(noh1[sc]['noh'])
		n_oh2.append(noh2[sc]['noh'])
		count = count + 1

	n_oh1 = np.asarray(n_oh1, dtype=np.float32)
	n_oh2 = np.asarray(n_oh2, dtype=np.float32)
	print n_oh1/n_oh2

	## Plot histogram ##
	size = 1.
	fig  = cplot()
	trace1 = fig.hist(n_oh1/n_oh2,label='[N(OH1665)/N(OH1667)]',autobinx=False,
	                  xbins=dict(start=-5.0, end=5.0, size=size),
	                  opacity=1.0,
	                  histtype='step',
	                  marker=dict(
	                    color = 'r',
	                    linewidth=2                   
	                    )
	                 )
	
	data   = [trace1]

	layout = dict(title  = 'Histogram of N(OH1665)/N(OH1667)',
					  title_fontsize=30,
		              grid   = True,
		              legend = dict(loc='upper right', fontsize=18),
		              xaxis  = dict(label='N(OH1665)/N(OH1667)',tick_size=18,fontsize=35),
		              yaxis  = dict(label='Counts',tick_size=18,fontsize=35),
		              # text   = [dict(loc=[0.2,0.4],text='a = '+str(m)+'$\pm$'+str(ea) +',  b = '+str(b)+'$\pm$'+str(eb),color='blue',fontsize=17),
		              # 			dict(loc=[0.2,0.3],text='(Available sources with the presence of OH line are shown)',color='red',fontsize=19)
		              # 		   ],
		             )
	fig.iplot(data,layout)

## plot N(OH) vs N(HI) #
 #
 # params dict dat Data on N(HI)
 # params dict noh1 Data on N(OH1)
 # params dict noh2 Data on N(OH2)
 # return void
 # 
 # version 09/2016 
 # author Nguyen Van Hiep ##
def plot_noh_diff(noh1,noh2):
	n_oh1 = []
	n_oh2 = []
	ind   = []
	src   = []
	count = 0
	diff  = []
	ratio = []
	for sc in noh1:
		print count, sc, noh1[sc]['noh'], noh2[sc]['noh']
		ind.append(count)
		src.append(sc)
		n_oh1.append(noh1[sc]['noh'])
		n_oh2.append(noh2[sc]['noh'])
		diff.append(noh2[sc]['noh']-noh1[sc]['noh'])
		ratio.append(noh2[sc]['noh']/noh1[sc]['noh'])
		count = count + 1

	n_oh1 = np.asarray(n_oh1, dtype=np.float32)
	n_oh2 = np.asarray(n_oh2, dtype=np.float32)
	diff  = np.asarray(diff, dtype=np.float32)
	ratio = n_oh1/n_oh2

	# plt.plot(ind, n_oh2-n_oh1)
	# plt.show()

	# Error bar for x-axis and y-axis
	xer = np.zeros(n_oh1.shape, dtype=np.float32)
	yer = np.zeros(n_oh1.shape, dtype=np.float32)

	fig    = cplot()
	trace1 = fig.error_bar(ind,ratio,label='N(OH1665)/N(OH1667)',
			err_x=dict(val=xer),
			err_y=dict(val=yer),			
			prop=dict(fmt='ro',
				markersize=8,
				markeredgecolor='b',
				markeredgewidth=1)
			)

	data   = [trace1]
	layout = dict(title  = 'Ratio N(OH1665)/N(OH1667)',
				  title_fontsize=30,
	              grid   = True,
	              legend = dict(loc='upper left', fontsize=18),
	              xaxis  = dict(label='Source indexes',tick_size=18,fontsize=35,xlim=[-1.,22.]),
	              yaxis  = dict(label='N(OH1665)/N(OH1667)',tick_size=18,fontsize=35,ylim=[-.1,3.]),
	              # text   = [dict(loc=[0.2,0.4],text='a = '+str(m)+'$\pm$'+str(ea) +',  b = '+str(b)+'$\pm$'+str(eb),color='blue',fontsize=17),
	              # 			dict(loc=[0.2,0.3],text='(Available sources with the presence of OH line are shown)',color='red',fontsize=19)
	              # 		   ],
	             )

	for i in range(len(src)):
		fig.annotation(dict(label='('+src[i]+')',
			x=ind[i],y=ratio[i],
			fontsize=12,
			xycoords='data',
			xoff=-50.,yoff=20.,
			textcoords='offset points',
			arrowprops=dict(arrowstyle="->") )
			)

	fig.iplot(data,layout)

	#########################
	trace1 = fig.error_bar(n_oh1,ratio,label='N(OH1665)/N(OH1667)',
			err_x=dict(val=xer),
			err_y=dict(val=yer),			
			prop=dict(fmt='ro',
				markersize=8,
				markeredgecolor='b',
				markeredgewidth=1)
			)

	data   = [trace1]
	layout = dict(title  = 'Ratio N(OH1665)/N(OH1667)',
				  title_fontsize=30,
	              grid   = True,
	              legend = dict(loc='upper left', fontsize=18),
	              xaxis  = dict(label='N(OH1665)',tick_size=18,fontsize=35,xlim=[-1.,5.]),
	              yaxis  = dict(label='N(OH1665)/N(OH1667)',tick_size=18,fontsize=35,ylim=[-.1,3.]),
	              # text   = [dict(loc=[0.2,0.4],text='a = '+str(m)+'$\pm$'+str(ea) +',  b = '+str(b)+'$\pm$'+str(eb),color='blue',fontsize=17),
	              # 			dict(loc=[0.2,0.3],text='(Available sources with the presence of OH line are shown)',color='red',fontsize=19)
	              # 		   ],
	             )

	for i in range(len(src)):
		fig.annotation(dict(label='('+src[i]+')',
			x=n_oh1[i],y=ratio[i],
			fontsize=12,
			xycoords='data',
			xoff=-50.,yoff=20.,
			textcoords='offset points',
			arrowprops=dict(arrowstyle="->") )
			)

	fig.iplot(data,layout)

	#########################
	trace1 = fig.error_bar(n_oh2,ratio,label='N(OH1665)/N(OH1667)',
			err_x=dict(val=xer),
			err_y=dict(val=yer),			
			prop=dict(fmt='ro',
				markersize=8,
				markeredgecolor='b',
				markeredgewidth=1)
			)

	data   = [trace1]
	layout = dict(title  = 'Ratio N(OH1665)/N(OH1667)',
				  title_fontsize=30,
	              grid   = True,
	              legend = dict(loc='upper left', fontsize=18),
	              xaxis  = dict(label='N(OH1667)',tick_size=18,fontsize=35,xlim=[-1.,6.]),
	              yaxis  = dict(label='N(OH1665)/N(OH1667)',tick_size=18,fontsize=35,ylim=[-.1,3.]),
	              # text   = [dict(loc=[0.2,0.4],text='a = '+str(m)+'$\pm$'+str(ea) +',  b = '+str(b)+'$\pm$'+str(eb),color='blue',fontsize=17),
	              # 			dict(loc=[0.2,0.3],text='(Available sources with the presence of OH line are shown)',color='red',fontsize=19)
	              # 		   ],
	             )

	for i in range(len(src)):
		fig.annotation(dict(label='('+src[i]+')',
			x=n_oh2[i],y=ratio[i],
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

hist_noh_diff(noh1,noh2)
plot_noh_diff(noh1,noh2)

sys.exit()