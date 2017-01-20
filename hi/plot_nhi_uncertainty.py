import os, sys, shutil
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt
import pylab             as pl
import copy

from numpy               import array
from restore             import restore
import operator

## Read N(HI) of 78 sources and its uncertainty #
 #
 # params string fname Filename
 # return dict info of N(H) and N(HI) & uncertaintes
 # 
 # version 10/2016
 # Author Van Hiep##
def read_nhi_uncertainty(fname = ''):
	cols  = ['indx','src','nhi','nhi_er','cnm','cnm_er','wnm','wnm_er']
	fmt   = ['f',    's',  'f',  'f',      'f', 'f',     'f',    'f']
	dat   = restore(fname, 3, cols, fmt)
	ret   = dat.read()
	return ret

## Plot histogram of N(HI) uncertainties - 78 sources #
 #
 # params dict data data of uncertainties
 # return void
 #
 # version 10/2016
 # Author Van Hiep ##
def plot_hist_nhi_err_78src(data):
	src    = data['src']
	nhi    = np.asarray(data['nhi'])
	nhi_er = np.asarray(data['nhi_er'])
	a      = 100.*nhi_er/nhi

	b = list(np.where(a>50.) ) ## 3 sources with high Uncertainties, in Galactic plane
	b = b[0]
	string = ''
	for i in b:
		string = string + src[i] + ', '

	print string

	bins = np.histogram(np.hstack((a)), bins=25)[1] #get the bin edges

	plt.hist(a, bins, label='Uncertainties of N(HI)', histtype='step', linewidth=3)

	plt.xlabel('Uncertainties (%)', fontsize = 35)
	plt.ylabel('Counts', fontsize = 35)
	plt.title('Histogram of total $N_{HI}$ uncertainties for 78 lines-of-sight', fontsize = 35)
	plt.grid(True)
	plt.tick_params(axis='x', labelsize=20)
	plt.tick_params(axis='y', labelsize=20)
	plt.text(55, 5, string, color='b', fontsize=14)
	plt.text(55, 35, 'Mean = ' + str(round(a.mean(),2) ) + '%',color='b',fontsize=14)

	plt.legend(loc='upper right', fontsize = 18)
	plt.show()

## Plot histogram of N(HI) uncertainties - 26 no-CO sources #
 #
 # params dict data data of uncertainties
 # return void
 #
 # version 10/2016
 # Author Van Hiep ##
def plot_hist_nhi_err_26src_no_co(data):
	src    = data['src']
	nhi    = np.asarray(data['nhi'])
	nhi_er = np.asarray(data['nhi_er'])
	a      = 100.*nhi_er/nhi

	b = list(np.where(a>50.) )
	b = b[0]
	string = ''
	for i in b:
		string = string + src[i] + ', '

	print string

	bins = np.histogram(np.hstack((a)), bins=25)[1] #get the bin edges

	plt.hist(a, bins, label='Uncertainties of N(HI)', histtype='step', linewidth=2)

	plt.xlabel('Uncertainties (%)', fontsize = 35)
	plt.ylabel('Counts', fontsize = 35)
	plt.title('Histogram of $N_{HI}$ uncertainties for 26 sources without CO', fontsize = 35)
	plt.grid(True)
	plt.tick_params(axis='x', labelsize=20)
	plt.tick_params(axis='y', labelsize=20)
	# plt.text(55, 5,string,color='b',fontsize=14)
	# plt.text(55, 35, 'Mean = ' + str(round(a.mean(),2) ) + '%',color='b',fontsize=14)

	plt.legend(loc='upper right', fontsize = 18)
	plt.show()

# plot the N(HI) uncertainties of 26 sources  #
#
# params dict data data to plot
#
# return void
# 
# Author Van Hiep
##
def plot_nhi_uncertainty(data):
	x    = data['id']
	#y    = map(abs, data['nhi_diff'])
	y    = data['err']
	plt.plot(x, y, 'r*-')
	plt.grid()
	plt.title('N(HI) uncertainties')
	plt.ylabel('Uncertainties (%)')
	plt.xlabel('index')
	plt.xlim(-0., 26)
	plt.ylim(-1, 25)
	for i,j in zip(x,y):
		k = j
		# if(i==12):
		# 	k = k-0.5
		if (j>-1.0) :
			plt.annotate(str(j)+'% ('+str(data['src'][i])+')',xy=(i,k)) # with offset plt.annotate(str(j),xy=(i,j+0.5))
	plt.show()

#================= MAIN ========================#

data = read_nhi_uncertainty('result/78src_nhi_cnm_wnm_with_err.txt')
plot_hist_nhi_err_78src(data)

sys.exit()

data = read_nhi_uncertainty('result/26src_no_co_cnm_uncertainties_arcb.txt')
plot_hist_nhi_err_26src_no_co(data)
#plot_nhi_uncertainty(data)