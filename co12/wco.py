import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import operator

from numpy    import array
from restore  import restore
from plotting import cplot

# Read Ts #
 #
 # params string fname Filename
 #
 # return dict info of Ts
 # 
 # Author Van Hiep
 ##
def read_vrange(fname = 'data/co_src_peak_vrange.txt'):
	cols = ['file','v1','v2', 'src']
	fmt  = ['s',    'f', 'f',  's' ]
	data = restore(fname, 2, cols, fmt)
	return data.read()

## Read and plot the reduced-12co(10) spectrum #
 #
 # params string fname Filename
 #
 # return void
 # 
 # Author Van Hiep
 ##
def read_co_spec(fname):
	cols = ['v','T']
	fmt  = ['f','f']
	data = restore(fname, 0, cols, fmt)
	dat  = data.read()

	return dat['v'],dat['T']

## Get the index of a given velocity #
 #
 # params list v-axis List of Velocity_axis
 # params float vel Velocity
 #
 # return int idx Index of vel in List of velocity_axis
 # 
 # Author Van Hiep
 ##
def get_vel_index(v_axis, vel):
    idx = (np.abs(np.array(v_axis)-vel)).argmin()
    return idx

## Plot spectra of CO(1-0) #
 #
 # params 
 # return void
 # 
 # version 12/2016
 # Author Van Hiep  ##
def plot_spec():
	info = read_vrange()
	src  = info['src']
	fil  = info['file']
	v1   = info['v1']
	v2   = info['v2']

	for i in range(len(src)):

		# if (src[i] != 'T0629+10'):
		# 	continue

		v,T = read_co_spec('data/' + fil[i])
		
		dv   = round((v[-1]-v[0])/len(v), 2)
		dv   = (v[-1]-v[0])/len(v)

		n1 = get_vel_index(v, v1[i])
		n2 = get_vel_index(v, v2[i])

		n1 = n1 - 40
		n2 = n2 + 40

		plt.plot(v[n1:n2], T[n1:n2], 'b-', lw=2)
		plt.grid()
		plt.title(src[i], fontsize = 35)
		plt.ylabel('$T_{b} (K)$', fontsize = 35)
		plt.xlabel('VLSR (km/s)', fontsize = 35)
		plt.tick_params(axis='x', labelsize=20)
		plt.tick_params(axis='y', labelsize=20)
		# plt.text(55, 5,string,color='b',fontsize=14)
		# plt.legend(loc='upper right', fontsize = 18)
		# plt.axvline(x=60., lw=4)
		# plt.axvline(x=-90., lw=4)
		plt.show()

		print src[i], v1[i], v2[i]

## Plot spectra of CO(1-0) #
 #
 # params 
 # return void
 # 
 # Author Van Hiep ##
def cal_wco():
	info = read_vrange()
	src  = info['src']
	fil  = info['file']
	v1   = info['v1']
	v2   = info['v2']

	for i in range(len(src)):

		v,T = read_co_spec('data/' + fil[i])
		
		dv  = round((v[-1]-v[0])/len(v), 2)
		dv  = (v[-1]-v[0])/len(v)

		n1  = get_vel_index(v, v1[i])
		n2  = get_vel_index(v, v2[i])

		t   = np.array(T[n1:n2], dtype=np.float64)
		wco = t.sum()*dv
		wco = round(wco, 10)

		x   = []
		y   = []
		for j in range(len(v)):
			if( ((v[j]>v1[i]-50.) and (v[j]<v1[i])) or ((v[j]>v2[i]) and (v[j]<v2[i]+50.)) ):
				x.append(v[j])
				y.append(T[j])

		x = np.asarray(x, dtype=np.float32)
		y = np.asarray(y, dtype=np.float32)

		ery  = np.std(y)
		n    = n2-n1
		erw2 = n*(ery**2)*(dv**2)
		erw  = np.sqrt(erw2)

		# print np.std(y), ery, wco, n, wco, erw/wco, dv, erw

		print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(i, src[i], v1[i], v2[i], wco, erw ))

## Plot spectra of CO(1-0) #
 #
 # params 
 # return void
 # 
 # Author Van Hiep ##
def plot_full_spec():
	info = read_vrange()
	src  = info['src']
	fil  = info['file']
	v1   = info['v1']
	v2   = info['v2']

	for i in range(len(src)):

		# if (src[i] != 'T0629+10'):
		# 	continue

		v,T = read_co_spec('data/' + fil[i])
		
		dv   = round((v[-1]-v[0])/len(v), 2)
		dv   = (v[-1]-v[0])/len(v)

		n1 = get_vel_index(v, v1[i])
		n2 = get_vel_index(v, v2[i])

		plt.plot(v, T, 'r')
		plt.ylim(-5., 10.)
		plt.xlim(-500., 500.)
		plt.grid()
		plt.title('Spectrum: ' + src[i])
		plt.ylabel('T(K)')
		plt.xlabel('v(km/s')
		plt.show()

		print src[i], v1[i], v2[i]

## Plot spectra of CO(1-2) #
 #
 # params 
 # return void
 # 
 # Author Van Hiep
 ##
def get_1sigma_co_spec():
	info = read_vrange()
	src  = info['src']
	v1   = info['v1']
	v2   = info['v2']

	for i in range(len(src)):
		v,T = read_co_spec('data/' + src[i])

		x   = []
		y   = []
		for j in range(len(v)):
			if( ((v[j]>v1[i]-50.) and (v[j]<v1[i])) or ((v[j]>v2[i]) and (v[j]<v2[i]+50.)) ):
				x.append(v[j])
				y.append(T[j])

		x = np.asarray(x, dtype=np.float32)
		y = np.asarray(y, dtype=np.float32)

		print np.std(y)

		plt.plot(v,T, 'b.')
		plt.plot(x,y,'r.')
		plt.xlim(v1[i]-60., v2[i]+60.)
		plt.ylim(-5.,10.)
		plt.show()

############## MAIN ###########
plot_spec()
# cal_wco()
# plot_full_spec()
# get_1sigma_co_spec()
sys.exit()