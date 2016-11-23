import numpy             as np
import matplotlib.pyplot as plt
import operator
from scipy.stats import norm
import matplotlib.mlab as mlab

from copy import copy

# =========== Define Functions =============== #
## Read and plot the difference between the paper's N(HI) and obtained-N(HI) #
 #
 # params string fname Filename
 #
 # return void
 #  
 # Author Van Hiep
 ##
def read_comp_nhi(fname = "result/nhi2comp_with_paper_20160316.txt"):
	ret = {}

	ret['idx']      = []
	ret['vs']       = []
	ret['ve']       = []
	ret['vs_id']    = []
	ret['ve_id']    = []

	ret['nhi_i']     = []
	ret['sources']   = []
	ret['ma_nhi']    = []
	ret['nhi_diff']  = []

	file    = open (fname,"r")
	file.readline() # comment line
	file.readline() # comment line
	for line in file:
	    line    = line.strip()
	    columns = line.split()

	    ret['idx'].append(int(columns[0]))
	    ret['vs'].append(float(columns[1]))
	    ret['ve'].append(float(columns[2]))
	    ret['vs_id'].append(float(columns[3]))
	    ret['ve_id'].append(float(columns[4]))
	    ret['nhi_i'].append(float(columns[5]))
	    ret['ma_nhi'].append(float(columns[6]))
	    ret['nhi_diff'].append(float(columns[7]))
	    ret['sources'].append(columns[8])

	file.close()

	return ret

## Read data from file #
 #
 # params string file filename
 #
 # return dict ret data
 # 
 # Author Van Hiep
 ##
def read_data_from_file(file = "data/specs.txt"):
	ret  = {}

	name = []
	v    = []
	t    = []

	f    = open (file,"r")
	#read line into array 
	for line in f:
	    line    = line.strip()
	    columns = line.split()

	    name.append(columns[0])
	    v.append(float(columns[1]))
	    t.append(float(columns[2]))

	f.close()

	ret['name'] = name
	ret['v']    = v
	ret['t']    = t

	return ret

## Read spectra #
 #
 # params string file filename
 #
 # return dict specs spectra
 # 
 # Author Van Hiep
 ##
def read_specs_data():
	a     = read_data_from_file()
	specs = {}

	name  = a['name']
	v     = a['v']
	t     = 0.5*np.array(a['t'])

	for i in range(0,79):
		specs[i] = {}

		start = i*2048
		end   = start + 2047
		
		specs[i]['name'] = name[end]
		specs[i]['v']    = v[start:end+1]
		specs[i]['t']    = t[start:end+1]

	return specs

## plot spectrum and calculate RMS of Gaussian noise #
 #
 # params dict data data
 # params string src source-name
 #
 # return rms
 # 
 # Author Van Hiep
 ##
def plot_specs(data, src):
	# Get key from Source_name #
	get_keys(data, src)
	i = result[0][0]

	nbin = len(data[i]['v'])
	t = []
	for k in range(0, nbin) :
		if (data[i]['v'][k] > 100. or data[i]['v'][k]<-100.) :
			t.append(data[i]['t'][k])

	(mu, sigma) = norm.fit(t)

	n, bins, patches = plt.hist(t, 60, normed=1, facecolor='green', alpha=0.75)

	y = mlab.normpdf( bins, mu, sigma)
	l = plt.plot(bins, y, 'r--', linewidth=2)

	print sigma, rms

	plt.grid()
	plt.title('Source ' + src)
	plt.ylabel('Tb(K)')
	plt.xlabel('v(km/s)')
	plt.show()

## plot spectrum and calculate RMS of Gaussian noise #
 #
 # params dict data data
 #
 # return rms
 # 
 # Author Van Hiep
 ##
def get_rms(data):
	ret = read_comp_nhi() # Nothing but just to get more infor from 79 src

	nsrc = 79
	nbin = 2048
	dv   = 0.161 # km/s
	for i in range(0,nsrc) :
		t = []
		for k in range(0, nbin) :
			if ((data[i]['v'][k]>ret['ve'][i]) or (data[i]['v'][k]<ret['vs'][i])) :
				t.append(data[i]['t'][k])

		(mu, sigma) = norm.fit(t)

		n, bins, patches = plt.hist(t, 60, normed=1, facecolor='green', alpha=0.75)

		y = mlab.normpdf(bins, mu, sigma)
		l = plt.plot(bins, y, 'r--', linewidth=2)

		sigma = 1.8224*sigma/100. # 10^18/10^20
		sigma = round(sigma,4)
		print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t\t{9}'.
			format(i, ret['vs'][i], ret['ve'][i], ret['vs_id'][i], ret['ve_id'][i], ret['nhi_i'][i], ret['ma_nhi'][i], ret['nhi_diff'][i],ret['sources'][i], sigma))

		## Uncomment these lines to see plots ##
		# plt.grid()
		# plt.title('Source ' + data[i]['name'])
		# plt.ylabel('Tb(K)')
		# plt.xlabel('v(km/s)')
		# plt.show()

## Get keys from value - for dictionary, in General #
 #
 # params dict d a dictionary
 # params '' target value to find correspondent key
 #
 # return ...
 # 
 # Author Van Hiep
 ##
def get_keys(d, target):	
    for k, v in d.iteritems():
        path.append(k)
        if isinstance(v, dict):
            get_keys(v, target)
        if v == target:
            result.append(copy(path))
        path.pop()

## Get channel-bin-index from velocity #
#
# params list v_axis list of channel
# params float vel velocity
#
# return index of channel
# 
# Author Van Hiep
##
def get_vel_index(v_axis, vel):
    idx = (np.abs(np.array(v_axis)-vel)).argmin()
    return idx

# =========== End - Define Functions =============== #

# =========== MAIN =============== #
specs = read_specs_data()
path   = []
result = []
get_rms(specs)
# plot_specs(specs, '3C223')
# print specs