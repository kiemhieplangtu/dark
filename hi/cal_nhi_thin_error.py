import sys
sys.path.insert(0, r'/home/vnguyen/dark/common') # add folder of Class
import numpy             as np
import matplotlib.pyplot as plt
import operator

from copy                import copy
from restore             import restore

# =========== Define Functions =============== #

## Read data from file into columns #
 #
 # params string fname Filename
 # return dict data columns
 # 
 # Author Van Hiep ##
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

## Read data into structure #
 #
 # params none
 # return dict Structure of Specs data
 # 
 # Author Van Hiep ##
def read_specs_data():
	a     = read_data_from_file()
	specs = {}

	name  = a['name']
	v     = a['v']
	t     = 0.5*np.array(a['t'])

	for i in range(len(v)):
		if (name[i] not in specs):
			specs[name[i]]      = {}
			specs[name[i]]['v'] = [v[i]]
			specs[name[i]]['t'] = [t[i]]
		else:
			specs[name[i]]['v'] = specs[name[i]]['v'] + [v[i]]
			specs[name[i]]['t'] = specs[name[i]]['t'] + [t[i]]

	return specs

## Read 78 sources with l,b and ra, dec
 # Note: 3C223 with No N(HI) and N(HI)_uncertainty #
 #
 # params string fname Filename
 # return dict info
 # 
 # version 04/10/2016
 # Author Van Hiep ##
def read_lb(fname = ''):
	cols = ['l', 'b', 'src', 'ra', 'dec']
	fmt  = ['f', 'f',  's',   'f',  'f']
	data = restore(fname, 3, cols, fmt)
	return data.read()

## Read Vrange for RMS
 # Note: 3C223 with No N(HI) and N(HI)_uncertainty #
 #
 # params string fname Filename
 # return dict info
 # 
 # version 12/2016
 # Author Van Hiep ##
def read_vrange4rms(fname = ''):
	cols = [ 'idx', 'v1', 'v2', 'src']
	fmt  = [ 'i',   'f',  'f',  's'  ]
	data = restore(fname, 0, cols, fmt)
	return data.read()

## Get the index of a given velocity #
 #
 # params list v-axis List of Velocity_axis
 # params float vel Velocity
 #
 # return int idx Index of vel in List of velocity_axis
 # 
 # Author Van Hiep ##
def get_vel_index(v_axis, vel):
    idx = (np.abs(np.array(v_axis)-vel)).argmin()
    return idx

## Read velocity_range from file -- k, vstart, vend, source_name #
 #
 # params dict specs_data Structure data of Spectra
 # params string fname Filename
 #
 # return dict info
 # 
 # Author Van Hiep ##
def read_vrange(fname = "data/vrange_78src.txt"):
	cols = [ 'idx', 'v1', 'v2', 'src']
	fmt  = [ 'i',   'f',  'f',  's'  ]
	data = restore(fname, 0, cols, fmt)
	return data.read()


## Get average of channel_width, vaxis, vstart_index, vend_index #
 #
 # params list vaxis List of Velocity_axis
 #
 # return float average of velocity_channel_width (about 0.16km/s) 0.160947265625km/s
 # 
 # Author Van Hiep ##
def av_channel_width(vaxis):
	vi = vaxis[0]
	vf = vaxis[-1]
	dv = (vf-vi)/len(vaxis)
	return np.abs(dv)

## Plot spectra #
 #
 # params dict data Structure data of Spectra
 # params string src Source_name
 #
 # return void
 # 
 # Author Van Hiep ##
def plot_specs(data):
	lb  = read_lb('data/78src_radec_lb.txt')
	sc  = lb['src']
	l   = lb['l']
	b   = lb['b']
	vr  = read_vrange4rms('data/vrange4rms.txt')

	for i in range(len(sc)):
		src = sc[i]
		if(src != '3C234'):
			continue
		v   = data[src]['v']
		t   = data[src]['t']

		plt.plot(v, t, 'r-', lw=2)
		plt.grid()
		plt.title(src, fontsize = 35)
		plt.ylabel('$T_{b} (K)$', fontsize = 35)
		plt.xlabel('VLSR (km/s)', fontsize = 35)
		plt.tick_params(axis='x', labelsize=20)
		plt.tick_params(axis='y', labelsize=20)
		plt.axvline(vr['v1'][i], ymin=-100., ymax = 100., linewidth=2, color='k')
		plt.axvline(vr['v2'][i], ymin=-100., ymax = 100., linewidth=2, color='k')
		plt.show()

## Get 1-sigma STD from Spectra #
 #
 # params dict data Structure data of Spectra
 # return void
 # 
 # Author Van Hiep ##
def get_std(data):
	lb  = read_lb('data/78src_radec_lb.txt')
	sc  = lb['src']
	vr  = read_vrange4rms('data/vrange4rms.txt')

	ret = {}
	for i in range(len(sc)):
		src = sc[i]
		v   = data[src]['v']
		t   = data[src]['t']

		rt = []
		for k in range(len(v)):
			if( (v[k] > vr['v1'][i]) and (v[k] < vr['v2'][i]) ):
				rt.append(t[k])

		rt       = np.asarray(rt)
		rms      = np.std(rt)
		ret[src] = rms
	return ret

## Get NHI Thin #
 #
 # params dict data Structure data of Spectra #
 # return void
 # 
 # Author Van Hiep ##
def get_nhi_thin(data):
	lb  = read_lb('data/78src_radec_lb.txt')
	sc  = lb['src']
	l   = lb['l']
	b   = lb['b']
	vr  = read_vrange(fname = "data/vrange_78src.txt")
	sg  = get_std(specs) ## sg[src]

	ret = {}
	for i in range(len(sc)):
		src = sc[i]
		v   = data[src]['v']
		t   = data[src]['t']

		vi_idx = get_vel_index(v, vr['v2'][i])
		vf_idx = get_vel_index(v, vr['v1'][i])

		dv    = av_channel_width(v)
		nhi_i = sum( t[ vi_idx:vf_idx ] )*dv*0.018224		
		nhi_i = round(nhi_i, 4)

		## Uncertainties
		nchnl = vf_idx-vi_idx+1
		sigma = sg[src]
		err   = np.sqrt(nchnl) * sigma * dv
		err   = round(err,4)

		# print src, nhi_i, err, 100*err/nhi_i, sigma
		xl = round(l[i],4)
		xb = round(b[i],4)
		print('{}    {}\t{}    {}    {}   {}'.format(i, src, xl,xb, nhi_i, err  ))


# =========== MAIN =============== #
specs = read_specs_data()
# plot_specs(specs)
# print get_std(specs)

get_nhi_thin(specs)
