import numpy             as np
import matplotlib.pyplot as plt
import operator

from copy import copy

# =========== Define Functions =============== #

## Read Spectra data from file ##
 #
 # params string file Filename
 # return dict ret Name, vel, Temperature
 #
 # version 03/2016 
 # author Nguyen Van Hiep
 # DON'T CHANGE ##
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

## Read Spectra data ##
 #
 # params NONE
 # return dict specs Name, vel, Temperature
 #
 # version 03/2016 
 # author Nguyen Van Hiep ##
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

## Plot spectra ##
 #
 # params dict data Spectra of sources
 # params string src Source-name
 # params list vrange Velocity-range
 # params list vpixout Do Nothing :) ^^
 # return None, just plot a spectrum
 #
 # version 03/2016 
 # author Nguyen Van Hiep ##
def plot_specs(data, src, vrange = [], vpixout = [1,2]):
	# Get key from Source_name #
	get_keys(data, src)
	i = result[0][0]

	# vstart = get_vel_index(data[i]['v'], 100)
	# print vstart, data[i]['v'][410]

	if (len(vrange) == 2) :
		vstart = get_vel_index(data[i]['v'], vrange[0])
		vend   = get_vel_index(data[i]['v'], vrange[1])
	else :
		vstart = 0
		vend   = 2047

	plt.plot(data[i]['v'][vstart:vend], data[i]['t'][vstart:vend], 'r-', lw=3)
	plt.grid()
	plt.title(src, fontsize = 35)
	plt.ylabel('$T_{b} (K)$', fontsize = 35)
	plt.xlabel('VLSR (km/s)', fontsize = 35)
	plt.tick_params(axis='x', labelsize=20)
	plt.tick_params(axis='y', labelsize=20)
	# plt.text(55, 5,string,color='b',fontsize=14)
	# plt.legend(loc='upper right', fontsize = 18)
	# plt.axvline(x=60., lw=4)
	# plt.axvline(x=-90., lw=4)
	plt.show()

	if (len(vpixout) != 0) :
		vpixout = [vstart, vend]

## Get the index of a given velocity ##
 #
 # params list v_axis Velocity, vlsr
 # params float vel: (One) Velocity value
 # return int idx Index of velocity in List
 #
 # version 03/2016 
 # author Nguyen Van Hiep ##
def get_vel_index(v_axis, vel):
    idx = (np.abs(np.array(v_axis)-vel)).argmin()
    return idx

## Get keys from value - In General ##
 #
 # params dict d A dictionary
 # params ? target Value
 # return ? key Key of value
 #
 # version 03/2016 
 # author Nguyen Van Hiep ##
def get_keys(d, target):	
    for k, v in d.iteritems():
        path.append(k)
        if isinstance(v, dict):
            get_keys(v, target)
        if v == target:
            result.append(copy(path))
        path.pop()

## Read velocity_range from file -- k, vstart, vend, source_name ##
 #
 # params dict specs_data Spectra-data
 # params string fname Filename
 # return dict ret Info on vel-range
 #
 # version 03/2016 
 # author Nguyen Van Hiep ##
def read_vrange_from_file(specs_data, fname = "vrange4sources.txt"):
	ret  = {}

	idx    = []
	vstart = []
	vend   = []
	src    = []

	f = open (fname,"r")
	for line in f:
	    line    = line.strip()
	    columns = line.split()

	    idx.append(int(columns[0]))
	    vstart.append(float(columns[1]))
	    vend.append(float(columns[2]))
	    src.append(columns[3])

	f.close()

	for i in range(len(idx)):
		ret[i] = {}
		
		ret[i]['idx']    = idx[i]
		ret[i]['vstart'] = vstart[i]
		ret[i]['vend']   = vend[i]
		ret[i]['src']    = src[i]

		ret[i]['ve_id']  = get_vel_index(specs_data[i]['v'], vstart[i])
		ret[i]['vs_id']  = get_vel_index(specs_data[i]['v'], vend[i])

	return ret

## Get Channel-Width ##
 #
 # params list vaxis Velocity Axis
 # params int vs_id Starting-ID
 # params int ve_id End-ID
 # return float ret Channel-Width
 #
 # version 03/2016 
 # author Nguyen Van Hiep ##
def av_channel_width(vaxis, vs_id, ve_id):
	res = map(operator.sub, vaxis[vs_id:(ve_id-1)], vaxis[(vs_id+1):ve_id])
	return round(sum(res)/float(len(res)), 2)

##  Cal. N(HI) in opt.thin assumption ##
 #
 # params list vaxis Velocity Axis
 # return float ret N(HI)
 #
 # version 03/2016 
 # author Nguyen Van Hiep ##
def get_nhi(data):
	vrange = read_vrange_from_file(data)
	print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}'.format('#', 'Vel_start', 'Vel_end', 'Vel_start_index', 'Vel_end_index', 'N(HI) (1e20)', 'Source'))
	for i in range(0,79):
		av_ch_width = av_channel_width(data[i]['v'], vrange[i]['vs_id'], vrange[i]['ve_id'])
		nhi_i       = sum(data[i]['t'][vrange[i]['vs_id']:vrange[i]['ve_id']])*av_ch_width*0.018
		nhi_i       = round(nhi_i, 2)
		print('{0}\t{1}\t\t{2}\t\t{3}\t\t{4}\t\t{5}\t\t{6}'.format(i, vrange[i]['vstart'], vrange[i]['vend'], vrange[i]['ve_id'], vrange[i]['vs_id'], nhi_i, data[i]['name']))

## =========== MAIN =============== ##
specs = read_specs_data()

# Plot spectra #
k = 0
while k < 79:
	src    = raw_input('Enter source: ')
	src    = src.upper()
	path   = []
	result = []

	plot_specs(specs, src)
	k = k +1