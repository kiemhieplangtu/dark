import numpy             as np
import matplotlib.pyplot as plt
import operator

from copy import copy

# =========== Define Functions =============== #

# Read data from file into columns #
#
# params string fname Filename
#
# return dict data columns
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

# Read data into structure #
#
# params none
#
# return dict Structure of Specs data
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

# Plot spectra #
#
# params dict data Structure data of Spectra
# params string src Source_name
#
# return void
# 
# Author Van Hiep
##
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

	plt.plot(data[i]['v'][vstart:vend], data[i]['t'][vstart:vend], 'r-')
	plt.grid()
	plt.title('Source ' + src)
	plt.ylabel('Tb(K)')
	plt.xlabel('v(km/s)')
	plt.show()

	if (len(vpixout) != 0) :
		vpixout = [vstart, vend]

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

# Get keys from value - In General #
#
# params ndarray d A N-Dimension array
# params float target Value in array
#
# return list Keys in ndarray of target value
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

# Read velocity_range from file -- k, vstart, vend, source_name #
#
# params dict specs_data Structure data of Spectra
# params string fname Filename
#
# return dict velocity_range, corresponding velocity_index
# 
# Author Van Hiep
##
def read_vrange_from_file(specs_data, fname = "data/vrange4sources.txt"):
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

## Get average of channel_width, vaxis, vstart_index, vend_index #
 #
 # params list vaxis List of Velocity_axis
 # params integer vs_id Velocity_start_index
 # params integer ve_id Velocity_end_index
 #
 # return float average of velocity_channel_width (about 0.16km/s)
 # 
 # Author Van Hiep ##
def av_channel_width(vaxis, vs_id, ve_id):
	res = map(operator.sub, vaxis[vs_id:(ve_id-1)], vaxis[(vs_id+1):ve_id])
	return round(sum(res)/float(len(res)), 2)

## Get HI Column Density under optically thin assumption N(HI) #
 #
 # params dict data Structure data of Spectra
 #
 # return NHI_i
 # Author Van Hiep ##
def get_nhi_thin(data):
	vrange = read_vrange_from_file(data)
	print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}'.format('#', 'Vel_start', 'Vel_end', 'Vel_start_index', 'Vel_end_index', 'N(HI) (1e20)', 'Source'))
	for i in range(0,79):
		av_ch_width = av_channel_width(data[i]['v'], vrange[i]['vs_id'], vrange[i]['ve_id'])
		nhi_i       = sum(data[i]['t'][vrange[i]['vs_id']:vrange[i]['ve_id']])*av_ch_width*0.018224		
		nhi_i       = round(nhi_i, 2)
		print('{0}\t{1}\t\t{2}\t\t{3}\t\t{4}\t\t{5}\t\t{6}'.format(i, vrange[i]['vstart'], vrange[i]['vend'], vrange[i]['ve_id'], vrange[i]['vs_id'], nhi_i, data[i]['name']))



# =========== MAIN =============== #
specs = read_specs_data()
get_nhi_thin(specs)