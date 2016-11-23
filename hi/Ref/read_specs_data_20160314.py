import numpy             as np
import matplotlib.pyplot as plt

from copy import copy

# =========== Define Functions =============== #

# Read data from file into columns #
def read_data_from_file():
	ret  = {}

	name = []
	v    = []
	t    = []

	f    = open ("specs.txt","r")
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

# Get the index of a given velocity #
def get_vel_index(v_axis, vel):
    idx = (np.abs(np.array(v_axis)-vel)).argmin()
    return idx

# Get keys from value - In General #
def get_keys(d, target):	
    for k, v in d.iteritems():
        path.append(k)
        if isinstance(v, dict):
            get_keys(v, target)
        if v == target:
            result.append(copy(path))
        path.pop()

#def get_nhi(data, vrange_arr = []):
	# vrange_arr = np.array([3,3])
	# if (len(vrange_arr) == 0):
	# 	vrange_arr    = np.empty((2,79,))
	# 	vrange_arr[:] = np.NAN
	# elif (len(vrange_arr[0]) != 79):
	# 	print 'Wrong size vrange_arr'

	# print vrange_arr




# =========== End - Define Functions =============== #

# =========== MAIN =============== #
k = 0
specs  = read_specs_data()
while k < 79:
	src = raw_input('Enter source: ')
	path   = []
	result = []

	
	plot_specs(specs, src)
	k = k +1
#get_nhi(specs)
