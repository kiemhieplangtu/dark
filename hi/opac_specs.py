import numpy             as np
import matplotlib.pyplot as plt
import operator

from copy import copy

# =========== Define Functions =============== #

# Read data from file into columns #
def read_data_from_file(file = "nhi_opac_specs.txt"):
	ret  = {}

	name = []
	v    = []
	t    = []
	opac = []

	f    = open (file,"r")
	f.readline() # comment line
	f.readline() # comment line
	f.readline() # comment line

	#read line into array 
	for line in f:
	    line    = line.strip()
	    columns = line.split()

	    name.append(columns[0])
	    v.append(float(columns[1]))
	    t.append(float(columns[2]))
	    opac.append((-1.0)*np.log(float(columns[3])))

	f.close()

	ret['name'] = name
	ret['v']    = v
	ret['t']    = t
	ret['opac'] = opac

	return ret

# Read data into structure #
def read_specs_data():
	a     = read_data_from_file()
	specs = {}

	name  = a['name']
	v     = a['v']
	t     = 0.5*np.array(a['t'])
	opac  = a['opac']

	for i in range(0,79):
		specs[i] = {}

		start = i*2048
		end   = start + 2047
		
		specs[i]['name'] = name[end]
		specs[i]['v']    = v[start:end+1]
		specs[i]['t']    = t[start:end+1]
		specs[i]['opac'] = opac[start:end+1]

	return specs


# Plot opacity spectra #
def plot_opac(data, src, vrange = [], vpixout = [1,2]):
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

	plt.plot(data[i]['v'][vstart:vend], data[i]['opac'][vstart:vend], 'r-')
	plt.grid()
	plt.title('Source ' + src)
	plt.ylabel('1-exp(-tau)')
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

# Read velocity_range from file -- k, vstart, vend, source_name #
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
# Get average of channel_width, vaxis, vstart_index, vend_index #
def av_channel_width(vaxis, vs_id, ve_id):
	res = map(operator.sub, vaxis[vs_id:(ve_id-1)], vaxis[(vs_id+1):ve_id])
	return round(sum(res)/float(len(res)), 2)

# Get N(HI) #
def get_nhi(data):
	vrange = read_vrange_from_file(data)
	print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}'.format('#', 'Vel_start', 'Vel_end', 'Vel_start_index', 'Vel_end_index', 'N(HI) (1e20)', 'Source'))
	for i in range(0,79):
		av_ch_width = av_channel_width(data[i]['v'], vrange[i]['vs_id'], vrange[i]['ve_id'])
		nhi_i       = sum(data[i]['t'][vrange[i]['vs_id']:vrange[i]['ve_id']])*av_ch_width*0.018
		nhi_i       = round(nhi_i, 2)
		print('{0}\t{1}\t\t{2}\t\t{3}\t\t{4}\t\t{5}\t\t{6}'.format(i, vrange[i]['vstart'], vrange[i]['vend'], vrange[i]['ve_id'], vrange[i]['vs_id'], nhi_i, data[i]['name']))


# Get N(HI) #
def get_tau_integration(data):
	vrange = read_vrange_from_file(data)
	print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}'.format('#', 'Vel_start', 'Vel_end', 'Vel_start_index', 'Vel_end_index', 'N(HI) (1e20)', 'Source'))
	for i in range(0,79):
		av_ch_width = av_channel_width(data[i]['v'], vrange[i]['vs_id'], vrange[i]['ve_id'])
		tau_intg_i       = sum(data[i]['opac'][vrange[i]['vs_id']:vrange[i]['ve_id']])*av_ch_width
		tau_intg_i       = round(tau_intg_i, 2)
		print('{0}\t{1}\t\t{2}\t\t{3}\t\t{4}\t\t{5}\t\t{6}'.format(i, vrange[i]['vstart'], vrange[i]['vend'], vrange[i]['ve_id'], vrange[i]['vs_id'], tau_intg_i, data[i]['name']))

# =========== End - Define Functions =============== #



# =========== MAIN =============== #
specs = read_specs_data()

# Plot spectra #
k = 0
while k < 79:
	src = raw_input('Enter source: ')
	path   = []
	result = []

	
	#plot_specs(specs, src)
	plot_opac(specs, src)
	k = k +1

#get_tau_integration(specs)