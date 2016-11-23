import matplotlib.pyplot as plt
import numpy             as np
import healpy            as hp
import operator

# Read info of 79 sources #
#
# params string fname Filename
#
# return void
# 
# Author Van Hiep
##
def read_info(fname = '79src_info.txt'):
	ret = {}

	ret['src'] = []
	ret['yn']  = []
	ret['l']  = []
	ret['b']  = []
	ret['ra_icrs']  = []
	ret['de_icrs']  = []
	ret['ra_j']  = []
	ret['de_j']  = []

	file = open (fname,'r')
	file.readline() # comment line
	file.readline() # comment line

	for line in file:
	    line    = line.strip()
	    columns = line.split()

	    ret['src'].append(columns[1])
	    ret['yn'].append(int(columns[2]))
	    ret['l'].append(float(columns[3]))
	    ret['b'].append(float(columns[4]))
	    ret['ra_icrs'].append(float(columns[5]))
	    ret['de_icrs'].append(float(columns[6]))
	    ret['ra_j'].append(str(columns[7]))
	    ret['de_j'].append(str(columns[8]))

	file.close()

	return ret

# Find 26 sources with no CO #
#
# params string fname Filename
#
# return void
# 
# Author Van Hiep
##
def read_info_no_co(fname = '79src_info.txt'):
	info = read_info()
	#print map(operator.sub, info['b'], info['bc'])

	j=0
	for i in range(0,79):
		if (info['yn'][i] == 0) :
			print('{0}\t{1}\t\t{2}\t\t{3}\t\t{4}\t\t{5}\t\t{6}\t\t{6}'.format(j, info['src'][i], info['l'][i], info['b'][i], info['ra_icrs'][i], info['de_icrs'][i], info['ra_j'][i], info['de_j'][i]))
			j = j + 1


# Read info of 79 sources #
#
# params string fname Filename
#
# return void
# 
# Author Van Hiep
##
def read_info_no_co(fname = '26src_no_co_info.dat'):
	ret = {}

	ret['src'] = []
	ret['l']  = []
	ret['b']  = []
	ret['ra_icrs']  = []
	ret['de_icrs']  = []
	ret['ra_j']  = []
	ret['de_j']  = []

	file = open (fname,'r')
	file.readline() # comment line
	file.readline() # comment line
	for line in file:
	    line    = line.strip()
	    columns = line.split()

	    ret['src'].append(columns[1])
	    ret['l'].append(float(columns[2]))
	    ret['b'].append(float(columns[3]))
	    ret['ra_icrs'].append(float(columns[4]))
	    ret['de_icrs'].append(float(columns[5]))
	    ret['ra_j'].append(str(columns[6]))
	    ret['de_j'].append(str(columns[7]))

	file.close()

	return ret	

#================= MAIN ========================#	
map_file = 'data/HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
quantity = 'tau353'
map_unit = ''

# map_file = 'data/HFI_SkyMap_353_2048_R2.02_full.fits'
# quantity = 'w353'
# map_unit = 'K_CMB'


info  = read_info_no_co()
beam  = 3.5
dbeam = beam/120.0 # Beam = 3.5' -> dbeam = beam/60/2
dbeam = 0.5
edge  = dbeam/40. # To plot a rectangle about the source

# Define constants #
nside     = 2048
deg2rad   = np.pi/180.

test_map = hp.read_map(map_file, field = 0)

for i in range(0,26):
	
	if (i != 16):
		continue

	# Define constants #
	deg2rad   = np.pi/180.

	sl    = info['l'][i]
	sb    = info['b'][i]

	theta = (90.0 - sb)*deg2rad
	phi   = sl*deg2rad
	pix   = hp.ang2pix(nside, theta, phi, nest=False)

	val   = test_map[pix]

	if (sl > 180.) :
		sl = sl - 360.

	long1 = sl - dbeam
	long2 = sl + dbeam

	lat1  = sb - dbeam
	lat2  = sb + dbeam

	
	tmap  = hp.cartview(test_map, title=info['src'][i]+'('+str(info['l'][i])+','+str(info['b'][i])+') - '+map_file, coord='G', unit=map_unit, rot=[sl,sb],
				norm='hist', xsize=800, lonra=[-dbeam,dbeam], latra=[-dbeam,dbeam],
				return_projected_map=True)

	# tmap  = hp.cartview(test_map, title=info['src'][i]+'('+str(info['l'][i])+','+str(info['b'][i])+') - '+map_file, coord='G', unit=map_unit,
	# 			norm='hist', xsize=800, lonra=[sl-dbeam,sl+dbeam], latra=[sb-dbeam,sb+dbeam],
	# 			return_projected_map=True)

	# gmap = hp.gnomview(test_map, rot=[sl,sb], coord=None, unit='K', title=info['src'][i], nest=False, min=1e-7, max=1e-3,
	# 	return_projected_map=True)

	#t353 = np.average(tmap)

	l_edge = edge/np.cos(sb*deg2rad)

	equateur_lon = [sl-l_edge, sl+l_edge, sl+l_edge, sl-l_edge, sl-l_edge]
	equateur_lat = [sb+edge, sb+edge, sb-edge, sb-edge, sb+edge]
	hp.projplot(equateur_lon, equateur_lat, lonlat=True, coord='G')

	hp.projtext(sl, sb, str("{0:.4e}".format(val)), lonlat=True, coord='G') 

	#plt.savefig('s'+str(i) + '_'+ info['src'][i]+'('+str(info['l'][i])+','+str(info['b'][i])+')_' + quantity + '_1o_wdth_accurate_lb.png')

	plt.show()