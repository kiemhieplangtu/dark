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

## Read info of 26 no-CO sources #
 # l,b, nhi, and nhi_error
 #
 # params string fname Filename
 # return dict infocd 
 # 
 # version 1/2017
 # Author Van Hiep ##
def read_info_no_co(fname = '../../co12/result/26src_no_co_with_sponge.dat'):
	cols = ['idx','src','l','b','ra_icrs','de_icrs','ra_j','de_j', 'oh', 'nhi','nhi_er','thin','thin_er', 'cnm','cnm_er','wnm','wnm_er']
	fmt  = ['i',  's',  'f','f', 'f',    'f',       's',    's',    'i', 'f',   'f',     'f',    'f'    , 'f',   'f',     'f',  'f'    ]
	data = restore(fname, 2, cols, fmt)
	dat  = data.read()
	return dat

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

tau  = hp.read_map(map_file, field = 0)
r    = hp.read_map(map_file, field = 3)
t    = hp.read_map(map_file, field = 4)
beta = hp.read_map(map_file, field = 6)

# hp.cartview(dmap, title=map_file, coord='G', unit=map_unit, rot=[90.,-80., 0.], nest=False, flip=None,
# 				norm='log', xsize=800, lonra=[-6.,6.], latra=[-6.,6.], min=0.3e-6,max=1.7e-6, hold=True, cbar=True,
# 				return_projected_map=True)

# hp.cartview(dmap, title=map_file, coord='C', unit=map_unit, rot=[90.,-80., -30.],
# 				norm='log', xsize=800, lonra=[-6.,6.], latra=[-6.,6.], min=2.7e-8,max=8.6e-8,
# 				return_projected_map=True)
# hp.orthview(map=dmap, fig=None, rot=[90.,-80.,0.], coord='G', unit='', xsize=800, 
# 	half_sky=False, title='Orthographic view', nest=False, min=None, max=None, 
# 	flip='astro', remove_dip=False, remove_mono=False, gal_cut=0, format='%g', format2='%g', 
# 	cbar=True, cmap=None, notext=False, norm='hist', hold=False, margins=None, sub=None, return_projected_map=False)

# hp.mollview(dmap, title=map_file,
# 			coord='G', unit='', rot=[0,0,0], norm='hist', min=7e-10,max=0.025, xsize=800)

# hp.mollview(np.log10(dmap), title=map_file,
# 			coord='G', unit='', rot=[0,0,0], norm='hist', min=-6,max=-3, xsize=800)

# hp.mollview(t, title=map_file,
# 			coord='G', unit='', rot=[0,0,0], norm='hist', min=15,max=27, xsize=800)

# hp.mollview(beta, title=map_file,
# 			coord='G', unit='', rot=[0,0,0], norm='hist', min=1,max=2.2, xsize=800)

hp.mollview(np.log10(r), title=map_file,
			coord='G', unit='', rot=[0,0,0], norm='hist', min=-8,max=-5, xsize=800)

# for i in range(0,26):
	
# 	if (i != 16):
# 		continue

# 	# Define constants #
# 	deg2rad   = np.pi/180.

# 	sl    = info['l'][i]
# 	sb    = info['b'][i]

# 	theta = (90.0 - sb)*deg2rad
# 	phi   = sl*deg2rad
# 	pix   = hp.ang2pix(nside, theta, phi, nest=False)

# 	val   = test_map[pix]

# 	if (sl > 180.) :
# 		sl = sl - 360.

# 	long1 = sl - dbeam
# 	long2 = sl + dbeam

# 	lat1  = sb - dbeam
# 	lat2  = sb + dbeam

	
# 	tmap  = hp.cartview(test_map, title=info['src'][i]+'('+str(info['l'][i])+','+str(info['b'][i])+') - '+map_file, coord='G', unit=map_unit, rot=[sl,sb],
# 				norm='hist', xsize=800, lonra=[-dbeam,dbeam], latra=[-dbeam,dbeam],
# 				return_projected_map=True)

# 	# tmap  = hp.cartview(test_map, title=info['src'][i]+'('+str(info['l'][i])+','+str(info['b'][i])+') - '+map_file, coord='G', unit=map_unit,
# 	# 			norm='hist', xsize=800, lonra=[sl-dbeam,sl+dbeam], latra=[sb-dbeam,sb+dbeam],
# 	# 			return_projected_map=True)

# 	# gmap = hp.gnomview(test_map, rot=[sl,sb], coord=None, unit='K', title=info['src'][i], nest=False, min=1e-7, max=1e-3,
# 	# 	return_projected_map=True)

# 	#t353 = np.average(tmap)

# 	l_edge = edge/np.cos(sb*deg2rad)

# 	equateur_lon = [sl-l_edge, sl+l_edge, sl+l_edge, sl-l_edge, sl-l_edge]
# 	equateur_lat = [sb+edge, sb+edge, sb-edge, sb-edge, sb+edge]
# 	hp.projplot(equateur_lon, equateur_lat, lonlat=True, coord='G')

# 	hp.projtext(sl, sb, str("{0:.4e}".format(val)), lonlat=True, coord='G') 

# 	#plt.savefig('s'+str(i) + '_'+ info['src'][i]+'('+str(info['l'][i])+','+str(info['b'][i])+')_' + quantity + '_1o_wdth_accurate_lb.png')

hp.graticule()
plt.show()