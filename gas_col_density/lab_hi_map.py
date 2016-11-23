import os, sys
import matplotlib.pyplot as plt
import numpy             as np
import healpy            as hp
import operator
import matplotlib.cm as cm

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
# map_file = 'data/LAB_Healpix/LAB_cut+440.fits' #LAB_cut+30.fits
map_file = 'data/LAB_Healpix/LAB_IV.fits' #
# map_file = 'data/LAB_Healpix/LAB_fullvel.fits' #
whi = hp.read_map(map_file, field = 0, h=False)
whi  = whi*1.8224e18
# whii = hp.read_map('data/LAB_Healpix/LAB_IV.fits', field = 0, h=False)
# whii  = whii*1.8224e18
title = map_file

# Open a file
path = "/home/vnguyen/dark/gas_col_density/data/LAB_Healpix/"
dirs = os.listdir( path )

# # This would print all the files and directories
# namelist = []
# vlist    = []
# for names in dirs:
#     if (names.endswith(".fits") and names.startswith('LAB_cut')) :
#         namelist.append(names)
#         names = names.replace("LAB_cut","")
#         names = names.replace(".fits","")
#         vlist.append(float(names))

# # print len(namelist)

# whi = 0.
# title = ''
# for i in range(len(namelist)):
# 	# print i, namelist[i], vlist[i]
# 	# if (np.abs(vlist[i]) < 35.):
# 	# if ( (np.abs(vlist[i]) < 75.) and (np.abs(vlist[i]) > 30.) ):
# 	if ( (vlist[i] < -30.) and (vlist[i] > -100.) ):
# 		map_file = 'data/LAB_Healpix/' + namelist[i]
# 		title = title + str(vlist[i])+', '
# 		whi = whi + hp.read_map(map_file, field = 0, h=False)

# whi  = whi*1.8224e18
# wmax = np.max(whi)
# wmin = np.min(whi)
# dw   = (wmax - wmin)/100.
# wmin1 = wmin+dw


# whii = 0.
# title = ''
# for i in range(len(namelist)):
# 	if ( (np.abs(vlist[i]) < 105.) and (np.abs(vlist[i]) > 30.) ):
# 		map_file = 'data/LAB_Healpix/' + namelist[i]
# 		title = title + str(vlist[i])+', '
# 		whii = whii + hp.read_map(map_file, field = 0, h=False)

# whii  = whii*1.8224e18

# npix = 0
# for i in range(len(whi)):
# 	# if( (whi[i] > 2.0e20) or (whi[i] > 1.0e19) ):
# 	# if( whii[i] > 1.0e19 ):
# 	if( whi[i] > wmin ):
# 		whi[i] = hp.UNSEEN
# 	else:
# 		npix = npix + 1
# 		whi[i] = 1

# print 100.*npix/len(whi)

# whi  = np.log10(whi)
whi = whi/1e20
cmap = plt.cm.get_cmap('jet')
# hp.mollview(whi, title=title, coord='G', unit='', rot=[0,0,0], norm=None,  xsize=800)
# hp.mollview(whi, title=map_file, coord='G', unit='', rot=[0,0,0], norm=None, min=0.,max=2.e20, xsize=800, cmap = cmap)
hp.mollview(whi, title=title, coord='G', unit='', rot=[0,0,0], norm=None, min=0., max=2.)

# hp.cartview(whi, title=title, coord='G', unit='', rot=[0.,0., 0.], nest=False, flip=None,
# 				norm=None, return_projected_map=True)
# hp.graticule()
plt.show()
# hp.write_map("my_map.fits", whi)
sys.exit()

# map_file = 'data/lambda_combined_nh.fits'
# quantity = 'tau353'
# map_unit = ''

# # map_file = 'data/HFI_SkyMap_353_2048_R2.02_full.fits'
# # quantity = 'w353'
# # map_unit = 'K_CMB'


# info  = read_info_no_co()
# beam  = 3.5
# dbeam = beam/120.0 # Beam = 3.5' -> dbeam = beam/60/2
# dbeam = 0.5
# edge  = dbeam/40. # To plot a rectangle about the source

# # Define constants #
# nside     = 2048
# deg2rad   = np.pi/180.

whi = hp.read_map(map_file, field = 0, h=False)

# print whi
# whi = whi*1.8224e18
# whi = np.log10(whi)

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

hp.mollview(whi, title=map_file, coord='G', unit='', rot=[0,0,0], norm='hist', xsize=800)
# hp.mollview(whi, title=map_file, coord='G', unit='', rot=[0,0,0], norm='hist', min=1.e9,max=2.e22, xsize=800)
# hp.mollview(whi, title=map_file, coord='G', unit='', rot=[0,0,0], norm='hist', min=19.7,max=22., xsize=800)
# hp.mollview(fullvel, title=map_file, coord='G', unit='', rot=[0,0,0], norm='hist', xsize=800)

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