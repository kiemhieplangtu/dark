import sys
sys.path.insert(0, r'/home/vnguyen/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import healpy            as hp
import pylab             as pl
import operator

from restore              import restore

## Read Infor of 26 sources with Tex values ##
 #
 # params string fname File-name
 #
 # return dict ret Info of Vel-range
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def read_src_lb(fname='result/tbg_408_to_compare.txt'):
	cols     = ['src','l','b', 'il', 'ib', 'tbg', 'l-idx','b-idx','tbg1','tbg_hpy']
	fmt      = ['s','f','f','f','f','f','f','f','f','f']
	src      = restore(fname, 2, cols, fmt)
	info     = src.read()

	return info

## Calculate Continuum Temperature of 408MHz from Haslam ##
 # PLot Galatic map
 #
 # params string fname File-name
 #
 # return dict ret Info of Vel-range
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def plot_haslam_408mhz_map():
	# map_file = 'data/haslam408_dsds_Remazeilles2014.fits'
	# lambda_haslam408_nofilt.fits
	# lambda_haslam408_dsds.fits
	# lambda_zea_haslam408_dsds.fits
	# haslam408_dsds_Remazeilles2014.fits
	# haslam408_ds_Remazeilles2014.fits
	# lambda_chipass_healpix_r10.fits

	# Define constants #
	map_file = 'data/haslam408_dsds_Remazeilles2014.fits'
	deg2rad  = np.pi/180.
	factor   = (408./1665.)**2.8

	# Define the width of area #
	beam   = 3.5            # Beam = 3.5'
	dbeam  = beam/120.0     # Beam = 3.5' -> dbeam = beam/60/2 in degree
	offset = dbeam          # degree
	edge   = dbeam/40. # To plot a rectangle about the source

	info = read_src_lb()
	src = info['src']
	gl  = info['l']
	gb  = info['b']
	il  = info['il']
	ib  = info['ib']
	tbg = info['tbg']
	lid = info['l-idx']
	bid = info['b-idx']
	bg1 = info['tbg1']
	bgh = info['tbg_hpy']

	tb408 = hp.read_map(map_file,field=0, nest=False, hdu=1, h=False, verbose=False)
	nside = hp.get_nside(tb408)
	res   = hp.nside2resol(nside, arcmin=False)
	dd    = res/deg2rad/10.0
	hp.cartview(tb408, title='Continuum background at 408 MHz from Haslam et al. Survey', coord='G', unit='K',\
		norm='hist', xsize=800) #min=19.9,max=20.6

	for i in range(len(src)):
		sl = gl[i]
		sb = gb[i]

		if (sl > 180.) :
			sl = sl - 360.

		## Plot cartview a/o mollview #
		# ll = sl
		# if (sl>180):
		# 	ll = ll-360.

		# hp.cartview(tb408, title=src[i]+': ' + str(sl) + ',' + str(sb), coord='G', unit='',
		# 		norm='hist', xsize=800, lonra=[ll-offset-0.1*offset,ll+offset+0.1*offset], latra=[sb-offset-0.1*offset,sb+offset+0.1*offset],
		# 		return_projected_map=True)
		## End Plot cartview a/o mollview ##

		theta = (90.0 - sb)*deg2rad
		phi   = sl*deg2rad
		pix   = hp.ang2pix(nside, theta, phi, nest=False)

		tbg_i = []
		if (tb408[pix] > 0.) : # Some pixels not defined
			tbg_i.append(tb408[pix])

		for x in pl.frange(sl-offset, sl+offset, dd):
				for y in pl.frange(sb-offset, sb+offset, dd):
					cosb = np.cos(sb*deg2rad)
					cosy = np.cos(y*deg2rad)

					if ( ((x-sl)**2 + (y-sb)**2) <= offset**2 ):
						# hp.projtext(x, y, '.', lonlat=True, coord='G')
						theta = (90.0 - y)*deg2rad
						phi   = x*deg2rad
						pix   = hp.ang2pix(nside, theta, phi, nest=False)

						if (tb408[pix] > 0.) :
							tbg_i.append(tb408[pix])

		# plt.show()
		# tbg_i = list(set(tbg_i))
		tbg_i = np.asarray(tbg_i, dtype = np.float32)	
		val   = np.mean(tbg_i)

		if(i<18):
			print '{0}\t\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}'\
			.format(src[i], gl[i], gb[i], il[i], ib[i], tbg[i], lid[i], bid[i], bg1[i], val)
		else:
			print '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}'\
			.format(src[i], gl[i], gb[i], il[i], ib[i], tbg[i], lid[i], bid[i], bg1[i], val)

		b_edge = 300.*edge
		l_edge = b_edge/np.cos(sb*deg2rad)

		equateur_lon = [sl-l_edge, sl+l_edge, sl+l_edge, sl-l_edge, sl-l_edge]
		equateur_lat = [sb+b_edge, sb+b_edge, sb-b_edge, sb-b_edge, sb+b_edge]
		hp.projplot(equateur_lon, equateur_lat, lonlat=True, coord='G')
		hp.projplot(sl, sb, 'kx', lonlat=True, coord='G')

		# hp.projtext(sl, sb, src[i] +' ('+ str("{0:.4e}".format(val))+')', lonlat=True, coord='G') 
		hp.projtext(sl, sb, src[i], lonlat=True, coord='G') 

	plt.show()

## Calculate Continuum Temperature of 408MHz from Haslam ##
 # Plot patches around each line-of-sight
 #
 # params string fname File-name
 #
 # return dict ret Info of Vel-range
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def plot_haslam_408mhz_patch():
	# map_file = 'data/haslam408_dsds_Remazeilles2014.fits'
	# lambda_haslam408_nofilt.fits
	# lambda_haslam408_dsds.fits
	# lambda_zea_haslam408_dsds.fits
	# haslam408_dsds_Remazeilles2014.fits
	# haslam408_ds_Remazeilles2014.fits
	# lambda_chipass_healpix_r10.fits

	# Define constants #
	map_file = 'data/haslam408_dsds_Remazeilles2014.fits'
	deg2rad  = np.pi/180.
	factor   = (408./1665.)**2.8

	# Define the width of area #
	beam   = 3.5            # Beam = 3.5'
	dbeam  = beam/120.0     # Beam = 3.5' -> dbeam = beam/60/2 in degree
	offset = dbeam          # degree
	edge   = dbeam/40. # To plot a rectangle about the source

	info = read_src_lb()
	src = info['src']
	gl  = info['l']
	gb  = info['b']
	il  = info['il']
	ib  = info['ib']
	tbg = info['tbg']
	lid = info['l-idx']
	bid = info['b-idx']
	bg1 = info['tbg1']
	bgh = info['tbg_hpy']

	tb408 = hp.read_map(map_file,field=0, nest=False, hdu=1, h=False, verbose=False)
	nside = hp.get_nside(tb408)
	res   = hp.nside2resol(nside, arcmin=False)
	dd    = res/deg2rad/10.0
	# hp.cartview(tb408, title='Continuum background at 408 MHz from Haslam et al. Survey', coord='G', unit='K',norm='hist', xsize=800)

	for i in range(len(src)):
		if(i != 8):
			continue
		sl = gl[i]
		sb = gb[i]

		if (sl > 180.) :
			sl = sl - 360.

		# Plot cartview a/o mollview #
		ll = sl
		if (sl>180):
			ll = ll-360.

		margin = 100.
		hp.cartview(tb408, title=src[i]+' (l=' + str(round(sl,2)) + ', b=' + str(round(sb,2)) + ')', coord='G', unit='K',
				norm='hist', xsize=800, lonra=[ll-offset-margin*offset,ll+offset+margin*offset], latra=[sb-offset-margin*offset,sb+offset+margin*offset],
				return_projected_map=True, min=11.2, max=4.04e3)
		# End Plot cartview a/o mollview ##

		theta = (90.0 - sb)*deg2rad
		phi   = sl*deg2rad
		pix   = hp.ang2pix(nside, theta, phi, nest=False)

		tbg_i = []
		if (tb408[pix] > 0.) : # Some pixels not defined
			tbg_i.append(tb408[pix])

		for x in pl.frange(sl-offset, sl+offset, dd):
				for y in pl.frange(sb-offset, sb+offset, dd):
					cosb = np.cos(sb*deg2rad)
					cosy = np.cos(y*deg2rad)

					if ( ((x-sl)**2 + (y-sb)**2) <= offset**2 ):
						# hp.projtext(x, y, '*', lonlat=True, coord='G')
						hp.projplot(x, y, 'kx', lonlat=True, coord='G')
						theta = (90.0 - y)*deg2rad
						phi   = x*deg2rad
						pix   = hp.ang2pix(nside, theta, phi, nest=False)

						if (tb408[pix] > 0.) :
							tbg_i.append(tb408[pix])

		plt.show()
		# tbg_i = list(set(tbg_i))
		tbg_i = np.asarray(tbg_i, dtype = np.float32)	
		val   = np.mean(tbg_i)

		if(i<18):
			print '{0}\t\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}'\
			.format(src[i], gl[i], gb[i], il[i], ib[i], tbg[i], lid[i], bid[i], bg1[i], val)
		else:
			print '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}'\
			.format(src[i], gl[i], gb[i], il[i], ib[i], tbg[i], lid[i], bid[i], bg1[i], val)


## Calculate Continuum Temperature of 408MHz from Haslam ##
## ======================= MAIN ========================== ##
# plot_haslam_408mhz_map()
plot_haslam_408mhz_patch()