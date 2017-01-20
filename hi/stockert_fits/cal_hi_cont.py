import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class
import matplotlib.pyplot as plt
import numpy             as np
import healpy            as hp
import pylab             as pl

import operator
from restore             import restore
from astropy.io          import fits

## Read info of 26 no-CO sources #
 # l,b, nhi, and nhi_error
 #
 # params string fname Filename
 # return dict info
 # 
 # version 11/2016
 # Author Van Hiep ##
def read_23oh_src(fname = '../../oh/result/23src_with_oh.txt'):
	cols = ['src','l','b']
	fmt  = ['s',  'f','f']
	data = restore(fname, 2, cols, fmt)
	dat  = data.read()
	return dat

## Get continuum temperature from Stockert Telescope #
 #
 # params string src Source-name
 # return float tb Continuum Temperature
 # 
 # version 11/2016
 # Author Van Hiep ##
def get_cont(src):
	src     = src.lower()
	hdulist = fits.open('data/fits'+src+'_17min.bin')
	hdu     = hdulist[0]
	tb      = np.mean(hdu.data)

	return tb

## Read Tbg of 408MHz from Healpy map ##
 #
 # params string fname File-name
 #
 # return dict ret Info of Vel-range
 #
 # version 09/2016 
 # author Nguyen Van Hiep ##
def read_tbg408_healpy(fname='../../oh/result/bg408_to_compare.txt'):
	cols     = ['src','l','b', 'il', 'ib', 'tbg', 'l-idx','b-idx','tbg1','tbg_hpy']
	fmt      = ['s','f','f','f','f','f','f','f','f','f']
	src      = restore(fname, 2, cols, fmt)
	info     = src.read()

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

	ret = {}
	for i in range(len(src)):
		ret[src[i]] = bgh[i]

	return ret

## Get 21cm continuum temperature values #
#
# params str map_file File of maps
# params int src_num Number of sources
# params dict info Information of sources
#
# return void
# 
# Author Van Hiep
##	
def plot_patches(src_num, info):
	# Define constants #
	deg2rad   = np.pi/180.
	# Define the width of area #
	beam   = 14.4           # Beam = 3.5'
	dbeam  = beam/60./2.     # Beam = 3.5' -> dbeam = beam/60/2 in degree
	offset = dbeam          # degree

	# HI continuum map and resolution #
	cont  = hp.read_map(os.getenv("HOME")+'/hdata/hi/lambda_chipass_healpix_r10.fits', field = 0, h=False)
	nside = hp.get_nside(cont)
	res   = hp.nside2resol(nside, arcmin=False)
	dd    = res/deg2rad/5.0

	# OK - Go #
	tb    = {}
	for i in range(0,src_num):
		if(i ==2): continue
		if(i ==3): continue
		if(i ==6): continue
		if(i ==11): continue
		## Find the values of Continuum temperature #
		tb[i] = []
		
		## Longitude and latitude ##
		l     = info['l'][i]
		b     = info['b'][i]

		# if(i != 14): continue

		# Plot cartview a/o mollview #
		ll = l
		if (l>180):
			ll = ll-360.

		f = 10.
		# if (i == sr):
		proj = hp.cartview(cont, title=info['src'][i]+'('+str(info['l'][i])+','+str(info['b'][i])+')', coord='G', unit='',
				norm=None, xsize=1920, lonra=[ll-0.5,ll+0.5], latra=[b-0.5,b+0.5],
				return_projected_map=True)

		# hp.cartview(cont, title=info['src'][i]+'('+str(info['l'][i])+','+str(info['b'][i])+')', coord='G', unit='',
		# 		norm='hist', xsize=800, lonra=[ll-offset-f*offset,ll+offset+f*offset], latra=[b-offset-f*offset,b+offset+f*offset],
		# 		return_projected_map=True)

		# hp.orthview(cont, title=info['src'][i]+'('+str(info['l'][i])+','+str(info['b'][i])+')', coord='G', unit='',
		# 		norm='hist', xsize=800, return_projected_map=True)

		# hp.mollview(cont, title=info['src'][i]+'('+str(info['l'][i])+','+str(info['b'][i])+')',
		# 	coord='G', unit='', rot=[0,0,0], norm='hist')

		# hp.mollzoom(cont, title=info['src'][i]+'('+str(info['l'][i])+','+str(info['b'][i])+')',
		# 	coord='G', unit='', rot=[0,0,0], norm=None, min=4599., max=4600.)

		print proj
		print proj.shape

		# Cal. #
		theta = (90.0-b)*deg2rad
		phi   = l*deg2rad
		pix   = hp.ang2pix(nside, theta, phi, nest=False)

		if (cont[pix] > -1.0e30) : # Some pixels not defined
			tb[i].append(cont[pix])

		for x in pl.frange(l-offset, l+offset, dd):
			for y in pl.frange(b-offset, b+offset, dd):
				cosb = np.cos(b*deg2rad)
				cosy = np.cos(y*deg2rad)
				if ( (((x-l)**2 + (y-b)**2) <= offset**2) ):
					theta = (90.0 - y)*deg2rad
					phi   = x*deg2rad
					pix   = hp.ang2pix(nside, theta, phi, nest=False)
					# hp.projtext(x, y, '.'+str(pix)+str(cont[pix]), lonlat=True, coord='G')
					# hp.projtext(x, y, '.'+str(cont[pix]), lonlat=True, coord='G')
					hp.projtext(x, y, '.', lonlat=True, coord='G')

		plt.show()

##================= MAIN ========================##
## Read Tbg408 from Haslam, Read infor of 23 OH src ##
bg408 = read_tbg408_healpy()
info  = read_23oh_src('../../oh/result/23src_with_oh.txt')

# Define constants #
deg2rad = np.pi/180.

# Define the width of area #
beam    = 14.4           # Beam = 3.5'
dbeam   = beam/60./2.    # Beam = 3.5' -> dbeam = beam/60/2 in degree
offset  = dbeam          # degree

# HI continuum map and resolution #
cont  = hp.read_map(os.getenv("HOME")+'/hdata/hi/lambda_chipass_healpix_r10.fits', field = 0, h=False)
nside = hp.get_nside(cont)
res   = hp.nside2resol(nside, arcmin=False)
dd    = res/deg2rad/5.0

num_of_src = len(info['src'])
plot_patches(num_of_src, info)

sys.exit()

# Product Name
# CHIPASS 1.4 GHz Continuum Map
# Coord. System
# Galactic
# Projection Type
# HEALPix, nested, res 10 (Nside=1024)
# Resolution
# 14.4 arcmin
# Original Data Source ATNF

tb = {}
for i in range(len(info['src'])):
	# if(i != 14): continue

	## Find the values of Tau353 and Err_tau353 in small area #
	tb[i] = []

	## Longitude and latitude ##
	l = info['l'][i]
	b = info['b'][i]
	s = info['src'][i]

	# # Plot cartview a/o mollview #
	# ll = l
	# if (l>180):
	# 	ll = ll-360.

	# hp.cartview(tau_map, title=info['src'][i]+'('+str(info['l'][i])+','+str(info['b'][i])+') - '+map_file, coord='G', unit='',
	# 		norm='hist', xsize=800, lonra=[ll-offset-0.1*offset,ll+offset+0.1*offset], latra=[b-offset-0.1*offset,b+offset+0.1*offset],
	# 		return_projected_map=True)
	# hp.cartview(err_map, title='Err', coord='G', unit='',
	# 		norm='hist', xsize=800, lonra=[ll-offset-0.1*offset,ll+offset+0.1*offset], latra=[b-offset-0.1*offset,b+offset+0.1*offset],
	# 		return_projected_map=True)
	# # End Plot cartview a/o mollview ##

	# Cal. #
	theta = (90.0-b)*deg2rad
	phi   = l*deg2rad
	pix   = hp.ang2pix(nside, theta, phi, nest=False)

	# if ( tau_map[pix] > -1.0e30 ): # Checked Invalid error & Some pixels not defined
	tb[i].append(cont[pix])

	for x in pl.frange(l-offset, l+offset, dd):
		for y in pl.frange(b-offset, b+offset, dd):
			cosb = np.cos(b*deg2rad)
			cosy = np.cos(y*deg2rad)

			if ( ((x-l)**2 + (y-b)**2) <= offset**2 ):
				# hp.projplot(x, y, 'kx', lonlat=True, coord='G')
				theta = (90.0 - y)*deg2rad
				phi   = x*deg2rad
				pix   = hp.ang2pix(nside, theta, phi, nest=False)

				# if (tau_map[pix] > -1.0e30): # Checked Invalid error & Some pixels not defined
				tb[i].append(cont[pix])  # in mK (mili Kelvin)

	tb[i] = list(set(tb[i]))
	tb[i] = np.asarray(tb[i])
	tb[i] = np.mean( tb[i] )
	if(tb[i]<0.):
		tb[i] = get_cont(s)

	tb[i]   = tb[i]/1000.
	tc1     = 2.8+tb[i]*(1420./1665.402)**2.8
	tc2     = 2.8+tb[i]*(1420./1667.359)**2.8
	tbg1665 = 2.8+bg408[s]*(408./1665.402)**2.8 # Tbg from 408MHz
	print( '{}    {}    {}    {}    {}    {}    {}    {}'\
		.format(i, s, l, b, tb[i], tc1, tc2, tbg1665) )