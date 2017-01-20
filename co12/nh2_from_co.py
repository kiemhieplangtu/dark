import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import operator

from numpy    import array
from restore  import restore
from plotting import cplot

## Read Wco from Dame to compare #
 #
 # params string fname Filename
 #
 # return dict info of Wco
 # 
 # version 9/2016
 # Author Van Hiep ##
def read_wco_dame(fname = 'result/18src_wco12_comp_dame.dat'):
	cols = ['indx', 'src','tele','l', 'b', 'xplane', 'yplane', 'v1', 'v2', 'vx1', 'vx2', 'wco', 'wco_dame', 'diff']
	fmt  = ['i','s','s','f','f','i','i','f','f','i','i','f','f','s']
	data = restore(fname, 2, cols, fmt)
	return data.read()

## Read Wco and Error #
 #
 # params string fname Filename
 #
 # return dict info of Wco
 # 
 # version 9/2016
 # Author Van Hiep ##
def read_wco_18src(fname = 'result/wco12_18src.txt'):
	cols = ['indx', 'src', 'v1', 'v2', 'wco', 'wco_er']
	fmt  = ['i','s','f','f','f','f']
	data = restore(fname, 2, cols, fmt)
	return data.read()

## Read Wco and NHI and Error from common 16src #
 #
 # params string fname Filename
 #
 # return dict info of Wco
 # 
 # version 9/2016
 # Author Van Hiep ##
def read_wco_16src(fname = 'result/wco12_nhi_16src.txt'):
	cols = ['indx', 'src','file', 'v1', 'v2', 'wco', 'wco_er', 'nhi', 'nhi_er']
	fmt  = ['i','s','s','f','f','f','f','f','f']
	data = restore(fname, 2, cols, fmt)
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

## Calculate N(H2) from X-factor from 18 src #
 #
 # params 
 # return void
 # 
 # Author Van Hiep ##
def cal_nh2_18src():
	xfact = 2.54e20 # H2 cm-2/ (K km/s)
	xferr = 0.13e20

	info  = read_wco_dame()

	src = info['src']
	wco = info['wco']
	tel = info['tele']
	l   = info['l']
	b   = info['b']
	xp  = info['xplane']
	yp  = info['yplane']
	v1  = info['v1']
	v2  = info['v2']

	vx1 = info['vx1']
	vx2 = info['vx2']

	wco      = info['wco']
	wco_dame = info['wco_dame']
	diff_p   = info['diff']

	## Read Wco and error
	inf = read_wco_18src()
	wco = inf['wco']
	wer = inf['wco_er']


	for i in range(len(src)):
		nh2    = 2.*xfact*wco[i]/1.e21
		nh2_er = nh2 * np.sqrt( (xferr/xfact)**2 + (wer[i]/wco[i])**2 )

		# print wer[i]/wco[i], nh2_er/nh2
		print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}'\
			.format(i, src[i], tel[i], l[i],b[i], xp[i], yp[i],v1[i],v2[i],vx1[i],vx2[i],wco[i], nh2, wco_dame[i], diff_p[i]))

		# print('{0}\t{1}\t{2}\t{3}\t{4}'.format(i, src[i], v1[i], v2[i], wco))

## Calculate N(H2) from X-factor from 16 src #
 # X = [2.54+/-0.13]e20
 #
 # params 
 # return void
 # 
 # Author Van Hiep ##
def cal_nh2_16src():
	xfact = 2.54e20 # H2 cm-2/ (K km/s)
	xferr = 0.13e20

	## Read Wco and error
	inf    = read_wco_16src()
	wco    = inf['wco']
	wer    = inf['wco_er']
	src    = inf['src']
	file   = inf['file']
	v1     = inf['v1']
	v2     = inf['v2']
	nhi    = inf['nhi']
	nhi_er = inf['nhi_er']


	for i in range(len(src)):
		nh2    = 2.*xfact*wco[i]/1.e20
		nh2_er = nh2 * np.sqrt( (xferr/xfact)**2 + (wer[i]/wco[i])**2 )

		nh     = nhi[i] + nh2
		nh_er  = np.sqrt( (nhi_er[i])**2 + (nh2_er)**2 )

		# print nh, nh_er/nh, nhi_er[i]/nhi[i], nh2_er/nh2
		if(i<=12):
			print('{0}\t{1}\t\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}'\
				.format(i, src[i],file[i],v1[i],v2[i],wco[i], wer[i], nhi[i], nhi_er[i], nh2, nh2_er, nh, nh_er))
		else:
			print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}'\
				.format(i, src[i],file[i],v1[i],v2[i],wco[i], wer[i], nhi[i], nhi_er[i], nh2, nh2_er, nh, nh_er))

## Calculate N(H2) from X-factor from 16 src #
 # X = [2.0+/-30%]e20 from Alberto D.Bolatto 2013
 #
 # params 
 # return void
 # 
 # Author Van Hiep ##
def cal_nh2_16src_f2():
	xfact = 2.00e20 # H2 cm-2/ (K km/s)
	xferr = 0.6e20

	## Read Wco and error
	inf    = read_wco_16src()
	wco    = inf['wco']
	wer    = inf['wco_er']
	src    = inf['src']
	file   = inf['file']
	v1     = inf['v1']
	v2     = inf['v2']
	nhi    = inf['nhi']
	nhi_er = inf['nhi_er']


	for i in range(len(src)):
		nh2    = 2.*xfact*wco[i]/1.e20
		nh2_er = nh2 * np.sqrt( (xferr/xfact)**2 + (wer[i]/wco[i])**2 )

		nh     = nhi[i] + nh2
		nh_er  = np.sqrt( (nhi_er[i])**2 + (nh2_er)**2 )

		# print nh, nh_er/nh, nhi_er[i]/nhi[i], nh2_er/nh2
		if(i<=12):
			print('{0}\t{1}\t\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}'\
				.format(i, src[i],file[i],v1[i],v2[i],wco[i], wer[i], nhi[i], nhi_er[i], nh2, nh2_er, nh, nh_er))
		else:
			print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}'\
				.format(i, src[i],file[i],v1[i],v2[i],wco[i], wer[i], nhi[i], nhi_er[i], nh2, nh2_er, nh, nh_er))			

############## MAIN ###########
# cal_nh2_18src()
# cal_nh2_16src()
cal_nh2_16src_f2()
sys.exit()