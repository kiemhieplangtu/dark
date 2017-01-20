import os, sys, shutil
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt
import pylab             as pl
import copy

from scipy.io.idl        import readsav
from restore             import restore
from mpfit               import mpfit

## Get name of sources ##
 # params string code Code of source
 # return string src  Name of source
 #
 # version 10/2016
 # Author Van Hiep ##
def get_name(code):
	src = '-'
	if(code == 'SRC0'): src='J0022'
	if(code == 'SRC1'): src='3C018'
	if(code == 'SRC2'): src='3C019'
	if(code == 'SRC3'): src='3C041'
	if(code == 'SRC4'): src='3C48'
	if(code == 'SRC5'): src='4C15.05'
	if(code == 'SRC6'): src='3C78'
	if(code == 'SRC7'): src='4C16.09'
	if(code == 'SRC8'): src='J0407'
	if(code == 'SRC9'): src='3C120'
	if(code == 'SRC10'): src='3C123'
	if(code == 'SRC11'): src='3C132'
	if(code == 'SRC12'): src='3C133'
	if(code == 'SRC13'): src='3C138'
	if(code == 'SRC14'): src='J05344A'
	if(code == 'SRC15'): src='J05344B'
	if(code == 'SRC16'): src='PKS0531'
	if(code == 'SRC17'): src='3C154'
	if(code == 'SRC18'): src='PKS0742'
	if(code == 'SRC19'): src='3C225A'
	if(code == 'SRC20'): src='3C236'
	if(code == 'SRC21'): src='3C237'
	if(code == 'SRC22'): src='3C245'
	if(code == 'SRC23'): src='1055+018'
	if(code == 'SRC24'): src='3C263.1'
	if(code == 'SRC25'): src='3C273'
	if(code == 'SRC26'): src='4C32.44'
	if(code == 'SRC27'): src='4C25.43'
	if(code == 'SRC28'): src='3C286'
	if(code == 'SRC29'): src='4C12.50'
	if(code == 'SRC30'): src='3C298'
	if(code == 'SRC31'): src='UGC09799'
	if(code == 'SRC32'): src='4C04.51'
	if(code == 'SRC33'): src='3C327.1'
	if(code == 'SRC34'): src='PKS1607'
	if(code == 'SRC35'): src='J1613'
	if(code == 'SRC36'): src='3C346'
	if(code == 'SRC37'): src='J1651'
	if(code == 'SRC38'): src='3C390'
	if(code == 'SRC39'): src='J1853'
	if(code == 'SRC40'): src='3C395'
	if(code == 'SRC41'): src='4C33.48'
	if(code == 'SRC42'): src='PKS1944'
	if(code == 'SRC43'): src='3C409'
	if(code == 'SRC44'): src='3C410'
	if(code == 'SRC45'): src='4C31.56'
	if(code == 'SRC46'): src='B2050'
	if(code == 'SRC47'): src='3C433'
	if(code == 'SRC48'): src='PKS2127'
	if(code == 'SRC49'): src='J2136'
	if(code == 'SRC50'): src='J2232'
	if(code == 'SRC51'): src='3C454.3'
	if(code == 'SRC52'): src='3C459'

	return src

## Get sources & number of sources ##
 # params string srcs All 127 lines-of-sight
 # return number-of-src & list of sources
 #
 # version 10/2016
 # Author Van Hiep ##
def get_src_and_nr(srcs):
	all_src = np.unique(srcs)

	name = []
	for sc in all_src:
		name.append(get_name(sc))

	return len(all_src), name

## Read Claire's results on 31 src #
 #
 # params string fname Filename
 # return dict info of Ts
 # 
 # version 10/2016
 # Author Van Hiep ##
def read_claire_res(fname = 'sub_data/column_densities.txt'):
	cols = ['src','cnm','cnm_er','wnm','wnm_er','nhi','nhi_er','thin','thin_er']
	fmt  = ['s',  'f',   'f',     'f',  'f',     'f',   'f'     , 'f', 'f'   ]
	data = restore(fname, 2, cols, fmt)
	return data.read()

## Read 23 OH sources #
 #
 # params string fname Filename
 # return dict info of 23-OH sources
 # 
 # version 10/2016
 # Author Van Hiep ##
def read_oh_src(fname = '../sub_data/nh_nhi_noh_19src_er.txt'):
	cols = ['idx','src','noh','noh_er','nhi','nhi_er','nh','nh_er']
	fmt  = ['i','s',  'f','f',  'f','f',  'f','f']
	data = restore(fname, 3, cols, fmt)
	return data.read()

## in 51 src, get common src and different src with Claire's results ##
 # params list src51 51-src from data
 # params list src30 30-src from Claire's results
 #
 # return list common Common src 
 #        list diff   different src 
 #
 # version 10/2016
 # Author Van Hiep ##
def get_common_diff_src(src51,src30):
	common = []
	diff   = []
	for sc in src51:
		if(sc in src30):
			common.append(sc)
		else:
			diff.append(sc)

	return common,diff

## Different src from Claire's results that not in Data##
 # params list src51 51-src from data
 # params list src30 30-src from Claire's results
 #
 # return list common Common src 
 #        list diff   different src 
 #
 # version 10/2016
 # Author Van Hiep ##
def get_diff_src(src51,src30):
	common = []
	diff   = []
	for sc in src30:
		if(sc in src51):
			common.append(sc)
		else:
			diff.append(sc)

	return common,diff

## get common src and different src between Claire's results and 19-OH src ##
 #
 # params list src19 19-OH sources
 # params list src30 30-src from Claire's results
 #
 # return list common Common src 
 #        list diff   different src 
 #
 # version 10/2016
 # Author Van Hiep ##
def get_common_diff_ohsrc_with_claire(src19,src30):
	common = []
	diff   = []
	for sc in src30:
		if(sc in src19):
			common.append(sc)
		else:
			diff.append(sc)

	return common,diff

## Get indexes of a src in list-of-sources ##
 # 
 # params string src Name of src
 # params list sr_names List of sources
 #
 # return list idx List of Indexes
 #
 # version 10/2016
 # Author Van Hiep ##
def get_idx(src,src_names):
	idx = []
	for i in range(len(src_names)):
		if(src == src_names[i]):
			idx.append(i)
	return idx

## Get VLSR ##
 # 
 # params dict hdr1info Header-Info
 # params int npatt Index of Observation
 #
 # return 1D-array vlsr 
 #
 # version 10/2016
 # Author Van Hiep ##
def get_vlsr(hdr1info,freq,npatt = 0):
	#GET THE CHANNEL SEPARATION in HZ...
    chnlsep  = 1e6*hdr1info[npatt,4]/2048.
    vlsrcntr = hdr1info[npatt,7]
    cntrfreq = hdr1info[npatt,5]
    nch      = (freq - 1e6*cntrfreq)/chnlsep

    if (cntrfreq < 1500.):
    	cntrchnl = 1023.
    	vlsr = vlsrcntr - (np.arange(2048)-cntrchnl)*chnlsep*2.99792458e5/(1e6*hdr1info[npatt,5])
    	dn   = nch
    	vlsr = vlsr + dn*chnlsep*2.99792458e5/(1e6*hdr1info[npatt,5])
    else:
    	cntrchnl = 1024.
    	vlsr = vlsrcntr - (np.arange(2048)-cntrchnl)*chnlsep*2.99792458e5/(1e6*hdr1info[npatt,5])
    	dn   = nch
    	vlsr = vlsr + dn*chnlsep*2.99792458e5/(1e6*hdr1info[npatt,5])

    return vlsr


##============== MAIN ==============##
## 172 observations, 0->20 and 21-171
## p1_claire_16feb2013-23jun_2013_a2770_bd1.sav
## 0->20: 1612.23 MHz, Bandwidth: 1.5625MHz
## 21->171: 1666.38 MHz, Bandwidth: 3.125MHz (1667MHz in the range, 1665MHz in Noise)

## p1_claire_16feb2013-23jun_2013_a2770_bd2.sav
## 0->20: 1665.40 MHz, Bandwidth: 1.5625MHz
## 21->171: 1720.53 MHz, Bandwidth: 3.125MHz 

## p1_claire_16feb2013-23jun_2013_a2770_bd3.sav
## 0->20: 1667.36 MHz, Bandwidth: 1.5625MHz
## 21->171: 1300.00 MHz, Bandwidth: 3.125MHz 

## For OH, I will analyse:
## BD1: 21->171 in 1666.38MHz
## BD3: 0->20 in 1667.36MHz

## This script to test 3C18 in SPONGE and MS.
## In MS, Absorption line at -7.8 (km/s)
## Now check in SPONGE.

data     = readsav('../data/p1_claire_16feb2013-23jun_2013_a2770_bd2.sav')
srcs     = np.asarray(data['HDRSRCNAME'])
nr,src51 = get_src_and_nr(srcs) ## 51 sources observed

sc_names = []
for i in range(len(srcs)):
	sc_names.append(get_name(srcs[i]))

# ## Results from Claire, nhi, cnm, wnm, nhi_thin and uncertainties ##
# res      = read_claire_res('sub_data/column_densities.txt')
# src30    = res['src'] ## 30-src from Claire's results

# ## Results for OH: nh, nhi, noh and uncertainties ##
# oh       = read_oh_src(fname = '../sub_data/nh_nhi_noh_19src_er.txt')
# oh_src   = oh['src']

# # common,diff = get_common_diff_src(src51,src30)
# # common1,diff1 = get_diff_src(src51,src30)
# comsrc,diffsrc = get_common_diff_ohsrc_with_claire(oh_src,src30)

## Constants ##
nch = 2048
c   = 2.99792458e5
pi  = 2.*np.arccos(0.)

f1 = 1612.231*1e6
f2 = 1665.402*1e6
f3 = 1667.359*1e6
f4 = 1720.530*1e6

## Data ##
stkon  = data['stkon']
stkoff = data['stkoff']
hdr1   = data['hdr1info']

src = '4C16.09'
idx = get_idx(src,sc_names)
print idx
rc = idx[0]

f0  = 1e6*hdr1[rc,5] # Center Freq
bw  = 1e6*hdr1[rc,4] # Bandwidth
df  = 1e6*bw/nch # Channel Separation

print f0, f2-f0,f3-f0

v0 = hdr1[rc,7]
a  = bw/nch
b  = f0-bw/2

# f = np.zeros(nch, dtype=np.float32)
f = []
for i in range(nch):
	f.append(a*i+b-f0)

v = get_vlsr(hdr1,f2,npatt=rc)

t_on  = 0.
for i in idx:
	t_on = t_on + stkon[i,0,:]

t_on = t_on/len(idx)

t_off = stkoff[rc,0,:]

plt.plot(v,t_on)
plt.grid()
# plt.axvline(x=f2-f0, ymin=-1000., ymax = 1000., linewidth=2, color='k')
# plt.axvline(x=f3-f0, ymin=-1000., ymax = 1000., linewidth=2, color='k')
plt.title(src + ', f0 = ' +str(f0) )
plt.show()