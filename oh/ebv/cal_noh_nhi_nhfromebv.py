import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import operator

from numpy    import array
from restore  import restore
from plotting import cplot

## Compute N(OH) ##
 #
 # params dict dat Data on N(OH)
 # return dict ret N(OH) of each source
 #
 # version 09/2016 
 # author Nguyen Van Hiep ##
def get_noh_for_each_src(fname):
	cols = ['idx','src','tau','v0','wid','tex','tex_er','noh','noh_er']
	fmt  = ['i','s','f','f','f','f','f','f','f']
	data = restore(fname, 2, cols, fmt)
	dat  = data.read()
	noh  = dat['noh']
	er2  = dat['noh_er']
	src  = dat['src']

	ret  = {}
	for i in range(0,len(dat['src'])):
		if dat['src'][i] not in ret.keys():
			ret[src[i]] = {}
			ret[src[i]]['noh'] = noh[i]
			ret[src[i]]['er2'] = er2[i]**2
		else:
			ret[src[i]]['noh'] = ret[src[i]]['noh'] + noh[i]
			ret[src[i]]['er2'] = ret[src[i]]['er2'] + er2[i]**2

	return ret

##================= MAIN ========================##

## N(HI) and N(H) from Dust ##
cols = ['idx','src', 'l', 'b', 'nhi','nhi_er', 'nh','nh_er']
fmt  = ['i','s','f','f','f','f','f','f']
data = restore('../../dust/ebv2nh/result/nh_nhi_uncert_78src_from_ebv.txt', 3, cols, fmt)
dat  = data.read()

nhi    = dat['nhi']
nhi_er = dat['nhi_er']
nh     = dat['nh']
nh_er  = dat['nh_er']
src    = dat['src']
idx    = dat['idx']

## N(OH) ##
# noh1  = get_noh_for_each_src(fname = 'result/noh1_src96_er.txt')
noh2    = get_noh_for_each_src(fname = '../result/noh2_src96_er.txt')
nhi2    = []
nhi_er2 = [] # WARNING: Here Sigma^2

n       = 0
for sc in noh2:
	if((sc in src) and (noh2[sc] > 0.)):
		nhi2.append(nhi[src.index(sc)])
		nhi_er2.append(nhi_er[src.index(sc)])

		print("{}   {}\t{}\t{}\t{}\t{}\t{}    {}" \
			.format(n, sc, round(noh2[sc]['noh'],8), round(np.sqrt(noh2[sc]['er2']),8), nhi[src.index(sc)],nhi_er[src.index(sc)], nh[src.index(sc)],nh_er[src.index(sc)]) )
		n = n + 1