import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import operator

from numpy    import array
from restore  import restore
from plotting import cplot

fname = 'test.txt'
cols = ['xl','xm','wi']
fmt  = ['f','f','f']
data = restore(fname, 2, cols, fmt)
dat  = data.read()
xl   = dat['xl']
xm   = dat['xm']
wi   = dat['wi']

xl = np.asarray(xl)
xm = np.asarray(xm)
wi = np.asarray(wi)

print len(wi)

tmin = np.min(wi)
tmax = np.max(wi)

print tmin
print tmax

lv = np.arange(tmin, tmax ,0.1)

xl = np.unique(xl)
xm = np.unique(xm)

xlen = len(xl)
ylen = len(xm)
w    = wi.reshape((xlen,ylen))
# w    = np.rot90(w,k=1)

plt.figure()
ax = plt.gca()
# cs = plt.contour(xl, xb, tb, lv)
# plt.clabel(cs, inline=1, fontsize=10, colors=('r', 'g', 'b'),
	# 	                  origin='lower',
	# 	                  extend='both')
plt.contourf(xm, xl, w, lv, cmap=plt.cm.jet)
plt.colorbar()

# for i in range(len(l_src)):
# 	if(b_src[i] < 10. and b_src[i] > -10.):
# 		c = plt.Circle((l_src[i], b_src[i]), 1., color='magenta', linewidth=2, fill=False)
# 		ax.add_artist(c)
# 		plt.text(l_src[i], b_src[i], srcs[i], fontsize=10, color='white')

# plt.hist(tb)
plt.title('WI , Min:'+str(tmin)+', Max:'+str(tmax))
plt.xlabel('Galatic l')
plt.ylabel('Galatic b')
# plt.xlim(-0.1, 0.1)
# plt.ylim(-0.1, 0.1)
plt.grid()
plt.show()


