import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import operator

from numpy    import array
from restore  import restore
from plotting import cplot

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import scipy.interpolate as interp

fname = 'config1.txt'
cols = ['xl','xm','wi']
fmt  = ['f','f','f']
data = restore(fname, 0, cols, fmt)
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

## PLOT CONTOURS ##
lv = np.arange(tmin, tmax ,0.1)

xl1 = np.unique(xl)
xm1 = np.unique(xm)

xlen = len(xl1)
ylen = len(xm1)
w   = wi.reshape((xlen,ylen))
plt.figure()
ax = plt.gca()
# cs = plt.contour(xl, xb, tb, lv)
# plt.clabel(cs, inline=1, fontsize=10, colors=('r', 'g', 'b'),
	# 	                  origin='lower',
	# 	                  extend='both')
plt.contourf(xm1, xl1, w, lv, cmap=plt.cm.jet)
plt.colorbar()

# for i in range(len(l_src)):
# 	if(b_src[i] < 10. and b_src[i] > -10.):
# 		c = plt.Circle((l_src[i], b_src[i]), 1., color='magenta', linewidth=2, fill=False)
# 		ax.add_artist(c)
# 		plt.text(l_src[i], b_src[i], srcs[i], fontsize=10, color='white')

# plt.hist(tb)
plt.title('Dirty Map')
plt.xlabel('m (arcsec)')
plt.ylabel('l (arcsec)')
# plt.xlim(-0.1, 0.1)
# plt.ylim(-0.1, 0.1)
plt.grid()
plt.show()


## PLOT 3D ##
plotx,ploty, = np.meshgrid(np.linspace(np.min(xl),np.max(xl),100),\
                           np.linspace(np.min(xm),np.max(xm),100))
plotz = interp.griddata((xl,xm),wi,(plotx,ploty),method='linear')


fig = plt.figure()
ax = fig.gca(projection='3d')
# surf = ax.plot_surface(plotx,ploty,plotz, rstride=8, cstride=8, alpha=0.3)
surf = ax.plot_surface(plotx,ploty,plotz, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# ax.set_xlim3d(-0.1, 0.1)
# ax.set_ylim3d(-0.1, 0.1)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()