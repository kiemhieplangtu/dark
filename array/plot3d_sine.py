from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

pi = 2.*np.arccos(0.)
X = np.arange(-2, 2, 1./16.)
Y = np.arange(-2, 2, 1./16.)
X, Y = np.meshgrid(X, Y)
Z = np.sin(2.*pi*X+2.*pi*1.5*Y)

# fig = plt.figure()
# ax = fig.gca(projection='3d')
# surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False)
# # ax.set_zlim(-1.01, 1.01)

# ax.zaxis.set_major_locator(LinearLocator(10))
# ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# fig.colorbar(surf, shrink=0.5, aspect=5)

# plt.show()

tmin = np.min(Z)
tmax = np.max(Z)

print tmin
print tmax

lv = np.arange(tmin, tmax ,0.1)

plt.figure()
ax = plt.gca()
# cs = plt.contour(xl, xb, tb, lv)
# plt.clabel(cs, inline=1, fontsize=10, colors=('r', 'g', 'b'),
	# 	                  origin='lower',
	# 	                  extend='both')
plt.contourf(X, Y, Z, lv, cmap=plt.cm.jet)
plt.colorbar()

# for i in range(len(l_src)):
# 	if(b_src[i] < 10. and b_src[i] > -10.):
# 		c = plt.Circle((l_src[i], b_src[i]), 1., color='magenta', linewidth=2, fill=False)
# 		ax.add_artist(c)
# 		plt.text(l_src[i], b_src[i], srcs[i], fontsize=10, color='white')

# plt.hist(tb)
plt.title('Sin')
plt.xlabel('X')
plt.ylabel('Y')
# plt.ylim(-10., 10.)
plt.grid()
plt.show()