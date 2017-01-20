import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class
import matplotlib.pyplot as plt
import numpy             as np
import healpy            as hp
import pylab             as pl
import matplotlib        as mpl
import copy
import operator


from matplotlib           import cm
from matplotlib.ticker    import LinearLocator, FormatStrFormatter
from restore              import restore
from astropy.io           import fits
from mpfit                import mpfit
from mpl_toolkits.mplot3d import Axes3D

## Read info of 26 no-CO sources #
 # l,b, nhi, and nhi_error
 #
 # params string fname Filename
 # return dict info
 # 
 # version 11/2016
 # Author Van Hiep ##
def read_23oh_src(fname = '../oh/result/23src_with_oh.txt'):
	cols = ['src','l','b']
	fmt  = ['s',  'f','f']
	data = restore(fname, 2, cols, fmt)
	dat  = data.read()
	return dat

## 2D gaussian function ##
 #
 # params 
 # params 
 #
 # return 
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def myfunc(p, fjac=None, x=None, y=None, z=None, err=None, fit=True):
	model  = p[0] + p[1]*np.exp( -0.5*( ( ( (x-p[2])/p[3] )**2 ) + ( ( (y-p[4])/p[5] )**2 ) ) )
	status = 0
	
	if(fit):
		return [status, (z-model)/err]
	else:
		return model

#================= MAIN ========================#
## Define constants ##
deg2rad = np.pi/180.
beam    = 3.5
dbeam   = beam/120.0 # Beam = 3.5' -> dbeam = beam/60/2
dbeam   = 0.5
edge    = dbeam/40. # To plot a rectangle about the source
# info    = read_23oh_src('../oh/result/23src_with_oh.txt')

print 'mmm: ', 2.8+6.32*(1420./1665.402)**2.8
print 'mmm: ', 2.8+22.6*(408./1665.402)**2.8

file = 'fits3c123_3deg.bin'
file = 'fits3c131_1.5deg.bin'
file = 'fits3c154_3deg.bin'
file = 'fits3c410_3deg.bin'
file = 'fitsp042820.bin'
hdulist = fits.open('data/'+file)
hdu = hdulist[0]

shape = hdu.data.shape
xlen = shape[1]
ylen = shape[2]

# print hdu.data[0,:,:]
print 'Mean: ',np.mean(hdu.data)
# print hdu.header
plt.imshow(hdu.data[0,:,:], origin='lower')
plt.title(file)
plt.show()



# Data & Errors
X = np.arange(0., xlen, 1.)
Y = np.arange(0., ylen, 1.)

xdat, ydat = np.meshgrid(X, Y)
xdata = xdat.ravel()
ydata = ydat.ravel()
zdat  = hdu.data[0,:,:]
zdata = zdat.ravel()
err  = 150.*np.ones(zdat.shape)
err = err.ravel()

## Fit ##
lguess  = [3750., 2500.,6.,1.,6.,1.]  ## Guess the initial values
npar    = len(lguess)    ## Number of params
guessp  = np.array(lguess, dtype='float64') ## Nothing, put initial values into array
plimd   = [[False,False]]*npar  ## Limit the Params or not
plims   = [[0.,0.]]*npar        ## If limit the Params, define the range of params, if not just put Zeros
parbase = {'value': 0., 'fixed': 0, 'parname': '', 'limited': [0, 0], 'limits': [0., 0.]}
pname   = ['offset','amp','x0','xwid','y0','ywid']  ## Names of params
pfix    = [False]*npar  ## Fix the params or not

parinfo = []
for i in range(len(guessp)):
	parinfo.append(copy.deepcopy(parbase))

for i in range(len(guessp)):
	parinfo[i]['value']   = guessp[i]
	parinfo[i]['fixed']   = pfix[i]
	parinfo[i]['parname'] = pname[i]
	parinfo[i]['limited'] = plimd[i]

##  1665 ###
x    = xdata.astype(np.float64)
y    = ydata.astype(np.float64)
z    = zdata.astype(np.float64)
er   = err.astype(np.float64)

fa   = {'x':x, 'y':y, 'z':z, 'err':er} ## Arrange the params before fitting
mp   = mpfit(myfunc, guessp, parinfo=parinfo, functkw=fa, quiet=True)  ## Do the fit

## ********* Results ********* ##
print '********* Results *********'
p   = mp.params
per = mp.perror
for i in range(len(parinfo)):
	print "%s = %f +/- %f" % (parinfo[i]['parname'],p[i],per[i])  ## Print out the results


print '##===##'
print 'Offset: ', "%f +/- %f" % (p[0]/1000.,per[0]/1000.)
print 'Amplitude: ', "%f +/- %f" % (p[1]/1000.,per[1]/1000.)


# sys.exit()


# Plot the original, fit & residual
fig = pl.figure(figsize=(18,4.3))

ax1 = fig.add_subplot(1,3,1)
xyData = hdu.data[0,:,:]
cax1 = ax1.imshow(xyData, origin='lower',cmap=mpl.cm.jet)
cbar1=fig.colorbar(cax1, pad=0.0)
# ax1.scatter(X, Y, c=Z, s=40, cmap=mpl.cm.jet)
ax1.set_title("Data: " + file)
# ax1.set_xlim(0, shape[-1]-1)
# ax1.set_ylim(0, shape[-2]-1)
ax1.set_aspect('equal')

ax2 = fig.add_subplot(1,3,2)
xyDataFit = myfunc(p, x=xdata, y=ydata, fit=False)
xyDataFit.shape = zdat.shape
cax2 = ax2.imshow(xyDataFit, origin='lower', cmap=mpl.cm.jet)
cbar2=fig.colorbar(cax2, pad=0.0)
ax2.set_title("2D Gaussian Fit")

ax3 = fig.add_subplot(1,3,3)
xyDataRes = xyData - xyDataFit
cax3 = ax3.imshow(xyDataRes, origin='lower', cmap=mpl.cm.jet)
cbar2=fig.colorbar(cax3, pad=0.0)
ax3.set_title("Residual")

pl.show()

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(xdat, ydat, zdat, rstride=1, cstride=1, cmap=mpl.cm.jet,
                       linewidth=0, antialiased=False)
# ax.set_zlim(-1.01, 1.01)

# ax.zaxis.set_major_locator(LinearLocator(10))
# ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()

xf = np.arange(0., 13., 0.2)
yf = np.arange(0., 13., 0.2)
xdat, ydat = np.meshgrid(xf, yf)
zf = myfunc(p, x=xdat.ravel(), y=ydat.ravel(), fit=False)
zf.shape = xdat.shape

fig = plt.figure()
ax = fig.gca(projection='3d')

surf = ax.plot_surface(xdat, ydat, zf, rstride=1, cstride=1, cmap=mpl.cm.jet,
                       linewidth=0, antialiased=False)
# ax.set_zlim(-1.01, 1.01)

# ax.zaxis.set_major_locator(LinearLocator(10))
# ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()

