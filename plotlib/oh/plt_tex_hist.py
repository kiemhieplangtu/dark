import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt
import matplotlib
import pandas            as pd

from scipy.optimize import curve_fit
from scipy.io.idl   import readsav
from numpy          import array
from restore        import restore
from plotting       import cplot
from gauss_fit      import gfit

## Create a line ##
 #
 # params list x x-data
 # params list y y-data
 # params string label Label of line
 # params dict prop Properties of line
 # return dict ret All infor of line
 #
 # version 07/2016 
 # author Nguyen Van Hiep ##
def daddsada(tb):
	return ''

#============== MAIN ==============#
dat                 = readsav('data/doall_fitgauss_stokesiover2.sav')

## For OH1665 Absorption 
src1                = dat.amg_1665.sname
ell1                = dat.amg_1665.ell
bee1                = dat.amg_1665.bee
ng1                 = dat.amg_1665.nrgauss
zro1                = dat.amg_1665.zrocnm
taucnm1             = dat.amg_1665.taucnm
cencnm1             = dat.amg_1665.cencnm
widcnm1             = dat.amg_1665.widcnm
tspincnm1           = dat.amg_1665.tspincnm
sigzrocnm1          = dat.amg_1665.sigzrocnm
sigtaucnm1          = dat.amg_1665.sigtaucnm
sigcencnm1          = dat.amg_1665.sigcencnm
sigwidcnm1          = dat.amg_1665.sigwidcnm
sigtspincnm1        = dat.amg_1665.sigtspincnm
continuum_defln1    = dat.amg_1665.continuum_defln
sigcontinuum_defln1 = dat.amg_1665.sigcontinuum_defln
continuum_em1       = dat.amg_1665.continuum_em


## For OH1667 Absorption 
src2                = dat.amg_1667.sname
ng2                 = dat.amg_1667.nrgauss
ell2                = dat.amg_1667.ell
bee2                = dat.amg_1667.bee
zro2                = dat.amg_1667.zrocnm
taucnm2             = dat.amg_1667.taucnm
cencnm2             = dat.amg_1667.cencnm
widcnm2             = dat.amg_1667.widcnm
tspincnm2           = dat.amg_1667.tspincnm
sigzrocnm2          = dat.amg_1667.sigzrocnm
sigtaucnm2          = dat.amg_1667.sigtaucnm
sigcencnm2          = dat.amg_1667.sigcencnm
sigwidcnm2          = dat.amg_1667.sigwidcnm
sigtspincnm2        = dat.amg_1667.sigtspincnm
continuum_defln2    = dat.amg_1667.continuum_defln
sigcontinuum_defln2 = dat.amg_1667.sigcontinuum_defln
continuum_em2       = dat.amg_1667.continuum_em

## For OH1665 Emmission
esrc1                = dat.emg_1665.sname
ell2                 = dat.emg_1665.ell
bee2                 = dat.emg_1665.bee
eng1                 = dat.emg_1665.nrgauss
bsl1                 = dat.emg_1665.tbaseline
cont_em1             = dat.emg_1665.continuum_em
hgtwnm1              = dat.emg_1665.hgtwnm
cenwnm1              = dat.emg_1665.cenwnm
widwnm1              = dat.emg_1665.widwnm
fwnm1                = dat.emg_1665.fwnm
sigcont1             = dat.emg_1665.sigcontinuum
sighgtwnm1           = dat.emg_1665.sighgtwnm
sigcenwnm1           = dat.emg_1665.sigcenwnm
sigwidwnm1           = dat.emg_1665.sigwidwnm
sigfwnm1             = dat.emg_1665.sigfwnm

## For OH1667 Emmission 
esrc2                = dat.emg_1667.sname
ell2                 = dat.emg_1667.ell
bee2                 = dat.emg_1667.bee
eng2                 = dat.emg_1667.nrgauss
bsl2                 = dat.emg_1667.tbaseline
cont_em2             = dat.emg_1667.continuum_em
hgtwnm2              = dat.emg_1667.hgtwnm
cenwnm2              = dat.emg_1667.cenwnm
widwnm2              = dat.emg_1667.widwnm
fwnm2                = dat.emg_1667.fwnm
sigcont2             = dat.emg_1667.sigcontinuum
sighgtwnm2           = dat.emg_1667.sighgtwnm
sigcenwnm2           = dat.emg_1667.sigcenwnm
sigwidwnm2           = dat.emg_1667.sigwidwnm
sigfwnm2             = dat.emg_1667.sigfwnm

# for i in range(len(src1)):
# 	# print src1[i], taucnm1[i], cencnm1[i], widcnm1[i], tspincnm1[i]
# 	print i, src1[i], tspincnm1[i], tspincnm2[i], continuum_em1[i], continuum_em2[i]

## Plot

matplotlib.style.use('grayscale')
tbg = list(set(continuum_em1))

print len(tbg)

## Classic StepHist
fig         = plt.figure(figsize=(12,16))
ax          = fig.add_subplot(111)                                 
major_ticks = np.arange(0, 100, 5)                                              
minor_ticks = np.arange(0, 100, 1) 
kwargs      = dict(histtype='step', alpha=0.7, stacked=False, fill=False, range=(0.0,26.0), normed=False, bins=13, lw=3)
plt.hist(tspincnm1, hatch='', label='$T^{OH^{1665}}_{ex}$', color='r', **kwargs)
plt.hist(tspincnm2,hatch='', label='$T^{OH^{1667}}_{ex}$', color='k', **kwargs)
plt.hist(tbg, hatch='', label='$T_{bg}$', color='g', **kwargs)

bins             = np.arange(0.,20.1,2.)
tex1, bin_edges1 = np.histogram(tspincnm1, bins=bins, range=(0.0,20.0), normed=False)
tex2, bin_edges2 = np.histogram(tspincnm2, bins=bins, range=(0.0,20.0), normed=False)
bins             = np.arange(0.,26.1,2.)
tbg,  bin_edges3 = np.histogram(tbg,       bins=bins, range=(0.0,26.0), normed=False)

bincen1 = 0.5*(bin_edges1[1:] + bin_edges1[:-1])
bincen2 = 0.5*(bin_edges2[1:] + bin_edges2[:-1])
bincen3 = 0.5*(bin_edges3[1:] + bin_edges3[:-1])

plot1   = ax.plot(bincen1, tex1, ls='', color='r', marker='d',  markersize=8)
plot2   = ax.plot(bincen2, tex2, ls='',  color='k', marker='d', markersize=8)
plot3   = ax.plot(bincen3, tbg,  ls='',  color='g', marker='d', markersize=8)

plt.title('', fontsize=30)
plt.ylabel('Numbers', fontsize=35)
plt.xlabel('Temperature (K)', fontsize=35)
ax.set_xticks(major_ticks)                                                       
ax.set_xticks(minor_ticks, minor=True)                                           
ax.set_yticks(major_ticks)                                                       
ax.set_yticks(minor_ticks, minor=True)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)
plt.legend(loc='upper right', fontsize=18)
plt.grid(False)
plt.xlim(0.,26.)
plt.ylim(0.,20.)
plt.show()

## Use numpy histogram, no-color
tbg              = list(set(continuum_em1))
fig              = plt.figure(figsize=(12,16))
ax               = fig.add_subplot(111)
major_ticks      = np.arange(0, 100, 5)                                              
minor_ticks      = np.arange(0, 100, 1) 

bins             = np.arange(0.,20.1,2.)
tex1, bin_edges1 = np.histogram(tspincnm1, bins=bins, range=(0.0,20.0), normed=False)
tex2, bin_edges2 = np.histogram(tspincnm2, bins=bins, range=(0.0,20.0), normed=False)
bins             = np.arange(0.,26.1,2.)
tbg,  bin_edges3 = np.histogram(tbg,       bins=bins, range=(0.0,26.0), normed=False)

bincen1 = 0.5*(bin_edges1[1:] + bin_edges1[:-1])
bincen2 = 0.5*(bin_edges2[1:] + bin_edges2[:-1])
bincen3 = 0.5*(bin_edges3[1:] + bin_edges3[:-1])

plot1   = ax.plot(bincen1, tex1, ls='-', color='k', label='$T^{OH^{1665}}_{ex}$', marker='s', markersize=8,  lw=2)
plot2   = ax.plot(bincen2, tex2, ls=':',  color='k', label='$T^{OH^{1667}}_{ex}$', marker='*', markersize=13, lw=2)
plot3   = ax.plot(bincen3, tbg,  ls='--',  color='k', label='$T_{bg}$',             marker='D', markersize=9,  lw=1)

plt.title('', fontsize=30)
plt.ylabel('Numbers', fontsize=35)
plt.xlabel('Temperature (K)', fontsize=35)
ax.set_xticks(major_ticks)                                                       
ax.set_xticks(minor_ticks, minor=True)                                           
ax.set_yticks(major_ticks)                                                       
ax.set_yticks(minor_ticks, minor=True)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)
plt.legend(loc='upper right', fontsize=18)
plt.grid(False)
plt.xlim(0.,26.)
plt.ylim(0.,20.)
plt.show()

## Use numpy histogram, with colors
fig              = plt.figure(figsize=(12,16))
ax               = fig.add_subplot(111)
major_ticks      = np.arange(0, 100, 5)                                              
minor_ticks      = np.arange(0, 100, 1)

ax.plot(bincen1, tex1, ls='-', color='k', label='$T^{OH^{1665}}_{ex}$', marker='s', markersize=8, lw=2)
ax.plot(bincen2, tex2, ls=':',  color='r', label='$T^{OH^{1667}}_{ex}$', marker='*', markersize=13, lw=2)
ax.plot(bincen3, tbg,  ls='--',  color='g', label='$T_{bg}$',             marker='D', markersize=9, lw=1)

plt.title('', fontsize=30)
plt.ylabel('Numbers', fontsize=35)
plt.xlabel('Temperature (K)', fontsize=35)
ax.set_xticks(major_ticks)                                                       
ax.set_xticks(minor_ticks, minor=True)                                           
ax.set_yticks(major_ticks)                                                       
ax.set_yticks(minor_ticks, minor=True)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)
plt.legend(loc='upper right', fontsize=18)
plt.grid(False)
plt.xlim(0.,26.)
plt.ylim(0.,20.)
plt.show()