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

## Calculate the uncertainties of ratio a/b #
 #
 # params float a
 # params float b
 # params float aer
 # params float ber
 #
 # return float ret Uncertainty of a/b
 # 
 # Author Van Hiep ##
def uncertainty_of_ratio(a, b, aer, ber):
	r  = np.abs(a/b)
	d1 = aer/a
	d1 = d1*d1
	d2 = ber/b
	d2 = d2*d2
	d  = r * np.sqrt(d1+d2)
	return d

## Calculate the uncertainties of a-b #
 #
 # params float aer
 # params float ber
 #
 # return float ret Uncertainty of a-b
 # 
 # Author Van Hiep ##
def uncertainty_of_diff(aer, ber):
	d1 = aer**2
	d2 = ber**2
	d  = np.sqrt(d1+d2)
	return d	

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

r = taucnm2/taucnm1
d = tspincnm2-tspincnm1

dr = uncertainty_of_ratio(taucnm2, taucnm1, sigtaucnm2, sigtaucnm1)
dd = uncertainty_of_diff(sigtspincnm2, sigtspincnm1)

## Plot
fig          = plt.figure(figsize=(16,16))
ax           = fig.add_subplot(111)                                 
major_xticks = np.arange(-20., 20., 5.)                                              
minor_xticks = np.arange(-20., 20., 1.)
major_yticks = np.arange(0.5, 5., 0.5)                                              
minor_yticks = np.arange(0.5, 5., 0.1)

# ax.plot(d, r, 'k.', markersize=16)
plt.errorbar(d, r,xerr=dd, yerr=dr, color='k', marker='o', ls='None', markersize=8, \
	markeredgecolor='k', markeredgewidth=1, capsize=0)
plt.axhline(y=1.8, xmin=-20., xmax=20., color='k', ls=':', lw=3)
ax.axvspan(-2., 2., ymin=0., ymax=10., alpha=0.3, color='grey')

plt.title('', fontsize=30)
plt.ylabel(r'$R_{67/65}=\tau_{1667}/\tau_{1665}$', fontsize=35)
plt.xlabel(r'$\Delta T_{EX}=T_{ex}(1667) - T_{ex}(1665)$', fontsize=35)
ax.set_xticks(major_xticks)                                                       
ax.set_xticks(minor_xticks, minor=True)                                           
ax.set_yticks(major_yticks)                                                       
ax.set_yticks(minor_yticks, minor=True)
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=18)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)
# plt.legend(loc='upper right', fontsize=18)
plt.grid(False)
plt.xlim(-10.,20.)
plt.ylim(0.5,3.5)
plt.show()