'''
cmbmap.py
Plot intensity skymap from Planck, and get power spectrum.
NOTE that I have not been able to install healpy in Canopy on Windows :(
'''

import numpy as np
import matplotlib.pyplot as plt
# from pylab import * 
import healpy as hp

# here are a bunch of different maps; the last one seems to work best
# http://healpy.readthedocs.org/en/latest/tutorial.html
# http://irsa.ipac.caltech.edu/data/Planck/release_1/all-sky-maps/previews/COM_CompMap_CMB-commrul_2048_R1.00/index.html
fn = '/Users/ajw/Downloads/COM_CompMap_CMB-commrul_2048_R1.00.fits'
# http://irsa.ipac.caltech.edu/data/Planck/release_1/all-sky-maps/previews/COM_CompMap_CMB-smica_2048_R1.20/index.html
fn = '/Users/ajw/Downloads/COM_CompMap_CMB-smica_2048_R1.20.fits'
# http://irsa.ipac.caltech.edu/data/Planck/release_1/all-sky-maps/previews/COM_CompMap_CMB-nilc_2048_R1.20/index.html
fn = '/Users/ajw/Downloads/COM_CompMap_CMB-nilc_2048_R1.20.fits'
map_I = hp.read_map(fn)

# plot the Planck T map
#figure()
hp.mollview(map_I, coord='G', title='Planck R1, Histogram equalized Galactic', unit='mK', norm='hist')
#hp.graticule(dpar=30,dmer=30)

# analyze the map, get Cl power spectrum:
LMAX = 2048
cl = hp.anafast(map_I, lmax=LMAX)
#hp.write_cl('cl_Planck.fits', cl)
ell = np.arange(len(cl))
ecl = ell * (ell+1) * cl/(2.*np.pi)

# Get the Planck fit result: http://pla.esac.esa.int/pla/aio/planckResults.jsp?
# http://irsa.ipac.caltech.edu/data/Planck/release_1/ancillary-data/previews/COM_PowerSpect_CMB_R1.10/index.html
fn = 'COM_PowerSpect_CMB_R1.10.txt'
f = open(fn)
lines = f.readlines()
f.close()

# extract low-ell and high-ell power spectra
a = np.genfromtxt(lines[10:58],delimiter=",").transpose()
b = np.genfromtxt(lines[68:],delimiter=",").transpose()

# I don't know how to get the 2nd (high-ell) table from thre fits file
#fn = 'COM_PowerSpect_CMB_R1.10.fits'
#cls = hp.read_cl(fn)

# "fit" from CAMB: http://lambda.gsfc.nasa.gov/toolbox/tb_camb_form.cfm
fn = 'test_scalCls.dat'
fn = 'camb_53319956_scalcls.dat.txt'
f = open(fn)
lines = f.readlines()
f.close()
c = np.genfromtxt(lines).transpose()

# here we see a significant unexplained discrepancy mid-ell.
plt.figure()
ax = plt.subplot()
plt.plot(ell,ecl,'k',label='anafast fit')
plt.errorbar(a[0],a[1],yerr=a[3],fmt='bo',label='COM_PowerSpect_CMB_R1.10')
plt.errorbar(b[0],b[3],xerr=b[2]-b[0],yerr=b[4],fmt='bo',label='COM_PowerSpect_CMB_R1.10')
plt.plot(c[0],c[1],'r',label='camb_53319956_scalcls')
#ax.set_xscale("log", nonposx='clip')
plt.xlim([2.,2500.])
plt.xlabel(r'$\ell$')
plt.ylabel(r'$\ell(\ell+1) C_\ell/2\pi$')
plt.title('Data from Planck, prediction from CAMB')
plt.legend()
plt.grid()

# part of the discrepancy come from masking, so reanalyze map with mask:
mask = hp.read_map('HFI_PowerSpect_Mask_2048_R1.10.fits').astype(np.bool)
map_I_masked = hp.ma(map_I)
map_I_masked.mask = np.logical_not(mask)
clm = hp.anafast(map_I_masked, lmax=LMAX)
eclm = ell * (ell+1) * clm/(2.*np.pi) * 2.5   # note the norm fudge-factor; should come from the mask

# plot results:
plt.figure()
ax = plt.subplot()
plt.plot(ell,eclm,'k',label='anafast fit - masked')
plt.errorbar(a[0],a[1],yerr=a[3],fmt='bo',label='COM_PowerSpect_CMB_R1.10')
plt.errorbar(b[0],b[3],xerr=b[2]-b[0],yerr=b[4],fmt='bo',label='COM_PowerSpect_CMB_R1.10')
plt.plot(c[0],c[1],'r',label='camb_53319956_scalcls')
#ax.set_xscale("log", nonposx='clip')
plt.xlim([2.,2500.])
plt.xlabel(r'$\ell$')
plt.ylabel(r'$\ell(\ell+1) C_\ell/2\pi$')
plt.title('Data from Planck - masked')
plt.legend()
plt.grid()
plt.show()


'''
# for fun, make a simulated map
nside = hp.get_nside(pl)
smap = hp.synfast(cl,nside)
hp.mollview(smap, coord='G', title='simulated data', unit='mK', norm='hist')
hp.graticule(dpar=30,dmer=30)
LMAX = 2048
cl1 = hp.anafast(smap, lmax=LMAX)
ell = arange(len(cl))
ecl1 = ell * (ell+1) * cl1/(2.*np.pi)
plt.figure()
ax = plt.subplot()
plt.plot(ell,ecl,'b',label="into to synfast")
plt.plot(ell,ecl1,'r',label="output from synfast")
plt.xlim([2.,2500.])
plt.xlabel(r'$\ell$')
plt.ylabel(r'$\ell(\ell+1) C_\ell/2\pi$')
plt.legend()
plt.grid()
plt.show()
'''
