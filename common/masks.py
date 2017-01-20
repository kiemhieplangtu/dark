## Import libs ##
import sys, os
import numpy             as     np
import matplotlib.pyplot as     plt
import healpy            as     hp

## Class - Get the Dust opacity and it Standard Deviation for each pixel (l,b) ##
 # 
 # Eg:
 #  opc           = masks()
 # 	opacity,error = opc.get_dust_opacity(msk,l,b)
 # All Masks:
 #   lowest  weight = 7
 #   low_hi  weight = 6
 #   sth_cap weight = 5
 #   g35     weight = 4
 #   g45     weight = 3
 #   gt15    weight = 2
 #   g56     weight = 1
 #   whole   weight = 0
 #
 # version 11/2016
 # author Nguyen Van Hiep
 ##
class masks:

	## Initiate function ##
	 #
	 # params 
	 # return void
	 #
	 # version 011/2016 
	 # author Nguyen Van Hiep
	 ##
	def __init__(self):
		self.ret = ''

	## Get the Dust opacity and it Standard Deviation for each pixel (l,b) #
	 #
	 # params float l Gal-longitude
	 # params float b Gal-latitude
	 # params array msk Map of masks
	 #
	 # return list [f,er] Dust opacity and it Standard Deviation 
	 # 
	 # version 11/2016
	 # Author Van Hiep ##
	def get_dust_opacity(self,msk,l,b):	
		deg2rad = np.pi/180.
		nside   = hp.get_nside(msk)
		opac    = [8.4e-27, 7.1e-27, 7.0e-27, 6.8e-27, 6.5e-27, 6.5e-27, 6.6e-27, 7.9e-27 ]
		err     = [3.0e-27, 1.9e-27, 2.0e-27, 1.8e-27, 1.8e-27, 1.9e-27, 1.7e-27, 1.9e-27 ]

		theta   = (90.0 - b)*deg2rad
		phi     = l*deg2rad
		pix     = hp.ang2pix(nside, theta, phi, nest=False)
		val     = msk[pix]

		return [ opac[int(val)], err[int(val)] ]