import sys, os
import numpy             as np
import matplotlib.pyplot as plt

## Class - fit gauss ##
 # 
 # Using, eg:
 # cols = ['index, 'source', 'temperature']                  # Columns
 # fmt  = ['i', 's', 'f']                                    # Format of each column (eg: ['s', 'i', 'f'])
 # x    = restore('filename.txt', skip_li!=s=4, cols, fmt)
 # dat  = x.read()
 #
 # version 07/2016 
 # author Nguyen Van Hiep ##
class gfit:
	# Formats of data columns #
	fmts = ['s', 'i', 'f']

	## Initiate function ##
	 #
	 # params 
	 # return void
	 #
	 # version 07/2016 
	 # author Nguyen Van Hiep ##
	def __init__(self, file=''):
		self.ret  = {}

	##  ##
	 #
	 # params 
	 # return 	
	 #
	 # version 08/2016 
	 # author Nguyen Van Hiep ##
	def lsum(self,vect):
		if(type(vect) is int):
			return vect
		else:
			return sum(vect)

	## Convert angle to a specified range by finding the angle modulo the extent of the range. ##
	 #
	 # params 
	 # params 
	 #
	 # return 
	 #
	 # version 08/2016 
	 # author Nguyen Van Hiep ##
	def tb_exp(self,xdata, zrocnm, hgtcnm, cencnm, widcnm, tspincnm, ordercnm, continuum, hgtwnm, cenwnm, widwnm, fwnm):
		#ZEROTH STEP IS TO REARRANGE CLOUDS IN ORDER OF 'ORDER'.
		# zro1 = zrocnm
		# hgt1 = hgtcnm[ordercnm]  # od = [0,1,2,3],  q1 = a[od], q1 = [0,2,3,4]
		# cen1 = cencnm[ordercnm]
		# wid1 = widcnm[ordercnm]
		# tspin1 = tspincnm[ordercnm]
		zro1   = zrocnm
		hgt1   = list(hgtcnm)
		cen1   = list(cencnm)
		wid1   = list(widcnm)
		tspin1 = list(tspincnm)

		#FIRST STEP IS TO CALCULATE THE OPACITY OF EACH COLD CLOUD...
		nchnl  = len(xdata)
		nrcnm  = len(hgt1)
		taucnm = np.zeros( (nchnl, nrcnm) )

		for nrc in range(nrcnm):
			tau1nrc = self.gcurv(xdata, zro1, hgt1[nrc], cen1[nrc], wid1[nrc])
			taucnm[:, nrc] = tau1nrc

		if (len(ordercnm) != 1):
			tausum = taucnm.sum(1)
		else:
			tausum = taucnm

		exp_tausum = np.exp(-tausum)
		exp_tausum = exp_tausum.reshape(nchnl,)

		##********** CALCULATE THE WNM CONTRIBUTION ********************
		##  EXPRESS THE WNM CONTRIBUTION AS A SUM OF GAUSSIANS:
		##	FWNM, ZROWNM, HGTWNM, CENWNM, WIDWNM
		tb_cont    = continuum* exp_tausum

		tb_wnm_tot = np.zeros(nchnl)
		nrwnm      = len(hgtwnm)
		for nrw in range(nrwnm):
			tb_wnm_nrw = self.gcurv(xdata, 0., hgtwnm[nrw], cenwnm[nrw], widwnm[nrw])
			tb_wnm_tot = tb_wnm_tot + tb_wnm_nrw*(fwnm[nrw] + (1.-fwnm[nrw])*exp_tausum)

		#*************** NEXT CALCULATE THE CNM CONTRIBUTION ****************

		tb_cnm_tot = np.zeros(nchnl)

		# BRIGHT TEMP OF EACH CNM CLUMP:
		tbclump = np.zeros((nchnl, nrcnm))
		for nrc in range(nrcnm):
			tbclump[:,nrc] = tspin1[nrc] * (1. - np.exp(-taucnm[:,nrc]))

		## Cho nay dung' ro`i, cong lai ca 1->m roi Tru` di tau_m. Lam` the' nay` de? tranh' di "[]" luc' nrc=0
		for nrc in range(nrcnm):
			temp        = np.reshape(taucnm[:, 0:nrc+1], (nchnl, nrc+1))
			tausum_nrc  = temp.sum(1)
			exp_tau_nrc = np.exp(taucnm[:, nrc] - tausum_nrc)
			tb_cnm_tot  = tb_cnm_tot + tspin1[nrc] * (1. - np.exp(-taucnm[:,nrc]) ) * exp_tau_nrc

		tb_tot = tb_cont+ tb_cnm_tot + tb_wnm_tot

		# print tb_wnm_tot, max(tb_wnm_tot), min(tb_wnm_tot)

		return tb_cont, tb_wnm_tot, tb_cnm_tot, tb_tot, exp_tausum
		
	## Multiple (N) Gaussians + offset. ##
	 #
	 # params list  v    VLSR
	 # params float zr   estimated constant zero offset of the data points.
	 # params list  h    the array of N estimated heights of the Gaussians.
	 # params list  v0   the array of N estimated centers of the Gaussians.
	 # params list  w    the array of N estimated halfwidths of the Gaussians.
	 #
	 # return 1-D-array  tf  The calculated points.
	 #
	 # version 01/2017 
	 # author Nguyen Van Hiep ##
	def gcurv(self, v, zr, h, v0, w):
		dp600 = np.float64(0.60056120)
		if(np.isscalar(v)):
			v  = np.array([v], dtype=np.float64)
		if(np.isscalar(h)):
			h  = np.array([h], dtype=np.float64)
		if(np.isscalar(v0)):
			v0 = np.array([v0], dtype=np.float64)
		if(np.isscalar(w)):
			w  = np.array([w], dtype=np.float64)

		#DETERMINE NR OF GAUSSIANS...
		ng = len(h)
		
		tf = 0.*v + zr
		for i in range(ng):
			if (w[i] > 0.):
				tf = tf + h[i]*np.exp(- ( (v-v0[i])/(dp600*w[i]))**2) # 0.6005612 - 1/e width

		return tf


	## Divide list by list ##
	 #
	 # params 
	 # params 
	 #
	 # return 
	 #
	 # version 08/2016 
	 # author Nguyen Van Hiep ##
	def ldiv(self,a,b):
		ret = []
		for i in range(len(a)):
			ret.append(a[i]/b[i])

		return ret

	## Take ABS of a list ##
	 #
	 # params 
	 # params 
	 #
	 # return 
	 #
	 # version 08/2016 
	 # author Nguyen Van Hiep ##
	def labs(self,a):
		ret = []
		for i in range(len(a)):
			ret.append(abs(a[i]))

		return ret

	## Multiply list by number ##
	 #
	 # params 
	 # params 
	 #
	 # return 
	 #
	 # version 08/2016 
	 # author Nguyen Van Hiep ##
	def listxnum(self,a,x):
		ret = []
		for i in range(len(a)):
			ret.append(a[i]*x)

		return ret

	## Multiply list by list ##
	 #
	 # params 
	 # params 
	 #
	 # return 
	 #
	 # version 08/2016 
	 # author Nguyen Van Hiep ##
	def lxl(self,a,b):
		ret = []
		for i in range(len(a)):
			ret.append(a[i]*b[i])

		return ret

	## Add list by list ##
	 #
	 # params 
	 # params 
	 #
	 # return 
	 #
	 # version 08/2016 
	 # author Nguyen Van Hiep ##
	def laddl(self,a,b):
		ret = []
		for i in range(len(a)):
			ret.append(a[i] + b[i])

		return ret


	##  Fit Abs spectrum: tau ##
	 #
	 # params...
	 # return... 	
	 #
	 # version 01/2017 
	 # author Nguyen Van Hiep ##
	def abfit(self,look, xdataa, tdataa, xindxrange, \
			zro0,   hgt0,   cen0,   wid0,\
			zro0yn, hgt0yn, cen0yn, wid0yn):

		dp600    = np.float64(0.60056120)
		xdataa   = xdataa.astype(np.float64)
		tdataa   = tdataa.astype(np.float64)

		nr_of_ns = len(xindxrange)/2
		datasize = 0L
		for nnr in range(nr_of_ns):
			datasize = datasize + xindxrange[2*nnr+1]-xindxrange[2*nnr]+1L

		xdata = np.zeros(datasize, dtype=np.float64)
		tdata = np.zeros(datasize, dtype=np.float64)

		## Data ##
		dtsiz = 0L
		for nnr in range(nr_of_ns):
			dtsiz1              = dtsiz + xindxrange[2*nnr+1]-xindxrange[2*nnr] +1L
			xdata[dtsiz:dtsiz1] = xdataa[xindxrange[2*nnr]:xindxrange[2*nnr+1]+1]
			tdata[dtsiz:dtsiz1] = tdataa[xindxrange[2*nnr]:xindxrange[2*nnr+1]+1]
			dtsiz               = dtsiz1

		# AX1 IS THE PERCENTAGE OF CHANGE THAT WE ALLOW# 1% IS THE DEFAULT...
		ax1 = np.float64(0.01)

		# HALFASSED IS THE MULTIPLIER FOR THE CORRECTIONS IN NONLI!=AR REGIME.
		# if (halfasseduse == No!=):
		halfassed = np.float64(0.5)

		#if (nloopmax == No!=):
		nloopmax = 350

		#A NONZERO PROBLEM INDICATES A PROBLEM...
		problem = 0

		#DFSTOP IS THE MAXIMUM WIDTH WE ALLOW, = 80% of the total window...
		dfstop = 0.8*abs(xdata[datasize-1]-xdata[0])

		# THE OUTPUT GAUSSIAN PARAMETERS# SCALE WID FROM FWHM TO 1/E...
		# THESE ARE THE SAME AS THE PARAMETERS THAT ARE ITERATED.
		zro1   = np.float64(zro0) # Scalar
		hgt1   = np.asarray(hgt0, dtype=np.float64)
		cen1   = np.asarray(cen0, dtype=np.float64)
		wid1   = dp600*np.asarray(wid0, dtype=np.float64)

		nloop     = 0
		nloopn    = 0

		# get Number OF GAUSSIANS TO FIT...
		ngaussians = len(hgt0)

		## get NR OF PARAMETERS TO FIT...
		nparams = self.lsum(zro0yn) + self.lsum(hgt0yn) + self.lsum(cen0yn) + self.lsum(wid0yn)


		# TOTAL NR OF PARAMGERS THAT WOULD BE FIT IF ALL YesNo'S = 1...
		nparams_max = 1 + 3*ngaussians

		## EQUATION-OF-CONDITION ARRAY, S AND ITS COUNTERPART SFULL...
		s          = np.zeros((datasize,nparams), dtype=np.float64)
		sfull      = np.zeros((datasize,nparams_max), dtype=np.float64)
		afull      = np.zeros(nparams_max, dtype=np.float64)
		sfull_to_s = [0]*nparams
		s_to_sfull = [0]*nparams_max
		sigarraya  = np.zeros(nparams_max, dtype=np.float64)

		# RELATIONSHIP BETWEEN COLS IN S AND SFULL...
		scol     = 0
		sfullcol = 0

		if (zro0yn != 0):
			s_to_sfull[0]    = int(scol)
			sfull_to_s[scol] = 0
			scol             = scol + 1

		for ng in range(ngaussians):
			if (hgt0yn[ng] != 0):
				s_to_sfull[3*ng+1] = int(scol)
				sfull_to_s[scol]   = 3*ng + 1
				scol               = scol + 1

			if (cen0yn[ng] != 0):
				s_to_sfull[3*ng+2] = int(scol)
				sfull_to_s[scol]   = 3*ng + 2
				scol               = scol + 1

			if (wid0yn[ng] != 0):
				s_to_sfull[3*ng+3] = int(scol)
				sfull_to_s[scol]   = 3*ng + 3
				scol               = scol + 1

		### BAT DAU VONG LAP ###
		redoit = 1
		while (redoit == 1):
			nloop  = nloop  + 1
			nloopn = nloopn + 1

			## FIRST DEFINE SFULL...
			sum_of_gaussians = self.gcurv(xdata, zro1, hgt1, cen1, np.array(wid1)/dp600)
			t_predicted      = np.exp(-sum_of_gaussians)
			expfactor        = -t_predicted
			sfull[:,0]       = expfactor 			## THE CONSTANT

			for ng in range(ngaussians):
				xdel                 = (xdata - cen1[ng])/wid1[ng]
				edel                 = np.exp(-xdel*xdel)
				sum1                 = edel
				sum2                 = edel*xdel
				sum3                 = sum2*xdel
				sum6                 = 2.0*hgt1[ng]/wid1[ng]
				sfull[ :, (3*ng+1) ] = expfactor*sum1          ## HGT
				sfull[ :, (3*ng+2) ] = expfactor*sum2*sum6     ## CNTR
				sfull[ :, (3*ng+3) ] = expfactor*sum3*sum6     ## WIDTH

			s = sfull[:,sfull_to_s]

			## CREATE AND SOLVE THE NORMAL EQUATION MATRICES...
			t   = tdata-t_predicted
			ss  = np.dot(np.transpose(s),s)
			st  = np.dot(np.transpose(s), np.transpose(t))
			ssi = np.linalg.inv(ss)
			a   = np.dot(ssi,st)
			afull[sfull_to_s] = a

			## CHECK THE DERIVED CNM PARAMETERS...
			## THE AMPLITUDES...
			delt  = afull[ [x+1 for x in (x*3 for x in list(range(ngaussians)))]  ]
			adelt = [abs(x) for x in delt]
			hgt1  = [abs(x) for x in hgt1]

			for i in range(len(adelt)):
				if(0.2*hgt1[i] < adelt[i]):
					adelt[i] = 0.2*hgt1[i]

			delthgt = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					delthgt.append(-adelt[i])
				else:
					delthgt.append(adelt[i])

			## CENTERS
			delt  = afull[ [x+2 for x in (x*3 for x in list(range(ngaussians)))]  ]
			adelt = [abs(x) for x in delt]
			wid1  = [abs(x) for x in wid1]

			for i in range(len(adelt)):
				if(0.2*wid1[i] < adelt[i]):
					adelt[i] = 0.2*wid1[i]

			deltcen = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					deltcen.append(-adelt[i])
				else:
					deltcen.append(adelt[i])

			## WIDTHS
			delt  = afull[ [x+3 for x in (x*3 for x in list(range(ngaussians)))]  ]
			adelt = [abs(x) for x in delt]
			wid1  = [abs(x) for x in wid1]

			for i in range(len(adelt)):
				if(0.2*wid1[i] < adelt[i]):
					adelt[i] = 0.2*wid1[i]

			deltwid = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					deltwid.append(-adelt[i])
				else:
					deltwid.append(adelt[i])

			## CHECK FOR CONVERGENCE AND IF CNM PARAMETERS are REASONABLE ##
			hgtf   = self.labs(self.ldiv(delthgt,hgt1))
			cenf   = self.labs(self.ldiv(deltcen,wid1))
			widf   = self.labs(self.ldiv(deltwid,wid1))

			redoit = 0
			if (max(hgtf) > ax1):
				redoit = 1
			if (max(cenf) > ax1):
				redoit = 1
			if (max(widf) > ax1):
				redoit = 1

			# INCREMENT THE PARAMETERS...
			if(redoit == 0):
				halfassed = 1.0

			zro1   = zro1 + halfassed * afull[0]
			hgt1   = self.laddl(hgt1, self.listxnum(delthgt,halfassed))
			cen1   = self.laddl(cen1, self.listxnum(deltcen,halfassed))
			wid1   = self.laddl(wid1, self.listxnum(deltwid,halfassed))

			# CHECK TO SEE IF WIDTH IS TOO BIG..but ignore if these params are fixed.
			if ((max(self.lxl(wid0yn,wid1)) > dfstop)  or \
				(min(self.lxl(wid0yn,wid1)) < 0.)  ):
			    problem = -1
			    break

			if (nloop >= nloopmax-1):
				problem = -2
				break
		
		## CONVERT THE 1/E WIDTHS TO HALFWIDTHS...
		wid1 = self.listxnum(wid1,1./dp600)

		## DERIVE THE FITTED POINTS, RESIDUALS, THE ERRORS IN DERIVED COEFFICIENTS...
		## NOTE THAT THE WIDTHS HAVE BEEN CONVERTED TO HALFWIDTHS HERE, SO THE
		## 0.6 FACTORS ARE NOT REQUIRED...
		sum_of_gaussians = self.gcurv(xdata, zro1, hgt1, cen1, np.array(wid1)/dp600)
		t_predicted      = np.exp(-sum_of_gaussians)
		        
		resid  = tdata - t_predicted
		resid2 = resid**2
		ressum = resid2.sum()
		sigsq  = ressum/(datasize - nparams)
		sigma  = sigsq**0.5

		plt.plot(xdata, t_predicted)
		plt.show()

		print nloop
		print resid
		for i in range(256):
			print i, resid[i], t_predicted[i]

		ltemp  = list(range(nparams))
		ltemp  = [x*(nparams+1) for x in ltemp]
		ssi_temp = ssi.ravel()
		sigarray = sigsq*ssi_temp[ltemp]

		countsqrt = 0
		indxsqrt  = []
		jj        = 0
		for x in np.nditer(sigarray):
			if (x<0.):
				countsqrt = countsqrt + 1
				indxsqrt.append(jj)
			jj = jj + 1

		sigarray = np.sqrt( abs(sigarray))

		## TEST FOR NEG SQRTS...
		if (countsqrt != 0):
			sigarray[indxsqrt] = -sigarray[indxsqrt]
			problem = -3

		## TEST FOR INFINITIES, ETC...
		countbad = 0
		indxbad  = []
		kk       = 0
		for x in np.nditer(a):
			if (np.isfinite(x) == False):
				countbad = countbad + 1
				indxbad.append(kk)

			kk = kk + 1

		if (countbad != 0):
			problem = -4

		sigarraya[sfull_to_s] = sigarray
		sigzro1               = sigarraya[0]
		temp_list             = [x*3 for x in list(range(ngaussians))]
		sighgt1               = sigarraya[ [x+1 for x in temp_list] ]
		sigcen1               = sigarraya[ [x+2 for x in temp_list] ]
		sigwid1               = self.listxnum(sigarraya[ [x+3 for x in temp_list] ], 1./dp600)

		## DERIVE THE NORMALIZED COVARIANCE ARRAY
		temp_list = [x*(nparams+1) for x in list(range(nparams))]
		ssi_temp  = ssi.ravel()
		doug      = ssi_temp[temp_list]
		doug_temp = doug[np.newaxis]
		doug_t    = doug[np.newaxis].T
		doug      = np.dot(doug_t,doug_temp)
		cov       = ssi/np.sqrt(doug)

		sum_of_gaussians = self.gcurv(xdata, zro1, hgt1, cen1, wid1)
		tfita            = np.exp(-sum_of_gaussians)

		# ## Absoprtion line
		# plt.plot(xdata,tdata, 'b-', linewidth=1, label='data')
		# plt.plot(xdata,tfita, 'r-', linewidth=1, label='fit')
		# plt.title('3C98', fontsize=30)
		# plt.ylabel('$T_{b} [K]$', fontsize=35)
		# plt.xlabel('$V_{lsr} (km/s)$', fontsize=35)
		# # plt.xlim(0.0, 2.0)
		# # plt.xlim(-1.0, 6.0)
		# plt.grid(True)
		# plt.tick_params(axis='x', labelsize=18)
		# plt.tick_params(axis='y', labelsize=15)

		# # plt.text(0.0, 3.2, r'$f = [0.32\pm0.06]\cdot log_{10}(N^*_{HI}/10^{20}) + [0.81\pm0.05]$, Lee et al.', color='blue', fontsize=20)
		# plt.legend(loc='upper left', fontsize=18)
		# plt.show()
		        
		return tfita, sigma, \
				zro1, hgt1, cen1, wid1,\
				sigzro1, sighgt1, sigcen1, sigwid1,\
				cov, problem,\
				nparams

	##  ##
	 #
	 # params 
	 # return 	
	 #
	 # version 08/2016 
	 # author Nguyen Van Hiep ##
	def fit(self,look, xdataa, tdataa, xindxrange, \
			zrocnm, hgtcnm, cencnm, widcnm, tspincnm, ordercnm, \
			zrocnmyn, hgtcnmyn, cencnmyn, widcnmyn, tspincnmyn, \
			zrownm, hgtwnm, cenwnm, widwnm, fwnm, \
			zrownmyn, hgtwnmyn, cenwnmyn, widwnmyn, fwnmyn, halfasseduse=0.5):

		dp600    = np.float64(0.60056120)
		xdataa   = xdataa.astype(np.float64)
		tdataa   = tdataa.astype(np.float64)

		nr_of_ns = len(xindxrange)/2
		datasize = 0L
		for nnr in range(nr_of_ns):
			datasize = datasize + xindxrange[2*nnr+1]-xindxrange[2*nnr]+1L

		xdata = np.zeros(datasize, dtype=np.float64)
		tdata = np.zeros(datasize, dtype=np.float64)

		## Data ##
		dtsiz = 0L
		for nnr in range(nr_of_ns):
			dtsiz1              = dtsiz + xindxrange[2*nnr+1]-xindxrange[2*nnr] +1L
			xdata[dtsiz:dtsiz1] = xdataa[xindxrange[2*nnr]:xindxrange[2*nnr+1]+1]
			tdata[dtsiz:dtsiz1] = tdataa[xindxrange[2*nnr]:xindxrange[2*nnr+1]+1]
			dtsiz               = dtsiz1

		# AX1 IS THE PERCENTAGE OF CHANGE THAT WE ALLOW# 1% IS THE DEFAULT...
		ax1 = 0.01
		ax1 = np.float64(0.003)

		# HALFASSED IS THE MULTIPLIER FOR THE CORRECTIONS IN NONLI!=AR REGIME.
		# if (halfasseduse == No!=):
		halfassed = np.float64(0.5)

		#if (nloopmax == No!=):
		nloopmax = 350

		#A NONZERO PROBLEM INDICATES A PROBLEM...
		problem = 0

		#DFSTOP IS THE MAXIMUM WIDTH WE ALLOW, = 80% of the total window...
		dfstop = 0.8*abs(xdata[datasize-1]-xdata[0])

		# THE OUTPUT GAUSSIAN PARAMETERS# SCALE WID FROM FWHM TO 1/E...
		# THESE ARE THE SAME AS THE PARAMETERS THAT ARE ITERATED.
		zrocnm1   = np.float64(zrocnm) # Scalar
		hgtcnm1   = np.asarray(hgtcnm, dtype=np.float64)
		cencnm1   = np.asarray(cencnm, dtype=np.float64)
		widcnm1   = dp600*np.asarray(widcnm, dtype=np.float64)
		tspincnm1 = np.asarray(tspincnm, dtype=np.float64)

		zrownm1   = np.float64(zrownm) # Scalar 
		hgtwnm1   = np.asarray(hgtwnm, dtype=np.float64)
		cenwnm1   = np.asarray(cenwnm, dtype=np.float64)
		widwnm1   = dp600*np.asarray(widwnm, dtype=np.float64)
		fwnm1     = np.asarray(fwnm, dtype=np.float64)

		nloop     = 0
		nloopn    = 0

		# get Number OF CNM AND WNM GAUSSIANS TO FIT...
		ngaussians_cnm = len(hgtcnm)
		ngaussians_wnm = len(hgtwnm)

		## get NR OF PARAMETERS TO FIT...
		nparams = self.lsum(zrocnmyn) + self.lsum(hgtcnmyn) + self.lsum(cencnmyn) + self.lsum(widcnmyn) + self.lsum(tspincnmyn) + \
			      self.lsum(zrownmyn) + self.lsum(hgtwnmyn) + self.lsum(cenwnmyn) + self.lsum(widwnmyn) + self.lsum(fwnmyn) 


		# TOTAL NR OF PARAMGERS THAT WOULD BE FIT IF ALL YesNo'S = 1...
		nparams_max = 2 + 4*(ngaussians_cnm + ngaussians_wnm)

		# print 'N params',nparams, nparams_max

		## EQUATION-OF-CONDITION ARRAY, S AND ITS COUNTERPART SFULL...
		s          = np.zeros((datasize,nparams), dtype=np.float64)
		sfull      = np.zeros((datasize,nparams_max), dtype=np.float64)
		afull      = np.zeros(nparams_max, dtype=np.float64)
		sfull_to_s = [0]*nparams
		sigarraya  = np.zeros(nparams_max, dtype=np.float64)

		# RELATIONSHIP BETWEEN COLS IN S AND SFULL...
		# BEGIN WITH THE CNM PARAMETERS...
		scol     = 0
		sfullcol = 0

		if (zrocnmyn != 0):
			sfull_to_s[scol] = int(sfullcol)
			scol             = scol + 1 

		sfullcol = sfullcol + 1

		for ng in range(ngaussians_cnm):
			if (hgtcnmyn[ng] != 0):
				sfull_to_s[scol] = int(sfullcol)
				scol             = scol + 1 

			sfullcol = sfullcol + 1

			if (cencnmyn[ng] != 0):
				sfull_to_s[scol] = int(sfullcol)
				scol             = scol + 1 

			sfullcol = sfullcol + 1

			if (widcnmyn[ng] != 0):
				sfull_to_s[scol] = int(sfullcol)
				scol             = scol + 1 

			sfullcol = sfullcol + 1

			if (tspincnmyn[ng] != 0):
				sfull_to_s[scol] = int(sfullcol)
				scol             = scol + 1 

			sfullcol = sfullcol + 1

		## THEN THE WNM PARAMETERS...
		if (zrownmyn != 0):
			sfull_to_s[scol] = int(sfullcol)
			scol             = scol + 1 
		
		sfullcol = sfullcol + 1

		if(type(hgtwnmyn) is int):
			hgtwnmyn = [hgtwnmyn]
			cenwnmyn = [cenwnmyn]
			widwnmyn = [widwnmyn]
			fwnmyn   = [fwnmyn]

		for ng in range(ngaussians_wnm): 
			if (hgtwnmyn[ng] != 0):
				sfull_to_s[scol] = int(sfullcol)
				scol             = scol + 1 
			
			sfullcol = sfullcol + 1

			if (cenwnmyn[ng] != 0):
				sfull_to_s[scol] = int(sfullcol)
				scol             = scol + 1 
			
			sfullcol = sfullcol + 1

			if (widwnmyn[ng] != 0):
				sfull_to_s[scol] = int(sfullcol)
				scol             = scol + 1 
			
			sfullcol = sfullcol + 1

			if (fwnmyn[ng] != 0):
				sfull_to_s[scol] = int(sfullcol)
				scol             = scol + 1 
			
			sfullcol = sfullcol + 1

		### BAT DAU VONG LAP ###
		redoit = 1
		while (redoit == 1):
		# for j in range(1):
			nloop  = nloop  + 1
			nloopn = nloopn + 1

			sfullcol = 0

			# EVALUATE CONSTANT DERIVATIVE FOR CNM:
			xdel    = np.float64(0.0000025)
			zrocnm1 = zrocnm1 + xdel
			tb_cont, tb_wnm_tot, tb_cnm_tot, tb_totplus, exp_tau_sum = self.tb_exp(xdata, \
				zrocnm1, hgtcnm1, cencnm1, self.listxnum(widcnm1,1./dp600), tspincnm1, ordercnm, \
			        zrownm1, hgtwnm1, cenwnm1, self.listxnum(widwnm1,1./dp600), fwnm1)
			        

			zrocnm1 = zrocnm1 - 2.*xdel
			tb_cont, tb_wnm_tot, tb_cnm_tot, tb_totminus, exp_tau_sum = self.tb_exp(xdata, \
				zrocnm1, hgtcnm1, cencnm1, self.listxnum(widcnm1,1./dp600), tspincnm1, ordercnm, \
			        zrownm1, hgtwnm1, cenwnm1, self.listxnum(widwnm1,1./dp600), fwnm1)
			        

			zrocnm1   = zrocnm1 + xdel
			zrocnmder = (tb_totplus - tb_totminus)/(2.*xdel)

			sfull[:,sfullcol] = zrocnmder				#THE CONSTANT
			sfullcol          = sfullcol + 1

			# WORK THROUGH CNM GAUSSIANS...
			for ng in range(ngaussians_cnm):
				#EVALUATE HGT DERIVATIVE:
				xdel = np.float64(0.0000025)
				hgtcnm1[ng] = hgtcnm1[ng] + xdel
				tb_cont, tb_wnm_tot, tb_cnm_tot, tb_totplus, exp_tau_sum = self.tb_exp(xdata, \
					zrocnm1, hgtcnm1, cencnm1, self.listxnum(widcnm1,1./dp600), tspincnm1, ordercnm, \
				        zrownm1, hgtwnm1, cenwnm1, self.listxnum(widwnm1,1./dp600), fwnm1)
				        
				hgtcnm1[ng] = hgtcnm1[ng] -2.* xdel
				tb_cont, tb_wnm_tot, tb_cnm_tot, tb_totminus, exp_tau_sum = self.tb_exp(xdata, \
					zrocnm1, hgtcnm1, cencnm1, self.listxnum(widcnm1,1./dp600), tspincnm1, ordercnm, \
				        zrownm1, hgtwnm1, cenwnm1, self.listxnum(widwnm1,1./dp600), fwnm1)
				        
				hgtcnm1[ng] = hgtcnm1[ng] + xdel
				hgtder      = (tb_totplus - tb_totminus)/(2.*xdel)

				## EVALUATE CEN DERIVATIVE:
				xdel        = 0.0000025*widcnm1[ng]
				cencnm1[ng] = cencnm1[ng] + xdel
				tb_cont, tb_wnm_tot, tb_cnm_tot, tb_totplus, exp_tau_sum = self.tb_exp(xdata, \
					zrocnm1, hgtcnm1, cencnm1, self.listxnum(widcnm1,1./dp600), tspincnm1, ordercnm, \
				        zrownm1, hgtwnm1, cenwnm1, self.listxnum(widwnm1,1./dp600), fwnm1)
				        
				cencnm1[ng] = cencnm1[ng] - 2.*xdel
				tb_cont, tb_wnm_tot, tb_cnm_tot, tb_totminus, exp_tau_sum = self.tb_exp(xdata, \
					zrocnm1, hgtcnm1, cencnm1, self.listxnum(widcnm1,1./dp600), tspincnm1, ordercnm, \
				        zrownm1, hgtwnm1, cenwnm1, self.listxnum(widwnm1,1./dp600), fwnm1)
				        
				cencnm1[ng] = cencnm1[ng] + xdel
				cender      = (tb_totplus - tb_totminus)/(2.*xdel)

				## EVALUATE WID DERIVATIVE:
				xdel        = np.float64(0.0000025)*widcnm1[ng]
				widcnm1[ng] = widcnm1[ng] + xdel
				tb_cont, tb_wnm_tot, tb_cnm_tot, tb_totplus, exp_tau_sum = self.tb_exp(xdata, \
					zrocnm1, hgtcnm1, cencnm1, self.listxnum(widcnm1,1./dp600), tspincnm1, ordercnm, \
				        zrownm1, hgtwnm1, cenwnm1, self.listxnum(widwnm1,1./dp600), fwnm1)
				        
				widcnm1[ng] = widcnm1[ng] - 2.*xdel
				tb_cont, tb_wnm_tot, tb_cnm_tot, tb_totminus, exp_tau_sum = self.tb_exp(xdata, \
					zrocnm1, hgtcnm1, cencnm1, self.listxnum(widcnm1,1./dp600), tspincnm1, ordercnm, \
				        zrownm1, hgtwnm1, cenwnm1, self.listxnum(widwnm1,1./dp600), fwnm1)
				        
				widcnm1[ng] = widcnm1[ng] + xdel
				widder      = (tb_totplus - tb_totminus)/(2.*xdel)

				## EVALUATE TSPIN DERIVATIVE:
				xdel          = np.float64(0.0000025)*tspincnm1[ng]
				tspincnm1[ng] = tspincnm1[ng] + xdel
				tb_cont, tb_wnm_tot, tb_cnm_tot, tb_totplus, exp_tau_sum = self.tb_exp(xdata, \
					zrocnm1, hgtcnm1, cencnm1, self.listxnum(widcnm1,1./dp600), tspincnm1, ordercnm, \
				        zrownm1, hgtwnm1, cenwnm1, self.listxnum(widwnm1,1./dp600), fwnm1)
				        
				tspincnm1[ng] = tspincnm1[ng] -2.*xdel
				tb_cont, tb_wnm_tot, tb_cnm_tot, tb_totminus, exp_tau_sum = self.tb_exp(xdata, \
					zrocnm1, hgtcnm1, cencnm1, self.listxnum(widcnm1,1./dp600), tspincnm1, ordercnm, \
				        zrownm1, hgtwnm1, cenwnm1, self.listxnum(widwnm1,1./dp600), fwnm1)
				        
				tspincnm1[ng] = tspincnm1[ng] + xdel
				tspinder      = (tb_totplus - tb_totminus)/(2.*xdel)
				sfull[:, sfullcol] = hgtder     #HGT-*/41-*/7
				sfullcol           = sfullcol + 1

				sfull[:, sfullcol] = cender     #CNTR
				sfullcol           = sfullcol + 1

				sfull[:, sfullcol] = widder     #WIDTH
				sfullcol           = sfullcol + 1

				sfull[:, sfullcol] = tspinder   #TSPIN
				sfullcol           = sfullcol + 1

			## EVALUATE CONSTANT DERIVATIVE FOR WNM:
			xdel = np.float64(0.0000025)
			zrownm1 = zrownm1 + xdel
			tb_cont, tb_wnm_tot, tb_cnm_tot, tb_totplus, exp_tau_sum = self.tb_exp(xdata, \
				zrocnm1, hgtcnm1, cencnm1, self.listxnum(widcnm1,1./dp600), tspincnm1, ordercnm, \
			        zrownm1, hgtwnm1, cenwnm1, self.listxnum(widwnm1,1./dp600), fwnm1)
			        
			zrownm1 = zrownm1 - 2.*xdel
			tb_cont, tb_wnm_tot, tb_cnm_tot, tb_totminus, exp_tau_sum = self.tb_exp(xdata, \
				zrocnm1, hgtcnm1, cencnm1, self.listxnum(widcnm1,1./dp600), tspincnm1, ordercnm, \
			        zrownm1, hgtwnm1, cenwnm1, self.listxnum(widwnm1,1./dp600), fwnm1)
			        
			zrownm1   = zrownm1 + xdel
			zrownmder = (tb_totplus - tb_totminus)/(2.*xdel)

			sfull[:, sfullcol] = zrownmder				#THE CONSTANT
			sfullcol           = sfullcol + 1

			# WORK THROUGH wnm GAUSSIANS...
			for ng in range(ngaussians_wnm):
				## EVALUATE HGT DERIVATIVE:
				xdel = np.float64(0.0000025)
				hgtwnm1[ng] = hgtwnm1[ng] + xdel
				tb_cont, tb_wnm_tot, tb_cnm_tot, tb_totplus, exp_tau_sum = self.tb_exp(xdata, \
					zrocnm1, hgtcnm1, cencnm1, self.listxnum(widcnm1,1./dp600), tspincnm1, ordercnm, \
				        zrownm1, hgtwnm1, cenwnm1, self.listxnum(widwnm1,1./dp600), fwnm1)
				        
				hgtwnm1[ng] = hgtwnm1[ng] -2.* xdel
				tb_cont, tb_wnm_tot, tb_cnm_tot, tb_totminus, exp_tau_sum = self.tb_exp(xdata, \
					zrocnm1, hgtcnm1, cencnm1, self.listxnum(widcnm1,1./dp600), tspincnm1, ordercnm, \
				        zrownm1, hgtwnm1, cenwnm1, self.listxnum(widwnm1,1./dp600), fwnm1)
				        
				hgtwnm1[ng] = hgtwnm1[ng] + xdel
				hgtder      = (tb_totplus - tb_totminus)/(2.*xdel)

				# EVALUATE CEN DERIVATIVE:
				xdel = 0.0000025*widwnm1[ng]
				cenwnm1[ng] = cenwnm1[ng] + xdel
				tb_cont, tb_wnm_tot, tb_cnm_tot, tb_totplus, exp_tau_sum = self.tb_exp(xdata, \
					zrocnm1, hgtcnm1, cencnm1, self.listxnum(widcnm1,1./dp600), tspincnm1, ordercnm, \
				        zrownm1, hgtwnm1, cenwnm1, self.listxnum(widwnm1,1./dp600), fwnm1)
				        
				cenwnm1[ng] = cenwnm1[ng] - 2.*xdel
				tb_cont, tb_wnm_tot, tb_cnm_tot, tb_totminus, exp_tau_sum = self.tb_exp(xdata, \
					zrocnm1, hgtcnm1, cencnm1, self.listxnum(widcnm1,1./dp600), tspincnm1, ordercnm, \
				        zrownm1, hgtwnm1, cenwnm1, self.listxnum(widwnm1,1./dp600), fwnm1)
				        
				cenwnm1[ng] = cenwnm1[ng] + xdel
				cender      = (tb_totplus - tb_totminus)/(2.*xdel)

				# EVALUATE WID DERIVATIVE:
				xdel = np.float64(0.0000025)*widwnm1[ng]
				widwnm1[ng] = widwnm1[ng] + xdel
				tb_cont, tb_wnm_tot, tb_cnm_tot, tb_totplus, exp_tau_sum = self.tb_exp(xdata, \
					zrocnm1, hgtcnm1, cencnm1, self.listxnum(widcnm1,1./dp600), tspincnm1, ordercnm, \
				        zrownm1, hgtwnm1, cenwnm1, self.listxnum(widwnm1,1./dp600), fwnm1)
				        
				widwnm1[ng] = widwnm1[ng] - 2.*xdel
				tb_cont, tb_wnm_tot, tb_cnm_tot, tb_totminus, exp_tau_sum = self.tb_exp(xdata, \
					zrocnm1, hgtcnm1, cencnm1, self.listxnum(widcnm1,1./dp600), tspincnm1, ordercnm, \
				        zrownm1, hgtwnm1, cenwnm1, self.listxnum(widwnm1,1./dp600), fwnm1)
				        
				widwnm1[ng] = widwnm1[ng] + xdel
				widder = (tb_totplus - tb_totminus)/(2.*xdel)

				#EVALUATE F DERIVATIVE:
				xdel = np.float64(0.0000025) * fwnm1[ng]
				fwnm1[ng] = fwnm1[ng] + xdel
				tb_cont, tb_wnm_tot, tb_cnm_tot, tb_totplus, exp_tau_sum = self.tb_exp(xdata, \
					zrocnm1, hgtcnm1, cencnm1, self.listxnum(widcnm1,1./dp600), tspincnm1, ordercnm, \
				        zrownm1, hgtwnm1, cenwnm1, self.listxnum(widwnm1,1./dp600), fwnm1)
				        
				fwnm1[ng] = fwnm1[ng] -2.*xdel
				tb_cont, tb_wnm_tot, tb_cnm_tot, tb_totminus, exp_tau_sum = self.tb_exp(xdata, \
					zrocnm1, hgtcnm1, cencnm1, self.listxnum(widcnm1,1./dp600), tspincnm1, ordercnm, \
				        zrownm1, hgtwnm1, cenwnm1, self.listxnum(widwnm1,1./dp600), fwnm1)
				        
				fwnm1[ng] = fwnm1[ng] + xdel
				fder = (tb_totplus - tb_totminus)/(2.*xdel)
				sfull[ :, sfullcol] = hgtder #HGT
				sfullcol            = sfullcol + 1
				sfull[ :, sfullcol] = cender #CNTR
				sfullcol            = sfullcol + 1
				sfull[ :, sfullcol] = widder #WIDTH
				sfullcol            = sfullcol + 1
				sfull[ :, sfullcol] = fder      #f
				sfullcol            = sfullcol + 1

			s = sfull[:,sfull_to_s]

			# print 'S shape',s.shape

			# CALCULATE T_PREDICTED...
			tb_cont, tb_wnm_tot, tb_cnm_tot, t_predicted, exp_tau_sum = self.tb_exp(xdata, \
				zrocnm1, hgtcnm1, cencnm1, self.listxnum(widcnm1,1./dp600), tspincnm1, ordercnm, \
			        zrownm1, hgtwnm1, cenwnm1, self.listxnum(widwnm1,1./dp600), fwnm1)
			        
			## CREATE AND SOLVE THE NORMAL EQUATION MATRICES...
			t   = tdata-t_predicted
			ss  = np.dot(np.transpose(s),s)
			st  = np.dot(np.transpose(s), np.transpose(t))
			ssi = np.linalg.inv(ss)
			a   = np.dot(ssi,st)
			afull[sfull_to_s] = a

			## CHECK THE DERIVED CNM PARAMETERS...
			## THE CNM AMPLITUDES...
			delt    = afull[ [x+1 for x in (x*4 for x in list(range(ngaussians_cnm)))]  ]
			adelt   = [abs(x) for x in delt]
			hgtcnm1 = [abs(x) for x in hgtcnm1]

			for i in range(len(adelt)):
				if(0.2*hgtcnm1[i] < adelt[i]):
					adelt[i] = 0.2*hgtcnm1[i]

			delthgtcnm = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					delthgtcnm.append(-adelt[i])
				else:
					delthgtcnm.append(adelt[i])

			## CNM CENTERS
			delt    = afull[ [x+2 for x in (x*4 for x in list(range(ngaussians_cnm)))]  ]
			adelt   = [abs(x) for x in delt]
			widcnm1 = [abs(x) for x in widcnm1]

			for i in range(len(adelt)):
				if(0.2*widcnm1[i] < adelt[i]):
					adelt[i] = 0.2*widcnm1[i]

			deltcencnm = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					deltcencnm.append(-adelt[i])
				else:
					deltcencnm.append(adelt[i])

			## CNM WIDTHS
			delt    = afull[ [x+3 for x in (x*4 for x in list(range(ngaussians_cnm)))]  ]
			adelt   = [abs(x) for x in delt]
			widcnm1 = [abs(x) for x in widcnm1]

			for i in range(len(adelt)):
				if(0.2*widcnm1[i] < adelt[i]):
					adelt[i] = 0.2*widcnm1[i]

			deltwidcnm = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					deltwidcnm.append(-adelt[i])
				else:
					deltwidcnm.append(adelt[i])

			## CNM Tex
			delt      = afull[ [x+4 for x in (x*4 for x in list(range(ngaussians_cnm)))]  ]
			adelt     = [abs(x) for x in delt]
			tspincnm1 = [abs(x) for x in tspincnm1]

			for i in range(len(adelt)):
				if(0.2*tspincnm1[i] < adelt[i]):
					adelt[i] = 0.2*tspincnm1[i]

			delttspincnm = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					delttspincnm.append(-adelt[i])
				else:
					delttspincnm.append(adelt[i])

			## CHECK FOR CONVERGENCE AND IF CNM PARAMETERS are REASONABLE ##
			hgtf   = self.labs(self.ldiv(delthgtcnm,hgtcnm1))
			cenf   = self.labs(self.ldiv(deltcencnm,widcnm1))
			widf   = self.labs(self.ldiv(deltwidcnm,widcnm1))
			tspinf = self.labs(self.ldiv(delttspincnm,tspincnm1))

			redoit = 0
			if (max(hgtf) > ax1):
				redoit = 1
			if (max(cenf) > ax1):
				redoit = 1
			if (max(widf) > ax1):
				redoit = 1
			if (max(tspinf) > ax1):
				redoit = 1


			## CHECK THE DERIVED CNM PARAMETERS...
			## THE WNM AMPLITUDES...
			delt    = afull[ [x+4*ngaussians_cnm+2 for x in (x*4 for x in list(range(ngaussians_wnm)))]  ]
			adelt   = [abs(x) for x in delt]
			hgtwnm1 = [abs(x) for x in hgtwnm1]

			for i in range(len(adelt)):
				if(0.2*hgtwnm1[i] < adelt[i]):
					adelt[i] = 0.2*hgtwnm1[i]

			delthgtwnm = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					delthgtwnm.append(-adelt[i])
				else:
					delthgtwnm.append(adelt[i])

			## THE WNM CENTERS
			delt    = afull[ [x+4*ngaussians_cnm+3 for x in (x*4 for x in list(range(ngaussians_wnm)))]  ]
			adelt   = [abs(x) for x in delt]
			widwnm1 = [abs(x) for x in widwnm1]

			for i in range(len(adelt)):
				if(0.2*widwnm1[i] < adelt[i]):
					adelt[i] = 0.2*widwnm1[i]

			deltcenwnm = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					deltcenwnm.append(-adelt[i])
				else:
					deltcenwnm.append(adelt[i])

			## THE WNM WIDTHS
			delt    = afull[ [x+4*ngaussians_cnm+4 for x in (x*4 for x in list(range(ngaussians_wnm)))]  ]
			adelt   = [abs(x) for x in delt]
			widwnm1 = [abs(x) for x in widwnm1]

			for i in range(len(adelt)):
				if(0.2*widwnm1[i] < adelt[i]):
					adelt[i] = 0.2*widwnm1[i]

			deltwidwnm = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					deltwidwnm.append(-adelt[i])
				else:
					deltwidwnm.append(adelt[i])

			## THE WNM FRACTIONs
			delt    = afull[ [x+4*ngaussians_cnm+5 for x in (x*4 for x in list(range(ngaussians_wnm)))]  ]
			adelt   = [abs(x) for x in delt]
			fwnm1   = [abs(x) for x in fwnm1]

			for i in range(len(adelt)):
				if(0.2*fwnm1[i] < adelt[i]):
					adelt[i] = 0.2*fwnm1[i]

			deltfwnm = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					deltfwnm.append(-adelt[i])
				else:
					deltfwnm.append(adelt[i])

			## CHECK FOR CONVERGENCE AND IF WNM PARAMETERS are REASONABLE ##
			hgtwnmf   = self.labs(self.ldiv(delthgtwnm,hgtwnm1))
			cenwnmf   = self.labs(self.ldiv(deltcenwnm,widwnm1))
			widwnmf   = self.labs(self.ldiv(deltwidwnm,widwnm1))
			fwnmf     = self.labs(self.ldiv(deltfwnm,fwnm1))

			if (max(hgtwnmf) > ax1):
				redoit = 1
			if (max(cenwnmf) > ax1):
				redoit = 1
			if (max(widwnmf) > ax1):
				redoit = 1
			if (max(fwnmf) > ax1):
				redoit = 1

			# INCREMENT THE PARAMETERS...
			# halfassed = 0.5
			# halfassed = 0.4
			if(redoit == 0):
				halfassed = 1.0

			zrocnm1   = zrocnm1   + halfassed * afull[0]
			hgtcnm1   = self.laddl(hgtcnm1, self.listxnum(delthgtcnm,halfassed))
			cencnm1   = self.laddl(cencnm1, self.listxnum(deltcencnm,halfassed))
			widcnm1   = self.laddl(widcnm1, self.listxnum(deltwidcnm,halfassed))
			tspincnm1 = self.laddl(tspincnm1, self.listxnum(delttspincnm,halfassed))

			zrownm1 = zrownm1 + halfassed * afull[4*ngaussians_cnm + 1]
			hgtwnm1 = self.laddl(hgtwnm1, self.listxnum(delthgtwnm,halfassed))
			cenwnm1 = self.laddl(cenwnm1, self.listxnum(deltcenwnm,halfassed))
			widwnm1 = self.laddl(widwnm1, self.listxnum(deltwidwnm,halfassed))
			fwnm1   = self.laddl(fwnm1, self.listxnum(deltfwnm,halfassed))

			# CHECK TO SEE IF WIDTH IS TOO BIG..but ignore if these params are fixed.
			if ((max(self.lxl(widcnmyn,widcnm1)) > dfstop)  or \
				(min(self.lxl(widcnmyn,widcnm1)) < 0.)     or \
				(max(self.lxl(widwnmyn,widwnm1)) > dfstop) or \
				(min(self.lxl(widwnmyn,widwnm1)) < 0.) ):
			    problem = -1
			    break

		   	# nloop = 200
			if (nloop >= nloopmax-1):
				problem = -2
				break
		
		## CONVERT THE 1/E WIDTHS TO HALFWIDTHS...
		widcnm1 = self.listxnum(widcnm1,1./dp600)
		widwnm1 = self.listxnum(widwnm1,1./dp600)

		## DERIVE THE FITTED POINTS, RESIDUALS, THE ERRORS IN DERIVED COEFFICIENTS...
		## NOTE THAT THE WIDTHS HAVE BEEN CONVERTED TO HALFWIDTHS HERE, SO THE
		## 0.6 FACTORS ARE NOT REQUIRED...
		tb_cont, tb_wnm_tot, tb_cnm_tot, t_predicted, exp_tau_sum = self.tb_exp(xdata, \
			zrocnm1, hgtcnm1, cencnm1, widcnm1, tspincnm1, ordercnm, \
		        zrownm1, hgtwnm1, cenwnm1, widwnm1, fwnm1)
		        
		resid  = tdata - t_predicted
		resid2 = np.square(resid)
		sigsq  = resid2.sum()/(datasize - nparams)
		sigma  = sigsq**0.5

		ltemp  = list(range(nparams))
		ltemp  = [x*(nparams+1) for x in ltemp]
		ssi_temp = ssi.ravel()
		sigarray = sigsq*ssi_temp[ltemp]

		countsqrt = 0
		indxsqrt  = []
		jj        = 0
		for x in np.nditer(sigarray):
			if (x<0.):
				countsqrt = countsqrt + 1
				indxsqrt.append(jj)

			jj = jj + 1

		sigarray = np.sqrt( abs(sigarray))

		## TEST FOR NEG SQRTS...
		if (countsqrt != 0):
			sigarray[indxsqrt] = -sigarray[indxsqrt]
			problem = -3

		## TEST FOR INFINITIES, ETC...
		countbad = 0
		indxbad  = []
		kk       = 0
		for x in np.nditer(a):
			if (np.isfinite(x) == False):
				countbad = countbad + 1
				indxbad.append(kk)

			kk = kk + 1

		if (countbad != 0):
			problem = -4

		sigarraya[sfull_to_s] = sigarray
		sigzrocnm1            = sigarraya[0]
		temp_list             = [x*4 for x in list(range(ngaussians_cnm))]
		sighgtcnm1            = sigarraya[ [x+1 for x in temp_list] ]
		sigcencnm1            = sigarraya[ [x+2 for x in temp_list] ]
		sigwidcnm1            = self.listxnum(sigarraya[ [x+3 for x in temp_list] ], 1./dp600)
		sigtspincnm1          = sigarraya[ [x+4 for x in temp_list] ]

		temp_list             = [x*4 for x in list(range(ngaussians_wnm))]
		sigzrownm1            = sigarraya[ 4*ngaussians_cnm + 1]
		sighgtwnm1            = sigarraya[ [x+4*ngaussians_cnm+2 for x in temp_list] ]
		sigcenwnm1            = sigarraya[ [x+4*ngaussians_cnm+3 for x in temp_list] ]
		sigwidwnm1            = self.listxnum(sigarraya[ [x+4*ngaussians_cnm+4 for x in temp_list] ], 1./dp600)
		sigfwnm1              = sigarraya[ [x+4*ngaussians_cnm+5 for x in temp_list] ]


		## DERIVE THE NORMALIZED COVARIANCE ARRAY
		temp_list = [x*(nparams+1) for x in list(range(nparams))]
		ssi_temp  = ssi.ravel()
		doug      = ssi_temp[temp_list]
		doug_temp = doug[np.newaxis]
		doug_t    = doug[np.newaxis].T
		doug      = np.dot(doug_t,doug_temp)
		cov       = ssi/np.sqrt(doug)

		tb_cont, tb_wnm_tot, tb_cnm_tot, tfita, exp_tau_sum = self.tb_exp(xdata, \
			zrocnm1, hgtcnm1, cencnm1, widcnm1, tspincnm1, ordercnm, \
		        zrownm1, hgtwnm1, cenwnm1, widwnm1, fwnm1)

		# ## Absoprtion line
		# plt.plot(xdata,tdata, 'b-', linewidth=1, label='data')
		# plt.plot(xdata,tfita, 'r-', linewidth=1, label='fit')
		# plt.title('3C98', fontsize=30)
		# plt.ylabel('$T_{b} [K]$', fontsize=35)
		# plt.xlabel('$V_{lsr} (km/s)$', fontsize=35)
		# # plt.xlim(0.0, 2.0)
		# # plt.xlim(-1.0, 6.0)
		# plt.grid(True)
		# plt.tick_params(axis='x', labelsize=18)
		# plt.tick_params(axis='y', labelsize=15)

		# # plt.text(0.0, 3.2, r'$f = [0.32\pm0.06]\cdot log_{10}(N^*_{HI}/10^{20}) + [0.81\pm0.05]$, Lee et al.', color='blue', fontsize=20)
		# plt.legend(loc='upper left', fontsize=18)
		# plt.show()
		        
		return tfita, sigma, \
				zrocnm1, hgtcnm1, cencnm1, widcnm1, tspincnm1, \
				sigzrocnm1, sighgtcnm1, sigcencnm1, sigwidcnm1, sigtspincnm1, \
				zrownm1, hgtwnm1, cenwnm1, widwnm1, fwnm1, \
				sigzrownm1, sighgtwnm1, sigcenwnm1, sigwidwnm1, sigfwnm1, \
				cov, problem, nloop, \
		        tb_cont, tb_wnm_tot, tb_cnm_tot, \
		        exp_tau_sum, nloopmax, halfasseduse
		