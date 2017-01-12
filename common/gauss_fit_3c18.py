import sys
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

	## Model of Absorption line and Emission line ##
	 #
	 # params 
	 # params 
	 #
	 # return 
	 #
	 # version 08/2016 
	 # author Nguyen Van Hiep ##
	def tb_exp(self,x,bg,tau,v0,wid,tex,cont):
		bg      = bg
		amp     = list(tau) # copy
		cntr    = list(v0)  # copy
		width   = list(wid) # copy
		ts      = list(tex) # copy

		# CAL. THE OPACITY TAU
		lenx   = len(x) # sometimes = 2048
		ng     = len(amp)
		arrtau = np.zeros((lenx, ng),dtype=np.float64)

		for i in range(ng):
			arrtau[:, i] = self.gcurv(x,bg,amp[i],cntr[i],width[i])

		sumtau = arrtau.sum(1) #2048

		exp_sumtau = np.exp(-sumtau)
		tb_cont    = cont*exp_sumtau
		tb_peaks   = np.zeros(lenx)

		# BRIGHT TEMP OF EACH PEAK:
		for i in range(ng):
			temp        = np.reshape(arrtau[:, 0:i+1], (lenx, i+1))
			sumtau_i    = temp.sum(1)
			exp_tau_i   = np.exp(arrtau[:, i] - sumtau_i)
			tb_peaks    = tb_peaks + ts[i] * (1. - np.exp(-arrtau[:,i]) ) * exp_tau_i ## From Carl's paper

		tb_tot = tb_cont + tb_peaks # 2048

		return tb_tot
		
	## Calculate multiple (N) Gaussians + offset. ##
	 #
	 # params 
	 # params 
	 #
	 # return 
	 #
	 # version 08/2016 
	 # author Nguyen Van Hiep ##
	def gcurv(self, xdata, zro1, hgt1, cen1, widfit):
		tfit = 0.*xdata + zro1
		if (widfit > 0.):
			tfit = tfit + hgt1*np.exp(- ( (xdata-cen1)/(0.6005612*widfit))**2)

		return tfit


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

	##  ##
	 #
	 # params 
	 # return 	
	 #
	 # version 08/2016 
	 # author Nguyen Van Hiep ##
	def fit(self,vdata, tbdata, vrangeid, \
			zrolv, tau, v0, wid, tex, \
			zrolvyn, tauyn, v0yn, widyn, texyn, \
			tcont,tcontyn):

		dp600   = np.float64(0.60056120)
		vdata   = vdata.astype(np.float64)
		tbdata  = tbdata.astype(np.float64)

		nr_of_ns = len(vrangeid)/2
		datasize = 0L
		for nnr in range(nr_of_ns):
			datasize = datasize + vrangeid[2*nnr+1]-vrangeid[2*nnr]+1L

		xdata = np.zeros(datasize, dtype=np.float64)
		tdata = np.zeros(datasize, dtype=np.float64)

		## Data ##
		dtsiz = 0L
		for nnr in range(nr_of_ns):
			dtsiz1              = dtsiz + vrangeid[2*nnr+1]-vrangeid[2*nnr] +1L
			xdata[dtsiz:dtsiz1] = vdata[vrangeid[2*nnr]:vrangeid[2*nnr+1]+1]
			tdata[dtsiz:dtsiz1] = tbdata[vrangeid[2*nnr]:vrangeid[2*nnr+1]+1]

		# Criterion
		criterion = np.float64(0.001)

		# xfact 
		xfact    = np.float64(0.5)
		nloopmax = 210

		#DFSTOP IS THE MAXIMUM WIDTH ALLOWED
		dfstop = abs(xdata[datasize-1]-xdata[0])

		# THE OUTPUT GAUSSIAN PARAMETERS# SCALE WID FROM FWHM TO 1/E
		zrolvfit = np.float64(zrolv) # Scalar
		taufit   = np.asarray(tau, dtype=np.float64)
		v0fit    = np.asarray(v0, dtype=np.float64)
		widfit   = dp600*np.asarray(wid, dtype=np.float64)
		texfit   = np.asarray(tex, dtype=np.float64)
		tcontfit = np.float64(tcont) # Scalar 

		nloop    = 0

		# get Number OF GAUSSIANS TO FIT
		ngaussians = len(tau)

		## NR OF PARAMETERS
		nparams = self.lsum(zrolvyn) + self.lsum(tauyn) + self.lsum(v0yn) + self.lsum(widyn) + self.lsum(texyn) + \
			      self.lsum(tcontyn)

		# Max Number of Params
		nparams_max = 2 + 4*(ngaussians)

		## EQUATION-OF-CONDITION
		jacob        = np.zeros((datasize,nparams), dtype=np.float64)
		jacobfull    = np.zeros((datasize,nparams_max), dtype=np.float64)
		parfull      = np.zeros(nparams_max, dtype=np.float64)
		jfull_to_jcb = [0]*nparams
		erarray      = np.zeros(nparams_max, dtype=np.float64)

		# RELATIONSHIP BETWEEN COLS IN Jacobian AND jacobianfull
		jcol     = 0
		jfullcol = 0

		if (zrolvyn != 0):
			jfull_to_jcb[jcol] = int(jfullcol)
			jcol               = jcol + 1 

		jfullcol = jfullcol + 1

		for ng in range(ngaussians):
			if (tauyn[ng] != 0):
				jfull_to_jcb[jcol] = int(jfullcol)
				jcol               = jcol + 1 

			jfullcol = jfullcol + 1

			if (v0yn[ng] != 0):
				jfull_to_jcb[jcol] = int(jfullcol)
				jcol               = jcol + 1 

			jfullcol = jfullcol + 1

			if (widyn[ng] != 0):
				jfull_to_jcb[jcol] = int(jfullcol)
				jcol               = jcol + 1 

			jfullcol = jfullcol + 1

			if (texyn[ng] != 0):
				jfull_to_jcb[jcol] = int(jfullcol)
				jcol               = jcol + 1 

			jfullcol = jfullcol + 1

		## THEN THE WNM PARAMETERS...
		if (tcontyn != 0):
			jfull_to_jcb[jcol] = int(jfullcol)
			jcol               = jcol + 1 
		
		jfullcol = jfullcol + 1

		### BAT DAU VONG LAP ###
		redoit = 1
		while (redoit == 1):
			nloop  = nloop  + 1

			jfullcol = 0

			# Tb continuum derivative
			xdel       = np.float64(0.0000025)
			zrolvfit   = zrolvfit + xdel
			tb_totplus = self.tb_exp(xdata, \
				zrolvfit, taufit, v0fit, np.array(widfit)/dp600, texfit, \
			        tcontfit)
			        

			zrolvfit = zrolvfit - 2.*xdel
			tb_totminus = self.tb_exp(xdata, \
				zrolvfit, taufit, v0fit, np.array(widfit)/dp600, texfit, \
			        tcontfit)
			        

			zrolvfit = zrolvfit + xdel
			zrolvder = (tb_totplus - tb_totminus)/(2.*xdel)

			jacobfull[:,jfullcol] = zrolvder				#THE CONSTANT
			jfullcol              = jfullcol + 1

			# WORK THROUGH GAUSSIANS...
			for ng in range(ngaussians):
				# Tau derivative:
				xdel = np.float64(0.0000025)
				taufit[ng] = taufit[ng] + xdel
				tb_totplus = self.tb_exp(xdata, \
					zrolvfit, taufit, v0fit, np.array(widfit)/dp600, texfit, \
				        tcontfit)
				        
				taufit[ng] = taufit[ng] -2.* xdel
				tb_totminus = self.tb_exp(xdata, \
					zrolvfit, taufit, v0fit, np.array(widfit)/dp600, texfit, \
				        tcontfit)
				        
				taufit[ng] = taufit[ng] + xdel
				hgtder      = (tb_totplus - tb_totminus)/(2.*xdel)

				## V0 derivative:
				xdel       = 0.0000025*widfit[ng]
				v0fit[ng]  = v0fit[ng] + xdel
				tb_totplus = self.tb_exp(xdata, \
					zrolvfit, taufit, v0fit, np.array(widfit)/dp600, texfit, \
				        tcontfit)
				        
				v0fit[ng]   = v0fit[ng] - 2.*xdel
				tb_totminus = self.tb_exp(xdata, \
					zrolvfit, taufit, v0fit, np.array(widfit)/dp600, texfit, \
				        tcontfit)
				        
				v0fit[ng] = v0fit[ng] + xdel
				cender    = (tb_totplus - tb_totminus)/(2.*xdel)

				## WID derivative:
				xdel       = np.float64(0.0000025)*widfit[ng]
				widfit[ng] = widfit[ng] + xdel
				tb_totplus = self.tb_exp(xdata, \
					zrolvfit, taufit, v0fit, np.array(widfit)/dp600, texfit, \
				        tcontfit)
				        
				widfit[ng]  = widfit[ng] - 2.*xdel
				tb_totminus = self.tb_exp(xdata, \
					zrolvfit, taufit, v0fit, np.array(widfit)/dp600, texfit, \
				        tcontfit)
				        
				widfit[ng] = widfit[ng] + xdel
				widder     = (tb_totplus - tb_totminus)/(2.*xdel)

				## Tex derivative:
				xdel       = np.float64(0.0000025)*texfit[ng]
				texfit[ng] = texfit[ng] + xdel
				tb_totplus = self.tb_exp(xdata, \
					zrolvfit, taufit, v0fit, np.array(widfit)/dp600, texfit, \
				        tcontfit)
				        
				texfit[ng]  = texfit[ng] -2.*xdel
				tb_totminus = self.tb_exp(xdata, \
					zrolvfit, taufit, v0fit, np.array(widfit)/dp600, texfit, \
				        tcontfit)
				        
				texfit[ng]             = texfit[ng] + xdel
				tspinder               = (tb_totplus - tb_totminus)/(2.*xdel)
				jacobfull[:, jfullcol] = hgtder     #HGT-*/41-*/7
				jfullcol               = jfullcol + 1

				jacobfull[:, jfullcol] = cender     #CNTR
				jfullcol               = jfullcol + 1

				jacobfull[:, jfullcol] = widder     #WIDTH
				jfullcol               = jfullcol + 1

				jacobfull[:, jfullcol] = tspinder   #TSPIN
				jfullcol               = jfullcol + 1

			## Emission Offset derivative:
			xdel       = np.float64(0.0000025)
			tcontfit   = tcontfit + xdel
			tb_totplus = self.tb_exp(xdata, \
				zrolvfit, taufit, v0fit, np.array(widfit)/dp600, texfit, \
			        tcontfit)
			        
			tcontfit    = tcontfit - 2.*xdel
			tb_totminus = self.tb_exp(xdata, \
				zrolvfit, taufit, v0fit, np.array(widfit)/dp600, texfit, \
			        tcontfit)
			        
			tcontfit = tcontfit + xdel
			tcontder = (tb_totplus - tb_totminus)/(2.*xdel)

			jacobfull[:, jfullcol] = tcontder				#THE CONSTANT
			jfullcol               = jfullcol + 1

			jacob = jacobfull[:,jfull_to_jcb]

			# CALCULATE T_PREDICTED...
			t_predicted = self.tb_exp(xdata, \
				zrolvfit, taufit, v0fit, np.array(widfit)/dp600, texfit, \
			        tcontfit)
			        
			## CREATE AND SOLVE THE NORMAL EQUATION MATRICES...
			y    = tdata - t_predicted
			jtj  = np.dot(np.transpose(jacob),jacob)
			jy   = np.dot(np.transpose(jacob), np.transpose(y))
			c    = np.linalg.inv(jtj)
			parfull[jfull_to_jcb] = np.dot(c,jy)

			## CHECK THE DERIVED CNM PARAMETERS...
			## THE CNM AMPLITUDES...
			delt   = parfull[ [x+1 for x in (x*4 for x in list(range(ngaussians)))]  ]
			adelt  = [abs(x) for x in delt]
			taufit = [abs(x) for x in taufit]

			for i in range(len(adelt)):
				if(0.2*taufit[i] < adelt[i]):
					adelt[i] = 0.2*taufit[i]

			delttau = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					delttau.append(-adelt[i])
				else:
					delttau.append(adelt[i])

			## CNM CENTERS
			delt   = parfull[ [x+2 for x in (x*4 for x in list(range(ngaussians)))]  ]
			adelt  = [abs(x) for x in delt]
			widfit = [abs(x) for x in widfit]

			for i in range(len(adelt)):
				if(0.2*widfit[i] < adelt[i]):
					adelt[i] = 0.2*widfit[i]

			deltv0 = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					deltv0.append(-adelt[i])
				else:
					deltv0.append(adelt[i])

			## CNM WIDTHS
			delt   = parfull[ [x+3 for x in (x*4 for x in list(range(ngaussians)))]  ]
			adelt  = [abs(x) for x in delt]
			widfit = [abs(x) for x in widfit]

			for i in range(len(adelt)):
				if(0.2*widfit[i] < adelt[i]):
					adelt[i] = 0.2*widfit[i]

			deltwid = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					deltwid.append(-adelt[i])
				else:
					deltwid.append(adelt[i])

			## CNM Tex
			delt   = parfull[ [x+4 for x in (x*4 for x in list(range(ngaussians)))]  ]
			adelt  = [abs(x) for x in delt]
			texfit = [abs(x) for x in texfit]

			for i in range(len(adelt)):
				if(0.2*texfit[i] < adelt[i]):
					adelt[i] = 0.2*texfit[i]

			delttex = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					delttex.append(-adelt[i])
				else:
					delttex.append(adelt[i])

			## CHECK FOR CONVERGENCE AND IF CNM PARAMETERS are REASONABLE ##
			hgtf   = self.labs(self.ldiv(delttau,taufit))
			cenf   = self.labs(self.ldiv(deltv0,widfit))
			widf   = self.labs(self.ldiv(deltwid,widfit))
			tspinf = self.labs(self.ldiv(delttex,texfit))

			redoit = 0
			if (max(hgtf) > criterion):
				redoit = 1
			if (max(cenf) > criterion):
				redoit = 1
			if (max(widf) > criterion):
				redoit = 1
			if (max(tspinf) > criterion):
				redoit = 1

			# INCREASE THE PARAMETERS
			if(redoit == 0):
				xfact = 1.0

			zrolvfit = zrolvfit   + xfact * parfull[0]
			taufit   = self.laddl(taufit, self.listxnum(delttau,xfact))
			v0fit    = self.laddl(v0fit, self.listxnum(deltv0,xfact))
			widfit   = self.laddl(widfit, self.listxnum(deltwid,xfact))
			texfit   = self.laddl(texfit, self.listxnum(delttex,xfact))
			tcontfit = tcontfit + xfact * parfull[4*ngaussians + 1]

			# CHECK TO SEE IF WIDTH IS TOO BIG..but ignore if these params are fixed.
			if ((max(self.lxl(widyn,widfit)) > dfstop)  or \
				(min(self.lxl(widyn,widfit)) < 0.)):
			    break

		   	# nloop = 200
			if (nloop >= nloopmax-1):
				break
		
		## CONVERT THE 1/E WIDTHS TO HALFWIDTHS...
		widfit = np.array(widfit)/dp600

		## DERIVE THE FITTED POINTS, RESIDUALS, THE ERRORS of Params
		t_predicted = self.tb_exp(xdata, \
			zrolvfit, taufit, v0fit, widfit, texfit, \
		        tcontfit)
		        
		resid  = tdata - t_predicted
		resid2 = np.square(resid)
		sigsq  = resid2.sum()/(datasize - nparams)
		sigma  = sigsq**0.5

		ltemp  = list(range(nparams))
		ltemp  = [x*(nparams+1) for x in ltemp]
		c_temp = c.ravel()
		sigarray = sigsq*c_temp[ltemp]

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

		## TEST FOR INFINITIES
		countbad = 0
		indxbad  = []
		kk       = 0
		for x in np.nditer(parfull):
			if (np.isfinite(x) == False):
				countbad = countbad + 1
				indxbad.append(kk)

			kk = kk + 1

		if (countbad != 0):
			print 'Count Bad, Inf:', countbad

		erarray[jfull_to_jcb] = sigarray
		sigzrolvfit           = erarray[0]
		temp_list             = [x*4 for x in list(range(ngaussians))]
		sigtaufit             = erarray[ [x+1 for x in temp_list] ]
		sigv0fit              = erarray[ [x+2 for x in temp_list] ]
		sigwidfit             = np.array(erarray[ [x+3 for x in temp_list] ])/dp600
		sigtexfit             = erarray[ [x+4 for x in temp_list] ]
		sigtcontfit           = erarray[ 4*ngaussians + 1]


		## DERIVE THE NORMALIZED COVARIANCE ARRAY
		temp_list = [x*(nparams+1) for x in list(range(nparams))]
		c_temp    = c.ravel()
		doug      = c_temp[temp_list]
		doug_temp = doug[np.newaxis]
		doug_t    = doug[np.newaxis].T
		doug      = np.dot(doug_t,doug_temp)
		cov       = c/np.sqrt(doug)

		tfita = self.tb_exp(xdata, \
			zrolvfit, taufit, v0fit, widfit, texfit, \
		        tcontfit)
		        
		return tfita, sigma, \
				zrolvfit, taufit, v0fit, widfit, texfit, \
				sigzrolvfit, sigtaufit, sigv0fit, sigwidfit, sigtexfit, \
				tcontfit, \
				sigtcontfit, \
				cov, nloop, nloopmax
		