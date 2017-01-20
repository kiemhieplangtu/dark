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

	## Initiate function ##
	 #
	 # params 
	 # return void
	 #
	 # version 07/2016 
	 # author Nguyen Van Hiep ##
	def __init__(self, file=''):
		self.ret  = {}

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

	##  Sum of elements in a list ##
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

	## Get expeted Tb for Absorption line##
	 #
	 # params 1D array  x     v-data
	 # params float     bg    Backgrnd Baseline
	 # params 1D array  tau   amplitude of Tau   
	 # params 1D array  v0
	 # params 1D array  wid   width
	 # params 1D array  tex   Ex-temperature
	 # params float     cont  Tc - Continuum
	 #
	 # return 
	 #
	 # version 08/2016 
	 # author Nguyen Van Hiep ##
	def get_abs_gfuncs(self,x,offset,tau,v0,wid,cont):
		t_offset = offset
		amp      = list(tau) # copy
		cntr     = list(v0)  # copy
		width    = list(wid) # copy

		# OPACITY TAU
		lenx   = len(x) # sometimes = 2048
		ng     = len(amp)
		arrtau = np.zeros((lenx, ng),dtype=np.float64) # Initialize array for each Gaussian

		for i in range(ng):
			arrtau[:, i] = self.gfuncs(x,t_offset,amp[i],cntr[i],width[i])

		sumtau = arrtau.sum(1) #2048 sum of N gaucans

		exp_sumtau   = np.exp(-sumtau) # exp(-tau) from the [sum of N gaucans]
		spec_xpected = cont*exp_sumtau
		# TO HERE: 1st fitting, opacity spectrum e(-tau) where, tau = sum of Gaucans

		return spec_xpected

	## Get expeted Tb ##
	 #
	 # params 1D array  x     v-data
	 # params float     bg    Backgrnd Baseline
	 # params 1D array  tau   amplitude of Tau   
	 # params 1D array  v0    Center
	 # params 1D array  wid   width
	 # params 1D array  tex   Ex-temperature
	 # params float     cont  Tc - Continuum
	 #
	 # return Expected temperature
	 #
	 # version 08/2016 
	 # author Nguyen Van Hiep ##
	def get_tb(self,x,bg,tau,v0,wid,tex,cont):
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
			arrtau[:, i] = self.gfuncs(x,bg,amp[i],cntr[i],width[i])

		sumtau = arrtau.sum(1) #2048

		exp_sumtau = np.exp(-sumtau)
		tb_cont    = cont*exp_sumtau

		tb_peaks   = np.zeros(lenx)
		# BRIGHT TEMP OF EACH PEAK:
		for i in range(ng):
			temp        = np.reshape(arrtau[:, 0:i+1], (lenx, i+1))
			temp        = arrtau[:, 0:i+1]
			sumtau_i    = temp.sum(1)
			exp_tau_i   = np.exp(-sumtau_i)
			# tb_peaks    = tb_peaks + ts[i] * (1. - np.exp(-arrtau[:,i]) ) * exp_tau_i ## From Carl's paper
			tb_peaks    = tb_peaks + ts[i] * (1. - np.exp(-arrtau[:,i]) ) * exp_tau_i ## 

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
	def gfuncs(self, x, bg, amp, x0, wid):
		beam_fact = np.float64(0.60056120439) # w(@1/e) = FWHM/2sqrt(ln2)
		fit       = 0.*x + bg
		if (wid > 0.):
			fit = fit + amp*np.exp(-( (x-x0)/(beam_fact*wid))**2) # Beam = Width at 1/e

		return fit

	##  ##
	 #
	 # params 
	 # return 	
	 #
	 # version 08/2016 
	 # author Nguyen Van Hiep ##
	def fit(self,vdata, tdata, vrangeid, bg, tau, v0, wid, ts, cont,
			bgyn, ampyn,v0yn,widyn,tsyn,contyn):

		beam_fact = np.float64(0.60056120439)
		vdata     = vdata.astype(np.float64)
		tdata     = tdata.astype(np.float64)

		nrange = len(vrangeid)/2
		datasize = 0L
		for i in range(nrange):
			datasize = datasize + vrangeid[2*i+1]-vrangeid[2*i]+1L

		xdat = np.zeros(datasize, dtype=np.float64)
		tdat = np.zeros(datasize, dtype=np.float64)

		## Get each Segment of Data considered to Fit
		## bini: initial Channel, binf: final Channal for each Segment ##
		bini = 0L
		for i in range(nrange):
			binf            = bini + vrangeid[2*i+1]-vrangeid[2*i] +1L
			xdat[bini:binf] = vdata[vrangeid[2*i]:vrangeid[2*i+1]+1]
			tdat[bini:binf] = tdata[vrangeid[2*i]:vrangeid[2*i+1]+1]

		# CONVERGENCE criterion
		criterion = np.float64(0.001)

		# MULTIPLIER FOR THE CORRECTIONS
		xfact     = np.float64(0.5)
		nloopmax  = 100

		# VRange LIMIT
		dv_lim = 0.8*abs(xdat[datasize-1]-xdat[0])

		# Fit GAUSSIAN PARAMETERS
		amp_fit   = np.asarray(tau, dtype=np.float64)
		v0_fit    = np.asarray(v0, dtype=np.float64)
		wid_fit   = beam_fact*np.asarray(wid, dtype=np.float64)
		ts_fit    = np.asarray(ts, dtype=np.float64)
		cont_fit  = np.float64(cont) # Scalar 
		bg_fit    = np.float64(bg) # Scalar 

		nloop  = 0
		nloopn = 0

		# Get Number OF GAUSSIANS TO FIT
		ngfit = len(tau)

		## get NR OF PARAMETERS TO FIT
		nparams = self.lsum(bgyn) + self.lsum(ampyn) + self.lsum(v0yn) + self.lsum(widyn) + self.lsum(tsyn) + self.lsum(contyn)

		# TOTAL NR OF PARAMGERS THAT WOULD BE FIT
		nparams_max = 2 + 4*(ngfit)

		## EQUATION-OF-CONDITION ARRAY
		jacobfit      = np.zeros((datasize,nparams), dtype=np.float64)
		jacobfull     = np.zeros((datasize,nparams_max), dtype=np.float64)
		delta_params  = np.zeros(nparams_max, dtype=np.float64)
		jfull_to_jfit = [0]*nparams
		err_array     = np.zeros(nparams_max, dtype=np.float64)

		# RELATIONSHIP BETWEEN COLS IN Jacobian-Fit AND Jacobian-Full
		jcol     = 0
		jfullcol = 0

		if (bgyn != 0):
			jfull_to_jfit[jcol] = int(jfullcol)
			jcol                = jcol + 1

		jfullcol = jfullcol + 1

		for ng in range(ngfit):
			if (ampyn[ng] != 0):
				jfull_to_jfit[jcol] = int(jfullcol)
				jcol                = jcol + 1 

			jfullcol = jfullcol + 1

			if (v0yn[ng] != 0):
				jfull_to_jfit[jcol] = int(jfullcol)
				jcol             = jcol + 1 

			jfullcol = jfullcol + 1

			if (widyn[ng] != 0):
				jfull_to_jfit[jcol] = int(jfullcol)
				jcol             = jcol + 1 

			jfullcol = jfullcol + 1

			if (tsyn[ng] != 0):
				jfull_to_jfit[jcol] = int(jfullcol)
				jcol             = jcol + 1 

			jfullcol = jfullcol + 1

		## THE CONTINUUM BG_OFF PARAMETER
		if (contyn != 0):
			jfull_to_jfit[jcol] = int(jfullcol)
			jcol             = jcol + 1 
		
		jfullcol = jfullcol + 1

		### BAT DAU VONG LAP ###
		repeat = 1
		while (repeat == 1):
			nloop  = nloop  + 1
			nloopn = nloopn + 1

			jfullcol = 0

			# BG DERIVATIVE:
			xdel = np.float64(0.000001)
			bg_fit     = bg_fit + xdel
			tb_totplus = self.get_tb(xdat, \
				bg_fit,amp_fit, v0_fit, self.listxnum(wid_fit,1./beam_fact), ts_fit, cont_fit)
			        
			bg_fit      = bg_fit -2.*xdel
			tb_totminus = self.get_tb(xdat, \
				bg_fit,amp_fit, v0_fit, self.listxnum(wid_fit,1./beam_fact), ts_fit, cont_fit)
			        
			bg_fit    = bg_fit + xdel
			bgder     = (tb_totplus - tb_totminus)/(2.*xdel)

			jacobfull[:, jfullcol] = bgder     # Tau AMP
			jfullcol               = jfullcol + 1

			# WORK THROUGH GAUSSIANS...
			for ng in range(ngfit):
				# TAU AMP DERIVATIVE:
				xdel = np.float64(0.000001)
				amp_fit[ng] = amp_fit[ng] + xdel
				tb_totplus = self.get_tb(xdat, \
					bg_fit,amp_fit, v0_fit, self.listxnum(wid_fit,1./beam_fact), ts_fit, cont_fit)
				        
				amp_fit[ng] = amp_fit[ng] -2.*xdel
				tb_totminus = self.get_tb(xdat, \
					bg_fit,amp_fit, v0_fit, self.listxnum(wid_fit,1./beam_fact), ts_fit, cont_fit)
				        
				amp_fit[ng] = amp_fit[ng] + xdel
				tauder   = (tb_totplus - tb_totminus)/(2.*xdel)

				## CEN DERIVATIVE:
				xdel    = np.float64(0.000001)*wid_fit[ng]
				v0_fit[ng] = v0_fit[ng] + xdel
				tb_totplus = self.get_tb(xdat, \
					bg_fit,amp_fit, v0_fit, self.listxnum(wid_fit,1./beam_fact), ts_fit, cont_fit)
				        
				v0_fit[ng] = v0_fit[ng] - 2.*xdel
				tb_totminus = self.get_tb(xdat, \
					bg_fit,amp_fit, v0_fit, self.listxnum(wid_fit,1./beam_fact), ts_fit, cont_fit)
				        
				v0_fit[ng] = v0_fit[ng] + xdel
				v0der   = (tb_totplus - tb_totminus)/(2.*xdel)

				## WID DERIVATIVE:
				xdel     = np.float64(0.000001)*wid_fit[ng]
				wid_fit[ng] = wid_fit[ng] + xdel
				tb_totplus = self.get_tb(xdat, \
					bg_fit,amp_fit, v0_fit, self.listxnum(wid_fit,1./beam_fact), ts_fit, cont_fit)
				        
				wid_fit[ng] = wid_fit[ng] - 2.*xdel
				tb_totminus = self.get_tb(xdat, \
					bg_fit,amp_fit, v0_fit, self.listxnum(wid_fit,1./beam_fact), ts_fit, cont_fit)
				        
				wid_fit[ng] = wid_fit[ng] + xdel
				widder   = (tb_totplus - tb_totminus)/(2.*xdel)

				## TSPIN DERIVATIVE:
				xdel    = np.float64(0.000001)*ts_fit[ng]
				ts_fit[ng] = ts_fit[ng] + xdel
				tb_totplus = self.get_tb(xdat, \
					bg_fit,amp_fit, v0_fit, self.listxnum(wid_fit,1./beam_fact), ts_fit, cont_fit)
				        
				ts_fit[ng] = ts_fit[ng] - 2.*xdel
				tb_totminus = self.get_tb(xdat, \
					bg_fit,amp_fit, v0_fit, self.listxnum(wid_fit,1./beam_fact), ts_fit, cont_fit)
				        
				ts_fit[ng] = ts_fit[ng] + xdel
				tsder   = (tb_totplus - tb_totminus)/(2.*xdel)

				jacobfull[:, jfullcol] = tauder     # Tau AMP
				jfullcol               = jfullcol + 1

				jacobfull[:, jfullcol] = v0der      # CNTR
				jfullcol               = jfullcol + 1

				jacobfull[:, jfullcol] = widder     # WIDTH
				jfullcol               = jfullcol + 1

				jacobfull[:, jfullcol] = tsder      # TSPIN
				jfullcol               = jfullcol + 1

			## CONSTANT DERIVATIVE:
			xdel = np.float64(0.000001)
			cont_fit = cont_fit + xdel
			tb_totplus = self.get_tb(xdat, \
					bg_fit,amp_fit, v0_fit, self.listxnum(wid_fit,1./beam_fact), ts_fit, cont_fit)
			        
			cont_fit = cont_fit - 2.*xdel
			tb_totminus = self.get_tb(xdat, \
					bg_fit,amp_fit, v0_fit, self.listxnum(wid_fit,1./beam_fact), ts_fit, cont_fit)
			        
			cont_fit   = cont_fit + xdel
			contder = (tb_totplus - tb_totminus)/(2.*xdel)

			jacobfull[:, jfullcol] = contder				# Scalar
			jfullcol           = jfullcol + 1

			jacobfit = jacobfull[:,jfull_to_jfit]

			# CALCULATE T_PREDICTED #
			t_predicted = self.get_tb(xdat, \
					bg_fit,amp_fit, v0_fit, self.listxnum(wid_fit,1./beam_fact), ts_fit, cont_fit)
			        
			## NORMAL EQUATION MATRICES [(J^T)*J]*deltaParams = (J^T)*DeltaY
			y   = tdat - t_predicted                              # deltaY
			ss  = np.dot(np.transpose(jacobfit),jacobfit)         # ss = (J^T)*J
			st  = np.dot(np.transpose(jacobfit), np.transpose(y)) # st = (J^T)*deltaY
			ssi = np.linalg.inv(ss)                               # ssi = C = h^-1
			delta_params[jfull_to_jfit] = np.dot(ssi,st)          # delta_params, it's a vector

			## CHECK THE PARAMETERS
			## AMP
			delt    = delta_params[ [x+1 for x in (x*4 for x in list(range(ngfit)))]  ]
			adelt   = [abs(x) for x in delt]
			amp_fit = [abs(x) for x in amp_fit]

			for i in range(len(adelt)):
				if(0.2*amp_fit[i] < adelt[i]):
					adelt[i] = 0.2*amp_fit[i]

			delttau = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					delttau.append(-adelt[i])
				else:
					delttau.append(adelt[i])

			## CENTERS
			delt    = delta_params[ [x+2 for x in (x*4 for x in list(range(ngfit)))]  ]
			adelt   = [abs(x) for x in delt]
			wid_fit = [abs(x) for x in wid_fit]

			for i in range(len(adelt)):
				if(0.2*wid_fit[i] < adelt[i]):
					adelt[i] = 0.2*wid_fit[i]

			deltv0 = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					deltv0.append(-adelt[i])
				else:
					deltv0.append(adelt[i])

			## WIDTHS
			delt    = delta_params[ [x+3 for x in (x*4 for x in list(range(ngfit)))]  ]
			adelt   = [abs(x) for x in delt]
			wid_fit = [abs(x) for x in wid_fit]

			for i in range(len(adelt)):
				if(0.2*wid_fit[i] < adelt[i]):
					adelt[i] = 0.2*wid_fit[i]

			deltwid = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					deltwid.append(-adelt[i])
				else:
					deltwid.append(adelt[i])

			## Tex
			delt      = delta_params[ [x+4 for x in (x*4 for x in list(range(ngfit)))]  ]
			adelt     = [abs(x) for x in delt]
			ts_fit    = [abs(x) for x in ts_fit]

			for i in range(len(adelt)):
				if(0.2*ts_fit[i] < adelt[i]):
					adelt[i] = 0.2*ts_fit[i]

			deltts = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					deltts.append(-adelt[i])
				else:
					deltts.append(adelt[i])

			## CHECK FOR CONVERGENCE AND IF PARAMETERS are REASONABLE ##
			hgtf   = self.labs(self.ldiv(delttau,amp_fit))
			cenf   = self.labs(self.ldiv(deltv0,wid_fit))
			widf   = self.labs(self.ldiv(deltwid,wid_fit))
			tspinf = self.labs(self.ldiv(deltts,ts_fit))

			repeat = 0
			if (max(hgtf) > criterion):
				repeat = 1
			if (max(cenf) > criterion):
				repeat = 1
			if (max(widf) > criterion):
				repeat = 1
			if (max(tspinf) > criterion):
				repeat = 1

			# INCREMENT THE PARAMETERS...
			if(repeat == 0):
				xfact = 1.0

			bg_fit    = bg_fit + xfact * delta_params[0]
			amp_fit   = self.laddl(amp_fit, self.listxnum(delttau,xfact))
			v0_fit    = self.laddl(v0_fit, self.listxnum(deltv0,xfact))
			wid_fit   = self.laddl(wid_fit, self.listxnum(deltwid,xfact))
			ts_fit    = self.laddl(ts_fit, self.listxnum(deltts,xfact))
			cont_fit  = cont_fit + xfact * delta_params[4*ngfit+1]

			# CHECK TO SEE IF WIDTH IS TOO BIG..but ignore if these params are fixed.
			if (min(self.lxl(widyn,wid_fit)) < 0.):
			    break

		   	# nloop = 200
			if (nloop >= nloopmax-1):
				break
		
		## CONVERT THE 1/E WIDTHS TO HALFWIDTHS...
		wid_fit = self.listxnum(wid_fit,1./beam_fact)

		## DERIVE THE FITTED POINTS, RESIDUALS, ERRORS OF DERIVED COEFFICIENTS
		## WIDTHS = FWHM
		t_predicted = self.get_tb(xdat, bg_fit,amp_fit, v0_fit, wid_fit, ts_fit, cont_fit)
		        
		resid  = tdat - t_predicted
		resid2 = np.square(resid)
		sigsq  = resid2.sum()/(datasize - nparams)
		error  = sigsq**0.5

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

		## CHECK FOR NEG SQRTS
		if (countsqrt != 0):
			sigarray[indxsqrt] = -sigarray[indxsqrt]
			problem = -3

		## CHECK FOR INFINITIES
		countbad = 0
		indxbad  = []
		kk       = 0
		for x in np.nditer(delta_params):
			if (np.isfinite(x) == False):
				countbad = countbad + 1
				indxbad.append(kk)

			kk = kk + 1

		if (countbad != 0):
			problem = -4

		err_array[jfull_to_jfit] = sigarray
		bg_er                 = err_array[0]
		temp_list             = [x*4 for x in list(range(ngfit))]
		tau_er                = err_array[ [x+1 for x in temp_list] ]
		v0_er                 = err_array[ [x+2 for x in temp_list] ]
		wid_er                = self.listxnum(err_array[ [x+3 for x in temp_list] ], 1./beam_fact)
		ts_er                 = err_array[ [x+4 for x in temp_list] ]
		cont_er               = err_array[ 4*ngfit+1]

		## DERIVE THE NORMALIZED COVARIANCE ARRAY
		temp_list = [x*(nparams+1) for x in list(range(nparams))]
		ssi_temp  = ssi.ravel()
		doug      = ssi_temp[temp_list]
		doug_temp = doug[np.newaxis]
		doug_t    = doug[np.newaxis].T
		doug      = np.dot(doug_t,doug_temp)
		cov       = ssi/np.sqrt(doug)

		tfit = self.get_tb(xdat, bg_fit,amp_fit, v0_fit, wid_fit, ts_fit, cont_fit)
		        
		return tfit, error, \
				xdat,tdat,\
				bg_fit,amp_fit, v0_fit, wid_fit, ts_fit, \
				bg_er, tau_er, v0_er, wid_er, ts_er, \
				cont_fit, cont_er, cov, nloop, \
		        nloopmax

    ##  ##
	 #
	 # params 
	 # return 	
	 #
	 # version 08/2016 
	 # author Nguyen Van Hiep ##
	def abfit(self,vdata, tdata, vrangeid, cont, 
		specbg_offset, tau, v0, wid):

		beam_fact   = np.float64(0.60056120)
		vdata   = vdata.astype(np.float64)
		tdata   = tdata.astype(np.float64)

		nrange   = len(vrangeid)/2
		datasize = 0L
		for i in range(nrange):
			datasize = datasize + vrangeid[2*i+1]-vrangeid[2*i]+1L

		xdat = np.zeros(datasize, dtype=np.float64)
		tdat = np.zeros(datasize, dtype=np.float64)

		## Get each Segment of Data considered to Fit
		## bini: initial Channel, binf: final Channal for each Segment ##
		bini = 0L
		for i in range(nrange):
			binf            = bini + vrangeid[2*i+1]-vrangeid[2*i] +1L
			xdat[bini:binf] = vdata[vrangeid[2*i]:vrangeid[2*i+1]+1]
			tdat[bini:binf] = tdata[vrangeid[2*i]:vrangeid[2*i+1]+1]

		# CONVERGENCE CRITERION: ABS(Delta(Params_i)/Params_i < 0.001)
		criterion = np.float64(0.001)

		# xfact IS THE MULTIPLIER FOR THE CORRECTIONS IN NONLINEAR
		xfact     = np.float64(0.5)
		nloopmax  = 100

		# THE OUTPUT GAUcAN PARAMETERS# SCALE WID FROM FWHM TO 1/E...
		# THESE ARE THE SAME AS THE PARAMETERS THAT ARE ITERATED.
		toffset  = np.float64(specbg_offset) # Scalar
		fit_tau  = np.asarray(tau, dtype=np.float64)
		fit_v0   = np.asarray(v0, dtype=np.float64)
		fit_wid  = beam_fact*np.asarray(wid, dtype=np.float64)
		fit_cont = np.float64(cont) # Scalar 

		nloop  = 0
		nloopn = 0

		# Get Number OF GAUcANS TO FIT
		ngfit = len(tau)

		# TOTAL NR OF PARAMs
		nparams = 3*(ngfit) + 1

		## EQUATION-OF-CONDITION ARRAY
		jacob     = np.zeros((datasize,nparams), dtype=np.float64)
		dt_params = np.zeros(nparams, dtype=np.float64)
		err_arr   = np.zeros(nparams, dtype=np.float64)

		### BAT DAU VONG LAP ###
		repeat = 1
		while (repeat == 1):
			nloop  = nloop  + 1
			nloopn = nloopn + 1

			jacobcol = 0

			# WORK THROUGH GAUcANS...
			for ng in range(ngfit):
				#TAU AMP DERIVATIVE:
				xdel        = np.float64(0.0000001)
				fit_tau[ng] = fit_tau[ng] + xdel
				tb_totplus  = self.get_abs_gfuncs(xdat, \
					toffset, fit_tau, fit_v0, self.listxnum(fit_wid,1./beam_fact),fit_cont)
				        
				fit_tau[ng] = fit_tau[ng] - xdel
				tb_totminus = self.get_abs_gfuncs(xdat, \
					toffset, fit_tau, fit_v0, self.listxnum(fit_wid,1./beam_fact),fit_cont)
				        
				fit_tau[ng] = fit_tau[ng] + xdel  # back to fit_tau[ng]
				tauder      = (tb_totplus - tb_totminus)/(2.*xdel)

				## CEN DERIVATIVE:
				xdel       = np.float64(0.0000001)*fit_wid[ng]
				fit_v0[ng] = fit_v0[ng] + xdel
				tb_totplus = self.get_abs_gfuncs(xdat, \
					toffset, fit_tau, fit_v0, self.listxnum(fit_wid,1./beam_fact),fit_cont)
				        
				fit_v0[ng]  = fit_v0[ng] - xdel
				tb_totminus = self.get_abs_gfuncs(xdat, \
					toffset, fit_tau, fit_v0, self.listxnum(fit_wid,1./beam_fact),fit_cont)
				        
				fit_v0[ng] = fit_v0[ng] + xdel
				v0der      = (tb_totplus - tb_totminus)/(2.*xdel)

				## WID DERIVATIVE:
				xdel        = np.float64(0.0000001)*fit_wid[ng]
				fit_wid[ng] = fit_wid[ng] + xdel
				tb_totplus  = self.get_abs_gfuncs(xdat, \
					toffset, fit_tau, fit_v0, self.listxnum(fit_wid,1./beam_fact),fit_cont)
				        
				fit_wid[ng] = fit_wid[ng] - xdel
				tb_totminus = self.get_abs_gfuncs(xdat, \
					toffset, fit_tau, fit_v0, self.listxnum(fit_wid,1./beam_fact),fit_cont)
				        
				fit_wid[ng] = fit_wid[ng] + xdel
				widder      = (tb_totplus - tb_totminus)/(2.*xdel)

				jacob[:, jacobcol] = tauder     # Tau AMP
				jacobcol           = jacobcol + 1

				jacob[:, jacobcol] = v0der      # CNTR
				jacobcol           = jacobcol + 1

				jacob[:, jacobcol] = widder     # WIDTH
				jacobcol           = jacobcol + 1

			## Tc DERIVATIVE:
			xdel       = np.float64(0.0000001)
			fit_cont   = fit_cont + xdel
			tb_totplus = self.get_abs_gfuncs(xdat,toffset,fit_tau, fit_v0,self.listxnum(fit_wid,1./beam_fact),fit_cont)
			        
			fit_cont    = fit_cont - xdel
			tb_totminus = self.get_abs_gfuncs(xdat,toffset,fit_tau,fit_v0,self.listxnum(fit_wid,1./beam_fact),fit_cont)
			        
			fit_cont = fit_cont + xdel
			contder  = (tb_totplus - tb_totminus)/(2.*xdel)

			jacob[:, jacobcol] = contder				# Scalar
			jacobcol           = jacobcol + 1

			# COMPUTE T_EXPECTED
			t_xpected = self.get_abs_gfuncs(xdat,toffset, fit_tau, fit_v0, self.listxnum(fit_wid,1./beam_fact),fit_cont)
			        
			## CREATE AND SOLVE THE NORMAL EQUATION MATRICES, [(J^T)*J]*deltaParams = (J^T)*DeltaY
			y         = tdat - t_xpected                              # deltaY
			h         = np.dot(np.transpose(jacob),jacob)             # h = (J^T)*J
			ht        = np.dot(np.transpose(jacob), np.transpose(y))  # ht = (J^T)*deltaY
			c         = np.linalg.inv(h)                              # c = h^-1
			dt_params = np.dot(c,ht)                                  # delta_params, it's a vector

			## CHECK THE OBTAINED PARAMETERS
			## AMP
			delt    = dt_params[ [x for x in (x*3 for x in list(range(ngfit)))]  ]
			adelt   = [abs(x) for x in delt]
			fit_tau = [abs(x) for x in fit_tau]

			for i in range(len(adelt)):
				if(0.1*fit_tau[i] < adelt[i]):
					adelt[i] = 0.1*fit_tau[i]

			delttau = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					delttau.append(-adelt[i])
				else:
					delttau.append(adelt[i])

			## V0
			delt    = dt_params[ [x+1 for x in (x*3 for x in list(range(ngfit)))]  ]
			adelt   = [abs(x) for x in delt]
			fit_wid = [abs(x) for x in fit_wid]

			for i in range(len(adelt)):
				if(0.1*fit_wid[i] < adelt[i]):
					adelt[i] = 0.1*fit_wid[i]

			deltv0 = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					deltv0.append(-adelt[i])
				else:
					deltv0.append(adelt[i])

			## WIDTH
			delt    = dt_params[ [x+2 for x in (x*3 for x in list(range(ngfit)))]  ]
			adelt   = [abs(x) for x in delt]
			fit_wid = [abs(x) for x in fit_wid]

			for i in range(len(adelt)):
				if(0.1*fit_wid[i] < adelt[i]):
					adelt[i] = 0.1*fit_wid[i]

			deltwid = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					deltwid.append(-adelt[i])
				else:
					deltwid.append(adelt[i])

			## CHECK FOR CONVERGENCE AND IF PARAMETERS are REASONABLE ##
			# CONVERGENCE CRITERION: ABS(Delta(Params_i)/Params_i < 0.001)
			tauf = self.labs(self.ldiv(delttau,fit_tau))
			cenf = self.labs(self.ldiv(deltv0,fit_wid))
			widf = self.labs(self.ldiv(deltwid,fit_wid))

			repeat = 0
			if (max(tauf) > criterion): # criterion = 0.001
				repeat = 1
			if (max(cenf) > criterion):
				repeat = 1
			if (max(widf) > criterion):
				repeat = 1

			# INCREMENT THE PARAMETERS
			if(repeat == 0):
				xfact = 1.0

			fit_tau   = self.laddl(fit_tau, self.listxnum(delttau,xfact)) # tau(k+1) = tau(k) + dt_tau*xfact
			fit_v0    = self.laddl(fit_v0, self.listxnum(deltv0,xfact))
			fit_wid   = self.laddl(fit_wid, self.listxnum(deltwid,xfact))
			fit_cont  = fit_cont + xfact * dt_params[3*ngfit]

			# CHECK WIDTH
			if (min(fit_wid) < 0.):
			    break

		   	# nloopmax = 100
			if (nloop >= nloopmax-1):
				break
		
		## OK OK OK OK
		## CONVERT THE 1/E WIDTHS TO FWHM...
		fit_wid = self.listxnum(fit_wid,1./beam_fact)

		## DERIVE THE FITTED POINTS, RESIDUALS, THE ERRORS 
		## WIDTH = FWHM HERE, 0.6 FACTORS NOT REQUIRED
		t_fit = self.get_abs_gfuncs(xdat,toffset, fit_tau, fit_v0, fit_wid,fit_cont)
		        
		resid  = tdat - t_fit
		resid2 = np.square(resid)
		sig2   = resid2.sum()/(datasize - nparams) # reduced Sigma2
		error  = sig2**0.5

		ltemp   = list(range(nparams))
		ltemp   = [x*(nparams+1) for x in ltemp]
		c_temp  = c.ravel()
		err_arr = sig2*c_temp[ltemp]

		temp_list  = [x*3 for x in list(range(ngfit))] # ngfit=2, [0,3]
		tau_er     = err_arr[ [x for x in temp_list] ] # tau = col0, col3
		v0_er      = err_arr[ [x+1 for x in temp_list] ] #v0 = col1, col4
		wid_er     = self.listxnum(err_arr[ [x+2 for x in temp_list] ], 1./beam_fact) # wid = col2,col5
		cont_er    = err_arr[3*ngfit] # col6

		## DERIVE THE NORMALIZED COVARIANCE ARRAY
		temp_list = [x*(nparams+1) for x in list(range(nparams))]
		c_temp    = c.ravel()
		doug      = c_temp[temp_list]
		doug_temp = doug[np.newaxis]
		doug_t    = doug[np.newaxis].T
		doug      = np.dot(doug_t,doug_temp)
		cov       = c/np.sqrt(doug)

		tfit = self.get_abs_gfuncs(xdat,toffset,fit_tau,fit_v0,fit_wid,fit_cont)

		plt.plot(xdat,tdat)
		plt.plot(xdat,tfit)

		print 'tdat - average', np.average(tdat)

		plt.grid()

		plt.show()

        ## Width here is FWHM
		return tfit, error, \
				toffset, fit_tau, fit_v0, fit_wid, \
				tau_er, v0_er, wid_er, \
				fit_cont, cont_er, \
				cov, nloop, \
		        nloopmax