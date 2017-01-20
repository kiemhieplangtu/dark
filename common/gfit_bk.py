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

	## Get expeted Tb ##
	 #
	 # params 1D array  x     v-data
	 # params float     bg    Backgrnd Baseline
	 # params 1D array  tau   amplitude of Tau   
	 # params 1D array  v0
	 # params 1D array  wid   width
	 # params 1D array  tex   Ex-temperature
	 # params 1D array  order Order of peaks
	 # params float     cont  Tc - Continuum
	 #
	 # return 
	 #
	 # version 08/2016 
	 # author Nguyen Van Hiep ##
	def get_tb(self,x,bg,tau,v0,wid,tex,order,cont):
		# REARRANGE CLOUDS IN ORDER OF 'ORDER'.
		t_bg    = bg
		amp     = list(tau) # copy
		cntr    = list(v0)  # copy
		width   = list(wid) # copy
		ts      = list(tex) # copy

		# CALCULATE THE OPACITY TAU OF EACH COLD CLOUD
		lenx   = len(x) # sometimes = 2048
		ng     = len(amp)
		arrtau = np.zeros((lenx, ng),dtype=np.float64)

		for i in range(ng):
			arrtau[:, i] = self.gcurv(x,t_bg,amp[i],cntr[i],width[i])

		sumtau = arrtau.sum(1) #2048

		exp_sumtau = np.exp(-sumtau)
		tb_cont    = cont*exp_sumtau
		tb_peaks   = np.zeros(lenx)

		# plt.plot(x,exp_sumtau)
		# plt.show()
		# sys.exit()

		# BRIGHT TEMP OF EACH PEAK:
		for i in range(ng):
			temp        = np.reshape(arrtau[:, 0:i+1], (lenx, i+1))
			sumtau_i    = temp.sum(1)
			exp_tau_i   = np.exp(arrtau[:, i] - sumtau_i)
			tb_peaks    = tb_peaks + ts[i] * (1. - np.exp(-arrtau[:,i]) ) * exp_tau_i

		tb_tot = tb_cont + tb_peaks # 2048

		return tb_cont, tb_peaks, tb_tot, exp_sumtau
		
	## Calculate multiple (N) Gaussians + offset. ##
	 #
	 # params 
	 # params 
	 #
	 # return 
	 #
	 # version 08/2016 
	 # author Nguyen Van Hiep ##
	def gcurv(self, x,bg,amp,x0,wid):
		# NR OF GAUSSIANS...
		if((type(amp) == list) or (type(amp) == np.ndarray)):
			ng  = len(amp)
			arr = True
		else:
			ng = 1
			arr = False

		fit = 0.*x + bg
		for i in range(ng):
			if(not arr):
				widi = wid
				ampi = amp
				x0i  = x0
			else:
				widi = wid
				ampi = amp
				x0i  = x0

			if (widi > 0.):
				fit = fit + ampi*np.exp(-( (x-x0i)/(0.6005612*widi))**2)

		return fit


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
	def fit(self,vdata, tdata, vrangeid, cont, tbg, tau, v0, wid, ts, order,
			tbgyn,tauyn,v0yn,widyn,tsyn,contyn, halfasseduse=0.5):

		dp600   = np.float64(0.60056120)
		vdata   = vdata.astype(np.float64)
		tdata   = tdata.astype(np.float64)

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

		# AX1 IS THE PERCENTAGE OF CHANGE THAT WE ALLOW# 1% IS THE DEFAULT...
		ax1 = 0.01
		ax1 = np.float64(0.003)

		# HALFASSED IS THE MULTIPLIER FOR THE CORRECTIONS IN NONLI!=AR REGIME.
		halfassed = np.float64(0.5)
		nloopmax  = 50

		#A NONZERO PROBLEM INDICATES A PROBLEM...
		problem = 0

		#DFSTOP IS THE MAXIMUM WIDTH WE ALLOW, = 80% of the total window...
		dfstop = 0.8*abs(xdat[datasize-1]-xdat[0])

		# THE OUTPUT GAUSSIAN PARAMETERS# SCALE WID FROM FWHM TO 1/E...
		# THESE ARE THE SAME AS THE PARAMETERS THAT ARE ITERATED.
		tbg1   = np.float64(tbg) # Scalar
		tau1   = np.asarray(tau, dtype=np.float64)
		v01    = np.asarray(v0, dtype=np.float64)
		wid1   = dp600*np.asarray(wid, dtype=np.float64)
		ts1    = np.asarray(ts, dtype=np.float64)
		cont1  = np.float64(cont) # Scalar 

		nloop  = 0
		nloopn = 0

		# Get Number OF GAUSSIANS TO FIT
		ngfit = len(tau)

		## get NR OF PARAMETERS TO FIT
		nparams = self.lsum(tbgyn) + self.lsum(tauyn) + self.lsum(v0yn) + self.lsum(widyn) + self.lsum(tsyn) + self.lsum(contyn)

		# TOTAL NR OF PARAMGERS THAT WOULD BE FIT IF ALL YesNo'S = 1...
		nparams_max = 2 + 4*(ngfit)

		## EQUATION-OF-CONDITION ARRAY, S AND ITS COUNTERPART SFULL...
		s          = np.zeros((datasize,nparams), dtype=np.float64)
		sfull      = np.zeros((datasize,nparams_max), dtype=np.float64)
		afull      = np.zeros(nparams_max, dtype=np.float64)
		sfull_to_s = [0]*nparams
		sigarraya  = np.zeros(nparams_max, dtype=np.float64)

		# RELATIONSHIP BETWEEN COLS IN S AND SFULL
		scol     = 0
		sfullcol = 0

		if (tbg != 0):
			sfull_to_s[scol] = int(sfullcol)
			scol             = scol + 1 

		sfullcol = sfullcol + 1

		for ng in range(ngfit):
			if (tauyn[ng] != 0):
				sfull_to_s[scol] = int(sfullcol)
				scol             = scol + 1 

			sfullcol = sfullcol + 1

			if (v0yn[ng] != 0):
				sfull_to_s[scol] = int(sfullcol)
				scol             = scol + 1 

			sfullcol = sfullcol + 1

			if (widyn[ng] != 0):
				sfull_to_s[scol] = int(sfullcol)
				scol             = scol + 1 

			sfullcol = sfullcol + 1

			if (tsyn[ng] != 0):
				sfull_to_s[scol] = int(sfullcol)
				scol             = scol + 1 

			sfullcol = sfullcol + 1

		## THEN THE BG_OFF PARAMETER
		if (contyn != 0):
			sfull_to_s[scol] = int(sfullcol)
			scol             = scol + 1 
		
		sfullcol = sfullcol + 1

		### BAT DAU VONG LAP ###
		repeat = 1
		while (repeat == 1):
			nloop  = nloop  + 1
			nloopn = nloopn + 1

			sfullcol = 0

			# EVALUATE CONSTANT DERIVATIVE:
			xdel = np.float64(0.0000025)
			tbg1 = tbg1 + xdel
			tb_cont, tb_peaks, tb_totplus, exp_tau_sum = self.get_tb(xdat, \
				tbg1, tau1, v01, self.listxnum(wid1,1./dp600), ts1, order, cont1)

			tbg1 = tbg1 - 2.*xdel
			tb_cont, tb_peaks, tb_totminus, exp_tau_sum = self.get_tb(xdat, \
				tbg1, tau1, v01, self.listxnum(wid1,1./dp600), ts1, order, cont1)			        

			tbg1   = tbg1 + xdel  # back to init-tbg
			tbgder = (tb_totplus - tb_totminus)/(2.*xdel)

			sfull[:,sfullcol] = tbgder				#THIS IS THE CONSTANT
			sfullcol          = sfullcol + 1

			# WORK THROUGH GAUSSIANS...
			for ng in range(ngfit):
				#EVALUATE TAU AMP DERIVATIVE:
				xdel = np.float64(0.0000025)
				tau1[ng] = tau1[ng] + xdel
				tb_cont, tb_peaks, tb_totplus, exp_tau_sum = self.get_tb(xdat, \
					tbg1, tau1, v01, self.listxnum(wid1,1./dp600), ts1, order, cont1)
				        
				tau1[ng] = tau1[ng] -2.* xdel
				tb_cont, tb_peaks, tb_totminus, exp_tau_sum = self.get_tb(xdat, \
					tbg1, tau1, v01, self.listxnum(wid1,1./dp600), ts1, order, cont1)
				        
				tau1[ng] = tau1[ng] + xdel
				tauder   = (tb_totplus - tb_totminus)/(2.*xdel)

				## EVALUATE CEN DERIVATIVE:
				xdel    = np.float64(0.0000025)*wid1[ng]
				v01[ng] = v01[ng] + xdel
				tb_cont, tb_peaks, tb_totplus, exp_tau_sum = self.get_tb(xdat, \
					tbg1, tau1, v01, self.listxnum(wid1,1./dp600), ts1, order, cont1)
				        
				v01[ng] = v01[ng] - 2.*xdel
				tb_cont, tb_peaks, tb_totminus, exp_tau_sum = self.get_tb(xdat, \
					tbg1, tau1, v01, self.listxnum(wid1,1./dp600), ts1, order, cont1)
				        
				v01[ng] = v01[ng] + xdel
				v0der   = (tb_totplus - tb_totminus)/(2.*xdel)

				## EVALUATE WID DERIVATIVE:
				xdel     = np.float64(0.0000025)*wid1[ng]
				wid1[ng] = wid1[ng] + xdel
				tb_cont, tb_peaks, tb_totplus, exp_tau_sum = self.get_tb(xdat, \
					tbg1, tau1, v01, self.listxnum(wid1,1./dp600), ts1, order, cont1)
				        
				wid1[ng] = wid1[ng] - 2.*xdel
				tb_cont, tb_peaks, tb_totminus, exp_tau_sum = self.get_tb(xdat, \
					tbg1, tau1, v01, self.listxnum(wid1,1./dp600), ts1, order, cont1)
				        
				wid1[ng] = wid1[ng] + xdel
				widder   = (tb_totplus - tb_totminus)/(2.*xdel)

				## EVALUATE TSPIN DERIVATIVE:
				xdel    = np.float64(0.0000025)*ts1[ng]
				ts1[ng] = ts1[ng] + xdel
				tb_cont, tb_peaks, tb_totplus, exp_tau_sum = self.get_tb(xdat, \
					tbg1, tau1, v01, self.listxnum(wid1,1./dp600), ts1, order, cont1)
				        
				ts1[ng] = ts1[ng] -2.*xdel
				tb_cont, tb_peaks, tb_totminus, exp_tau_sum = self.get_tb(xdat, \
					tbg1, tau1, v01, self.listxnum(wid1,1./dp600), ts1, order, cont1)
				        
				ts1[ng] = ts1[ng] + xdel
				tsder   = (tb_totplus - tb_totminus)/(2.*xdel)

				sfull[:, sfullcol] = tauder     # Tau AMP
				sfullcol           = sfullcol + 1

				sfull[:, sfullcol] = v0der      # CNTR
				sfullcol           = sfullcol + 1

				sfull[:, sfullcol] = widder     # WIDTH
				sfullcol           = sfullcol + 1

				sfull[:, sfullcol] = tsder      # TSPIN
				sfullcol           = sfullcol + 1

			## EVALUATE CONSTANT DERIVATIVE:
			xdel = np.float64(0.0000025)
			cont1 = cont1 + xdel
			tb_cont, tb_peaks, tb_totplus, exp_tau_sum = self.get_tb(xdat, \
					tbg1, tau1, v01, self.listxnum(wid1,1./dp600), ts1, order, cont1)
			        
			cont1 = cont1 - 2.*xdel
			tb_cont, tb_peaks, tb_totminus, exp_tau_sum = self.get_tb(xdat, \
					tbg1, tau1, v01, self.listxnum(wid1,1./dp600), ts1, order, cont1)
			        
			cont1   = cont1 + xdel
			contder = (tb_totplus - tb_totminus)/(2.*xdel)

			sfull[:, sfullcol] = contder				# Scalar
			sfullcol           = sfullcol + 1

			s = sfull[:,sfull_to_s]

			# CALCULATE T_PREDICTED...
			tb_cont, tb_peaks, t_predicted, exp_tau_sum = self.get_tb(xdat, \
					tbg1, tau1, v01, self.listxnum(wid1,1./dp600), ts1, order, cont1)
			        
			## CREATE AND SOLVE THE NORMAL EQUATION MATRICES...
			t   = tdat - t_predicted
			ss  = np.dot(np.transpose(s),s)
			st  = np.dot(np.transpose(s), np.transpose(t))
			ssi = np.linalg.inv(ss)
			a   = np.dot(ssi,st)
			afull[sfull_to_s] = a

			## CHECK THE DERIVED PARAMETERS...
			## AMP
			delt    = afull[ [x+1 for x in (x*4 for x in list(range(ngfit)))]  ]
			adelt   = [abs(x) for x in delt]
			tau1    = [abs(x) for x in tau1]

			for i in range(len(adelt)):
				if(0.2*tau1[i] < adelt[i]):
					adelt[i] = 0.2*tau1[i]

			delttau = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					delttau.append(-adelt[i])
				else:
					delttau.append(adelt[i])

			## CENTERS
			delt    = afull[ [x+2 for x in (x*4 for x in list(range(ngfit)))]  ]
			adelt   = [abs(x) for x in delt]
			wid1    = [abs(x) for x in wid1]

			for i in range(len(adelt)):
				if(0.2*wid1[i] < adelt[i]):
					adelt[i] = 0.2*wid1[i]

			deltv0 = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					deltv0.append(-adelt[i])
				else:
					deltv0.append(adelt[i])

			## WIDTHS
			delt    = afull[ [x+3 for x in (x*4 for x in list(range(ngfit)))]  ]
			adelt   = [abs(x) for x in delt]
			wid1    = [abs(x) for x in wid1]

			for i in range(len(adelt)):
				if(0.2*wid1[i] < adelt[i]):
					adelt[i] = 0.2*wid1[i]

			deltwid = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					deltwid.append(-adelt[i])
				else:
					deltwid.append(adelt[i])

			## Tex
			delt      = afull[ [x+4 for x in (x*4 for x in list(range(ngfit)))]  ]
			adelt     = [abs(x) for x in delt]
			ts1       = [abs(x) for x in ts1]

			for i in range(len(adelt)):
				if(0.2*ts1[i] < adelt[i]):
					adelt[i] = 0.2*ts1[i]

			deltts = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					deltts.append(-adelt[i])
				else:
					deltts.append(adelt[i])

			## CHECK FOR CONVERGENCE AND IF PARAMETERS are REASONABLE ##
			hgtf   = self.labs(self.ldiv(delttau,tau1))
			cenf   = self.labs(self.ldiv(deltv0,wid1))
			widf   = self.labs(self.ldiv(deltwid,wid1))
			tspinf = self.labs(self.ldiv(deltts,ts1))

			repeat = 0
			if (max(hgtf) > ax1):
				repeat = 1
			if (max(cenf) > ax1):
				repeat = 1
			if (max(widf) > ax1):
				repeat = 1
			if (max(tspinf) > ax1):
				repeat = 1

			# INCREMENT THE PARAMETERS...
			if(repeat == 0):
				halfassed = 1.0

			tbg1   = tbg1   + halfassed * afull[0]
			tau1   = self.laddl(tau1, self.listxnum(delttau,halfassed))
			v01    = self.laddl(v01, self.listxnum(deltv0,halfassed))
			wid1   = self.laddl(wid1, self.listxnum(deltwid,halfassed))
			ts1    = self.laddl(ts1, self.listxnum(deltts,halfassed))
			cont1  = cont1 + halfassed * afull[4*ngfit + 1]

			# CHECK TO SEE IF WIDTH IS TOO BIG..but ignore if these params are fixed.
			if ((max(self.lxl(widyn,wid1)) > dfstop)  or \
				(min(self.lxl(widyn,wid1)) < 0.)):
			    problem = -1
			    break

		   	# nloop = 200
			if (nloop >= nloopmax-1):
				problem = -2
				break
		
		## CONVERT THE 1/E WIDTHS TO HALFWIDTHS...
		wid1 = self.listxnum(wid1,1./dp600)

		## DERIVE THE FITTED POINTS, RESIDUALS, THE ERRORS IN DERIVED COEFFICIENTS...
		## NOTE THAT THE WIDTHS HAVE BEEN CONVERTED TO HALFWIDTHS HERE, SO THE
		## 0.6 FACTORS ARE NOT REQUIRED...
		tb_cont, tb_peaks, t_predicted, exp_tau_sum = self.get_tb(xdat, \
			tbg1, tau1, v01, wid1, ts1, order, cont1)
		        
		resid  = tdat - t_predicted
		resid2 = np.square(resid)
		sigsq  = resid2.sum()/(datasize - nparams)
		error  = sigsq**0.5

		print 'sig2',error

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
		tbg_er                = sigarraya[0]
		temp_list             = [x*4 for x in list(range(ngfit))]
		tau_er                = sigarraya[ [x+1 for x in temp_list] ]
		v0_er                 = sigarraya[ [x+2 for x in temp_list] ]
		wid_er                = self.listxnum(sigarraya[ [x+3 for x in temp_list] ], 1./dp600)
		ts1_er                = sigarraya[ [x+4 for x in temp_list] ]
		cont_er               = sigarraya[ 4*ngfit + 1]

		## DERIVE THE NORMALIZED COVARIANCE ARRAY
		temp_list = [x*(nparams+1) for x in list(range(nparams))]
		ssi_temp  = ssi.ravel()
		doug      = ssi_temp[temp_list]
		doug_temp = doug[np.newaxis]
		doug_t    = doug[np.newaxis].T
		doug      = np.dot(doug_t,doug_temp)
		cov       = ssi/np.sqrt(doug)

		tb_cont, tb_peaks, tfit, exp_tau_sum = self.get_tb(xdat, \
			tbg1, tau1, v01, wid1, ts1, order, cont1)

		plt.plot(xdat,tdat)
		plt.plot(xdat,tfit)

		plt.grid()

		plt.show()
		        
		return tfit, error, \
				tbg1, tau1, v01, wid1, ts1, \
				tbg_er, tau_er, v0_er, wid_er, ts1_er, \
				cont1, cont_er, cov, problem, nloop, \
		        tb_cont, tb_peaks, \
		        exp_tau_sum, nloopmax, halfasseduse