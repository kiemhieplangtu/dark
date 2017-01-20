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

	## Get expeted Tb ##
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
	def get_em_gfuncs0(self,x,amp,v0,wid,cont,tex):
		amp     = list(amp) # copy
		cntr    = list(v0)  # copy
		width   = list(wid) # copy
		ts      = list(tex) # copy

		# OPACITY TAU 
		lenx   = len(x) # sometimes = 2048
		ng     = len(amp)
		arrtau = np.zeros((lenx, ng),dtype=np.float64)

		for i in range(ng):
			arrtau[:, i] = self.gfuncs(x,0.,amp[i],cntr[i],width[i])

		sumtau = arrtau.sum(1) #2048

		exp_sumtau = np.exp(-sumtau)
		tb_cont    = cont*exp_sumtau

		tb_sum   = np.zeros(lenx)
		# BRIGHT TEMP OF EACH PEAK:
		for i in range(ng):
			temp        = np.reshape(arrtau[:, 0:i+1], (lenx, i+1))
			sumtau_i    = temp.sum(1)
			exp_tau_i   = np.exp(arrtau[:, i] - sumtau_i)
			tb_sum      = tb_sum + ts[i] * (1. - np.exp(-arrtau[:,i]) ) * exp_tau_i

		tb_tot = cont + tb_sum # 2048

		return tb_tot
	
	## Get expeted Tb ##
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
	def get_em_gfuncs(self,x,amp,v0,wid,baseline):
		amp     = list(amp) # copy
		cntr    = list(v0)  # copy
		width   = list(wid) # copy

		# Amplitude
		lenx   = len(x) # sometimes = 2048
		ng     = len(amp)
		arrtau = np.zeros((lenx, ng),dtype=np.float64)

		for i in range(ng):
			arrtau[:, i] = self.gfuncs(x,0.,amp[i],cntr[i],width[i])

		sumtau = arrtau.sum(1) #2048

		tb_tot = baseline + sumtau # 2048

		return tb_tot
		
	## Get expeted Tb ##
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

	## Calculate multiple (N) Gaucans + offset. ##
	 #
	 # params 
	 # params 
	 #
	 # return 
	 #
	 # version 08/2016 
	 # author Nguyen Van Hiep ##
	def gfuncs(self,x,offset,amp,x0,wid):
		beam_fact = np.float64(0.60056120439) # w(@1/e) = FWHM/2sqrt(ln2)
		fit       = 0.*x + offset
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
	def emfit(self,vdata, tdata, vrangeid, tbaseline, amp, v0, wid):

		beam_fact = np.float64(0.60056120439) # w(@1/e) = FWHM/2sqrt(ln2)
		vdata     = vdata.astype(np.float64)
		tdata     = tdata.astype(np.float64)

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
		fit_amp  = np.asarray(amp, dtype=np.float64)
		fit_v0   = np.asarray(v0, dtype=np.float64)
		fit_wid  = beam_fact*np.asarray(wid, dtype=np.float64)
		fit_base = np.float64(tbaseline) # Scalar 

		nloop  = 0
		nloopn = 0

		# Get Number OF GAUcANS TO FIT
		ngfit = len(amp)

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
				# AMP DERIVATIVE:
				xdel        = np.float64(0.0000001)
				fit_amp[ng] = fit_amp[ng] + xdel
				tb_totplus  = self.get_em_gfuncs(xdat, \
					fit_amp, fit_v0, self.listxnum(fit_wid,1./beam_fact),fit_base)
				        
				fit_amp[ng] = fit_amp[ng] - xdel
				tb_totminus = self.get_em_gfuncs(xdat, \
					fit_amp, fit_v0, self.listxnum(fit_wid,1./beam_fact),fit_base)
				        
				fit_amp[ng] = fit_amp[ng] + xdel  # back to fit_amp[ng]
				ampder      = (tb_totplus - tb_totminus)/(2.*xdel)

				## CEN DERIVATIVE:
				xdel       = np.float64(0.0000001)*fit_wid[ng]
				fit_v0[ng] = fit_v0[ng] + xdel
				tb_totplus = self.get_em_gfuncs(xdat, \
					fit_amp, fit_v0, self.listxnum(fit_wid,1./beam_fact),fit_base)
				        
				fit_v0[ng]  = fit_v0[ng] - xdel
				tb_totminus = self.get_em_gfuncs(xdat, \
					fit_amp, fit_v0, self.listxnum(fit_wid,1./beam_fact),fit_base)
				        
				fit_v0[ng] = fit_v0[ng] + xdel
				v0der      = (tb_totplus - tb_totminus)/(2.*xdel)

				## WID DERIVATIVE:
				xdel        = np.float64(0.0000001)*fit_wid[ng]
				fit_wid[ng] = fit_wid[ng] + xdel
				tb_totplus  = self.get_em_gfuncs(xdat, \
					fit_amp, fit_v0, self.listxnum(fit_wid,1./beam_fact),fit_base)
				        
				fit_wid[ng] = fit_wid[ng] - xdel
				tb_totminus = self.get_em_gfuncs(xdat, \
					fit_amp, fit_v0, self.listxnum(fit_wid,1./beam_fact),fit_base)
				        
				fit_wid[ng] = fit_wid[ng] + xdel
				widder      = (tb_totplus - tb_totminus)/(2.*xdel)

				jacob[:, jacobcol] = ampder     # amp AMP
				jacobcol           = jacobcol + 1

				jacob[:, jacobcol] = v0der      # CNTR
				jacobcol           = jacobcol + 1

				jacob[:, jacobcol] = widder     # WIDTH
				jacobcol           = jacobcol + 1

			## Baseline DERIVATIVE:
			xdel       = np.float64(0.0000001)
			fit_base   = fit_base + xdel
			tb_totplus = self.get_em_gfuncs(xdat,fit_amp, fit_v0,self.listxnum(fit_wid,1./beam_fact),fit_base)
			        
			fit_base    = fit_base - xdel
			tb_totminus = self.get_em_gfuncs(xdat,fit_amp,fit_v0,self.listxnum(fit_wid,1./beam_fact),fit_base)
			        
			fit_base = fit_base + xdel
			tbaseder = (tb_totplus - tb_totminus)/(2.*xdel)

			jacob[:, jacobcol] = tbaseder				# Scalar
			jacobcol           = jacobcol + 1

			# COMPUTE T_EXPECTED
			t_xpected = self.get_em_gfuncs(xdat,fit_amp, fit_v0, self.listxnum(fit_wid,1./beam_fact),fit_base)
			        
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
			fit_amp = [abs(x) for x in fit_amp]

			for i in range(len(adelt)):
				if(0.1*fit_amp[i] < adelt[i]):
					adelt[i] = 0.1*fit_amp[i]

			deltamp = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					deltamp.append(-adelt[i])
				else:
					deltamp.append(adelt[i])

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
			ampf = self.labs(self.ldiv(deltamp,fit_amp))
			cenf = self.labs(self.ldiv(deltv0,fit_wid))
			widf = self.labs(self.ldiv(deltwid,fit_wid))

			repeat = 0
			if (max(ampf) > criterion): # criterion = 0.001
				repeat = 1
			if (max(cenf) > criterion):
				repeat = 1
			if (max(widf) > criterion):
				repeat = 1

			# INCREMENT THE PARAMETERS
			if(repeat == 0):
				xfact = 1.0

			fit_amp   = self.laddl(fit_amp, self.listxnum(deltamp,xfact)) # amp(k+1) = amp(k) + dt_amp*xfact
			fit_v0    = self.laddl(fit_v0, self.listxnum(deltv0,xfact))
			fit_wid   = self.laddl(fit_wid, self.listxnum(deltwid,xfact))
			fit_base  = fit_base + xfact * dt_params[3*ngfit]

			# CHECK WIDTH
			if (min(fit_wid) < 0.):
			    break

		   	# nloopmax = 100
			if (nloop >= nloopmax-1):
				break
		
		## OK OK OK OK
		## CONVERT THE 1/E WIDTHS TO FWHM
		fit_wid = self.listxnum(fit_wid,1./beam_fact)

		## DERIVE THE FITTED POINTS, RESIDUALS, THE ERRORS 
		## WIDTH = FWHM , 0.6 FACTORS NOT REQUIRED
		t_fit = self.get_em_gfuncs(xdat,fit_amp, fit_v0, fit_wid,fit_base)
		        
		resid  = tdat - t_fit
		resid2 = np.square(resid)
		sig2   = resid2.sum()/(datasize - nparams) # reduced Sigma2
		error  = sig2**0.5

		ltemp   = list(range(nparams))
		ltemp   = [x*(nparams+1) for x in ltemp]
		c_temp  = c.ravel()
		err_arr = sig2*c_temp[ltemp]

		temp_list  = [x*3 for x in list(range(ngfit))] # ngfit=2, [0,3]
		amp_er     = err_arr[ [x for x in temp_list] ] # amp = col0, col3
		v0_er      = err_arr[ [x+1 for x in temp_list] ] #v0 = col1, col4
		wid_er     = self.listxnum(err_arr[ [x+2 for x in temp_list] ], 1./beam_fact) # wid = col2,col5
		base_er    = err_arr[3*ngfit] # col6

		## DERIVE THE NORMALIZED COVARIANCE ARRAY
		temp_list = [x*(nparams+1) for x in list(range(nparams))]
		c_temp    = c.ravel()
		doug      = c_temp[temp_list]
		doug_temp = doug[np.newaxis]
		doug_t    = doug[np.newaxis].T
		doug      = np.dot(doug_t,doug_temp)
		cov       = c/np.sqrt(doug)

		tfit = self.get_em_gfuncs(xdat,fit_amp,fit_v0,fit_wid,fit_base)
		        
		return tfit, error, \
				xdat, tdat,\
				fit_base, fit_amp, fit_v0, fit_wid, \
				base_er, amp_er, v0_er, wid_er,\
				cov, nloop, nloopmax

	##  ##
	 #
	 # params 
	 # return 	
	 #
	 # version 08/2016 
	 # author Nguyen Van Hiep ## xd,tdat,vrangid,tbaseline,amp,v0,wid,ts
	def emfit_carl(self,vdata, tdata, vrangeid, tbaseline, amp, v0, wid, ts):

		beam_fact = np.float64(0.60056120439) # w(@1/e) = FWHM/2sqrt(ln2)
		vdata     = vdata.astype(np.float64)
		tdata     = tdata.astype(np.float64)

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
		fit_amp  = np.asarray(amp, dtype=np.float64)
		fit_v0   = np.asarray(v0, dtype=np.float64)
		fit_wid  = beam_fact*np.asarray(wid, dtype=np.float64)
		fit_ts   = np.asarray(ts, dtype=np.float64)
		fit_base = np.float64(tbaseline) # Scalar 

		nloop  = 0
		nloopn = 0

		# Get Number OF GAUcANS TO FIT
		ngfit = len(amp)

		# TOTAL NR OF PARAMs
		nparams = 4*(ngfit) + 1

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
				# AMP DERIVATIVE:
				xdel        = np.float64(0.0000001)
				fit_amp[ng] = fit_amp[ng] + xdel
				tb_totplus  = self.get_em_gfuncs0(xdat, \
					fit_amp, fit_v0, self.listxnum(fit_wid,1./beam_fact),fit_base, fit_ts)
				        
				fit_amp[ng] = fit_amp[ng] - xdel
				tb_totminus = self.get_em_gfuncs0(xdat, \
					fit_amp, fit_v0, self.listxnum(fit_wid,1./beam_fact),fit_base, fit_ts)
				        
				fit_amp[ng] = fit_amp[ng] + xdel  # back to fit_amp[ng]
				ampder      = (tb_totplus - tb_totminus)/(2.*xdel)

				## CEN DERIVATIVE:
				xdel       = np.float64(0.0000001)*fit_wid[ng]
				fit_v0[ng] = fit_v0[ng] + xdel
				tb_totplus = self.get_em_gfuncs0(xdat, \
					fit_amp, fit_v0, self.listxnum(fit_wid,1./beam_fact),fit_base, fit_ts)
				        
				fit_v0[ng]  = fit_v0[ng] - xdel
				tb_totminus = self.get_em_gfuncs0(xdat, \
					fit_amp, fit_v0, self.listxnum(fit_wid,1./beam_fact),fit_base, fit_ts)
				        
				fit_v0[ng] = fit_v0[ng] + xdel
				v0der      = (tb_totplus - tb_totminus)/(2.*xdel)

				## WID DERIVATIVE:
				xdel        = np.float64(0.0000001)*fit_wid[ng]
				fit_wid[ng] = fit_wid[ng] + xdel
				tb_totplus  = self.get_em_gfuncs0(xdat, \
					fit_amp, fit_v0, self.listxnum(fit_wid,1./beam_fact),fit_base, fit_ts)
				        
				fit_wid[ng] = fit_wid[ng] - xdel
				tb_totminus = self.get_em_gfuncs0(xdat, \
					fit_amp, fit_v0, self.listxnum(fit_wid,1./beam_fact),fit_base, fit_ts)
				        
				fit_wid[ng] = fit_wid[ng] + xdel
				widder      = (tb_totplus - tb_totminus)/(2.*xdel)

				# Ts DERIVATIVE:
				xdel        = np.float64(0.0000001)
				fit_ts[ng]  = fit_ts[ng] + xdel
				tb_totplus  = self.get_em_gfuncs0(xdat, \
					fit_amp, fit_v0, self.listxnum(fit_wid,1./beam_fact),fit_base, fit_ts)
				        
				fit_ts[ng]  = fit_ts[ng] - xdel
				tb_totminus = self.get_em_gfuncs0(xdat, \
					fit_amp, fit_v0, self.listxnum(fit_wid,1./beam_fact),fit_base, fit_ts)
				        
				fit_ts[ng]  = fit_ts[ng] + xdel  # back to fit_ts[ng]
				tsder       = (tb_totplus - tb_totminus)/(2.*xdel)

				jacob[:, jacobcol] = ampder     # AMP
				jacobcol           = jacobcol + 1

				jacob[:, jacobcol] = v0der      # CNTR
				jacobcol           = jacobcol + 1

				jacob[:, jacobcol] = widder     # WIDTH
				jacobcol           = jacobcol + 1

				jacob[:, jacobcol] = tsder     # Ts
				jacobcol           = jacobcol + 1

			## Baseline DERIVATIVE:
			xdel       = np.float64(0.0000001)
			fit_base   = fit_base + xdel
			tb_totplus = self.get_em_gfuncs0(xdat,fit_amp, fit_v0,self.listxnum(fit_wid,1./beam_fact),fit_base, fit_ts)
			        
			fit_base    = fit_base - xdel
			tb_totminus = self.get_em_gfuncs0(xdat,fit_amp,fit_v0,self.listxnum(fit_wid,1./beam_fact),fit_base, fit_ts)
			        
			fit_base = fit_base + xdel
			tbaseder = (tb_totplus - tb_totminus)/(2.*xdel)

			jacob[:, jacobcol] = tbaseder				# Scalar
			jacobcol           = jacobcol + 1

			# COMPUTE T_EXPECTED
			t_xpected = self.get_em_gfuncs0(xdat,fit_amp, fit_v0, self.listxnum(fit_wid,1./beam_fact),fit_base, fit_ts)
			        
			## CREATE AND SOLVE THE NORMAL EQUATION MATRICES, [(J^T)*J]*deltaParams = (J^T)*DeltaY
			y         = tdat - t_xpected                              # deltaY
			h         = np.dot(np.transpose(jacob),jacob)             # h = (J^T)*J
			ht        = np.dot(np.transpose(jacob), np.transpose(y))  # ht = (J^T)*deltaY
			c         = np.linalg.inv(h)                              # c = h^-1
			dt_params = np.dot(c,ht)                                  # delta_params, it's a vector

			## CHECK THE OBTAINED PARAMETERS
			## AMP
			delt    = dt_params[ [x for x in (x*4 for x in list(range(ngfit)))]  ]
			adelt   = [abs(x) for x in delt]
			fit_amp = [abs(x) for x in fit_amp]

			for i in range(len(adelt)):
				if(0.1*fit_amp[i] < adelt[i]):
					adelt[i] = 0.1*fit_amp[i]

			deltamp = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					deltamp.append(-adelt[i])
				else:
					deltamp.append(adelt[i])

			## V0
			delt    = dt_params[ [x+1 for x in (x*4 for x in list(range(ngfit)))]  ]
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
			delt    = dt_params[ [x+2 for x in (x*4 for x in list(range(ngfit)))]  ]
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

			## Ts
			delt    = dt_params[ [x+3 for x in (x*4 for x in list(range(ngfit)))]  ]
			adelt   = [abs(x) for x in delt]
			fit_ts  = [abs(x) for x in fit_ts]

			for i in range(len(adelt)):
				if(0.1*fit_ts[i] < adelt[i]):
					adelt[i] = 0.1*fit_ts[i]

			deltts = []
			for i in range(len(delt)):
				if(delt[i] < 0.):
					deltts.append(-adelt[i])
				else:
					deltts.append(adelt[i])

			## CHECK FOR CONVERGENCE AND IF PARAMETERS are REASONABLE ##
			# CONVERGENCE CRITERION: ABS(Delta(Params_i)/Params_i < 0.001)
			ampf = self.labs(self.ldiv(deltamp,fit_amp))
			cenf = self.labs(self.ldiv(deltv0,fit_wid))
			widf = self.labs(self.ldiv(deltwid,fit_wid))
			tsf  = self.labs(self.ldiv(deltts,fit_ts))

			repeat = 0
			if (max(ampf) > criterion): # criterion = 0.001
				repeat = 1
			if (max(cenf) > criterion):
				repeat = 1
			if (max(widf) > criterion):
				repeat = 1
			if (max(tsf) > criterion):
				repeat = 1

			# INCREMENT THE PARAMETERS
			if(repeat == 0):
				xfact = 1.0

			fit_amp   = self.laddl(fit_amp, self.listxnum(deltamp,xfact)) # amp(k+1) = amp(k) + dt_amp*xfact
			fit_v0    = self.laddl(fit_v0, self.listxnum(deltv0,xfact))
			fit_wid   = self.laddl(fit_wid, self.listxnum(deltwid,xfact))
			fit_ts    = self.laddl(fit_ts, self.listxnum(deltts,xfact))
			fit_base  = fit_base + xfact * dt_params[4*ngfit]

			# CHECK WIDTH
			if (min(fit_wid) < 0.):
			    break

		   	# nloopmax = 100
			if (nloop >= nloopmax-1):
				break
		
		## OK OK OK OK
		## CONVERT THE 1/E WIDTHS TO FWHM
		fit_wid = self.listxnum(fit_wid,1./beam_fact)

		## DERIVE THE FITTED POINTS, RESIDUALS, THE ERRORS 
		## WIDTH = FWHM , 0.6 FACTORS NOT REQUIRED
		t_fit = self.get_em_gfuncs0(xdat,fit_amp, fit_v0, fit_wid,fit_base, fit_ts)
		        
		resid  = tdat - t_fit
		resid2 = np.square(resid)
		sig2   = resid2.sum()/(datasize - nparams) # reduced Sigma2
		error  = sig2**0.5

		ltemp   = list(range(nparams))
		ltemp   = [x*(nparams+1) for x in ltemp]
		c_temp  = c.ravel()
		err_arr = sig2*c_temp[ltemp]

		temp_list  = [x*4 for x in list(range(ngfit))] # ngfit=2, [0,3]
		amp_er     = err_arr[ [x for x in temp_list] ] # amp = col0, col3
		v0_er      = err_arr[ [x+1 for x in temp_list] ] #v0 = col1, col4
		wid_er     = self.listxnum(err_arr[ [x+2 for x in temp_list] ], 1./beam_fact) # wid = col2,col5
		ts_er      = err_arr[ [x+3 for x in temp_list] ] #v0 = col1, col4
		base_er    = err_arr[4*ngfit] # col6

		## DERIVE THE NORMALIZED COVARIANCE ARRAY
		temp_list = [x*(nparams+1) for x in list(range(nparams))]
		c_temp    = c.ravel()
		doug      = c_temp[temp_list]
		doug_temp = doug[np.newaxis]
		doug_t    = doug[np.newaxis].T
		doug      = np.dot(doug_t,doug_temp)
		cov       = c/np.sqrt(doug)

		tfit = self.get_em_gfuncs0(xdat,fit_amp,fit_v0,fit_wid,fit_base, fit_ts)

		plt.plot(xdat, tfit)
		plt.plot(xdat, tdat)
		# plt.xlim(5.,12.)
		# plt.ylim(-1.,30.)
		plt.show()
		        
		return tfit, error, \
				xdat, tdat,\
				fit_base, fit_amp, fit_v0, fit_wid, fit_ts,\
				base_er, amp_er, v0_er, wid_er,ts_er,\
				cov, nloop, nloopmax

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

        ## Width here is FWHM
		return tfit, error, \
				xdat, tdat,\
				toffset, fit_tau, fit_v0, fit_wid, \
				tau_er, v0_er, wid_er, \
				fit_cont, cont_er, \
				cov, nloop, \
		        nloopmax