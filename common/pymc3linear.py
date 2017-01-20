## Import libs ##
import sys, os
import pymc
import numpy             as     np
import matplotlib.pyplot as     plt
import pymc3             as     pm

from   pprint            import pprint

## Class - linear fit using MCMC with uniform Distribution of priors (alpha and beta) ##
 # Fit: y = alpha + beta*x
 #          beta is uniformly distributed
 #          alpha is uniformly distributed
 # 
 # Using, eg:
 # 		lfit        = pymc3lfit()
 # 		alpha,beta  = lfit.fit(xdata,ydata,yerr,arange=[-100.,100.],brange=[-100.,100] )
 # 
 # Note: 'arange' is the range of slope, 'brange' is offset, 
 #        depending to your Slope's value and offset's value, you should change 'arange' and 'brange' appropriately.
 #
 # version 11/2016
 # author Nguyen Van Hiep
 # Ref: http://jakevdp.github.io/blog/2014/06/14/frequentism-and-bayesianism-4-bayesian-in-python/
 #      http://exordio.qfb.umich.mx/archivos%20pdf%20de%20trabajo%20umsnh/aphilosofia/bayesian%20importantes/leapz.pdf
 #      http://stackoverflow.com/questions/24804298/fit-a-non-linear-function-to-data-observations-with-pymcmc-pymc
 #      https://people.duke.edu/~ccc14/sta-663/PyMC2.html
 ##
class pymc3lfit:

	## Initiate function ##
	 #
	 # params 
	 # return void
	 #
	 # version 011/2016 
	 # author Nguyen Van Hiep
	 ##
	def __init__(self):
		self.ret  = ''

	## MCMC linear fit ##
	 #
	 # params 1D-array x x-data
	 # params 1D-array y y-data
	 # params 1D-array yerr y-error
	 # params list arange Range of Slope default = [-100.,100]
	 # params list brange Range of Offset default = [-100.,100]
	 #
	 # return list [alpha[:], beta[:]]
	 #
	 # version 11/2016 
	 # author Nguyen Van Hiep
	 ##
	def fit(self,xdata,ydata,yerr,arange=[-100.,100],brange=[-100.,100]):
		trace = None
		with pm.Model() as model:
		    # alpha = pm.Normal('alpha', mu=1.0e7, sd=1.0e6)
		    # beta  = pm.Normal('beta', mu=1.0e7, sd=1.0e6)
		    # sigma = pm.Uniform('sigma', lower=0, upper=20)
		    alpha = pm.Uniform('alpha', lower=arange[0], upper=arange[1])
		    beta  = pm.Uniform('beta',  lower=brange[0], upper=brange[1])
		    sigma = yerr
		    
		    y_est = alpha + beta * xdata
		    
		    likelihood = pm.Normal('y', mu=y_est, sd=sigma, observed=ydata)
		    
		    # obtain starting values via MAP
		    start = pm.find_MAP()
		    step  = pm.NUTS(state=start)
		    trace = pm.sample(2000, step, start=start, progressbar=False)
		    
		    # pm.traceplot(trace)

		# plt.show()
		# pprint(trace['alpha'].mean())
		# pprint(trace['alpha'].std())
		# print pm.summary(trace)
		# print pm.summary(trace, ['alpha'])
		# print pm.stats()
		# print(trace.__dict__)

		# Return the traces
		return [trace['alpha'], trace['beta']]