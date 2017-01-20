## Import libs ##
import sys, os
import pymc
import numpy             as     np
import matplotlib.pyplot as     plt
from   pprint            import pprint

## Class - linear fit using MCMC ##
 # Fit: y = alpha * x + beta
 #          beta is uniformly distributed
 #          alpha distribution: ~ (1+beta^2)^(-3/2), see the references below
 # 
 # 
 # 
 # Using, eg:
 # 		lfit        = mclfit()
 # 		pymc_trace  = lfit.fit(xdata,ydata,yerr,brange=[-100.,100.],print_stats = False )
 # 		alpha_stats = pymc_trace['alpha_stats']
 # 		beta_stats  = pymc_trace['beta_stats']
 # 		alpha       = pymc_trace['alpha']
 # 		beta        = pymc_trace['beta']
 # 
 # Note: 'brange' is the range of Offset, depending to your Slope's value, you should change 'arange' appropriately.
 #
 # version 11/2016
 # author Nguyen Van Hiep
 # Ref: http://jakevdp.github.io/blog/2014/06/14/frequentism-and-bayesianism-4-bayesian-in-python/
 #      http://exordio.qfb.umich.mx/archivos%20pdf%20de%20trabajo%20umsnh/aphilosofia/bayesian%20importantes/leapz.pdf
 #      http://stackoverflow.com/questions/24804298/fit-a-non-linear-function-to-data-observations-with-pymcmc-pymc
 #      https://people.duke.edu/~ccc14/sta-663/PyMC2.html
 ##
class mclfit:

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
	 # params list arange Range of Offset default = [-100.,100]
	 # params Boolean print_stats Print params-Stats or not
	 #
	 # return dict {
	 #	        'alpha_stats': S.stats()['alpha'],
	 #	        'beta_stats' : S.stats()['beta'],
     #	        'alpha'      : S.trace('alpha')[:],
	 #	        'beta'       : S.trace('beta')[:]
	 #	       }
	 #
	 # version 11/2016 
	 # author Nguyen Van Hiep
	 ##
	def fit(self,x,y,yerr,brange=[-100.,100],print_stats = False):
		S = pymc.MCMC(self.model(x,y,yerr,brange))
		S.sample(iter=1e5, burn=5e4, thin=20)

		if(print_stats):
			print('alpha.stats():')
			pprint(S.stats()['alpha'])
			print('beta.stats():')
			pprint(S.stats()['beta'])

		# print S.stats()
		# pymc.Matplot.plot(S)

		# print S.alpha.trace()

		# S.alpha.summary() ## Infor for alpha
		# S.beta.summary() ## Infor for beta

		# Return the traces
		return {
		        'alpha_stats': S.stats()['alpha'],
		        'beta_stats' : S.stats()['beta'],
		        'alpha'      : S.trace('alpha')[:],
		        'beta'       : S.trace('beta')[:]
		       }

	## Define the model/function to be fitted. ##
	 #
	 # params 1D-array x x-data
	 # params 1D-array y y-data
	 # params 1D-array yerr y-error
	 # params list brange Range of Offset default = [-100.,100]
	 #
	 # return All local-variables
	 #        eg: alpha, beta, func (yfit)
	 #
	 # version 11/2016 
	 # author Nguyen Van Hiep
	 ##
	def model(self,x,y,yerr,brange):
		##
		#  func = alpha*x + beta
		##

		## beta is uniformly distributed ##
		beta = pymc.Uniform('beta', brange[0], brange[1])

		## P(alpha,beta) ~ (1+beta^2)^(-3/2) ##
		@pymc.stochastic(observed=False)
		def alpha(value=1.0e7):
			return -1.5 * np.log(1.0 + value ** 2)

		## Linear fir function: y = alpha*x + beta ##
		@pymc.deterministic
		def func(x=x, alpha=alpha, beta=beta):
			return alpha*x + beta

		f = pymc.Normal('y', mu=func, tau=1. / yerr ** 2, observed=True, value=y)

		## package the full model in a dictionary
		model1 = dict(alpha=alpha, beta=beta, func=func, y=f)
		return model1 #locals()