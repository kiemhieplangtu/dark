import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import operator
import pymc

from numpy    import array
from restore  import restore
from mclinear import mclfit
from plotting import cplot
# from pprint   import pprint

## Read data #
 # Author Van Hiep
 ##
cols  = ['idx','src','nhi','xerr','ratio','yerr','oh']
fmt   = ['i',   's',  'f',  'f',   'f',    'f',   'f']
data  = restore('26src_no_co_xy_err.txt', 2, cols, fmt)
dat   = data.read()

src   = dat['src']
nhi   = dat['nhi']
xerr  = dat['xerr']
ratio = dat['ratio']
yerr  = np.asarray(dat['yerr'])
oh    = dat['oh']

xdata = np.asarray(nhi)
ydata = np.asarray(ratio)

lfit        = mclfit()
pymc_trace  = lfit.fit(nhi,ratio,yerr,brange=[0.,2.],print_stats = False)
alpha_stats = pymc_trace['alpha_stats']
beta_stats  = pymc_trace['beta_stats']
alpha       = pymc_trace['alpha']
beta        = pymc_trace['beta']

xfit = np.linspace(xdata.min(), xdata.max(), 20)
yfit = alpha[:, None] * xfit + beta[:, None]
mu   = yfit.mean(0)
sig  = 1.0*yfit.std(0)

a  = round(111, 2)
b  = round(222, 2)
ea = round(333, 2)
eb = round(444, 2)

a  = round(alpha_stats['mean'], 2)
b  = round(beta_stats['mean'], 2)
ea = round(alpha_stats['standard deviation'], 2)
eb = round(beta_stats['standard deviation'], 2)

plt.plot(xdata, ydata, 'ok', ls='None', marker='.', lw=1, label='Factor $f = N_{HI}$/$N^*_{HI}$')
plt.errorbar(xdata,ydata,xerr=xerr, yerr=yerr, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='Observed')
plt.plot(xfit, mu, '-b', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='MCMC linear fit')
plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)
# plt.plot(xfit, alpha*xfit[:, None] + beta, c='blue', alpha=0.01)
# plt.plot(xfit, alpha_stats['mean']*xfit + beta_stats['mean'], linewidth=2, c='red')

plt.title('Correlation between Total HI column densities $N_{HI}$ and \n HI optically thin column densities $N^*_{HI}$ \
	along 79 Millennium Survey lines-of-sight', fontsize=30)
plt.ylabel('$Factor f = N_{HI}$/$N^*_{HI}$', fontsize=35)
plt.xlabel('log$_{10}$($N^*_{HI}/10^{20}$ cm$^{-2}$)', fontsize=35)
plt.xlim(0.,1.6)
plt.ylim(-0.5,4.)
plt.grid(True)
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=15)

plt.text(0., 0., 'a = '+str(a)+'$\pm$'+str(ea) +',  b = '+str(b)+'$\pm$'+str(eb), color='blue', fontsize=17)
plt.legend(loc='upper right', fontsize=18)
# plt.savefig("test.png",bbox_inches='tight')
plt.show()






# # define the model/function to be fitted.
# def model(x,f):
# 	alpha = pymc.Uniform('alpha', 0.37, 0.69)
# 	@pymc.stochastic(observed=False)
# 	def beta(value=0):
# 		return -1.5 * np.log(1 + value ** 2)
# 	@pymc.deterministic
# 	def y_model(x=xdata, alpha=alpha, beta=beta):
# 		return alpha + beta * x

# 	y = pymc.Normal('y', mu=y_model, tau=1. / yerr ** 2, observed=True, value=f)
# 	print locals()
# 	return locals()

# S = pymc.MCMC(model(xdata,ydata))
# S.sample(iter=5e4, burn=1e4, thin=20)

# print('alpha.stats():')
# pprint(S.stats()['alpha'])
# print('beta.stats():')
# pprint(S.stats()['beta'])


# # extract the traces and plot the results
# pymc_trace = [S.trace('alpha')[:],
#               S.trace('beta')[:] ]