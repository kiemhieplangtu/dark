import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import operator
import pymc

from numpy    import array
from restore  import restore
from plotting import cplot
from pprint   import pprint

## Read data #
 # Author Van Hiep
 ##
cols = ['idx','src','nhi','xerr','ratio','yerr','oh']
fmt  = ['i',   's',  'f',  'f',   'f',    'f',   'f']
data = restore('26src_no_co_xy_err.txt', 2, cols, fmt)
dat  = data.read()

src   = dat['src']
nhi   = dat['nhi']
xerr  = dat['xerr']
ratio = dat['ratio']
yerr  = np.asarray(dat['yerr'])
oh    = dat['oh']


xdata = nhi
ydata = ratio

# Define the variables needed for the routine, with their prior distributions
alpha = pymc.Uniform('alpha', 0.37, 0.69)

@pymc.stochastic(observed=False)
def beta(value=0):
    return -1.5 * np.log(1 + value ** 2)

@pymc.stochastic(observed=False)
def sigma(value=1):
    return -np.log(abs(value))

# Define the form of the model and likelihood
@pymc.deterministic
def y_model(x=xdata, alpha=alpha, beta=beta):
    return alpha + beta * x

y = pymc.Normal('y', mu=y_model, tau=1. / yerr ** 2, observed=True, value=ydata)

# package the full model in a dictionary
model1 = dict(alpha=alpha, beta=beta,
              y_model=y_model, y=y)

# run the basic MCMC: we'll do 100000 iterations
S = pymc.MCMC(model1)
S.sample(iter=100000, burn=50000)

print('alpha.stats():')
pprint(alpha.stats())
print('\nbeta.stats():')
pprint(beta.stats())

# extract the traces and plot the results
pymc_trace = [S.trace('alpha')[:],
              S.trace('beta')[:] ]
              # S.trace('sigma')[:]]

"""Plot the linear model and 1sigma contours"""
alpha, beta = pymc_trace[:2]
xfit = np.linspace(-20, 120, 10)
yfit = alpha[:, None] + beta[:, None] * xfit
mu = yfit.mean(0)
sig = 1. * yfit.std(0)

# plt.plot(xdata, ydata, 'ok')
# plt.plot(xfit, mu, '-b')
# plt.fill_between(xfit, mu - sig, mu + sig, color='lightgray')

# plt.xlabel('x')
# plt.ylabel('y')
# # plt.plot(nhi,ratio,'r.')
# plt.xlim(0.,1.6)
# plt.ylim(0.,3.)
# # plt.show()
# plt.show()

plt.plot(xdata, ydata, 'ok', ls='None', marker='.', lw=1, label='True')
plt.errorbar(xdata,ydata,yerr=yerr, color='r', marker='.', ls='None', label='Observed')
plt.plot(xfit, mu, '-b', marker='None', ls='-', ms=5, mew=2, label='Fit')
plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)
plt.xlim(0.,1.6)
plt.ylim(0.,3.)
plt.legend()
plt.show()