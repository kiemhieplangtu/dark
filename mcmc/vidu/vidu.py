import numpy as np
import pymc as pm
import matplotlib.pyplot as plt

from pymc.Matplot import plot as mcplot
from pprint import pprint


# Generate noise-free data
nSamples = 101
a_true = 1 # Inclination
b_true = 0 # Offset
xs = np.linspace(0, 10, nSamples)
ys = a_true * xs + b_true 

# print xs

# plt.plot(xs,ys)
# plt.show()

# Uninformative priors
a = pm.Uniform('a', -10, 10)
b = pm.Uniform('b', -10, 10)
tau = pm.Uniform('tau', 0, 1./10**2)

print tau

@pm.deterministic
def mu(a=a, b=b):
    return a * xs + b

ys_obs = pm.Normal('ys_obs', mu, tau, value=ys, observed=True)

model = pm.Model([a, b, tau, mu, ys_obs])
mcmc = pm.MCMC(model)
mcmc.sample(iter=5e4, burn=1e4, thin=20)

print('a.stats():')
pprint(a.stats())
print('\nb.stats():')
pprint(b.stats())
