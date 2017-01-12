import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from mpfit import mpfit

def myfunc(p, fjac=None, x=None, y=None, err=None):

	status  = 0

	# CAL. THE OPACITY TAU
	lenx   = len(x) # sometimes = 2048
	ng     = (len(p)-1)/3
	arrtau = np.zeros((lenx, ng),dtype=np.float64)

	yy = p[0]
	for i in range(1, len(p), 3):
		amp = p[i]
		ctr = p[i+1]
		wid = p[i+2]
		yy   = yy + amp * np.exp(- ( (x-ctr)/(0.6005612*wid))**2)

	return [status, (y - yy) / err]

# Generate fake data
np.random.seed(0)
x = np.linspace(-5., 5., 200)
y = 3 * np.exp(-0.5 * (x - 1.3)**2 / 0.8**2)
y += np.random.normal(0., 0.2, x.shape)

err = 0.01 * np.ones(y.shape, dtype='float64')
fa = {'x':x, 'y':y, 'err':err}
guessp = np.array([0.2,3.5,1.2,2.], dtype='float64')

# Initial model parameters
parinfo=[

      {'value': 0.2,
       'fixed': False,
       'parname': 'cont',
       'limited': [False, False]},

      {'value': 3.5,
       'fixed': False,
       'parname': 'amp1',
       'limited': [False, False]},

      {'value': 1.2,
       'fixed': False,
       'parname': 'v01',
       'limited': [False, False]},

      {'value': 2.0,
       'fixed': False,
       'parname': 'wid1',
       'limited': [False, False]}

       ]

mp = mpfit(myfunc, guessp, parinfo=parinfo, functkw=fa)
print mp
p = mp.params
yy = p[0] + p[1] * np.exp(- ( (x-p[2])/(0.6005612*p[3]))**2)
yi = guessp[0] + guessp[1] * np.exp(- ( (x-guessp[2])/(0.6005612*guessp[3]))**2)

# Plot the data with the best-fit model
plt.figure(figsize=(8,5))
plt.plot(x, y, 'ko')
plt.plot(x, yy, 'r-')
plt.plot(x, yi, 'g-')
plt.xlabel('Position')
plt.ylabel('Flux')
plt.show()