import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as mpl

# Let's create a function to model and create data
def func(x, m, n, a, x0, sigma):
    return m*x+n+a*np.exp(-(x-x0)**2/(2*sigma**2))

# Read data
a = np.loadtxt('millennium_specs.txt', dtype='S', skiprows=54)

name_list = list(set(a[:,0]))

print len(name_list)
print name_list
print a[:,1]
print a.shape

name = a[:,0]
x = a[:,1]
y = a[:,2]

# Plot out the current state of the data and model
fig = mpl.figure()
ax = fig.add_subplot(111)

filter = (name == '3C318')
ax.plot(x[filter], y[filter], c='k', label='Function')

#ax.plot(x, y, c='k', label='Function')
#ax.scatter(x, yn)
fig.savefig('d1.png')

# Executing curve_fit on noisy data
popt, pcov = curve_fit(func, x, y)

#popt returns the best fit values for parameters of the given model (func)
print popt

ym = func(x, popt[0], popt[1], popt[2], popt[3], popt[4])
ax.plot(x, ym, c='r', label='Best fit')
ax.legend()
fig.savefig('d2.png')
