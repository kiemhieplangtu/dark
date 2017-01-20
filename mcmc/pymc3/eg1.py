import sys, os
import scipy
import pymc3 as pm
import numpy as np
import matplotlib.pyplot as plt

x = np.array(range(0,50))
y = np.random.uniform(low=0.0, high=40.0, size=200)
y = map((lambda a: a[0] + a[1]), zip(x,y))

print zip(x,y)

# print y
# print len(y)

sys.exit()

# plt.scatter(x,y)

trace = None
with pm.Model() as model:
    alpha = pm.Normal('alpha', mu=0, sd=20)
    beta = pm.Normal('beta', mu=0, sd=20)
    sigma = pm.Uniform('sigma', lower=0, upper=20)
    
    y_est = alpha + beta * x
    
    likelihood = pm.Normal('y', mu=y_est, sd=sigma, observed=y)
    
    # obtain starting values via MAP
    start = pm.find_MAP()
    step = pm.NUTS(state=start)
    trace = pm.sample(2000, step, start=start, progressbar=False)
    
    # pm.traceplot(trace)

# plt.show()
a = trace['alpha']
b = trace['beta']
xfit = np.linspace(x.min(), x.max(), 20)
yfit = b[:, None] * xfit + a[:, None]
mu   = yfit.mean(0)
sig  = 1.0*yfit.std(0)
# sys.exit()

def graph(formula, x_range, color='black', alpha=1):  
    x = np.array(x_range)  
    y = eval(formula)
    plt.plot(x, y, color=color, alpha=alpha)
    
plt.scatter(x,y)

for i in xrange(0,2000):
    point = trace.point(i)
    graph('{0} + {1}*x'.format(point['alpha'], point['beta']), range(0,50), color='black', alpha=.0098035)

# graph('20 + .95*x', range(0,50), color='red', alpha=1.0)
# graph(xfit, mu, range(0,50), color='red', alpha=1.0)
plt.plot(xfit, mu, '-b', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='MCMC linear fit')
plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)

plt.show()