import sys, os
import random
from pymc import stochastic, deterministic,Normal, Uniform
import numpy as np
from scipy.stats import multivariate_normal
import matplotlib.pyplot as plt

TRUE_SLOPE = 10
TRUE_INTERCEPT = 5.6
TRUE_OUTPUT_NOISE_PRECISION = 0.1
NUMBER_OF_SAMPLES = 1000

# Input and output data points
X_data = np.array( [random.random()*500 for i in range(NUMBER_OF_SAMPLES)] )

OUTPUT_NOISE = np.random.normal(0, 1./TRUE_OUTPUT_NOISE_PRECISION, NUMBER_OF_SAMPLES)
Y_data = TRUE_SLOPE*X_data + TRUE_INTERCEPT + OUTPUT_NOISE

# add scatter to points
X_data = np.random.normal(X_data, 20)
Y_data = np.random.normal(Y_data, 20)

# print OUTPUT_NOISE
# print len(OUTPUT_NOISE)

# hist, bins = np.histogram(OUTPUT_NOISE, bins=50)
# width = 0.7 * (bins[1] - bins[0])
# center = (bins[:-1] + bins[1:]) / 2
# plt.bar(center, hist, align='center', width=width)

# # plt.plot(X_data,Y_data,'r.')
# plt.show()

# sys.exit()

# Vector representing slope and intercept; prior on both is relatively flat Normal centered at zero
@stochastic
def w(value=np.array([1,1])):
    var = multivariate_normal(mean=[0,0], cov=[[10000,0],[0,10000]])
    return var.pdf(value)

# Input variable. Probability dist is not used, this is only here so X is observed
x = Normal('x', mu=0, tau=1, value=X_data, observed=True)

# Function mapping inputs to most likely outputs; y = mx + b
@deterministic(plot=False)
def f(x=x, w=w):
    m,b = w
    return m*x + b

# Output variable; normally distributed with mean f(x,w) and precision beta
ALPHA = 0.0001
beta = Uniform('beta', lower=0.00001, upper=1.0)
y = Normal('y', mu=f, tau=beta, value=Y_data, observed=True)