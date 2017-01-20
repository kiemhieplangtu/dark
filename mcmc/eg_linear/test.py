import model
from pymc import MCMC
import pprint
import sys, os

# Run sampling for 40000 iterations, with a burn-in of 2000 iterations and thinning for every 10 iterations.
M = MCMC(model)
print M

sys.exit()
M.sample(iter=40000, burn=5000, thin=10)

# Refer to sample_output.txt for example of posterior sampling summary.
pprint.pprint(M.stats())