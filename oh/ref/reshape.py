import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy.io.idl   import readsav
from numpy          import array
from restore        import restore
from plotting       import cplot


a = np.zeros((2048,3))
# b = a[:, 0:1]
# # c = np.reshape(b, (2048, 2))
# print a.shape
# print b.shape
# print b
for nrc in range(3):
	b  = a[:, 0:(nrc+1)]
	print b
	print b.shape
	temp        = np.reshape(b, (2048, nrc+1))
	tausum_nrc  = temp.sum(1)

	print nrc, temp.shape