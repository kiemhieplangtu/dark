from scipy import signal
import numpy as np
import matplotlib.pyplot as plt

xs = np.arange(0, np.pi, 0.05)
data = np.sin(xs)
peakind = signal.find_peaks_cwt(data, np.arange(1,2))

print peakind, xs[peakind], data[peakind]
print np.arange(1,10)
plt.plot(xs,data)
plt.show()