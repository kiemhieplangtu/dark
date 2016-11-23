import numpy             as np
import matplotlib.pyplot as plt
import matplotlib.cm     as cm
import matplotlib.image  as mpimg

from scipy.fftpack import ifftn

N = 100
f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')

img = mpimg.imread('a1.png')

xf = np.zeros((N,N))
xf[51, 41] = 1
xf[61, 61] = 1
Z = np.fft.fft2(xf)
ax1.imshow(xf, cmap=None)
ax4.imshow(np.real(Z), cmap=None)

# xf = np.zeros((N, N))
# xf[5, 0] = 1
# xf[N-5, 0] = 1
# Z = ifftn(xf)
# ax2.imshow(xf, cmap=cm.Reds)
# ax5.imshow(np.real(Z), cmap=cm.gray)

# xf = np.zeros((N, N))
# xf[5, 10] = 1
# xf[N-5, N-10] = 1
# Z = ifftn(xf)
# ax3.imshow(xf, cmap=cm.Reds)
# ax6.imshow(np.real(Z), cmap=cm.gray)

plt.show()

