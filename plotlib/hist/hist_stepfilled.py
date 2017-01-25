import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn-white')

data = np.random.randn(1000)
# plt.hist(data)
# plt.show()


# plt.hist(data, bins=30, normed=True, alpha=0.5,
#                 histtype='stepfilled', color='steelblue',
#                 edgecolor='none')

# plt.show()

x1 = np.random.normal(0, 0.8, 1000)
x2 = np.random.normal(-2, 1, 1000)
x3 = np.random.normal(3, 2, 1000)

kwargs = dict(histtype='stepfilled', alpha=0.3, normed=True, bins=40)

plt.hist(x1, hatch='/', **kwargs)
plt.hist(x2, **kwargs)
plt.hist(x3, **kwargs)
plt.show()