import numpy             as np    
from operator import truediv
a=[3.,6.,8.,65.,3.]
b=[3.,2.,5.,5.,5.]
#aa = map(truediv, a, b)
aa = np.array(a)/np.array(b)
bb = np.array(a)/np.array(b)
aa = np.square(aa)
bb = np.square(bb)
cc = np.sqrt(aa + bb)*np.array(aa)
print aa
print cc
print type(cc.tolist())