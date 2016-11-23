import matplotlib.pyplot as plt
import numpy             as np

from numpy  import array
import operator

# tau = 0.669
# e_tau = 0.011
# dv = 2.43
# e_dv = 0.03

# ts = 34.25
# e_ts = 7.89

tau = 0.077
e_tau = 0.019
dv = 8.62
e_dv = 0.7

ts = 267.80
e_ts = 8.28

fct = 2.0*((np.log(2))**0.5)
const = 1.8224

dv = dv/fct
e_dv = e_dv/fct

d1 = dv*tau*e_ts
d1 = d1**2

d2 = dv*ts*e_tau
d2 = d2**2

d3 = ts*tau*e_dv
d3 = d3**2

temp = (d1+d2+d3)*np.pi
dy = const*temp**0.5

print dy/5.98
print np.pi
