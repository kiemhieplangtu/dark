import numpy as np

pi = np.pi
c  = 3.e5 #km/s tau spectrum: tau vs v (km/s), 1km2 = 10^10cm2

# a0 = 1.29e-11
# a1 = 7.11e-11
# a2 = 7.71e-11
# a3 = 0.94e-11

# f1 = 1665.401e6 # MHz
# f2 = 1667.358e6 # MHz

a0 = 1.302e-11
a1 = 7.177e-11
a2 = 7.778e-11
a3 = 9.496e-12

f1 = 1665.4017e6 # MHz
f2 = 1667.3589e6 # MHz

k  = 1.38064852e-23 # x10-23 m2 kg s-2 K-1
h  = 6.62607004e-34 # x10-34 m2 kg / s

# Constant for 1667 MHz#
print 'OH1667'
c2 = 8.*16.*pi*f2*f2*k/(c*c*c*a2*h)/5./1.e10

print c2 # 2.21841824609e+14

# Constant for 1665 MHz#
print 'OH1665'
c1 = 8.*16.*pi*f1*f1*k/(c*c*c*a1*h)/3./1.e10

print c1 # 3.99757843817e+14
