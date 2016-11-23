import numpy as np
import matplotlib.pyplot as plt

# millennium_specs.txt
# asu.tsv 

# data = np.genfromtxt('aaa.txt',
# 	                  names='column_1,column_2,column_3,column_4,column_5,column_6, col7, col8, col9',
# 	                  dtype=None,
# 	                  unpack=True, 
# 	                  skiprows=54)
# print data[:,0]
# print data.shape

a = np.loadtxt('millennium_specs.txt', dtype='S', skiprows=54)

name_list = list(set(a[:,0]))

print len(name_list)
print name_list
print a[:,1]
print a.shape

name = a[:,0]
x = a[:,1]
y = a[:,2]


input_name = raw_input('Please enter name: ')
filter = (name == input_name)
plt.plot(x[filter], y[filter])
plt.show()

exit()
