import numpy             as np
import matplotlib.pyplot as plt

import numpy as np

# Read and plot the difference between the paper's N(HI) and obtained-N(HI) #
#
# params string fname Filename
#
# return void
# 
# Author Van Hiep
##
def read_comp_nhi(fname = "nhi2comp_with_paper_20160316.txt"):
	ret = {}

	ret['idx']      = []
	ret['vs']       = []
	ret['ve']       = []
	ret['vs_id']    = []
	ret['ve_id']    = []

	ret['nhi_i']     = []
	ret['sources']   = []
	ret['ma_nhi']    = []
	ret['nhi_diff']  = []

	file    = open (fname,"r")
	file.readline() # comment line
	file.readline() # comment line
	for line in file:
	    line    = line.strip()
	    columns = line.split()

	    ret['idx'].append(int(columns[0]))
	    ret['vs'].append(float(columns[1]))
	    ret['ve'].append(float(columns[2]))
	    ret['vs_id'].append(float(columns[3]))
	    ret['ve_id'].append(float(columns[4]))
	    ret['nhi_i'].append(float(columns[5]))
	    ret['ma_nhi'].append(float(columns[6]))
	    ret['nhi_diff'].append(float(columns[7]))
	    ret['sources'].append(columns[8])

	file.close()

	return ret

# Read 48 sources from China #
#
# params string fname Filename
#
# return list-structured data
# 
# Author Van Hiep
##
def src_from_china(fname = "sources_from_china.txt"):

	src = []

	file = open (fname,"r")
	for line in file:
	    line    = line.strip()
	    columns = line.split()

	    src.append(columns[0])

	file.close()

	return src	


############## MAIN ###########
data = read_comp_nhi()
x    = data['idx']

cina_src  = src_from_china();

## ======= Difference ====== ##
# y    = map(abs, data['nhi_diff'])
# #y    = data['nhi_diff']
# plt.plot(x, y, 'r-')
# plt.grid()
# plt.title('difference')
# plt.ylabel('diff')
# plt.xlabel('index')
# for i,j in zip(x,y):
# 	if (j>9.0) :
# 		plt.annotate(str(j)+' ('+str(data['sources'][i])+')',xy=(i,j)) # with offset plt.annotate(str(j),xy=(i,j+0.5))

# plt.show()

## ======= N(HI) of each source ====== ##
y  = data['nhi_i'] # nhi_i # ma_nhi

plt.plot(x, y, 'r-')
plt.grid()
plt.title('Calculated N(HI)- Show sources not from China')
plt.ylabel('N(HI)')
plt.xlabel('source-index')
for i,j in zip(x,y):
	 if (data['sources'][i] not in cina_src) :
		plt.annotate(str(data['sources'][i]),xy=(i,j)) # with offset plt.annotate(str(j),xy=(i,j+0.5))

plt.show()