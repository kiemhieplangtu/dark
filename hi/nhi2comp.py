import numpy as np
##
# Read N(HI) from file  and calculate the N(HI) from the paper Heiles 2003 #
#
# params string fname Filename
#
# return N(HI) from the paper Heiles 2003 and the difference to obtained-N(HI)
# 
# Author Van Hiep
##
def read__comp_nhi(fname = "millennium_2compare.txt"):
	ret = {}

	idx     = []
	vs      = []
	ve      = []
	vs_id   = []
	ve_id   = []
	nhi_i   = []
	sources = []

	ma_nhi  = {}

	file    = open ("nhi_79src.txt","r")
	for line in file:
	    line    = line.strip()
	    columns = line.split()

	    idx.append(int(columns[0]))
	    vs.append(float(columns[1]))
	    ve.append(float(columns[2]))
	    vs_id.append(float(columns[3]))
	    ve_id.append(float(columns[4]))
	    nhi_i.append(float(columns[5]))
	    sources.append(columns[6])

	    ma_nhi[columns[6]] =[]

	file.close()

	nhi = []
	src = []

	f = open (fname,"r")
	for line in f:
	    line    = line.strip()
	    columns = line.split()

	    nhi.append(float(columns[0]))
	    src.append(columns[1])

	    ma_nhi[columns[1]].append(float(columns[0]))

	f.close()

	for i in range(len(sources)):
		ma_nhi[sources[i]] = sum(ma_nhi[sources[i]])
		print('{0}\t{1}'.format(idx[i], sources[i]))


############## MAIN ###########
read__comp_nhi()	