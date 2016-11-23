import matplotlib.pyplot as plt
import numpy             as np

from numpy  import array
import operator

# Read glon & glat #
#
# params string fname Filename
#
# return dict info of src
# 
# Author Van Hiep
##
def read_lb(fname = 'data/79src_radec_lb.txt'):
	ret = {}

	ret['l0']   = [] # 0 - no CO
	ret['b0']   = []
	ret['l1']   = [] # 1 - with CO
	ret['b1']   = []
	ret['l2']   = [] # 2 - Heiles only, not considered.
	ret['b2']   = []

	ret['src']  = []
	ret['ra']   = []
	ret['dec']  = []

	file = open (fname,'r')
	file.readline() # comment line
	file.readline() # comment line
	file.readline() # comment line

	for line in file:
	    line    = line.strip()
	    columns = line.split()

	    l = float(columns[0])
	    if (l > 180.) :
	    	l = l - 360.

	    typ = int(columns[5])
	    if (typ == 0) :
	    	ret['l0'].append(l)
	    	ret['b0'].append(float(columns[1]))
	    elif (typ == 1) :
	    	ret['l1'].append(l)
	    	ret['b1'].append(float(columns[1]))
	    else :
	    	ret['l2'].append(l)
	    	ret['b2'].append(float(columns[1]))

	    ret['src'].append(columns[2])
	    ret['ra'].append(float(columns[3]))
	    ret['dec'].append(float(columns[4]))

	file.close()

	return ret

# plot sorces in MW map #
#
# params dict data data to plot
#
# return void
# 
# Author Van Hiep
##
def map_src_in_mw(data):

	l0 = data['l0'] # 
	b0 = data['b0'] # 

	l1 = data['l1'] # 
	b1 = data['b1'] # 

	l2 = data['l2'] # 
	b2 = data['b2'] # 

	t = np.arange(-180., 180.0, 0.1)

	plt.plot(l0,b0, 'r^', label='without CO', markersize=10)
	plt.plot(l1,b1, 'b^', label='with CO', markersize=10)
	plt.plot(l2,b2, 'ko', label='no CO data', markersize=6)
	plt.plot(t,0.*t, 'k--', label='', linewidth=0.6)

	plt.xlabel('Galactic longitude $(^{o})$', fontsize=35)
	plt.ylabel('Galactic latitude $(^{o})$', fontsize=35)
	plt.title('Locations of 79 Sources', fontsize=30)
	plt.xlim(-180., 180.)
	plt.ylim(-90., 90.)
	# plt.grid(True)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=18)

	# plt.text(0.21, 1.31, 'Arecibo beam at 1.4GHz = 3.5\'', color='blue', fontsize=12)
	# plt.text(0.21, 0.92, 'a = '+str(a)+'  b = '+str(b), color='blue', fontsize=12)

	plt.legend(loc='upper right', fontsize=18)
	plt.show()

#================= MAIN ========================#

data = read_lb()
map_src_in_mw(data)