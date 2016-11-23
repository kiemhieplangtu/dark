import matplotlib.pyplot as plt
import numpy             as np

from numpy  import array
import operator

# Read N(HI)_fukui N(H)_Planck N(HI)_Heiles #
#
# params string fname Filename
#
# return dict info of N(H) and N(HI)
# 
# Author Van Hiep
##
def read_nhi_fukui_nh_planck(fname = 'result/26nhi_nh_halfbeam_accurate_lb.txt'):
	ret = {}

	ret['src']    = []
	ret['nhi_fk'] = []
	ret['err_fk'] = []

	ret['nh_pl']  = []
	ret['err_pl'] = []

	ret['nhi_hl'] = []

	ret['fk_hl']  = []
	ret['pl_hl']  = []

	file = open (fname,'r')
	file.readline() # comment line
	file.readline() # comment line
	file.readline() # comment line
	file.readline() # comment line

	for line in file:
	    line    = line.strip()
	    columns = line.split()

	    ret['src'].append(columns[1])
	    ret['nhi_fk'].append(float(columns[2]))
	    ret['err_fk'].append(float(columns[3]))

	    ret['nh_pl'].append(float(columns[4]))
	    ret['err_pl'].append(float(columns[5]))

	    ret['nhi_hl'].append(float(columns[6]))

	    ret['fk_hl'].append(float(columns[7]))
	    ret['pl_hl'].append(float(columns[8]))

	file.close()

	return ret

#================= MAIN ========================#

col_density = read_nhi_fukui_nh_planck()
diff_fk_hl  = col_density['fk_hl']
diff_pl_hl  = col_density['pl_hl']

fk_mean = sum(diff_fk_hl) / float(len(diff_fk_hl))
pl_mean = sum(diff_pl_hl) / float(len(diff_pl_hl))

fk_mean = round(fk_mean, 2)
pl_mean = round(pl_mean, 2)

a = np.array(diff_fk_hl)
b = np.array(diff_pl_hl)

bins=np.histogram(np.hstack((a,b)), bins=32)[1] #get the bin edges

plt.hist(a, bins, label='N(HI_Fukui)/N(HI_Heiles) - '+str(fk_mean))
plt.hist(b, bins, label='N(HI_Planck)/N(HI_Heiles) - '+str(pl_mean))

plt.xlabel('Factor - N(H)/N(HI)')
plt.ylabel('Counts')
plt.title('Histogram of Factors')
plt.grid(True)

plt.legend(loc='upper right')
plt.show()