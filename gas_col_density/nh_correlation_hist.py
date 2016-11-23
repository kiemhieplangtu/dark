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
offsets =[0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 10.0, 15.0] # beam eg: 0.5beam, 1beam, 1.5beam ...
files = ['result/26nhi_nh_halfbeam_accurate_lb.txt',
		 'result/26nhi_nh_1beam_accurate_lb.txt',
		 'result/26nhi_nh_1beam5_accurate_lb.txt',	
		 'result/26nhi_nh_2beam_accurate_lb.txt',	
		 'result/26nhi_nh_2beam5_accurate_lb.txt',	
		 'result/26nhi_nh_3beam_accurate_lb.txt',	
		 'result/26nhi_nh_3beam5_accurate_lb.txt',	
		 'result/26nhi_nh_4beam_accurate_lb.txt',	
		 'result/26nhi_nh_4beam5_accurate_lb.txt',			 
		 'result/26nhi_nh_5beam_accurate_lb.txt',			 
		 'result/26nhi_nh_1degwidth_accurate_lb.txt',			 
		 'result/26nhi_nh_1deg5width_accurate_lb.txt'
		]

diff_fk_hl = []
diff_pl_hl = []
for i in range(0, len(files)):
	file_name = files[i]

	col_density = read_nhi_fukui_nh_planck(file_name)
	# nhi_hl  = col_density['nhi_hl']
	# nh_pl   = col_density['nh_pl']
	# nhi_fk  = col_density['nhi_fk']

	# diff_fk_hl  = col_density['fk_hl']
	# diff_pl_hl  = col_density['pl_hl']

	diff_fk_hl.extend(col_density['fk_hl'])
	diff_pl_hl.extend(col_density['pl_hl'])

	# x = [1] * len(nhi_hl)

	fk_mean = sum(diff_fk_hl) / float(len(diff_fk_hl))
	pl_mean = sum(diff_pl_hl) / float(len(diff_pl_hl))

	fk_mean = round(fk_mean, 2)
	pl_mean = round(pl_mean, 2)

	

	# plt.plot(nhi_hl, nhi_fk, 'r.')
	# plt.plot(nhi_hl, nh_pl, 'b.')

	# plt.plot(diff_fk_hl, x, 'r.')
	#plt.plot(diff_pl_hl, 'b.')	


a = np.array(diff_fk_hl)
b = np.array(diff_pl_hl)

bins=np.histogram(np.hstack((a,b)), bins=100)[1] #get the bin edges

plt.hist(a, bins, label='N(HI_Fukui)/N(HI_Heiles)')
plt.hist(b, bins, label='N(HI_Planck)/N(HI_Heiles)')

plt.xlabel('Factor - N(H)/N(HI)')
plt.ylabel('Counts')
plt.title('Histogram of Factors for all offsets')
plt.grid(True)

plt.legend(loc='upper right')
plt.show()