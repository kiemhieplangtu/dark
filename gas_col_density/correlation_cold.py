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

	ret['nhi_hl'] 	= []
	ret['nhi_warm'] = []
	ret['nhi_cold'] = []

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
	    ret['nhi_warm'].append(float(columns[7]))
	    ret['nhi_cold'].append(float(columns[8]))

	    ret['fk_hl'].append(float(columns[9]))
	    ret['pl_hl'].append(float(columns[10]))

	file.close()

	return ret

#================= MAIN ========================#
offsets =[0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 10.0, 15.0] # beam eg: 0.5beam, 1beam, 1.5beam ...
files = ['result/res_26nhi_nh_halfbeam.txt',
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

for i in range(0, len(files)):
	if (i > 0) :
		continue

	file_name   = files[i]
	col_density = read_nhi_fukui_nh_planck(file_name)

	nhi_hl  = col_density['nhi_hl']
	nhi_w   = col_density['nhi_warm']
	nhi_c   = col_density['nhi_cold']

	nh_pl   = col_density['nh_pl']
	nhi_fk  = col_density['nhi_fk']

	diff_fk_hl  = col_density['fk_hl']
	diff_pl_hl  = col_density['pl_hl']

	for k in range(0, len(nhi_hl)):
		nhi_hl[k] = np.log10(nhi_hl[k])
		nhi_w[k]  = np.log10(nhi_w[k])
		nhi_c[k]  = np.log10(nhi_c[k])

	nhi_plot = nhi_c

	m, b = np.polyfit(nhi_plot,diff_pl_hl,1)	
	n, c = np.polyfit(nhi_plot,diff_fk_hl,1)	

	plt.plot(nhi_plot, diff_fk_hl, 'r.', label='N(HI_Fukui)/N(HI)_CNM')
	plt.plot(nhi_plot, diff_pl_hl, 'b.', label='N(HI_Planck)/N(HI)_CNM')
	plt.plot(nhi_plot, m*np.array(nhi_plot) + b, 'b-')
	plt.plot(nhi_plot, n*np.array(nhi_plot) + c, 'r-')

	plt.text(0.2, 2.7, 'a='+str(round(n,2))+', b='+str(round(c,2)), color='red', fontsize=12)
	plt.text(0.2, 0.3, 'a='+str(round(m,2))+', b='+str(round(b,2)), color='blue', fontsize=12)


	plt.xlabel('log[N(HI)/1e20]_CNM')
	plt.ylabel('Factor')
	plt.title('Correlation - Factor vs N(HI)_CNM')
	plt.grid(True)


# a = np.array(diff_fk_hl)
# b = np.array(diff_pl_hl)
# bins=np.histogram(np.hstack((a,b)), bins=100)[1] #get the bin edges
# plt.hist(a, bins, label='N(HI_Fukui)/N(HI_Heiles)')
# plt.hist(b, bins, label='N(HI_Planck)/N(HI_Heiles)')

plt.xlim(0, 1.6)
plt.ylim(0, 5)
plt.legend(loc='upper right')
plt.show()