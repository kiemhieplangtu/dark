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


# Plot the correlation between Factors and log(N(HI)) #
#
# params list files filenames
# params list offsets offsets - width of areas
#
# return Void
# 
# Author Van Hiep
##
def plot_correlation(files, offsets):

	for i in range(0, len(files)):
		if (i > 1) :
			continue

		file_name = files[i]

		col_density = read_nhi_fukui_nh_planck(file_name)
		nhi_hl      = col_density['nhi_hl']
		nh_pl       = col_density['nh_pl']
		nhi_fk      = col_density['nhi_fk']

		diff_fk_hl  = col_density['fk_hl']
		diff_pl_hl  = col_density['pl_hl']

		for k in range(0, len(nhi_hl)):
			nhi_hl[k] = np.log10(nhi_hl[k])

		# Fit and Plot #

		# m, b = np.polyfit(nhi_hl,diff_pl_hl,1)	
		# n, c = np.polyfit(nhi_hl,diff_fk_hl,1)	

		plt.plot(nhi_hl, diff_fk_hl, 'r.', label='')
		plt.plot(nhi_hl, diff_pl_hl, 'b.', label='')
		# plt.plot(nhi_hl, m*np.array(nhi_hl) + b, 'b-')
		# plt.plot(nhi_hl, n*np.array(nhi_hl) + c, 'r-')
		
	plt.xlabel('log(N(HI)_heiles)')
	plt.ylabel('Factor')
	plt.title('Histogram of Factors for all offsets')
	plt.grid(True)
	plt.xlim(0, 1.6)
	plt.ylim(0, 5)
	#plt.legend(loc='upper right')
	plt.show()

# Plot the correlation between Factors and  N(HI)_fukui N(H)_Planck N(HI)_Heiles #
#
# params list files filenames
# params list offsets offsets - width of areas
#
# return Void
# 
# Author Van Hiep
##
def plot_factors_all_offsets(files, offsets):	

	first_factor_pl = []
	first_factor_fk = []

	for i in range(0, len(files)):
		# if (i > 1) :
		# 	continue

		file_name = files[i]

		col_density = read_nhi_fukui_nh_planck(file_name)

		nhi_hl      = col_density['nhi_hl']
		diff_fk_hl  = col_density['fk_hl']
		diff_pl_hl  = col_density['pl_hl']

		if (i == 0) :
			first_factor_pl = list(diff_pl_hl)
			first_factor_fk = list(diff_fk_hl)

		for k in range(0, len(diff_pl_hl)):
			nhi_hl[k]     = np.log10(nhi_hl[k])
			diff_pl_hl[k] = diff_pl_hl[k]/first_factor_pl[k]
			diff_fk_hl[k] = diff_fk_hl[k]/first_factor_fk[k]

			if (diff_pl_hl[k] > 1.4) :
				print k, files[i], diff_pl_hl[k]

		b    = np.array(diff_fk_hl)
		a    = np.array(diff_pl_hl)
		bins = np.histogram(np.hstack((a,b)), bins=20)[1] #get the bin edges
		plt.hist(a, bins, label='Width='+str(offsets[i]*2.) + ' beams')
		#plt.hist(b, bins, label='N(HI_Planck)/N(HI_Heiles)')

		# plt.plot(nhi_hl, diff_pl_hl, 'b.', label='')


		plt.xlabel('[Factor/Factor_at_1beam]')
		plt.ylabel('Counts')
		plt.title('Histogram of Factor-Ratios for all offsets')
		plt.grid(True)

	# plt.xlim(0.6, 1.4)
	# plt.ylim(0, 5)
	plt.text(0.72, 25.1, 'Arecibo beam at 1.4GHz = 3.5\'', color='blue', fontsize=12)
	plt.legend(loc='upper right')
	plt.show()

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

plot_correlation(files, offsets)
#plot_factors_all_offsets(files, offsets)
