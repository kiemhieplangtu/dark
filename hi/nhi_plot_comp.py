import numpy             as np
import matplotlib.pyplot as plt

## Class for Style #
 # 
 # Author Van Hiep
 ##
class style:
   BOLD = '\033[1m'
   END = '\033[0m'

## Read and plot the difference between the paper's N(HI) and obtained-N(HI) #
 #
 # params string fname Filename
 #
 # return void
 #  
 # Author Van Hiep
 ##
def read_comp_nhi(fname = "result/nhi2comp_with_paper_20160316.txt"):
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

	file = open (fname,"r")
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

## plot the difference between the paper's N(HI) and obtained-N(HI) #
 #
 # params dict data data to plot
 #
 # return void
 # 
 # Author Van Hiep
 ##
def plot_comp_nhi(data):
	x    = data['idx']
	#y    = map(abs, data['nhi_diff'])
	y    = data['nhi_diff']
	plt.plot(x, y, 'r.')
	plt.grid()
	plt.title('difference')
	plt.ylabel('diff')
	plt.xlabel('index')
	for i,j in zip(x,y):
		if (j<-9.0) :
			plt.annotate(str(j)+' ('+str(data['sources'][i])+')',xy=(i,j)) # with offset plt.annotate(str(j),xy=(i,j+0.5))
	plt.show()

## plot the difference between the paper's N(HI) and obtained-N(HI) #
 #
 # params dict data data to plot
 #
 # return void
 # 
 # Author Van Hiep
 ##
def plot_comp_nhi_percentage(data):
	x    = data['idx']
	y    = []

	nhi        = data['nhi_i']
	nhi_heiles = data['ma_nhi']

	for i in range(0, len(nhi)):
		per = round(100.*np.absolute(nhi[i]-nhi_heiles[i])/nhi_heiles[i],2)
		y.append(per)
		print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t\t{6}\t\t{7}\t\t{8}'.
			format(i, data['vs'][i],data['ve'][i],data['vs_id'][i], data['ve_id'][i], nhi[i], nhi_heiles[i], per, data['sources'][i]))

	plt.plot(x, y, 'r-')
	plt.grid()
	plt.title('aaa')
	plt.ylabel('(%)')
	plt.xlabel('Sources')
	for i,j in zip(x,y):
		if (j>10.0) :
			plt.annotate(str(j)+' ('+str(data['sources'][i])+')',xy=(i,j)) # with offset plt.annotate(str(j),xy=(i,j+0.5))
	plt.show()

## plot histogram for the difference between the paper's N(HI) and obtained-N(HI) #
 #
 # params dict data data to plot
 #
 # return void
 # 
 # Author Van Hiep
 ##
def plot_hist_comp_nhi(data):
	diff       = []
	nhi        = data['nhi_i']
	nhi_heiles = data['ma_nhi']

	count = 0
	for i in range(0, len(nhi)):
		delta = 100.*np.absolute(nhi[i]-nhi_heiles[i])/nhi_heiles[i]
		diff.append(delta)
		if (delta > 10.):
			count = count + 1

	print count, 79-count, np.mean(diff)
	a = np.array(diff)

	bins = np.histogram(np.hstack((a)), bins=35)[1] #get the bin edges

	plt.hist(a, bins, label='$N_{HI}$/$N^*_{HI}$' + ', mean value = ' + str(round(np.mean(diff),1) )+ '%', histtype='step', linewidth=4)
	plt.xlabel('Percentage of $N_{HI}$/$N^*_{HI} (\%)$',fontsize=35) # style.BOLD + 'This is my text string.' + style.END
	plt.ylabel('Counts',fontsize=35)
	plt.title('Ratio $N_{HI}$/$N^*_{HI}$', fontsize=30)

	#\n

	plt.grid(True)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=18)
	# plt.text(0.21, 0.8, 'a = '+str(a)+'$\pm$'+str(ea) +',  b = '+str(b)+'$\pm$'+str(eb), color='blue', fontsize=17)

	# plt.text(55., 14., 'Arecibo beam at 1.4GHz = 3.5\'', color='blue', fontsize=18)
	# plt.text(55., 13.5, 'Mean = ' + str(round(np.mean(diff),1) ) + '%', color='blue', fontsize=18)
	plt.text(42., 2., '$N_{HI}$: Atomic hydrogen column densities\n$N^*_{HI}$: Atomic hydrogen column densities under optically thin assumption', color='k', fontsize=14)
	plt.legend(loc='upper right', fontsize=18)
	# plt.savefig("test.png",bbox_inches='tight')
	plt.show()

## plot N(HI)/N(HI)_thin vs N(HI)_thin #
 #
 # params dict data data to plot
 #
 # return void
 # 
 # Author Van Hiep
 ##
def plot_factor_vs_nhi(data):
	fact       = []
	lognhi     = []
	nhi        = data['nhi_i'] # Optically-thin assumption
	nhi_heiles = data['ma_nhi'] # N(HI) from Carl Heiles Paper

	for i in range(0, len(nhi)):
		temp   = round(nhi_heiles[i]/nhi[i], 3)
		lognhi.append(np.log10(nhi[i]))
		fact.append(temp)

	# Fit and Plot #
	params = linear_fit(lognhi,fact)
	a      = params['a']
	b      = params['b']
	ea     = params['ea']
	eb     = params['eb']

	plt.plot(lognhi, fact, 'b^', label='Ratio $f = N_{HI}$/$N^*_{HI}$', markersize=10)
	plt.plot(lognhi, a*np.array(lognhi) + b, 'r-', linewidth=4, label='Best linear fit')

	a  = round(a, 2)
	b  = round(b, 2)
	ea = round(ea, 2)
	eb = round(eb, 2)

	plt.ylabel('$Ratio f = N_{HI}$/$N^*_{HI}$', fontsize=35)
	plt.xlabel('$log_{10}(N^*_{HI}/10^{20} cm^{-2}$)', fontsize=35)
	plt.title('Correlation between Total HI column densities $N_{HI}$ and \n HI optically thin column densities $N^*_{HI}$ along 79 Millennium Survey lines-of-sight', fontsize=30)
	plt.grid(True)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=15)

	plt.text(0.21, 1.5, '$f = ['+str(a)+'\pm'+str(ea) +']\cdot log_{10}(N^*_{HI}/10^{20}) + ['+str(b)+'\pm'+str(eb)+']$', color='blue', fontsize=20)

	plt.legend(loc='upper right', fontsize=18)
	# plt.savefig("test.png",bbox_inches='tight')
	plt.show()

## plot N(HI) vs N(HI)_thin #
 #
 # params dict data data to plot
 #
 # return void
 # 
 # Author Van Hiep
 ##
def plot_nhi_vs_thin(data):
	nhi        = data['nhi_i'] # Optically-thin assumption
	nhi_heiles = data['ma_nhi'] # N(HI) from Carl Heiles Paper

	# Fit and Plot #
	params = linear_fit(nhi,nhi_heiles)
	a = params['a']
	b = params['b']
	ea = params['ea']
	eb = params['eb']

	plt.plot(nhi, nhi_heiles, 'b^', label='data', markersize=10)
	plt.plot(nhi, a*np.array(nhi) + b, 'r-', linewidth=4, label='Best linear fit')

	a = round(a, 2)
	b = round(b, 2)
	ea = round(ea, 2)
	eb = round(eb, 2)

	plt.ylabel('$N_{HI} [10^{20} cm^{-2}]$', fontsize=35)
	plt.xlabel('$N^*_{HI} [10^{20} cm^{-2}]$)', fontsize=35)
	plt.title('Correlation between Total HI column densities $N_{HI}$ and \n HI optically thin column densities $N^*_{HI}$ along 79 Millennium Survey lines-of-sight', fontsize=30)
	plt.grid(True)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=15)

	plt.text(25., 5., '$N_{HI}[10^{20} cm^{-2}] = ('+str(a)+'\pm'+str(ea) +')\cdot N^*_{HI}[10^{20} cm^{-2}] + ('+str(b)+'\pm'+str(eb)+')$', color='blue', fontsize=20)

	plt.legend(loc='upper right', fontsize=18)
	# plt.savefig("test.png",bbox_inches='tight')
	plt.show()

## linear fit #
 #
 # params x list x-data
 # params y list y-data
 #
 # return fit parameters and errors
 # 
 # Author Van Hiep ##
def linear_fit(x,y):
	ret ={}

	sxy = 0.
	sx  = 0.
	sy  = 0.
	sx2 = 0.
	n   = len(x)
	for i in range(0,n) :
		sxy = sxy + x[i]*y[i]
		sx  = sx  + x[i]
		sy  = sy  + y[i]
		sx2 = sx2 + x[i]**2

	denom = (n*sx2 - sx**2)
	a = (n*sxy - sx*sy)/denom
	b = (sx2*sy - sx*sxy)/denom

	t    = n*sx2 - sx**2
	er_a = np.sqrt(n/t) 
	er_b = np.sqrt(sx2/t) 

	chi2 = 0.
	for i in range(0,n) :
		chi2 = chi2 + (y[i]-a*x[i]-b)**2

	chi2 = np.sqrt(chi2/(n-2))
	er_a = chi2*er_a
	er_b = chi2*er_b

	ret['a']  = a
	ret['b']  = b
	ret['ea'] = er_a
	ret['eb'] = er_b

	return ret

############## MAIN ###########

data = read_comp_nhi()
#plot_comp_nhi(data)
# plot_hist_comp_nhi(data)
# plot_comp_nhi_percentage(data)
# plot_factor_vs_nhi(data)

plot_nhi_vs_thin(data)
