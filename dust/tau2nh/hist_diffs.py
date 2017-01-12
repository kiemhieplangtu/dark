import sys
sys.path.insert(0, r'/home/vnguyen/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import operator

from numpy    import array
from restore  import restore
from plotting import cplot

# Read N(HI)_fukui N(H)_Planck N(HI)_Heiles #
#
# params string fname Filename
#
# return dict info of N(H) and N(HI)
# 
# Author Van Hiep
##
def read_nhi_fukui_nh_planck(fname = '../result/26src_no_co_nhi_and_uncertainties_full.txt'):
	cols = ['idx','src','nhi_fk','err_fk','nh_pl','err_pl','nhi_hl','err_hl','wnm','cnm','fk_hl','pl_hl','oh','nhi_thin']
	fmt  = ['i','s','f','f','f','f','f','f','f','f','f','f','f','f']
	data = restore(fname, 4, cols, fmt)
	return data.read()

#================= MAIN ========================#

info       = read_nhi_fukui_nh_planck()
rat_fk_hl  = info['fk_hl']
rat_pl_hl  = info['pl_hl']
oh         = info['oh']
src        = info['src']

fk_mean = sum(rat_fk_hl) / float(len(rat_fk_hl))
pl_mean = sum(rat_pl_hl) / float(len(rat_pl_hl))

fk_mean = round(fk_mean, 2)
pl_mean = round(pl_mean, 2)

a = np.array(rat_fk_hl)
b = np.array(rat_pl_hl)

# Plot histogram #
# size = 0.25
# fig  = cplot()
# trace1 = fig.hist(diff_fk_hl,label='N(HI_Fukui)/N(HI_Heiles) - '+str(fk_mean),autobinx=False,
#                   xbins=dict(start=0.0, end=4.0, size=size),
#                   opacity=1.0,
#                   histtype='step',
#                   marker=dict(
#                     color = 'b',                    
#                     )
#                  )
# trace2 = fig.hist(diff_pl_hl,label='N(HI_Planck)/N(HI_Heiles) - '+str(pl_mean),autobinx=False,
#                   xbins=dict(start=0.0, end=4.0, size=size),
#                   opacity=1.0,
#                   histtype='step',
#                   marker=dict(
#                     color = 'g',                    
#                     )
#                  )
# data   = [trace1, trace2]

# layout = dict(title  = 'Histogram of Factors',
# 				  title_fontsize=30,
# 	              grid   = True,
# 	              legend = dict(loc='upper right', fontsize=18),
# 	              xaxis  = dict(label='Factor - N(H)/N(HI)',tick_size=18,fontsize=35),
# 	              yaxis  = dict(label='Counts',tick_size=18,fontsize=35),
# 	              # text   = [dict(loc=[0.2,0.4],text='a = '+str(m)+'$\pm$'+str(ea) +',  b = '+str(b)+'$\pm$'+str(eb),color='blue',fontsize=17),
# 	              # 			dict(loc=[0.2,0.3],text='(Available sources with the presence of OH line are shown)',color='red',fontsize=19)
# 	              # 		   ],
# 	             )
# fig.iplot(data,layout)

nbins = 12

bins = np.histogram(np.hstack((a)), bins=nbins)[1] #get the bin edges
plt.hist(a, bins, label='$N_{H}^{Fukui}/N_{HI}, mean = '+str(fk_mean)+'$', histtype='step', linewidth=4)

print bins

bins = np.histogram(np.hstack((b)), bins=nbins)[1] #get the bin edges
plt.hist(b, bins, label='$N_{H}^{Dust}/N_{HI}, mean = '+str(pl_mean)+'$', histtype='step', linewidth=4)

print bins

plt.xlabel('$N_{H}/N_{HI}$',fontsize=35) # style.BOLD + 'This is my text string.' + style.END
plt.ylabel('Counts',fontsize=35)
plt.title('$N_{H}/N_{HI}$ along 26 sightlines without CO', fontsize=30)

plt.grid(True)
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=18)

# plt.text(55., 14., 'Arecibo beam at 1.4GHz = 3.5\'', color='blue', fontsize=18)
plt.legend(loc='upper right', fontsize=18)
# plt.savefig("test.png",bbox_inches='tight')

for i in range(len(src)):
	if (oh[i] > 0) :
		plt.annotate('('+str(src[i])+')', xy=(a[i], 0.), xycoords='data',
               xytext=(-5.,5.+a[i]**3.5), textcoords='offset points',
               arrowprops=dict(arrowstyle="->", color='b', lw=3),fontsize=15, color='b',
               )

		plt.annotate('('+str(src[i])+')', xy=(b[i], 0.), xycoords='data',
	               xytext=(-90.+a[i]**3.5,40.+10*a[i]), textcoords='offset points',
	               arrowprops=dict(arrowstyle="->", color='g', lw=3),fontsize=15, color='r',
	               )
plt.show()