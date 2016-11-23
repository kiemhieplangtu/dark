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
def read_nhi_fukui_nh_planck(fname = 'result/26src_no_co_nhi_and_uncertainties_full.txt'):

	cols = ['idx','src','nhi_fk','err_fk','nh_pl','err_pl','nhi_hl','er_hl','err_hl','wnm','cnm','fk_hl','pl_hl','oh','nhi_thin']
	fmt  = ['i','s','f','f','f','f','f','f','f','f','f','f','f','f','f']
	data = restore(fname, 4, cols, fmt)
	return data.read()

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

# Plot histogram #
size = 0.125
fig  = cplot()
trace1 = fig.hist(diff_fk_hl,label='N(HI_Fukui)/N(HI_Heiles) - '+str(fk_mean),autobinx=False,
                  xbins=dict(start=0.0, end=4.0, size=size),
                  opacity=1.0,
                  marker=dict(
                    color = 'b',                    
                    )
                 )
trace2 = fig.hist(diff_pl_hl,label='N(HI_Planck)/N(HI_Heiles) - '+str(pl_mean),autobinx=False,
                  xbins=dict(start=0.0, end=4.0, size=size),
                  opacity=1.0,
                  marker=dict(
                    color = 'g',                    
                    )
                 )
data   = [trace1, trace2]

layout = dict(title  = 'Histogram of Factors',
				  title_fontsize=30,
	              grid   = True,
	              legend = dict(loc='upper right', fontsize=18),
	              xaxis  = dict(label='Factor - N(H)/N(HI)',tick_size=18,fontsize=35),
	              yaxis  = dict(label='Counts',tick_size=18,fontsize=35),
	              # text   = [dict(loc=[0.2,0.4],text='a = '+str(m)+'$\pm$'+str(ea) +',  b = '+str(b)+'$\pm$'+str(eb),color='blue',fontsize=17),
	              # 			dict(loc=[0.2,0.3],text='(Available sources with the presence of OH line are shown)',color='red',fontsize=19)
	              # 		   ],
	             )
fig.iplot(data,layout)