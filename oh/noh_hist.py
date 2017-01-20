import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import operator

from numpy    import array
from restore  import restore
from plotting import cplot

# Read Ts #
#
# params string fname Filename
#
# return dict info of Ts
# 
# Author Van Hiep
##
def read_ts(fname = 'result/noh1_src96_er.txt'):
	cols = ['idx','src','tau','v0','wid','tau_again','tex','noh']
	fmt  = ['i','s','f','f','f','f','f','f']
	data = restore(fname, 2, cols, fmt)
	return data.read()

# Read Tbg of 1666 from 408MHz #
#
# params string fname Filename
#
# return dict info of Tbg
# 
# Author Van Hiep
##
def read_tbg(fname = 'result/tbg1666from408.txt'):

	cols = ['idx','src','amp','v0','wid','ts1','er1','ts2','er2','tbg']
	fmt  = ['i','s','f','f','f','f','f','f','f','f']
	data = restore(fname, 3, cols, fmt)
	return data.read()	

#================= MAIN ========================#

ts   = read_ts()
ts1  = ts['tex']
noh1 = ts['noh']

ts   = read_ts(fname = 'result/noh2_src96_er.txt')
ts2  = ts['tex']
noh2 = ts['noh']

tbg = read_tbg()
tbg = tbg['tbg']

print tbg

noh65 = []
noh67 = []
for i in range(len(ts1)):
	if (noh1[i] > 0):
		noh65.append(noh1[i])
	if(noh2[i] > 0):
		noh67.append(noh2[i])

# Plot histogram #
size = 0.05
fig  = cplot()
trace1 = fig.hist(noh65,label='N(OH) - OH1665',autobinx=False,
                  xbins=dict(start=0.0, end=2.4, size=size),
                  opacity=1.0,
                  histtype='step',
                  marker=dict(
                    color = 'r',
                    linewidth=2                   
                    )
                 )
trace2 = fig.hist(noh67,label='N(OH) - OH1667',autobinx=False,
                  xbins=dict(start=0.0, end=4.0, size=size),
                  opacity=1.0,
                  histtype='step',
                  marker=dict(
                    color = 'k',
                    linewidth=2                   
                    )
                 )
data   = [trace1, trace2]

layout = dict(title  = 'Histogram of N(OH)',
				  title_fontsize=30,
	              grid   = True,
	              legend = dict(loc='upper right', fontsize=18),
	              xaxis  = dict(label='N(OH)/1e14',tick_size=18,fontsize=35),
	              yaxis  = dict(label='Counts',tick_size=18,fontsize=35),
	              # text   = [dict(loc=[0.2,0.4],text='a = '+str(m)+'$\pm$'+str(ea) +',  b = '+str(b)+'$\pm$'+str(eb),color='blue',fontsize=17),
	              # 			dict(loc=[0.2,0.3],text='(Available sources with the presence of OH line are shown)',color='red',fontsize=19)
	              # 		   ],
	             )
fig.iplot(data,layout)