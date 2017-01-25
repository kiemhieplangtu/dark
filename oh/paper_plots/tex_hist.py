import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import operator
import matplotlib

from numpy    import array
from restore  import restore
from plotting import cplot

matplotlib.style.use('grayscale')

## Read Ts #
 #
 # params string fname Filename
 #
 # return dict info of Ts
 # 
 # Author Van Hiep ##
def read_ts(fname = ''):
  cols = ['idx','src','tau','v0','wid','tex','tex_er','noh','noh_er']
  fmt  = ['i','s','f','f','f','f','f','f','f']
  data = restore(fname, 2, cols, fmt)
  return data.read()

## Read Tbg of 1666 from 408MHz #
 #
 # params string fname Filename
 #
 # return dict info of Tbg
 # 
 # Author Van Hiep ##
def read_tbg(fname = '../result/tbg1666from408.txt'):

	cols = ['idx','src','amp','v0','wid','ts1','er1','ts2','er2','tbg']
	fmt  = ['i','s','f','f','f','f','f','f','f','f']
	data = restore(fname, 3, cols, fmt)
	return data.read()

## Plot histogram of Tex for 1665 and 1667 lines #
 #
 # params None
 # return None
 #
 # version 10/2016 
 # Author Van Hiep ##
def tex_hist():
  ts   = read_ts(fname = '../result/noh1_src96_er.txt')
  ts1  = ts['tex']
  noh1 = ts['noh']

  ts   = read_ts(fname = '../result/noh2_src96_er.txt')
  ts2  = ts['tex']
  noh2 = ts['noh']

  tbg = read_tbg()
  tbg = tbg['tbg']

  print tbg

  ts65 = []
  ts67 = []
  for i in range(len(ts1)):
    if (ts1[i] > 0):
      ts65.append(ts1[i])
    if(ts2[i] > 0):
      ts67.append(ts2[i])

  # Plot histogram #

  fig = plt.figure()
  ax1 = fig.add_subplot(111)
  kwargs = dict(histtype='stepfilled', alpha=0.999, normed=True, bins=40)

  ax1.hist(ts65, hatch='/', label='aaa', **kwargs)
  ax1.hist(ts67, hatch='x', label='bbb', **kwargs)

  plt.title('Hist', fontsize=30)
  plt.ylabel('$Ratio f = N_{HI}$/$N^*_{HI}$', fontsize=35)
  plt.xlabel('$log_{10}(N^*_{HI}/10^{20} cm^{-2}$)', fontsize=35)
  # plt.xlim(0.0, 2.0)
  # plt.ylim(-1.0, 6.0)
  plt.grid(True)
  plt.tick_params(axis='x', labelsize=18)
  plt.tick_params(axis='y', labelsize=15)
  # plt.text(0.0, 3.2, 'adsaadas', color='blue', fontsize=20)
  plt.legend(loc='upper right', fontsize=18)
  # plt.savefig("test.png",bbox_inches='tight')
  # for i in range(len(sc)):
  #   plt.annotate('('+str(sc[i])+')', xy=(xdata[i], ydata[i]), xycoords='data',
  #             xytext=(-50.,30.), textcoords='offset points',
  #             arrowprops=dict(arrowstyle="->"),fontsize=18,
  #             )
  plt.show()



  # size = 2.
  # fig  = cplot()
  # trace1 = fig.hist(ts65,label='Tex - OH1665',autobinx=False,
  #                   xbins=dict(start=0.0, end=50.0, size=size),
  #                   opacity=1.0,
  #                   histtype='step',
  #                   marker=dict(
  #                     color = 'r',
  #                     linewidth=2                   
  #                     )
  #                  )
  # trace2 = fig.hist(ts67,label='Tex - OH1667',autobinx=False,
  #                   xbins=dict(start=0.0, end=50.0, size=size),
  #                   opacity=1.0,
  #                   histtype='step',
  #                   marker=dict(
  #                     color = 'k',
  #                     linewidth=2                   
  #                     )
  #                  )
  # trace3 = fig.hist(tbg,label='Tbg of 1666 from 408',autobinx=False,
  #                   xbins=dict(start=0.0, end=6.0, size=size),
  #                   opacity=1.0,
  #                   histtype='step',
  #                   marker=dict(
  #                     color = 'g',
  #                     linewidth=2                   
  #                     )
  #                  )
  # data   = [trace1, trace2, trace3]

  # layout = dict(title  = 'Histogram of Tex and Tbg',
  #           title_fontsize=30,
  #                 grid   = True,
  #                 legend = dict(loc='upper right', fontsize=18),
  #                 xaxis  = dict(label='T(K)',tick_size=18,fontsize=35),
  #                 yaxis  = dict(label='Counts',tick_size=18,fontsize=35),
  #                 # text   = [dict(loc=[0.2,0.4],text='a = '+str(m)+'$\pm$'+str(ea) +',  b = '+str(b)+'$\pm$'+str(eb),color='blue',fontsize=17),
  #                 #       dict(loc=[0.2,0.3],text='(Available sources with the presence of OH line are shown)',color='red',fontsize=19)
  #                 #        ],
  #                )
  # fig.iplot(data,layout) 

#================= MAIN ========================#
tex_hist()
