import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import operator

from numpy    import array
from restore  import restore
from plotting import cplot

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
def read_tbg(fname = 'result/tbg1666from408.txt'):

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
  ts   = read_ts(fname = 'result/noh1_src96_er.txt')
  ts1  = ts['tex']
  noh1 = ts['noh']

  ts   = read_ts(fname = 'result/noh2_src96_er.txt')
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
  size = 2.
  fig  = cplot()
  trace1 = fig.hist(ts65,label='Tex - OH1665',autobinx=False,
                    xbins=dict(start=0.0, end=50.0, size=size),
                    opacity=1.0,
                    histtype='step',
                    marker=dict(
                      color = 'r',
                      linewidth=2                   
                      )
                   )
  trace2 = fig.hist(ts67,label='Tex - OH1667',autobinx=False,
                    xbins=dict(start=0.0, end=50.0, size=size),
                    opacity=1.0,
                    histtype='step',
                    marker=dict(
                      color = 'k',
                      linewidth=2                   
                      )
                   )
  trace3 = fig.hist(tbg,label='Tbg of 1666 from 408',autobinx=False,
                    xbins=dict(start=0.0, end=6.0, size=size),
                    opacity=1.0,
                    histtype='step',
                    marker=dict(
                      color = 'g',
                      linewidth=2                   
                      )
                   )
  data   = [trace1, trace2, trace3]

  layout = dict(title  = 'Histogram of Tex and Tbg',
            title_fontsize=30,
                  grid   = True,
                  legend = dict(loc='upper right', fontsize=18),
                  xaxis  = dict(label='T(K)',tick_size=18,fontsize=35),
                  yaxis  = dict(label='Counts',tick_size=18,fontsize=35),
                  # text   = [dict(loc=[0.2,0.4],text='a = '+str(m)+'$\pm$'+str(ea) +',  b = '+str(b)+'$\pm$'+str(eb),color='blue',fontsize=17),
                  #       dict(loc=[0.2,0.3],text='(Available sources with the presence of OH line are shown)',color='red',fontsize=19)
                  #        ],
                 )
  fig.iplot(data,layout) 

## Plot difference between Tex1665 and Tex 1667 #
 #
 # params None
 # return None
 #
 # version 10/2016 
 # Author Van Hiep ##
def component_tex_diff():
  ts     = read_ts(fname = 'result/noh1_src96_er.txt')
  ts1    = ts['tex']
  noh1   = ts['noh']
  ts1_er = ts['tex_er']

  indx   = ts['idx']
  print len(set(indx))

  idx = []
  xer = []
  for i in range(len(indx)):
    idx.append(i)
    xer.append(0.)

  ts     = read_ts(fname = 'result/noh2_src96_er.txt')
  ts2    = ts['tex']
  noh2   = ts['noh']
  ts2_er = ts['tex_er']

  tbg = read_tbg()
  tbg = tbg['tbg']

  # plt.plot(idx,ts1,'b.')
  # plt.plot(idx,ts2,'r.')
  # plt.grid()
  # plt.show()

  fig    = cplot()
  trace1 = fig.error_bar(idx,ts1,label='Tex1665',
      err_x=dict(val=xer),
      err_y=dict(val=ts1_er),      
      prop=dict(fmt='bo',
        markersize=8,
        markeredgecolor='b',
        markeredgewidth=1)
      )

  trace2 = fig.error_bar(idx,ts2,label='Tex1667',
      err_x=dict(val=xer),
      err_y=dict(val=ts2_er),      
      prop=dict(fmt='ro',
        markersize=8,
        markeredgecolor='r',
        markeredgewidth=1)
      )

  data   = [trace1, trace2]
  layout = dict(title  = 'Tex1667 and Tex1665',
          title_fontsize=30,
                grid   = True,
                legend = dict(loc='upper right', fontsize=18),
                xaxis  = dict(label='Index',tick_size=18,fontsize=35,xlim=[-1.,50.]),
                yaxis  = dict(label='Tex (K)',tick_size=18,fontsize=35,ylim=[-10.,35.]),
                # text   = [dict(loc=[0.2,0.4],text='a = '+str(m)+'$\pm$'+str(ea) +',  b = '+str(b)+'$\pm$'+str(eb),color='blue',fontsize=17),
                #       dict(loc=[0.2,0.3],text='(Available sources with the presence of OH line are shown)',color='red',fontsize=19)
                #        ],
               )

  # for i in range(len(src)):
  #   fig.annotation(dict(label='('+src[i]+')',
  #     x=n_oh2[i],y=ratio[i],
  #     fontsize=12,
  #     xycoords='data',
  #     xoff=-50.,yoff=20.,
  #     textcoords='offset points',
  #     arrowprops=dict(arrowstyle="->") )
  #     )

  fig.iplot(data,layout)

  ts_diff = np.abs(np.asarray(ts2) - np.asarray(ts1))
  print ts_diff.mean()

  #------------------------------#
  # Plot ts1-ts2 #
  trace3 = fig.error_bar(idx,ts_diff,label='|$Tex_{1667}-Tex{1665}$|',
      err_x=dict(val=xer),
      err_y=dict(val= np.sqrt(np.asarray(ts1_er)**2 + np.asarray(ts2_er)**2 ) ) ,      
      prop=dict(fmt='bo',
        markersize=8,
        markeredgecolor='b',
        markeredgewidth=1)
      )

  data   = [trace3]
  layout = dict(title  = 'Tex anomalies |$Tex_{1667}-Tex{1665}$|',
          title_fontsize=30,
                grid   = True,
                legend = dict(loc='upper right', fontsize=18),
                xaxis  = dict(label='Index',tick_size=18,fontsize=35,xlim=[-1.,50.]),
                yaxis  = dict(label='Tex (K)',tick_size=18,fontsize=35,ylim=[-20.,20.]),
                # text   = [dict(loc=[0.2,0.4],text='a = '+str(m)+'$\pm$'+str(ea) +',  b = '+str(b)+'$\pm$'+str(eb),color='blue',fontsize=17),
                #       dict(loc=[0.2,0.3],text='(Available sources with the presence of OH line are shown)',color='red',fontsize=19)
                #        ],
               )

  # for i in range(len(src)):
  #   fig.annotation(dict(label='('+src[i]+')',
  #     x=n_oh2[i],y=ratio[i],
  #     fontsize=12,
  #     xycoords='data',
  #     xoff=-50.,yoff=20.,
  #     textcoords='offset points',
  #     arrowprops=dict(arrowstyle="->") )
  #     )

  fig.iplot(data,layout)



  #------------------------------#
  # Plot histogram ts1/ts2 #
  ts1    = np.asarray(ts1, dtype=np.float32)
  ts2    = np.asarray(ts2, dtype=np.float32)
  size   = 0.5
  trace1 = fig.hist(ts1/ts2,label='Tex1665/Tex1667',autobinx=False,
                    xbins=dict(start=0.0, end=10.0, size=size),
                    opacity=1.0,
                    histtype='step',
                    marker=dict(
                      color = 'r',
                      linewidth=2                   
                      )
                   )
  data   = [trace1]

  layout = dict(title  = 'Histogram of Tex1665/Tex1667',
            title_fontsize=30,
                  grid   = True,
                  legend = dict(loc='upper right', fontsize=18),
                  xaxis  = dict(label='Tex1665/Tex1667',tick_size=18,fontsize=35),
                  yaxis  = dict(label='Counts',tick_size=18,fontsize=35),
                  # text   = [dict(loc=[0.2,0.4],text='a = '+str(m)+'$\pm$'+str(ea) +',  b = '+str(b)+'$\pm$'+str(eb),color='blue',fontsize=17),
                  #       dict(loc=[0.2,0.3],text='(Available sources with the presence of OH line are shown)',color='red',fontsize=19)
                  #        ],
                 )
  fig.iplot(data,layout)

  #------------------------------#
  # Plot histogram ts1-ts2 #
  ts1    = np.asarray(ts1, dtype=np.float32)
  ts2    = np.asarray(ts2, dtype=np.float32)
  size   = 0.5
  trace1 = fig.hist(np.abs(-ts1+ts2),label='|$T_{ex_{1667}}-T_{ex_{1665}}$|, Mean value: ' + str(round(ts_diff.mean(),1 ))+'K',
                    xbins=dict(start=0.0, end=10.0, size=size),
                    opacity=1.0,
                    histtype='step',
                    marker=dict(
                      color = 'r',
                      linewidth=4                   
                      )
                   )
  data   = [trace1]

  layout = dict(title  = 'OH main-line Excitation temperature anomalies |$T_{ex_{1667}}-T_{ex_{1665}}$|',
            title_fontsize=30,
                  grid   = True,
                  legend = dict(loc='upper right', fontsize=18),
                  xaxis  = dict(label='|$T_{ex_{1667}}-T_{ex_{1665}}$| [K]',tick_size=18,fontsize=35),
                  yaxis  = dict(label='Counts',tick_size=18,fontsize=35),
                  # text   = [dict(loc=[17.07,18.07],text='Mean value: ' + str(round(ts_diff.mean(),1 )),color='blue',fontsize=17),
                        # dict(loc=[0.2,0.3],text='(Available sources with the presence of OH line are shown)',color='red',fontsize=19)
                         # ],
                 )
  fig.iplot(data,layout) 

#================= MAIN ========================#
tex_hist()
# component_tex_diff()
