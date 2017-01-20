import sys, os
import numpy             as np
import matplotlib.pyplot as plt

## Class - Read data from file to array, column by column ##
 # 
 # Using, eg:
 # cols = ['index, 'source', 'temperature']                  # Columns
 # fmt  = ['i', 's', 'f']                                    # Format of each column (eg: ['s', 'i', 'f'])
 # x    = restore('filename.txt', skip_lines=4, cols, fmt)
 # dat  = x.read()
 #
 # version 07/2016 
 # author Nguyen Van Hiep ##
class cplot:

    # Define Constants #
    line     = 0
    errbar   = 1
    hist     = 2
    bar      = 3
    contour  = 4
    box      = 5
    vline    = 6
    vlines   = 7

    solid      = 'solid'
    circl_mkr  = 'o'

    black      = 'k'
    red        = 'r'

    red_circle = 'ro'

    size_zero  = 0
    er_mk_size = 8

    lw_default = 1
    dft_nbin   = 50

    loc_deflt  = 'upper right'
    lgnd_fs    = 12

    ## Initiate function ##
     #
     # params string filename Save figure to filename
     # return void
     #
     # version 07/2016 
     # author Nguyen Van Hiep ##
    def __init__(self, filename=''):
        self.filename = filename

    ## Check x-data and y-data ##
     #
     # params list x x-axis data
     # params list y y-axis data
     # return boolean
     #
     # version 07/2016 
     # author Nguyen Van Hiep ##
    def check(self,x,y):
        if (len(x) != len(y)):
            sys.exit('X-data and Y-data do not have the same length.')
            return False

        return True

    ## Set arguments of line ##
     #
     # params dict trace Infor of a line
     # return properties of line
     #
     # version 07/2016 
     # author Nguyen Van Hiep ##
    def set_line_args(self,trace):
        # Set line-label #
        label = ''
        if (('label' in trace) and (trace['label'] != None)):
            label = trace['label']

        # Set line-color #
        color = self.black
        if (('prop' in trace) and ('color' in trace['prop'])):
            color = trace['prop']['color']

        # Set line-style #
        ls = self.solid
        if (('prop' in trace) and ('linestyle' in trace['prop'])):
            ls = trace['prop']['linestyle']

        # Set marker #
        mkr = self.circl_mkr
        if (('prop' in trace) and ('marker' in trace['prop'])):
            mkr = trace['prop']['marker']

        # Set marker-face-color #
        mfc = self.red
        if (('prop' in trace) and ('markerfacecolor' in trace['prop'])):
            mfc = trace['prop']['markerfacecolor']

        # Set markersize #
        mks = self.size_zero
        if (('prop' in trace) and ('markersize' in trace['prop'])):
            mks = trace['prop']['markersize']

        # Set line-width #
        lw = self.lw_default
        if (('prop' in trace) and ('linewidth' in trace['prop'])):
            lw = trace['prop']['linewidth']

        return label,color,ls,mkr,mfc,mks,lw

    ## set arguments for error-bar ##
     #
     # params dict trace Infor of errorbar-line
     # return properties of line
     #
     # version 07/2016 
     # author Nguyen Van Hiep ##
    def set_erbar_args(self, trace):
        # Set line-label #
        label = ''
        if (('label' in trace) and (trace['label'] != None)):
            label = trace['label']

        # Set color #
        fmt = self.red_circle
        if (('fmt' in trace['prop']) and (trace['prop']['fmt'] != None)):
            fmt = trace['prop']['fmt']

        # Set size #
        size = self.er_mk_size
        if (('markersize' in trace['prop']) and (trace['prop']['markersize'] != None)):
            size = trace['prop']['markersize']

        # Set markeredgecolor #
        mec = self.black
        if (('markeredgecolor' in trace['prop']) and (trace['prop']['markeredgecolor'] != None)):
            mec = trace['prop']['markeredgecolor']

        # Set markeredgewidth #
        mew = self.lw_default
        if (('markeredgewidth' in trace['prop']) and (trace['prop']['markeredgewidth'] != None)):
            mew = trace['prop']['markeredgewidth']

        return label,fmt,size,mec,mew

    ## Set arguments of vertical line ##
     #
     # params dict trace Infor of line
     # return properties of line
     #
     # version 07/2016 
     # author Nguyen Van Hiep ##
    def set_vline_args(self,trace):
        # Set line-label #
        label = ''
        if (('label' in trace) and (trace['label'] != None)):
            label = trace['label']

        # Set line-color #
        color = self.black
        if (('prop' in trace) and ('color' in trace['prop'])):
            color = trace['prop']['color']

        # Set line-style #
        ls = self.solid
        if (('prop' in trace) and ('linestyle' in trace['prop'])):
            ls = trace['prop']['linestyle']

        # Set marker #
        mkr = self.circl_mkr
        if (('prop' in trace) and ('marker' in trace['prop'])):
            mkr = trace['prop']['marker']

        # Set marker-face-color #
        mfc = self.red
        if (('prop' in trace) and ('markerfacecolor' in trace['prop'])):
            mfc = trace['prop']['markerfacecolor']

        # Set markersize #
        mks = self.size_zero
        if (('prop' in trace) and ('markersize' in trace['prop'])):
            mks = trace['prop']['markersize']

        # Set line-width #
        lw = self.lw_default
        if (('prop' in trace) and ('linewidth' in trace['prop'])):
            lw = trace['prop']['linewidth']

        return label,color,ls,mkr,mfc,mks,lw

    ## Create a line ##
     #
     # params list x x-data
     # params list y y-data
     # params string label Label of line
     # params dict prop Properties of line
     # return dict ret All infor of line
     #
     # version 07/2016 
     # author Nguyen Van Hiep ##
    def lines(self,x,y,label='',prop={}):
        # Check the x-data and y-data
        self.check(x,y)

        ret = {}
        ret['chart']  = self.line

        ret['x']      = x
        ret['y']      = y

        ret['label']  = label
        ret['prop']   = prop

        return ret

    ## Plot line ##
     #
     # params dict trace Infor of line
     # return void
     #
     # version 07/2016 
     # author Nguyen Van Hiep ##
    def pl_line(self, trace):
        x = trace['x']
        y = trace['y']

        # Set line-properties #
        label,color,ls,mkr,mfc,mks,lw = self.set_line_args(trace)
        # plot #
        plt.plot(x,y,label=label,color=color,linestyle=ls,marker=mkr,markerfacecolor=mfc,markersize=mks,linewidth=lw)

    ## Create a histogram ##
     #
     # params list x Data to plot
     # params string label Label of histogram
     # params boolean autobinx Auto-bin or not
     # params dict xbins Detail for binning
     # params dict marker Properties of histogram
     # params float opacity Opacities of histogram
     # return dict ret All infor of Histogram
     #
     # version 07/2016 
     # author Nguyen Van Hiep ##
    def hist(self,x,label='',autobinx=False,xbins={},marker={},opacity=0.75, histtype='bar'):
        ret = {}

        ret['chart']      = self.hist

        ret['x']          = x

        ret['label']      = label
        ret['autobinx']   = autobinx
        ret['xbins']      = xbins
        ret['marker']     = marker
        ret['opacity']    = opacity
        ret['histtype']   = histtype

        return ret

    ## Plot histogram ##
     #
     # params dict trace Infor of Histogram
     # return void
     #
     # version 07/2016 
     # author Nguyen Van Hiep ##
    def pl_hist(self,trace):
        x        = trace['x']
        label    = trace['label']
        autobinx = trace['autobinx']
        xbins    = trace['xbins']
        marker   = trace['marker']
        opacity  = trace['opacity']
        htype    = trace['histtype']

        nbins = self.dft_nbin # 50
        if (autobinx == False):
            nbins = (xbins['end']-xbins['start'])/xbins['size']

        clor = self.black 
        if (('color' in marker) and (marker['color'] != None)):
            clor = marker['color']

        lw = 1.
        if (('linewidth' in marker) and (marker['linewidth'] != None)):
            lw = marker['linewidth']

        plt.hist(x,nbins,color=clor,alpha=opacity,label=label, histtype=htype, linewidth=lw)

    ## Create a error-bar trace ##
     #
     # params list x x-data
     # params list y y-data
     # params dict err_x error-data on X-axis
     # params dict err_y error-data on Y-axis
     # params string label Label of line
     # params dict prop Properties of line
     # return dict ret All infor of Error-bar
     #
     # version 07/2016 
     # author Nguyen Van Hiep ##
    def error_bar(self,x,y,err_x={},err_y={},label='',prop={}):
        # Check the x-data and y-data
        self.check(x,y)

        ret = {}
        ret['chart'] = self.errbar

        ret['x'] = x
        ret['y'] = y

        ret['label']   = label
        ret['err_y']   = err_y
        ret['err_x']   = err_x
        ret['prop']    = prop

        return ret

    ## Plot error-bar ##
     #
     # params dict trace Infor of Error-bar
     # return void
     #
     # version 07/2016 
     # author Nguyen Van Hiep ##
    def pl_errbar(self, trace):
        x = trace['x']
        y = trace['y']

        xerr = trace['err_x']['val']
        yerr = trace['err_y']['val']

        label,fmt,size,mec,mew = self.set_erbar_args(trace)        

        plt.errorbar(x, y,yerr=yerr, xerr=xerr,fmt=fmt,label=label,markersize=size,markeredgecolor=mec,markeredgewidth=mew)

    ## Create a vline ##
     #
     # params float x x-data
     # params float ymin y-min limit
     # params float ymax y-max limit
     # params string label Label of line
     # params dict prop Properties of line
     #
     # return dict ret All info of line
     #
     # version 07/2016 
     # author Nguyen Van Hiep ##
    def vline(self,x,ymin=0,ymax=1,label='',prop={}):
        ret = {}
        ret['chart']  = self.vline

        ret['x']      = x
        ret['label']  = label
        ret['ymin']   = ymin
        ret['ymax']   = ymax
        ret['prop']   = prop

        return ret

    ## Plot vertical line ##
     #
     # params dict trace Infor of Vertical line
     # return void
     #
     # version 07/2016 
     # author Nguyen Van Hiep ##
    def pl_vline(self, trace):
        x    = trace['x']
        ymin = trace['ymin']
        ymax = trace['ymax']

        # Set line-properties #
        label,color,ls,mkr,mfc,mks,lw = self.set_vline_args(trace)
        # plot #
        plt.axvline(x,ymin=ymin,ymax=ymax,label=label,color=color,linestyle=ls,marker=mkr,markerfacecolor=mfc,markersize=mks,linewidth=lw)

    ## Create vlines ##
     #
     # params list x x-data
     # params float ymin y-min limit
     # params float ymax y-max limit
     # params string label Label of line
     # params dict prop Properties of line
     #
     # return dict ret All infor of line
     #
     # version 07/2016 
     # author Nguyen Van Hiep ##
    def vlines(self,x=[],ymin=0,ymax=1,label='',prop={}):
        ret = {}
        ret['chart']  = self.vlines

        ret['x']      = x
        ret['label']  = label
        ret['ymin']   = ymin
        ret['ymax']   = ymax
        ret['prop']   = prop

        return ret

    ## Plot vertical lines ##
     #
     # params dict trace Infor of Vertical lines
     # return void
     #
     # version 07/2016 
     # author Nguyen Van Hiep ##
    def pl_vlines(self, trace):
        x    = trace['x']
        ymin = trace['ymin']
        ymax = trace['ymax']

        # Set line-properties #
        label,color,ls,mkr,mfc,mks,lw = self.set_vline_args(trace)
        # plot vlines #
        for xc in x:
            plt.axvline(xc,ymin=ymin,ymax=ymax,label=label,color=color,linestyle=ls,marker=mkr,markerfacecolor=mfc,markersize=mks,linewidth=lw)

    ## Set layout properties of plot ##
     #
     # params dict layout Detail of layout
     # return Properties of Layout
     #
     # version 07/2016 
     # author Nguyen Van Hiep ##
    def set_layout_args(self, layout={}):
        title = ''
        if (('title' in layout) and (layout['title'] != None)):
            title = layout['title']

        tfs = ''
        if (('title_fontsize' in layout) and (layout['title_fontsize'] != None)):
            tfs = layout['title_fontsize']

        grid = False
        if (('grid' in layout) and (layout['grid'] == True)):
            grid = True

        xlbl = ''
        if (('label' in layout['xaxis']) and (layout['xaxis']['label'] != None)):
            xlbl = layout['xaxis']['label']

        xlb_fs = ''
        if (('fontsize' in layout['xaxis']) and (layout['xaxis']['fontsize'] != None)):
            xlb_fs = layout['xaxis']['fontsize']

        ylbl = ''
        if (('label' in layout['yaxis']) and (layout['yaxis']['label'] != None)):
            ylbl = layout['yaxis']['label']

        ylb_fs = ''
        if (('fontsize' in layout['yaxis']) and (layout['yaxis']['fontsize'] != None)):
            ylb_fs = layout['xaxis']['fontsize']

        #========= Tick params =============#
        xts = yts = 15
        if (('tick_size' in layout['xaxis']) and (layout['xaxis']['tick_size'] != None)):
            xts = layout['xaxis']['tick_size']

        if (('tick_size' in layout['yaxis']) and (layout['yaxis']['tick_size'] != None)):
            yts = layout['yaxis']['tick_size']

        #========= Legend params =============#
        lgd_loc = self.loc_deflt # upper right
        lgd_fs  = self.lgnd_fs   # 12
        if (('legend' in layout) and (layout['legend'] != False)):
            if(('loc' in layout['legend']) and (layout['legend']['loc'] != None)):
                lgd_loc = layout['legend']['loc']
            if(('fontsize' in layout['legend']) and (layout['legend']['fontsize'] != None)):
                lgd_fs = layout['legend']['fontsize']

        return title,tfs,grid,xlbl,xlb_fs,ylbl,ylb_fs,xts,yts,lgd_loc,lgd_fs

    ## Write text in plot ##
     #
     # params list text List of text-infor
     # return void
     #
     # version 07/2016 
     # author Nguyen Van Hiep ##
    def write_text(self,text=[]):
        for txt in text:
            plt.text(txt['loc'][0], txt['loc'][1],txt['text'],color=txt['color'],fontsize=txt['fontsize'])

    ## Write annotation ##
     #
     # params dict inf Infor of text
     # return void
     #
     # version 07/2016 
     # author Nguyen Van Hiep ##
    def annotation(self,inf={}):
        plt.annotate(inf['label'], xy=(inf['x'],inf['y']), xycoords=inf['xycoords'],
                xytext=(inf['xoff'],inf['yoff']), textcoords=inf['textcoords'],
                arrowprops=inf['arrowprops'],fontsize=inf['fontsize']
                )

    ## Save figure to file ##
     #
     # params none
     # return void
     #
     # version 07/2016 
     # author Nguyen Van Hiep ##
    def save_figure(self):
        # manager = plt.get_current_fig_manager()
        # manager.resize(*manager.window.maxsize())
        # # Save figure using [number] dots per inch
        plt.savefig(self.filename, dpi=72)

    ## Plot layout ##
     #
     # params dict layout Detail of Layout
     # return void
     #
     # version 07/2016 
     # author Nguyen Van Hiep ##
    def layout(self,layout={}):

        title,tfs,grid,xlbl,xlb_fs,ylbl,ylb_fs,xts,yts,lgd_loc,lgd_fs = self.set_layout_args(layout)

        plt.title(title,fontsize=tfs)
        plt.grid(grid)
        plt.xlabel(xlbl, fontsize=xlb_fs)
        plt.ylabel(ylbl, fontsize=ylb_fs)

        plt.tick_params(axis='x', labelsize=xts)
        plt.tick_params(axis='y', labelsize=yts)

        if (('xlim' in layout['xaxis']) and (len(layout['xaxis']['xlim']) == 2) ):
            plt.xlim(layout['xaxis']['xlim'][0],layout['xaxis']['xlim'][1])
        if (('ylim' in layout['yaxis']) and (len(layout['yaxis']['ylim']) == 2) ):
            plt.ylim(layout['yaxis']['ylim'][0],layout['yaxis']['ylim'][1])
            
        plt.legend(loc=lgd_loc, fontsize=lgd_fs)

        # Write text #
        if(('text' in layout) and (type(layout['text']) is list) and (len(layout['text']) >= 1) ):
            self.write_text(layout['text'])

    ## Plot ##
     #
     # params list data List of Traces to plot (lines or hists or error-bars)
     # params dict layout Detail of Layout
     # params bool fullscreen Fullscreen or not
     # return void
     #
     # version 07/2016 
     # author Nguyen Van Hiep ##
    def iplot(self,data=[],layout={}, fullscreen=False):

        if(fullscreen):
            plt.switch_backend('QT4Agg') #default on my system
        for trace in data:
            if(trace['chart'] == self.line):
                self.pl_line(trace)
            elif(trace['chart'] == self.hist):
                self.pl_hist(trace)
            elif(trace['chart'] == self.errbar):
                self.pl_errbar(trace)
            elif(trace['chart'] == self.vline):
                self.pl_vline(trace)
            elif(trace['chart'] == self.vlines):
                self.pl_vlines(trace)

        # set layout and plot
        self.layout(layout)
        if(fullscreen):
            mng = plt.get_current_fig_manager()
            mng.window.showMaximized()            

        plt.show()
