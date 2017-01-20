    !p.multi=[0,1,4]
    !x.margin = [8,3]
    !y.margin = [0,0]
    !y.omargin = [4,2]
    !x.omargin = [4,2]
	!p.thick=5
	!x.thick=5
	!y.thick=5
!p.charthick=2.2
!x.charsize=1
!y.charsize=1
!p.charsize=3.2
plotthick=5

!x.style=1
datacolor=cgcolor('black')
datathick=2
modelthick=2
cize=3.2
datastyle=0
xplotrange=[min(xdata),max(xdata)]

gclr = [cgcolor('red'),cgcolor('green'), cgcolor('blue'), cgcolor('cyan'), cgcolor('magenta'), cgcolor('yellow'), cgcolor('orange'),cgcolor('red'),cgcolor('magenta'),cgcolor('blue')]
gclrnames=['red','green','blue','cyan','magenta','yellow','orange','red','magenta','blue']

;;;;;;;;;;;;;;;;; Plot 1/4: EMISSION SPECTRA ;;;;;;;;;;;;;;;;;;
q=where((xdataem gt -50) and (xdataem lt 50))
yplotrange=[min(tdata[q])-0.1*(max(tdata[q])),max(tdata[q])+0.1*max(tdata[q])]

name1=name

plot, xdataem, tb_tot_fit, linestyle=0, /xstyle, /ystyle, ytitle=textoidl('T_{exp} (K)'), title=name1, yrange=yplotrange, charsize=cize, xrange=[min(xdata), max(xdata)], thick=modelthick, position=[0.12,0.675,0.95,0.95], xcharsize=0.000001
hline, 0, thick=0.5

oplot, xdataem, tdata, thick=datathick, linestyle=datastyle
oplot, xdataem, tb_wnm_tot, thick=modelthick, linestyle=4

oplot, xdataem, tb_cnm_tot, linestyle=2, thick=5

for numb=0,n_elements(hgt1)-1 do begin
	gcurv,xdata,zro1,hgt1[numb],cen1[numb],wid1[numb],tfitg
	oplot, xdata, tspin2[numb]*(1-exp(-tfitg)), thick=6, color=gclr[numb], linestyle=2
endfor

al_legend, ['EmissionT_exp (K)', 'Total Fit', 'T_B, ABS(v)','T_B, EM(v)', '1 sigma error'], lines=[datastyle,0,2,4,1], thick=[2,2,2,3,4], charsize=1.1, /top, /right

;;;;;;;;;;;;;;;;; Plot 2/4: EMISSION RESIDUALS ;;;;;;;;;;;;;;;;;;
plot, xdataem, (tdata-tb_tot_fit), position=[0.12,0.525,0.95,0.675], xrange=xplotrange, ytitle='Emission Residuals'

hline, sigmaw, linestyle=1
hline, (-1.)*sigmaw, linestyle=1

;;;;;;;;;;;;;;;;; Plot 3/4: ABSORPTION SPECTRA ;;;;;;;;;;;;;;;;;;
val=(max(taudata)-min(taudata))
yplotrange2=[min(taudata)-val, max(taudata)+val/3]

plot, xdata, taudata, xrange=xplotrange, /xstyle, /ystyle, ytitle='exp(-'+textoidl('\tau')+')', charsize=cize, thick=modelthick, position=[0.12,0.25,0.95,0.525], yrange=yplotrange2
oplot, xdata, tfit

for numb=0,n_elements(hgt1)-1 do begin
	gcurv,xdata,zro1,hgt1[numb],cen1[numb],wid1[numb],tfitg
	oplot, xdata, exp(-tfitg), thick=6, color=gclr[numb], linestyle=2
endfor

al_legend, ['VLA absorption', 'Total fit', 'CNM, WNM components','1 sigma error'], lines=[datastyle,0,2,1], thick=[datathick,modelthick,modelthick,4], charsize=1.1, /bottom, /right

;;;;;;;;;;;;;;;;; Plot 4/4: ABSORPTION RESIDUALS ;;;;;;;;;;;;;;;;;;
plot, xdata, taudata-tfit, xrange=xplotrange, position=[0.12,0.1,0.95,0.25], xtitle='Velocity (km/s)', ytitle='Absorption Residuals'
hline, sigma, linestyle=1
hline, sigma*(-1.), linestyle=1



