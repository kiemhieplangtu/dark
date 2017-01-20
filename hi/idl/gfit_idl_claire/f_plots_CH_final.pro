pro f_plots_CH_final 

; *** All double-checked (as of 01/15/2015) *** 

; *** Method 1 (Carl) *** 
; f_plots_CH.pro 
; This code is based on /d/bip3/sstanimi/roadrunner/Min/Jesse/NHICorrection/correct_carl_17jan14.pro. 
; See Snez's email on 10/15/2014. 
; Nexp, sig_Nexp, Ncarl, f = derived from observations 
; sig_Ncarl, sig_f = derived from MC simulation 
; MC-based errors in Ncnm and Nwnm are propagated for sig_f.  
 
; CHMC_Min.pro 
; This code performs a MC simulation to estimate errors in f. 
; See Snez's email on 12/18/2014.   

; f_plots_CH_newMC.pro 
; f_plots_CH.pro is updated with CHMC_Min.pro. 
; Nexp, sig_Nexp, Ncarl, f = derived from observations 
; sig_f = derived from MC simulation
; MC-based errors in f are calculated from the beginning. 

; f_plots_CH_final.pro 
; This code is essentially f_plots_CH.pro but uses CHMC_Min.pro to calculate MC-based errors in Ncnm and Nwnm. 

sources=['NV0157+28','4C+29.05','4C+27.07','5C06.237','B20218+35','3C067',$
         '4C+34.07','NV0232+34','3C068.2','4C+28.06','4C+28.07','4C+34.09',$
         '4C+30.04','B20326+27','3C092','3C093.1','4C+26.12','B20400+25','3C108',$
         'B20411+34','4C+25.14','4C+33.10','3C131','3C132','4C+27.14','3C133']

; *** Nexp vs Ncube *** 
; Read the large GALFA-HI cube. 
; FWHM = 4' and pixel size = 1'
; Array size = [3431,1051,121] 
; Velocity range = from -100.20 km/sec to +98.55 km/sec 

filesS='/d/bip3/sstanimi/roadrunner/Min/Jesse/NHICorrection/'
cube=readfits(filesS+'GALFA_Large_Cube/Perseus_Taurus_large_cube.fits',hc)
; 1.1 = correction to match LAB data (See Section 2.2 of Stanimirovic+14)
cube=cube*1.1 
; Get velocity information. 
n=sxpar(hc,'NAXIS3')
crval=sxpar(hc,'CRVAL3')/1.e3
crpix=sxpar(hc,'CRPIX3')
cdelt=sxpar(hc,'CDELT3')/1.e3
cube_vel=fltarr(n)
for i = 0, n_elements(cube_vel)-1 do cube_vel[i]=crval+((i+1)-crpix)*cdelt

; Get source coordinates (in degrees). 
data=read_ascii(filesS+'a2644_src.list')
ra1950=15.*(data.field1[1,*]+(data.field1[2,*]/60.)+(data.field1[3,*]/3600.))
dec1950=data.field1[4,*]+(data.field1[5,*]/60.)+(data.field1[6,*]/3600.)
jprecess, ra1950, dec1950, ra2000, dec2000

; Calculate Tavg and Ncube. 
expSpec=fltarr(n_elements(sources),n)
sig_expSpec=fltarr(n_elements(sources),n)
Ncube=fltarr(n_elements(sources))
sig_Ncube=fltarr(n_elements(sources))
Ncube_vel=fltarr(n_elements(sources),2)

; Extract 72 (= (9 * 9) - (3 * 3)) spectra.  
; (9' by 9' region) - (3' by 3' region) = 
; (9 by 9 pixel region) - (3 by 3 pixel region) 
grid=72
ii=[-4,-3,-2,-1,0,1,2,3,4,$
    -4,-3,-2,-1,0,1,2,3,4,$
    -4,-3,-2,-1,0,1,2,3,4,$
    -4,-3,-2,2,3,4,$
    -4,-3,-2,2,3,4,$
    -4,-3,-2,2,3,4,$
    -4,-3,-2,-1,0,1,2,3,4,$
    -4,-3,-2,-1,0,1,2,3,4,$
    -4,-3,-2,-1,0,1,2,3,4]
jj=[4,4,4,4,4,4,4,4,4,$
    3,3,3,3,3,3,3,3,3,$
    2,2,2,2,2,2,2,2,2,$
    1,1,1,1,1,1,$
    0,0,0,0,0,0,$
    -1,-1,-1,-1,-1,-1,$
    -2,-2,-2,-2,-2,-2,-2,-2,-2,$
    -3,-3,-3,-3,-3,-3,-3,-3,-3,$
    -4,-4,-4,-4,-4,-4,-4,-4,-4]
spectra=fltarr(grid,n)
for i = 0, n_elements(sources)-1 do begin
; Get source locations. 
  adxy, hc, ra2000[i], dec2000[i], x, y

; Calculate Tavg. 
  for j = 0, grid-1 do begin
    spectra[j,*]=cube[x+ii[j],y+jj[j],*]
    expSpec[i,*]=expSpec[i,*]+cube[x+ii[j],y+jj[j],*]
  endfor 
; expSpec[i,*] = Tavg for each source
  expSpec[i,*]=expSpec[i,*]/float(grid)
  temp=expSpec[i,*]

; Calculate the uncertainty in Tavg. 
  for k = 0, grid-1 do sig_expSpec[i,*]=sig_expSpec[i,*]+(spectra[k,*]-temp)^2
; sig_expSpec[i,*] = uncertainty in Tavg for each source 
  sig_expSpec[i,*]=sqrt(sig_expSpec[i,*]/float(grid))
  
; Calculate Ncube. 
  range=where(expSpec[i,*] gt (sig_expSpec[i,*]*3.))
  Ncube_vel[i,0]=min(range)
  Ncube_vel[i,1]=max(range)
  Ncube[i]=1.823e18*abs(cube_vel[1]-cube_vel[0])*total(expSpec[i,range]) 

; Calculate the uncertainty in Ncube. 
  for k = 0, n_elements(range)-1 do sig_Ncube[i]=sig_Ncube[i]+sig_expSpec[i,range[k]]^2
  sig_Ncube[i]=1.823e18*abs(cube_vel[1]-cube_vel[0])*sqrt(sig_Ncube[i])
endfor  

; *** Correction factor ***
Nexp=fltarr(n_elements(sources))
sig_Nexp=fltarr(n_elements(sources))
Nexp_MC_mean=fltarr(n_elements(sources))
sig_Nexp_MC=fltarr(n_elements(sources))
Ncarl=fltarr(n_elements(sources))
Ncarl_MC_mean=fltarr(n_elements(sources))
sig_Ncarl_MC=fltarr(n_elements(sources))
sig_Ncarl_MC_easy=fltarr(n_elements(sources))
Fjd=fltarr(n_elements(sources)) 
Fjd_MC_mean=fltarr(n_elements(sources)) 
sig_Fjd_MC=fltarr(n_elements(sources)) 
sig_Fjd_MC_easy=fltarr(n_elements(sources))
int_tau=fltarr(n_elements(sources))
sig_int_tau=fltarr(n_elements(sources))
Nexp_vel=fltarr(n_elements(sources),2)

absorption=fltarr(n_elements(sources),1180)
sig_absorption=fltarr(n_elements(sources),1180)
emission=fltarr(n_elements(sources),1180)
sig_emission=fltarr(n_elements(sources),1180)

; xdata = LSR velocity array  
; taudata = exp(-tau) array 
; sigtaudata = exp(-tau) uncertainty array 
; cen1 = CNM central velocity array 
; wid1 = CNM FWHM array 
; hgt1 = CNM peak optical depth array 
; tspin2 = CNM spin temperature array  
; tdata = Texp profile 
; tb_tot_fit = total CNM and WNM fit (in brightness temperatures) 
; zrownm2 = ? 
; emsigmadata = Texp uncertainty profile 
; cenwnm2 = WNM central velocity array 
; widwnm2 = WNM FWHM array 
; hgtwnm2 = WNM peak brightness temperature array 

filesC='/d/leffe2/cmurray/IDL/jesse/claireplots/Correction/'
for i = 0, n_elements(sources)-1 do begin
; Read absorption data. 
  name=filesC+sources[i]+'_ABS_params_LOS_2.sav'
  restore, /ver, name 
  absorption[i,*]=-alog(taudata) 
  sig_absorption[i,*]=sigtaudata/taudata

; Calculate Tcnm profile.
  cfit=xdata*0.
  for ng = 0, n_elements(cen1)-1 do begin 
    if(wid1[ng] gt 0.) then begin
     p=hgt1[ng]*exp(-((xdata-cen1[ng])/(0.6005612*wid1[ng]))^2)
     cfit=cfit+tspin2[ng]*p
    endif
  endfor
; cfits = Tcnm profile (including all CNM components; Equation 5 of Lee+15)

; Read emission data. 
  name=filesC+sources[i]+'_EM_params_LOS_2.sav'
  restore, name, /ver 
  emission[i,*]=tdata
  sig_emission[i,*]=emsigmadata

; Calculate Ncnm. 
; check1 = to find channels with Texp larger than 3sigma 
  check1=where(tdata gt (emsigmadata*3.))
; NH_cnm = CNM column density (including all CNM components; in 1e20)
  NH_cnm=0.01823*abs(xdata[2]-xdata[1])*total(cfit[check1])

; Calculate Twnm profile. 
  wfit=xdata*0. 
  for ng = 0, n_elements(cenwnm2)-1 do begin
    if(widwnm2[ng] gt 0.) then begin 
     wfit=wfit+hgtwnm2[ng]*exp(-((xdata-cenwnm2[ng])/(0.6005612*widwnm2[ng]))^2)
    endif 
  endfor
; wfit = Twnm profile (including all WNM components; Equation 5 of Lee+15) 
; NH_wnm = WNM column density (including all WNM components; in 1e20)
  NH_wnm=0.01823*abs(xdata[2]-xdata[1])*total(wfit[check1])

; Ncarl[i] = Ntrue for each source (in 1e20)
  Ncarl[i]=NH_cnm+NH_wnm

; Nexp[i] = HI column density from Texp profile for each source (in 1e20) 
  Nexp[i]=0.01823*abs(xdata[2]-xdata[1])*total(tdata[check1])
; sig_Nexp[i] = uncertainty in Nexp for each source (in 1e20)
  sig_Nexp[i]=0.01823*abs(xdata[2]-xdata[1])*sqrt(total(emsigmadata[check1]^2))
  Nexp_vel[i,0]=min(check1)
  Nexp_vel[i,1]=max(check1)

  CHMC_Min, xdata, tdata, emsigmadata, hgt1, sighgt1, cen1, sigcen1, wid1, sigwid1, $
  tspin2, sigtspin2, hgtwnm2, sighgtwnm2, cenwnm2, sigcenwnm2, widwnm2, sigwidwnm2, $
  NH_cnm_mean, NH_cnm_sd, NH_wnm_mean, NH_wnm_sd, Ncarl_mean, Ncarl_sd, Nexp_mean, Nexp_sd, Fjd_mean, Fjd_sd  
  Nexp_MC_mean[i]=Nexp_mean
  sig_Nexp_MC[i]=Nexp_sd
  Ncarl_MC_mean[i]=Ncarl_mean
  sig_Ncarl_MC[i]=Ncarl_sd
  sig_Ncarl_MC_easy[i]=sqrt(NH_cnm_sd^2+NH_wnm_sd^2)

; Calculate the correction factor and its uncertainty. 
; Fjd[i] = correction factor for each source 
  Fjd[i]=Ncarl[i]/Nexp[i]
  Fjd_MC_mean[i]=Fjd_mean
  sig_Fjd_MC[i]=Fjd_sd
  sig_Fjd_MC_easy[i]=Fjd[i]*sqrt((sig_Ncarl_MC_easy[i]/Ncarl[i])^2+(sig_Nexp[i]/Nexp[i])^2)

; Calculate the integrated tau and its uncertainty. 
  int_tau[i]=abs(xdata[2]-xdata[1])*total(absorption[i,check1])
  sig_int_tau[i]=abs(xdata[2]-xdata[1])*sqrt(total(sig_absorption[i,check1]^2))
endfor

;save, Ncube, sig_Ncube, Ncube_vel, Nexp, sig_Nexp, Nexp_vel, Nexp_MC_mean, sig_Nexp_MC, $
;      Ncarl, Ncarl_MC_mean, sig_Ncarl_MC, sig_Ncarl_MC_easy, Fjd, Fjd_MC_mean, sig_Fjd_MC, sig_Fjd_MC_easy, $
;      int_tau, sig_int_tau, filename='f_plots_CH_final.sav' 

device, decomposed=0
tvlct, [[0],[0],[255]], 100 
tvlct, [[34],[139],[34]], 101 
tvlct, [[255],[0],[0]], 102

xlabel1 = textoidl('log_{10}(N_{exp}/10^{20} cm^{-2})')
ylabel1 = textoidl('Correction Factor f')

xlabel2 = textoidl('Integrated \tau (km s^{-1})')
ylabel2 = textoidl('Correction Factor f')

xlabel3 = textoidl('N_{cube} (10^{20} cm^{-2})') 
ylabel3 = textoidl('N_{exp} (10^{20} cm^{-2})')

set_plot, 'ps'
device, filename='f_plots_CH_final.eps', xsize=30, ysize=11, /color, /encapsulated 

!p.multi = [0,2,1]
!x.margin = [6.5,1.5]
!y.margin = [3.5,0.5]

; Correction factor vs Nexp 
plot, alog10(Nexp), Fjd, xtit=xlabel1, ytit=ylabel1, tit='!6', yrange=[0.6,2.2], ystyle=1, psym=7, $
charsize=1.6, thick=5, xthick=5, ythick=5, charthick=5
; 0.43 = 1/ln(10)
oploterror, alog10(Nexp), Fjd, (0.43*sig_Nexp/Nexp), sig_Fjd_MC_easy, psym=3, errthick=3

; Calculate 1/sigma^2 weighted mean values. 
mean_logN=fltarr(4)
sig_mean_logN=fltarr(4)
mean_Fjd1=fltarr(4)
sig_mean_Fjd1=fltarr(4)

logN=alog10(Nexp)
sig_logN=0.43*sig_Nexp/Nexp

; See http://ned.ipac.caltech.edu/level5/Leo/Stats4_5.html for the weighted mean. 
for jjj = 0, 3 do begin
  mmm=where((logN ge (min(logN)+jjj*0.2)) and (logN lt (min(logN)+(jjj+1)*0.2)))

  mean_logN[jjj]=total(logN[mmm]/sig_logN[mmm]^2)/total(1./sig_logN[mmm]^2)
  sig_mean_logN[jjj]=sqrt(1./total(1./sig_logN[mmm]^2))

  mean_Fjd1[jjj]=total(Fjd[mmm]/sig_Fjd_MC_easy[mmm]^2)/total(1./sig_Fjd_MC_easy[mmm]^2)
  sig_mean_Fjd1[jjj]=sqrt(1./total(1./sig_Fjd_MC_easy[mmm]^2))
endfor 

plotsym, 8, 1.5, color=100, /fill
oplot, mean_logN, mean_Fjd1, psym=8
oploterror, mean_logN, mean_Fjd1, sig_mean_logN, sig_mean_Fjd1, psym=3, errthick=5, errcolor=100
oplot, alog10(Nexp), Fjd, psym=7 
oploterror, alog10(Nexp), Fjd, (0.43*sig_Nexp/Nexp), sig_Fjd_MC_easy, psym=3, errthick=3
legend, pos=[1.1,2.1], psym=[8], [' Weighted mean'], charsize=1.4, charthick=5, box=0

; Perform linear fitting. 
res1=svdfit(logN, Fjd, 2, measure_errors=sig_Fjd_MC_easy, chisq=chi1, covar=cov1, sigma=sig1)

; These are the coefficients we used to rederive N(H2) (as of 10/17/2014). 
; Slope = 0.32 
; Y-intercept = 0.81 
x1=[0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8]
oplot, x1, x1*0.32+0.81, thick=6, color=101 

; Correction factor vs Integrated tau 
plot, int_tau, Fjd, xtit=xlabel2, ytit=ylabel2, tit='!6', yrange=[0.6,2.2], ystyle=1, psym=7, $
charsize=1.6, thick=5, xthick=5, ythick=5, charthick=5
oploterror, int_tau, Fjd, sig_int_tau, sig_Fjd_MC_easy, psym=3, errthick=3

; Calculate 1/sigma^2 weighted mean values. 
mean_int_tau=fltarr(5)
sig_mean_int_tau=fltarr(5)
mean_Fjd2=fltarr(5)
sig_mean_Fjd2=fltarr(5)

for jjj = 0, 4 do begin
  mmm=where((int_tau ge (min(int_tau)+jjj*3.35)) and (int_tau lt (min(int_tau)+(jjj+1)*3.35)))

  mean_int_tau[jjj]=total(int_tau[mmm]/sig_int_tau[mmm]^2)/total(1./sig_int_tau[mmm]^2)
  sig_mean_int_tau[jjj]=sqrt(1./total(1./sig_int_tau[mmm]^2))

  mean_Fjd2[jjj]=total(Fjd[mmm]/sig_Fjd_MC_easy[mmm]^2)/total(1./sig_Fjd_MC_easy[mmm]^2)
  sig_mean_Fjd2[jjj]=sqrt(1./total(1./sig_Fjd_MC_easy[mmm]^2))
endfor 

plotsym, 8, 1.5, color=100, /fill
oplot, mean_int_tau, mean_Fjd2, psym=8
oploterror, mean_int_tau, mean_Fjd2, sig_mean_int_tau, sig_mean_Fjd2, psym=3, errthick=5, errcolor=100
oplot, int_tau, Fjd, psym=7 
oploterror, int_tau, Fjd, sig_int_tau, sig_Fjd_MC_easy, psym=3, errthick=3

; Perform linear fitting. 
res2=svdfit(int_tau, Fjd, 2, measure_errors=sig_Fjd_MC_easy, chisq=chi2, covar=cov2, sigma=sig2)

x2=[0.0,5.0,10.0,15.0,20.0,25.0]
oplot, x2, x2*res2[1]+res2[0], thick=6, color=101

;save, res1, chi1, cov1, sig1, mean_logN, mean_Fjd1, sig_mean_logN, sig_mean_Fjd1, $
;      res2, chi2, cov2, sig2, mean_int_tau, mean_Fjd2, sig_mean_int_tau, sig_mean_Fjd2, filename='f_plots_CH_final_fitting.sav'

device, /close 
set_plot, 'x' 

;set_plot, 'ps' 
;device, filename='Ncube_Nexp_compare_final.eps', xsize=15, ysize=11, /color, /encapsulated

;!p.multi = [0,1,1]
;!x.margin = [6,1.5] 
;!y.margin = [3.5,0.5] 

;plot, Ncube/1.e20, Nexp, xtit=xlabel3, ytit=ylabel3, tit='!6 ', xrange=[0,35], yrange=[0,35], xstyle=1, ystyle=1, $
;psym=7, symsize=1.4, charsize=1.6, charthick=5, thick=5, xthick=5, ythick=5, /nodata 
;oplot, Ncube/1.e20, Nexp, psym=7, symsize=1.4, thick=5, color=102
;oploterror, Ncube/1.e20, Nexp, sig_Ncube/1.e20, sig_Nexp, psym=3, errthick=3, errcolor=102

;x3=findgen(50)
;y3=x3
;oplot, x3, y3, thick=5

;device, /close
;set_plot, 'x'

stop
end
