PRO fit

;---***---
;---***Code for fitting Gaussians to HI emission/absorption pairs
;---***Uses codes by Carl Heiles, modified from S.S. and A.P. versions
;---***Edited by Claire Murray 08/15
;---***---

;---Set directory where the output files will go:
dir='~/Documents/HI_absorption/Test/'

;---Running LOS iterations? If yes, simple=1. If no, simple=0
simple=1

;---Fit for the emission spectrum too? If yes, wfit=1. If no, wfit=0.
wfit=1 

;---List of sources to fit
sourcelist=['3C286']

!p.font=0
for t=0, n_elements(sourcelist)-1 do begin
	name=sourcelist[t]
	
	print, 'Fitting...'+name

	;---Restore the emission and absorption arrays
	restore, dir+name+'.sav'

	;---Re-name the absorption arrays 
	xdata=vlsr
	taudata=spec1
	sigtaudata=sigma

	;---Set the velocity range for the emission spectra (to trim Carl's data)
	xdataem=vlsrem
	m=where( (xdataem gt -100) and (xdataem lt 100) )
	;---Rename the emission arrays within this velocity range
	tdata=em_spec[m]
	emsigmadata=emsigma[m]
	xdataem=xdataem[m]

	;---Retrieve initial Gaussian guesses for the absorption spectrum
	@cnm.pro

	;---Fit these guesses
	gcurv, xdata, zro0, hgt0, cen0, wid0, tau0
	tfit0 = exp(-tau0)


	gfitflex_exp, -1, xdata, taudata, [0,n_elements(taudata)-1], $
		zro0, hgt0, cen0, wid0, $      
		zro0yn, hgt0yn, cen0yn, wid0yn, $
		tfit, sigma, zro1, hgt1, cen1, wid1, $       
		sigzro1, sighgt1, sigcen1, sigwid1, cov, problem
	if (problem ne 0) then print, 'PROBLEM!!! = ', problem

	;---Retrieve initial Gaussian guesses for the emission spectrum
	@wnm.pro

	;---Start Tspin guesses at 30 and hold all absorption line fits fixed
	fwnm1=intarr(nrgwnm)
	order1=indgen(nrg)
	tspin1 = 30.*(fltarr(nrg)+1.)
	tspin1yn=1+intarr(nrg)

	zrocnm1 = 0. & hgtcnm1 = hgt1 & cencnm1 = cen1
	widcnm1 = wid1
	zrocnm1yn = 0
	hgtcnm1yn = 0+intarr(nrg)
	cencnm1yn = 0+intarr(nrg)
	widcnm1yn = 0+intarr(nrg)

	look=-1
	xindxrange=[0,n_elements(xdataem)-1]

	;---Parameters within tbgfitflex_exp.pro, sets number of loops (nloopmax)
	;---and the fractional change in each parameter per loop iteration
	nloopmax=100
	halfasseduse=0.2

	;---Compute Tsky at the source position [**predict_sync method currently not working, 
        ;---so I just set a generic value here, will fix in the future]
	;@galactic_coordinates.pro
	;print, 'gl gb', gl, ' ', gb
	;tsky=predict_sync(gl,gb, nu=1.4, /interp)+2.725 
	tsky=2.8

	tdata=tdata+tsky
	if simple eq 1 then begin
		tbgfitflex_exp, look, xdataem, tdata, xindxrange, $
			zrocnm1, hgtcnm1, cencnm1, widcnm1, tspin1, order1, $
			zrocnm1yn, hgtcnm1yn, cencnm1yn, widcnm1yn, tspin1yn, $
			cont, hgtwnm1, cenwnm1, widwnm1, fwnm1, $
			1, hgtwnm1yn, cenwnm1yn, widwnm1yn, fwnm1yn, $
			tfita, sigmaw, $
			zrocnm2, hgtcnm2, cencnm2, widcnm2, tspin2, $
			sigzrocnm2, sighgtcnm2, sigcencnm2, sigwidcnm2, sigtspin2, $
			zrownm2, hgtwnm2, cenwnm2, widwnm2, fwnm2, $
			sigzrownm2, sighgtwnm2, sigcenwnm2, sigwidwnm2, sigfwnm2, $
			cov, problem, nloop,$
        		tb_cont=tb_cont, tb_wnm_tot=tb_wnm_tot, tb_cnm_tot=tb_cnm_tot, $
        		exp_tau_sum=exp_tau_sum, nloopmax=nloopmax, halfasseduse=halfasseduse

		tb_tot_fit= tb_cnm_tot+tb_wnm_tot+tb_cont-tsky
		tb_tot_fit=tfita-tsky
		tdata=tdata-tsky

	endif else begin
		@losfit.pro
	endelse

	;;;;Calculating column densities:;;;;
	NH_cnm =  1.064467 * 0.0183 * 1. * hgt1 * wid1 * tspin2
	pom=sighgt1^2*(0.0194797*wid1*tspin2)^2 + $
		sigwid1^2*(0.0194797*hgt1*tspin2)^2 +$
		sigtspin2*(0.0194797*hgt1*wid1)^2
	delta_NH_cnm=sqrt(pom)

	print, 'NH_cnm: ', NH_cnm, ' ', delta_NH_cnm

	NH_wnm = 1.064467 * 0.0183 * 1. * hgtwnm2 * widwnm2
	pomw=NH_wnm^2*((sighgtwnm2/hgtwnm2)^2 + (sigwidwnm2/widwnm2)^2)
	delta_NH_wnm=sqrt(pomw)

	print, 'NH_wnm: ', NH_wnm, ' ', delta_NH_wnm

	;;;solving for Tkmax
	Tkmax=21.855*(wid1)^2
	sigTkmax=(sigwid1/wid1)*2*Tkmax
	Tkmaxw=21.855*(widwnm2)^2

	print, 'Sigma for the fit: ', sigmaw
	for nr=0,nrg-1 do begin 
		print, name, ' ', nr, ' Tkmax:', Tkmax[nr], ' ', sigTkmax[nr], ' Tspin:', tspin2[nr], ' ', sigtspin2[nr]
	endfor
	print
	print, cenwnm2
	print, widwnm2
	print, hgtwnm2 
	print

	;;;Saving parameters from the fits
	if simple eq 1 then begin
		save,filename=dir+name+'_ABS_params.sav',$
		xdata, taudata, tfit, hgt1, sighgt1, cen1, sigcen1, wid1, sigwid1, zro1, sigma, Tkmax, tspin2, sigtspin2, NH_cnm, delta_NH_cnm, name, order1, tsky

		save,filename=dir+name+'_EM_params.sav',$
		xdataem, tdata, tb_tot_fit, tb_cnm_tot, tb_wnm_tot, hgtwnm2, sighgtwnm2, cenwnm2, sigcenwnm2, widwnm2, sigwidwnm2, zrownm2, sigmaw, Tkmaxw, NH_wnm, delta_NH_wnm, fwnm2, cov, name  
	endif 

	if simple eq 0 then begin
		save,filename=dir+name+'_ABS_params_LOS.sav',$
		xdata, taudata, tfit, hgt1, sighgt1, cen1, sigcen1, wid1, sigwid1, zro1, sigma, Tkmax, tspin2, sigtspin2, NH_cnm, name, order1, tsky

		save,filename=dir+name+'_EM_params_LOS.sav',$
		xdataem, tdata, tb_tot_fit, tb_cnm_tot, tb_wnm_tot, hgtwnm2, sighgtwnm2, cenwnm2, sigcenwnm2, widwnm2, sigwidwnm2, zrownm2, sigmaw, Tkmaxw, NH_wnm, fwnm1, cov, name  
	endif 

	@plotter.pro
	;endelse
endfor


end
