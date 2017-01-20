nrgcnm=nrg
@orders.pro

f_pom = fltarr(nrgwnm,3^nrgwnm)
sigmas= fltarr(3^nrgwnm*factorial(nrgcnm))

Tspins=fltarr(nrgcnm,3^nrgwnm*factorial(nrgcnm))
Tspins_err=fltarr(nrgcnm,3^nrgwnm*factorial(nrgcnm))

sigmatemp=fltarr(factorial(nrgcnm))
orders_all=fltarr(nrgcnm,3^nrgwnm*factorial(nrgcnm))

for k = 0, nrgwnm-1 do f_pom[k,*] = indgen(3^nrgwnm)/3^k mod 3
f_pom=f_pom*0.5

for j=0, 3^nrgwnm-1 do begin
  fwnm1=f_pom[*,j]
	print, 'Starting orders loop: ', j 
	loop_progress, j, 0, nrgwnm-1

for oval=0, (factorial(nrgcnm)-1) do begin
	order1=orders[*,oval]

look=-1
tbgfitflex_exp, look, xdataem, tdata, [0,n_elements(xdataem)-1], $
	zrocnm1, hgtcnm1, cencnm1, widcnm1, tspin1, order1, $
	zrocnm1yn, hgtcnm1yn, cencnm1yn, widcnm1yn, tspin1yn, $
	zrownm1, hgtwnm1, cenwnm1, widwnm1, fwnm1, $
	zrownm1yn, hgtwnm1yn, cenwnm1yn, widwnm1yn, fwnm1yn, $
	tfita, sigmaw, $
	zrocnm2, hgtcnm2, cencnm2, widcnm2, tspin2, $
	sigzrocnm2, sighgtcnm2, sigcencnm2, sigwidcnm2, sigtspin2, $
	zrownm2, hgtwnm2, cenwnm2, widwnm2, fwnm2, $
	sigzrownm2, sighgtwnm2, sigcenwnm2, sigwidwnm2, sigfwnm2, $
	cov, problem, nloop,$
        tb_cont=tb_cont, tb_wnm_tot=tb_wnm_tot, tb_cnm_tot=tb_cnm_tot, $
        exp_tau_sum=exp_tau_sum, nloopmax=nloopmax, halfasseduse=halfasseduse

	print, 'fwnm is: ', f_pom[*,j]
	print, 'order1 is:', order1
	print, 'tspin is: ', tspin2, ' p/m:', sigtspin2

  sigmas[j*factorial(nrgcnm)+oval]=sigmaw
  Tspins[*,j*factorial(nrgcnm)+oval]= tspin2
  Tspins_err[*,j*factorial(nrgcnm)+oval]= sigtspin2
 
endfor ; for oval

endfor ; for j   

tspin_final=fltarr(nrgcnm)
tspin_err_final=fltarr(nrgcnm)
F=(factorial(nrgcnm))*(3^nrgwnm)
for j=0, nrgcnm-1 do begin

	w=(1/sigmas)^2

	tspin_val=total(w*Tspins[j,*])/total(w)

	tspin_err=(total(w*((Tspins[j,*]-tspin_val)^2+Tspins_err[j,*]^2))/total(w))*(F/(F-1))
	tspin_err=sqrt(Tspin_err)

	tspin_final[j]=tspin_val
	tspin_err_final[j]=tspin_err
endfor

print, 'Final:'
for jj=0, nrgcnm-1 do print, $
'Tspin is:', tspin_final[jj], '+/-', tspin_err_final[jj]

m=where(sigmas eq min(sigmas))
print, 'where sigmas is min:', m

fwnm1=f_pom[*,(floor(m[0]/factorial(nrgcnm)))]
order1=orders[*,(m[0] mod factorial(nrgcnm))]
fwnm1yn=intarr(nrgwnm)

print, 'fwnm: ', fwnm1
print, 'order1: ', order1


look=-1
tbgfitflex_exp, look, xdataem, tdata, [0,n_elements(xdataem)-1], $
	zrocnm1, hgtcnm1, cencnm1, widcnm1, tspin1, order1, $
	zrocnm1yn, hgtcnm1yn, cencnm1yn, widcnm1yn, tspin1yn, $
	zrownm1, hgtwnm1, cenwnm1, widwnm1, fwnm1, $
	zrownm1yn, hgtwnm1yn, cenwnm1yn, widwnm1yn, fwnm1yn, $
	tfita, sigmaw, $
	zrocnm2, hgtcnm2, cencnm2, widcnm2, tspin2, $
	sigzrocnm2, sighgtcnm2, sigcencnm2, sigwidcnm2, sigtspin2, $
	zrownm2, hgtwnm2, cenwnm2, widwnm2, fwnm2, $
	sigzrownm2, sighgtwnm2, sigcenwnm2, sigwidwnm2, sigfwnm2, $
	cov, problem, nloop,$
        tb_cont=tb_cont, tb_wnm_tot=tb_wnm_tot, tb_cnm_tot=tb_cnm_tot, $
        exp_tau_sum=exp_tau_sum, nloopmax=nloopmax, halfasseduse=halfasseduse

tb_tot_fit=tfita-tsky
tdata=tdata-tsky


save, filename=dir+name+'_losparamaters.sav', $
Tspins, Tspins_err, sigmas, f_pom, orders, NH_cnm, delta_NH_cnm, NH_wnm, delta_NH_wnm
