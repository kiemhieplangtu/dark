; File containing perameters for absorption spectra, for fitting CNM
if name eq '3C286' then begin
   zro0=0.0
	hgt0=[1,1,1]
	cen0=[-30,-15,-5]
	wid0=[1,1,1]
   look=0
   nrg=n_elements(hgt0)
   zro0yn=0
   hgt0yn=1+intarr(nrg)
   cen0yn=1+intarr(nrg)
   wid0yn=1+intarr(nrg)
 	corder='no'
endif else if name eq '4C12.50' then begin
   zro0=0.0
	hgt0=[1]
	cen0=[0]
	wid0=[1]
   look=0
   nrg=n_elements(hgt0)
   zro0yn=0
   hgt0yn=1+intarr(nrg)
   cen0yn=1+intarr(nrg)
   wid0yn=1+intarr(nrg)
 	corder='no'
endif else if name eq '3C433' then begin
   zro0=0.0
	hgt0=[1]
	cen0=[0]
	wid0=[1]
   look=0
   nrg=n_elements(hgt0)
   zro0yn=0
   hgt0yn=1+intarr(nrg)
   cen0yn=1+intarr(nrg)
   wid0yn=1+intarr(nrg)
 	corder='no'
endif 
