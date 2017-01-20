; File containing WNM parameters for all sources
if name eq '3C286' then begin
  	zrownm1=0.0
  	hgtwnm1=[1,1]
  	cenwnm1=[-5,-20]
   	widwnm1=[10,10]
   look=0
   nrgwnm=n_elements(hgtwnm1)
   zrownm1yn=1
   hgtwnm1yn=1+intarr(nrgwnm)
   cenwnm1yn=1+intarr(nrgwnm)
   widwnm1yn=1+intarr(nrgwnm)
   fwnm1=fwnm1
   fwnm1yn=intarr(nrgwnm)
endif else if name eq '4C12.50' then begin
  	zrownm1=0.0
  	hgtwnm1=[1]
  	cenwnm1=[0]
   	widwnm1=[10]
   look=0
   nrgwnm=n_elements(hgtwnm1)
   zrownm1yn=1
   hgtwnm1yn=1+intarr(nrgwnm)
   cenwnm1yn=1+intarr(nrgwnm)
   widwnm1yn=1+intarr(nrgwnm)
   fwnm1=fwnm1
   fwnm1yn=intarr(nrgwnm)
endif else if name eq '3C433' then begin
  	zrownm1=0.0
  	hgtwnm1=[1]
  	cenwnm1=[0]
   	widwnm1=[10]
   look=0
   nrgwnm=n_elements(hgtwnm1)
   zrownm1yn=1
   hgtwnm1yn=1+intarr(nrgwnm)
   cenwnm1yn=1+intarr(nrgwnm)
   widwnm1yn=1+intarr(nrgwnm)
   fwnm1=fwnm1
   fwnm1yn=intarr(nrgwnm)
endif
