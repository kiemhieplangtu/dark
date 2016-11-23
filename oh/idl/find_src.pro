PRO test1
 
fmt = 'A,F,F' 
filename = 'data/26src_radec_no_co.dat'
readcol,filename,f=fmt,name,ra0,dec0

restore, 'data/p1_claire_16feb2013-23jun_2013_a2770_bd3.sav', /ve
ra  = 15.*hdr2info[*,3,*]
dec = hdr2info[*,4,*]

src = hdrsrcname

for i=0,171 do begin
 	for j=0,21 do begin
		for k=0,25 do begin
			if (abs(ra0[k]-ra[j,*,i]) lt 0.1) && (abs(dec0[k]-dec[j,*,i]) lt 0.1) then begin
				print,i,' ',k,' ',name[k],' ', src[i]
			endif
		endfor
	endfor
endfor
 
END