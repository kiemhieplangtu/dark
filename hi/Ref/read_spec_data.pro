; read downloaded millennium survey data into structure

function read_spec_data, filename, n

fmt='A,F,F,F,F,F,F,F,F' 
readcol,filename,f=fmt,name,v,offspec,dta1,dta2,d2ta1,d2ta2,d2ta3,opspec

specs = create_struct('name', '', 'vaxis', fltarr(2048), 'offspec', fltarr(2048), 'opspec', fltarr(2048))
specs = replicate(specs,n)

; Uncomment commented stuff for thing that accepts flexible spectral lengths

;line=long(0)
;endline=long(n_elements(name))-1
for i=0,n-1 do begin
 ;count=0 
 ;endloop=0
 ;while endloop eq 0 do begin
 ; line=line+1
 ; count=count+1
 ; if (line ne endline) then begin
 ;   if name[line] ne name[line+1] then endloop=1
 ; endif else endloop=1
 ;endwhile   
 s=i*2048.0
 e=s+2047.0
 ;s=line-count
 ;e=line
 specs[i].name=name[e]
 specs[i].vaxis=v[s:e]
 specs[i].offspec=0.5*offspec[s:e]
 specs[i].opspec=opspec[s:e]
 ;print,'Source ',name[line],'    Speclen',count+1
 ;line=line+1
endfor

return,specs

end
