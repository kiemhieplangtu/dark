; millennium survey plot by source, velocity range & spectra type
; assumes you've used read_spec_data and have a structure called specs
; hardcodes in the fact that velocity axis is reversed
; optional output for pixel start/end points

pro specplot, data, sourcename, vrange=vrange, vpixout=vpixout

i=where(data[*].name eq sourcename)

if n_elements(vrange) ne 0 then begin
 vstart=where(min(abs(data[i].vaxis - vrange[0])) eq (abs(data[i].vaxis - vrange[0])))
 vend=where(min(abs(data[i].vaxis - vrange[1])) eq (abs(data[i].vaxis - vrange[1])))
endif else begin
 vstart=2047
 vend=0
endelse

plot,data[i].vaxis[vend:vstart],data[i].offspec[vend:vstart],title=data[i].name,xstyle=1

if n_elements(vpixout) ne 0 then vpixout=[vend,vstart]

end
