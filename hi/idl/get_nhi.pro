; output array of optically-thin column densities for millennium
; survey off-source spectra, in units of 1e20 cm-2
; optionally takes array of velocity endpoints, otherwise user inputs
; assumes the usual structure format built by read_spec_data.pro 
; hardcoded channel width of 0.16 km/s

function get_nhi, data, vrange_array=vrange_array

nhistr=create_struct('name', '', 'vrange', fltarr(2), 'nhi', 0.0)
nhistr=replicate(nhistr,n_elements(data[*]))

if n_elements(vrange_array) eq 0 then begin
 vrange_array=fltarr(2,79)
 vrange_array[*]=-!values.f_nan
endif else begin
 if n_elements(vrange_array[0,*]) ne 79 then print, 'Wrong size vrange_array!'
endelse 

for i=0,78 do begin

 vmin=min(data[i].vaxis) & vmax=max(data[i].vaxis)

 if finite(vrange_array[0,i]) eq 0 then begin
  ok=0
  vpix=0.0

  while ok eq 0 do begin 
   ;print,'Enter velocity range (two values separated by space):'
   ;read,s,e

   window,0
   loadct,0
   specplot, data, data[i].name
   imagesysvar = {xsysvar:!x,ysysvar:!y} 

   wset,0
   print, 'click on line edges'  
   lineok=0
   while lineok eq 0 do begin  
    cursor, x1, y1, 3
    cursor, x2, y2, 3
    if (x1 lt vmin) or (x1 gt vmax) or (x2 lt vmin) or (x2 gt vmax) or (x1 eq x2) then begin
     print, 'bad points, click again'
    endif else lineok=1 
   endwhile
 
   temp=0.0
   if (x1 gt x2) then begin
    temp=x1
    x1=x2
    x2=temp
   endif   
 
   window,1
   specplot,data, data[i].name, vrange=[x1,x2], vpixout=vpix
   wset,0
   !x = imagesysvar.xsysvar & !y = imagesysvar.ysysvar
   loadct,13
   oplot,data[i].vaxis[vpix[0]:vpix[1]],data[i].offspec[vpix[0]:vpix[1]],color=90
   loadct,0
 
   av_ch_width=mean(data[0].vaxis[vpix[0]:vpix[1]-1]-data[0].vaxis[vpix[0]+1:vpix[1]])
   nhi_i=total(data[i].offspec[vpix[0]:vpix[1]])*av_ch_width*0.018   
   print,'nhi:',nhi_i,'(1e20)'
   print,'This OK? y/n?'
   yninput, answer
   if answer eq 'y' then ok=1 else ok=0
  endwhile 
 endif else begin 
  x1=vrange_array[0,i]
  x2=vrange_array[1,i]
  specplot,data, data[i].name, vrange=[x1,x2], vpixout=vpix
  av_ch_width=mean(data[0].vaxis[vpix[0]:vpix[1]-1]-data[0].vaxis[vpix[0]+1:vpix[1]])
  nhi_i=total(data[i].offspec[vpix[0]:vpix[1]])*av_ch_width*0.018    
 endelse
 nhistr[i].name=data[i].name
 nhistr[i].vrange=[x1,x2]
 nhistr[i].nhi=nhi_i
endfor

return,nhistr

end
