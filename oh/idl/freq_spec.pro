PRO freq_spec, file_id, src_id
 
fmt = 'A,F,F' 
filename = 'data/26src_radec_no_co.dat'
readcol,filename,f=fmt,name,ra0,dec0

; SRC1, 3C18 @1665.4
files = ['p1_claire_16feb2013-23jun_2013_a2770_bd0.sav', 'p1_claire_16feb2013-23jun_2013_a2770_bd1.sav', 'p1_claire_16feb2013-23jun_2013_a2770_bd2.sav', 'p1_claire_16feb2013-23jun_2013_a2770_bd3.sav']

files = ['16feb13.a2770.1_bd0.sav', '16feb13.a2770.1_bd1.sav', '16feb13.a2770.1_bd2.sav', '16feb13.a2770.1_bd3.sav']

fname = 'data/'+files[file_id]
restore, fname, /ve
ra  = 15.*hdr2info[*,3,*]
dec = hdr2info[*,4,*]
src = hdrsrcname

c  = 3.e5
bw = hdr1info[4,src_id]
f0 = hdr1info[5,src_id]
a  = bw/2047.
b  = f0-a*(2047./2)

print, 'Source: ', src[src_id], ' f0: ', f0

x = fltarr(2048)
for i=0,2047 do begin
	x[i] = a*i+b - f0
endfor

v = fltarr(2048)
for i=0,2047 do begin
	v[i] = -c*x[i]/f0
endfor

t = stkon(*,0,src_id)-stkoff(*,0,src_id)
t = stkon(*,0,src_id)

ylim_min = min(t[600:1500])-2.
ylim_max = max(t[600:1500])+2.

plot, v, t, title='Spectrum', XTITLE='v(km/s)', YTITLE='T(K)', YRANGE = [ylim_min,ylim_max], XRANGE = [-50,50]

 
END
