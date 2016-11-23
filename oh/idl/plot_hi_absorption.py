PRO plot_hi_absorption
 
fname = 'data/16feb13.a2770.1_bd0.sav'
restore, fname, /ve
ra  = 15.*hdr2info[*,3,*]
dec = hdr2info[*,4,*]
src = hdrsrcname

c  = 3.e5
bw = hdr1info[4,0]
f0 = hdr1info[5,0]
a  = bw/2047.
b  = f0-a*(2047./2)

x = fltarr(2048)
for i=0,2047 do begin
	x[i] = a*i+b - f0
endfor

v = fltarr(2048)
for i=0,2047 do begin
	v[i] = -c*x[i]/f0
endfor

aon  = stkon(*,0,0)-171.
aoff = stkoff(*,0,0)-76.2

plot, v, aon-aoff, title='Spectrum', XTITLE='v(km/s)', YTITLE='T(K)', XRANGE = [-50,50]

 
END
