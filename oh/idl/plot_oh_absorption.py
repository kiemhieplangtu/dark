PRO plot_oh_absorption
 
fname = 'data/p1_claire_16feb2013-23jun_2013_a2770_bd2.sav'
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

on1 = stkon(*,0,0)-151.7
on2 = stkon(*,0,1)-149.
on3 = stkon(*,0,2)-150.

off1 = stkoff(*,0,0)-71.2
off2 = stkoff(*,0,1)-70.1
off3 = stkoff(*,0,2)-71.2

aon  = (on1+on2+on3)/3.
aoff = (off1+off3+off2)/3.

plot, v, aon-aoff, title='Spectrum', XTITLE='v(km/s)', YTITLE='T(K)', XRANGE = [-50,50]

 
END
