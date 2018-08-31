pro specfilt,data,lon,lat,nr,nc
;
;  4000 km fourier filter
;
a=6.37e6
lmax=4.e6
circum=2.*!pi*a*cos(lat*!pi/180.)
wmax=circum/lmax
pindex=where(wmax mod 1 ge .5)
mindex=where(wmax mod 1 lt .5)
wint=0*wmax
wint(pindex)=1+wmax(pindex)-wmax(pindex) mod 1
wint(mindex)=wmax(mindex)-wmax(mindex) mod 1
;
for j=0,nr-1 do begin
dummy=reform(data(*,j),nc)
coef=fft(dummy,-1)
nmax=wint(j)
fcoef=coef
for n=nmax+1,nc/2 do begin
fcoef(n)=0.
fcoef(nc-n)=0.
endfor
data(*,j)=fft(fcoef,1)
endfor
return
end
