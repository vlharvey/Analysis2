idir='/usr34/users/ukmo/Datfiles/'
idir='' 
yearfile='files.fil'
nl=22
nr=72
nc=96
      t=fltarr(nc,nr,nl)
icount=0
openr,2,yearfile 
ntime=0
readf,2,ntime
print,ntime
daysave=fltarr(ntime)

;
uday=36
for m=0,ntime-1 do begin
ifile='' 
readf,2,ifile
ifile=strmid(ifile,0,12)
close,10
openr,10,idir+ifile,/f77
nlm1=nl-1
umax=-999999.
umin= 999999.
vmax=-999999.
vmin= 999999.
tmax=-999999.
tmin= 999999.
zmax=-999999.
zmin= 999999.
for l=0,nl-1 do begin
plev=0.
ul=fltarr(nc,nr)
vl=fltarr(nc,nr)
tl=fltarr(nc,nr+1)
zl=fltarr(nc,nr+1)
forrd,10,plev
;
; reverse order of levels
;
forrd,10,ul,vl,tl,zl
if max(ul) ge umax then umax=max(ul)
if min(ul) le umin then umin=min(ul)
if max(vl) ge vmax then vmax=max(vl)
if min(vl) le vmin then vmin=min(vl)
if max(tl) ge tmax then tmax=max(tl)
if min(tl) le tmin then tmin=min(tl)
if max(zl) ge zmax then zmax=max(zl)
if min(zl) le zmin then zmin=min(zl)
endfor
print,ifile
print,umax,vmax,tmax,zmax
print,umin,vmin,tmin,zmin
endfor
end
