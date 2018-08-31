
; Plot markers.  cyclones = 1;  anticyclones = -1

@stddat
@kgmt
@ckday
@kdate

num=0L
nlvls=20
loadct,38
mcolor=!p.color
icolmax=byte(!p.color)
icmm1=icolmax-1B
icmm2=icolmax-2B
col1=1+indgen(nlvls)*icolmax/nlvls
!P.FONT=0
SETPLOT='ps'
read,'setplot',setplot
nxdim=750
nydim=750
xorig=[0.10,0.10]
yorig=[0.60,0.20]
xlen=0.8
ylen=0.3
cbaryoff=0.03
cbarydel=0.01
mon=strarr(12)*4
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
lstmn=0l
lstdy=0l
lstyr=0l
ledmn=0l
leddy=0l
ledyr=0l
thlev=0l
lstday=0l
ledday=0l
date=''
print, ' '
print, '      UKMO Version '
print, ' '
read,'Enter date (month, day, year) ',lstmn,lstdy,lstyr
read,'Enter ending date (month, day, year) ',ledmn,leddy,ledyr
print, ' '
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1991 then stop,'Year out of range '
if ledyr lt 1991 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '

   lc=icolmax
if setplot eq 'ps' then begin
   lc=0
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   wdelete,1
   set_plot,'ps'
   device,/color,/landscape,bits=8,filename='ukmo_marks.ps'
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
   xsize=xsize,ysize=ysize
endif

; Compute initial Julian date
iyr = lstyr
idy = lstdy
imn = lstmn
z=kgmt(imn,idy,iyr,iday)
iday = iday - 1

; loop over days
jump: iday = iday + 1
kdate,float(iday),iyr,imn,idy
ckday,iday,iyr
iyr1=iyr-1900
print,imn,idy,iyr1

; Test for end condition and close windows.
z=stddat(imn,idy,iyr,ndays)
if ndays lt lstday then stop,' starting day outside range '
if ndays gt ledday then stop,' Normal termination condition '

date=strcompress(string(FORMAT='(A3,A1,I2,A2,I4)',$
     month(imn-1),' ',idy,', ',iyr))
if idy ge 10 then datelab=strcompress(string(FORMAT='(A4,I2,A1,I4)',$
     mon(imn-1),idy,'_',iyr))
if idy lt 10 then datelab=strcompress(string(FORMAT='(A4,A1,I1,A1,I4)',$
     mon(imn-1),'0',idy,'_',iyr))
;
; read input
;
nc=0l
nr=0l 
nth=0l 
close,10
openr,10,datelab+'_rhj_div.dat'
readu,10,nr,nc,nth
xlat=fltarr(nr)
xlon=fltarr(nc)
th=fltarr(nth)
readu,10,xlat,xlon,th
div2=fltarr(nr,nc,nth)
rhj2=fltarr(nr,nc,nth)
q2=fltarr(nr,nc,nth)
p2=fltarr(nr,nc,nth)
readu,10,div2,rhj2,q2,p2
;goto,jumpsmooth1
;
; smooth in vertical
;
;for j=0,nr-1 do begin
;for i=0,nc-1 do begin
;dummy1=reform(div2(j,i,*),nth)
;dummy2=reform(rhj2(j,i,*),nth)
;for n=0,4 do begin
;dummy1=smooth(dummy1,3)
;dummy2=smooth(dummy2,3)
;endfor
;div2(j,i,*)=dummy1
;rhj2(j,i,*)=dummy2
;endfor
;endfor
;
; smooth in latitude
;
for l=0,nth-1 do begin
for i=0,nc-1 do begin
dummy1=reform(div2(*,i,l),nr)
;dummy2=reform(rhj2(*,i,l),nr)
for n=0,4 do begin
dummy1=smooth(dummy1,6)
;dummy2=smooth(dummy2,3)
endfor
div2(*,i,l)=dummy1
;rhj2(*,i,l)=dummy2
endfor
endfor
;
; smooth in longitude
;
for l=0,nth-1 do begin
for j=0,nr-1 do begin
dummy1=reform(div2(j,*,l),nc)
;dummy2=reform(rhj2(j,*,l),nc)
for n=0,4 do begin
dummy1=smooth(dummy1,6)
;dummy2=smooth(dummy2,3)
endfor
div2(j,*,l)=dummy1
;rhj2(j,*,l)=dummy2
endfor
endfor
jumpsmooth1:
if ndays eq lstday then begin
rhjnm1=rhj2
divnm1=div2
qnm1=q2
pnm1=p2
datelabnm1=datelab
goto, jump
endif
if ndays eq lstday+1 then begin
divn=div2
rhjn=rhj2
qn=q2
pn=p2
datelabn=datelab
goto, jump
endif
if ndays ge lstday+2 then begin
rhjnp1=rhj2
divnp1=div2
qnp1=q2
pnp1=p2
datelabnp1=datelab
endif
print,'n-1=',datelabnm1,' n=',datelabn,' n+1=',datelabnp1
thd2=rhjn*qn/(24.*60.*60.)
dt=2.*24.*60.*60.
drhjdt=(rhjnp1-rhjnm1)/dt
for i=0,nc-1 do begin
for j=0,nr-1 do begin
pres=reform(pn(j,i,*),nth)
index=where(pres ge 70.)
bth=index(0)-1
;
; loop over theta from top down
FOR thlev=bth,nth-1 DO BEGIN
dth=th(thlev)-th(thlev-1)
divbar=.5*(divn(j,i,thlev)+divn(j,i,thlev-1))
drhjdtbar=.5*(drhjdt(j,i,thlev)+drhjdt(j,i,thlev-1)) 
thd2(j,i,thlev)=thd2(j,i,thlev-1)-dth*(divn(j,i,thlev)+drhjdt(j,i,thlev))
ENDFOR  ; loop over theta
endfor
endfor
;
; loop over theta from top down
FOR thlev=0,nth-1 DO BEGIN
dummy1=reform(thd2(*,*,thlev),nr,nc)
dummy2=reform(rhjn(*,*,thlev),nr,nc)
index=where(dummy2 gt 0.)
if index(0) ne -1 then dummy1(index)=dummy1(index)/dummy2(index)
index=where(dummy2 eq 0.)
if index(0) ne -1 then dummy1(index)=0. 
thd2(*,*,thlev)=24.*60.*60.*dummy1 
ENDFOR  ; loop over theta
goto, jumpsmooth2
;
; smooth in vertical
;
for j=0,nr-1 do begin
for i=0,nc-1 do begin
dummy1=reform(thd2(j,i,*),nth)
for n=0,4 do begin
dummy1=smooth(dummy1,3)
endfor
thd2(j,i,*)=dummy1
endfor
endfor
;
; smooth in latitude
;
for l=0,nth-1 do begin
for i=0,nc-1 do begin
dummy1=reform(thd2(*,i,l),nr)
for n=0,4 do begin
dummy1=smooth(dummy1,3)
endfor
thd2(*,i,l)=dummy1
endfor
endfor
;
; smooth in longitude
;
for l=0,nth-1 do begin
for j=0,nr-1 do begin
dummy1=reform(thd2(j,*,l),nc)
for n=0,4 do begin
dummy1=smooth(dummy1,3)
endfor
thd2(j,*,l)=dummy1
endfor
endfor
jumpsmooth2:


thdbar=fltarr(nr,nth)
qbar=fltarr(nr,nth)
pbar=fltarr(nr,nth)
ybar=fltarr(nr,nth)
for j=0,nr-1 do begin
for l=0,nth-1 do begin
thdbar(j,l)=total(thd2(j,*,l))/nc
qbar(j,l)=total(qn(j,*,l))/nc
pbar(j,l)=total(pn(j,*,l))/nc
ybar(j,l)=xlat(j)
endfor
endfor
level=-3.0+.2*findgen(31)
!p.multi=[0,3,1]
plot_io,ybar,pbar,xrange=[-90,90],yrange=[1000.,10.],/nodata,$
xtitle='Latitude',ytitle='Pressure'
contour,qbar,ybar,pbar,level=level,$
c_linestyle = level lt 0.,/overplot
plot_io,ybar,pbar,xrange=[-90,90],yrange=[1000.,10.],/nodata,$
xtitle='Latitude',ytitle='Pressure'
contour,thdbar,ybar,pbar,level=level,$ 
c_linestyle = level lt 0.,/overplot
plot_io,ybar,pbar,xrange=[-90,90],yrange=[1000.,10.],/nodata,$
xtitle='Latitude',ytitle='Pressure'
contour,qbar-thdbar,ybar,pbar,level=level,$ 
c_linestyle = level lt 0.,/overplot
;
; write output
;
close,10
openw,10,datelabn+'_theta_dot.dat'
writeu,10,nr,nc,nth
writeu,10,xlat,xlon,th
writeu,10,thd2,qn,pn
;
rhjnm1=rhjn
divnm1=divn
qnm1=qn
datelabnm1=datelabn
;
rhjn=rhjnp1
divn=divnp1
qn=qnp1
datelabn=datelabnp1
goto, jump	; loop over days
end
