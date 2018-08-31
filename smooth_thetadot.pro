

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
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1991 or lstyr gt 1999 then stop,'Year out of range '
if ledyr lt 1991 or ledyr gt 1999 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '

if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
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
if ndays gt ledday then goto, jumpout 

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
openr,10,datelab+'_theta_dot.dat'
readu,10,nr,nc,nth
xlat=fltarr(nr)
xlon=fltarr(nc)
th=fltarr(nth)
readu,10,xlat,xlon,th
thd2=fltarr(nr,nc,nth)
q2=fltarr(nr,nc,nth)
p2=fltarr(nr,nc,nth)
readu,10,thd2,q2,p2
if ndays eq lstday then begin
icount=0
thdbar=fltarr(nr,nc,nth)
qbar=fltarr(nr,nc,nth)
pbar=fltarr(nr,nc,nth)
endif
thdbar=thdbar+thd2
qbar=qbar+q2
pbar=pbar+p2
icount=icount+1
goto, jump	; loop over days
jumpout:
thdbar=thdbar/icount
qbar=qbar/icount
pbar=pbar/icount
;
; zonal average
;
thdz=fltarr(nr,nth)
qz=fltarr(nr,nth)
pz=fltarr(nr,nth)
yz=fltarr(nr,nth)
for l=0,nth-1 do begin
yz(*,l)=xlat
for j=0,nr-1 do begin
thdz(j,l)=total(thdbar(j,*,l))/nc
qz(j,l)=total(qbar(j,*,l))/nc
pz(j,l)=total(pbar(j,*,l))/nc
endfor
endfor
!p.multi=[0,2,1]
;
level=-20+2*findgen(21)
plot_io,yz,pz,xrange=[-90,90],yrange=[1000.,1.],$
xtitle='Latitude',ytitle='Pressure',title='Shine Heating (K/day)',/nodata
contour,qz,yz,pz,/overplot,level=level,c_linestyle=level lt 0.,/follow
level=-2+.2*findgen(21)
index=where(pz le 70.)
qz(index)=9999.
contour,qz,yz,pz,max_val=9999.,/overplot
;
level=-20+2*findgen(21)
plot_io,yz,pz,xrange=[-90,90],yrange=[1000.,1.],$
xtitle='Latitude',ytitle='Pressure',title='Div Heating (K/day)',/nodata
contour,thdz,yz,pz,/overplot,level=level,c_linestyle=level lt 0.,/follow
level=-2+.2*findgen(21)
index=where(pz le 70.)
thdz(index)=9999.
contour,thdz,yz,pz,max_val=9999.,/overplot
if setplot eq 'x' then stop

xlon0=xlon
xlon(0:nc/2-1)=xlon(nc/2:nc-1)
xlon(nc/2:nc-1)=xlon0(0:nc/2-1)
index=where(xlon ge 180.)
xlon(index)=xlon(index)-360.
;
; layer by layer
;
for l=0,nth-1 do begin
!p.multi=[0,2,1]
map_set,0,0,0,/continents,title='Shine Heating '+string(th(l))+'K'
q0=transpose(reform(qbar(*,*,l),nr,nc))
q=q0
q(0:nc/2-1,*)=q(nc/2:nc-1,*)
q(nc/2:nc-1,*)=q0(0:nc/2-1,*)
if max(abs(q)) ge 10 then level=-20+2*findgen(21)
if max(abs(q)) ge 5 and max(abs(q)) lt 10 then level=-10+1*findgen(21)
if max(abs(q)) lt 5 then level=-5+.5*findgen(21)
if max(abs(q)) lt 2 then level=-2+.2*findgen(21)
!p.multi=[0,2,1]
contour,q,xlon,xlat,level=level,/overplot,c_linestyle=level lt 0.,/follow
!p.multi=[1,2,1]
map_set,0,0,0,/continents,/noeras,title='Div Heating '+string(th(l))+'K'
q0=transpose(reform(thdbar(*,*,l),nr,nc))
q=q0
q(0:nc/2-1,*)=q(nc/2:nc-1,*)
q(nc/2:nc-1,*)=q0(0:nc/2-1,*)
!p.multi=[1,2,1]
contour,q,xlon,xlat,level=level,/overplot,c_linestyle=level lt 0.,/follow
if setplot eq 'x' then stop
endfor
end
