;
; zonal mean wind from MERRA pressure data
;
@stddat
@kgmt
@ckday
@kdate

loadct,39
device,decompose=0
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
icmm1=icolmax-1
icmm2=icolmax-2
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
!noeras=1
nxdim=750
nydim=750
xorig=[0.15,0.15]
yorig=[0.60,0.15]
xlen=0.7
ylen=0.3
cbaryoff=0.08
cbarydel=0.01
setplot='x'
read,'setplot=',setplot
if setplot ne 'ps' then begin
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=icolmax
endif
dir='/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_press_'
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
lstmn=1
lstdy=1
lstyr=2008
ledmn=12
leddy=31
ledyr=2008
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
;print, ' '
;print, '      MERRA Version '
;print, ' '
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1991 then stop,'Year out of range '
if ledyr lt 1991 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
;
; Compute initial Julian date
;
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then stop,' normal termination condition '
      if iyr ge 2000L then iyr1=iyr-2000L
      if iyr lt 2000L then iyr1=iyr-1900L
;
;***Read data
;
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
;
; LATITUDE_WACCM  FLOAT     = Array[96]
; LONGITUDE_WACCM FLOAT     = Array[144]
; PRESSURE        FLOAT     = Array[41]
; PSGRD           FLOAT     = Array[144, 96]
; QVGRD           FLOAT     = Array[144, 96, 41]
; TGRD            FLOAT     = Array[144, 96, 41]
; UGRD            FLOAT     = Array[144, 96, 41]
; VGRD            FLOAT     = Array[144, 96, 41]
; ZGRD            FLOAT     = Array[144, 96, 41]
;
      restore,dir+sdate+'.sav'
      alat=LATITUDE_WACCM
;
; zonal mean temperature and zonal wind
;
      tzm=mean(tgrd,dim=1)
      uzm=mean(ugrd,dim=1)
;
; read mark
;
   ncid=ncdf_open('/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_theta_'+sdate+'.nc3')
   ncdf_diminq,ncid,0,name,nlat
   ncdf_diminq,ncid,1,name,nlg
   ncdf_diminq,ncid,2,name,nth
   mlon=fltarr(nlg)
   mlat=fltarr(nlat)
   thlev=fltarr(nth)
   p2=fltarr(nlat,nlg,nth)
   mark2=fltarr(nlat,nlg,nth)
   ncdf_varget,ncid,0,mlon
   ncdf_varget,ncid,1,mlat
   ncdf_varget,ncid,2,thlev
   ncdf_varget,ncid,4,p2
   ncdf_varget,ncid,8,mark2
   ncdf_close,ncid

mzm=mean(mark2,dim=2)
pzm=mean(p2,dim=2)

if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='yz_merra_ubarp_'+sdate+'.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
   !p.thick=2.0                   ;Plotted lines twice as thick
   !p.charsize=1.0
endif
;
; plot zonal mean zonal wind
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7
nlvls=26
col1=1+indgen(nlvls)*icolmax/nlvls
level=170.+5.*findgen(nlvls)
contour,tzm,alat,pressure,/ylog,/noeras,xrange=[-90.,90.],yrange=[10.,min(pressure)],charsize=1.5,color=0,$
      ytitle='Pressure (hPa)',title='MERRA Zonal Mean Temperature '+sdate,xticks=6,/fill,c_color=col1,$
      levels=level
index=where(level mod 10. eq 0)
contour,tzm,alat,pressure,levels=level(index),color=0,/follow,/overplot,c_labels=1+0*index
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(K)'
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor

xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7
nlvls=25
col1=1+indgen(nlvls)*icolmax/nlvls
level=-120.+10.*findgen(nlvls)
contour,uzm,alat,pressure,/ylog,/noeras,xrange=[-90.,90.],yrange=[10.,min(pressure)],charsize=1.5,color=0,$
      ytitle='Pressure (hPa)`',title='Zonal Mean Zonal Wind '+sdate,xticks=6,/fill,c_color=col1,$
      levels=level
index=where(level gt 0.)
contour,uzm,alat,pressure,levels=level(index),color=mcolor,/follow,/overplot
index=where(level lt 0.)
contour,uzm,alat,pressure,levels=level(index),color=0,/follow,/overplot
contour,uzm,alat,pressure,levels=[0],color=0,/follow,/overplot,thick=3

contour,mzm,mlat,pzm,levels=0.1+0.1*findgen(9),/follow,thick=5,color=0,/overplot

imin=min(level)
imax=max(level)
ymnb=yorig(1) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(m/s)'
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor

if setplot ne 'ps' then wait,1
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim yz_merra_ubarp_'+sdate+'.ps -rotate -90 yz_merra_ubarp_'+sdate+'.jpg'
endif
goto,jump
end
