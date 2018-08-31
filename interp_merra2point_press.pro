;
; input x,y,date
; ouput MERRA pressure data at Poker Flat 65.1200° N, 147.4700W (212.53E)
;
@stddat
@kgmt
@ckday
@kdate

a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
loadct,39
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
device,decompose=0
setplot='ps'
read,'setplot=',setplot
nxdim=700
nydim=700
xorig=[0.1,0.1]
yorig=[0.6,0.15]
xlen=0.8
ylen=0.3
cbaryoff=0.03
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=mcolor
endif
diru='/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_press_'
pthout='/Volumes/Data/MERRA_data/Datfiles/'
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
lstmn=1L & lstdy=1L & lstyr=1998L
ledmn=9L & leddy=30L & ledyr=2014L
lstday=0L & ledday=0L
;
; Ask interactive questions- get starting/ending date and p surface
;
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
;read,' Enter starting year ',lstyr
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1979 then stop,'Year out of range '
if ledyr lt 1979 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L
kcount=0L
slat=-67.57 & slon=291.88		; Rothera
slat=65.12 & slon=212.53			; Poker Flat
;read,' Enter longitude, latitude, theta ',slon,slat
sloc='('+string(FORMAT='(f6.2)',slon)+','+string(FORMAT='(f6.2)',slat)+')'
re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd
ks=1.931853d-3
ecc=0.081819
gamma45=9.80

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,plotit
;
; read MERRA data 
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
sdate=string(FORMAT='(i4,i2.2,i2.2)',iyr,imn,idy)
ifile=diru+sdate+'.sav'
print,ifile
dum=findfile(ifile)
if dum(0) ne '' then begin
   restore,ifile
   alon=LONGITUDE_WACCM
   alat=LATITUDE_WACCM
   ncw=n_elements(alon)
   nrw=n_elements(alat)
   nthw=n_elements(pressure)
   ggrd=zgrd
;
; on the first day declare profile arrays
;
   if kcount eq 0L then begin
      dates=lonarr(kday)
      u_profiles=fltarr(kday,nthw)
      v_profiles=fltarr(kday,nthw)
      t_profiles=fltarr(kday,nthw)
      z_profiles=fltarr(kday,nthw)
      kcount=1
   endif
endif
dates(icount)=long(sdate)
;
; check that point is within latitude range
;
if slat lt min(alat) then stop,'slat lt min lat'   ; slat=min(alat)
if slat gt max(alat) then stop,'slat gt max lat'   ; slat=max(alat)
;
; convert geopotential height to geometric height
;
zgrd=0.*tgrd
for k=0L,nthw-1L do begin
    for j=0L,nrw-1L do begin
        sin2=sin( (alat(j)*dtr)^2.0 )
        numerator=1.0+ks*sin2
        denominator=sqrt( 1.0 - (ecc^2.0)*sin2 )
        gammas=gamma45*(numerator/denominator)
        r=6378.137/(1.006803-(0.006706*sin2))
        zgrd(*,j,k)=(r*ggrd(*,j,k))/ ( (gammas/gamma45)*r - ggrd(*,j,k) )
    endfor
endfor
;
; interpolate MERRA to point location
;
if slon lt alon(0) then slon=slon+360.
for i=0L,ncw-1L do begin
    ip1=i+1
    if i eq ncw-1L then ip1=0L
    xlon=alon(i)
    xlonp1=alon(ip1)
    if i eq ncw-1L then xlonp1=360.+alon(ip1)
    if slon ge xlon and slon le xlonp1 then begin
       xscale=(slon-xlon)/(xlonp1-xlon)
       goto,jumpx
    endif
endfor
jumpx:
for j=0L,nrw-2L do begin
    jp1=j+1
    xlat=alat(j)
    xlatp1=alat(jp1)
    if slat ge xlat and slat le xlatp1 then begin
        yscale=(slat-xlat)/(xlatp1-xlat)
        goto,jumpy
    endif
endfor
jumpy:
;
; interpolate gridded data to profile lon/lat
;
pj1=ugrd(i,j,*)+xscale*(ugrd(ip1,j,*)-ugrd(i,j,*))
pjp1=ugrd(i,jp1,*)+xscale*(ugrd(ip1,jp1,*)-ugrd(i,jp1,*))
u_profiles(icount,*)=reform(pj1+yscale*(pjp1-pj1))

pj1=vgrd(i,j,*)+xscale*(vgrd(ip1,j,*)-vgrd(i,j,*))
pjp1=vgrd(i,jp1,*)+xscale*(vgrd(ip1,jp1,*)-vgrd(i,jp1,*))
v_profiles(icount,*)=reform(pj1+yscale*(pjp1-pj1))

pj1=tgrd(i,j,*)+xscale*(tgrd(ip1,j,*)-tgrd(i,j,*))
pjp1=tgrd(i,jp1,*)+xscale*(tgrd(ip1,jp1,*)-tgrd(i,jp1,*))
t_profiles(icount,*)=reform(pj1+yscale*(pjp1-pj1))

pj1=zgrd(i,j,*)+xscale*(zgrd(ip1,j,*)-zgrd(i,j,*))
pjp1=zgrd(i,jp1,*)+xscale*(zgrd(ip1,jp1,*)-zgrd(i,jp1,*))
z_profiles(icount,*)=reform(pj1+yscale*(pjp1-pj1))

icount=icount+1L
goto,jump

plotit:
;
; plot for posterity
;
if setplot eq 'ps' then begin
   xsize=nxdim/100.
   ysize=nydim/100.
   set_plot,'ps'
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename=pthout+'merra2poker_'+sdate+'_press.ps'
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=5
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
erase
!type=2^2+2^3
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
tmax=max(t_profiles)+5.
tmin=min(t_profiles)-5.
nlvls=20
level=tmin+((tmax-tmin)/float(nlvls))*findgen(nlvls+1)
nlvls=nlvls+1
col1=1+indgen(nlvls)*icolmax/float(nlvls)
contour,t_profiles,findgen(kday),pressure,levels=level,c_color=col1,/cell_fill,/noeras,title='T at Poker Flat '+sloc,xrange=[0.,kday-1],color=0,yrange=[max(pressure),min(pressure)],/ylog
;contour,t_profiles,findgen(kday),pressure,/overplot,levels=level,/follow,c_labels=0*level,/noeras,color=0
imin=min(level)
imax=max(level)
ymnb=ymn -cbaryoff
ymxb=ymnb+cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,xtitle='Temperature (K)'
ybox=[0,10,10,0,0]
x2=imin
dx=(imax-imin)/(float(nlvls)-1)
for j=1,nlvls-1 do begin
    xbox=[x2,x2,x2+dx,x2+dx,x2]
    polyfill,xbox,ybox,color=col1(j)
    x2=x2+dx
endfor

!type=2^2+2^3
xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
tmax=max(u_profiles)+5.
tmin=min(u_profiles)-5.
nlvls=20
level=tmin+((tmax-tmin)/float(nlvls))*findgen(nlvls+1)
nlvls=nlvls+1
col1=1+indgen(nlvls)*icolmax/float(nlvls)
contour,u_profiles,findgen(kday),pressure,levels=level,c_color=col1,/cell_fill,/noeras,title='U at Poker Flat '+sloc,xrange=[0.,kday-1],color=0,yrange=[max(pressure),min(pressure)],/ylog
;contour,u_profiles,findgen(kday),pressure,/overplot,levels=level,/follow,c_labels=0*level,/noeras,color=0
imin=min(level)
imax=max(level)
ymnb=ymn -cbaryoff
ymxb=ymnb+cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,xtitle='Zonal Wind (m/s)'
ybox=[0,10,10,0,0]
x2=imin
dx=(imax-imin)/(float(nlvls)-1)
for j=1,nlvls-1 do begin
    xbox=[x2,x2,x2+dx,x2+dx,x2]
    polyfill,xbox,ybox,color=col1(j)
    x2=x2+dx
endfor
;
; save all profiles
;
longitude=slon
latitude=slat
date=long(sdate)
comment=strarr(9)
comment(0)='Profiles of U, V, T, Z based on MERRA daily averaged data'
comment(1)='dates in YYYYMMDD'
comment(2)='longitude in degrees east'
comment(3)='latitude in degrees'
comment(4)='pressure = MERRA pressure levels (hPa)'
comment(5)='t_profiles = Temperature profiles (K)'
comment(6)='z_profiles = Geometric Altitude profiles (km)'
comment(7)='u_profiles = Zonal Wind profiles (m/s)'
comment(8)='v_profiles = Meridional Wind profiles (m/s)'

save,file=pthout+'merra_at_poker_press.sav',dates,longitude,latitude,$
     pressure,t_profiles,z_profiles,u_profiles,v_profiles,comment
print,'saved DMP '+sdate

if setplot eq 'x' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim '+pthout+'merra2poker_'+sdate+'_press.ps -rotate -90 '+pthout+'merra2poker_'+sdate+'_press.jpg'
endif
end
