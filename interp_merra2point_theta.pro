;
; input x,y,date,time
; ouput MERRA theta data at Poker Flat 65.1200° N, 147.4700W (212.53E)
;
@calcelat2d
@stddat
@kgmt
@ckday
@kdate
@rd_merra_nc3

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
diru='/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_theta_'
pthout='/Volumes/Data/MERRA_data/Datfiles/'
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
lstmn=1L & lstdy=1L & lstyr=1998L
ledmn=12L & leddy=31L & ledyr=1998L
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
time=12.
;read,' Enter UT time (0-24 hours) ',time
stime=strcompress(string(time))
time=time/24.   ; convert to fractional day
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
nr=91L
yeq=findgen(nr)
latcircle=fltarr(nr)
hem_frac=fltarr(nr)
for j=0,nr-1 do begin
    hy=re*dtr
    dx=re*cos(yeq(j)*dtr)*360.*dtr
    latcircle(j)=dx*hy
endfor
;
; need latcircle fully initialized
;
for j=0,nr-1 do begin
    index=where(yeq ge yeq(j))
    if index(0) ne -1 then hem_frac(j)=total(latcircle(index))/hem_area
    if yeq(j) eq 0. then hem_frac(j)=1.
endfor

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,plotit
;
; calculate day of year
;
z = kgmt(imn,idy,iyr,jday)
;
; determine 2 bounding dates in (month, day, year) format
; based on fractional year and day and analyses valid at 12Z
;
if time lt 0.5 then begin
    jday0=jday-1.0
    jday1=jday
    tscale=time+0.5
endif
if time ge 0.5 then begin
    jday0=jday
    jday1=jday+1.0
    tscale=time-0.5
endif
iyr0=iyr
iyr1=iyr
kdate,float(jday0),iyr0,imn0,idy0
ckday,jday0,iyr0
kdate,float(jday1),iyr1,imn1,idy1
ckday,jday1,iyr1
if iyr0 lt 2000 then iyr0=iyr0-1900
if iyr0 ge 2000 then iyr0=iyr0-2000
if iyr1 lt 2000 then iyr1=iyr1-1900
if iyr1 ge 2000 then iyr1=iyr1-2000
;
; read MERRA data on day 0
;
newday=0L
sdate=string(FORMAT='(i4,i2.2,i2.2)',iyr,imn,idy)
ifile=diru+sdate+'.nc3'
print,ifile
dum=findfile(ifile)
if dum(0) ne '' then begin
;   rd_merra_nc3,ifile,ncw,nrw,nthw,alon,alat,th,pvold,pold,$
;               uold,vold,qdfold,markold,qvold,gold,sfold,qold,iflag
    rd_merra_nc3,ifile,ncw,nrw,nthw,alon,alat,th,pvgrd,pgrd,$
                ugrd,vgrd,qdfgrd,markgrd,qvgrd,ggrd,sfgrd,qgrd,iflag
   index=where(markgrd lt -1.)
   if index(0) ne -1L then markgrd(index)=-1.0*markgrd(index)/markgrd(index)	; anticyclones= -1
   markgrd=smooth(markgrd,3,/edge_truncate)
;
; on the first day declare profile arrays
;
   if kcount eq 0L then begin
      dates=lonarr(kday)
      zonal_wind_profiles=fltarr(kday,nthw)
      meridional_wind_profiles=fltarr(kday,nthw)
      vortex_marker_profiles=fltarr(kday,nthw)
      pressure_profiles=fltarr(kday,nthw)
      temperature_profiles=fltarr(kday,nthw)
      altitude_profiles=fltarr(kday,nthw)
      kcount=1
   endif
endif
dates(icount)=long(sdate)
;;
;; read MERRA data on day 1
;;
;ifile=diru+sdate+'.nc3'
;print,'file 2 ',ifile
;dum=findfile(ifile)
;if dum(0) ne '' then begin
;   rd_merra_nc3,ifile,ncw,nrw,nthw,alon,alat,th,pvnew,pnew,$
;               unew,vnew,qdfnew,marknew,qvnew,gnew,sfnew,qnew,iflag
;   index=where(marknew lt -1.)
;   if index(0) ne -1L then marknew(index)=-1.0*marknew(index)/marknew(index)    ; anticyclones= -1
;   marknew=smooth(marknew,3,/edge_truncate)
;endif
;
; check that point is within latitude range
;
if slat lt min(alat) then stop,'slat lt min lat'   ; slat=min(alat)
if slat gt max(alat) then stop,'slat gt max lat'   ; slat=max(alat)
;;
;; perform time interpolation
;;
;ugrd=uold+TSCALE*(unew-uold)
;vgrd=vold+TSCALE*(vnew-vold)
;pgrd=pold+TSCALE*(pnew-pold)
;ggrd=gold+TSCALE*(gnew-gold)
;markgrd=markold+TSCALE*(marknew-markold)
;
; calculate geopotential height of isentropic surface = (msf - cp*T)/g
; where T = theta* (p/po)^R/cp and divide by 1000 for km
;
tgrd=0.*pgrd
zgrd=0.*pgrd
for k=0L,nthw-1L do begin
    tgrd(0:nrw-1,0:ncw-1,k)=th(k)*( (pgrd(0:nrw-1,0:ncw-1,k)/1000.)^(.286) )
;
; convert geopotential height to geometric height
;
    for j=0L,nrw-1L do begin
        sin2=sin( (alat(j)*dtr)^2.0 )
        numerator=1.0+ks*sin2
        denominator=sqrt( 1.0 - (ecc^2.0)*sin2 )
        gammas=gamma45*(numerator/denominator)
        r=6378.137/(1.006803-(0.006706*sin2))
        zgrd(j,*,k)=(r*ggrd(j,*,k))/ ( (gammas/gamma45)*r - ggrd(j,*,k) )
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
pj1=ugrd(j,i,*)+xscale*(ugrd(j,ip1,*)-ugrd(j,i,*))
pjp1=ugrd(jp1,i,*)+xscale*(ugrd(jp1,ip1,*)-ugrd(jp1,i,*))
zonal_wind_profiles(icount,*)=reform(pj1+yscale*(pjp1-pj1))

pj1=vgrd(j,i,*)+xscale*(vgrd(j,ip1,*)-vgrd(j,i,*))
pjp1=vgrd(jp1,i,*)+xscale*(vgrd(jp1,ip1,*)-vgrd(jp1,i,*))
meridional_wind_profiles(icount,*)=reform(pj1+yscale*(pjp1-pj1))

pj1=markgrd(j,i,*)+xscale*(markgrd(j,ip1,*)-markgrd(j,i,*))
pjp1=markgrd(jp1,i,*)+xscale*(markgrd(jp1,ip1,*)-markgrd(jp1,i,*))
vortex_marker_profiles(icount,*)=reform(pj1+yscale*(pjp1-pj1))

pj1=pgrd(j,i,*)+xscale*(pgrd(j,ip1,*)-pgrd(j,i,*))
pjp1=pgrd(jp1,i,*)+xscale*(pgrd(jp1,ip1,*)-pgrd(jp1,i,*))
pressure_profiles(icount,*)=reform(pj1+yscale*(pjp1-pj1))

pj1=tgrd(j,i,*)+xscale*(tgrd(j,ip1,*)-tgrd(j,i,*))
pjp1=tgrd(jp1,i,*)+xscale*(tgrd(jp1,ip1,*)-tgrd(jp1,i,*))
temperature_profiles(icount,*)=reform(pj1+yscale*(pjp1-pj1))

pj1=zgrd(j,i,*)+xscale*(zgrd(j,ip1,*)-zgrd(j,i,*))
pjp1=zgrd(jp1,i,*)+xscale*(zgrd(jp1,ip1,*)-zgrd(jp1,i,*))
altitude_profiles(icount,*)=reform(pj1+yscale*(pjp1-pj1))

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
          /bold,/color,bits_per_pixel=8,/helvetica,filename=pthout+'merra2poker_'+sdate+'.ps'
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
tmax=max(temperature_profiles)+5.
tmin=min(temperature_profiles)-5.
nlvls=20
level=tmin+((tmax-tmin)/float(nlvls))*findgen(nlvls+1)
nlvls=nlvls+1
col1=1+indgen(nlvls)*icolmax/float(nlvls)
contour,temperature_profiles,findgen(kday),th,levels=level,c_color=col1,/cell_fill,/noeras,title='T at Poker Flat '+sloc,xrange=[0.,kday-1],color=0,yrange=[min(th),max(th)]
;contour,temperature_profiles,findgen(kday),th,/overplot,levels=level,/follow,c_labels=0*level,/noeras,color=0
;contour,vortex_marker_profiles,findgen(kday),th,/overplot,levels=[0.1],thick=5,color=0
;contour,vortex_marker_profiles,findgen(kday),th,/overplot,levels=[-0.1],thick=5,color=mcolor
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
tmax=max(zonal_wind_profiles)+5.
tmin=min(zonal_wind_profiles)-5.
nlvls=20
level=tmin+((tmax-tmin)/float(nlvls))*findgen(nlvls+1)
nlvls=nlvls+1
col1=1+indgen(nlvls)*icolmax/float(nlvls)
contour,zonal_wind_profiles,findgen(kday),th,levels=level,c_color=col1,/cell_fill,/noeras,title='U at Poker Flat '+sloc,xrange=[0.,kday-1],color=0,yrange=[min(th),max(th)]
;contour,zonal_wind_profiles,findgen(kday),th,/overplot,levels=level,/follow,c_labels=0*level,/noeras,color=0
;contour,vortex_marker_profiles,findgen(kday),th,/overplot,levels=[0.1],thick=5,color=0
;contour,vortex_marker_profiles,findgen(kday),th,/overplot,levels=[-0.1],thick=5,color=mcolor
xyouts,-.98,450.,'Anticyclone',/data,color=0
xyouts,0.25,450.,'Vortex',/data,color=0
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
potential_temperature_profile=th
date=long(sdate)
comment=strarr(11)
comment(0)='Profiles based on MERRA daily averaged data'
comment(1)='dates in YYYYMMDD'
comment(2)='longitude in degrees east'
comment(3)='latitude in degrees'
comment(4)='potential_temperature = potential temperature profile (K)'
comment(5)='vortex_marker_profiles = positive (negative) values in vortex (anticyclones)'
comment(6)='pressure_profiles = Pressure profiles (hPa)'
comment(7)='temperature_profiles = Temperature profiles (K)'
comment(8)='altitude_profiles = Geometric Altitude profiles (km)'
comment(9)='zonal_wind_profiles = Zonal Wind profiles (km)'
comment(10)='meridional_wind_profiles = Meridional Wind profiles (km)'

save,file=pthout+'merra_at_poker_'+sdate+'.sav',dates,time,longitude,latitude,$
     potential_temperature_profile,vortex_marker_profiles,pressure_profiles,$
     temperature_profiles,altitude_profiles,zonal_wind_profiles,meridional_wind_profiles,comment
print,'saved DMP '+sdate

if setplot eq 'x' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim '+pthout+'merra2poker_'+sdate+'.ps -rotate -90 '+pthout+'merra2poker_'+sdate+'.jpg'
   spawn,'/usr/bin/rm '+pthout+'merra2poker_'+sdate+'.ps'
endif
end
