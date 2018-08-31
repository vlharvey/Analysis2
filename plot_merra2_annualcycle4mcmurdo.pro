;
; read timeseries of MERRA2 at McMurdo (and vortex distance and vortex edge speed) and plot 5-year average annual cycles
;
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
xorig=[0.10,0.10,0.10,0.6,0.6,0.6]
yorig=[0.75,0.425,0.1,0.75,0.425,0.1]
xlen=0.35
ylen=0.2
cbaryoff=0.05
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=mcolor
endif
diru='/atmos/harvey/MERRA2_data/Datfiles/MERRA2-on-WACCM_theta_'
pthout='/atmos/harvey/MERRA2_data/Pre_process/McMurdo/'
mon=['J','F','M','A','M','J','J','A','S','O','N','D']
nmonth=n_elements(mon)

;comment(0)='Profiles based on MERRA-2 0Z data'
;comment(1)='dates in YYYYMMDD'
;comment(2)='longitude in degrees east'
;comment(3)='latitude in degrees'
;comment(4)='altitude = Geometric altitude profile (km)'
;comment(5)='vortex_marker_profiles = positive values in vortex'
;comment(6)='pressure_profiles = Pressure profiles (hPa)'
;comment(7)='temperature_profiles = Temperature profiles (K)'
;comment(8)='zonal_wind_profiles = Zonal Wind profiles (km)'
;comment(9)='meridional_wind_profiles = Meridional Wind profiles (km)'
;comment(10)='vortex_distance_profiles = km distance McMurdo Station is from the Antarctic vortex edge (negative if station is inside vortex, positive if station is outside vortex, /Nan if no vortex is present)'
;comment(11)='vortex_edge_wind_profiles = mean wind speed at the edge of the polar vortex (not at McMurdo)'
;
; kday,nz,dates,time,longitude,latitude,altitude,$
; vortex_marker_profiles_z,pressure_profiles_z,temperature_profiles_z,$
; zonal_wind_profiles_z,meridional_wind_profiles_z,vortex_distance_profiles_z,vortex_edge_wind_profiles_z,comment

restore,pthout+'merra2_at_mcmurdo_20110101-20151231.sav'	

plotit:
datelab=strcompress(dates(0),/r)+'-'+strcompress(dates(-1),/r)
slat=-77.85 & slon=166.67                       ; McMurdo is 77°51'S, 166°40'E
sloc='('+string(FORMAT='(f6.2)',slon)+','+string(FORMAT='(f6.2)',slat)+')'
;
; monthly average of variables between 30-50 km
;
u_avg=fltarr(nmonth,nz)
v_avg=fltarr(nmonth,nz)
t_avg=fltarr(nmonth,nz)
dist_avg=fltarr(nmonth,nz)
mark_avg=fltarr(nmonth,nz)
edgespeed_avg=fltarr(nmonth,nz)
sdates=strcompress(dates,/r)
smon=strmid(sdates,4,2)
sday=strmid(sdates,6,2)
syear=strmid(sdates,0,4)
;
; loop over months
;
for imonth=0L,nmonth-1L do begin
    thismonth=where(imonth+1 eq long(smon),npts)
    u_avg(imonth,*)=mean(ZONAL_WIND_PROFILES_Z(thismonth,*),dim=1,/Nan)
    v_avg(imonth,*)=mean(MERIDIONAL_WIND_PROFILES_Z(thismonth,*),dim=1,/Nan)
    t_avg(imonth,*)=mean(TEMPERATURE_PROFILES_Z(thismonth,*),dim=1,/Nan)
    dist_avg(imonth,*)=mean(VORTEX_DISTANCE_PROFILES_Z(thismonth,*),dim=1,/Nan)
    mark_avg(imonth,*)=mean(VORTEX_MARKER_PROFILES_Z(thismonth,*),dim=1,/Nan)
    edgespeed_avg(imonth,*)=mean(VORTEX_EDGE_WIND_PROFILES_Z(thismonth,*),dim=1,/Nan)
endfor
zindex=where(altitude eq 30)
uavg30=reform(u_avg(*,zindex))
vavg30=reform(v_avg(*,zindex))
tavg30=reform(t_avg(*,zindex))
distavg30=reform(dist_avg(*,zindex))
markavg30=reform(mark_avg(*,zindex))
edgespeedavg30=reform(edgespeed_avg(*,zindex))
index=where(markavg30 eq 0.)
if index(0) ne -1L then edgespeedavg30(index)=0./0.
if index(0) ne -1L then distavg30(index)=0./0.
if index(0) ne -1L then markavg30(index)=0./0.

zindex=where(altitude eq 40)
uavg40=reform(u_avg(*,zindex))
vavg40=reform(v_avg(*,zindex))
tavg40=reform(t_avg(*,zindex))
distavg40=reform(dist_avg(*,zindex))
markavg40=reform(mark_avg(*,zindex))
edgespeedavg40=reform(edgespeed_avg(*,zindex))
index=where(markavg40 eq 0.)
if index(0) ne -1L then edgespeedavg40(index)=0./0.
if index(0) ne -1L then distavg40(index)=0./0.
if index(0) ne -1L then markavg40(index)=0./0.

zindex=where(altitude eq 50)
uavg50=reform(u_avg(*,zindex))
vavg50=reform(v_avg(*,zindex))
tavg50=reform(t_avg(*,zindex))
distavg50=reform(dist_avg(*,zindex))
markavg50=reform(mark_avg(*,zindex))
edgespeedavg50=reform(edgespeed_avg(*,zindex))
index=where(markavg50 eq 0.)
if index(0) ne -1L then edgespeedavg50(index)=0./0.
if index(0) ne -1L then distavg50(index)=0./0.
if index(0) ne -1L then markavg50(index)=0./0.
;
; plot for posterity
;
if setplot eq 'ps' then begin
   xsize=nxdim/100.
   ysize=nydim/100.
   set_plot,'ps'
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename=pthout+'annualcycle_merra2mcmurdo.ps'
   !p.font=0
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
;
; xticks
;
erase
!type=2^2+2^3
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
plot,1+findgen(nmonth),uavg30,thick=8,title='U at McMurdo',color=0,xtickname=mon,xticks=nmonth-1,psym=-2,yrange=[-20,70],charsize=1.25,ytitle='m/s'
oplot,1+findgen(nmonth),0*uavg30,color=0
oplot,1+findgen(nmonth),uavg40,thick=8,color=100,psym=-2
oplot,1+findgen(nmonth),uavg50,thick=8,color=250,psym=-2

xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
plot,1+findgen(nmonth),vavg30,thick=8,title='V at McMurdo',color=0,xtickname=mon,xticks=nmonth-1,psym=-2,charsize=1.25,yrange=[-10,25],ytitle='m/s'
oplot,1+findgen(nmonth),0*vavg30,color=0
oplot,1+findgen(nmonth),vavg40,thick=8,color=100,psym=-2
oplot,1+findgen(nmonth),vavg50,thick=8,color=250,psym=-2

xmn=xorig(2)
xmx=xorig(2)+xlen
ymn=yorig(2)
ymx=yorig(2)+ylen
set_viewport,xmn,xmx,ymn,ymx
plot,1+findgen(nmonth),tavg30,thick=8,title='T at McMurdo',color=0,xtickname=mon,xticks=nmonth-1,psym=-2,charsize=1.25,yrange=[195,300],ytitle='K'
oplot,1+findgen(nmonth),tavg40,thick=8,color=100,psym=-2
oplot,1+findgen(nmonth),tavg50,thick=8,color=250,psym=-2

xmn=xorig(3)
xmx=xorig(3)+xlen
ymn=yorig(3)
ymx=yorig(3)+ylen
set_viewport,xmn,xmx,ymn,ymx
plot,1+findgen(nmonth),edgespeedavg30,thick=8,title='EDGE SPEED',color=0,xtickname=mon,xticks=nmonth-1,psym=-2,charsize=1.25,yrange=[0,140],ytitle='m/s'
oplot,1+findgen(nmonth),edgespeedavg40,thick=8,color=100,psym=-2
oplot,1+findgen(nmonth),edgespeedavg50,thick=8,color=250,psym=-2

xmn=xorig(4)
xmx=xorig(4)+xlen
ymn=yorig(4)
ymx=yorig(4)+ylen
set_viewport,xmn,xmx,ymn,ymx
plot,1+findgen(nmonth),markavg30,thick=8,title='VORTEX MARK at McMurdo',color=0,xtickname=mon,xticks=nmonth-1,psym=-2,charsize=1.25,yrange=[0,1.1]
oplot,1+findgen(nmonth),markavg40,thick=8,color=100,psym=-2
oplot,1+findgen(nmonth),markavg50,thick=8,color=250,psym=-2

xmn=xorig(5)
xmx=xorig(5)+xlen
ymn=yorig(5)
ymx=yorig(5)+ylen
set_viewport,xmn,xmx,ymn,ymx
plot,1+findgen(nmonth),distavg30,thick=8,title='DIST FROM EDGE',color=0,xtickname=mon,xticks=nmonth-1,psym=-2,charsize=1.25,yrange=[-4000,2000],ytitle='km'
oplot,1+findgen(nmonth),0.*tavg30,color=0
oplot,1+findgen(nmonth),distavg40,thick=8,color=100,psym=-2
oplot,1+findgen(nmonth),distavg50,thick=8,color=250,psym=-2

xyouts,0.15+min(xorig),0.675,'--',color=0,charsize=2,charthick=5,/normal
xyouts,0.15+min(xorig)+0.05,0.675,'30 km',color=0,charsize=2,charthick=2,/normal
xyouts,0.35+min(xorig),0.675,'--',color=100,charsize=2,charthick=5,/normal
xyouts,0.35+min(xorig)+0.05,0.675,'40 km',color=100,charsize=2,charthick=2,/normal
xyouts,0.55+min(xorig),0.675,'--',color=250,charsize=2,charthick=5,/normal
xyouts,0.55+min(xorig)+0.05,0.675,'50 km',color=250,charsize=2,charthick=2,/normal

if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim '+pthout+'annualcycle_merra2mcmurdo.ps -rotate -90 '+pthout+'annualcycle_merra2mcmurdo.jpg'
;  spawn,'/usr/bin/rm '+pthout+'annualcycle_merra2mcmurdo.ps'
endif
end
