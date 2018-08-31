;
; plot MERRA yearly files
;
loadct,0
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
device,decompose=0
!p.background=icolmax
a=findgen(8)*(2*!pi/8.)
usersym,1.5*cos(a),1.5*sin(a),/fill
nxdim=700
nydim=700
xorig=[0.15,.15]
yorig=[.55,.1]
xlen=0.7
ylen=0.35
cbaryoff=0.075
cbarydel=0.01
set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif
; save postscript version
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       device,font_size=9
       device,/landscape,bits=8,filename='merra_t_u_clim0+2012-2013.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
       !p.thick=2.0                   ;Plotted lines twice as thick
       !p.charsize=2.0
    endif

dirh='/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_press_'
model_years=1979+indgen(37)
model_years=string(FORMAT='(i4)',long(model_years))
nyears=n_elements(model_years)
;
; plot timeseries of polar cap temperature
;
erase
!type=2^2+2^3
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
nlvls=20
col1=1+mcolor*findgen(20)/nlvls
;
; restore first year
;
ifile=dirh+model_years(0)+'_daily_zm.sav'
restore,ifile   ;,sdate_all,latbin,pressure,UGRD_AVG,VGRD_AVG,GPGRD_AVG,H2OGRD_AVG,TGRD_AVG
nday=n_elements(sdate_all)
tsave=fltarr(nday,nyears)
usave=fltarr(nday,nyears)
datesave=strarr(nday,nyears)
rpress=10.
rlat=70.
ilev=where(pressure eq rpress)
ilat=where(min(abs(latbin-rlat)) eq abs(latbin-rlat))
tsave(*,0)=mean(reform(TGRD_AVG(ilat(0):-1,ilev,*)),dim=1)
rlat2=60.
ilat2=where(min(abs(latbin-rlat2)) eq abs(latbin-rlat2))
usave(*,0)=mean(reform(UGRD_AVG(ilat2(0),ilev,*)),dim=1)
for iyear=1L,nyears-1L do begin
    ifile=dirh+model_years(iyear)+'_daily_zm.sav'
    restore,ifile	;,sdate_all,latbin,pressure,UGRD_AVG,VGRD_AVG,GPGRD_AVG,H2OGRD_AVG,TGRD_AVG
    plotarray=mean(reform(TGRD_AVG(ilat(0):-1,ilev,*)),dim=1)
    index=where(plotarray eq 0.)
    if index(0) ne -1L then plotarray(index)=0./0.
    tsave(*,iyear)=plotarray
    plotarray=reform(UGRD_AVG(ilat2(0),ilev,*))
    index=where(plotarray eq 0.)
    if index(0) ne -1L then plotarray(index)=0./0.
    usave(*,iyear)=plotarray
    datesave(*,iyear)=model_years(iyear)+sdate_all
endfor	; loop over years
tmin=fltarr(nday)
tmax=fltarr(nday)
tmean=fltarr(nday)
tsig=fltarr(nday)
umin=fltarr(nday)
umax=fltarr(nday)
umean=fltarr(nday)
usig=fltarr(nday)
for i=0,nday-1 do begin
    tvals=reform(tsave(i,*))
    index=where(tvals ne 0. and finite(tvals) eq 1)
    if index(0) ne -1L then tmin(i)=min(tvals(index))
    if index(0) ne -1L then tmax(i)=max(tvals(index))
    if index(0) ne -1L then tmean(i)=mean(tvals(index))
    if index(0) ne -1L then tsig(i)=stdev(tvals(index))
    uvals=reform(usave(i,*))
    index=where(uvals ne 0. and finite(uvals) eq 1)
    if index(0) ne -1L then umin(i)=min(uvals(index))
    if index(0) ne -1L then umax(i)=max(uvals(index))
    if index(0) ne -1L then umean(i)=mean(uvals(index))
    if index(0) ne -1L then usig(i)=stdev(uvals(index))
endfor
;
; extract DJF climo and D2012 and JF2013
;
nday=90
xindex=where(strmid(sdate_all,0,2) eq '12')	; Dec
tmin=[tmin(xindex),tmin(0:nday-31-1)]
tmax=[tmax(xindex),tmax(0:nday-31-1)]
tmean=[tmean(xindex),tmean(0:nday-31-1)]
tsig=[tsig(xindex),tsig(0:nday-31-1)]
umin=[umin(xindex),umin(0:nday-31-1)]
umax=[umax(xindex),umax(0:nday-31-1)]
umean=[umean(xindex),umean(0:nday-31-1)]
usig=[usig(xindex),usig(0:nday-31-1)]
plot,1+findgen(nday),tmean,xrange=[1,nday],yrange=[190,265],color=0,charsize=2,charthick=2,/noeras,ytitle='>70N Temperature (K)',thick=3,title='MERRA '+strcompress(long(rpress),/r)+' hPa',$
     xticks=4,xtickv=[15,31,15+31,31+31,15+31+28],xtickname=['Dec',' ','Jan',' ','Feb']
tmax=smooth(tmax,3,/edge_truncate)
tmin=smooth(tmin,3,/edge_truncate)
for i=1,nday-1 do begin
    if i lt nday-1 then begin
    plots,i+1,tmin(i)
    plots,i+1,tmax(i),color=200,/continue,thick=15
    plots,i+1,tmean(i)-tsig(i)
    plots,i+1,tmean(i)+tsig(i),color=100,/continue,thick=15
    endif
endfor
oplot,1+findgen(nday),tmean,color=0,thick=8
oplot,1+findgen(nday),tmin,color=0,thick=3
oplot,1+findgen(nday),tmax,color=0,thick=3
loadct,39
ifile=dirh+'2012_daily_zm.sav'
restore,ifile   ;,sdate_all,latbin,pressure,UGRD_AVG,VGRD_AVG,GPGRD_AVG,H2OGRD_AVG,TGRD_AVG
plotarray=reform(TGRD_AVG(ilat(0),ilev,*))
plotarray=plotarray(xindex)	; Dec 2012
ifile=dirh+'2013_daily_zm.sav'
restore,ifile   ;,sdate_all,latbin,pressure,UGRD_AVG,VGRD_AVG,GPGRD_AVG,H2OGRD_AVG,TGRD_AVG
plotarray2=reform(TGRD_AVG(ilat(0),ilev,*))
plotarray=[plotarray,plotarray2(0:nday-31-1)]
oplot,1+findgen(nday),plotarray,color=0.9*mcolor,thick=10
xyouts,xmn+0.02,ymx-0.03,'2012-2013',charsize=2,charthick=2,color=.9*mcolor,/normal
xyouts,xmx-0.2,ymx-0.03,'1979-2015',charsize=2,charthick=2,color=0,/normal

!type=2^2+2^3
xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
nlvls=20
col1=1+mcolor*findgen(20)/nlvls
loadct,0
plot,1+findgen(nday),umean,xrange=[1,nday],yrange=[-30,80],color=0,charsize=2,charthick=2,/noeras,ytitle='60N Ubar (m/s)',thick=3,$
     xticks=4,xtickv=[15,31,15+31,31+31,15+31+28],xtickname=['Dec',' ','Jan',' ','Feb']
umax=smooth(umax,3,/edge_truncate)
umin=smooth(umin,3,/edge_truncate)
for i=1,nday-1 do begin
    if i lt nday-1 then begin
    plots,i+1,umin(i)
    plots,i+1,umax(i),color=200,/continue,thick=15
    plots,i+1,umean(i)-usig(i)
    plots,i+1,umean(i)+usig(i),color=100,/continue,thick=15
    endif
endfor
oplot,1+findgen(nday),umean,color=0,thick=8
oplot,1+findgen(nday),umin,color=0,thick=3
oplot,1+findgen(nday),umax,color=0,thick=3
loadct,39
ifile=dirh+'2012_daily_zm.sav'
restore,ifile   ;,sdate_all,latbin,pressure,UGRD_AVG,VGRD_AVG,GPGRD_AVG,H2OGRD_AVG,TGRD_AVG
plotarray=reform(UGRD_AVG(ilat2(0),ilev,*))
plotarray=plotarray(xindex)     ; Dec 2012
ifile=dirh+'2013_daily_zm.sav'
restore,ifile   ;,sdate_all,latbin,pressure,UGRD_AVG,VGRD_AVG,GPGRD_AVG,H2OGRD_AVG,TGRD_AVG
plotarray2=reform(UGRD_AVG(ilat2(0),ilev,*))
plotarray=[plotarray,plotarray2(0:nday-31-1)]
oplot,1+findgen(nday),plotarray,color=0.9*mcolor,thick=10
xyouts,xmn+0.02,ymx-0.03,'2012-2013',charsize=2,charthick=2,color=.9*mcolor,/normal
xyouts,xmx-0.2,ymx-0.03,'1979-2015',charsize=2,charthick=2,color=0,/normal

    if setplot eq 'ps' then begin
       device,/close
       spawn,'convert -trim merra_t_u_clim0+2012-2013.ps -rotate -90 merra_t_u_clim0+2012-2013.png'
    endif

end
