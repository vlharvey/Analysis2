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
tsave=fltarr(n_elements(sdate_all),n_elements(model_years))
rpress=10.
rlat=70.
ilev=where(pressure eq rpress)
ilat=where(min(abs(latbin-rlat)) eq abs(latbin-rlat))
nday=n_elements(sdate_all)
plotarray=mean(reform(TGRD_AVG(ilat(0):-1,ilev,*)),dim=1)
tsave(*,0)=plotarray
plot,1+findgen(nday),plotarray,xrange=[1,92],yrange=[190,265],color=0,charsize=2,charthick=2,/noeras,ytitle='>70N Temperature (K)',xtitle='DOY',thick=3,title='MERRA '+strcompress(long(rpress),/r)+' hPa'
for iyear=1L,nyears-1L do begin
    ifile=dirh+model_years(iyear)+'_daily_zm.sav'
    restore,ifile	;,sdate_all,latbin,pressure,UGRD_AVG,VGRD_AVG,GPGRD_AVG,H2OGRD_AVG,TGRD_AVG
    plotarray=mean(reform(TGRD_AVG(ilat(0):-1,ilev,*)),dim=1)
    oplot,findgen(nday),plotarray,color=0,thick=3	;(float(iyear)/float(nyears))*mcolor,thick=3
    tsave(*,iyear)=plotarray
endfor	; loop over years
loadct,39
ifile=dirh+'2013_daily_zm.sav'
restore,ifile   ;,sdate_all,latbin,pressure,UGRD_AVG,VGRD_AVG,GPGRD_AVG,H2OGRD_AVG,TGRD_AVG
plotarray=mean(reform(TGRD_AVG(ilat(0):-1,ilev,*)),dim=1)
loadct,0
tmin=fltarr(nday)
tmax=fltarr(nday)
for i=0,nday-1 do begin
    tmin(i)=min(tsave(i,*))
    tmax(i)=max(tsave(i,*))
endfor
for i=0,91 do begin
    plots,i+1,min(tsave(i,*))
    plots,i+1,max(tsave(i,*)),color=200,/continue,thick=8
endfor
oplot,findgen(nday),tmin,color=0,thick=5
oplot,findgen(nday),tmax,color=0,thick=5
loadct,39
oplot,findgen(nday),plotarray,color=0.9*mcolor,thick=10
xyouts,xmx-0.12,ymn+0.01,'2013',charsize=3,charthick=2,color=.9*mcolor,/normal
xyouts,xmx-0.2,ymx-0.03,'1979-2015',charsize=2,charthick=2,color=0,/normal

!type=2^2+2^3
xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
nlvls=20
col1=1+mcolor*findgen(20)/nlvls
;
; restore first year
;
ifile=dirh+model_years(0)+'_daily_zm.sav'
restore,ifile   ;,sdate_all,latbin,pressure,UGRD_AVG,VGRD_AVG,GPGRD_AVG,H2OGRD_AVG,TGRD_AVG
usave=fltarr(n_elements(sdate_all),n_elements(model_years))
rpress=10.
rlat=60.
ilev=where(pressure eq rpress)
ilat=where(min(abs(latbin-rlat)) eq abs(latbin-rlat))
nday=n_elements(sdate_all)
plotarray=reform(UGRD_AVG(ilat(0),ilev,*))
usave(*,0)=plotarray
plot,1+findgen(nday),plotarray,xrange=[1,92],yrange=[-40,80],color=0,charsize=2,charthick=2,/noeras,ytitle='60N Ubar (K)',xtitle='DOY',thick=3
for iyear=1L,nyears-1L do begin
    ifile=dirh+model_years(iyear)+'_daily_zm.sav'
    restore,ifile       ;,sdate_all,latbin,pressure,UGRD_AVG,VGRD_AVG,GPGRD_AVG,H2OGRD_AVG,TGRD_AVG
    plotarray=reform(UGRD_AVG(ilat(0),ilev,*))
    oplot,findgen(nday),plotarray,color=0,thick=3       ;(float(iyear)/float(nyears))*mcolor,thick=3
    usave(*,iyear)=plotarray
endfor  ; loop over years
loadct,39
ifile=dirh+'2013_daily_zm.sav'
restore,ifile   ;,sdate_all,latbin,pressure,UGRD_AVG,VGRD_AVG,GPGRD_AVG,H2OGRD_AVG,TGRD_AVG
plotarray=reform(UGRD_AVG(ilat(0),ilev,*))
loadct,0
umin=fltarr(nday)
umax=fltarr(nday)
for i=0,nday-1 do begin
    umin(i)=min(usave(i,*))
    umax(i)=max(usave(i,*))
endfor
for i=0,91 do begin
    plots,i+1,min(usave(i,*))
    plots,i+1,max(usave(i,*)),color=200,/continue,thick=8
endfor
oplot,findgen(nday),umin,color=0,thick=5
oplot,findgen(nday),umax,color=0,thick=5
loadct,39
oplot,findgen(nday),plotarray,color=0.9*mcolor,thick=10
xyouts,xmx-0.12,ymn+0.01,'2013',charsize=3,charthick=2,color=.9*mcolor,/normal
xyouts,xmx-0.2,ymx-0.03,'1979-2015',charsize=2,charthick=2,color=0,/normal

end
