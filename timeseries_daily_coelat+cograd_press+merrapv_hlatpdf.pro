;
; DJF
; timeseries of Elat edge values based on CO and PV
;
re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd
nrr=91L
yeq=findgen(nrr)
;index=where(yeq mod 2 eq 0,nrr)
;yeq=yeq(index)				; 2 degree bins

latcircle=fltarr(nrr)
hem_frac=fltarr(nrr)
for j=0,nrr-2 do begin
    hy=re*dtr
    dx=re*cos(yeq(j)*dtr)*360.*dtr
    latcircle(j)=dx*hy
endfor
for j=0,nrr-1 do begin
    if yeq(j) ge 0. then index=where(yeq ge yeq(j))
    if index(0) ne -1 then hem_frac(j)=100.*total(latcircle(index))/hem_area
    if yeq(j) eq 0. then hem_frac(j)=100.
endfor

loadct,39
mcolor=byte(!p.color)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
nxdim=700
nydim=700
xorig=[0.15]
yorig=[0.35]
cbaryoff=0.02
cbarydel=0.01
xlen=0.8
ylen=0.4
PI2=6.2831853071796
DTR=PI2/360.
RADEA=6.37E6
!NOERAS=-1
syear=['2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014']
nyear=n_elements(syear)
smon=['01','02','03','04','05','06','07','08','09','10','11','12']
nmon=n_elements(smon)
set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif
;
; get file listing
;
spawn,'ls daily_mls_coelatedge+merra_pvelatedge_??????_2.15443hPa.sav',ifiles
restore,ifiles(0)       ;'daily_mls_coelatedge+merra_pvelatedge_'+syear(iyear)+smon(imon)+'_2.15443hPa.sav'     ;,elatedge_time,pvelatedge_time,nashedge_time,sdate_time,spress,sth,hlatpdf_time,llatpdf_time
hlatpdf_time_all=hlatpdf_time
llatpdf_time_all=llatpdf_time
elatedge_time_all=elatedge_time
lowlat_elatedge_time_all=LOWLAT_ELATEDGE_TIME
pvelatedge_time_all=pvelatedge_time
nashedge_time_all=nashedge_time
sdate_time_all=sdate_time

for ifile=1L,n_elements(ifiles)-1L do begin
    restore,ifiles(ifile)
    hlatpdf_time_all=[hlatpdf_time_all,hlatpdf_time]
    llatpdf_time_all=[llatpdf_time_all,llatpdf_time]
    elatedge_time_all=[elatedge_time_all,elatedge_time]
    LOWLAT_ELATEDGE_TIME_all=[LOWLAT_ELATEDGE_TIME_all,LOWLAT_ELATEDGE_TIME]
    pvelatedge_time_all=[pvelatedge_time_all,pvelatedge_time]
    nashedge_time_all=[nashedge_time_all,nashedge_time]
    sdate_time_all=[sdate_time_all,sdate_time]
endfor

; postscript file
;
if setplot eq 'ps' then begin
   lc=0
   xsize=nxdim/100.
   ysize=nydim/100.
   set_plot,'ps'
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='timeseries_daily_coelat+cograd_press_merrapv_'+spress+'_hlatpdf.ps'
   !p.charsize=1.2
   !p.thick=2
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
;
; 2004-2005
;
index=where(sdate_time_all ge '20041101' and sdate_time_all le '20050331')
if index(0) ne -1L then begin
   sdate_time_all=sdate_time_all(index)
   hlatpdf_time_all=hlatpdf_time_all(index,*)
   llatpdf_time_all=llatpdf_time_all(index,*)
   LOWLAT_ELATEDGE_TIME_all=LOWLAT_ELATEDGE_TIME_all(index)
endif
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
xindex=where(strmid(sdate_time_all,6,2) eq '15',nxticks)
xlabs=strmid(sdate_time_all(xindex),4,2)
contour,hlatpdf_time_all,findgen(n_elements(sdate_time_all)),yeq,color=0,/noeras,charsize=1.5,ytitle='Geographic Latitude',yrange=[0,90],xtickname=xlabs,xtickv=xindex,xticks=nxticks-1,/nodata
for i=0L,n_elements(sdate_time_all)-1L do begin
    hlatday=reform(hlatpdf_time_all(i,*))
    index=where(hlatday gt 0.)
    if index(0) ne -1L then oplot,i+0.*index,yeq(index),psym=8,color=0
endfor
oplot,findgen(n_elements(sdate_time_all)),LOWLAT_ELATEDGE_TIME_all,psym=8,color=mcolor*.9,symsize=1.5
for i=0L,n_elements(sdate_time_all)-1L do begin
    llatday=reform(llatpdf_time_all(i,*))
    index=where(llatday gt 0.)
;    if index(0) ne -1L then oplot,i+0.*index,yeq(index),psym=8,color=mcolor*.3,symsize=0.5
endfor
;
; indicate new December
;
xindex=where(strmid(sdate_time_all,4,4) eq '1201',nxticks)
;for i=0L,nxticks-1L do begin
;    plots,xindex(i),20
;    plots,xindex(i),90,/continue,color=0
;endfor
;
; Difference between CO Elat and PV Elat
;
;xmn=xorig(1)
;xmx=xorig(1)+xlen
;ymn=yorig(1)
;ymx=yorig(1)+ylen
;set_viewport,xmn,xmx,ymn,ymx
;!type=2^2+2^3
;xindex=where(strmid(sdate_time_all,4,4) eq '1201',nxticks)
;xlabs0=strmid(sdate_time_all(xindex),2,2)
;xindex=where(strmid(sdate_time_all,4,4) eq '0101',nxticks)
;xlabs1=strmid(sdate_time_all(xindex),2,2)
;xlabs=xlabs0+xlabs1
;diff=LOWLAT_ELATEDGE_TIME_ALL-pvelatedge_time_all
;ymax=60
;ymin=-20
;plot,findgen(n_elements(sdate_time_all)),diff,color=0,/noeras,charsize=1.5,ytitle='CO Elat Edge minus PV Elat Edge',yrange=[ymin,ymax],psym=8,xtickname=xlabs,xtickv=xindex,xticks=nxticks-1,symsize=1.25
;oplot,findgen(n_elements(sdate_time_all)),0.*diff,color=0
;
; indicate new December
;
xindex=where(strmid(sdate_time_all,4,4) eq '1201',nxticks)
;for i=0L,nxticks-1L do begin
;    plots,xindex(i),ymin
;    plots,xindex(i),ymax,/continue,color=0,/data
;endfor

; Close PostScript file and return control to X-windows
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim timeseries_daily_coelat+cograd_press_merrapv_'+spress+'_hlatpdf.ps -rotate -90 '+$
                       'timeseries_daily_coelat+cograd_press_merrapv_'+spress+'_hlatpdf.jpg'
;  spawn,'rm -f timeseries_daily_coelat+cograd_press_merrapv_'+spress+'_hlatpdf.ps'
endif
end
