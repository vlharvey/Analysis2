;
; timeseries of PDF peak geographical latitudes of high CO elats and low CO elats
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
yorig=[0.40]
cbaryoff=0.02
cbarydel=0.01
xlen=0.75
ylen=0.5
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
spawn,'ls daily_mls_coelatedge+merra_sfelatedge_??????_1800K.sav',ifiles
restore,ifiles(0)
hlat_time_all=highlat_time			; PDF Peak Geographical Lat of High Elats
llat_time_all=lowlat_time			; PDF Peak Geographical Lat of Low Elats
sfelatedge_time_all=sfelatedge_time		; Elat of poleward-most edge if there are two co-located QDF nodes and local max in speed (nc5)
elatedge_time_all=elatedge_time			; CO Elat Edge (absolute maximum in 1st derivative)
markcoelatedge_time_all=markcoelatedge_time	; CO Elat of Marker Edge
marksfelatedge_time_all=marksfelatedge_time	; SF Elat of Marker Edge
coedge_time_all=coedge_time			; CO Value at Edge
lowlat_elatedge_time_all=lowlat_elatedge_time	; CO Elat Edge (equatorward-most maximum in 1st derivative that is half the magnitude of absolute maximum) *** Use this as CO Edge ***
lowlat_coedge_time_all=lowlat_coedge_time	; CO Value at Lowlat Edge
pvelatedge_time_all=pvelatedge_time		; PV Value at Edge
nashedge_time_all=nashedge_time			; Elat of Nash Edge
hlatpdf_time_all=hlatpdf_time			; Full PDFs of Geographical Lats of High Elats > 80
llatpdf_time_all=llatpdf_time			; Full PDFs of Geographical Lats of Low Elats < 60
sdate_time_all=sdate_time

index=where(yeq mod 2 eq 0,nrr2)
yeq2=yeq(index)                         ; 2 degree bins

for ifile=1L,n_elements(ifiles)-1L do begin
    restore,ifiles(ifile)
    hlat_time_all=[hlat_time_all,highlat_time]
    llat_time_all=[llat_time_all,lowlat_time]
    sfelatedge_time_all=[sfelatedge_time_all,sfelatedge_time]
    elatedge_time_all=[elatedge_time_all,elatedge_time]
    markcoelatedge_time_all=[markcoelatedge_time_all,markcoelatedge_time]
    marksfelatedge_time_all=[marksfelatedge_time_all,marksfelatedge_time]
    coedge_time_all=[coedge_time_all,coedge_time]
    lowlat_elatedge_time_all=[lowlat_elatedge_time_all,lowlat_elatedge_time]
    lowlat_coedge_time_all=[lowlat_coedge_time_all,lowlat_coedge_time]
    pvelatedge_time_all=[pvelatedge_time_all,pvelatedge_time]
    nashedge_time_all=[nashedge_time_all,nashedge_time]
    hlatpdf_time_all=[hlatpdf_time_all,hlatpdf_time]
    llatpdf_time_all=[llatpdf_time_all,llatpdf_time]
    sdate_time_all=[sdate_time_all,sdate_time]
endfor
nday=n_elements(sdate_time_all)

; postscript file
;
if setplot eq 'ps' then begin
   lc=0
   xsize=nxdim/100.
   ysize=nydim/100.
   set_plot,'ps'
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='timeseries_daily_coelat_latpdf_'+spress+'.ps'
   !p.charsize=1.2
   !p.thick=2
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
;
; retain DJF
; retain Jan
;
index=where(strmid(sdate_time_all,4,2) eq '12' or strmid(sdate_time_all,4,2) eq '01' or strmid(sdate_time_all,4,2) eq '02',nday)
index=where(strmid(sdate_time_all,4,2) eq '01',nday)
hlat_time_all=hlat_time_all(index)
llat_time_all=llat_time_all(index)
sfelatedge_time_all=sfelatedge_time_all(index)
elatedge_time_all=elatedge_time_all(index)
markcoelatedge_time_all=markcoelatedge_time_all(index)
marksfelatedge_time_all=marksfelatedge_time_all(index)
coedge_time_all=coedge_time_all(index)
lowlat_elatedge_time_all=lowlat_elatedge_time_all(index)
lowlat_coedge_time_all=lowlat_coedge_time_all(index)
pvelatedge_time_all=pvelatedge_time_all(index)
nashedge_time_all=nashedge_time_all(index)
hlatpdf_time_all=hlatpdf_time_all(index,*)
llatpdf_time_all=llatpdf_time_all(index,*)
sdate_time_all=sdate_time_all(index)
;
; plot PDF maximum geographical Lat and Latitude range for High Elat values
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
xindex=where(strmid(sdate_time_all,4,4) eq '1201',nxticks)
xlabs0=strmid(sdate_time_all(xindex),2,2)
xindex=where(strmid(sdate_time_all,4,4) eq '0101',nxticks)
xlabs0=strarr(nxticks)
xlabs1=strmid(sdate_time_all(xindex),2,2)
xlabs=xlabs0+xlabs1
plot,findgen(nday),hlat_time_all,color=0,/noeras,charsize=1.5,ytitle='Geographical Latitude of COElat>80',title='January '+sth,yrange=[30,90],xtickname=xlabs,xtickv=xindex,xticks=nxticks-1,/nodata
oplot,findgen(nday),hlat_time_all,color=90,psym=8
goodflag=-99.+0.*fltarr(nday)
for i=0L,nday-1L do begin
    today=hlatpdf_time_all(i,*)
    index=where(today gt 0.)
    if min(yeq2(index)) gt 30. then plots,i,min(yeq2(index))
    if min(yeq2(index)) le 30. then plots,i,30.
    plots,i,max(yeq2(index)),/continue,color=90
;
; connect days off equator
;
    if i gt 0L then begin
       index=where(yesterday gt 0.)
print,'yesterday min ',min(yeq2(index))
       if min(yeq2(index)) gt 0. then begin
          index=where(today gt 0.)
print,'today min ',min(yeq(index))
          if min(yeq(index)) gt 0. then oplot,[i-1,i],[hlat_time_all(i-1),hlat_time_all(i)],color=60,thick=3
          if min(yeq(index)) gt 0. then goodflag(i)=1.
;if min(yeq(index)) gt 0. then stop
       endif  
    endif
    yesterday=today
endfor
oplot,findgen(nday),hlat_time_all,color=60,psym=8
;
; what fraction is good?
;
index=where(goodflag eq 1.)
print,100.*n_elements(index)/float(nday)
sval=strcompress(long(round(100.*n_elements(index)/float(nday))),/remove_all)
xyouts,xmx-0.075,ymn+0.02,sval+'%',/normal,color=0,charsize=2,charthick=2
;
; indicate new season
;
;xindex=where(strmid(sdate_time_all,4,4) eq '1201',nxticks)
xindex=where(strmid(sdate_time_all,4,4) eq '0101',nxticks)
for i=0L,nxticks-1L do begin
    plots,xindex(i),30
    plots,xindex(i),90,/continue,color=0
endfor
;
; superimpose Elat edges
;
axis,/yax,yrange=[30,90],/save,ytitle='Equivalent Latitude',charsize=1.5,color=mcolor*.9
;hlat_time_all=highlat_time                      ; PDF Peak Geographical Lat of High Elats
;llat_time_all=lowlat_time                       ; PDF Peak Geographical Lat of Low Elats
;sfelatedge_time_all=sfelatedge_time             ; Elat of poleward-most edge if there are two co-located QDF nodes and local max in speed (nc5)
;elatedge_time_all=elatedge_time                 ; CO Elat Edge (absolute maximum in 1st derivative)
;markcoelatedge_time_all=markcoelatedge_time     ; CO Elat of Marker Edge
;marksfelatedge_time_all=marksfelatedge_time     ; SF Elat of Marker Edge
;coedge_time_all=coedge_time                     ; CO Value at Edge
;lowlat_elatedge_time_all=lowlat_elatedge_time   ; CO Elat Edge (equatorward-most maximum in 1st derivative that is half the magnitude of absolute maximum) *** Use this as CO Edge ***
;lowlat_coedge_time_all=lowlat_coedge_time       ; CO Value at Lowlat Edge
;pvelatedge_time_all=pvelatedge_time             ; PV Value at Edge
;nashedge_time_all=nashedge_time                 ; Elat of Nash Edge
;hlatpdf_time_all=hlatpdf_time                   ; Full PDFs of Geographical Lats of High Elats > 80
;llatpdf_time_all=llatpdf_time                   ; Full PDFs of Geographical Lats of Low Elats < 60
oplot,findgen(nday),lowlat_elatedge_time_all,psym=8,color=250
oplot,findgen(nday),marksfelatedge_time_all,psym=8,color=30
;index=where(sfelatedge_time_all ne 0.)
;oplot,index,sfelatedge_time_all(index),psym=8,color=0		; nc5 fix for double jets

set_viewport,xmn,xmx,0.1,0.3
plot,findgen(nday),lowlat_elatedge_time_all,psym=8,color=0,ytitle='Elat',/noeras,xtickname=xlabs,xtickv=xindex,xticks=nxticks-1,/nodata,yrange=[0,90],charsize=1.5
oplot,findgen(nday),lowlat_elatedge_time_all,psym=8,color=250
oplot,findgen(nday),marksfelatedge_time_all,psym=8,color=30
;index=where(sfelatedge_time_all ne 0.)
;oplot,index,sfelatedge_time_all(index),psym=8,color=0
;
; what fraction is good?
;
goodflag=-99.+0.*fltarr(nday)
diff=abs(lowlat_elatedge_time_all-marksfelatedge_time_all)
index=where(diff le 5.)
print,100.*n_elements(index)/float(nday)
sval=strcompress(long(round(100.*n_elements(index)/float(nday))),/remove_all)
xyouts,xmx-0.75,0.1+0.02,sval+'% within 5!uo!n',/normal,color=0,charsize=1.5,charthick=2
index=where(diff le 10.)
print,100.*n_elements(index)/float(nday)
sval=strcompress(long(round(100.*n_elements(index)/float(nday))),/remove_all)
xyouts,xmx-0.55,0.1+0.02,sval+'% within 10!uo!n',/normal,color=0,charsize=1.5,charthick=2
index=where(diff le 15.)
print,100.*n_elements(index)/float(nday)
sval=strcompress(long(round(100.*n_elements(index)/float(nday))),/remove_all)
xyouts,xmx-0.35,0.1+0.02,sval+'% within 15!uo!n',/normal,color=0,charsize=1.5,charthick=2
index=where(diff le 20.)
print,100.*n_elements(index)/float(nday)
sval=strcompress(long(round(100.*n_elements(index)/float(nday))),/remove_all)
xyouts,xmx-0.15,0.1+0.02,sval+'% within 20!uo!n',/normal,color=0,charsize=1.5,charthick=2

;
; indicate new season
;
xindex=where(strmid(sdate_time_all,4,4) eq '0101',nxticks)
for i=0L,nxticks-1L do begin
    plots,xindex(i),0
    plots,xindex(i),90,/continue,color=0
endfor
;
; Close PostScript file and return control to X-windows
;
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim timeseries_daily_coelat_latpdf_'+spress+'.ps -rotate -90 timeseries_daily_coelat_latpdf_'+spress+'.jpg'
;  spawn,'rm -f timeseries_daily_coelat_latpdf_'+spress+'.ps'
endif
end
