;
; timeseries of daily correlation between GEOS and MLS and SABER sabertratopause temperaure
; separate inside vortex from inside anticyclones
;
@sabertddat
@kgeosmarklsabert
@ckday
@kdate

loadct,39
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
device,decompose=0
!p.background=icolmax
setplot='ps'
read,'setplot=',setplot
nxdim=600
nydim=600
xorig=[0.20,0.2,0.2]
yorig=[0.7,0.4,0.1]
xlen=0.65
ylen=0.175
cbaryoff=0.05
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
smonth=['J','A','S','O','N','D','J','F','M','A','M','J']
;
; restore save file
;
; DATES           STRING    = Array[1800]
; GEOS5HEIGHTS    FLOAT     = Array[1800, 96, 72]
; GEOS5LAT        FLOAT     = Array[72]
; GEOS5LON        FLOAT     = Array[96]
; GEOS5MARK       FLOAT     = Array[96, 1800]
; GEOS5MARKERS    FLOAT     = Array[1800, 96, 72]
; GEOS5PTS        FLOAT     = Array[96, 1800]
; GEOS5TEMPS      FLOAT     = Array[1800, 96, 72]
; GEOS5TS         FLOAT     = Array[96, 1800]
; MLSHEIGHTS      FLOAT     = Array[1800, 96, 72]
; MLSHEIGHTS1     FLOAT     = Array[1800, 20, 40]
; MLSLAT          FLOAT     = Array[40]
; MLSLON          FLOAT     = Array[20]
; MLSPTS          FLOAT     = Array[1800, 96, 72]
; MLSTEMPS        FLOAT     = Array[1800, 96, 72]
; MLSTEMPS1       FLOAT     = Array[1800, 20, 40]
; MLSTS           FLOAT     = Array[1800, 96, 72]
; SABERHEIGHTS    FLOAT     = Array[1800, 96, 72]
; SABERHEIGHTS1   FLOAT     = Array[1800, 20, 40]
; SABERLAT        FLOAT     = Array[40]
; SABERLON        FLOAT     = Array[20]
; SABERPTS        FLOAT     = Array[1800, 96, 72]
; SABERTEMPS      FLOAT     = Array[1800, 96, 72]
; SABERTEMPS1     FLOAT     = Array[1800, 20, 40]
; SABERTS         FLOAT     = Array[1800, 96, 72]
; 
restore,'strat_Z_T_MLS_GEOS_SABER.sav'
nday=n_elements(dates)
;nday=100	; testing
nc=n_elements(geos5lon)
nr=n_elements(geos5lat)
geos_mls_vcorr_nh=-99.+0.*fltarr(nday)
geos_mls_anticorr_nh=-99.+0.*fltarr(nday)
geos_mls_ambcorr_nh=-99.+0.*fltarr(nday)
geos_mls_vcorr_sh=-99.+0.*fltarr(nday)
geos_mls_anticorr_sh=-99.+0.*fltarr(nday)
geos_mls_ambcorr_sh=-99.+0.*fltarr(nday)
geos_saber_vcorr_nh=-99.+0.*fltarr(nday)
geos_saber_anticorr_nh=-99.+0.*fltarr(nday)
geos_saber_ambcorr_nh=-99.+0.*fltarr(nday)
geos_saber_vcorr_sh=-99.+0.*fltarr(nday)
geos_saber_anticorr_sh=-99.+0.*fltarr(nday)
geos_saber_ambcorr_sh=-99.+0.*fltarr(nday)
mls_saber_vcorr_nh=-99.+0.*fltarr(nday)
mls_saber_anticorr_nh=-99.+0.*fltarr(nday)
mls_saber_ambcorr_nh=-99.+0.*fltarr(nday)
mls_saber_vcorr_sh=-99.+0.*fltarr(nday)
mls_saber_anticorr_sh=-99.+0.*fltarr(nday)
mls_saber_ambcorr_sh=-99.+0.*fltarr(nday)
lat2d=fltarr(nc,nr)
for i=0L,nc-1L do lat2d(i,*)=GEOS5LAT
;
; loop over days
;
for iday=0,nday-1 do begin
;
; GEOS MLS correlations
;
geost=reform(geos5temps(iday,*,*))
mlst=reform(mlstemps(iday,*,*))
sabert=reform(sabertemps(iday,*,*))
geosmark=reform(geos5markers(iday,*,*))
index=where(geost ne 0. and finite(mlst) eq 1 and geosmark eq 1. and lat2d gt 0.)
if index(0) ne -1L then begin
result=correlate(mlst(index),geost(index))
geos_mls_vcorr_nh(iday)=result(0)
if result(0) gt 1. or result(0) lt -1. then stop
endif

index=where(geost ne 0. and finite(mlst) eq 1 and geosmark lt 0. and lat2d gt 0.)
if index(0) ne -1L then begin
result=correlate(mlst(index),geost(index))                            
geos_mls_anticorr_nh(iday)=result(0)
if result(0) gt 1. or result(0) lt -1. then stop
endif

index=where(geost ne 0. and finite(mlst) eq 1 and geosmark eq 0. and lat2d gt 30.)
if index(0) ne -1L then begin
result=correlate(mlst(index),geost(index))
geos_mls_ambcorr_nh(iday)=result(0)
if result(0) gt 1. or result(0) lt -1. then stop
endif

index=where(geost ne 0. and finite(mlst) eq 1 and geosmark eq 1. and lat2d lt 0.)
if index(0) ne -1L then begin
result=correlate(mlst(index),geost(index))
geos_mls_vcorr_sh(iday)=result(0)
if result(0) gt 1. or result(0) lt -1. then stop
endif

index=where(geost ne 0. and finite(mlst) eq 1 and geosmark lt 0. and lat2d lt 0.)
if index(0) ne -1L then begin
result=correlate(mlst(index),geost(index))
geos_mls_anticorr_sh(iday)=result(0)
if result(0) gt 1. or result(0) lt -1. then stop
endif

index=where(geost ne 0. and finite(mlst) eq 1 and geosmark eq 0. and lat2d lt -30.)
if index(0) ne -1L then begin
result=correlate(mlst(index),geost(index))
geos_mls_ambcorr_sh(iday)=result(0)
if result(0) gt 1. or result(0) lt -1. then stop
endif
print,dates(iday),geos_mls_vcorr_nh(iday),geos_mls_anticorr_nh(iday),geos_mls_ambcorr_nh(iday),geos_mls_vcorr_sh(iday),geos_mls_anticorr_sh(iday),geos_mls_ambcorr_sh(iday)
;
; GEOS SABER correlations
;
index=where(geost ne 0. and finite(sabert) eq 1 and geosmark eq 1. and lat2d gt 0.)
if index(0) ne -1L then begin
result=correlate(sabert(index),geost(index))
geos_saber_vcorr_nh(iday)=result(0)
if result(0) gt 1. or result(0) lt -1. then stop
endif

index=where(geost ne 0. and finite(sabert) eq 1 and geosmark lt 0. and lat2d gt 0.)
if index(0) ne -1L then begin
result=correlate(sabert(index),geost(index))
geos_saber_anticorr_nh(iday)=result(0)
if result(0) gt 1. or result(0) lt -1. then stop
endif

index=where(geost ne 0. and finite(sabert) eq 1 and geosmark eq 0. and lat2d gt 30.)
if index(0) ne -1L then begin
result=correlate(sabert(index),geost(index))
geos_saber_ambcorr_nh(iday)=result(0)
if result(0) gt 1. or result(0) lt -1. then stop
endif

index=where(geost ne 0. and finite(sabert) eq 1 and geosmark eq 1. and lat2d lt 0.)
if index(0) ne -1L then begin
result=correlate(sabert(index),geost(index))
geos_saber_vcorr_sh(iday)=result(0)
if result(0) gt 1. or result(0) lt -1. then stop
endif

index=where(geost ne 0. and finite(sabert) eq 1 and geosmark lt 0. and lat2d lt 0.)
if index(0) ne -1L then begin
result=correlate(sabert(index),geost(index))
geos_saber_anticorr_sh(iday)=result(0)
if result(0) gt 1. or result(0) lt -1. then stop
endif

index=where(geost ne 0. and finite(sabert) eq 1 and geosmark eq 0. and lat2d lt -30.)
if index(0) ne -1L then begin
result=correlate(sabert(index),geost(index))
geos_saber_ambcorr_sh(iday)=result(0)
if result(0) gt 1. or result(0) lt -1. then stop
endif
;
; MLS SABER correlations
;
index=where(finite(mlst) eq 1 and finite(sabert) eq 1 and geosmark eq 1. and lat2d gt 0.)
if index(0) ne -1L then begin
result=correlate(mlst(index),sabert(index))
mls_saber_vcorr_nh(iday)=result(0)
if result(0) gt 1. or result(0) lt -1. then stop
endif

index=where(finite(mlst) eq 1 and finite(sabert) eq 1 and geosmark lt 0. and lat2d gt 0.)
if index(0) ne -1L then begin
result=correlate(mlst(index),sabert(index))
mls_saber_anticorr_nh(iday)=result(0)
if result(0) gt 1. or result(0) lt -1. then stop
endif

index=where(finite(mlst) eq 1 and finite(sabert) eq 1 and geosmark eq 0. and lat2d gt 30.)
if index(0) ne -1L then begin
result=correlate(mlst(index),sabert(index))
mls_saber_ambcorr_nh(iday)=result(0)
if result(0) gt 1. or result(0) lt -1. then stop
endif

index=where(finite(mlst) eq 1 and finite(sabert) eq 1 and geosmark eq 1. and lat2d lt 0.)
if index(0) ne -1L then begin
result=correlate(mlst(index),sabert(index))
mls_saber_vcorr_sh(iday)=result(0)
if result(0) gt 1. or result(0) lt -1. then stop
endif

index=where(finite(mlst) eq 1 and finite(sabert) eq 1 and geosmark lt 0. and lat2d lt 0.)
if index(0) ne -1L then begin
result=correlate(mlst(index),sabert(index))
mls_saber_anticorr_sh(iday)=result(0)
if result(0) gt 1. or result(0) lt -1. then stop
endif

index=where(finite(mlst) eq 1 and finite(sabert) eq 1 and geosmark eq 0. and lat2d lt -30.)
if index(0) ne -1L then begin
result=correlate(mlst(index),sabert(index))
mls_saber_ambcorr_sh(iday)=result(0)
if result(0) gt 1. or result(0) lt -1. then stop
endif

endfor
;
; plot timeseries of correlation coefficients
;
if setplot eq 'ps' then begin
   xsize=nxdim/100.
   ysize=nydim/100.
   set_plot,'ps'
   !p.font=0
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='timeseries_sabertratopause_corr.ps'
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=5
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
xyouts,.1,.95,'Arctic Stratopause Temperature',/normal,color=0,charsize=2
imin=-1.
imax=1.
index=where(strmid(dates,4,4) eq '0101' or strmid(dates,4,4) eq '0701',nxticks)
xlabs=strmid(dates(index),4,2)
plot,findgen(nday),geos_mls_ambcorr_nh,title='GEOS vs. MLS',yrange=[imin,imax],psym=8,$
     color=0,symsize=1.5,xticks=nxticks-1,xtickv=index,xtickname=xlabs
plots,0,0
plots,nday-1,0,/continue,color=0
oplot,findgen(nday),geos_mls_vcorr_nh,psym=8,color=0.3*mcolor
oplot,findgen(nday),geos_mls_anticorr_nh,psym=8,color=0.9*mcolor,symsize=0.75
oplot,findgen(nday),geos_mls_ambcorr_nh,psym=8,color=0,symsize=0.25

;oplot,findgen(nday),geos_mls_ambcorr_sh,psym=8,color=0.2*mcolor
;oplot,findgen(nday),geos_mls_vcorr_sh,psym=4,color=0.3*mcolor
;oplot,findgen(nday),geos_mls_anticorr_sh,psym=4,color=0.9*mcolor

xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
imin=-1.
imax=1.
plot,findgen(nday),geos_saber_ambcorr_nh,title='GEOS vs. SABER',yrange=[imin,imax],psym=8,$
     ytitle='Correlation Coefficient',color=0,symsize=1.5,xticks=nxticks-1,xtickv=index,xtickname=xlabs
plots,0,0
plots,nday-1,0,/continue,color=0
oplot,findgen(nday),geos_saber_vcorr_nh,psym=8,color=0.3*mcolor
oplot,findgen(nday),geos_saber_anticorr_nh,psym=8,color=0.9*mcolor,symsize=0.75
oplot,findgen(nday),geos_saber_ambcorr_nh,psym=8,color=0,symsize=0.25

;oplot,findgen(nday),geos_saber_ambcorr_sh,psym=8,color=0.2*mcolor
;oplot,findgen(nday),geos_saber_vcorr_sh,psym=4,color=0.3*mcolor
;oplot,findgen(nday),geos_saber_anticorr_sh,psym=4,color=0.9*mcolor

xmn=xorig(2)
xmx=xorig(2)+xlen
ymn=yorig(2)
ymx=yorig(2)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
imin=-1.
imax=1.
plot,findgen(nday),mls_saber_ambcorr_nh,title='MLS vs. SABER',yrange=[imin,imax],psym=8,$
     color=0,symsize=1.5,xticks=nxticks-1,xtickv=index,xtickname=xlabs
plots,0,0
plots,nday-1,0,/continue,color=0
oplot,findgen(nday),mls_saber_vcorr_nh,psym=8,color=0.3*mcolor
oplot,findgen(nday),mls_saber_anticorr_nh,psym=8,color=0.9*mcolor,symsize=0.75
oplot,findgen(nday),mls_saber_ambcorr_nh,psym=8,color=0,symsize=0.25

;oplot,findgen(nday),mls_saber_ambcorr_sh,psym=8,color=0.2*mcolor
;oplot,findgen(nday),mls_saber_vcorr_sh,psym=4,color=0.3*mcolor
;oplot,findgen(nday),mls_saber_anticorr_sh,psym=4,color=0.9*mcolor

index=where(strmid(dates,4,4) eq '0201',nxticks)
syears=strmid(dates(index),0,4)
for i=0,nxticks-1L do xyouts,index(i),-2,syears(i),/data,color=0,alignment=0.

if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim timeseries_sabertratopause_corr.ps -rotate -90 timeseries_sabertratopause_corr.jpg'
endif
end
