;
; timeseries of MLS IWC derived from T and H2O using Hervig's equation
; NH version
;
device,decompose=0
!p.background=255
setplot='ps'
read,'setplot=',setplot
nxdim=750
nydim=750
xorig=[0.2]
yorig=[0.2]
xlen=0.7
ylen=0.5
cbaryoff=0.12
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=255
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smonth=['J','F','M','A','M','J','J','A','S','O','N','D']
col1=[0,1,2,3,4,5,10]
for lstyr=2007,2013 do begin
slstyr=strcompress(lstyr,/remove_all)
kday=92+31L
restore,'avg_mls_iwc_may_'+slstyr+'.sav'
icemay=ice
restore,'avg_mls_iwc_jun_'+slstyr+'.sav'
icejun=ice
if lstyr ne 2013 then begin
restore,'avg_mls_iwc_jul_'+slstyr+'.sav'
icejul=ice
restore,'avg_mls_iwc_aug_'+slstyr+'.sav'
iceaug=ice
endif
if lstyr eq 2013 then begin
icejul=0.*icejul
iceaug=0.*iceaug
endif
;
; create NDJF ICE
;
nr=n_elements(latbin)
nc=n_elements(lonbin)
icemjja=fltarr(nc,nr,kday)
icemjja(*,*,0:30)=icemay
icemjja(*,*,31:31+30-1)=icejun
icemjja(*,*,31+30:31+30+31-1)=icejul
icemjja(*,*,31+30+31:31+30+31+31-1)=iceaug
;
; convert sdate to DFS
;
doy=121+findgen(kday)
jdaysol=172.
dfs=doy-jdaysol
;
; latitude
;
if lstyr eq 2007 then begin
rlat=75.
print,latbin
read,'Enter latitude ',rlat
slat=strcompress(long(rlat),/remove_all)
index=where(latbin eq rlat)
ilat=index(0)
endif
xtice=reform(icemjja(*,ilat,*))
;
; zonal mean
;
zmice=fltarr(kday)
for k=0L,kday-1L do begin
    good=where(xtice(*,k) gt 0.,nn)
    if good(0) ne -1L then zmice(k)=mean(xtice(good,k))
endfor
;
if lstyr eq 2007 then begin
if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !p.font=0
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='timeseries_mls_iwc_'+slat+'.ps'
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
;
restore,'c11.tbl'
tvlct,c1,c2,c3
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7	; ticks outward
plot,dfs,zmice,/nodata,color=0,charsize=2,charthick=2,xrange=[-50.,80.],$
     ytitle='IWC (ug/m2)',yrange=[0.,210.],xtitle='DFS'
xyouts,xmn+0.02,ymx-0.04,slat+' Latitude',/normal,charsize=2,charthick=2,color=0
endif
print,max(zmice)
good=where(zmice ge 0.001)
zmice=smooth(zmice,5,/edge_truncate)
if lstyr eq 2013 then loadct,0 & col1(6)=150
oplot,dfs(good),zmice(good),color=col1(lstyr-2007),thick=10
xyouts,xmx-0.12,ymx-0.04-0.04*(lstyr-2007),slstyr,/normal,charsize=2,charthick=2,color=col1(lstyr-2007)
endfor

    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert -trim timeseries_mls_iwc_'+slat+'.ps -rotate -90 timeseries_mls_iwc_'+slat+'.png'
    endif
end
