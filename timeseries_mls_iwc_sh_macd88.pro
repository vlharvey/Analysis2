;
; timeseries of MLS IWC derived from T and H2O using Hervig's equation
; SH version
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
col1=[0,1,3,4,6,11]
for lstyr=2007,2012 do begin
slstyr=strcompress(lstyr,/remove_all)
ledyr=lstyr+1L
sledyr=strcompress(ledyr,/remove_all)
lday=28
if ledyr mod 4 eq 0 then lday=29
kday=92L+lday
restore,'avg_mls_iwc_nov_'+slstyr+'.sav'
icenov=ice
restore,'avg_mls_iwc_dec_'+slstyr+'.sav'
icedec=ice
restore,'avg_mls_iwc_jan_'+sledyr+'.sav'
icejan=ice
restore,'avg_mls_iwc_feb_'+sledyr+'.sav'
icefeb=ice
;
; create NDJF ICE
;
nr=n_elements(latbin)
nc=n_elements(lonbin)
icendjf=fltarr(nc,nr,kday)
icendjf(*,*,0:29)=icenov
icendjf(*,*,30:30+31-1)=icedec
icendjf(*,*,30+31:30+31+31-1)=icejan
icendjf(*,*,30+31+31:30+31+31+lday-1)=icefeb
;
; convert sdate to DFS
;
doy=305+findgen(121)
jdaysol=355.
dfs=doy-jdaysol
;
; latitude
;
if lstyr eq 2007 then begin
rlat=-75.
print,latbin
read,'Enter latitude ',rlat
slat=strcompress(long(rlat),/remove_all)
index=where(latbin eq rlat)
ilat=index(0)
endif
xtice=reform(icendjf(*,ilat,*))
;
; zonal mean
;
zmice=fltarr(kday)
for k=0L,kday-1L do begin
    good=where(xtice(*,k) ne 0.,nn)
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
zmice=smooth(zmice,5,/edge_truncate)
oplot,dfs,zmice,color=col1(lstyr-2007),thick=10
xyouts,xmx-0.12,ymx-0.04-0.04*(lstyr-2007),slstyr,/normal,charsize=2,charthick=2,color=col1(lstyr-2007)
endfor

    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert -trim timeseries_mls_iwc_'+slat+'.ps -rotate -90 timeseries_mls_iwc_'+slat+'.png'
    endif
end
