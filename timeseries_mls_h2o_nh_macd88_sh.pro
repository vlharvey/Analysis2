;
; timeseries of MLS Water Vapor
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
;col1=[0,1,2,3,4,5,10]
col1=[0,1,3,4,6,11]
for ilat=-80,-60,5 do begin
rlat=float(ilat)
slat=strcompress(long(rlat),/remove_all)

for lstyr=2007,2012 do begin
yearlab=strcompress(lstyr,/remove_all)+'-'+strcompress(lstyr+1,/remove_all)

restore,'xt_mls_h2o_'+yearlab+'_'+slat+'.sav'		; mlspolarh2o_xt,kday,longrid,sdate_all
;
; convert sdate to DFS
;
doy=305+findgen(121)
jdaysol=355.
dfs=doy-jdaysol
;
; zonal mean
;
zmh2o=fltarr(kday)
for k=0L,kday-1L do begin
    good=where(MLSPOLARH2O_XT(*,k) gt 0.,nn)
    if good(0) ne -1L then zmh2o(k)=mean(MLSPOLARH2O_XT(good,k))
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
          /bold,/color,bits_per_pixel=8,/helvetica,filename='timeseries_mls_h2o_'+slat+'.ps'
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
plot,dfs,zmh2o,/nodata,color=0,charsize=2,charthick=2,xrange=[-50.,80.],$
     ytitle='Zonal Mean MLS Water Vapor (ppmv)',yrange=[0.,8.],xtitle='DFS'
xyouts,xmn+0.25,ymx-0.04,slat+' Latitude',/normal,charsize=2,charthick=2,color=0
endif
print,min(zmh2o),max(zmh2o)
bad=where(zmh2o le 0.1)
zmh2o(bad)=0./0.
zmh2o=smooth(zmh2o,5,/edge_truncate,/NaN)
oplot,dfs,zmh2o,color=col1(lstyr-2007),thick=10
yearlab=strcompress(lstyr,/remove_all)
xyouts,xmn+0.05+0.1*(lstyr-2007),ymn+0.01,yearlab,/normal,charsize=1.75,charthick=2,color=col1(lstyr-2007)
endfor

    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert -trim timeseries_mls_h2o_'+slat+'.ps -rotate -90 timeseries_mls_h2o_'+slat+'.png'
    endif
endfor
end
