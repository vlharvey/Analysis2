;
; timeseries of MetO vortex area at different latitudes
; NH version (SH vortex area)
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
restore,'MetO_antic_area_alllat_1991-2013_nc4.sav
;print,th
for ktheta=400,1200,100 do begin
rtheta=float(ktheta)
;read,'Enter theta ',rtheta
ith=where(th eq rtheta)
stheta=strcompress(long(rtheta),/remove_all)
sYYYYMMDD_ALL=strcompress(YYYYMMDD_ALL,/remove_all)
syear=strmid(syyyymmdd_all,0,4)
smon=strmid(syyyymmdd_all,4,2)
for ilat=-80,-60,5 do begin
rlat=float(ilat)
slat=strcompress(long(rlat),/remove_all)

index=where(abs(alat-float(ilat)) eq min(abs(alat-float(ilat))))
print,alat(index)
ztvortex_area=fltarr(nday)
for n=0L,nday-1L do begin
    ztvortex_area(n)=mean(VMARKBAR_ALL(n,index,ith(0)))
endfor
for lstyr=2007,2013 do begin
yearlab=strcompress(lstyr,/remove_all)
;
; extract MJJA for each year
;
index=where(long(syear) eq lstyr and (smon eq '05' or smon eq '06' or smon eq '07' or smon eq '08'),kday)
ztvortex=reform(ztvortex_area(index))
;
; convert sdate to DFS
;
doy=121+findgen(kday)
jdaysol=172.
dfs=doy-jdaysol
;
; zonal mean
;
if lstyr eq 2007 then begin
if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !p.font=0
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='timeseries_meto_shvortex_area_mjja_'+slat+'_'+stheta+'K.ps'
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
plot,dfs,ztvortex,/nodata,color=0,charsize=2,charthick=2,xrange=[-50.,80.],$
     ytitle='MetO SH Vortex Area (% of Lat circle)',yrange=[0.,100.],xtitle='DFS in NH'
xyouts,xmn+0.2,ymn+0.075,slat+' Latitude '+stheta+' K',/normal,charsize=2,charthick=2,color=0
endif
print,min(ztvortex),max(ztvortex)
bad=where(ztvortex le 0.1)
ztvortex(bad)=0./0.
ztvortex=smooth(ztvortex,7,/edge_truncate,/NaN)
if lstyr eq 2013 then loadct,0 & col1(6)=150
oplot,dfs,ztvortex,color=col1(lstyr-2007),thick=10
xyouts,xmn+0.01+0.1*(lstyr-2007),ymn+0.01,yearlab,/normal,charsize=1.75,charthick=2,color=col1(lstyr-2007)
endfor

    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert -trim timeseries_meto_shvortex_area_mjja_'+slat+'_'+stheta+'K.ps -rotate -90 timeseries_meto_shvortex_area_mjja_'+slat+'_'+stheta+'K.png'
    endif
endfor
endfor	; loop over theta
end
