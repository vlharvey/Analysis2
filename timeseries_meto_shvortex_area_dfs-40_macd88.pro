;
; timeseries of MetO SH vortex area at 900 K and 70 S on DFS -40 and NH PMC onset date
;
device,decompose=0
!p.background=255
setplot='ps'
read,'setplot=',setplot
nxdim=750
nydim=750
xorig=[0.2]
yorig=[0.2]
xlen=0.6
ylen=0.5
cbaryoff=0.12
cbarydel=0.01
!NOERAS=-1
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
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
ktheta=900.
rtheta=float(ktheta)
;read,'Enter theta ',rtheta
ith=where(th eq rtheta)
stheta=strcompress(long(rtheta),/remove_all)
sYYYYMMDD_ALL=strcompress(YYYYMMDD_ALL,/remove_all)
syear=strmid(syyyymmdd_all,0,4)
smon=strmid(syyyymmdd_all,4,2)
;for ilat=-80,-60,5 do begin
ilat=-70
rlat=float(ilat)
slat=strcompress(long(rlat),/remove_all)

index=where(abs(alat-float(ilat)) eq min(abs(alat-float(ilat))))
print,alat(index)
ztvortex_area=fltarr(nday)
for n=0L,nday-1L do begin
    ztvortex_area(n)=mean(VMARKBAR_ALL(n,index,ith(0)))
endfor
area1d=fltarr(7)

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

index=where(dfs ge -45. and dfs le -35.)
ztvortex=mean(ztvortex(index))
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
          /bold,/color,bits_per_pixel=8,/helvetica,filename='timeseries_meto_shvortex_area_dfs-40_'+slat+'_'+stheta+'K.ps'
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
plot,2007+findgen(7),findgen(7),/nodata,color=0,charsize=1.5,charthick=2,xrange=[2007,2013],$
     ytitle='MetO SH Vortex Area on DFS -40',yrange=[60.,110.],xticks=6,title=slat+' Latitude '+stheta+' K'
endif
oplot,[lstyr,lstyr],[ztvortex,ztvortex],color=0,psym=8
area1d(lstyr-2007)=ztvortex
endfor
oplot,2007+findgen(7),area1d,color=0,thick=15
NHyear=[2007,2008,2009,2010,2011,2012,2013]
NHonset=[-27.,-20.,-24.,-24.,-26.,-27.,-34.]
axis,yaxis=1,yrange=[-35.,-20.],/save,color=11,ytitle='NH PMC Onset Date (DFS)',charsize=2,charthick=2
oplot,nhyear,nhonset,color=11,thick=15
result=correlate(area1d,nhonset)
print,result
xyouts,xmn+0.2,ymn+0.05,'r = '+string(format='(f5.2)',result(0)),/normal,charsize=2,charthick=2,color=0

    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert -trim timeseries_meto_shvortex_area_dfs-40_'+slat+'_'+stheta+'K.ps -rotate -90 timeseries_meto_shvortex_area_dfs-40_'+slat+'_'+stheta+'K.png'
    endif
;endfor
end
