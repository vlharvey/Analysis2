;
; Plot time series of GEOS-5 and MLS Temperature
;
@stddat
@kgmt
@ckday
@kdate

loadct,39
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
device,decompose=0
!p.background=icolmax
setplot='ps'
read,'setplot=',setplot
nxdim=750
nydim=750
xorig=[0.2]
yorig=[0.25]
xlen=0.7
ylen=0.5
cbaryoff=0.05
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smonth=['J','F','M','A','M','J','J','A','S','O','N','D']

restore,'zt_saber_tbar.sav
altitude_saber=altitude
TBAR_ZT_saber=TBAR_ZT
SDATE_ALL_saber=SDATE_ALL
kday=n_elements(sdate_all)
index=where(tbar_zt_saber eq 0.)
if index(0) ne -1L then begin
   tbar_zt_saber(index)=0./0.
endif
tbar_zt_saber=smooth(tbar_zt_saber,3,/NaN)
index=where(finite(tbar_zt_saber) ne 1L)
if index(0) ne -1L then tbar_zt_saber(index)=-999.

restore,'zt_mls_tbar.sav
altitude_mls=altitude
TBAR_ZT_mls=TBAR_ZT
SDATE_ALL_mls=SDATE_ALL
kday=n_elements(sdate_all)
;
; interpolate small gaps in time
;
index=where(tbar_zt_mls eq 0.)
if index(0) ne -1L then begin
   tbar_zt_mls(index)=0./0.
endif
tbar_zt_mls=smooth(tbar_zt_mls,3,/NaN)
index=where(finite(tbar_zt_mls) ne 1L)
if index(0) ne -1L then tbar_zt_mls(index)=-999.
nlat=n_elements(latbin)
restore,'zt_geos5_at_mls_tbar.sav
TBAR_ZT_geos5=TBAR_ZT
SDATE_ALL_geos5=SDATE_ALL
index=where(tbar_zt_geos5 eq 0.)
if index(0) ne -1L then begin
   tbar_zt_geos5(index)=0./0.
endif
tbar_zt_geos5=smooth(tbar_zt_geos5,3,/NaN)
index=where(finite(tbar_zt_geos5) ne 1L)
if index(0) ne -1L then tbar_zt_geos5(index)=-999.
;
; choose altitude
;
ralt=0.
print,altitude
read,'Enter altitude ',ralt
index=where(altitude eq ralt)
ialt=index(0)
salt=strcompress(long(ralt),/remove_all)
;
; loop over latitude
;
;for j=0L,nlat-1L do begin
for j=29L,34L do begin
    slat=strcompress(string(format='(f5.1)',latbin(j)),/remove_all)
    if setplot eq 'ps' then begin
       lc=0
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       !p.thick=2.0                   ;Plotted lines twice as thick
       device,font_size=9
       device,/landscape,bits=8,filename='timeseries_mls+geos5_at_mls_tbar_'+slat+'_'+salt+'km.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
    endif
;
; plot zonal mean temperature
; MLS
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7
tzlat=transpose(reform(tbar_zt_mls(j,ialt,*)))
tzlats=transpose(reform(tbar_zt_saber(j,ialt,*)))
tzlatg=transpose(reform(tbar_zt_geos5(j,ialt,*)))
;
; hardwire xticks
;
xindex=[1.,31.+1.,31.+29.+1.]	; 15th of each month as record goes Jan 15 to Mar 15
xlabs=['J','F','M']
nxticks=n_elements(xlabs)
index=where(tzlat le 0.)
if index(0) ne -1L then tzlat(index)=-99.
index=where(tzlatg le 0.)
if index(0) ne -1L then tzlatg(index)=-99.
index=where(tzlat gt 0.)
tmin=min(tzlat(index))-5.
tmax=max(tzlat(index))+5.
plot,1.+findgen(kday),tzlatg,/noeras,xrange=[1.,kday],yrange=[tmin,tmax],xtitle='Month',$
     title=salt+' km Zonal Mean Temperature at '+slat+' Latitude',ytitle='Temperature',$
     charsize=1.5,color=0,thick=7,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex,min_value=-99.
oplot,1.+findgen(kday),tzlat,thick=7,color=100
;oplot,1.+findgen(kday),tzlats,thick=7,color=200
xyouts,xmx-0.1,ymx-0.02,'GEOS-5',color=0,/normal,charsize=1.5,charthick=2
xyouts,xmx-0.1,ymx-0.05,'MLS',color=100,/normal,charsize=1.5,charthick=2
;xyouts,xmx-0.1,ymx-0.08,'SABER',color=200,/normal,charsize=1.5,charthick=2

    if max(tzlat) gt 0. and setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert -trim timeseries_mls+geos5_at_mls_tbar_'+slat+'_'+salt+'km.ps '+$
             ' -rotate -90 timeseries_mls+geos5_at_mls_tbar_'+slat+'_'+salt+'km.jpg'
    endif
endfor	; loop over latitude
end
