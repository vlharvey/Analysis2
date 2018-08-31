;
; Plot time-altitude Tbar 
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
xlen=0.6
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
; loop over latitude
;
;for j=0L,nlat-1L do begin
for j=33,33L do begin
    slat=strcompress(long(latbin(j)),/remove_all)
    if setplot eq 'ps' then begin
       lc=0
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       !p.thick=2.0                   ;Plotted lines twice as thick
       device,font_size=9
       device,/landscape,bits=8,filename='timeseries_mls+geos5_at_mls_tbar+diffs_'+slat+'.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
    endif
;
; temperature timeseries at 35 km
;
    ilev=40L
    slev=strcompress(long(altitude(ilev)),/remove_all)
    erase
    xmn=xorig(0)
    xmx=xorig(0)+xlen
    ymn=yorig(0)
    ymx=yorig(0)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    !type=2^2+2^3+2^7
    nlvls=19
    col1=1+indgen(nlvls)*icolmax/nlvls
    level=180.+5.*findgen(nlvls)
    tzlat=reform(tbar_zt_mls(j,ilev,*))
;
; hardwire xticks
;
    xindex=[1.,31.+1.,31.+29.+1.]	; 15th of each month as record goes Jan 15 to Mar 15
    xlabs=['J','F','M']
    nxticks=n_elements(xlabs)
    index=where(tzlat le 0.)
    if index(0) ne -1L then tzlat(index)=-99.
    plot,1.+findgen(kday),tzlat,/noeras,xrange=[1.,kday],yrange=[220.,290.],$
          charsize=1.5,color=0,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex,$
          min_value=-99.,title='Zonal Mean Temperature at '+slat+' N and '+slev+' km',$
          ytitle='Temperature (K)'

    tzlatg=reform(tbar_zt_geos5(j,ilev,*))
    oplot,1.+findgen(kday),tzlatg,color=mcolor*.9

    pdiff=-99.+0.*tzlatg
    index=where(tzlatg gt 0. and tzlat gt 0.)
    if index(0) eq -1L then goto,jump
    pdiff(index)=tzlatg(index)-tzlat(index)
    axis,xrange=[1.,kday],yrange=[-15,5],/save,yaxis=1,/data,charsize=1.5,$
         yticklen=-0.02,yticks=8,color=0,ytitle='GEOS-5 - MLS'
    oplot,1.+findgen(kday),pdiff,color=0,thick=3
    plots,1.,0.
    plots,kday,0.,/continue,color=0

    jump:
    if max(tzlat) gt 0. and setplot ne 'ps' then stop

    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert -trim timeseries_mls+geos5_at_mls_tbar+diffs_'+slat+'.ps -rotate -90 timeseries_mls+geos5_at_mls_tbar+diffs_'+slat+'.jpg'
    endif
endfor	; loop over latitude
end
