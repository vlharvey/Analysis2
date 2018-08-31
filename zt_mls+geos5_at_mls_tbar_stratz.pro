;
; superimpose height of the stratopause
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
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
setplot='ps'
read,'setplot=',setplot
nxdim=750
nydim=750
xorig=[0.1,0.4,0.7]
yorig=[0.4,0.4,0.4]
xlen=0.25
ylen=0.25
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

restore,'stratopause_height_for_SABER.sav
restore,'stratopause_height_for_MLS.sav
restore,'stratopause_height_for_GEOS5.sav

restore,'zt_mls_tbar.sav
altitude_mls=altitude
TBAR_ZT_saber=TBAR_ZT
SDATE_ALL_saber=SDATE_ALL
kday=n_elements(sdate_all)
;
; interpolate small gaps in time
;
index=where(tbar_zt_saber eq 0.)
if index(0) ne -1L then begin
   tbar_zt_saber(index)=0./0.
endif
tbar_zt_saber=smooth(tbar_zt_saber,3,/NaN)
index=where(finite(tbar_zt_saber) ne 1L)
if index(0) ne -1L then tbar_zt_saber(index)=-999.
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

gheight=fltarr(n_elements(sdate_all))
gheightmax=fltarr(n_elements(sdate_all))
gheightmin=fltarr(n_elements(sdate_all))
mheight=fltarr(n_elements(sdate_all))
mheightmax=fltarr(n_elements(sdate_all))
mheightmin=fltarr(n_elements(sdate_all))

xindex=where(long(DATES) ge long(SDATE_ALL(0)) and long(dates) le max(long(sdate_all)))
MLSstratheight=reform(MLSstratheight(xindex,*,*))
geos5stratheight=reform(geos5stratheight(xindex,*,*))
dates=dates(xindex)
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
       device,/landscape,bits=8,filename='zt_mls+geos5_at_mls_tbar_'+slat+'.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
    endif
;
; interpolate stratopause height to latbin(j)
;
yindex=where(geos5lat gt latbin(j))
ilatg=yindex(0)
for iday=0L,n_elements(SDATE_ALL)-1L do begin
gheight(iday)=total((reform(GEOS5STRATHEIGHT(iday,*,ilatg-1))+reform(GEOS5STRATHEIGHT(iday,*,ilatg)))/2.)/float(n_elements(geos5lon))
gheightmax(iday)=max((reform(GEOS5STRATHEIGHT(iday,*,ilatg-1))+reform(GEOS5STRATHEIGHT(iday,*,ilatg)))/2.)
gheightmin(iday)=min((reform(GEOS5STRATHEIGHT(iday,*,ilatg-1))+reform(GEOS5STRATHEIGHT(iday,*,ilatg)))/2.)

for jj=1L,n_elements(mlslat)-1L do begin
    if mlslat(jj-1) lt latbin(j) and mlslat(jj) ge latbin(j) then begin
       yscale=(mlslat(jj)-latbin(j))/(mlslat(jj)-mlslat(jj-1))
       mh0=reform(MLSSTRATHEIGHT(iday,*,jj-1))
       mh1=reform(MLSSTRATHEIGHT(iday,*,jj))
       mheight(iday)=total(mh1-yscale*(mh1-mh0))/float(n_elements(mlslon))
       mheightmax(iday)=max(mh1-yscale*(mh1-mh0))
       mheightmin(iday)=min(mh1-yscale*(mh1-mh0))
    endif
endfor
endfor
;
; interpolate across data voids in MLS
;
for iday=1L,n_elements(SDATE_ALL)-2L do begin
    if mheight(iday) lt 30. then begin
       mheight(iday)=(mheight(iday-1)+mheight(iday+1))/2.
    endif
endfor
;
; plot zonal mean temperature
; MLS
;
erase
xyouts,.25,.7,'Zonal Mean Temperature at '+slat+' Latitude',color=0,charsize=2,/normal
xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7
nlvls=19
col1=1+indgen(nlvls)*icolmax/nlvls
level=180.+5.*findgen(nlvls)
tzlat=transpose(reform(tbar_zt_saber(j,*,*)))
print,slat,' ',min(tbar_zt_saber),max(tbar_zt_saber)
;
; hardwire xticks
;
xindex=[1.,31.+1.,31.+29.+1.]	; 15th of each month as record goes Jan 15 to Mar 15
xlabs=['J','F','M']
nxticks=n_elements(xlabs)
index=where(tzlat le 0.)
if index(0) ne -1L then tzlat(index)=-99.
contour,tzlat,1.+findgen(kday),altitude_mls,/noeras,xrange=[1.,kday],yrange=[20.,60.],$
      charsize=1.5,color=0,title='MLS',/fill,c_color=col1,$
      levels=level,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex,min_value=-99.
index=where(level mod 10 eq 0.)
contour,tzlat,1.+findgen(kday),altitude_mls,levels=level(index),color=0,/follow,/overplot,c_labels=0*index,min_value=-99.
loadct,0
oplot,1.+findgen(kday),mheight,color=150,thick=5
index=where(mheightmax gt 60.)
if index(0) ne -1L then mheightmax(index)=60.
for i=0L,kday-1 do begin
    if mheightmin(i) gt 30. then begin
    plots,i+1,mheightmin(i)
    plots,i+1,mheightmax(i),/continue,color=150,thick=2
    endif
endfor
loadct,39
imin=min(level)
imax=max(level)
ymnb=yorig(1) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(K)'
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for jj=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(jj)
x1=x1+dx
endfor
;
; GEOS-5
;
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7
tzlatg=transpose(reform(tbar_zt_geos5(j,*,*)))
syear=strmid(sdate_all_geos5,0,4)
smon=strmid(sdate_all_geos5,4,2)
sday=strmid(sdate_all_geos5,6,2)
xindex=where(sday eq '15',nxticks)
index=where(tzlatg le 0.)
if index(0) ne -1L then tzlatg(index)=-99.
contour,tzlatg,1.+findgen(kday),altitude,/noeras,xrange=[1.,kday],yrange=[20.,60.],$
      charsize=1.5,color=0,title='GEOS-5',/fill,c_color=col1,ytitle='Altitude (km)',$
      levels=level,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex,min_value=-99.
index=where(level mod 10 eq 0.)
contour,tzlatg,1.+findgen(kday),altitude,levels=level(index),color=0,/follow,/overplot,c_labels=0*index,min_value=-99.
oplot,1.+findgen(kday),gheight,color=0,thick=5	;,psym=8
index=where(gheightmax gt 60.)
if index(0) ne -1L then gheightmax(index)=60.
for i=0L,kday-1 do begin
    if gheightmin(i) gt 20. then begin
    plots,i+1,gheightmin(i)
    plots,i+1,gheightmax(i),/continue,color=0
    endif
endfor
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(K)'
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for jj=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(jj)
x1=x1+dx
endfor

; superimpose pdiff
;
restore,'c11_rb.tbl'
tvlct,c1,c2,c3
col2=1+indgen(11)
pdiff=-99.+0.*tzlatg
index=where(tzlatg gt 0. and tzlat gt 0.)
if index(0) eq -1L then goto,jump
pdiff(index)=tzlatg(index)-tzlat(index)
level=-20.+4.*findgen(11)
!type=2^2+2^3
xmn=xorig(2)
xmx=xorig(2)+xlen
ymn=yorig(2)
ymx=yorig(2)+ylen
set_viewport,xmn,xmx,ymn,ymx
contour,pdiff,1.+findgen(kday),altitude,xrange=[1.,kday],yrange=[20.,60.],$
        charsize=1.5,levels=level,/cell_fill,$
        title='GEOS-5 - MLS',c_color=col2,color=0,min_value=-99.,$
        xticks=nxticks-1,xtickname=xlabs,xtickv=xindex
contour,pdiff,1.+findgen(kday),altitude,/overplot,levels=level,color=0,/follow,min_value=-99.,$
        c_labels=0*level
contour,pdiff,1.+findgen(kday),altitude,/overplot,levels=[0],color=0,thick=1,min_value=-99.
;contour,tzlat,1.+findgen(kday),altitude,/overplot,levels=[230.],color=0,thick=6
oplot,1.+findgen(kday),gheight,color=0,thick=5
index=where(gheightmax gt 60.)
if index(0) ne -1L then gheightmax(index)=60.
for i=0L,kday-1 do begin
    if gheightmin(i) gt 20. then begin
    plots,i+1,gheightmin(i)
    plots,i+1,gheightmax(i),/continue,color=0
    endif
endfor
;loadct,0
;oplot,1.+findgen(kday),mheight,color=150,thick=5
;restore,'c11_rb.tbl'
;tvlct,c1,c2,c3
imin=min(level)
imax=max(level)
ymnb=ymn -cbaryoff
ymxb=ymnb+cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,$
      xtitle='(K)',color=0,xticks=n_elements(level)/2
ybox=[0,10,10,0,0]
x2=imin
dx=(imax-imin)/(float(n_elements(col2)))
for jj=0L,n_elements(col2)-1 do begin
    xbox=[x2,x2,x2+dx,x2+dx,x2]
    polyfill,xbox,ybox,color=col2(jj)
    x2=x2+dx
endfor
loadct,39
jump:
    if max(tzlat) gt 0. and setplot ne 'ps' then stop

    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert -trim zt_mls+geos5_at_mls_tbar_'+slat+'.ps -rotate -90 zt_mls+geos5_at_mls_tbar_'+slat+'.jpg'
    endif
endfor	; loop over latitude
end
