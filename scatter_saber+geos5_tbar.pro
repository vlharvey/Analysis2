;
; Plot time-altitude Tbar plus scatter plots of Tbar and Tdiff
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
xorig=[0.1,0.4,0.7,0.1,0.4,0.7]
yorig=[0.5,0.5,0.5,0.1,0.1,0.1]
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
;
; ALAT            FLOAT     = Array[35]
; ALTITUDE        FLOAT     = Array[121]
; GLAT            FLOAT     = Array[72]
; SDATE_ALL       STRING    = Array[65]
; TBARZ_YT_GEOS5  FLOAT     = Array[72, 121, 65]
; TBARZ_YT_SABER  FLOAT     = Array[35, 121, 65]
; UBARZ_YT_GEOS5  FLOAT     = Array[72, 121, 65]
; UBARZ_YT_SABER  FLOAT     = Array[35, 121, 65]       
;
restore,'yt_saber+geos5_ztu_20080101-20080305.sav'
nlat=n_elements(alat)
kday=n_elements(sdate_all)
;
; loop over latitude
;
for j=28L,nlat-1L do begin
    slat=strcompress(string(format='(f5.1)',alat(j)),/remove_all)
    if setplot eq 'ps' then begin
       lc=0
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       !p.thick=2.0                   ;Plotted lines twice as thick
       device,font_size=9
       device,/landscape,bits=8,filename='scatter_saber+geos5_tbar_'+slat+'.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
    endif
;
; plot zonal mean temperature
;
erase
xyouts,.25,.8,'Zonal Mean Temperature at '+slat+' Latitude',color=0,charsize=2,/normal
xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7
nlvls=19
col1=1+indgen(nlvls)*icolmax/nlvls
level=180.+5.*findgen(nlvls)
tzlat=transpose(reform(TBARZ_YT_SABER(j,*,*)))
;
; hardwire xticks
;
xindex=[1.,31.+1.,31.+29.+1.]	; 15th of each month as record goes Jan 15 to Mar 15
xlabs=['J','F','M']
nxticks=n_elements(xlabs)
index=where(tzlat le 0.)
if index(0) ne -1L then tzlat(index)=-99.
contour,tzlat,1.+findgen(kday),altitude,/noeras,xrange=[1.,kday],yrange=[20.,68.],$
      charsize=1.5,color=0,ytitle='Altitude (km)',title='SABER',/fill,c_color=col1,$
      levels=level,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex,min_value=-99.
index=where(level mod 10 eq 0.)
contour,tzlat,1.+findgen(kday),altitude,levels=level(index),color=0,/follow,/overplot,c_labels=0*index,min_value=-99.
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
syear=strmid(sdate_all,0,4)
smon=strmid(sdate_all,4,2)
sday=strmid(sdate_all,6,2)
xindex=where(sday eq '15',nxticks)
index=where(abs(glat-alat(j)) eq min(abs(glat-alat(j))))
tzlatg=transpose(reform(TBARZ_YT_GEOS5(index(0),*,*)))
index=where(tzlatg le 0.)
if index(0) ne -1L then tzlatg(index)=-99.
contour,tzlatg,1.+findgen(kday),altitude,/noeras,xrange=[1.,kday],yrange=[20.,68.],$
      charsize=1.5,color=0,title='GEOS-5',/fill,c_color=col1,$
      levels=level,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex,min_value=-99.
index=where(level mod 10 eq 0.)
contour,tzlatg,1.+findgen(kday),altitude,levels=level(index),color=0,/follow,/overplot,c_labels=0*index,min_value=-99.
;
; superimpose pdiff
;
pdiff=-99.+0.*tzlatg
index=where(tzlatg gt 0. and tzlat gt 0.)
if index(0) eq -1L then goto,jump
pdiff(index)=tzlatg(index)-tzlat(index)
;contour,smooth(pdiff,3),1.+findgen(kday),altitude,/overplot,levels=[0],color=0,thick=1,min_value=-99.
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

restore,'c11_rb.tbl'
tvlct,c1,c2,c3
col2=1+indgen(11)
pdiff=-99.+0.*tzlatg
index=where(tzlatg gt 0. and tzlat gt 0.)
if index(0) eq -1L then goto,jump
pdiff(index)=tzlatg(index)-tzlat(index)
level=-25.+5.*findgen(11)
!type=2^2+2^3
xmn=xorig(2)
xmx=xorig(2)+xlen
ymn=yorig(2)
ymx=yorig(2)+ylen
set_viewport,xmn,xmx,ymn,ymx
contour,pdiff,1.+findgen(kday),altitude,xrange=[1.,kday],yrange=[20.,68.],$
        charsize=1.5,levels=level,/cell_fill,$
        title='GEOS-5 - SABER',c_color=col2,color=0,min_value=-99.,$
        xticks=nxticks-1,xtickname=xlabs,xtickv=xindex
contour,pdiff,1.+findgen(kday),altitude,/overplot,levels=level,color=0,/follow,min_value=-99.,$
        c_labels=0*level
contour,pdiff,1.+findgen(kday),altitude,/overplot,levels=[0],color=0,thick=2,min_value=-99.
contour,tzlat,1.+findgen(kday),altitude,levels=[250.],color=0,/follow,/overplot,thick=3
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
;
; scatter plot of T vs Tdiff
;
tzlatg1=reform(tzlatg(*,20:35))
tzlat1=reform(tzlat(*,20:35))
pdiff=0.*tzlat1
index=where(tzlatg1 gt 0. and tzlat1 gt 0.)
if index(0) eq -1L then goto,jump
pdiff(index)=tzlatg1(index)-tzlat1(index)
!type=2^2+2^3
xmn=xorig(3)
xmx=xorig(3)+xlen
ymn=yorig(3)
ymx=yorig(3)+ylen
set_viewport,xmn,xmx,ymn,ymx
plot,tzlat1(index),pdiff(index),xrange=[200.,270.],yrange=[-15.,5.],psym=8,$
     charsize=1.5,xtitle='SABER Temperature',color=0,ytitle='T Difference',$
     title='20-35 km'
;oplot,tzlatg1(index),pdiff(index),psym=8,color=90
for iday=0L,kday-1 do begin
    tprof=reform(tzlat1(iday,*))
    tdiffprof=reform(pdiff(iday,*))
    index=where(tprof gt 0.)
    if index(0) ne -1L then oplot,tprof(index),tdiffprof(index),psym=8,color=((iday+1)/(kday+2.))*mcolor
endfor
plots,200,0
plots,270,0,/continue,color=0

tzlatg1=reform(tzlatg(*,35:45))
tzlat1=reform(tzlat(*,35:45))
pdiff=0.*tzlat1
index=where(tzlatg1 gt 0. and tzlat1 gt 0.)
if index(0) eq -1L then goto,jump
pdiff(index)=tzlatg1(index)-tzlat1(index)
!type=2^2+2^3
xmn=xorig(4)
xmx=xorig(4)+xlen
ymn=yorig(4)
ymx=yorig(4)+ylen
set_viewport,xmn,xmx,ymn,ymx
plot,tzlat1(index),pdiff(index),xrange=[200.,270.],yrange=[-15.,5.],psym=8,$
     charsize=1.5,xtitle='SABER Temperature',color=0,$
     title='35-45 km'
;oplot,tzlatg1(index),pdiff(index),psym=8,color=90
for iday=0L,kday-1 do begin
    tprof=reform(tzlat1(iday,*))
    tdiffprof=reform(pdiff(iday,*))
    index=where(tprof gt 0.)
    if index(0) ne -1L then oplot,tprof(index),tdiffprof(index),psym=8,color=((iday+1)/(kday+2.))*mcolor
endfor
plots,200,0
plots,270,0,/continue,color=0

tzlatg1=reform(tzlatg(*,55:65))
tzlat1=reform(tzlat(*,55:65))
pdiff=0.*tzlat1
index=where(tzlatg1 gt 0. and tzlat1 gt 0.)
if index(0) eq -1L then goto,jump
pdiff(index)=tzlatg1(index)-tzlat1(index)
!type=2^2+2^3
xmn=xorig(5)
xmx=xorig(5)+xlen
ymn=yorig(5)
ymx=yorig(5)+ylen
set_viewport,xmn,xmx,ymn,ymx
plot,tzlat1(index),pdiff(index),xrange=[200.,270.],yrange=[-15.,25.],psym=8,$
     charsize=1.5,xtitle='SABER Temperature',color=0,$
     title='55-65 km'
;oplot,tzlatg1(index),pdiff(index),psym=8,color=90
for iday=0L,kday-1 do begin
    tprof=reform(tzlat1(iday,*))
    tdiffprof=reform(pdiff(iday,*))
    index=where(tprof gt 0.)
    if index(0) ne -1L then oplot,tprof(index),tdiffprof(index),psym=8,color=((iday+1)/(kday+2.))*mcolor
endfor
plots,200,0
plots,270,0,/continue,color=0

jump:
    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert -trim scatter_saber+geos5_tbar_'+slat+'.ps -rotate -90 scatter_saber+geos_tbar_'+slat+'.jpg'
    endif
endfor	; loop over latitude
end
