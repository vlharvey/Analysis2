;
; latitude time plot of monthly mean zonal means at each altitude
; monthly mean zonal mean WACCM "specified met" version and MLS temperature
;
@stddat
@kgmt
@ckday
@kdate

loadct,39
device,decompose=0
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
icmm1=icolmax-1
icmm2=icolmax-2
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
!noeras=1
nxdim=750
nydim=750
xorig=[0.1,0.4,0.7]
yorig=[0.4,0.4,0.4]
xlen=0.25
ylen=0.25
cbaryoff=0.075
cbarydel=0.01
setplot='x'
read,'setplot=',setplot
if setplot ne 'ps' then begin
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=icolmax
endif
mdir='/aura2/harvey/MLS_data/Datfiles_SOSST/'
mdir2='/aura2/randall/mls_data/'
month=['January','February','March','April','May','June','July','August','September','October','November','December']
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
mday=[31,28,31,30,31,30,31,31,30,31,30,31]
;
; restore monthly mean WACCM data
;
restore,'/aura2/harvey/WACCM_data/Datfiles/wcm_geos.cam2.h0.2007.sav
date=[$
20070101L,$
20070201L,$
20070301L,$
20070401L,$
20070501L,$
20070601L,$
20070701L,$
20070801L,$
20070901L,$
20071001L,$
20071101L,$
20071201L,$
20080101L]

wdate=date
nmonth=n_elements(date)
sdate=strcompress(date,/remove_all)
smon=strcompress(long(strmid(sdate,4,2)),/remove_all)
for imonth=0L,nmonth-1L do begin
    syear=strmid(sdate(imonth),0,4)
    smonth=strmid(sdate(imonth),4,2)
    print,syear+smonth
;
; restore monthly mean save file
;
;save,file=mdir+'waccm_mls_v2.2_'+syear+smonth+'_temp.sav',wlat,altitude,tbar,tbarz
    restore,mdir+'waccm_mls_v2.2_'+syear+smonth+'_temp.sav'
;
; 3-d arrays
;   
    if imonth eq 0L then begin
       nz=n_elements(altitude)
       nr=n_elements(wlat)
       mls3d=fltarr(nr,nz,nmonth)
       waccm3d=fltarr(nr,nz,nmonth)
    endif
;
; save each month's zonal means into 3-d array
;
    mls3d(*,*,imonth)=tbar
    waccm3d(*,*,imonth)=tbarz
endfor  ; loop over months

for ilev=50L,nz-1L do begin
    salt=strcompress(long(altitude(ilev)),/remove_all)
    mls2d=transpose(reform(mls3d(*,ilev,*)))
    waccm2d=transpose(reform(waccm3d(*,ilev,*)))

if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='tz_waccm_mls_temp_'+salt+'.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
endif
;
; plot temperature
;
erase
xyouts,.3,.7,'Temperature at '+salt+' km',charsize=2,charthick=1.5,/normal,color=0
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
nlvls=25
col1=1+indgen(nlvls)*icolmax/nlvls
level=110.+5.*findgen(nlvls)
index=where(waccm2d eq 0.)
if index(0) ne -1L then waccm2d(index)=-9999.
contour,waccm2d,1.+findgen(nmonth),wlat,/noeras,xrange=[1,nmonth],yrange=[-90.,90.],charsize=1.25,color=0,$
      ytitle='Latitude',xticks=nmonth-1,/fill,c_color=col1,levels=level,min_value=-9999.,$
      xtitle='Month',title='WACCM',charthick=1.5,xtickname=smon
contour,waccm2d,1.+findgen(nmonth),wlat,levels=level,color=0,/follow,/overplot,min_value=-9999.
contour,waccm2d,1.+findgen(nmonth),wlat,levels=[150.],color=mcolor,thick=5,/follow,/overplot,min_value=-9999.
contour,waccm2d,1.+findgen(nmonth),wlat,levels=[280.],color=0,thick=5,/follow,/overplot,min_value=-9999.
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(K)',charsize=1.25,charthick=1.5
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor
;
; interpolate MLS through October
;
for kk=0L,nr-1 do begin
    for ii=1L,nmonth-2L do begin
        if mls2d(ii,kk) lt 100. then begin
           mls2d(ii,kk)=(mls2d(ii-1,kk)+mls2d(ii+1,kk))/2.
        endif
    endfor 
endfor
xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
index=where(mls2d eq 0.)
if index(0) ne -1L then mls2d(index)=-9999.
contour,mls2d,1.+findgen(nmonth),wlat,/noeras,xrange=[1,nmonth],yrange=[-90.,90.],charsize=1.25,color=0,$
      xticks=nmonth-1,/fill,c_color=col1,levels=level,min_value=-9999.,$
      xtitle='Month',title='MLS',charthick=1.5,xtickname=smon
contour,mls2d,1.+findgen(nmonth),wlat,levels=level,color=0,/follow,/overplot,min_value=-9999.
contour,mls2d,1.+findgen(nmonth),wlat,levels=[150.],color=mcolor,thick=5,/follow,/overplot,min_value=-9999.
contour,mls2d,1.+findgen(nmonth),wlat,levels=[280.],color=0,thick=5,/follow,/overplot,min_value=-9999.
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(K)',charsize=1.25,charthick=1.5
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor

restore,'c11_rb.tbl'
tvlct,c1,c2,c3
col2=1+indgen(11)
pdiff=-99.+0.*waccm2d
index=where(waccm2d gt 0. and mls2d gt 0.)
if index(0) eq -1L then goto,jump
pdiff(index)=waccm2d(index)-mls2d(index)
level=-50.+10.*findgen(11)
!type=2^2+2^3
xmn=xorig(2)
xmx=xorig(2)+xlen
ymn=yorig(2)
ymx=yorig(2)+ylen
set_viewport,xmn,xmx,ymn,ymx
contour,pdiff,1.+findgen(nmonth),wlat,xrange=[1,nmonth],yrange=[-90.,90.],$
        xticks=nmonth-1,xtitle='Month',charsize=1.25,levels=level,/cell_fill,$
        title='WACCM - MLS',c_color=col2,color=0,min_value=-99.,charthick=1.5,xtickname=smon
index=where(level gt 0.)
contour,pdiff,1.+findgen(nmonth),wlat,/overplot,levels=level(index),color=0,/follow,min_value=-99.
index=where(level lt 0.)
contour,pdiff,1.+findgen(nmonth),wlat,/overplot,levels=level(index),color=mcolor,/follow,min_value=-99.
contour,pdiff,1.+findgen(nmonth),wlat,/overplot,levels=[0],color=0,thick=3,min_value=-99.
imin=min(level)
imax=max(level)
ymnb=ymn -cbaryoff
ymxb=ymnb+cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,$
      xtitle='(K)',color=0,xticks=n_elements(level)/2,charsize=1.25,charthick=1.5
ybox=[0,10,10,0,0]
x2=imin
dx=(imax-imin)/(float(n_elements(col2)))
for jj=0L,n_elements(col2)-1 do begin
    xbox=[x2,x2,x2+dx,x2+dx,x2]
    polyfill,xbox,ybox,color=col2(jj)
    x2=x2+dx
endfor

if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim tz_waccm_mls_temp_'+salt+'.ps -rotate -90 tz_waccm_mls_temp_'+salt+'.jpg'
endif
jump:
loadct,39

endfor	; loop over altitude
end
