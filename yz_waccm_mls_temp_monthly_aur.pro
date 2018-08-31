;
; multi-year monthly mean zonal mean WACCM "aurora run" and MLS temperature
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
;
; 20070201 is for January monthly average
;
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

nz=nz-1
wlat=latitude
dy=wlat(1)-wlat(0)
wlon=longitude
wdate=date
nmonth=n_elements(date)
sdate=strcompress(date,/remove_all)
for imonth=0L,nmonth-1L do begin
;
; WACCM zonal mean temperature and geopotential height
;
    t3d=reform(temp_4d(*,*,*,imonth))
    z3d=reform(z_4d(*,*,*,imonth))/1000.
    tbarg=fltarr(nr,nz)
    zbarg=fltarr(nr,nz)
    for k=0L,nz-1L do begin
    for j=0L,nr-1L do begin
        tpts=reform(t3d(*,j,k))
        zpts=reform(z3d(*,j,k))
        index=where(tpts ne 0.,nx)
        if index(0) ne -1L then begin
           tbarg(j,k)=total(tpts(index))/float(nx)
           zbarg(j,k)=total(zpts(index))/float(nx)
        endif
    endfor
    endfor
;
; loop over days in the month
;
    iyear=long(strmid(sdate(imonth),0,4))
    imon=long(strmid(sdate(imonth),4,2))
    syear=strmid(sdate(imonth),0,4)
    smonth=strmid(sdate(imonth),4,2)
;
; check for monthly mean save file
;
    dum=findfile(mdir+'waccm_mls_v2.2_'+syear+smonth+'_temp.sav')
;   if dum(0) ne '' then goto,plotit
    nday=mday(imon-1)
    if iyear mod 4 eq 0 and imon eq 1 then nday=29
    for iday=0L,nday-1L do begin
;
;***Read MLS data
;
      sday=string(FORMAT='(I2.2)',iday+1)
      yyyymmdd=syear+smonth+sday
print,yyyymmdd
;
; restore MLS tpd_mls_v2.2_20080614.sav file
;
; ALTITUDE        FLOAT     = Array[121]
; COMMENT         STRING    = Array[7]
; DENTOT          FLOAT     = Array[3494, 121]
; ID              STRING    = Array[3494]
; PRESSURE        FLOAT     = Array[3494, 121]
; TEMPERATURE     FLOAT     = Array[3494, 121]
; TEMPERATURE_ERROR
; TEMPERATURE_MASK
;
    dum=findfile(mdir+'tpd_mls_v2.2_'+yyyymmdd+'.sav')
    if dum(0) ne '' then begin
    restore,mdir+'tpd_mls_v2.2_'+yyyymmdd+'.sav'
    restore,mdir+'cat_mls_v2.2_'+yyyymmdd+'.sav'   ; latitude
    endif
    if dum(0) eq '' then begin
       dum=findfile(mdir2+'cat_mls_v2.2_'+yyyymmdd+'.sav')
       if dum(0) eq '' then goto,jump
       restore,mdir2+'tpd_mls_v2.2_'+yyyymmdd+'.sav'
       restore,mdir2+'cat_mls_v2.2_'+yyyymmdd+'.sav'   ; latitude
    endif
    nlv=n_elements(altitude)
;
; compute zonal mean T, Z
;
    tbar=fltarr(nr,nlv)
    nbar=lonarr(nr,nlv)
    for ii=0L,n_elements(id)-1L do begin
        tmask_prof=reform(TEMPERATURE_MASK(ii,*))
        good=where(tmask_prof ne -99.,ngood)
        if good(0) ne -1L then begin
           ymean=latitude(ii)
           for j=0L,nr-1L do begin
               if ymean ge wlat(j)-dy/2. and ymean lt wlat(j)+dy/2. then begin
                  tbar(j,good)=tbar(j,good)+temperature(ii,good)
                  nbar(j,good)=nbar(j,good)+1L
               endif
           endfor
        endif
    endfor
jump:
endfor	; loop over months
good=where(nbar gt 0L)
if good(0) ne -1L then tbar(good)=tbar(good)/float(nbar(good))
;
; interpolate WACCM temperature to height surfaces
;
tbarz=fltarr(nr,nlv)
for kk=0L,nlv-1L do begin
    zz=altitude(kk)
    for j=0L,nr-1L do begin
        zprof=reform(zbarg(j,*))
        if zz gt max(zprof) then goto,jumplev
        for k=1L,nz-1L do begin
            zup=zbarg(j,k-1) & zlw=zbarg(j,k)
            if zup ge zz and zlw le zz then begin
               zscale=(zup-zz)/(zup-zlw)
               tbarz(j,kk)=tbarg(j,k-1)+zscale*(tbarg(j,k)-tbarg(j,k-1))
;print,zlw,zz,zup,zscale
;print,tbarg(j,k),tbarz(j,kk),tbarg(j,k-1)
;stop
            endif
         endfor
      endfor
jumplev:
endfor
;
; save monthly mean gridded MLS data
;
save,file=mdir+'waccm_mls_v2.2_'+syear+smonth+'_temp_aur.sav',wlat,altitude,tbar,tbarz
plotit:
restore,mdir+'waccm_mls_v2.2_'+syear+smonth+'_temp_aur.sav'

if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='yz_waccm_mls_temp_'+syear+smonth+'_aur.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
          xsize=xsize,ysize=ysize
endif
;
; plot temperature
;
erase
xyouts,.3,.7,month(long(smonth)-1)+' '+syear+' Temperature',charsize=2,charthick=2,/normal,color=0
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
nlvls=19
col1=1+indgen(nlvls)*icolmax/nlvls
level=110.+10.*findgen(nlvls)
index=where(tbarz eq 0.)
if index(0) ne -1L then tbarz(index)=-9999.
contour,tbarz,wlat,altitude,/noeras,xrange=[-90.,90.],yrange=[0.,100.],charsize=1.5,color=0,$
      ytitle='Altitude (km)',xticks=6,/fill,c_color=col1,levels=level,min_value=-9999.,$
      xtitle='Latitude',title='WACCM',charthick=2
index=where(level gt 0.)
contour,tbarz,wlat,altitude,levels=level(index),color=0,/follow,/overplot,min_value=-9999.,c_labels=0*level(index)
contour,tbarz,wlat,altitude,levels=[150.],color=mcolor,/follow,/overplot,min_value=-9999.,c_labels=0*level(index)
contour,tbarz,wlat,altitude,levels=[280.],color=0,thick=3,/follow,/overplot,min_value=-9999.,c_labels=0*level(index)
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(K)',charsize=1.5,charthick=2
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor

xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
index=where(tbar eq 0.)
if index(0) ne -1L then tbar(index)=-9999.
contour,tbar,wlat,altitude,/noeras,xrange=[-90.,90.],yrange=[0.,100.],charsize=1.5,color=0,$
      xticks=6,/fill,c_color=col1,levels=level,min_value=-9999.,$
      xtitle='Latitude',title='MLS',charthick=2
index=where(level gt 0.)
contour,tbar,wlat,altitude,levels=level(index),color=0,/follow,/overplot,min_value=-9999.,c_labels=0*level(index)
contour,tbar,wlat,altitude,levels=[150.],color=mcolor,/follow,/overplot,min_value=-9999.,c_labels=0*level(index)
contour,tbar,wlat,altitude,levels=[280.],color=0,thick=3,/follow,/overplot,min_value=-9999.,c_labels=0*level(index)
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(K)',charsize=1.5,charthick=2
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
pdiff=-99.+0.*tbarz
index=where(tbarz gt 0. and tbar gt 0.)
if index(0) eq -1L then goto,jump
pdiff(index)=tbarz(index)-tbar(index)
level=-50.+10.*findgen(11)
!type=2^2+2^3
xmn=xorig(2)
xmx=xorig(2)+xlen
ymn=yorig(2)
ymx=yorig(2)+ylen
set_viewport,xmn,xmx,ymn,ymx
contour,pdiff,wlat,altitude,xrange=[-90.,90.],yrange=[0.,100.],$
        xticks=6,xtitle='Latitude',charsize=1.5,levels=level,/cell_fill,$
        title='WACCM - MLS',c_color=col2,color=0,min_value=-99.,charthick=2
index=where(level gt 0.)
contour,pdiff,wlat,altitude,/overplot,levels=level(index),color=0,/follow,min_value=-99.,c_labels=0*level(index)
index=where(level lt 0.)
contour,pdiff,wlat,altitude,/overplot,levels=level(index),color=mcolor,/follow,min_value=-99.,c_labels=0*level(index)
contour,pdiff,wlat,altitude,/overplot,levels=[0],color=0,thick=3,min_value=-99.
imin=min(level)
imax=max(level)
ymnb=ymn -cbaryoff
ymxb=ymnb+cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,$
      xtitle='(K)',color=0,xticks=n_elements(level)/2,charsize=1.5,charthick=2
ybox=[0,10,10,0,0]
x2=imin
dx=(imax-imin)/(float(n_elements(col2)))
for jj=0L,n_elements(col2)-1 do begin
    xbox=[x2,x2,x2+dx,x2+dx,x2]
    polyfill,xbox,ybox,color=col2(jj)
    x2=x2+dx
endfor
loadct,39

if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim yz_waccm_mls_temp_'+syear+smonth+'_aur.ps -rotate -90 yz_waccm_mls_temp_'+syear+smonth+'_aur.jpg'
endif
endfor	; loop over months
end
