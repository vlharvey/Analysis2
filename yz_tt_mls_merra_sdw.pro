;
; plot zonal mean daily temperature tendency from MLS, MERRA, and SD-WACCM
;
@stddat
@kgmt
@ckday
@kdate

loadct,39
mcolor=byte(!p.color)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
nxdim=700
nydim=700
xorig=[0.15,0.425,0.725]
yorig=[0.3,0.3,0.3]
cbaryoff=0.1
cbarydel=0.01
xlen=0.25
ylen=0.4
PI2=6.2831853071796
DTR=PI2/360.
RADEA=6.37E6
!NOERAS=-1
set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif
;
; read Tbar previously saved
;
; LAT             DOUBLE    = Array[96]
; LATITUDE_WACCM  FLOAT     = Array[96]
; LEV             DOUBLE    = Array[88]
; MERRATT         FLOAT     = Array[62, 96, 41]
; MLAT            DOUBLE    = Array[96]
; MLSTT           FLOAT     = Array[62, 96, 55]
; PMLS2           FLOAT     = Array[55]
; PRESSURE        FLOAT     = Array[41]
; SDATE_ALL       STRING    = Array[62]
; SDWTT           FLOAT     = Array[62, 96, 88]
;
restore,'tbar_mls_merra_sdw_20071201.sav'
nday=n_elements(SDATE_ALL)
for n=0L+7L,nday-1L-7L do begin
    sdate=SDATE_ALL(n)
    print,sdate
;
; zonal mean temperature the week prior to today
;
    mtbar0=fltarr(n_elements(mlat),n_elements(pmls2))
    metbar0=fltarr(n_elements(latitude_waccm),n_elements(pressure))
    wtbar0=fltarr(n_elements(lat),n_elements(lev))
    for nn=n-7L,n-1L do begin
        mtbar0=mtbar0+reform(mlstt(nn,*,*))
        metbar0=metbar0+reform(merratt(nn,*,*))
        wtbar0=wtbar0+reform(sdwtt(nn,*,*))
print,'week before ',sdate_all(nn)
    endfor
    mtbar0=mtbar0/7.
    metbar0=metbar0/7.
    wtbar0=wtbar0/7.
;
; zonal mean temperature the week after today
;
    mtbar1=fltarr(n_elements(mlat),n_elements(pmls2))
    metbar1=fltarr(n_elements(latitude_waccm),n_elements(pressure))
    wtbar1=fltarr(n_elements(lat),n_elements(lev))
    for nn=n,n+7L-1L do begin
        mtbar1=mtbar1+reform(mlstt(nn,*,*))
        metbar1=metbar1+reform(merratt(nn,*,*))
        wtbar1=wtbar1+reform(sdwtt(nn,*,*))
print,'week after ',sdate_all(nn)
    endfor
    mtbar1=mtbar1/7.
    metbar1=metbar1/7.
    wtbar1=wtbar1/7.
;
; today's temperature to identify summer mesopause
;
mlstbar0=reform(mlstt(n,*,*))
sdwtbar0=reform(sdwtt(n,*,*))
;
; postscript file
;
        if setplot eq 'ps' then begin
           lc=0
           xsize=nxdim/100.
           ysize=nydim/100.
           set_plot,'ps'
           device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
                  /bold,/color,bits_per_pixel=8,/helvetica,filename='yz_tt_mls_merra_sdw_'+sdate+'.ps'
           !p.charsize=1.25
           !p.thick=2
           !p.charthick=5
           !y.thick=2
           !x.thick=2
        endif
;
; plot
;
       imin=-20.
       level=imin+2.*findgen(21)
level=[-50.,-40.,-30.,-25.,level,25.,30.,40.,50.]
imin=min(level)
       imax=max(level)
       nlvls=n_elements(level)
       col1=1+(indgen(nlvls)/float(nlvls))*mcolor

       erase
       xyouts,.4,.8,sdate,color=0,charsize=2,charthick=2,/normal
       xmn=xorig(0)
       xmx=xorig(0)+xlen
       ymn=yorig(0)
       ymx=yorig(0)+ylen
       set_viewport,xmn,xmx,ymn,ymx
       !type=2^2+2^3
mlstbar=mtbar1-mtbar0
print,min(mlstbar),max(mlstbar)
       contour,mlstbar,mlat,pmls2,levels=level,/noeras,c_color=col1,/cell_fill,yrange=[100.,1.e-4],/ylog,color=0,title='MLS'
       index=where(level gt 0.)
       contour,mlstbar,mlat,pmls2,levels=level(index),/overplot,/noeras,c_color=0,/follow
       index=where(level lt 0.)
       contour,mlstbar,mlat,pmls2,levels=level(index),/overplot,/noeras,color=mcolor,/follow,c_linestyle=5
contour,mlstbar0,mlat,pmls2,levels=[130.,140.,150.,160.],/overplot,/noeras,color=0,thick=3

       xmn=xorig(1)
       xmx=xorig(1)+xlen
       ymn=yorig(1)
       ymx=yorig(1)+ylen
       set_viewport,xmn,xmx,ymn,ymx
       !type=2^2+2^3
merratbar=metbar1-metbar0
print,min(merratbar),max(merratbar)
       contour,merratbar,latitude_waccm,pressure,levels=level,/noeras,c_color=col1,/cell_fill,yrange=[100.,1.e-4],/ylog,color=0,title='MERRA',ytickname=[' ',' ',' ',' ',' ',' ',' ']
       index=where(level gt 0.)
       contour,merratbar,latitude_waccm,pressure,levels=level(index),/overplot,/noeras,c_color=0,/follow
       index=where(level lt 0.)
       contour,merratbar,latitude_waccm,pressure,levels=level(index),/overplot,/noeras,color=mcolor,/follow,c_linestyle=5

       xmn=xorig(2)
       xmx=xorig(2)+xlen
       ymn=yorig(2)
       ymx=yorig(2)+ylen
       set_viewport,xmn,xmx,ymn,ymx
       !type=2^2+2^3
sdwtbar=wtbar1-wtbar0
print,min(sdwtbar),max(sdwtbar)
       contour,sdwtbar,lat,lev,levels=level,/noeras,c_color=col1,/cell_fill,yrange=[100.,1.e-4],/ylog,color=0,title='SD-WACCM',ytickname=[' ',' ',' ',' ',' ',' ',' ']
       index=where(level gt 0.)
       contour,sdwtbar,lat,lev,levels=level(index),/overplot,/noeras,c_color=0,/follow
       index=where(level lt 0.)
       contour,sdwtbar,lat,lev,levels=level(index),/overplot,/noeras,color=mcolor,/follow,c_linestyle=5
contour,sdwtbar0,lat,lev,levels=[130.,140.,150.,160.],/overplot,/noeras,color=0,thick=3


       ymnb=ymn -cbaryoff
       ymxb=ymnb+cbarydel
       set_viewport,xorig(0),xorig(2)+xlen,ymnb,ymxb
       !type=2^2+2^3+2^6
       plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,xtitle='Temperature Tendency (K)
       ybox=[0,10,10,0,0]
       x2=imin
       dx=(imax-imin)/(float(nlvls)-1)
       for j=1,nlvls-1 do begin
           xbox=[x2,x2,x2+dx,x2+dx,x2]
           polyfill,xbox,ybox,color=col1(j)
           x2=x2+dx
       endfor
;
; Close PostScript file and return control to X-windows
       if setplot ne 'ps' then stop	;wait,1
       if setplot eq 'ps' then begin
          device, /close
          spawn,'convert -trim yz_tt_mls_merra_sdw_'+sdate+'.ps -rotate -90 '+$
                              'yz_tt_mls_merra_sdw_'+sdate+'.jpg'
          spawn,'rm -f yz_tt_mls_merra_sdw_'+sdate+'.ps'
       endif

endfor

end
