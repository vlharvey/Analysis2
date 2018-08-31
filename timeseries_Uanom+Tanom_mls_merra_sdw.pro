;
; plot timeseries of mean daily temperature and Ubar tendency from MLS, MERRA, and SD-WACCM
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
xorig=[0.1,0.6,0.1,0.6]
yorig=[0.6,0.6,0.1,0.1]
cbaryoff=0.1
cbarydel=0.01
xlen=0.3
ylen=0.3
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
; read Tbar and Ubar previously saved
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
restore,'tbar+ubar_mls_merra_sdw_20071201.sav'
nn=n_elements(sdate_all)
dfs=fltarr(nn)
mlsubarnh=fltarr(nn)
mlsubarsh=fltarr(nn)
mlstbarsh=fltarr(nn)
sdwubarnh=fltarr(nn)
sdwubarsh=fltarr(nn)
sdwtbarsh=fltarr(nn)
;
; compute DFS relative to SH for all dates
;
iyr0=long(strmid(sdate_all(0),0,4))
for i=0L,nn-1L do begin
    iyr=long(strmid(sdate_all(i),0,4))
    imn=long(strmid(sdate_all(i),4,2))
    idy=long(strmid(sdate_all(i),6,2))
    dfs(i)=julday(imn,idy,iyr)-julday(12,21,iyr0)
endfor
;
; extract Ubar at 60N, 10-5 hPa
;
zindex=where(lev le 11. and lev ge 4.)
yindex=where(abs(lat-60.) lt 2.)
for i=0L,nn-1L do sdwubarnh(i)=mean(SDWUT(i,yindex,zindex))
zindex=where(pmls2 le 11. and pmls2 ge 4.)
yindex=where(abs(lat-60.) lt 2.)
for i=0L,nn-1L do mlsubarnh(i)=mean(MLSUT(i,yindex,zindex))
;
; extract Ubar at 60S and 5 hPa
;
zindex=where(lev le 5. and lev ge 4.)
yindex=where(abs(lat-(-60.)) lt 2.)
for i=0L,nn-1L do sdwubarsh(i)=mean(SDWUT(i,yindex,zindex))
zindex=where(pmls2 le 5. and pmls2 ge 4.)
yindex=where(abs(lat-(-60.)) lt 2.)
for i=0L,nn-1L do mlsubarsh(i)=mean(MLSUT(i,yindex,zindex))
;
; read PMC frequencies (mesopause T proxy for now)
;
yindex=where(abs(lat-(-60.)) lt 2.)
for i=0L,nn-1L do begin
    tprof=reform(sdwtt(i,0,*))
    zindex=where(tprof eq min(tprof))
zindex=where(lev le 0.005 and lev ge 0.001)
;   print,lev(zindex)
    sdwtbarsh(i)=mean(SDWTT(i,yindex,zindex))

    tprof=reform(mlstt(i,0,*))
    index=where(finite(tprof) ne 1)
    if index(0) ne -1L then tprof(index)=9999.
    zindex=where(tprof eq min(tprof))
zindex=where(pmls2 le 0.005 and pmls2 ge 0.001)
    print,pmls2(zindex)
    mlstbarsh(i)=mean(MLSTT(i,yindex,zindex))
endfor
;
; postscript file
;
if setplot eq 'ps' then begin
   lc=0
   xsize=nxdim/100.
   ysize=nydim/100.
   set_plot,'ps'
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='timeseries_Uanom+Tanom_mls_merra_sdw_'+sdate_all(0)+'.ps'
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
;
; plot
;
erase
!type=2^2+2^3
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!p.linestyle=0
plot,dfs,mlsubarnh,/noeras,yrange=[-30,120],color=0,xrange=[-30,60],thick=3,ytitle='Zonal Mean Wind (m/s)',title='2007/2008'
xyouts,-20,110,'MLS',color=0,charsize=1,charthick=2,/data
xyouts,-20,100,'SD-WACCM',color=100,charsize=1,charthick=2,/data
oplot,[-30,60],[0,0],color=0
!p.linestyle=5
plots,-20,0
plots,-20,80,/continue,color=0
plots,20,0
plots,20,80,/continue,color=0
plots,-20,0
plots,20,0,/continue,color=0
plots,-20,80
plots,20,80,/continue,color=0
!p.linestyle=0
oplot,dfs,sdwubarnh,color=100,thick=3
!p.linestyle=5
oplot,[-30,60],[0,0],color=0
oplot,dfs,mlsubarsh,thick=3,color=0
oplot,dfs,sdwubarsh,thick=3,color=100
axis,/yax,yrange=[110,200],/save,color=0,ytitle='Mesopause Temperature (K)'
oplot,dfs,smooth(mlstbarsh,3,/edge_truncate),thick=3,color=0
oplot,dfs,smooth(sdwtbarsh,3,/edge_truncate),thick=3,color=100

xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
!p.linestyle=0
plot,dfs,mlsubarnh,/noeras,yrange=[-30,120],color=0,xrange=[-30,60],thick=3,ytitle='Zonal Mean Wind (m/s)',title='2007/2008'
xyouts,-20,110,'MLS',color=0,charsize=1,charthick=2,/data
xyouts,-20,100,'SD-WACCM',color=100,charsize=1,charthick=2,/data
oplot,[-30,60],[0,0],color=0
!p.linestyle=5
plots,20,0
plots,20,80,/continue,color=0
plots,50,0
plots,50,80,/continue,color=0
plots,20,0
plots,50,0,/continue,color=0
plots,20,80
plots,50,80,/continue,color=0
!p.linestyle=0
oplot,dfs,sdwubarnh,color=100,thick=3
!p.linestyle=5
oplot,[-30,60],[0,0],color=0
oplot,dfs,mlsubarsh,thick=3,color=0
oplot,dfs,sdwubarsh,thick=3,color=100
axis,/yax,yrange=[110,200],/save,color=0,ytitle='Mesopause Temperature (K)'
oplot,dfs,smooth(mlstbarsh,3,/edge_truncate),thick=3,color=0
oplot,dfs,smooth(sdwtbarsh,3,/edge_truncate),thick=3,color=100

xmn=xorig(2)
xmx=xorig(2)+xlen
ymn=yorig(2)
ymx=yorig(2)+ylen
set_viewport,xmn,xmx,ymn,ymx
!p.linestyle=0
mlsuanom=mlsubarnh-smooth(mlsubarnh,13,/edge_truncate)
;mlsuanom=smooth(mlsuanom,3,/edge_truncate)
sdwuanom=sdwubarnh-smooth(sdwubarnh,13,/edge_truncate)
;sdwuanom=smooth(sdwuanom,3,/edge_truncate)

mlstanom=mlstbarsh-smooth(mlstbarsh,13,/edge_truncate)
mlstanom=smooth(mlstanom,3,/edge_truncate)
sdwtanom=sdwtbarsh-smooth(sdwtbarsh,13,/edge_truncate)
sdwtanom=smooth(sdwtanom,3,/edge_truncate)

plot,dfs,mlsuanom,/noeras,yrange=[-20,20],color=0,xrange=[-20,20],thick=5,ytitle='Ubar Anomaly (m/s)',/nodata,title='SD-WACCM'
oplot,dfs,sdwuanom,color=0,thick=5
oplot,[-30,60],[0,0],color=0
axis,/yax,yrange=[4,-4],/save,color=0,ytitle='Mesopause Temperature Anomaly (K)'
!p.linestyle=5
oplot,dfs,mlstanom,thick=5,color=150
oplot,dfs-2.,sdwtanom,thick=5,color=100


xmn=xorig(3)
xmx=xorig(3)+xlen
ymn=yorig(3)
ymx=yorig(3)+ylen
set_viewport,xmn,xmx,ymn,ymx
!p.linestyle=0
plot,dfs,mlsuanom,/noeras,yrange=[-30,30],color=0,xrange=[20,50],thick=5,ytitle='Ubar Anomaly (m/s)',/nodata,title='SD-WACCM'
oplot,dfs,sdwuanom,color=0,thick=5
oplot,[-30,60],[0,0],color=0
axis,/yax,yrange=[4,-4],/save,color=0,ytitle='Mesopause Temperature Anomaly (K)'
!p.linestyle=5
oplot,dfs-3.,mlstanom,thick=5,color=150
oplot,dfs-5.,sdwtanom,thick=5,color=100

;
; Close PostScript file and return control to X-windows
if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim timeseries_Uanom+Tanom_mls_merra_sdw_'+sdate_all(0)+'.ps -rotate -90 timeseries_Uanom+Tanom_mls_merra_sdw_'+sdate_all(0)+'.jpg'
;  spawn,'rm -f timeseries_Uanom+Tanom_mls_merra_sdw_'+sdate_all(0)+'.ps'
endif
end
