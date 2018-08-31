;
; zt of Elat edge values based on CO and SF
;
loadct,39
mcolor=byte(!p.color)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
nxdim=700
nydim=700
xorig=[0.15,0.60]
yorig=[0.35,0.35]
cbaryoff=0.08
cbarydel=0.01
xlen=0.3
ylen=0.4
PI2=6.2831853071796
DTR=PI2/360.
RADEA=6.37E6
!NOERAS=-1
syear=['2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014']
nyear=n_elements(syear)
set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif
goto,quick
;
; get file listing
;
spawn,'ls daily_mls_coelatedge+merra_sfelatedge_??????_2d.sav',ifiles
;
; DELATLEVS3D     FLOAT     = Array[29, 91, 37]
; LOWLAT_ELATEDGE_2D FLOAT  = Array[29, 37]
; NOVORTEX_FLAG_2D FLOAT    = Array[29, 37]
; PMLS            FLOAT     = Array[37]
; SDATE_TIME      STRING    = Array[29]
; SFELATEDGE_2D   FLOAT     = Array[29, 37]
; SPBIN3D         FLOAT     = Array[29, 91, 37]
; YEQ             FLOAT     = Array[91]
;
restore,ifiles(0)
lowlat_elatedge_2d_all=LOWLAT_ELATEDGE_2d
NOVORTEX_FLAG_2D_all=NOVORTEX_FLAG_2D
sfelatedge_2d_all=SFELATEDGE_2d
sdate_time_all=sdate_time

for ifile=1L,n_elements(ifiles)-1L do begin
    restore,ifiles(ifile)
    print,ifiles(ifile)
    LOWLAT_ELATEDGE_2d_all=[LOWLAT_ELATEDGE_2d_all,LOWLAT_ELATEDGE_2d]
    NOVORTEX_FLAG_2D_all=[NOVORTEX_FLAG_2D_all,NOVORTEX_FLAG_2D]
    SFELATEDGE_2d_all=[SFELATEDGE_2d_all,sfELATEDGE_2d]
    sdate_time_all=[sdate_time_all,sdate_time]
endfor

save,file='zt_coelat_sfelat_2d_stats.sav',sdate_time_all,SFELATEDGE_2d_all,NOVORTEX_FLAG_2D_all,LOWLAT_ELATEDGE_2d_all,pmls
quick:
restore,'zt_coelat_sfelat_2d_stats.sav

; postscript file
;
if setplot eq 'ps' then begin
   lc=0
   xsize=nxdim/100.
   ysize=nydim/100.
   set_plot,'ps'
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='zt_coelat_sfelat_2d_stats.ps'
   !p.charsize=1.2
   !p.thick=2
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
;
; DJF
;
syear=strmid(sdate_time_all,0,4)
smon=strmid(sdate_time_all,4,2)
sday=strmid(sdate_time_all,6,2)
index=where(smon ne '11' and smon ne '03')
index=where(smon ne '99')
LOWLAT_ELATEDGE_2d_all=LOWLAT_ELATEDGE_2d_all(index,*)
SFELATEDGE_2d_all=SFELATEDGE_2d_all(index,*)
sdate_time_all=sdate_time_all(index)
syear=strmid(sdate_time_all,0,4)
smon=strmid(sdate_time_all,4,2)
sday=strmid(sdate_time_all,6,2)

erase
xindex=where(strmid(sdate_time_all,4,4) eq '1201',nxticks)
xlabs0=strmid(sdate_time_all(xindex),2,2)
xindex=where(strmid(sdate_time_all,4,4) eq '0101',nxticks)
xlabs1=strmid(sdate_time_all(xindex),2,2)
xlabs=xlabs0+xlabs1
nlvls=21
col1=reverse(1+indgen(nlvls)*mcolor/nlvls)
imin=25.
imax=85.
level=imin+((imax-imin)/(nlvls-1.))*findgen(nlvls)
;
index=where(NOVORTEX_FLAG_2D_all eq -9999. or LOWLAT_ELATEDGE_2d_ALL eq 0. or LOWLAT_ELATEDGE_2d_ALL eq -99.)
if index(0) ne -1L then LOWLAT_ELATEDGE_2d_ALL(index)=0./0.
LOWLAT_ELATEDGE_2d_SAVE=LOWLAT_ELATEDGE_2d_ALL
;
x2d=0.*LOWLAT_ELATEDGE_2d_ALL
y2d=0.*LOWLAT_ELATEDGE_2d_ALL
for i=0L,n_elements(sdate_time_all)-1L do y2d(i,*)=pmls
for j=0L,n_elements(pmls)-1L do x2d(*,j)=findgen(n_elements(sdate_time_all))
LOWLAT_ELATEDGE_2d_ALL=smooth(LOWLAT_ELATEDGE_2d_ALL,5,/Nan)
index=where(finite(LOWLAT_ELATEDGE_2d_SAVE) ne 1)
if index(0) ne -1L then LOWLAT_ELATEDGE_2d_ALL(index)=0./0.
;contour,LOWLAT_ELATEDGE_2d_ALL,findgen(n_elements(sdate_time_all)),pmls,/ylog,color=0,/noeras,charsize=1.25,ytitle='Pressure (hPa)',yrange=[10,0.01],$
;     title='CO Vortex Elat',xtickname=xlabs,xtickv=xindex,xticks=nxticks-1,levels=level,c_color=col1,/cell_fill,charthick=2
;index=where(NOVORTEX_FLAG_2D_all le 0.)
;;if index(0) ne -1L then oplot,x2d(index),y2d(index),psym=8,color=mcolor
;;contour,LOWLAT_ELATEDGE_2d_ALL,findgen(n_elements(sdate_time_all)),pmls,/ylog,/overplot,levels=[60],color=0,/follow,c_labels=[0],thick=2
;;
;; indicate new year
;;
;index=where(strmid(sdate_time_all,4,4) eq '1001',nxticks)
;for i=0L,nxticks-1L do begin
;    plots,index(i),0.01
;    plots,index(i),10,/continue,color=0,thick=3
;endfor

smonth=strmid(SDATE_TIME_ALL,4,2)
index=where(smon eq '10')
LOWLAT_ELATEDGE_2D_10=reform(LOWLAT_ELATEDGE_2D_ALL(index,*))
SFELATEDGE_2D_10=reform(SFELATEDGE_2D_ALL(index,*))
index=where(smon eq '11')                                    
LOWLAT_ELATEDGE_2D_11=reform(LOWLAT_ELATEDGE_2D_ALL(index,*))
SFELATEDGE_2D_11=reform(SFELATEDGE_2D_ALL(index,*))
index=where(smon eq '12')                                    
LOWLAT_ELATEDGE_2D_12=reform(LOWLAT_ELATEDGE_2D_ALL(index,*))
SFELATEDGE_2D_12=reform(SFELATEDGE_2D_ALL(index,*))
index=where(smon eq '01')                                    
LOWLAT_ELATEDGE_2D_01=reform(LOWLAT_ELATEDGE_2D_ALL(index,*))
SFELATEDGE_2D_01=reform(SFELATEDGE_2D_ALL(index,*))
index=where(smon eq '02')                                    
LOWLAT_ELATEDGE_2D_02=reform(LOWLAT_ELATEDGE_2D_ALL(index,*))
SFELATEDGE_2D_02=reform(SFELATEDGE_2D_ALL(index,*))
index=where(smon eq '03')                                    
LOWLAT_ELATEDGE_2D_03=reform(LOWLAT_ELATEDGE_2D_ALL(index,*))
SFELATEDGE_2D_03=reform(SFELATEDGE_2D_ALL(index,*))
nl=n_elements(pmls)
octprof=fltarr(nl)
novprof=fltarr(nl)
decprof=fltarr(nl)
janprof=fltarr(nl)
febprof=fltarr(nl)
marprof=fltarr(nl)
for k=0L,nl-1L do begin
    index=where(long(finite(LOWLAT_ELATEDGE_2D_10(*,k)) eq 1L))
    octprof(k)=float(n_elements(index))/float(n_elements(LOWLAT_ELATEDGE_2D_10(*,k)))
    index=where(long(finite(LOWLAT_ELATEDGE_2D_11(*,k)) eq 1L))
    novprof(k)=float(n_elements(index))/float(n_elements(LOWLAT_ELATEDGE_2D_11(*,k)))
    index=where(long(finite(LOWLAT_ELATEDGE_2D_12(*,k)) eq 1L))
    decprof(k)=float(n_elements(index))/float(n_elements(LOWLAT_ELATEDGE_2D_12(*,k)))
    index=where(long(finite(LOWLAT_ELATEDGE_2D_01(*,k)) eq 1L))
    janprof(k)=float(n_elements(index))/float(n_elements(LOWLAT_ELATEDGE_2D_01(*,k)))
    index=where(long(finite(LOWLAT_ELATEDGE_2D_02(*,k)) eq 1L))
    febprof(k)=float(n_elements(index))/float(n_elements(LOWLAT_ELATEDGE_2D_02(*,k)))
    index=where(long(finite(LOWLAT_ELATEDGE_2D_03(*,k)) eq 1L))
    marprof(k)=float(n_elements(index))/float(n_elements(LOWLAT_ELATEDGE_2D_03(*,k)))
endfor
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
plot,octprof,pmls,color=0,/ylog,yrange=[10.,0.01],/noeras,thick=5,charthick=2,ytitle='Pressure (hPa)',xrange=[0.,1.],charsize=1.25,xtitle='Fraction'
oplot,novprof,pmls,color=50,thick=5
oplot,decprof,pmls,color=90,thick=5
oplot,janprof,pmls,color=150,thick=5
oplot,febprof,pmls,color=200,thick=5
oplot,marprof,pmls,color=250,thick=5
xyouts,xmn+0.01,ymx-0.02,'Oct',/normal,color=0,charthick=2,charsize=1.2
xyouts,xmn+0.01,ymx-0.04,'Nov',/normal,color=50,charthick=2,charsize=1.2
xyouts,xmn+0.01,ymx-0.06,'Dec',/normal,color=90,charthick=2,charsize=1.2
xyouts,xmn+0.01,ymx-0.08,'Jan',/normal,color=150,charthick=2,charsize=1.2
xyouts,xmn+0.01,ymx-0.10,'Feb',/normal,color=200,charthick=2,charsize=1.2
xyouts,xmn+0.01,ymx-0.12,'Mar',/normal,color=250,charthick=2,charsize=1.2
;
; avg monthly elat as a function of altitude
;
octprof=fltarr(nl)
novprof=fltarr(nl)
decprof=fltarr(nl)
janprof=fltarr(nl)
febprof=fltarr(nl)
marprof=fltarr(nl)
for k=0L,nl-1L do begin
    index=where(long(finite(SFELATEDGE_2D_10(*,k))) eq 1L)
    if index(0) ne -1L then octprof(k)=mean(SFELATEDGE_2D_10(index,k))
    index=where(long(finite(SFELATEDGE_2D_11(*,k))) eq 1L)
    if index(0) ne -1L then novprof(k)=mean(SFELATEDGE_2D_11(index,k))
    index=where(long(finite(SFELATEDGE_2D_12(*,k))) eq 1L)
    if index(0) ne -1L then decprof(k)=mean(SFELATEDGE_2D_12(index,k))
    index=where(long(finite(SFELATEDGE_2D_01(*,k))) eq 1L)
    if index(0) ne -1L then janprof(k)=mean(SFELATEDGE_2D_01(index,k))
    index=where(long(finite(SFELATEDGE_2D_02(*,k))) eq 1L)
    if index(0) ne -1L then febprof(k)=mean(SFELATEDGE_2D_02(index,k))
    index=where(long(finite(SFELATEDGE_2D_03(*,k))) eq 1L)
    if index(0) ne -1L then marprof(k)=mean(SFELATEDGE_2D_03(index,k))
endfor
xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
plot,octprof,pmls,color=0,/ylog,yrange=[10.,0.01],ytitle='Pressure (hPa)',/noeras,thick=5,charthick=2,xrange=[30.,90.],charsize=1.25,xtitle='Mean Elat',linestyle=5
oplot,novprof,pmls,color=50,thick=5,linestyle=5
oplot,decprof,pmls,color=90,thick=5,linestyle=5
oplot,janprof,pmls,color=150,thick=5,linestyle=5
oplot,febprof,pmls,color=200,thick=5,linestyle=5
oplot,marprof,pmls,color=250,thick=5,linestyle=5

octprof=fltarr(nl)*0./0.
novprof=fltarr(nl)*0./0.
decprof=fltarr(nl)*0./0.
janprof=fltarr(nl)*0./0.
febprof=fltarr(nl)*0./0.
marprof=fltarr(nl)*0./0.
for k=0L,nl-1L do begin
    index=where(long(finite(LOWLAT_ELATEDGE_2D_10(*,k))) eq 1L)
    if index(0) ne -1L then octprof(k)=mean(LOWLAT_ELATEDGE_2D_10(index,k))
    index=where(long(finite(LOWLAT_ELATEDGE_2D_11(*,k))) eq 1L)
    if index(0) ne -1L then novprof(k)=mean(LOWLAT_ELATEDGE_2D_11(index,k))
    index=where(long(finite(LOWLAT_ELATEDGE_2D_12(*,k))) eq 1L)
    if index(0) ne -1L then decprof(k)=mean(LOWLAT_ELATEDGE_2D_12(index,k))
    index=where(long(finite(LOWLAT_ELATEDGE_2D_01(*,k))) eq 1L)
    if index(0) ne -1L then janprof(k)=mean(LOWLAT_ELATEDGE_2D_01(index,k))
    index=where(long(finite(LOWLAT_ELATEDGE_2D_02(*,k))) eq 1L)
    if index(0) ne -1L then febprof(k)=mean(LOWLAT_ELATEDGE_2D_02(index,k))
    index=where(long(finite(LOWLAT_ELATEDGE_2D_03(*,k))) eq 1L)
    if index(0) ne -1L then marprof(k)=mean(LOWLAT_ELATEDGE_2D_03(index,k))
endfor
oplot,octprof,pmls,color=0,thick=5
oplot,novprof,pmls,color=50,thick=5
oplot,decprof,pmls,color=90,thick=5
oplot,janprof,pmls,color=150,thick=5
oplot,febprof,pmls,color=200,thick=5
oplot,marprof,pmls,color=250,thick=5
xyouts,xmx-0.05,ymx-0.02,'Oct',/normal,color=0,charthick=2,charsize=1.2
xyouts,xmx-0.05,ymx-0.04,'Nov',/normal,color=50,charthick=2,charsize=1.2
xyouts,xmx-0.05,ymx-0.06,'Dec',/normal,color=90,charthick=2,charsize=1.2
xyouts,xmx-0.05,ymx-0.08,'Jan',/normal,color=150,charthick=2,charsize=1.2
xyouts,xmx-0.05,ymx-0.10,'Feb',/normal,color=200,charthick=2,charsize=1.2
xyouts,xmx-0.05,ymx-0.12,'Mar',/normal,color=250,charthick=2,charsize=1.2
;
;
; Close PostScript file and return control to X-windows
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim zt_coelat_sfelat_2d_stats.ps -rotate -90 zt_coelat_sfelat_2d_stats.jpg
;  spawn,'rm -f zt_coelat_sfelat_2d_stats.ps'
endif
end
