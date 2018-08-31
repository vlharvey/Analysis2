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
xorig=[0.15]
yorig=[0.35]
cbaryoff=0.08
cbarydel=0.01
xlen=0.8
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
HLATPDF_TIME_3D_all=HLATPDF_TIME_3D
index=where(yeq mod 2 eq 0,nrr2)
yeq2=yeq(index)                         ; 2 degree bins

for ifile=1L,n_elements(ifiles)-1L do begin
    restore,ifiles(ifile)
    print,ifiles(ifile)
    LOWLAT_ELATEDGE_2d_all=[LOWLAT_ELATEDGE_2d_all,LOWLAT_ELATEDGE_2d]
    NOVORTEX_FLAG_2D_all=[NOVORTEX_FLAG_2D_all,NOVORTEX_FLAG_2D]
    SFELATEDGE_2d_all=[SFELATEDGE_2d_all,sfELATEDGE_2d]
    sdate_time_all=[sdate_time_all,sdate_time]
    HLATPDF_TIME_3D_all=[HLATPDF_TIME_3D_all,HLATPDF_TIME_3D]
help,HLATPDF_TIME_3D_all
endfor
save,file='figure_7.sav',sdate_time_all,SFELATEDGE_2d_all,NOVORTEX_FLAG_2D_all,LOWLAT_ELATEDGE_2d_all,pmls,yeq2,HLATPDF_TIME_3D_all
quick:
restore,'figure_7.sav

HPDF_FLAG_2D_ALL=0.*NOVORTEX_FLAG_2D_all
nday=n_elements(sdate_time_all)
for i=0L,nday-1L do begin
for k=0L,n_elements(pmls)-1L do begin
    pdfarray=reform(HLATPDF_TIME_3D_all(i,*,k))
    index=where(pdfarray ne 0.)
    if index(0) ne -1L then begin
       if min(yeq2(index)) gt 20. then HPDF_FLAG_2D_ALL(i,k)=1.		; good confined to high lats
    endif
endfor
endfor

; postscript file
;
if setplot eq 'ps' then begin
   lc=0
   xsize=nxdim/100.
   ysize=nydim/100.
   set_plot,'ps'
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='zt_coelat_sfelat_2d.ps'
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
index=where(smon ne '09')
LOWLAT_ELATEDGE_2d_all=LOWLAT_ELATEDGE_2d_all(index,*)
SFELATEDGE_2d_all=SFELATEDGE_2d_all(index,*)
NOVORTEX_FLAG_2D_all=NOVORTEX_FLAG_2D_all(index,*)
HPDF_FLAG_2D_all=HPDF_FLAG_2D_all(index,*)
sdate_time_all=sdate_time_all(index)
syear=strmid(sdate_time_all,0,4)
smon=strmid(sdate_time_all,4,2)
sday=strmid(sdate_time_all,6,2)

erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
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
;index=where(NOVORTEX_FLAG_2D_all eq -9999. or LOWLAT_ELATEDGE_2d_ALL eq 0. or LOWLAT_ELATEDGE_2d_ALL eq -99.)
index=where(HPDF_FLAG_2D_all eq 0. or LOWLAT_ELATEDGE_2d_ALL eq 0. or LOWLAT_ELATEDGE_2d_ALL eq -99.)
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
contour,LOWLAT_ELATEDGE_2d_ALL,findgen(n_elements(sdate_time_all)),pmls,/ylog,color=0,/noeras,charsize=1.25,ytitle='Pressure (hPa)',yrange=[10,0.01],$
     title='CO Vortex Elat',xtickname=xlabs,xtickv=xindex,xticks=nxticks-1,levels=level,c_color=col1,/cell_fill,charthick=2
index=where(NOVORTEX_FLAG_2D_all le 0.)
;if index(0) ne -1L then oplot,x2d(index),y2d(index),psym=8,color=mcolor
;contour,LOWLAT_ELATEDGE_2d_ALL,findgen(n_elements(sdate_time_all)),pmls,/ylog,/overplot,levels=[60],color=0,/follow,c_labels=[0],thick=2
;
; indicate new year
;
index=where(strmid(sdate_time_all,4,4) eq '1001',nxticks)
for i=0L,nxticks-1L do begin
    plots,index(i),0.01
    plots,index(i),10,/continue,color=0,thick=3
endfor
;
; SF Elat
;
;xmn=xorig(1)
;xmx=xorig(1)+xlen
;ymn=yorig(1)
;ymx=yorig(1)+ylen
;set_viewport,xmn,xmx,ymn,ymx
;!type=2^2+2^3
;index=where(SFELATEDGE_2d_ALL eq 0.)
;if index(0) ne -1L then SFELATEDGE_2d_ALL(index)=0./0.
;SFELATEDGE_2d_ALL=smooth(SFELATEDGE_2d_ALL,5,/Nan)
;contour,SFELATEDGE_2d_ALL,findgen(n_elements(sdate_time_all)),pmls,/ylog,color=0,/noeras,charsize=1.25,ytitle='Pressure (hPa)',yrange=[10,0.01],$
;     title='SF Vortex Elat',xtickname=xlabs,xtickv=xindex,xticks=nxticks-1,levels=level,c_color=col1,/cell_fill,charthick=2
;;contour,SFELATEDGE_2d_ALL,findgen(n_elements(sdate_time_all)),pmls,/ylog,/overplot,levels=[60],color=0,/follow,c_labels=[0],thick=2
;
;; indicate new year
;;
;index=where(strmid(sdate_time_all,4,4) eq '1001',nxticks)
;for i=0L,nxticks-1L do begin
;    plots,index(i),0.01
;    plots,index(i),10,/continue,color=0,/data,thick=3
;endfor

       ymnb=ymn -cbaryoff
       ymxb=ymnb+cbarydel
       set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
       !type=2^2+2^3+2^6
       plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1.5,xtitle='Equivalent Latitude',charthick=2
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
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim zt_coelat_sfelat_2d.ps -rotate -90 zt_coelat_sfelat_2d.jpg
;  spawn,'rm -f zt_coelat_sfelat_2d.ps'
endif
end
