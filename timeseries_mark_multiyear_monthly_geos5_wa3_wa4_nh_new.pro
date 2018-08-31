;
; monthly mean timeseries of the latitude of the Arctic vortex at a specified altitude
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3_meto

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
xorig=[0.2]
yorig=[0.2]
xlen=0.7
ylen=0.7
cbaryoff=0.125
cbarydel=0.01
setplot='x'
read,'setplot=',setplot
if setplot ne 'ps' then begin
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=icolmax
endif
dir='/aura7/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS5?0.MetO.'
dir='/Volumes/earth/harvey/GEOS5_data/Datfiles/DAS.ops.asm.tavg3d_dyn_v.GEOS5?0.MetO.'
month=['January','February','March','April','May','June','July','August','September','October','November','December']
xlab=['J','F','M','A','M','J','J','A','S','O','N','D']
xlab=['J','A','S','O','N','D','J','F','M','A','M','J']
mday=[31,28,31,30,31,30,31,31,30,31,30,31]
jday=[15,46,74,105,135,166,196,227,258,288,319,349]	; 15th of each month
jday=[196,227,258,288,319,349,15,46,74,105,135,166]	; 15th of each month

for ii=0L,11L do begin
;
;***Read GEOS data
;
    smn=string(FORMAT='(I2.2)',ii+1)
    restore,'yz_geos5_ubarth+mark_multiyear_'+smn+'.sav'	;,alat,th,tzm,uzm,markzm
    zindex=where(th eq 3000. or th eq 1800. or th eq 700.,nth2)
    if ii eq 0L then begin
       nr=n_elements(alat)
       nth=n_elements(th)
       gmarklat=fltarr(12,nr,nth2)
       w3marklat=fltarr(12,nr,nth2)
       w4marklat=fltarr(12,nr,nth2)
       zangle=fltarr(12,nr)
       sth='700-3000K'
    endif
    for j=0,nr-1L do gmarklat(ii,j,*)=markzm(j,zindex)
    restore,'yz_wa3_ubarth+mark_'+month(ii)+'.sav'        ;,alat,th,tzm,uzm,markzm
    zindex=where(th eq 3000. or th eq 1800. or th eq 700.,nth2)
    for j=0,nr-1L do w3marklat(ii,j,*)=markzm(j,zindex)
    restore,'yz_wa4_ubarth+mark_'+month(ii)+'.sav'        ;,alat,th,tzm,uzm,markzm
    zindex=where(th eq 3000. or th eq 1800. or th eq 700.,nth2)
    for j=0,nr-1L do w4marklat(ii,j,*)=markzm(j,zindex)
;
; compute solar zenith angle
;
      doy=jday(ii)
      pi=3.14159265
      dtor=pi/180.
      earinc=23.5
      for j=0L,nr-1 do begin
          rlat=alat(j)
          rlon=0.
          gmt=12.
          sinlat=sin(rlat*dtor)
          coslat=sqrt(1.-sinlat^2.)
          sinlon=sin(rlon*dtor)
          coslon=cos(rlon*dtor)
          soya=(doy-81.25)*pi/182.5           ; day angle
          soha=2.*pi*(gmt-12.)/24.            ; hour angle
          soha=-soha
          sininc=sin(earinc*dtor)
          sindec=sininc*sin(soya)
          cosdec= sqrt(1.-sindec^2.)
          coszen=cos(soha)*coslon+sin(soha)*sinlon
          coszen=coszen*cosdec*coslat
          coszen=sindec*sinlat+coszen
          coszen=min([max([coszen,-1.]),1.])
          chi = acos(coszen)
          zangle(ii,j) = chi/dtor
      endfor
endfor  ; loop over months

if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='timeseries_mark_multiyear_monthly_'+sth+'_nh.ps'
   !p.charsize=1.5
   !p.thick=2
   !p.charthick=5
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
;
; flip from Jan-Dec to Jul-Jun
;
gmarklat2=0.*gmarklat
w3marklat2=0.*gmarklat
w4marklat2=0.*gmarklat
for j=0L,nr-1L do begin
    gmarklat2(0:5,j,*)=gmarklat(6:11,j,*)
    w3marklat2(0:5,j,*)=w3marklat(6:11,j,*)
    w4marklat2(0:5,j,*)=w4marklat(6:11,j,*)
    gmarklat2(6:11,j,*)=gmarklat(0:5,j,*)
    w3marklat2(6:11,j,*)=w3marklat(0:5,j,*)
    w4marklat2(6:11,j,*)=w4marklat(0:5,j,*)
endfor
gmarklat=gmarklat2
w3marklat=w3marklat2
w4marklat=w4marklat2
;
; plot 
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
nlvls=20
col1=1+indgen(nlvls)*icolmax/nlvls
level=0.05+0.05*findgen(nlvls)
index=where(gmarklat eq 0.)
if index(0) ne -1L then gmarklat(index)=0./0.
contour,gmarklat(*,*,0),1+findgen(12),alat,/nodata,/noeras,xrange=[1.,12.],yrange=[20.,90.],charsize=2,color=0,/cell_fill,c_color=col1,$
      ytitle='Latitude',xticks=n_elements(xlab)-1,xtickname=xlab,title='Arctic Vortex '+sth,levels=level
contour,w3marklat(*,*,0),1+findgen(12),alat,/overplot,/follow,levels=[0.5],color=250,thick=15,c_linestyle=0,c_labels=[0]
contour,w4marklat(*,*,0),1+findgen(12),alat,/overplot,/follow,levels=[0.5],color=70,thick=15,c_linestyle=0,c_labels=[0]
contour,gmarklat(*,*,0),1+findgen(12),alat,/overplot,/follow,levels=[0.5],color=0,thick=15,c_linestyle=0,c_labels=[0]
contour,w3marklat(*,*,1),1+findgen(12),alat,/overplot,/follow,levels=[0.5],color=250,thick=15,c_linestyle=5,c_labels=[0]
contour,w4marklat(*,*,1),1+findgen(12),alat,/overplot,/follow,levels=[0.5],color=70,thick=15,c_linestyle=5,c_labels=[0]
contour,gmarklat(*,*,1),1+findgen(12),alat,/overplot,/follow,levels=[0.5],color=0,thick=15,c_linestyle=5,c_labels=[0]
contour,w3marklat(*,*,2),1+findgen(12),alat,/overplot,/follow,levels=[0.5],color=250,thick=15,c_linestyle=1,c_labels=[0]
contour,w4marklat(*,*,2),1+findgen(12),alat,/overplot,/follow,levels=[0.5],color=70,thick=15,c_linestyle=1,c_labels=[0]
contour,gmarklat(*,*,2),1+findgen(12),alat,/overplot,/follow,levels=[0.5],color=0,thick=15,c_linestyle=1,c_labels=[0]
loadct,0
contour,zangle,1+findgen(12),alat,/overplot,/follow,levels=[90],color=0,thick=20,c_linestyle=0,c_labels=[0]
contour,zangle,1+findgen(12),alat,/overplot,/follow,levels=[90],color=150,thick=20,c_linestyle=2,c_labels=[0]
xyouts,5.,80,'SZA=90',/data,color=150,charsize=2,charthick=2
loadct,39
xyouts,1.2,22.,'GEOS',/data,color=0,charsize=2,charthick=2
xyouts,1.2,30.,'WA-3',/data,color=250,charsize=2,charthick=2
xyouts,1.2,26.,'WA-4',/data,color=70,charsize=2,charthick=2

if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim timeseries_mark_multiyear_monthly_'+sth+'_nh.ps -rotate -90 timeseries_mark_multiyear_monthly_'+sth+'_nh.jpg'
endif
end
