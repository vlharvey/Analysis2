;
; plot daily UKMO Ubar
;
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
loadct,39
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
device,decompose=0
;setplot='ps'
;read,'setplot=',setplot
nxdim=700
nydim=700
xorig=[0.15]
yorig=[0.2]
xlen=0.7
ylen=0.7
cbaryoff=0.1
cbarydel=0.01
!NOERAS=-1
;if setplot ne 'ps' then begin
;   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
;   !p.background=mcolor
;endif
ofile='/atmos/harvey/UKMO_data/Datfiles/ukmo_12Z_Ubar_Tbar_3D.sav'
restore,ofile
kday=n_elements(sdate)
latjet=-99.+0.*fltarr(kday)
pjet=-99.+0.*fltarr(kday)
zzero=-99.+0.*fltarr(kday)
nr=n_elements(wlat)
nl=n_elements(p)
wlat2d=fltarr(nr,nl)
p2d=fltarr(nr,nl)
for j=0L,nr-1L do p2d(j,*)=p
for k=0L,nl-1L do wlat2d(*,k)=wlat
;
; loop over days
;
for iday=0L,kday-1L do begin
;   erase
    if max(ubar(*,*,iday)) le 0. then goto,skipday
;   if setplot eq 'ps' then begin
;      set_plot,'ps'
;      xsize=nxdim/100.
;      ysize=nydim/100.
;      device,font_size=9
;      device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
;             /bold,/color,bits_per_pixel=8,/times,filename='meto_ubar_'+sdate(iday)+'.ps'
;      !p.charsize=1.25
;      !p.thick=2
;      !p.charthick=5
;      !y.thick=2
;      !x.thick=2
;   endif
;   !type=2^2+2^3
;   xmn=xorig(0)
;   xmx=xorig(0)+xlen
;   ymn=yorig(0)
;   ymx=yorig(0)+ylen
;   set_viewport,xmn,xmx,ymn,ymx
;   imn=long(strmid(sdate,4,2))
;   contour,ubar(*,*,iday),wlat,p,/ylog,xrange=[-90.,90.],yrange=[1000.,1.],title='Ubar '+sdate(iday),$
;           xtitle='Latitude',ytitle='Pressure (hPa)',charsize=2,charthick=2,xticks=6,/nodata,color=0
;   contour,ubar(*,*,iday),wlat,p,/ylog,/overplot,levels=10+10*findgen(20),color=0,thick=2
;   contour,ubar(*,*,iday),wlat,p,/ylog,/overplot,levels=-200+10*findgen(20),c_linestyle=5,color=0,thick=2
;   contour,ubar(*,*,iday),wlat,p,/ylog,/overplot,levels=[0],color=0,thick=6
;   oplot,[-65.,-65.],[50.,50.],psym=8,symsize=2,color=0
;   oplot,[-60.,-60.],[20.,20.],psym=8,symsize=2,color=0
;
; jet max
;
    index1=where(wlat2d le -40. and p2d le 100.)
    wlattemp=reform(wlat2d(index1))
    ptemp=reform(p2d(index1))
    u2d=reform(ubar(*,*,iday))
    ubartemp=reform(u2d(index1))
    index2=where(ubartemp eq max(ubartemp))
    if ptemp(index2(0)) ne 100. then begin
       oplot,wlattemp(index2),ptemp(index2),psym=8,symsize=3,color=mcolor*.9
       latjet(iday)=wlattemp(index2(0))
       pjet(iday)=ptemp(index2(0))
    endif
;
; pressure of zero wind line - SAVE AS A FUNCTION OF LATITUDE?
;
    pzero=-99.
    ilat=where(wlat gt -62. and wlat lt -57.)
    ubar1d=fltarr(nl)
    for k=0L,nl-1L do ubar1d(k)=mean(ubar(ilat,k,iday))
    index=where(p lt 300.,nl2)
    ptemp2=reform(p(index))
    ubartemp=reform(ubar1d(index))
    for k=nl2-1L,1,-1 do begin
        if ubartemp(k)*ubartemp(k-1) lt 0. then begin
           p1=ptemp2(k)
           p0=ptemp2(k-1)
           u1=ubartemp(k)
           u0=ubartemp(k-1)
           scale=abs(u0)/(abs(u0)+abs(u1))
           uzero=u0-scale*(abs(u1)-abs(u0))
           pzero=p0-scale*(p0-p1)
           goto,jumpout
        endif
    endfor
    jumpout:
;   oplot,[mean(wlat(ilat)),mean(wlat(ilat))],[pzero,pzero],psym=8,symsize=3,color=mcolor*.5
    print,sdate(iday),' Jet Max at ',wlattemp(index2(0)),ptemp(index2(0)),' Height zero wind ',pzero
    zzero(iday)=pzero
;   if setplot ne 'ps' then stop
;   if setplot eq 'ps' then begin
;      device, /close
;      spawn,'convert -trim meto_ubar_'+sdate(iday)+'.ps -rotate -90 meto_ubar_'+sdate(iday)+'.png'
;      spawn,'rm -f meto_ubar_'+sdate(iday)+'.ps'
;  endif
   skipday:
endfor
save,file='ukmo_12Z_SH_Ujet_lat_p_zero.sav',sdate,latjet,pjet,zzero
end
