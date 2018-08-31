;****************************************************************************************
; Read OMI he5 data and save column ozone at a specified latitude
; Programed by: V. Lynn Harvey  10/1/14
;               CU/LASP									*
;****************************************************************************************
@stddat
@kgmt
@ckday
@kdate
@readl3omi

loadct,39
mcolor=byte(!p.color)
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
device,decompose=0
nxdim=1000
nydim=1000
xorig=[0.15]
yorig=[0.25]
cbaryoff=0.08
cbarydel=0.01
xlen=0.8
ylen=0.6
set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif
;
; version
;
sver='v2.2'
sver='v3.3'
;
; loop over years
;
col1=[10,50,100,140,190,200,210,250]
for iyear=2007,2014 do begin
    slat='-60.0000'
    yearlab=strcompress(iyear,/remove_all)
    restore,'omi_mino3_'+yearlab+'_'+slat+'.sav'	;,kday,mino3time,maxo3time,meano3time,minlattime,sdates
    if iyear eq 2007 then begin
       mino3time_all=mino3time
       maxo3time_all=maxo3time
       meano3time_all=meano3time
       minlattime_all=minlattime
       sdates_all=sdates

if setplot eq 'ps' then begin
   lc=0
   xsize=nxdim/100.
   ysize=nydim/100.
   set_plot,'ps'
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='timeseries_omi_allyears_'+slat+'.ps'
   !p.charsize=1.2
   !p.thick=2
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
nn=n_elements(sdates_all)
fdoy=fltarr(nn)
for i=0L,nn-1L do begin
    iyr=long(strmid(sdates_all(i),0,4))
    imn=long(strmid(sdates_all(i),4,2))
    idy=long(strmid(sdates_all(i),6,2))
    z = kgmt(imn,idy,iyr,iday)
    fdoy(i)=1.0*iday
endfor
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
xindex=where(strmid(sdates,6,2) eq '15',nxticks)
xlabs=strmid(sdates(xindex),4,2)
index=where(mino3time eq 0.)
if index(0) ne -1L then mino3time(index)=0./0.
index=where(maxo3time eq 0.)
if index(0) ne -1L then maxo3time(index)=0./0.
index=where(meano3time eq 0.)
if index(0) ne -1L then meano3time(index)=0./0.
index=where(minlattime ge 0.)
if index(0) ne -1L then minlattime(index)=0./0.
plot,findgen(kday),smooth(meano3time,7,/NaN),thick=5,xtickname=xlabs,xtickv=xindex,xticks=nxticks-1,ytitle='Minimum Column Ozone (DU) '+slat,color=0,yrange=[100,550]
    endif
    if iyear gt 2007 then begin
       mino3time_all=mino3time
       maxo3time_all=maxo3time
       meano3time_all=meano3time
       minlattime_all=minlattime
       sdates_all=sdates
    endif
index=where(mino3time eq 0.)
if index(0) ne -1L then mino3time(index)=0./0.
index=where(maxo3time eq 0.)
if index(0) ne -1L then maxo3time(index)=0./0.
index=where(meano3time eq 0.)
if index(0) ne -1L then meano3time(index)=0./0.
index=where(minlattime ge 0.)
if index(0) ne -1L then minlattime(index)=0./0.
oplot,findgen(kday),smooth(meano3time,7,/NaN),thick=5,color=col1(iyear-2007)
oplot,findgen(kday),smooth(mino3time,7,/NaN),thick=3,color=col1(iyear-2007)
oplot,findgen(kday),smooth(maxo3time,7,/NaN),thick=3,color=col1(iyear-2007)
endfor

loadct,0
axis,yrange=[-90,-30],yaxis=1,/save,ytitle='Minimum Latitude',color=150
oplot,findgen(kday),smooth(minlattime,7,/NaN),thick=3,color=150
loadct,39

!type=2^2+2^3+2^6
set_viewport,xmn,max(xorig)+xlen,ymn-cbaryoff,ymn-cbaryoff+cbarydel
imin=2007
imax=2014
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,xstyle=1,charsize=1.2,color=0,charthick=2
ybox=[0,10,10,0,0]
x1=imin
nlvls=n_elements(col1)
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
  xbox=[x1,x1,x1+dx,x1+dx,x1]
  polyfill,xbox,ybox,color=col1(j)
  x1=x1+dx
endfor
;
;
; Close PostScript file and return control to X-windows
;
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim timeseries_omi_allyears_'+slat+'.ps -rotate -90 timeseries_omi_allyears_'+slat+'.ps'
;  spawn,'rm -f timeseries_omi_allyears_'+slat+'.ps'
endif

end
