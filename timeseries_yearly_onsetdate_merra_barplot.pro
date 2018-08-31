;
; MERRA
; timeseries of vortex area
; Sep barplot as in Burnett
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
xorig=[0.25]
yorig=[0.3]
cbaryoff=0.07
cbarydel=0.01
xlen=0.6
ylen=.5
device,decompose=0
mcolor=byte(!p.color)
nlvls=20
col1=1+(indgen(nlvls)/float(nlvls))*mcolor
PI2=6.2831853071796
DTR=PI2/360.
RADEA=6.37E6
dum=1979+indgen(36)	; 1979-2014
syear=strcompress(dum,/remove_all)     ; 1979-2014
lyear=long(syear)
nyear=n_elements(syear)
nlvls=nyear
ssw=0*indgen(nyear)
index=where(syear eq '1985' or syear eq '1987' or syear eq '2004' or syear eq '2006' or syear eq '2009' or syear eq '2013')	; major SSW 
ssw(index)=1
smon=['09','10','11','12','01','02','03','04']
smon=['09']
febflux=[204,200,205,214,123,141,74,84,72,105,222,178,243,232,143,100,86,72,74,93,143,173,147,205,125,107,97,77,78,71,70,85,95,107,104,170]
kmon=[30,31,30,31,31,28,31,30]
smonth='     '+['Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr',' ']
smonth='     '+['Oct','Nov','Dec','Jan','Feb','Mar',' ']
smonth='          '+['Nov','Dec','Jan','Feb','Mar',' ']
nmon=n_elements(smon)
set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif
;
; loop over years
;
for iyear=0L,nyear-1L do begin
for imon=0,0 do begin	;0L,nmon-1L do begin ; January
    icount=0
    dum=findfile('vortex_shape_merra_'+syear(iyear)+smon(imon)+'.sav')
    if dum(0) eq '' then goto,skipmon
    restore,dum(0)
    print,dum(0)
    if icount ne 0L then begin
       area1_all=[area1_all,area1]
       centroid_longitude1_all=[centroid_longitude1_all,centroid_longitude1]
       centroid_latitude1_all=[centroid_latitude1_all,centroid_latitude1]
       number_vortex_lobes1_all=[number_vortex_lobes1_all,number_vortex_lobes1]
       ellipticity1_all=[ellipticity1_all,ellipticity1]
       altitude_all=[altitude_all,altitude]
       sdate_tot=[sdate_tot,sdate_all]
    endif
    if icount eq 0L then begin
       if iyear eq 0L then begin
          print,th
          rth=2000.
;         rth=1000.
;read,'Enter desired theta ',rth
          index=where(th eq rth)
          ith=index(0)
          sth=strcompress(long(th(ith)),/remove_all)+'K'
          erase
          if setplot eq 'ps' then begin
             set_plot,'ps'
             xsize=nxdim/100.
             ysize=nydim/100.
             !psym=0
             !p.font=0
             device,font_size=9
             device,/landscape,bits=8,filename='barplot_area_merra_'+sth+'_sep.ps'
             device,/color
             device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
                    xsize=xsize,ysize=ysize
             !p.charsize=2.0
          endif
       endif

       area1_all=area1
       centroid_longitude1_all=centroid_longitude1
       centroid_latitude1_all=centroid_latitude1
       number_vortex_lobes1_all=number_vortex_lobes1
       ellipticity1_all=ellipticity1
       altitude_all=altitude
       sdate_tot=sdate_all
       icount=1L
    endif

skipmon:
endfor  ; loop over months

if iyear eq 0L then begin
   ztd_allyears=fltarr(nyear)
   area_allyears=fltarr(nyear,nth)
   areasig_allyears=fltarr(nyear)
   ellip_allyears=fltarr(nyear)
   y0_allyears=fltarr(nyear)
   x0_allyears=fltarr(nyear)
   nvort_allyears=fltarr(nyear)
endif
index=where(area1_all lt 0.)
if index(0) ne -1L then area1_all(index)=0./0.
for ith=0,nth-1L do begin
    area1_lev=reform(area1_all(*,ith))
    index=where(finite(area1_lev) eq 1 and area1_lev ge 5)
    area_allyears(iyear,ith)=index(0)
endfor

endfor	; loop over years
;index=where(lyear lt max(lyear),nyear)
;lyear=lyear(index)
;syear=syear(index)
;area_allyears=area_allyears(index)
;
; calculate means and sigmas for non ssw years
;
arrmean=median(area_allyears)
arrsig=stdev(area_allyears)

xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
;
; plot all years
;
;b = BARPLOT(area_allyears-arrmean, NBARS=nyear, FILL_COLOR=0,xtickname=syear)

index=where(area_allyears eq 0.)
if index(0) ne -1L then area_allyears(index)=0./0.
area_allyears=smooth(area_allyears,3,/nan,/edge_truncate)
kindex=where(th le 2600. and th ge 900.,nth)
area_allyears=reform(area_allyears(*,kindex))
th=reform(th(kindex))
plot,lyear,area_allyears(*,0),yrange=[1,30],color=0,charsize=2,charthick=2,ytitle='Vortex Onset Day (Sep)',xrange=[min(lyear),max(lyear)],psym=8,/noeras,symsize=1.5,/nodata
for ith=0L,nth-1L do begin
    oplot,lyear,area_allyears(*,ith),thick=5,color=((nth-float(ith))/float(nth))*mcolor
print,th(ith)
endfor

xyouts,xmn+0.01,ymx-0.03,'September',/normal,charsize=2,charthick=2,color=0
;
nlvls=nth
col1=1+indgen(nlvls)*mcolor/nlvls

imin=long(min(th))
imax=long(max(th))
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,charsize=1.5,charthick=2,/noeras
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for jj=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(jj)
x1=x1+dx
endfor


if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device,/close
   spawn,'convert -trim barplot_area_merra_'+sth+'_sep.ps -rotate -90 barplot_area_merra_'+sth+'_sep.jpg'
endif

end
