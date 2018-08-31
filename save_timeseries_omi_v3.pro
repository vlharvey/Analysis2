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
nxdim=700
nydim=700
xorig=[0.25]
yorig=[0.35]
cbaryoff=0.08
cbarydel=0.01
xlen=0.6
ylen=0.4
set_plot,'ps'
setplot='ps'
read,'setplot= ',setplot
if setplot ne 'ps' then begin
   set_plot,'x'
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=255
endif
re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd
;
; loop over years
;
for iyear=2007,2014 do begin
;
; restore SOSST altitude grid
;
restore,'/Volumes/atmos/aura6/data/MLS_data/Datfiles_SOSST/o3_mls_v3.3_20040921.sav
;restore,'/aura3/data/SAGE_II_data/Datfiles_SOSST/altitude.sav'
nz=n_elements(altitude)
;
; enter dates to convert MLS pressure data
;
lstmn=9L & lstdy=1L & lstyr=iyear
ledmn=12L & leddy=31L & ledyr=iyear
lstday=0L & ledday=0L 
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
mino3time=fltarr(kday)
maxo3time=fltarr(kday)
meano3time=fltarr(kday)
areao3time=fltarr(kday)
minlattime=fltarr(kday)
sdates=strarr(kday)
dir='/Volumes/atmos/aura6/data/OMI_data/Datfiles/'
odir='/Volumes/atmos/aura6/data/OMI_data/Datfiles/'
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
mno=[31,28,31,30,31,30,31,31,30,31,30,31]
icount=0L
;
; loop over dates
;
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr
      z = stddat(imn,idy,iyr,ndays)
      if ndays ge ledday then goto,plotit
      syr=string(FORMAT='(I4.4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
      sdates(icount)=sdate
;
; look for EOS-OMI data files for today
;
      spawn,'ls '+dir+'OMI-Aura_L3-OMDOAO3e_'+syr+'m'+smn+sdy+'_*.he5',o3files
      result=size(o3files)
;
; this logic will jump day if missing
;
      if result(0) eq 0L then begin
         print,'OMI data missing on '+sdate
         goto,jumpday
      endif
;
; read data today
;
      readl3omi,o3files(0),longitude,latitude,o3,o3p
;
; check
;
;map_set,0,0,0,/contin,/grid,/noeras
;contour,o3,longitude,latitude,c_color=(findgen(20)/21.)*mcolor,nlevels=20,/noeras,/cell_fill,/overplot
;;contour,o3p,longitude,latitude,color=mcolor,levels=2.*findgen(20),/overplot,/foll
;stop
;
; on first day
;
      if icount eq 0L then begin
         rlat=-80.
;        print,lat
;        read,'Enter Latitude ',rlat
         index=where(latitude eq rlat)
         ilat=index(0)
         slat=strcompress(latitude(ilat),/remove_all)
         yindex=where(abs(latitude-rlat) le 5.)

         lon=0.*o3
         lat=0.*o3
         nc=n_elements(longitude)
         nr=n_elements(latitude)
         for i=0,nc-1 do lat(i,*)=latitude
         for j=0,nr-1 do lon(*,j)=longitude
         area=0.*o3
         deltax=longitude(1)-longitude(0)
         deltay=latitude(1)-latitude(0)
         for j=0,nr-1 do begin
             hy=re*deltay*dtr
             dx=re*cos(latitude(j)*dtr)*deltax*dtr
             area(*,j)=dx*hy    ; area of each grid point
         endfor
      endif
;
; retain column ozone at rlat
;
      dum=0.*o3
      o3data=reform(o3(*,yindex))	; 60S +/- 5
      index=where(finite(o3data) eq 1L)
      if index(0) ne -1L then mino3time(icount)=min(o3data(index))
      if index(0) ne -1L then maxo3time(icount)=max(o3data(index))
      if index(0) ne -1L then meano3time(icount)=mean(o3data(index))
;
; check latitude range
;
      o3slice=reform(o3(0,*))
      index=where(finite(o3slice) eq 1L)
      if index(0) ne -1L then minlattime(icount)=min(latitude(index))

map_set,-90,0,0,/ortho,/contin,/grid,/noeras,color=0
contour,o3,longitude,latitude,c_color=(findgen(20)/(21))*mcolor,/cell_fill,levels=100+20*findgen(20),/overplot
map_set,-90,0,0,/ortho,/contin,/grid,/noeras,color=0
;
; area where column ozone < 220 DU
; fill dum array with 1 if inside ozone hole. accommodate data void region over pole and noise inside the hole
;
      for i=0L,nc-1L do begin
          o3slice=reform(o3(i,*))
          index=where(latitude lt -40. and o3slice le 220.,nn)
          if index(0) ne -1L then dum(i,0:index(nn-1))=1.
      endfor
      index=where(dum eq 1.)
      if index(0) ne -1L then areao3time(icount)=total(area(index))
      print,mino3time(icount),' '+sdate+' ',min(latitude(index)),areao3time(icount)

index=where(dum eq 1.)
oplot,lon(index),lat(index),psym=4,color=0
stop
jumpday:
icount=icount+1L
goto,jump
;
; plot timeseries
;
plotit:
; postscript file
;
yearlab=strcompress(lstyr,/remove_all)
if lstyr ne ledyr then yearlab=strcompress(lstyr,/remove_all)+'-'+strcompress(ledyr,/remove_all)
if setplot eq 'ps' then begin
   lc=0
   xsize=nxdim/100.
   ysize=nydim/100.
   set_plot,'ps'
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='timeseries_omi_'+yearlab+'_'+slat+'.ps'
   !p.charsize=1.2
   !p.thick=2
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif

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

plot,findgen(kday),smooth(meano3time,7,/NaN),thick=5,xtickname=xlabs,xtickv=xindex,xticks=nxticks-1,ytitle='Minimum Column Ozone (DU) > '+slat,color=0,yrange=[0,500]
oplot,findgen(kday),smooth(mino3time,7,/NaN),thick=5,color=0,linestyle=5
oplot,findgen(kday),smooth(maxo3time,7,/NaN),thick=5,color=0,linestyle=5
axis,yrange=[-90,-30],yaxis=1,/save,ytitle='Minimum Latitude',color=150
oplot,findgen(kday),smooth(minlattime,7,/NaN),thick=3,color=150
xyouts,xmn+0.01,ymn+0.03,yearlab,/normal,charsize=2,charthick=5,color=0
;
; save yearly file
;
save,file='omi_mino3_'+yearlab+'_'+slat+'.sav',kday,mino3time,maxo3time,meano3time,minlattime,sdates,areao3time
;
; Close PostScript file and return control to X-windows
;
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim timeseries_omi_'+yearlab+'_'+slat+'.ps -rotate -90 timeseries_omi_'+yearlab+'_'+slat+'.jpg'
;  spawn,'rm -f timeseries_omi_'+yearlab+'_'+slat+'.ps'
endif

endfor  ; loop over years

end
