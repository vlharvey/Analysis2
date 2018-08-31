;
; save monthly vortex shape diagnostics
; MERRA version
;
@stddat
@kgmt
@ckday
@kdate
@rd_merra_nc3
@vortexshape

loadct,39
mcolor=byte(!p.color)
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
nxdim=700
nydim=700
xorig=[0.20]
yorig=[0.3]
cbaryoff=0.1
cbarydel=0.01
xlen=0.6
ylen=0.4
device,decompose=0
!NOERAS=-1
nlvls=20L
col1=1+(indgen(nlvls)/float(nlvls))*mcolor
PI2=6.2831853071796
DTR=PI2/360.
RADEA=6.37E6
re=40000./2./!pi
earth_area=4.*!pi*re*re
hem_area=earth_area/2.0
rtd=double(180./!pi)
dtr=1./rtd
;syear=['2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014']
dum=1979+indgen(36)
dum=2008
syear=strcompress(dum,/remove_all)
nyear=n_elements(syear)
smon=['01','02','03','04','05','06','07','08','09','10','11','12']
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
; get file listing
;
dir='/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_theta_'

for iyear=0L,nyear-1L do begin

lstmn=1
lstdy=1
lstyr=long(syear(iyear))
ledmn=2
leddy=1
ledyr=long(syear(iyear))
lstday=0
ledday=0
if lstyr lt 79 then lstyr=lstyr+2000
if ledyr lt 79 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1979 then stop,'Year out of range '
if ledyr lt 1979 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
sdate_all=strarr(kday)
dayno=lonarr(kday)
;
; Compute initial Julian date
;
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr
; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,plotit
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy

    print,sdate
    sdate_all(icount)=sdate
    dayno(icount)=iday
    dum=findfile(dir+sdate+'.nc3')
    if dum ne '' then ncfile0=dir+sdate+'.nc3'
    rd_merra_nc3,ncfile0,nc,nr,nth,alon,alat,th,pv2,p2,$
       u2,v2,qdf2,mark2,qv2,z2,sf2,q2,iflag
    if iflag ne 0L then goto,jumpstep
    tmp2=0.*p2
    for k=0L,nth-1L do tmp2(*,*,k)=th(k)*(p2(*,*,k)/1000.)^0.286
;
; monthly files
;
    if icount eq 0L then begin
       rlat=60.
;      print,alat
;      read,'Enter desired latitude ',rlat
       index=where(abs(rlat-alat) eq min(abs(rlat-alat)))
       ilat=index(0)
       slat=strcompress(long(alat(ilat)),/remove_all)
       x=fltarr(nc+1)
       x(0:nc-1)=alon(0:nc-1)
       x(nc)=alon(0)+360.
       x2d=fltarr(nc+1,nr)
       y2d=fltarr(nc+1,nr)
       for i=0,nc do y2d(i,*)=alat
       for j=0,nr-1 do x2d(*,j)=x

         lon=fltarr(nc,nr)
         lat=fltarr(nc,nr)
         for i=0,nc-1 do lat(i,*)=alat
         for j=0,nr-1 do lon(*,j)=alon
         area=0.*lat
         deltax=alon(1)-alon(0)
         deltay=alat(1)-alat(0)
         for j=0,nr-1 do begin
             hy=re*deltay*dtr
             dx=re*cos(alat(j)*dtr)*deltax*dtr
             area(*,j)=dx*hy    ; area of each grid point
         endfor

       ztd=fltarr(kday,nth)
       pressure=fltarr(kday,nth)
       area1=fltarr(kday,nth)
       centroid_longitude1=fltarr(kday,nth)
       centroid_latitude1=fltarr(kday,nth)
       number_vortex_lobes1=fltarr(kday,nth)
       ellipticity1=fltarr(kday,nth)
       altitude=fltarr(kday,nth)
       ubar=fltarr(kday,nth)
    endif
    ubar(icount,*)=reform(mean(u2(ilat,*,*),dim=2))
;
; find vortex centroids, ellipticity_profile
;
      marker_USLM = make_array(nc,nr,nth)
      for k=0,nth-1 do marker_USLM(*,*,k) = transpose(mark2(*,*,k))
      shape = vortexshape(marker_USLM, alat, alon)
      centroid=shape.nhcentroid
      centroidx=reform(centroid(0,*))
      centroidy=reform(centroid(1,*))
      axis=shape.axis
      majoraxis=reform(axis(0,*))
      minoraxis=reform(axis(1,*))
      ellipticity_profile=minoraxis/majoraxis
      index=where(centroidx lt 0.)
      if index(0) ne -1L then centroidx(index)=centroidx(index)+360.
      centroid_longitude1(icount,*)=centroidx
      centroid_latitude1(icount,*)=centroidy
      ellipticity1(icount,*)=ellipticity_profile
;
; compute vortex area based on SF and CO
;
      for k = 0L, nth - 1L do begin
          mark1=transpose(mark2(*,*,k))
          index=where(lat gt 0. and mark1 eq 1.0,nn)
          if index(0) ne -1L then area1(icount,k)=100.*total(area(index))/hem_area
          z1=transpose(z2(ilat,*,k))
          p1=transpose(p2(ilat,*,k))
          if index(0) eq -1L then index=where(z1 ne 0.)
          if index(0) ne -1L then begin
             altitude(icount,k)=mean(z1(index))
             pressure(icount,k)=mean(p1(index))
          endif
      endfor
icount=icount+1L
jumpstep:
goto,jump

plotit:
if setplot eq 'ps' then begin
   lc=0
   xsize=nxdim/100.
   ysize=nydim/100.
   set_plot,'ps'
   device,/color,/landscape,bits=8,filename='zt_merra_ubar_'+syear(iyear)+'_'+slat+'.ps'
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
endif

erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
imin=-20.
imax=80.
int=5.
nlvls=21
col1=1+indgen(nlvls)*mcolor/nlvls
level=imin+int*findgen(nlvls)
index=where(ubar eq 0. or pressure eq 0.)
if index(0) ne -1L then ubar(index)=0./0.
if index(0) ne -1L then pressure(index)=0./0.
contour,ubar,dayno,pressure,/ylog,levels=level,/cell_fill,c_color=col1,color=0,title=syear(iyear),$
            ytitle='Pressure (hPa)',yrange=[100.,0.01],charsize=1.5,charthick=2
contour,ubar,dayno,pressure,/ylog,levels=10.+10.*findgen(10),/overplot,/follow,color=0,thick=2
contour,ubar,dayno,pressure,/ylog,levels=-100.+10.*findgen(10),/overplot,/follow,color=mcolor,c_linestyle=5,thick=2
;contour,smooth(ubar,3,/Nan),dayno,pressure,/ylog,levels=0,/overplot,/follow,color=0,thick=5
contour,smooth(area1,3),dayno,pressure,/ylog,levels=5+5.*findgen(7),/overplot,/follow,c_color=1+indgen(7)*mcolor/7.,thick=5

imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='MERRA Ubar at '+slat+' N (m/s)',/noeras,charsize=1.5,charthick=2
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for j=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(j)
x1=x1+dx
endfor
;save,file='merra_ubar_'+syear(iyear)+'_'+slat+'.sav',ubar,dayno,pressure,altitude,area1

; Close PostScript file and return control to X-windows
     if setplot ne 'ps' then stop
     if setplot eq 'ps' then begin
        device, /close
        spawn,'convert -trim zt_merra_ubar_'+syear(iyear)+'_'+slat+'.ps -rotate -90 '+$
                            'zt_merra_ubar_'+syear(iyear)+'_'+slat+'.jpg'
     endif

endfor	; loop over years
end
