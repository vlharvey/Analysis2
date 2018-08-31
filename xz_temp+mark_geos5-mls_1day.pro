;
; use MLS dmp file for GEOS-5 data interpolate to the measurement locations
;
; GEOS-5 and MLS longitude-altitude sections on 4 different days
; bin MLS temperature on lon/lat grid and look for geographical 
; pattern to T differences with G5
; G5 temperature and vortex edge
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3_meto

ipan=0
npp=1
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
loadct,39
device,decompose=0
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
setplot='ps'
read,'setplot=',setplot
nxdim=750
nydim=750
xorig=[0.15,0.4,0.65]
yorig=[0.4,.4,.4]
xlen=.2
ylen=.2
cbaryoff=0.06
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
dir='/aura7/harvey/GEOS5_data/Datfiles/'
sdir='/aura6/data/SABER_data/Datfiles/'
mdir='/aura6/data/MLS_data/Datfiles_SOSST/'
stimes=[$
'_AVG.V01.']
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
lstmn=2
lstdy=21
lstyr=2005
ledmn=2
leddy=21
ledyr=2005
lstday=0
ledday=0
nlat=35L
latbin=-85+5.*findgen(nlat)
dy=latbin(1)-latbin(0)
nlon=12L
lonbin=15.+30.*findgen(nlon)
dx=lonbin(1)-lonbin(0)
rlat=70.
;print,latbin
;read,' Enter latitude ',rlat
index=where(rlat eq latbin)
ilat=index(0)
slat=strcompress(long(rlat),/remove_all)
;
; Ask interactive questions- get starting/ending date and p surface
;
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
;read,' Enter ending date   (month, day, year) ',ledmn,leddy,ledyr
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1991 then stop,'Year out of range '
if ledyr lt 1991 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
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
      if ndays gt ledday then stop
      if iyr ge 2000L then iyr1=iyr-2000L
      if iyr lt 2000L then iyr1=iyr-1900L
;
;***Read GEOS-5 data
;
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
      print,sdate
      if setplot eq 'ps' then begin
         set_plot,'ps'
         xsize=nxdim/100.
         ysize=nydim/100.
         !psym=0
         !p.font=0
         device,font_size=9
         device,/landscape,bits=8,filename='xz_temp+mark_geos5-mls_'+slat+'_'+sdate+'.ps'
         device,/color
         device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
                xsize=xsize,ysize=ysize
      endif
      ifile='DAS.ops.asm.tavg3d_dyn_v.GEOS510.MetO.'+sdate+stimes(0)+'nc3'
;
; read GEOS-5 data
;
      rd_geos5_nc3_meto,dir+ifile,nc,nr,nthg,alon,alat,th,$
                  pv2,p2,msf2,u2,v2,q2,qdf2,mark2,sf2,vp2,iflag
      if iflag eq 1 then goto,jump
      x=fltarr(nc+1)
      x(0:nc-1)=alon(0:nc-1)
      x(nc)=alon(0)+360.
      t2=0.*pv2
      for k=0,nthg-1 do t2(*,*,k)=th(k)*((p2(*,*,k)/1000.)^(.286))
      z2=(msf2-1004.*t2)/(9.86*1000.)
      index=where(mark2 lt -1.)
      if index(0) ne -1L then mark2(index)=-1.*(abs(mark2(index))/abs(mark2(index)))
;
; restore MLS dmp_mls_v2.2.geos5.20080614.sav file
;
; ALTITUDE        FLOAT     = Array[121]
; COMMENT         STRING    = Array[7]
; DENTOT          FLOAT     = Array[3494, 121]
; ID              STRING    = Array[3494]
; PRESSURE        FLOAT     = Array[3494, 121]
; TEMPERATURE     FLOAT     = Array[3494, 121]
; TEMPERATURE_ERROR
; TEMPERATURE_MASK
;
      dum=findfile(mdir+'tpd_mls_v2.2_'+sdate+'.sav')
      if dum(0) eq '' then goto,jump
      restore,mdir+'tpd_mls_v2.2_'+sdate+'.sav'
      restore,mdir+'cat_mls_v2.2_'+sdate+'.sav'   ; latitude
      restore,mdir+'dmps_mls_v2.2.geos5.'+sdate+'.sav'
      nth=n_elements(thlev)
      nprof=n_elements(time)
      nlv=n_elements(altitude)
;
; interpolate GEOS marker to altitude surfaces
;
      mark2z=fltarr(nr,nc,nlv)
      for kk=20L,80L do begin
      zz=altitude(kk)
      if max(z2) lt zz then goto,jumplev
      for j=0L,nr-1L do begin
      for i=0L,nc-1L do begin
          for k=1L,nthg-1L do begin
              zup=z2(j,i,k-1) & zlw=z2(j,i,k)
              if zup ne 0. and zlw ne 0. then begin
              if zup ge zz and zlw le zz then begin
                 zscale=(zup-zz)/(zup-zlw)
                 mark2z(j,i,kk)=mark2(j,i,k-1)+zscale*(mark2(j,i,k)-mark2(j,i,k-1))
              endif
              endif
          endfor
      endfor
      endfor
      jumplev:
      endfor
;
; bin MLS and GEOS temperature by lonbin/latbin. interp GEOS to altitudes.
;
      tgeos3d=fltarr(nlon,nlat,nlv)
      ngeos3d=lonarr(nlon,nlat,nlv)
      mls3d=fltarr(nlon,nlat,nlv)
      nmls3d=lonarr(nlon,nlat,nlv)
      nprof=n_elements(fdoy)
      for iprof=0L,nprof-1L do begin
          zprof=reform(z_prof(iprof,*))
          tprof=reform(tp_prof(iprof,*))
          xlon=longitude(iprof)
          if xlon gt lonbin(nlon-1)+dx then xlon=xlon-360.
          xlat=latitude(iprof)
          if latitude(iprof) lt -90. then goto,skipprof
          for i=0L,nlon-1L do begin
              if xlon ge lonbin(i)-dx/2. and xlon lt lonbin(i)+dx/2. then begin
              for j=ilat,ilat do begin
                  if xlat ge latbin(j)-dy/2. and xlat lt latbin(j)+dy/2. then begin
                  for kk=20L,80 do begin
                      zz=altitude(kk)
                      mls3d(i,j,kk)=mls3d(i,j,kk)+TEMPERATURE(iprof,kk)
                      nmls3d(i,j,kk)=nmls3d(i,j,kk)+1L
                      for k=1L,nth-1L do begin
                          zup=zprof(k-1) & zlw=zprof(k)        ; profiles are top down
                          if zup ge zz and zlw le zz then begin
                             zscale=(zup-zz)/(zup-zlw)
                             tgeos3d(i,j,kk)=tgeos3d(i,j,kk)+(tprof(k-1)+zscale*(tprof(k)-tprof(k-1)))
                             ngeos3d(i,j,kk)=ngeos3d(i,j,kk)+1L
                          endif
                      endfor
                  endfor
                  endif
              endfor
              endif
          endfor
          skipprof:
      endfor
      index=where(nmls3d gt 0L)
      if index(0) ne -1L then mls3d(index)=mls3d(index)/float(nmls3d(index))
      index=where(ngeos3d gt 0L)
      if index(0) ne -1L then tgeos3d(index)=tgeos3d(index)/float(ngeos3d(index))
;
; interpolate mark2z to lonbin,latbin at rlat
;
      geosmark1=fltarr(nlon,nlv)
      nmark1=lonarr(nlon,nlv)
      for i=0L,nlon-1L do begin
          for ii=0L,nc-1L do begin
              xlon=alon(ii)
              if xlon gt lonbin(nlon-1)+dx then xlon=xlon-360.
              if xlon ge lonbin(i)-dx/2. and xlon lt lonbin(i)+dx/2. then begin
              for jj=0,nr-1L do begin
                  xlat=alat(jj)
                  if xlat ge latbin(ilat)-dy/2. and xlat lt latbin(ilat)+dy/2. then begin
                     geosmark1(i,*)=geosmark1(i,kk)+mark2z(jj,ii,*)
                     nmark1(i,kk)=nmark1(i,kk)+1L
                  endif
              endfor
              endif
          endfor
      endfor
      index=where(nmark1 gt 0.)
      if index(0) ne -1L then geosmark1(index)=geosmark1(index)/float(nmark1(index))
;
; extract 2d arrays of GEOS and MLS temperatures at ilat
;
      mlst1=reform(mls3d(*,ilat,*))
      geost1=reform(tgeos3d(*,ilat,*))

erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
!psym=0
xyouts,0.45,yorig(0)+ylen+0.05,sdate,charsize=2,charthick=2,/normal,color=0
nlvls=10
col1=1+indgen(nlvls)*icolmax/nlvls
level=180.+10.*findgen(nlvls)
contour,geost1,lonbin,altitude,levels=level,c_color=col1,thick=1,/cell_fill,/noeras,yrange=[20,80],$
        ytitle='Altitude (km)',xtitle='Longitude',title='GEOS-5',color=0
contour,geost1,lonbin,altitude,/overplot,levels=level,/follow,color=0,c_labels=0*level
contour,geosmark1,lonbin,altitude,/overplot,levels=[0.1],/follow,color=0,c_labels=0*level,thick=5
loadct,0
contour,geosmark1,lonbin,altitude,/overplot,levels=[-0.1],/follow,color=125,c_labels=0*level,thick=5
loadct,39
imin=min(level)
imax=max(level)
ymnb=ymn -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],xtitle=slat+' N Temperature (K)',color=0
ybox=[0,10,10,0,0]
x1=imin
dxx=(imax-imin)/float(nlvls)
for jj=0,nlvls-1 do begin
xbox=[x1,x1,x1+dxx,x1+dxx,x1]
polyfill,xbox,ybox,color=col1(jj)
x1=x1+dxx
endfor

xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
!psym=0
contour,mlst1,lonbin,altitude,levels=level,c_color=col1,thick=1,/cell_fill,/noeras,yrange=[20,80],xtitle='Longitude',$
        ytickname=[' ',' ',' ',' ',' ',' '],title='MLS',color=0
contour,mlst1,lonbin,altitude,/overplot,levels=level,/follow,color=0,c_labels=0*level
contour,geosmark1,lonbin,altitude,/overplot,levels=[0.1],/follow,color=0,c_labels=0*level,thick=5
loadct,0
contour,geosmark1,lonbin,altitude,/overplot,levels=[-0.1],/follow,color=125,c_labels=0*level,thick=5
loadct,39
imin=min(level)
imax=max(level)
ymnb=ymn -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],xtitle=slat+' N Temperature (K)',color=0
ybox=[0,10,10,0,0]
x1=imin
dxx=(imax-imin)/float(nlvls)
for jj=0,nlvls-1 do begin
xbox=[x1,x1,x1+dxx,x1+dxx,x1]
polyfill,xbox,ybox,color=col1(jj)
x1=x1+dxx
endfor
;
; superimpose pdiff
;
restore,'c11_rb.tbl'
tvlct,c1,c2,c3
col2=1+indgen(11)
pdiff=-99.+0.*mlst1
index=where(mlst1 gt 0. and geost1 gt 0.)
if index(0) eq -1L then goto,jump
pdiff(index)=geost1(index)-mlst1(index)
level=-25.+5.*findgen(11)
!type=2^2+2^3
xmn=xorig(2)
xmx=xorig(2)+xlen
ymn=yorig(2)
ymx=yorig(2)+ylen
set_viewport,xmn,xmx,ymn,ymx
contour,pdiff,lonbin,altitude,levels=level,/cell_fill,c_color=col2,min_value=-99.,yrange=[20,80],$
        xtitle='Longitude',ytickname=[' ',' ',' ',' ',' ',' '],title='GEOS-5 - MLS',color=0
contour,pdiff,lonbin,altitude,/overplot,levels=level,color=0,/follow,min_value=-99.,$
        c_labels=0*level
contour,pdiff,lonbin,altitude,/overplot,levels=[0],color=0,thick=1,min_value=-99.
contour,geosmark1,lonbin,altitude,/overplot,levels=[0.1],/follow,color=0,c_labels=0*level,thick=5
loadct,0
contour,geosmark1,lonbin,altitude,/overplot,levels=[-0.1],/follow,color=125,c_labels=0*level,thick=5
restore,'c11_rb.tbl'
tvlct,c1,c2,c3
imin=min(level)
imax=max(level)
ymnb=ymn -cbaryoff
ymxb=ymnb+cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,$
      xtitle=slat+' N Differences (K)',xticks=n_elements(level)/2,color=0
ybox=[0,10,10,0,0]
x2=imin
dxx=(imax-imin)/(float(n_elements(col2)))
for jj=0L,n_elements(col2)-1 do begin
    xbox=[x2,x2,x2+dxx,x2+dxx,x2]
    polyfill,xbox,ybox,color=col2(jj)
    x2=x2+dxx
endfor
loadct,39

   if setplot ne 'ps' then stop
   if setplot eq 'ps' then begin
      device, /close
      spawn,'convert -trim xz_temp+mark_geos5-mls_'+slat+'_'+sdate+'.ps '+$
            ' -rotate -90  xz_temp+mark_geos5-mls_'+slat+'_'+sdate+'.jpg'
   endif
goto,jump
end
