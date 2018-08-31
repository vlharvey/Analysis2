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
@rd_merra_nc3

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
xorig=[0.15]
yorig=[0.2]
xlen=.8
ylen=.6
cbaryoff=0.06
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
mdir='/atmos/aura6/data/MLS_data/Datfiles_SOSST/'
dir='/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_theta_'
sdir='/atmos/aura6/data/SABER_data/Datfiles/'
stimes=[$
'_AVG.V01.']
month=['Jan','Feb','Mar','Apr','May','Jun',$
       'Jul','Aug','Sep','Oct','Nov','Dec']
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
lstmn=6
lstdy=1
lstyr=2009
ledmn=8
leddy=1
ledyr=2009
lstday=0
ledday=0
nlat=35L
latbin=-85+5.*findgen(nlat)
dy=latbin(1)-latbin(0)
nlon=12L
lonbin=15.+30.*findgen(nlon)
dx=lonbin(1)-lonbin(0)
rlat=-60.
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
         device,/landscape,bits=8,filename='xz_temp+mark_merra-mls_'+slat+'_'+sdate+'.ps'
         device,/color
         device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
                xsize=xsize,ysize=ysize
      endif
;
; read met data
;
    ncfile0=dir+sdate+'.nc3'
    rd_merra_nc3,ncfile0,nc,nr,nth,alon,alat,th,pv2,p2,$
       u2,v2,qdf2,mark2,qv2,z2,sf2,q2,iflag
    if iflag eq 1 then goto,jump
    x=fltarr(nc+1)
    x(0:nc-1)=alon(0:nc-1)
    x(nc)=alon(0)+360.
    t2=0.*pv2
    for k=0,nth-1 do t2(*,*,k)=th(k)*((p2(*,*,k)/1000.)^(.286))
    index=where(mark2 lt -1.)
    if index(0) ne -1L then mark2(index)=-1.*(abs(mark2(index))/abs(mark2(index)))
;
; restore MLS on this day
;
; ALTITUDE        FLOAT     = Array[121]
; COMMENT         STRING    = Array[4]
; DATE            LONG      =     20070101
; ERR             FLOAT     = Array[3491, 121]
; FDOY            FLOAT     = Array[3491]
; ID              STRING    = Array[3491]
; LATITUDE        FLOAT     = Array[3491]
; LONGITUDE       FLOAT     = Array[3491]
; MASK            FLOAT     = Array[3491, 121]
; MIX             FLOAT     = Array[3491, 121]
; TIME            FLOAT     = Array[3491]
;
      dum=findfile(mdir+'cat_mls_v3.3_'+sdate+'.sav')
      if dum(0) eq '' then goto,jump
      restore,mdir+'cat_mls_v3.3_'+sdate+'.sav'
      restore,mdir+'h2o_mls_v3.3_'+sdate+'.sav'
      restore,mdir+'tpd_mls_v3.3_'+sdate+'.sav'
      print,sdate
      nth=n_elements(thlev)
      nprof=n_elements(time)
      nlv=n_elements(altitude)
;
; apply mask
;
      index=where(mask eq -99.)
      if index(0) ne -1L then mix(index)=-99.
      mlsh2o=mix
;
; bin MLS XZ at rlat
;
      h2oxz=fltarr(nlon,nlv)
      nh2oxz=lonarr(nlon,nlv)
      nprof=n_elements(fdoy)
      for iprof=0L,nprof-1L do begin
          if abs(latitude(iprof)-rlat) le 2.5 then begin

          xlon=longitude(iprof)
          if xlon gt lonbin(nlon-1)+dx then xlon=xlon-360.
          for i=0L,nlon-1L do begin
              if xlon ge lonbin(i)-dx/2. and xlon lt lonbin(i)+dx/2. then begin
                 h2oprof=reform(mlsh2o(iprof,*))
                 index=where(h2oprof ne -99.)
                 if index(0) ne -1L then begin
                    h2oxz(i,index)=h2oxz(i,index)+h2oprof(index)*1.e6
                    nh2oxz(i,index)=nh2oxz(i,index)+1.
                 endif

              endif
          endfor
          endif
          skipprof:
      endfor
      index=where(nh2oxz gt 0L)
      if index(0) ne -1L then h2oxz(index)=h2oxz(index)/float(nh2oxz(index))
;
; zonal mean anomaly
;
anomxz=0.*h2oxz
dum=mean(h2oxz,dim=1)
for i=0L,nlon-1L do anomxz(i,*)=h2oxz(i,*)-dum
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
!psym=0
nlvls=20
col1=1+indgen(nlvls)*icolmax/nlvls
level=0.5*findgen(nlvls)
contour,h2oxz,lonbin,altitude,levels=level,c_color=col1,thick=1,/cell_fill,/noeras,yrange=[60,90],$
        ytitle='Altitude (km)',xtitle='Longitude',title=sdate,color=0
contour,anomxz,lonbin,altitude,/overplot,levels=-10+0.5*findgen(20),/follow,color=mcolor,c_labels=0*findgen(20),c_linestyle=5
contour,anomxz,lonbin,altitude,/overplot,levels=0.5+0.5*findgen(20),/follow,color=0,c_labels=0*findgen(20)
imin=min(level)
imax=max(level)
ymnb=ymn -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],xtitle=slat+' N MLS H2O (ppmv)',color=0
ybox=[0,10,10,0,0]
x1=imin
dxx=(imax-imin)/float(nlvls)
for jj=0,nlvls-1 do begin
xbox=[x1,x1,x1+dxx,x1+dxx,x1]
polyfill,xbox,ybox,color=col1(jj)
x1=x1+dxx
endfor

   if setplot ne 'ps' then stop
   if setplot eq 'ps' then begin
      device, /close
      spawn,'convert -trim xz_temp+mark_merra-mls_'+slat+'_'+sdate+'.ps '+$
            ' -rotate -90  xz_temp+mark_merra-mls_'+slat+'_'+sdate+'.jpg'
   endif
goto,jump
end
