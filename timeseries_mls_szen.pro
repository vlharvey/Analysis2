;
; daily solar zenith angle coverage for MLS
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3_meto

a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
loadct,39
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
device,decompose=0
!p.background=icolmax
setplot='ps'
read,'setplot=',setplot
nxdim=750
nydim=750
xorig=[0.20,0.20]
yorig=[0.55,0.15]
xlen=0.7
ylen=0.3
cbaryoff=0.05
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smonth=['J','F','M','A','M','J','J','A','S','O','N','D']
mdir='/Volumes/earth/aura6/data/MLS_data/Datfiles_SOSST/'
lstmn=1
lstdy=1
lstyr=2006
ledmn=12
leddy=31
ledyr=2006
lstday=0
ledday=0
;goto,quick
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
kcount=0L

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
;
; restore MLS
;
;ALTITUDE        FLOAT     = Array[121]
;COMMENT         STRING    = Array[3]
;DATE            LONG      =     20070101
;FDOY            FLOAT     = Array[3494]
;ID              STRING    = Array[3494]
;LATITUDE        FLOAT     = Array[3494]
;LONGITUDE       FLOAT     = Array[3494]
;TIME            FLOAT     = Array[3494]
;
      dum=findfile(mdir+'cat_mls_v2.2_'+sdate+'.sav')
      if dum(0) eq '' then goto,skip
      restore,mdir+'cat_mls_v2.2_'+sdate+'.sav'
      print,sdate
      nlv=n_elements(altitude)
;
; compute solar zenith angle for each MLS profile
;
      doy=iday
      pi=3.14159265
      dtor=pi/180.
      earinc=23.5
      zangle=fltarr(n_elements(latitude))
      for ii=0L,n_elements(latitude)-1 do begin
          rlat=latitude(ii)
          rlon=longitude(ii)
          gmt=time(ii)
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
          zangle(ii) = chi/dtor
      endfor
;erase
index=where(zangle ne -99.)
;plot,zangle(index),psym=3,title=sdate,color=0,yrange=[0.,180.]
;print,sdate,' ',min(zangle(index)),mean(zangle(index)),max(zangle(index))
;
; declare time period arrays on first day
;
      if kcount eq 0L then begin
         ltime_mean=fltarr(kday)
         zen_mean=fltarr(kday)
         zen_min=fltarr(kday)
         zen_max=fltarr(kday)
         sdate_all=strarr(kday)
         kcount=1
      endif
      sdate_all(icount)=sdate
      zen_mean(icount)=mean(zangle(index))
      zen_min(icount)=min(zangle(index))
      zen_max(icount)=max(zangle(index))

skip:
      icount=icount+1L
goto,jump

plotit:

if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='timeseries_mls_szen.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
endif

syear=strmid(sdate_all,0,4)
smon=strmid(sdate_all,4,2)
sday=strmid(sdate_all,6,2)
xindex=where(sday eq '01',nxticks)
xlabs=smon(xindex)
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
nxticks=n_elements(xlabs)
plot,1.+findgen(kday),zen_mean,/noeras,xrange=[1.,kday],yrange=[0.,180.],$
      charsize=1.5,color=0,ytitle='Solar Zenith Angle',title='MLS Daily Min/Mean/Max',psym=8,$
      xticks=nxticks-1,xtickname=xlabs,xtickv=xindex
for i=0L,kday-1L do begin
    plots,i,zen_min(i)
    plots,i,zen_max(i),/continue,color=0
endfor
;
    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert -trim timeseries_mls_szen.ps -rotate -90 timeseries_mls_szen.jpg'
    endif
end
