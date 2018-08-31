;
; zonal mean tbar is from GEOS5 DMP files
; (construct from data interpolated to hirdls locations)
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3_meto

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
xorig=[0.15]
yorig=[0.25]
xlen=0.7
ylen=0.5
cbaryoff=0.08
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smonth=['J','F','M','A','M','J','J','A','S','O','N','D']
sdir='/aura6/data/MLS_data/Datfiles_SOSST/'
lstmn=1
lstdy=15
lstyr=2005
ledmn=3
leddy=15
ledyr=2005
lstday=0
ledday=0
;goto,plotit
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
nlat=35L
latbin=-85+5.*findgen(nlat)
dy=latbin(1)-latbin(0)

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then goto,saveit
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
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
    dum=findfile(sdir+'dmps_mls_v2.2.geos5.'+sdate+'.sav')
    if dum(0) eq '' then goto,skip
    restore,sdir+'dmps_mls_v2.2.geos5.'+sdate+'.sav'
    restore,sdir+'cat_mls_v2.2_'+sdate+'.sav'	; latitude
    print,sdir+'cat_mls_v2.2_'+sdate+'.sav'
    nth=n_elements(thlev)
;
; compute zonal mean T, Z
;
    tbar=fltarr(nlat,nth)
    zbar=fltarr(nlat,nth)
    nbar=lonarr(nlat,nth)
    for ii=0L,n_elements(id)-1L do begin
        ymean=latitude(ii)
        for j=0L,nlat-1L do begin
            if ymean ge latbin(j)-dy/2. and ymean lt latbin(j)+dy/2. then begin
               tbar(j,*)=tbar(j,*)+tp_prof(ii,*)
               zbar(j,*)=zbar(j,*)+z_prof(ii,*)
               nbar(j,*)=nbar(j,*)+1L
            endif
        endfor
    endfor
    good=where(nbar gt 0L)
    if good(0) ne -1L then tbar(good)=tbar(good)/float(nbar(good))
    if good(0) ne -1L then zbar(good)=zbar(good)/float(nbar(good))
;
; interpolate GEOS temperature to height surfaces
;
nlv=n_elements(altitude)
tbarz=fltarr(nlat,nlv)
for kk=0L,nlv-1L do begin
    zz=altitude(kk)
    for j=0L,nlat-1L do begin
        if max(zbar(j,*)) eq 0. then goto,jumplat
        for k=1L,nth-1L do begin
            zup=zbar(j,k-1) & zlw=zbar(j,k)
            if zup ne 0. and zlw ne 0. then begin
            if zup ge zz and zlw le zz then begin
               zscale=(zup-zz)/(zup-zlw)
               tbarz(j,kk)=tbar(j,k-1)+zscale*(tbar(j,k)-tbar(j,k-1))
;print,zlw,zz,zup,zscale
;print,tbar(j,k),tbarz(j,kk),tbar(j,k-1)
;stop
            endif
            endif
         endfor
jumplat:
      endfor
jumplev:
endfor

;nlvls=26
;col1=1+indgen(nlvls)*icolmax/nlvls
;erase
;!type=2^2+2^3
;set_viewport,.15,.85,.55,.9
;contour,tbar,latbin,thlev,c_color=col1,color=0,/fill,levels=170.+5.*findgen(nlvls),/noerase,$
;        title='GEOS-5 at MLS '+sdate,yrange=[350.,3000.],xtitle='Latitude',ytitle='Theta (K)'
;contour,tbar,latbin,thlev,color=0,/follow,levels=[270.],/overplot
;contour,tbar,latbin,thlev,color=mcolor,/follow,levels=[200.],/overplot
;set_viewport,.15,.85,.1,.45
;contour,tbarz,latbin,altitude,c_color=col1,color=0,/fill,levels=170.+5.*findgen(nlvls),/noerase,$
;        yrange=[10.,60.],xtitle='Latitude',ytitle='Altitude (km)'
;contour,tbarz,latbin,altitude,color=0,/follow,levels=[270.],/overplot
;contour,tbarz,latbin,altitude,color=mcolor,/follow,levels=[200.],/overplot
;stop
;
; declare time period arrays on first day
;
      if kcount eq 0L then begin
         tbar_zt=fltarr(nlat,nlv,kday)
         sdate_all=strarr(kday)
         kcount=1
      endif
      sdate_all(icount)=sdate
;
; save zonal means on each day
;
      tbar_zt(*,*,icount)=tbarz
skip:
      icount=icount+1L
goto,jump

saveit:
sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)
save,file='zt_geos5_at_mls_tbar.sav',tbar_zt,sdate_all,latbin,altitude

plotit:
restore,'zt_geos5_at_mls_tbar.sav
sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)
syear=strmid(sdate_all,0,4)
smon=strmid(sdate_all,4,2)
sday=strmid(sdate_all,6,2)
xindex=where(sday eq '15',nxticks)
xlabs=smon(xindex)
;
; loop over latitude
;
for j=nlat-6L,nlat-1L do begin
    slat=strcompress(string(format='(f5.1)',latbin(j)),/remove_all)
    if setplot eq 'ps' then begin
       lc=0
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       device,font_size=9
       device,/landscape,bits=8,filename='zt_geos5_at_mls_tbar_'+slat+'.ps'
       device,/color
       device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,$
              xsize=xsize,ysize=ysize
    endif
;
; plot zonal mean temperature and z'
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7
nlvls=19
col1=1+indgen(nlvls)*icolmax/nlvls
level=180.+5.*findgen(nlvls)
tzlat=transpose(reform(tbar_zt(j,*,*)))
print,slat,' ',min(tbar_zt),max(tbar_zt)
contour,tzlat,1.+findgen(kday),altitude,/noeras,xrange=[1.,kday],yrange=[20.,60.],$
      charsize=1.5,color=0,ytitle='Altitude (km)',title='GEOS-5 at MLS Tbar at '+slat+' Latitude',/fill,c_color=col1,$
      levels=level,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex
index=where(level gt 0)
contour,tzlat,1.+findgen(kday),altitude,levels=level(index),color=0,/follow,/overplot,c_labels=1+0*index,c_linestyle=5
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(K)'
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
       device, /close
       spawn,'convert -trim zt_geos5_at_mls_tbar_'+slat+'.ps -rotate -90 zt_geos5_at_mls_tbar_'+slat+'.jpg'
    endif
endfor	; loop over latitude
end
