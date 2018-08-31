;
; Plot time-altitude Tbar for G5 and HIRDLS
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
sdir='/aura6/data/HIRDLS_data/Datfiles_SOSST/'
dir='/aura7/harvey/GEOS5_data/Datfiles/'
stimes=[$
'_0000.V01.',$
'_0600.V01.',$
'_1200.V01.',$
'_1800.V01.']
print,dir
lstmn=1
lstdy=15
lstyr=2008
ledmn=3
leddy=15
ledyr=2008
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
      ifile='DAS.ops.asm.tavg3d_dyn_v.GEOS510.MetO.'+sdate+stimes(0)+'nc3'
;
; read GEOS-5 data
;
      rd_geos5_nc3_meto,dir+ifile,nc,nr,nth,alon,alat,th,$
                  pv2,p2,msf2,u2,v2,q2,qdf2,mark2,sf2,vp2,iflag
      if iflag eq 1 then goto,jump
      t2=0.*pv2
      for k=0,nth-1 do t2(*,*,k)=th(k)*((p2(*,*,k)/1000.)^(.286))
      z2=(msf2-1004.*t2)/(9.86*1000.)
;
; zonal mean temperature
;
      tbarg=fltarr(nr,nth)
      zbarg=fltarr(nr,nth)
      for k=0L,nth-1L do begin
      for j=0L,nr-1L do begin
          tpts=reform(t2(j,*,k))
          zpts=reform(z2(j,*,k))
          index=where(tpts ne 0.,nx)
          if index(0) ne -1L then begin
              tbarg(j,k)=total(tpts(index))/float(nx)
              zbarg(j,k)=total(zpts(index))/float(nx)
          endif
      endfor
      endfor
;
; restore tpd_hirdls_v2.04.19_20080614.sav file
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
    dum=findfile(sdir+'tpd_hirdls_v2.04.19_'+sdate+'.sav')
    if dum(0) eq '' then goto,skip
    restore,sdir+'tpd_hirdls_v2.04.19_'+sdate+'.sav'
    restore,sdir+'cat_hirdls_v2.04.19_'+sdate+'.sav'	; latitude
    nlv=n_elements(altitude)
;
; compute zonal mean T, Z
;
    tbar=fltarr(nlat,nlv)
    nbar=lonarr(nlat,nlv)
    for ii=0L,n_elements(id)-1L do begin
        tmask_prof=reform(TEMPERATURE_MASK(ii,*))
        good=where(tmask_prof ne -99.,ngood)
        if good(0) ne -1L then begin
           ymean=latitude(ii)
           for j=0L,nlat-1L do begin
               if ymean ge latbin(j)-dy/2. and ymean lt latbin(j)+dy/2. then begin
                  tbar(j,good)=tbar(j,good)+temperature(ii,good)
                  nbar(j,good)=nbar(j,good)+1L
               endif
           endfor
        endif
    endfor
    good=where(nbar gt 0L)
    if good(0) ne -1L then tbar(good)=tbar(good)/float(nbar(good))
;
; interpolate GEOS temperature to height surfaces
;
tbarz=fltarr(nr,nlv)
for kk=0L,nlv-1L do begin
    zz=altitude(kk)
    for j=0L,nr-1L do begin
zprof=reform(zbarg(j,*))
if zz gt max(zprof) then goto,jumplev
        for k=1L,nth-1L do begin
            zup=zbarg(j,k-1) & zlw=zbarg(j,k)
            if zup ge zz and zlw le zz then begin
               zscale=(zup-zz)/(zup-zlw)
               tbarz(j,kk)=tbarg(j,k-1)+zscale*(tbarg(j,k)-tbarg(j,k-1))
;print,zlw,zz,zup,zscale
;print,tbarg(j,k),tbarz(j,kk),tbarg(j,k-1)
;stop
            endif
         endfor
      endfor
jumplev:
endfor

;nlvls=26
;col1=1+indgen(nlvls)*icolmax/nlvls
;erase
;!type=2^2+2^3
;set_viewport,.15,.85,.55,.9
;contour,tbar,latbin,altitude,c_color=col1,color=0,/fill,levels=170.+5.*findgen(nlvls),/noerase,title='HIRDLS',yrange=[10.,68.]
;contour,tbar,latbin,altitude,color=0,/follow,levels=[270.],/overplot
;contour,tbar,latbin,altitude,color=mcolor,/follow,levels=[200.],/overplot
;set_viewport,.15,.85,.1,.45
;contour,tbarz,alat,altitude,c_color=col1,color=0,/fill,levels=170.+5.*findgen(nlvls),/noerase,title='GEOS-5',yrange=[10.,68.]
;contour,tbarz,alat,altitude,color=0,/follow,levels=[270.],/overplot
;contour,tbarz,alat,altitude,color=mcolor,/follow,levels=[200.],/overplot
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
      tbar_zt(*,*,icount)=tbar
skip:
      icount=icount+1L
goto,jump

saveit:
sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)
save,file='zt_hirdls_tbar.sav',tbar_zt,sdate_all,latbin,altitude

plotit:
;restore,'zt_hirdls_tbar.sav
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
       device,/landscape,bits=8,filename='zt_hirdls+geos5_tbar_'+slat+'.ps'
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
contour,tzlat,1.+findgen(kday),altitude,/noeras,xrange=[1.,kday],yrange=[20.,68.],$
      charsize=1.5,color=0,ytitle='Altitude (km)',title='HIRDLS Tbar at '+slat+' Latitude',/fill,c_color=col1,$
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
       spawn,'convert -trim zt_hirdls+geos5_tbar_'+slat+'.ps -rotate -90 zt_hirdls+geos5_tbar_'+slat+'.jpg'
    endif
endfor	; loop over latitude
end
