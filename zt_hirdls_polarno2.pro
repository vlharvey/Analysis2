;
; Plot time-altitude Arctic HIRDLS NO2
;
@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3_meto

loadct,38
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
xorig=[0.20]
yorig=[0.25]
xlen=0.7
ylen=0.5
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
sdir='/aura6/data/HIRDLS_data/Datfiles_SOSST/'
lstmn=1
lstdy=1
lstyr=2006
ledmn=4
leddy=1
ledyr=2006
lstday=0
ledday=0
goto,quick
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
      restore,sdir+'no2_hirdls_v2.04.19_'+sdate+'.sav'	; no2
      print,sdate
      nlv=n_elements(altitude)
      index=where(mask eq -99.)
      if index(0) ne -1L then mix(index)=-99.
      no2mix=mix 
;
; declare time period arrays on first day
;
      if kcount eq 0L then begin
         polart_zt=fltarr(kday,nlv)
         polarno2_zt=fltarr(kday,nlv)
         sdate_all=strarr(kday)
         kcount=1
      endif
      sdate_all(icount)=sdate
;
; compute polar T, Z
;
      polart=fltarr(nlv)
      polarno2=fltarr(nlv)
      ntprof=lonarr(nlv)
      nno2prof=lonarr(nlv)
      for ii=0L,n_elements(id)-1L do begin
          if latitude(ii) ge 60. then begin
             tmask_prof=reform(TEMPERATURE_MASK(ii,*))
             good=where(tmask_prof ne -99.,ngood)
             if good(0) ne -1L then begin
                polart(good)=polart(good)+reform(temperature(ii,good))
                ntprof(good)=ntprof(good)+1L
             endif
             no2_prof=reform(no2mix(ii,*))
             good=where(no2_prof ne -99. and no2_prof lt 1.e-7,ngood)
             if good(0) ne -1L then begin
                polarno2(good)=polarno2(good)+reform(no2_prof(good))
                nno2prof(good)=nno2prof(good)+1L
             endif
          endif
      endfor
      good=where(ntprof gt 0L)
      if good(0) ne -1L then polart(good)=polart(good)/float(ntprof(good))
      good=where(nno2prof gt 0L)
      if good(0) ne -1L then polarno2(good)=polarno2(good)/float(nno2prof(good))
;erase
;!type=2^2+2^3
;set_viewport,.15,.85,.55,.9
;plot,polart,altitude,color=0,/noerase,title='Temperature',xrange=[170.,330.]
;set_viewport,.15,.85,.1,.45
;plot,polarno2,altitude,color=0,/noerase,title='NO!l2!n',xrange=[1.e-9,1.e-6],/xlog
;stop
;
; save Arctic means on each day
;
      polart_zt(icount,*)=polart
      polarno2_zt(icount,*)=polarno2
skip:
      icount=icount+1L
goto,jump

plotit:
save,file='zt_hirdls_polarno2_2006.sav',polart_zt,polarno2_zt,kday,altitude,sdate_all
quick:
restore,'zt_hirdls_polarno2_2006.sav'
sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)
syear=strmid(sdate_all,0,4)
smon=strmid(sdate_all,4,2)
sday=strmid(sdate_all,6,2)
xindex=where(sday eq '15',nxticks)
xlabs=smon(xindex)
if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !psym=0
   !p.font=0
   device,font_size=9
   device,/landscape,bits=8,filename='zt_hirdls_polarno2_'+sdate0+'.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
endif
;
; plot Arctic mean temperature and NO2
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
nlvls=19
col1=1+indgen(nlvls)*icolmax/nlvls
;index=where(polarno2_zt ge 1.e-7 or polarno2_zt eq 0.)
;if index(0) ne -1L then polarno2_zt(index)=0./0.
;polarno2_zt=polarno2_zt*1.e9
;polarno2_zt=smooth(polarno2_zt,7,/NaN)
;index=where(finite(polarno2_zt) ne 1)
;polarno2_zt(index)=0.
index=where(polarno2_zt gt 50.)
if index(0) ne -1L then polarno2_zt(index)=0.
level=5.+2.5*findgen(nlvls)
contour,polarno2_zt,1.+findgen(kday),altitude,/noeras,xrange=[1.,kday],yrange=[25.,55.],$
      charsize=1.5,color=0,ytitle='Altitude (km)',title='HIRDLS Avg NO!l2!n > 60 N',/fill,c_color=col1,$
      nlevels=nlvls,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex,min_value=0.
contour,polarno2_zt,1.+findgen(kday),altitude,nlevels=nlvls,color=0,/follow,/overplot,c_labels=1+0*index,c_linestyle=5
xyouts,xmn+0.02,ymn+0.02,syear(0),/normal,color=0,charsize=3,charthick=3
;
; why can't I set my own contour levels?
;
;contour,polarno2_zt,1.+findgen(kday),altitude,/noeras,xrange=[1.,kday],yrange=[25.,55.],$
;      charsize=1.5,color=0,ytitle='Altitude (km)',title='HIRDLS Arctic NO!l2!n',/fill,c_color=col1,$
;      levels=levels,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex,min_value=0.
;contour,polarno2_zt,1.+findgen(kday),altitude,levels=levels,color=0,/follow,/overplot,c_labels=1+0*index,c_linestyle=5

imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(ppbv)'
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
       spawn,'convert -trim zt_hirdls_polarno2_'+sdate0+'.ps -rotate -90 zt_hirdls_polarno2_'+sdate0+'.jpg'
    endif
end
