;
; Plot time-altitude Arctic HIRDLS NO2
; Superimpose MLS CO 0.5 and 5 ppmv contours
; Superimpoase HIRDLS CH4 contours
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
mdir='/aura6/data/MLS_data/Datfiles_SOSST/'
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
      restore,sdir+'ch4_hirdls_v2.04.19_'+sdate+'.sav'  ; ch4
      index=where(mask eq -99.)
      if index(0) ne -1L then mix(index)=-99.
      n2omix=mix
;
; declare time period arrays on first day
;
      if kcount eq 0L then begin
         polart_zt=fltarr(kday,nlv)
         polarno2_zt=fltarr(kday,nlv)
         polarn2o_zt=fltarr(kday,nlv)
         mlspolarco_zt=fltarr(kday,nlv)
         sdate_all=strarr(kday)
         kcount=1
      endif
      sdate_all(icount)=sdate
;
; compute polar T, NO2
;
      polart=fltarr(nlv)
      polarno2=fltarr(nlv)
      polarn2o=fltarr(nlv)
      ntprof=lonarr(nlv)
      nno2prof=lonarr(nlv)
      nn2oprof=lonarr(nlv)
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
             n2o_prof=reform(n2omix(ii,*))
             good=where(n2o_prof ne -99. and n2o_prof lt 1.e-4,ngood)
             if good(0) ne -1L then begin
                polarn2o(good)=polarn2o(good)+reform(n2o_prof(good))
                nn2oprof(good)=nn2oprof(good)+1L
             endif
          endif
      endfor
      good=where(ntprof gt 0L)
      if good(0) ne -1L then polart(good)=polart(good)/float(ntprof(good))
      good=where(nno2prof gt 0L)
      if good(0) ne -1L then polarno2(good)=polarno2(good)/float(nno2prof(good))
      good=where(nn2oprof gt 0L)
      if good(0) ne -1L then polarn2o(good)=polarn2o(good)/float(nn2oprof(good))
;
; check HIRDLS
;
erase
!type=2^2+2^3
set_viewport,.15,.85,.55,.9
plot,polart,altitude,color=0,/noerase,title='Temperature',xrange=[170.,330.]
set_viewport,.15,.85,.1,.45
plot,polarn2o,altitude,color=0,/noerase,title='CH!l4!n',xrange=[1.e-7,1.e-5],/xlog
;
; save Arctic means on each day
;
      polart_zt(icount,*)=polart
      polarno2_zt(icount,*)=polarno2
      polarn2o_zt(icount,*)=polarn2o
skip:
;
; restore MLS CO on this day
; ALTITUDE        FLOAT     = Array[121]
; COMMENT         STRING    = Array[4]
; DATE            LONG      =     20060101
; ERR             FLOAT     = Array[3491, 121]
; FDOY            FLOAT     = Array[3491]
; ID              STRING    = Array[3491]
; LATITUDE        FLOAT     = Array[3491]
; LONGITUDE       FLOAT     = Array[3491]
; MASK            FLOAT     = Array[3491, 121]
; MIX             FLOAT     = Array[3491, 121]
; TIME            FLOAT     = Array[3491]
;
dum=findfile(mdir+'cat_mls_v2.2_'+sdate+'.sav')
if dum(0) eq '' then goto,skipmls
restore,mdir+'cat_mls_v2.2_'+sdate+'.sav'
restore,mdir+'co_mls_v2.2_'+sdate+'.sav'
index=where(mask eq -99.)
if index(0) ne -1L then mix(index)=-99.
mlscomix=mix
;
; compute polar CO
;
      mlspolarco=fltarr(nlv)
      mlsncoprof=lonarr(nlv)
      for ii=0L,n_elements(id)-1L do begin
          if latitude(ii) ge 60. then begin
             co_prof=reform(mlscomix(ii,*))
             good=where(co_prof ne -99.,ngood)
             if good(0) ne -1L then begin
                mlspolarco(good)=mlspolarco(good)+reform(co_prof(good))
                mlsncoprof(good)=mlsncoprof(good)+1L
             endif
          endif
      endfor
      good=where(mlsncoprof gt 0L)
      if good(0) ne -1L then mlspolarco(good)=mlspolarco(good)/float(mlsncoprof(good))
      mlspolarco_zt(icount,*)=mlspolarco
skipmls:
      icount=icount+1L
goto,jump

plotit:
index=where(polarno2_zt ge 1.e-7 or polarno2_zt eq 0.)
if index(0) ne -1L then polarno2_zt(index)=0./0.
polarno2_zt=polarno2_zt*1.e9
polarno2_zt=smooth(polarno2_zt,7,/NaN,/edge_truncate)
index=where(finite(polarno2_zt) ne 1)
polarno2_zt(index)=0.
index=where(polarno2_zt gt 50.)
if index(0) ne -1L then polarno2_zt(index)=0.

index=where(polarn2o_zt eq 0.)
if index(0) ne -1L then polarn2o_zt(index)=0./0.
polarn2o_zt=polarn2o_zt*1.e6
polarn2o_zt=smooth(polarn2o_zt,7,/NaN,/edge_truncate)
index=where(finite(polarn2o_zt) ne 1)
polarn2o_zt(index)=0.

mlspolarco_zt=mlspolarco_zt*1.e6
save,file='zt_hirdls_polarno2_ch4+mlsco_2006.sav',polart_zt,polarno2_zt,polarn2o_zt,$
     mlspolarco_zt,kday,altitude,sdate_all
quick:
restore,'zt_hirdls_polarno2_ch4+mlsco_2006.sav'
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
   device,/landscape,bits=8,filename='zt_hirdls_polarno2_ch4+mlsco_'+syear(0)+'.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
endif
;
; plot Arctic mean NO2, N2O, CH4, MLS CO
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
level=5.+2.5*findgen(nlvls)
contour,polarno2_zt,1.+findgen(kday),altitude,/noeras,xrange=[1.,kday],yrange=[25.,57.],$
      charsize=1.5,color=0,ytitle='Altitude (km)',title='HIRDLS Avg NO!l2!n and CH!l4!n and MLS CO > 60 N',/fill,c_color=col1,$
      nlevels=nlvls,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex,min_value=0.
contour,polarno2_zt,1.+findgen(kday),altitude,nlevels=nlvls,color=0,/follow,/overplot,c_labels=1+fltarr(nlvls)
index=where(mlspolarco_zt eq 0.)
if index(0) ne -1L then mlspolarco_zt(index)=0./0.
mlspolarco_zt=smooth(mlspolarco_zt,7,/NaN,/edge_truncate)
index=where(finite(mlspolarco_zt) ne 1)
if index(0) ne -1L then mlspolarco_zt(index)=0.
mlspolarco_zt(kday-1,*)=mlspolarco_zt(kday-2,*)
contour,mlspolarco_zt,1.+findgen(kday),altitude,levels=[0.05,0.08,0.1,0.2,0.5],color=0,/follow,/overplot,$
        min_value=0.,c_labels=[1,1,1],thick=5

contour,polarn2o_zt,1.+findgen(kday),altitude,levels=[0.1,1.,2.,3.],color=mcolor,/follow,/overplot,$
        min_value=0.,c_labels=[1,1,1,1,1],thick=5
xyouts,xmn+0.02,ymn+0.02,syear(0),/normal,color=0,charsize=3,charthick=3
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
       spawn,'convert -trim zt_hirdls_polarno2_ch4+mlsco_'+syear(0)+'.ps -rotate -90 zt_hirdls_polarno2_ch4+mlsco_'+syear(0)+'.jpg'
    endif
end
