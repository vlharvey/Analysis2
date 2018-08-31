;
; timeseries of daily average local time
;
@stddat
@kgmt
@ckday
@kdate

loadct,39
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
device,decompose=0
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
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
sdir='/Volumes/earth/harvey/SABER_data/Datfiles/'
lstmn=1
lstdy=1
lstyr=2006
ledmn=12
leddy=31
ledyr=2006
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
;print, ' '
;print, '      SABER Version '
;print, ' '
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
;goto,plotit

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
; restore SABER TPZ files
;
; ALTITUDE        FLOAT     = Array[201]
; COMMENT         STRING    = Array[22]
; DATE            LONG      = Array[1317]
; GPALTITUDE      FLOAT     = Array[1317, 201]
; IAURORA         LONG      = Array[1317]
; ID              STRING    = Array[1317]
; JDATE           LONG      = Array[1317]
; LATITUDE        FLOAT     = Array[1317, 201]
; LONGITUDE       FLOAT     = Array[1317, 201]
; MODE            INT       = Array[1317]
; NEVENT          LONG      =         1317
; PRESSURE        FLOAT     = Array[1317, 201]
; SOLAP           FLOAT     =       14.0000
; SOLF10P7DAILY   FLOAT     =       69.4000
; SOLKP           FLOAT     =       2.00000
; SOLSPOTNO       LONG      = Array[1317]
; TEMPERATURE     FLOAT     = Array[1317, 201]
; TIME            FLOAT     = Array[1317, 201]
; TPAD            INT       = Array[1317]
; TPDN            INT       = Array[1317]
; TPSOLARLT       FLOAT     = Array[1317, 201]
; TPSOLARZEN      FLOAT     = Array[1317, 201]
; VER             STRING    = '01.07'
;
    dum=findfile(sdir+'SABER_TPZ_'+sdate+'.sav')
    if dum(0) eq '' then goto,skip
    restore,sdir+'SABER_TPZ_'+sdate+'.sav'
    nlv=n_elements(altitude)
;
;erase
index=where(TPSOLARLT ne -999.)       
;plot,TPSOLARLT(index),psym=3,title=sdate,color=0
;erase
;index=where(TPSOLARZEN ne -999.)
;plot,TPSOLARZEN(index),psym=3,title=sdate,color=0,yrange=[0.,180.]
;print,sdate,' ',min(TPSOLARZEN(index)),mean(TPSOLARZEN(index)),max(TPSOLARZEN(index))
;wait,.1
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
;
; save zonal means on each day
;
      ltime_mean(icount)=mean(TPSOLARLT(index))
      zen_mean(icount)=mean(TPSOLARZEN(index))
      zen_min(icount)=min(TPSOLARZEN(index))
      zen_max(icount)=max(TPSOLARZEN(index))
skip:
      icount=icount+1L
goto,jump

saveit:
sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)
;save,file='zt_saber_tbar.sav',tbar_zt,sdate_all,latbin,altitude

plotit:
;restore,'zt_saber_tbar.sav
sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)
syear=strmid(sdate_all,0,4)
smon=strmid(sdate_all,4,2)
sday=strmid(sdate_all,6,2)
xindex=where(sday eq '01',nxticks)
xlabs=smon(xindex)
;
    if setplot eq 'ps' then begin
       lc=0
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !psym=0
       !p.font=0
       !p.thick=2.0                   ;Plotted lines twice as thick
       device,font_size=9
       device,/landscape,bits=8,filename='timeseries_ltime.ps'
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
;
; hardwire xticks
;
;xindex=[1.,31.+1.,31.+29.+1.]	; 15th of each month as record goes Jan 15 to Mar 15
;xlabs=['J','F','M']
nxticks=n_elements(xlabs)
plot,1.+findgen(kday),zen_mean,/noeras,xrange=[1.,kday],yrange=[0.,180.],$
      charsize=1.5,color=0,ytitle='Solar Zenith Angle',title='SABER Daily Min/Mean/Max',psym=8,$
      xticks=nxticks-1,xtickname=xlabs,xtickv=xindex
for i=0L,kday-1L do begin
    plots,i,zen_min(i)
    plots,i,zen_max(i),/continue,color=0
endfor
;imin=min(level)
;imax=max(level)
;ymnb=yorig(0) -cbaryoff
;ymxb=ymnb  +cbarydel
;set_viewport,xmn,xmx,ymnb,ymxb
;!type=2^2+2^3+2^6
;plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='(K)'
;ybox=[0,10,10,0,0]
;x1=imin
;dx=(imax-imin)/float(nlvls)
;for jj=0,nlvls-1 do begin
;xbox=[x1,x1,x1+dx,x1+dx,x1]
;polyfill,xbox,ybox,color=col1(jj)
;x1=x1+dx
;endfor

    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert -trim timeseries_ltime.ps -rotate -90 timeseries_ltime.jpg'
    endif
end
