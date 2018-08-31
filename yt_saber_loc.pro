;
; Plot time-altitude Tbar 
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
sdir='/aura6/data/SABER_data/Datfiles/'
goto,plotit
lstmn=1
lstdy=1
lstyr=2007
ledmn=12
leddy=31
ledyr=2007
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
latmin=fltarr(kday)
latmax=fltarr(kday)
doy=fltarr(kday)
sdate_all=strarr(kday)
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
      if ndays gt ledday then stop,'Normal Termination Condition'
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
print,sdate
doy(icount)=iday
;
; restore SABER TPZ files
;
    dum=findfile(sdir+'SABER_TPZ_'+sdate+'.sav')
    if dum(0) eq '' then goto,skip
    restore,sdir+'SABER_TPZ_'+sdate+'.sav'
    nlv=n_elements(altitude)
;
; declare time period arrays on first day
;
      sdate_all(icount)=sdate
;
; save min/max lat each day
;
      index=where(latitude gt -90.)
      latmin(icount)=min(latitude(index))
      latmax(icount)=max(latitude(index))
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
       !p.thick=2.0                   ;Plotted lines twice as thick
       device,font_size=9
       device,/landscape,bits=8,filename='yt_saber_loc.ps'
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
!type=2^2+2^3
nlvls=19
col1=1+indgen(nlvls)*icolmax/nlvls
level=180.+5.*findgen(nlvls)

syr=strmid(sdate_all,0,4)
smn=strmid(sdate_all,4,2)
sdy=strmid(sdate_all,6,2)
index=where(sdy eq '15',nxticks)
xlab=smn(index)
plot,doy,-90.+findgen(181),/nodata,/noeras,xrange=[1.,kday],yrange=[-90.,90.],$
     charsize=1.5,color=0,ytitle='Latitude',title='SABER',xticks=nxticks-1,xtickv=index,xtickname=xlab
for i=0L,kday-1L do begin
    plots,i+1,latmin(i)
    plots,i+1,latmax(i),/continue,color=0,thick=3
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
       spawn,'convert -trim yt_saber_loc.ps -rotate -90 yt_saber_loc.jpg'
    endif
end
