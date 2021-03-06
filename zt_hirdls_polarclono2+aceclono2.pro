;
; clono2
; Plot 7-day average time-altitude Arctic HIRDLS
; Superimpose ACE contours
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
dira='/aura3/data/ACE_data/Datfiles_SOSST/v2.2/'
lstmn=1
lstdy=1
lstyr=2007
ledmn=4
leddy=1
ledyr=2007
lstday=0
ledday=0
;goto,quick
;
; restore year of ACE data
;
restore,dira+'cat_ace_v2.2.2007'
restore,dira+'clono2_ace_v2.2.2007'
dateace_all=date
yace_all=latitude
xace_all=longitude
modea_all=sctype
clono2ace_all=mix
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
      restore,sdir+'clono2_hirdls_v2.04.19_'+sdate+'.sav'	; clono2
      print,sdate
      nlv=n_elements(altitude)
      index=where(mask eq -99.)
      if index(0) ne -1L then mix(index)=-99.
      hirclono2mix=mix 
;
; declare time period arrays on first day
;
      if kcount eq 0L then begin
         hirpolarclono2_zt=fltarr(kday,nlv)
         acepolarclono2_zt=fltarr(kday,nlv)
         sdate_all=strarr(kday)
         kcount=1
      endif
      sdate_all(icount)=sdate
;
; compute daily polar averages
;
      polarclono2=fltarr(nlv)
      nclono2prof=lonarr(nlv)
      for ii=0L,n_elements(id)-1L do begin
          if latitude(ii) ge 60. then begin
             clono2_prof=reform(hirclono2mix(ii,*))
             good=where(clono2_prof ne -99. and clono2_prof lt 1.e-7,ngood)
             if good(0) ne -1L then begin
                polarclono2(good)=polarclono2(good)+reform(clono2_prof(good))
                nclono2prof(good)=nclono2prof(good)+1L
             endif
          endif
      endfor
      good=where(nclono2prof gt 0L)
      if good(0) ne -1L then polarclono2(good)=polarclono2(good)/float(nclono2prof(good))
      hirpolarclono2_zt(icount,*)=polarclono2
;
; check 
;
erase
!type=2^2+2^3
set_viewport,.15,.85,.2,.85
plot,polarclono2,altitude,color=0,/noerase,title=sdate+' ClONO!l2!n',xrange=[1.e-10,1.e-7],/xlog,thick=3
skip:
;
; extract ACE data today
;
      index=where(dateace_all eq date,nace)
      if index(0) eq -1L then goto,skipace
      yace_day=yace_all(index)
      aceclono2mix=clono2ace_all(index,*)
;
; ACE polar
;
      polarclono2=fltarr(nlv)
      nclono2prof=lonarr(nlv)
      for ii=0L,nace-1L do begin
          if yace_day(ii) ge 50. then begin
             clono2_prof=reform(aceclono2mix(ii,*))
             good=where(clono2_prof ne -99.,ngood)
             if good(0) ne -1L then begin
                polarclono2(good)=polarclono2(good)+reform(clono2_prof(good))
                nclono2prof(good)=nclono2prof(good)+1L
             endif
          endif
      endfor
      good=where(nclono2prof gt 0L)
      if good(0) ne -1L then polarclono2(good)=polarclono2(good)/float(nclono2prof(good))
      if good(0) ne -1L then oplot,polarclono2(good),altitude(good),psym=1,color=0
      acepolarclono2_zt(icount,*)=polarclono2
skipace:

      icount=icount+1L
goto,jump

plotit:
index=where(hirpolarclono2_zt ge 1.e-7 or hirpolarclono2_zt eq 0.)
if index(0) ne -1L then hirpolarclono2_zt(index)=0./0.
hirpolarclono2_zt=hirpolarclono2_zt*1.e9
hirpolarclono2_zt=smooth(hirpolarclono2_zt,7,/NaN,/edge_truncate)
index=where(finite(hirpolarclono2_zt) ne 1)
hirpolarclono2_zt(index)=0.
index=where(hirpolarclono2_zt gt 100.)
if index(0) ne -1L then hirpolarclono2_zt(index)=0.

index=where(acepolarclono2_zt ge 1.e-7 or acepolarclono2_zt eq 0.)
if index(0) ne -1L then acepolarclono2_zt(index)=0./0.
acepolarclono2_zt=acepolarclono2_zt*1.e9
acepolarclono2_zt=smooth(acepolarclono2_zt,7,/NaN,/edge_truncate)
index=where(finite(acepolarclono2_zt) ne 1)
acepolarclono2_zt(index)=0.
index=where(acepolarclono2_zt gt 1000.)
if index(0) ne -1L then acepolarclono2_zt(index)=0.
save,file='zt_hirdls_polarclono2+aceclono2_2007.sav',hirpolarclono2_zt,$
     acepolarclono2_zt,kday,altitude,sdate_all
quick:
restore,'zt_hirdls_polarclono2+aceclono2_2007.sav'
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
   device,/landscape,bits=8,filename='zt_hirdls_polarclono2+aceclono2_'+syear(0)+'.ps'
   device,/color
   device,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize
endif
;
; plot Arctic means
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3
nlvls=17
col1=1+indgen(nlvls)*icolmax/nlvls
level=[0.01,0.1,0.25,0.5,0.75,1.+0.5*findgen(nlvls)]
nlvls=n_elements(level)
col1=1+indgen(nlvls)*icolmax/nlvls
contour,hirpolarclono2_zt,1.+findgen(kday),altitude,/noeras,xrange=[1.,kday],yrange=[15.,57.],$
      charsize=1.5,color=0,ytitle='Altitude (km)',title='HIRDLS (color); ACE (white)',/fill,c_color=col1,$
      levels=level,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex,min_value=0.
contour,hirpolarclono2_zt,1.+findgen(kday),altitude,levels=level,color=mcolor,/follow,/overplot,c_labels=1+fltarr(nlvls)
contour,acepolarclono2_zt,1.+findgen(kday),altitude,levels=level,color=mcolor,/follow,/overplot,$
        min_value=0.,c_labels=1+0*indgen(nlvls),thick=5
xyouts,xmn+0.02,ymn+0.02,syear(0),/normal,color=0,charsize=3,charthick=3
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='7-day Avg ClONO!l2!n > 60 N (ppbv)',charsize=1.5
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
       spawn,'convert -trim zt_hirdls_polarclono2+aceclono2_'+syear(0)+'.ps -rotate -90 zt_hirdls_polarclono2+aceclono2_'+syear(0)+'.jpg'
    endif
end
