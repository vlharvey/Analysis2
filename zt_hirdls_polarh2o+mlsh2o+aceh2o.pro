;
; water
; Plot 7-day average time-altitude Arctic HIRDLS
; Superimpose MLS contours
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
mdir='/aura6/data/MLS_data/Datfiles_SOSST/'
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
restore,dira+'h2o_ace_v2.2.2007'
dateace_all=date
yace_all=latitude
xace_all=longitude
modea_all=sctype
h2oace_all=mix
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
      restore,sdir+'h2o_hirdls_v2.04.19_'+sdate+'.sav'	; h2o
      print,sdate
      nlv=n_elements(altitude)
      index=where(mask eq -99.)
      if index(0) ne -1L then mix(index)=-99.
      hirh2omix=mix 
;
; declare time period arrays on first day
;
      if kcount eq 0L then begin
         hirpolarh2o_zt=fltarr(kday,nlv)
         mlspolarh2o_zt=fltarr(kday,nlv)
         acepolarh2o_zt=fltarr(kday,nlv)
         sdate_all=strarr(kday)
         kcount=1
      endif
      sdate_all(icount)=sdate
;
; compute daily polar averages
;
      polarh2o=fltarr(nlv)
      nh2oprof=lonarr(nlv)
      for ii=0L,n_elements(id)-1L do begin
          if latitude(ii) ge 60. then begin
             h2o_prof=reform(hirh2omix(ii,*))
             good=where(h2o_prof ne -99. and h2o_prof lt 1.e-5,ngood)
             if good(0) ne -1L then begin
                polarh2o(good)=polarh2o(good)+reform(h2o_prof(good))
                nh2oprof(good)=nh2oprof(good)+1L
             endif
          endif
      endfor
      good=where(nh2oprof gt 0L)
      if good(0) ne -1L then polarh2o(good)=polarh2o(good)/float(nh2oprof(good))
      hirpolarh2o_zt(icount,*)=polarh2o
;
; check 
;
erase
!type=2^2+2^3
set_viewport,.15,.85,.2,.85
plot,polarh2o,altitude,color=0,/noerase,title=sdate+' H!l2!nO',xrange=[2.e-6,1.e-5],/xlog,thick=3
skip:
;
; restore MLS H2O on this day
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
dum=findfile(mdir+'cat_mls_v2.2_'+sdate+'.sav')
if dum(0) eq '' then goto,skipmls
restore,mdir+'cat_mls_v2.2_'+sdate+'.sav'
restore,mdir+'h2o_mls_v2.2_'+sdate+'.sav'
index=where(mask eq -99.)
if index(0) ne -1L then mix(index)=-99.
mlsh2omix=mix
;
; MLS polar H2O
;
      polarh2o=fltarr(nlv)
      nh2oprof=lonarr(nlv)
      for ii=0L,n_elements(id)-1L do begin
          if latitude(ii) ge 60. then begin
             h2o_prof=reform(mlsh2omix(ii,*))
             good=where(h2o_prof ne -99.,ngood)
             if good(0) ne -1L then begin
                polarh2o(good)=polarh2o(good)+reform(h2o_prof(good))
                nh2oprof(good)=nh2oprof(good)+1L
             endif
          endif
      endfor
      good=where(nh2oprof gt 0L)
      if good(0) ne -1L then polarh2o(good)=polarh2o(good)/float(nh2oprof(good))
      if good(0) ne -1L then oplot,polarh2o(good),altitude(good),psym=2,color=0
      mlspolarh2o_zt(icount,*)=polarh2o
skipmls:
;
; extract ACE data today
;
      index=where(dateace_all eq date,nace)
      if index(0) eq -1L then goto,skipace
      yace_day=yace_all(index)
      aceh2omix=h2oace_all(index,*)
;
; ACE polar H2O
;
      polarh2o=fltarr(nlv)
      nh2oprof=lonarr(nlv)
      for ii=0L,nace-1L do begin
          if yace_day(ii) ge 50. then begin
             h2o_prof=reform(aceh2omix(ii,*))
             good=where(h2o_prof ne -99.,ngood)
             if good(0) ne -1L then begin
                polarh2o(good)=polarh2o(good)+reform(h2o_prof(good))
                nh2oprof(good)=nh2oprof(good)+1L
             endif
          endif
      endfor
      good=where(nh2oprof gt 0L)
      if good(0) ne -1L then polarh2o(good)=polarh2o(good)/float(nh2oprof(good))
      if good(0) ne -1L then oplot,polarh2o(good),altitude(good),psym=1,color=0
      acepolarh2o_zt(icount,*)=polarh2o
skipace:

      icount=icount+1L
goto,jump

plotit:
index=where(hirpolarh2o_zt ge 1.e-5 or hirpolarh2o_zt eq 0.)
if index(0) ne -1L then hirpolarh2o_zt(index)=0./0.
hirpolarh2o_zt=hirpolarh2o_zt*1.e6
hirpolarh2o_zt=smooth(hirpolarh2o_zt,7,/NaN,/edge_truncate)
index=where(finite(hirpolarh2o_zt) ne 1)
hirpolarh2o_zt(index)=0.
index=where(hirpolarh2o_zt gt 10.)
if index(0) ne -1L then hirpolarh2o_zt(index)=0.

mlspolarh2o_zt=mlspolarh2o_zt*1.e6

index=where(acepolarh2o_zt ge 1.e-5 or acepolarh2o_zt eq 0.)
if index(0) ne -1L then acepolarh2o_zt(index)=0./0.
acepolarh2o_zt=acepolarh2o_zt*1.e6
acepolarh2o_zt=smooth(acepolarh2o_zt,7,/NaN,/edge_truncate)
index=where(finite(acepolarh2o_zt) ne 1)
acepolarh2o_zt(index)=0.
index=where(acepolarh2o_zt gt 10.)
if index(0) ne -1L then acepolarh2o_zt(index)=0.
save,file='zt_hirdls_polarh2o+mlsh2o+aceh2o_2007.sav',hirpolarh2o_zt,mlspolarh2o_zt,$
     acepolarh2o_zt,kday,altitude,sdate_all
quick:
restore,'zt_hirdls_polarh2o+mlsh2o+aceh2o_2007.sav'
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
   device,/landscape,bits=8,filename='zt_hirdls_polarh2o+mlsh2o+aceh2o_'+syear(0)+'.ps'
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
level=1.+0.5*findgen(nlvls)
contour,hirpolarh2o_zt,1.+findgen(kday),altitude,/noeras,xrange=[1.,kday],yrange=[15.,57.],$
      charsize=1.5,color=0,ytitle='Altitude (km)',title='HIRDLS (color); MLS (black); ACE (white)',/fill,c_color=col1,$
      levels=level,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex,min_value=0.
contour,hirpolarh2o_zt,1.+findgen(kday),altitude,levels=level,color=mcolor,/follow,/overplot,c_labels=1+fltarr(nlvls)
index=where(mlspolarh2o_zt eq 0.)
if index(0) ne -1L then mlspolarh2o_zt(index)=0./0.
mlspolarh2o_zt=smooth(mlspolarh2o_zt,7,/NaN,/edge_truncate)
index=where(finite(mlspolarh2o_zt) ne 1)
if index(0) ne -1L then mlspolarh2o_zt(index)=0.
mlspolarh2o_zt(kday-1,*)=mlspolarh2o_zt(kday-2,*)
contour,mlspolarh2o_zt,1.+findgen(kday),altitude,levels=level,color=0,/follow,/overplot,$
        min_value=0.,c_labels=1+0*indgen(nlvls),thick=5
contour,acepolarh2o_zt,1.+findgen(kday),altitude,levels=level,color=mcolor,/follow,/overplot,$
        min_value=0.,c_labels=1+0*indgen(nlvls),thick=5
xyouts,xmn+0.02,ymn+0.02,syear(0),/normal,color=0,charsize=3,charthick=3
imin=min(level)
imax=max(level)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='7-day Avg H!l2!nO > 60 N (ppmv)',charsize=1.5
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
       spawn,'convert -trim zt_hirdls_polarh2o+mlsh2o+aceh2o_'+syear(0)+'.ps -rotate -90 zt_hirdls_polarh2o+mlsh2o+aceh2o_'+syear(0)+'.jpg'
    endif
end
