;
; altitude-time of MLS PW1 in temperature
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
xorig=[0.20,0.20]
yorig=[0.15,0.55]
xlen=0.7
ylen=0.3
cbaryoff=0.07
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smonth=['J','F','M','A','M','J','J','A','S','O','N','D']
mdir='/atmos/aura6/data/MLS_data/Datfiles_SOSST/'
mdir2='/atmos/france/MLS_data/Datfiles_PW/'
lstmn=9
lstdy=1
lstyr=2006
ledmn=3
leddy=1
ledyr=lstyr+1L
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
;read,' Enter starting year ',lstyr
;read,' Enter ending year ',ledyr
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1991 then stop,'Year out of range '
if ledyr lt 1991 then stop,'Year out of range '
minyear=lstyr
maxyear=ledyr
yearlab=strcompress(maxyear,/remove_all)
if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)
;goto,quick

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
; restore MLS PW amplitude data
; DATE            STRING    = '20070614'
; LAT             FLOAT     = Array[36]
; MLSWAVE1HEIGHTASC FLOAT     = Array[36, 55]
; MLSWAVE1HEIGHTDES FLOAT     = Array[36, 55]
; MLSWAVE1PHASEASC FLOAT     = Array[36, 55]
; MLSWAVE1PHASEDES FLOAT     = Array[36, 55]
; MLSWAVE2HEIGHTASC FLOAT     = Array[36, 55]
; MLSWAVE2HEIGHTDES FLOAT     = Array[36, 55]
; MLSWAVE2PHASEASC FLOAT     = Array[36, 55]
; MLSWAVE2PHASEDES FLOAT     = Array[36, 55]
; PRESSURE        FLOAT     = Array[55]
;
      dum=findfile(mdir2+'Temperature_Wave_amplitudes_MLS_ascending_descending_'+sdate+'.sav')
      if dum(0) eq '' then goto,skipmls
      restore,dum
      if kcount eq 0L then begin
         print,lat
         rlat=-77.5
;        read,'Enter desired latitude for PW amplitudes ',rlat
         index=where(abs(lat-rlat) eq min(abs(lat-rlat)))
         ilat=index(0)
         slat1=strcompress(lat(ilat),/r)
         pressure_altitude=alog(1000./pressure)*7.
         nlv2=n_elements(pressure)
      endif
;
; restore MLS CO on this day
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
      if dum(0) eq '' then goto,skipmls
      restore,mdir+'cat_mls_v3.3_'+sdate+'.sav'
      restore,mdir+'co_mls_v3.3_'+sdate+'.sav'
      restore,mdir+'tpd_mls_v3.3_'+sdate+'.sav'
      print,sdate
;
; apply mask
;
      index=where(mask eq -99.)
      if index(0) ne -1L then mix(index)=-99.
      mlscomix=mix
      index=where(temperature_mask eq -99.)
      if index(0) ne -1L then temperature(index)=-99.
      mlscomix=mix
      nlv=n_elements(altitude)
;
; declare time period arrays on first day
;
      if kcount eq 0L then begin
         mlspw1a_zt=fltarr(kday,nlv2)
         mlspw1d_zt=fltarr(kday,nlv2)
         mlspw1a_phase_zt=fltarr(kday,nlv2)
         mlspw1d_phase_zt=fltarr(kday,nlv2)
         mlspolarco_zt=fltarr(kday,nlv)
         mlspolartp_zt=fltarr(kday,nlv)
         sdate_all=strarr(kday)
         kcount=1
      endif
      sdate_all(icount)=sdate
;
; wave PW amplitudes at ilat
;
      mlspw1a_zt(icount,*)=reform(MLSWAVE1HEIGHTASC(ilat,*))
      mlspw1d_zt(icount,*)=reform(MLSWAVE1HEIGHTDES(ilat,*))
      mlspw1a_phase_zt(icount,*)=reform(MLSWAVE1PHASEASC(ilat,*))
      mlspw1d_phase_zt(icount,*)=reform(MLSWAVE1PHASEDES(ilat,*))
;
; compute polar CO and temp
;
      mlspolarco=fltarr(nlv)
      mlsncoprof=lonarr(nlv)
      mlspolartp=fltarr(nlv)
      mlsntpprof=lonarr(nlv)
      for ii=0L,n_elements(id)-1L do begin
          if latitude(ii) le -60. then begin
             co_prof=reform(mlscomix(ii,*))
             good=where(co_prof ne -99.,ngood)
             if good(0) ne -1L then begin
                mlspolarco(good)=mlspolarco(good)+reform(co_prof(good))
                mlsncoprof(good)=mlsncoprof(good)+1L
             endif
             tp_prof=reform(temperature(ii,*))
             good=where(tp_prof ne -99.,ngood)
             if good(0) ne -1L then begin
                mlspolartp(good)=mlspolartp(good)+reform(tp_prof(good))
                mlsntpprof(good)=mlsntpprof(good)+1L
             endif
          endif
      endfor
      good=where(mlsncoprof gt 0L)
      if good(0) ne -1L then mlspolarco(good)=mlspolarco(good)/float(mlsncoprof(good))
      good=where(mlsntpprof gt 0L)
      if good(0) ne -1L then mlspolartp(good)=mlspolartp(good)/float(mlsntpprof(good))
      mlspolarco_zt(icount,*)=mlspolarco
      mlspolartp_zt(icount,*)=mlspolartp
skipmls:
      icount=icount+1L
goto,jump

plotit:
;
; convert phase in radians to degrees longitude. phase is where sin=0 so ridge is 90 degrees to the east
;
      index=where(mlspw1a_phase_zt ne -999.)
      if index(0) ne -1 then mlspw1a_phase_zt(index)=mlspw1a_phase_zt(index)/!dtor mod 360. +90.
      index=where(mlspw1d_phase_zt ne -999.)
      if index(0) ne -1 then mlspw1d_phase_zt(index)=mlspw1d_phase_zt(index)/!dtor mod 360. +90.
;
; interpolate small gaps in time
;
for k=0,nlv-1 do begin
    dlev=reform(mlspolarco_zt(*,k))
    for i=1,kday-1 do begin
        if dlev(i) eq 0. and dlev(i-1) ne 0. then begin
           for ii=i+1,kday-1 do begin
               naway=float(ii-i)
               if naway le 5.0 and dlev(ii) ne 0. then begin
                  dlev(i)=(naway*dlev(i-1)+dlev(ii))/(naway+1.0)
                  goto,jump1
               endif
           endfor
jump1:
        endif
    endfor
    mlspolarco_zt(*,k)=dlev
    dlev=reform(mlspolartp_zt(*,k))
    for i=1,kday-1 do begin
        if dlev(i) eq 0. and dlev(i-1) ne 0. then begin
           for ii=i+1,kday-1 do begin
               naway=float(ii-i)
               if naway le 5.0 and dlev(ii) ne 0. then begin
                  dlev(i)=(naway*dlev(i-1)+dlev(ii))/(naway+1.0)
                  goto,jump2
               endif
           endfor
jump2:
        endif
    endfor
    mlspolartp_zt(*,k)=dlev
endfor
;
; year date label
;
syear=strmid(sdate_all,0,4)
good=where(long(syear) ne 0L)
minyear=long(min(long(syear(good))))
maxyear=long(max(long(syear)))
yearlab=strcompress(maxyear,/remove_all)
if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)
mlspolarco_zt=mlspolarco_zt*1.e6
;
; save temp, co, etc
;
;save,file='zt_mls_temp+co_'+yearlab+'_sh.sav',mlspolarco_zt,mlspolartp_zt,kday,altitude,sdate_all
quick:
;restore,'zt_mls_temp+co_'+yearlab+'_sh.sav'
sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)
syear=strmid(sdate_all,0,4)
smon=strmid(sdate_all,4,2)
sday=strmid(sdate_all,6,2)
xindex=where(sday eq '15',nxticks)
xlabs=smon(xindex)+'/'+sday(xindex)
good=where(long(syear) ne 0L)
minyear=long(min(long(syear(good))))
maxyear=long(max(long(syear)))
yearlab=strcompress(maxyear,/remove_all)
if minyear ne maxyear then yearlab=strcompress(minyear,/remove_all)+'-'+strcompress(maxyear,/remove_all)
;
if setplot eq 'ps' then begin
   lc=0
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !p.font=0
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='zt_mls_tppw_'+yearlab+'_sh.ps'
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=5
   !p.charthick=5
   !y.thick=2
   !x.thick=2
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
!type=2^2+2^3+2^7	; ticks outward
level=[0.001,0.01,0.025,0.05,0.1,0.25,0.5,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,15.]
index=where(mlspolarco_zt eq 0.)
if index(0) ne -1L then mlspolarco_zt(index)=0./0.
mlspolarco_zt=smooth(mlspolarco_zt,7,/NaN,/edge_truncate)
if index(0) ne -1L then mlspolarco_zt(index)=0./0.
index=where(mlspolartp_zt eq 0.)
if index(0) ne -1L then mlspolartp_zt(index)=0./0.
tlevel=120.+5.*findgen(35)
nlvls=n_elements(tlevel)
col1=1+indgen(nlvls)*icolmax/nlvls
;
; t-tendency
;
n0=findgen(kday)
n1=1+findgen(kday)
ttend=0.*mlspolartp_zt
for k=0,nlv-1L do ttend(*,k)=mlspolartp_zt(n1,k)-mlspolartp_zt(n0,k)
;mlspolartp_zt=ttend
;tlevel=-2.+0.2*findgen(21)
;nlvls=n_elements(tlevel)
;col1=1+indgen(nlvls)*icolmax/nlvls
loadct,0
contour,smooth(mlspolartp_zt,5,/edge_truncate,/Nan),1.+findgen(kday),altitude,/noeras,xrange=[1.,kday],yrange=[15.,90.],$
      charsize=1.5,color=0,ytitle='Altitude (km)',title=yearlab+' PW1 Amplitudes at '+slat1,/cell_fill,c_color=col1,$
      levels=tlevel,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex,min_value=-99.
index=where(tlevel gt 0.)
contour,smooth(mlspolartp_zt,5,/edge_truncate,/Nan),1.+findgen(kday),altitude,levels=tlevel(index),color=0,/follow,/overplot,c_labels=fltarr(nlvls)
index=where(tlevel lt 0.)
contour,smooth(mlspolartp_zt,5,/edge_truncate,/Nan),1.+findgen(kday),altitude,levels=tlevel(index),color=mcolor,c_linestyle=5,/follow,/overplot,c_labels=fltarr(nlvls)
zlevel=findgen(30)
nlvls2=n_elements(zlevel)
col2=70+indgen(nlvls2)*icolmax/(nlvls2+1)
loadct,39
contour,smooth(mlspw1a_zt,3,/edge_truncate,/Nan),1.+findgen(kday),pressure_altitude,/overplot,levels=zlevel,c_color=col2,thick=4
xyouts,xmn,ymn+0.02,'ASC',/normal,color=0,charsize=3,charthick=3
;
; print end date
;
maxdate=max(long(sdate_all))
smaxdate=strcompress(maxdate,/remove_all)
datelab=strmid(smaxdate,4,2)+'/'+strmid(smaxdate,6,2)
;xyouts,(xmx+xmn)/2.,ymn+.02,'Last Day '+datelab,color=0,/normal,charsize=2

xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7       ; ticks outward
level=[0.001,0.01,0.025,0.05,0.1,0.25,0.5,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,15.]
index=where(mlspolarco_zt eq 0.)
if index(0) ne -1L then mlspolarco_zt(index)=0./0.
mlspolarco_zt=smooth(mlspolarco_zt,7,/NaN,/edge_truncate)
if index(0) ne -1L then mlspolarco_zt(index)=0./0.
index=where(mlspolartp_zt eq 0.)
if index(0) ne -1L then mlspolartp_zt(index)=0./0.
tlevel=120.+5.*findgen(35)
nlvls=n_elements(tlevel)
col1=1+indgen(nlvls)*icolmax/nlvls
;
; t-tendency
;
n0=findgen(kday)
n1=1+findgen(kday)
ttend=0.*mlspolartp_zt
for k=0,nlv-1L do ttend(*,k)=mlspolartp_zt(n1,k)-mlspolartp_zt(n0,k)
;mlspolartp_zt=ttend
;tlevel=-2.+0.2*findgen(21)
;nlvls=n_elements(tlevel)
;col1=1+indgen(nlvls)*icolmax/nlvls
loadct,0
contour,smooth(mlspolartp_zt,5,/edge_truncate,/Nan),1.+findgen(kday),altitude,/noeras,xrange=[1.,kday],yrange=[15.,90.],$
      charsize=1.5,color=0,ytitle='Altitude (km)',title=yearlab+' PW1 Amplitudes at '+slat1,/cell_fill,c_color=col1,$
      levels=tlevel,xticks=nxticks-1,xtickname=xlabs,xtickv=xindex,min_value=-99.
index=where(tlevel gt 0.)
contour,smooth(mlspolartp_zt,5,/edge_truncate,/Nan),1.+findgen(kday),altitude,levels=tlevel(index),color=0,/follow,/overplot,c_labels=fltarr(nlvls)
index=where(tlevel lt 0.)
contour,smooth(mlspolartp_zt,5,/edge_truncate,/Nan),1.+findgen(kday),altitude,levels=tlevel(index),color=mcolor,c_linestyle=5,/follow,/overplot,c_labels=fltarr(nlvls)
zlevel=findgen(30)
nlvls2=n_elements(zlevel)
col2=70+indgen(nlvls2)*icolmax/(nlvls2+1)
loadct,39
contour,smooth(mlspw1d_zt,3,/edge_truncate,/Nan),1.+findgen(kday),pressure_altitude,/overplot,levels=zlevel,c_color=col2,thick=4
xyouts,xmn,ymn+0.02,'DSC',/normal,color=0,charsize=3,charthick=3

loadct,0
imin=min(tlevel)
imax=max(tlevel)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='MLS dT/dt < 60 S (K)'
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
       spawn,'convert -trim zt_mls_tppw_'+yearlab+'_sh.ps -rotate -90 zt_mls_tppw_'+yearlab+'_sh.jpg'
    endif
end
