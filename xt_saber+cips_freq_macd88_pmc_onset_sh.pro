;
; Hovmoller of SABER temperature 1 Nov to 1 Jan
; overplot CIPS frequency
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
xorig=[0.25]
yorig=[0.25]
xlen=0.5
ylen=0.7
cbaryoff=0.15
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smonth=['J','F','M','A','M','J','J','A','S','O','N','D']
dir='/Volumes/Data/earth/harvey/SABER_data/Datfiles/SABER_Temp_O3_'
start_year=[2007,2008,2009,2010,2011,2012,2013]
start_date=[-18,-18,-31,-3,-23,-27,-28]
end_date=[48,60,62,54,56,59,61]
nyear=n_elements(start_year)
for nn=start_year(0),start_year(nyear-1L) do begin

lstmn=11
lstdy=10
lstyr=nn
ledmn=12
leddy=20
ledyr=nn
lstday=0
ledday=0
;
; loop over years
;
for iyear=lstyr,ledyr do begin
for ilat=-80,-50,5 do begin

syr=string(format='(i2.2)',iyear-2000)
slt=strcompress(abs(ilat),/remove_all)
;restore,'/Volumes/Data/CIPS_data/Pre_process/xt_cips_3c_SH'+syr+'_DSC_'+slt+'Lat_2G_freq.sav
;decfreq=XTFREQ
;restore,'/Volumes/Data/CIPS_data/Pre_process/xt_cips_3c_SH'+syr+'_ASC_'+slt+'Lat_2G_freq.sav
;ascfreq=XTFREQ
;cipsfreq=0.*decfreq
;index=where(decfreq ne 0. and ascfreq eq 0.)
;if index(0) ne -1L then cipsfreq(index)=decfreq(index)
;index=where(decfreq eq 0. and ascfreq ne 0.)
;if index(0) ne -1L then cipsfreq(index)=ascfreq(index)
;index=where(decfreq ne 0. and ascfreq ne 0.)
;if index(0) ne -1L then cipsfreq(index)=(decfreq(index)+ascfreq(index))/2.

kcount=0
rlat=float(ilat)
slat=strcompress(long(rlat),/remove_all)

;goto,quick

z = stddat(lstmn,lstdy,iyear,lstday)
z = stddat(ledmn,leddy,iyear,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
sdate_all=strarr(kday)
;
; Compute initial Julian date
;
iyr = iyear
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L
;
; longitude grid
;
dx=15.
nc=long(360./dx)+1
longrid=dx*findgen(nc)

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
      sdate_all(icount)=sdate
      print,sdate
;
; restore SABER on this day
; ALTITUDE        FLOAT     = Array[121]
; COMMENT         STRING    = Array[15]
; GPALTITUDE      FLOAT     = Array[1434, 121]
; LATITUDE        FLOAT     = Array[1434, 121]
; LONGITUDE       FLOAT     = Array[1434, 121]
; NEVENT          LONG      =         1434
; O3_127          FLOAT     = Array[1434, 121]
; O3_96           FLOAT     = Array[1434, 121]
; PRESSURE        FLOAT     = Array[1434, 121]
; SYYYYDOY        STRING    = Array[1434]
; SYYYYMMDD       STRING    = Array[1434]
; TEMPERATURE     FLOAT     = Array[1434, 121]
; TIME            FLOAT     = Array[1434, 121]
; TPSOLARLT       FLOAT     = Array[1434, 121]
; TPSOLARZEN      FLOAT     = Array[1434, 121]
;
      dum=findfile(dir+sdate+'_v2.0.sav')
      if dum(0) eq '' then dum=findfile('/Volumes/earth/harvey/SABER_data/Datfiles/SABER_TPD_'+sdate+'.sav')
      if dum(0) eq '' then goto,skipsaber
      restore,dum(0)
;
; SABER latitude/longitude are 2d
;
      lat2d=latitude
      lon2d=longitude
      latitude=fltarr(nevent)
      longitude=fltarr(nevent)
      for i=0L,nevent-1L do begin
          latprof=reform(lat2d(i,*))
          lonprof=reform(lon2d(i,*))
          index=where(latprof ne -999.)
          latitude(i)=mean(latprof(index))
          index=where(lonprof ne -999.)
          longitude(i)=mean(lonprof(index))
      endfor
      nlv=n_elements(altitude)
;
; declare time period arrays on first day
;
      if kcount eq 0L then begin
;        print,altitude
         ralt=83.
;        read,'Enter desired altitude ',ralt
         index=where(abs(altitude-ralt) eq min(abs(altitude-ralt)))
         ilev=index(0)
         salt=strcompress(long(ralt),/remove_all)
;        rlat=75.
;        read,'Enter desired latitude',rlat
;        slat=strcompress(long(rlat),/remove_all)
         saberpolartp_xt=fltarr(nc,kday)
         kcount=1
      endif
;
; compute polar temp
;
      saberpolartp=fltarr(nc)
      saberntpprof=lonarr(nc)
      for ii=0L,nevent-1L do begin
          if latitude(ii) ge rlat-2. and latitude(ii) le rlat+2. then begin
;         if latitude(ii) ge 80. then begin
             for i=0L,nc-2L do begin
                 if longitude(ii) ge longrid(i) and longitude(ii) le longrid(i+1) and temperature(ii,ilev) gt 0. then begin
                    saberpolartp(i)=saberpolartp(i)+temperature(ii,ilev)
                    saberntpprof(i)=saberntpprof(i)+1L
                 endif
             endfor
          endif
      endfor
      good=where(saberntpprof gt 0L)
      if good(0) ne -1L then saberpolartp(good)=saberpolartp(good)/float(saberntpprof(good))
      saberpolartp_xt(*,icount)=saberpolartp
skipsaber:
      icount=icount+1L
goto,jump

plotit:
;
; wrap around point
;
saberpolartp_xt(nc-1,*)=saberpolartp_xt(0,*)
;
; interpolate small gaps in time
;
for k=0,nc-1 do begin
    dlev=reform(saberpolartp_xt(k,*))
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
    saberpolartp_xt(k,*)=dlev
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
;
; save temp
;
save,file='xt_saber_temp_'+yearlab+'_'+slat+'_pmc_onset.sav',saberpolartp_xt,kday,longrid,sdate_all
quick:
yearlab=string(format='(i4)',iyear)+'-'+strcompress(iyear+1,/remove_all)
;restore,'xt_saber_temp_'+yearlab+'_'+slat+'_pmc_onset.sav'
if iyear eq 2013 then restore,'xt_saber_temp_2013_'+slat+'_pmc_onset.sav'
sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)
syear=strmid(sdate_all,0,4)
smon=strmid(sdate_all,4,2)
sday=strmid(sdate_all,6,2)
xindex=where(sday eq '01' or sday eq '05' or sday eq '10' or sday eq '15' or sday eq '20' or sday eq '25',nxticks)
;xindex=where(sday eq '01' or sday eq '15',nxticks)
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
          /bold,/color,bits_per_pixel=8,/helvetica,filename='xt_saber_temp+cips_freq_'+yearlab+'_'+slat+'_pmc_onset.ps'
   !p.charsize=1.55
   !p.thick=2
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
;
; plot Arctic mean temperature and CO
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7	; ticks outward
level=[0.001,0.01,0.025,0.05,0.1,0.25,0.5,1.]	;,2.,3.,4.,5.,6.,7.,8.,9.,10.,15.]
index=where(saberpolartp_xt eq 0.)
if index(0) ne -1L then saberpolartp_xt(index)=0./0.
saberpolartp_xt=smooth(saberpolartp_xt,3,/NaN,/edge_truncate)
if index(0) ne -1L then saberpolartp_xt(index)=0./0.

tlevel=140.+2.5*findgen(19)
nlvls=n_elements(tlevel)
col1=1+indgen(nlvls)*icolmax/nlvls
contour,saberpolartp_xt,longrid,findgen(kday),/noeras,yrange=[kday-1,0.],xrange=[0.,360.],$
      charsize=2,charthick=2,color=0,xtitle='Longitude',/cell_fill,c_color=col1,$
      levels=tlevel,yticks=nxticks-1,ytickname=xlabs,ytickv=xindex,min_value=-99.
contour,saberpolartp_xt,longrid,findgen(kday),levels=tlevel,color=0,/follow,/overplot,c_labels=fltarr(nlvls)
contour,saberpolartp_xt,longrid,findgen(kday),levels=[150.],color=mcolor,/follow,/overplot,c_labels=[1],thick=4
axis,yaxis=1,xaxis=1,xrange=[0.,360.],yrange=[-1,-42.],/save,charsize=2,color=0,charthick=2		; 11/10 to 12/20
level=[1,5,10,20,30,40,50,60,70,80,90,100]
;contour,smooth(cipsfreq,5,/NaN,/edge_truncate),lonbin,ddd,levels=level,color=mcolor,/overplot,/noeras,$
;        c_charsize=2,c_labels=1+0*level,thick=10
index=where(iyr eq start_year)
loadct,0
plots,0.,start_date(index)
plots,360.,start_date(index),/continue,color=0,thick=20
xyouts,xmn+0.01,ymx-0.05,yearlab+'/'+slat+'S/'+salt+'km',/normal,charsize=2.5,color=0,charthick=2
loadct,39

imin=min(tlevel)
imax=max(tlevel)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='SABER Temperature (K)',charsize=1.5,charthick=2
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
       spawn,'convert -trim xt_saber_temp+cips_freq_'+yearlab+'_'+slat+'_pmc_onset.ps -rotate -90 xt_saber_temp+cips_freq_'+yearlab+'_'+slat+'_pmc_onset.png'
    endif
endfor  ; loop over latitudes
endfor  ; loop over years


endfor
end
