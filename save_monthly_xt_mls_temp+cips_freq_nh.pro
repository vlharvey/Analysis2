;
; save monthly Hovmollers of MLS temperature, water, and CIPS frequency
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
xorig=[0.15,0.5]
yorig=[0.25,0.25]
xlen=0.325
ylen=0.6
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
mdir='/atmos/aura6/data/MLS_data/Datfiles_SOSST/'
smon=['01','02','03','04','05','06','07','08','09','10','11','12']
kmon=[31,28,31,30,31,30,31,31,30,31,30,31]*1L
nmon=n_elements(smon)
start_year=[2007,2008,2009,2010,2011,2012,2013,2014]
start_date=[-27, -21, -24, -24, -26, -27, -34, -28]
end_date=[66, 65, 61, 61, 64, 61, 64, 80]
nyear=n_elements(start_year)
;
; loop over years
;
for nn=start_year(0),start_year(nyear-1L) do begin

;
; loop over latitudes
;
for ilat=75,75 do begin	;50,80,5 do begin

syr=string(format='(i2.2)',nn-2000)
slt=strcompress(ilat,/remove_all)
restore,'/Volumes/Data/CIPS_data/Pre_process/Save_files/xt_cips_3c_NH'+syr+'_DSC_'+slt+'Lat_2G_freq.sav
decfreq_all=XTFREQ
restore,'/Volumes/Data/CIPS_data/Pre_process/Save_files/xt_cips_3c_NH'+syr+'_ASC_'+slt+'Lat_2G_freq.sav
ascfreq_all=XTFREQ
ddd_all=ddd
cipsfreq_all=0.*decfreq_all
index=where(decfreq_all ne 0. and ascfreq_all eq 0.)
if index(0) ne -1L then cipsfreq_all(index)=decfreq_all(index)
index=where(decfreq_all eq 0. and ascfreq_all ne 0.)
if index(0) ne -1L then cipsfreq_all(index)=ascfreq_all(index)
index=where(decfreq_all ne 0. and ascfreq_all ne 0.)
if index(0) ne -1L then cipsfreq_all(index)=(decfreq_all(index)+ascfreq_all(index))/2.
;
; loop over months
;
for imon=0L,nmon-1L do begin
print,nn,' ',imon+1
;
; compute DFS start and end for this month and truncate cipsfreq, ascfreq, and descfreq accordingly
;
lstmn=imon+1L
lstdy=1L
lstyr=nn*1L
ledmn=imon+1L
leddy=kmon(imon)
ledyr=nn*1L
jday0=julday(lstmn,lstdy,nn)
jday1=julday(ledmn,leddy,nn)
jdaysol=julday(6,21,nn)		; NH
dfsjday0=jday0-jdaysol
dfsjday1=jday1-jdaysol

ddd=0 & decfreq=0 & ascfreq=0 & cipsfreq=0		; zero out from previous month in case this month is zero
index=where(ddd_all ge dfsjday0 and ddd_all le dfsjday1,cday)
if index(0) eq -1L then print,cday,' days of CIPS data this month'
if index(0) ne -1L then begin
ddd=reform(ddd_all(index))
cipsfreq=reform(cipsfreq_all(*,index))
decfreq=reform(decfreq_all(*,index))
ascfreq=reform(ascfreq_all(*,index))
print,dfsjday0,dfsjday1,min(ddd),max(ddd)
endif

kcount=0
rlat=float(ilat)
slat=strcompress(long(rlat),/remove_all)

z = stddat(lstmn,lstdy,nn,lstday)
z = stddat(ledmn,leddy,nn,ledday)
if ledday lt lstday then stop,' Wrong dates! '
kday=ledday-lstday+1L
sdate_all=strarr(kday)
;
; Compute initial Julian date
;
iyr = nn
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
; restore MLS on this day
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
      restore,mdir+'tpd_mls_v3.3_'+sdate+'.sav'
      restore,mdir+'h2o_mls_v3.3_'+sdate+'.sav'
;
; apply mask
;
      index=where(temperature_mask eq -99.)
      if index(0) ne -1L then temperature(index)=-99.
      index=where(mask eq -99.)
      if index(0) ne -1L then mix(index)=-99.
      mlsh2omix=mix
      good=where(mlsh2omix ne -99.)
      if good(0) ne -1L then mlsh2omix(good)=mlsh2omix(good)*1.e6
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
         mlspolartp_xt=fltarr(nc,kday)
         mlspolarh2o_xt=fltarr(nc,kday)
         kcount=1
;goto,quick
      endif
;
; compute polar temp
;
      mlspolartp=fltarr(nc)
      mlsntpprof=lonarr(nc)
      for ii=0L,n_elements(id)-1L do begin
          if latitude(ii) ge rlat-2. and latitude(ii) le rlat+2. then begin
;         if latitude(ii) ge 80. then begin
             for i=0L,nc-2L do begin
                 if longitude(ii) ge longrid(i) and longitude(ii) le longrid(i+1) and temperature(ii,ilev) ne -99. then begin
                    mlspolartp(i)=mlspolartp(i)+temperature(ii,ilev)
                    mlsntpprof(i)=mlsntpprof(i)+1L
                 endif
             endfor
          endif
      endfor
      good=where(mlsntpprof gt 0L)
      if good(0) ne -1L then mlspolartp(good)=mlspolartp(good)/float(mlsntpprof(good))
      mlspolartp_xt(*,icount)=mlspolartp
;
; compute xt h2o
;
      mlspolarh2o=fltarr(nc)
      mlsnh2oprof=lonarr(nc)
      for ii=0L,n_elements(id)-1L do begin
          if latitude(ii) ge rlat-2. and latitude(ii) le rlat+2. then begin
             for i=0L,nc-2L do begin
                 if longitude(ii) ge longrid(i) and longitude(ii) le longrid(i+1) and mlsh2omix(ii,ilev) ne -99. then begin
                    mlspolarh2o(i)=mlspolarh2o(i)+mlsh2omix(ii,ilev)
                    mlsnh2oprof(i)=mlsnh2oprof(i)+1L
                 endif
             endfor
          endif
      endfor
      good=where(mlsnh2oprof gt 0L)
      if good(0) ne -1L then mlspolarh2o(good)=mlspolarh2o(good)/float(mlsnh2oprof(good))
      mlspolarh2o_xt(*,icount)=mlspolarh2o

skipmls:
      icount=icount+1L
goto,jump

plotit:
;
; wrap around point
;
mlspolartp_xt(nc-1,*)=mlspolartp_xt(0,*)
mlspolarh2o_xt(nc-1,*)=mlspolarh2o_xt(0,*)
;
; interpolate small gaps in time
;
for k=0,nc-1 do begin
    dlev=reform(mlspolartp_xt(k,*))
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
    mlspolartp_xt(k,*)=dlev
endfor
for k=0,nc-1 do begin
    dlev=reform(mlspolarh2o_xt(k,*))
    for i=1,kday-1 do begin
        if dlev(i) eq 0. and dlev(i-1) ne 0. then begin
           for ii=i+1,kday-1 do begin
               naway=float(ii-i)
               if naway le 5.0 and dlev(ii) ne 0. then begin
                  dlev(i)=(naway*dlev(i-1)+dlev(ii))/(naway+1.0)
                  goto,jump3
               endif
           endfor
jump3:
        endif
    endfor
    mlspolarh2o_xt(k,*)=dlev
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
save,file='xt_mls_temp_h2o_freq_'+yearlab+smon(imon)+'_'+slat+'.sav',mlspolartp_xt,mlspolarh2o_xt,kday,longrid,sdate_all,cipsfreq,decfreq,ascfreq,lonbin,ddd,ralt,ilat
quick:
yearlab=string(format='(i4)',nn)
restore,'xt_mls_temp_h2o_freq_'+yearlab+smon(imon)+'_'+slat+'.sav'
sdate0=sdate_all(0)
sdate1=sdate_all(n_elements(sdate_all)-1)
syear=strmid(sdate_all,0,4)
ssmon=strmid(sdate_all,4,2)
sday=strmid(sdate_all,6,2)
xindex=where(sday eq '01' or sday eq '05' or sday eq '10' or sday eq '15' or sday eq '20' or sday eq '25' or sday eq '30',nxticks)
;xindex=where(sday eq '01' or sday eq '15',nxticks)
xlabs=ssmon(xindex)+'/'+sday(xindex)
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
          /bold,/color,bits_per_pixel=8,/helvetica,filename='xt_mls_temp+cips_freq_'+yearlab+smon(imon)+'_'+slat+'.ps'
   !p.charsize=1.55
   !p.thick=2
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
;
; plot temperature, water, clouds
;
erase
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7	; ticks outward
level=[0.001,0.01,0.025,0.05,0.1,0.25,0.5,1.]	;,2.,3.,4.,5.,6.,7.,8.,9.,10.,15.]
index=where(mlspolartp_xt eq 0.)
if index(0) ne -1L then mlspolartp_xt(index)=0./0.
mlspolartp_xt=smooth(mlspolartp_xt,3,/NaN,/edge_truncate)
if index(0) ne -1L then mlspolartp_xt(index)=0./0.

tlevel=130.+5.*findgen(21)	;2.5*findgen(19)
if imon ge 4 and imon lt 9 then tlevel=130.+2.5*findgen(19)		; MJJA temps
nlvls=n_elements(tlevel)
col1=1+indgen(nlvls)*icolmax/nlvls
contour,mlspolartp_xt,longrid,findgen(kday),/noeras,yrange=[kday-1,0.],xrange=[0.,360.],$
      charsize=2,color=0,xtitle='Longitude',/cell_fill,c_color=col1,charthick=2,$
      levels=tlevel,yticks=nxticks-1,ytickname=xlabs,ytickv=xindex,min_value=-99.
contour,mlspolartp_xt,longrid,findgen(kday),levels=tlevel,color=0,/follow,/overplot,c_labels=fltarr(nlvls)
print,min(mlspolartp_xt),max(mlspolartp_xt)
loadct,0
;contour,mlspolartp_xt,longrid,findgen(kday),levels=[150.],color=0.5*mcolor,/follow,/overplot,c_labels=[1],thick=4
loadct,39
axis,yaxis=1,xaxis=1,xrange=[0.,360.],yrange=[dfsjday1,dfsjday0],/save,charsize=2,color=0,charthick=2,yticks=1,ytickname=[' ',' ']
level=[1,5,10,20,30,40,50,60,70,80,90,100]

index=where(iyr eq start_year)
loadct,0
plots,0.,start_date(index)
plots,360.,start_date(index),/continue,color=0,thick=20
xyouts,0.4,0.9,yearlab+' '+smon(imon),/normal,charsize=3,color=0,charthick=2
loadct,39

if cday ge 2 then begin
contour,smooth(cipsfreq,5,/NaN,/edge_truncate),lonbin,ddd,levels=[1],color=0,/overplot,/noeras,$
        c_charsize=2,c_labels=[1],thick=10
loadct,0
contour,smooth(cipsfreq,5,/NaN,/edge_truncate),lonbin,ddd,levels=[10],color=0.1*mcolor,/overplot,/noeras,c_charsize=2,c_labels=[1],thick=10
contour,smooth(cipsfreq,5,/NaN,/edge_truncate),lonbin,ddd,levels=[30],color=0.3*mcolor,/overplot,/noeras,c_charsize=2,c_labels=[1],thick=10
contour,smooth(cipsfreq,5,/NaN,/edge_truncate),lonbin,ddd,levels=[50],color=0.5*mcolor,/overplot,/noeras,c_charsize=2,c_labels=[1],thick=10
contour,smooth(cipsfreq,5,/NaN,/edge_truncate),lonbin,ddd,levels=[70],color=0.7*mcolor,/overplot,/noeras,c_charsize=2,c_labels=[1],thick=10
contour,smooth(cipsfreq,5,/NaN,/edge_truncate),lonbin,ddd,levels=[90],color=0.9*mcolor,/overplot,/noeras,c_charsize=2,c_labels=[1],thick=10
contour,smooth(cipsfreq,5,/NaN,/edge_truncate),lonbin,ddd,levels=[100],color=mcolor,/overplot,/noeras,c_charsize=2,c_labels=[1],thick=10
loadct,39
endif

imin=min(tlevel)
imax=max(tlevel)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='MLS Temperature (K)',charsize=1.5,charthick=2
ybox=[0,10,10,0,0]
x1=imin
dx=(imax-imin)/float(nlvls)
for jj=0,nlvls-1 do begin
xbox=[x1,x1,x1+dx,x1+dx,x1]
polyfill,xbox,ybox,color=col1(jj)
x1=x1+dx
endfor

xmn=xorig(1)
xmx=xorig(1)+xlen
ymn=yorig(1)
ymx=yorig(1)+ylen
set_viewport,xmn,xmx,ymn,ymx
!type=2^2+2^3+2^7       ; ticks outward
index=where(mlspolarh2o_xt eq 0.)
if index(0) ne -1L then mlspolarh2o_xt(index)=0./0.
mlspolarh2o_xt=smooth(mlspolarh2o_xt,3,/NaN,/edge_truncate)
if index(0) ne -1L then mlspolarh2o_xt(index)=0./0.

tlevel=0.5+.25*findgen(27)
nlvls=n_elements(tlevel)
col1=1+indgen(nlvls)*icolmax/nlvls
contour,mlspolarh2o_xt,longrid,findgen(kday),/noeras,yrange=[kday-1,0.],xrange=[0.,360.],$
      charsize=2,color=0,xtitle='Longitude',/cell_fill,c_color=col1,charthick=2,$
      levels=tlevel,yticks=nxticks-1,ytickname=' '+strarr(nxticks),ytickv=xindex,min_value=-99.
contour,mlspolarh2o_xt,longrid,findgen(kday),levels=tlevel,color=mcolor,/follow,/overplot,c_labels=fltarr(nlvls)
print,min(mlspolarh2o_xt),max(mlspolarh2o_xt)
axis,yaxis=1,xaxis=1,xrange=[0.,360.],yrange=[dfsjday1,dfsjday0],/save,charsize=2,color=0,charthick=2,ytitle='DFS'
level=[1,5,10,20,30,40,50,60,70,80,90,100]

index=where(iyr eq start_year)
loadct,0
plots,0.,start_date(index)
plots,360.,start_date(index),/continue,color=0,thick=20
loadct,39

if cday ge 2 then begin
contour,smooth(cipsfreq,5,/NaN,/edge_truncate),lonbin,ddd,levels=[1],color=0,/overplot,/noeras,$
        c_charsize=2,c_labels=[1],thick=10
loadct,0
contour,smooth(cipsfreq,5,/NaN,/edge_truncate),lonbin,ddd,levels=[10],color=0.1*mcolor,/overplot,/noeras,c_charsize=2,c_labels=[1],thick=10
contour,smooth(cipsfreq,5,/NaN,/edge_truncate),lonbin,ddd,levels=[30],color=0.3*mcolor,/overplot,/noeras,c_charsize=2,c_labels=[1],thick=10
contour,smooth(cipsfreq,5,/NaN,/edge_truncate),lonbin,ddd,levels=[50],color=0.5*mcolor,/overplot,/noeras,c_charsize=2,c_labels=[1],thick=10
contour,smooth(cipsfreq,5,/NaN,/edge_truncate),lonbin,ddd,levels=[70],color=0.7*mcolor,/overplot,/noeras,c_charsize=2,c_labels=[1],thick=10
contour,smooth(cipsfreq,5,/NaN,/edge_truncate),lonbin,ddd,levels=[90],color=0.9*mcolor,/overplot,/noeras,c_charsize=2,c_labels=[1],thick=10
contour,smooth(cipsfreq,5,/NaN,/edge_truncate),lonbin,ddd,levels=[100],color=mcolor,/overplot,/noeras,c_charsize=2,c_labels=[1],thick=10
loadct,39
endif

imin=min(tlevel)
imax=max(tlevel)
ymnb=yorig(0) -cbaryoff
ymxb=ymnb  +cbarydel
set_viewport,xmn,xmx,ymnb,ymxb
!type=2^2+2^3+2^6
plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,xtitle='H!l2!nO (ppmv) '+slat+'N/'+salt+'km',charsize=1.5,charthick=2
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
       spawn,'convert -trim xt_mls_temp+cips_freq_'+yearlab+smon(imon)+'_'+slat+'.ps -rotate -90 xt_mls_temp+cips_freq_'+yearlab+smon(imon)+'_'+slat+'.png'
    endif
endfor  ; loop over months
endfor  ; loop over latitudes
endfor  ; loop over years
end
