;
; x-axis is solar flux. y-axis is equatorial zonal wind
; separate QBO phases based on radiosondes
; SSWs based on ERA and MetO
;
; VLH 9/24/09
;
@stddat
@kgmt
@ckday
@kdate
@rd_era40_nc

loadct,39
device,decompose=0
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
icmm1=icolmax-1
icmm2=icolmax-2
a=findgen(9)*(2*!pi/8.)
usersym,1.5*cos(a),1.5*sin(a),/fill
!noeras=1
nxdim=700
nydim=700
xorig=[0.15]
yorig=[0.25]
xlen=0.7
ylen=0.6
cbaryoff=0.08
cbarydel=0.02
setplot='ps'
read,'setplot=',setplot
if setplot ne 'ps' then begin
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
   !p.background=mcolor
endif
mon=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
smon=['01','02','03','04','05','06','07','08','09','10','11','12']
nmonth=n_elements(mon)
;
; read f10.7 monthly mean data
;
close,1
openr,1,'f10.7cm_solar_flux_monthly_data.txt'
dum=' '
for i=0,5 do readf,1,dum
iyear=0
nyear=63
years=lonarr(nyear)
f10_jan=fltarr(nyear)
f10_feb=fltarr(nyear)
f10_mar=fltarr(nyear)
f10_apr=fltarr(nyear)
f10_may=fltarr(nyear)
f10_jun=fltarr(nyear)
f10_jul=fltarr(nyear)
f10_aug=fltarr(nyear)
f10_sep=fltarr(nyear)
f10_oct=fltarr(nyear)
f10_nov=fltarr(nyear)
f10_dec=fltarr(nyear)
while not eof(1) do begin
      readf,1,dum
      vals=strsplit(dum,' ',/extract)
;
; save years and monthly average f10.7 each year
;
      years(iyear)=long(vals(0))
      f10_jan(iyear)=vals(1)/10.
      f10_feb(iyear)=vals(2)/10.
      f10_mar(iyear)=vals(3)/10.
      f10_apr(iyear)=vals(4)/10.
      f10_may(iyear)=vals(5)/10.
      f10_jun(iyear)=vals(6)/10.
      f10_jul(iyear)=vals(7)/10.
      f10_aug(iyear)=vals(8)/10.
      f10_sep(iyear)=vals(9)/10.
      f10_oct(iyear)=vals(10)/10.
      f10_nov(iyear)=vals(11)/10.
      f10_dec(iyear)=vals(12)/10.
      iyear=iyear+1L
endwhile
close,1
;
; read radiosonde wind for QBO proxy
;
close,1
openr,1,'/aura2/harvey/Qbo/qbo_data'
dummy=' '
station_id=0l
nlev=7L
mb70=0l & mb50=0l & mb40=0l & mb30=0l & mb20=0l & mb15=0l & mb10=0l
num=0L & yymm1=0L & nlines=1000L
yymm=9999+lonarr(nlines)
ftime=9999+fltarr(nlines)
stime=strarr(nlines)
u=-9999+fltarr(nlines,nlev)
p=[70.,50.,40.,30.,20.,15.,10.]
readf,1, dummy
readf,1, dummy
readf,1, dummy
readf,1, dummy
readf,1, dummy
while not eof(1) do begin
      readf,1, station_id, dummy, yymm1, dummy, mb70, dummy, mb50, dummy, $
         mb40, dummy, mb30, dummy, mb20, dummy, mb15, dummy, $
         mb10, dummy, format='(I5,A1,I6,A1,7(I5,A2))'
      yymm(num)=yymm1
      u(num,0)=float(mb70)/10.
      u(num,1)=float(mb50)/10.
      u(num,2)=float(mb40)/10.
      u(num,3)=float(mb30)/10.
      u(num,4)=float(mb20)/10.
      u(num,5)=float(mb15)/10.
      u(num,6)=float(mb10)/10.
      ftime(num)=float(yymm1)
      stime(num)=strcompress(string(yymm(num)),/remove_all)
      num=num+1
endwhile
close,1
yymm=yymm(0:num-1)
ftime=ftime(0:num-1)
stime=stime(0:num-1)
u=u(0:num-1,*)
bad=where(u eq 99.9)
if bad(0) ne -1 then u(bad)=-9999.
index=where(p eq 50.)
up=reform(u(*,index(0)))
radyear=strmid(strcompress(yymm,/remove_all),0,4)
radmon=strmid(strcompress(yymm,/remove_all),4,2)
slev=strcompress(long(p(index(0))),/remove_all)+'hPa'
;
; save postscript
;
if setplot eq 'ps' then begin
   set_plot,'ps'
   xsize=nxdim/100.
   ysize=nydim/100.
   !p.font=0
   device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
          /bold,/color,bits_per_pixel=8,/helvetica,filename='solar_flux_vs_singapore_winds+ssw_'+slev+'_antic.ps'
   !p.charsize=1.25
   !p.thick=2
   !p.charthick=5
   !y.thick=2
   !x.thick=2
endif
;
; restore MetO marker zonal mean counts
;
restore,'MetO_antic_area_alllat_1991-2009.sav
syyyymmdd_all=strcompress(yyyymmdd_all,/remove_all)
syear=strmid(syyyymmdd_all,0,4)
smon=strmid(syyyymmdd_all,4,2)
;
; compute Feb mean anticyclone counts > 40 N 
ilev50=where(th eq 500.)
ilev50=ilev50(0)
ilev30=where(th eq 600.)
ilev30=ilev30(0)
ilev10=where(th eq 700.)
ilev10=ilev10(0)
ilev1=where(th eq 2000.)
ilev1=ilev1(0)
ilat=where(alat eq 41.2500)
ilat=ilat(0)
markavg=0.*fltarr(nyear)
markavg2d=0.*fltarr(nyear,nth)
for i=44,nyear-1L do begin
;
; sum latitudes poleward ilat
for j=ilat,nr-1 do begin
    marklev=reform(markBAR_ALL(*,j,ilev50))
    index=where(marklev ne -9999. and syear eq strcompress(years(i),/remove_all) and smon eq '02',nn)
    if index(0) ne -1 then begin        ; there are years ini f10.7 record before and after ERA40
       markavg(i)=markavg(i)+total(marklev(index))
    endif
;
; function of altitude
;
;    for k=0,nth-1L do begin
;    marklev=reform(markBAR_ALL(*,j,k))
;    index=where(marklev ne -9999. and syear eq strcompress(years(i),/remove_all) and smon eq '11',nn)
;;               (smon eq '09' or smon eq '10'),nn)
;    if index(0) ne -1 then begin        ; there are years ini f10.7 record before and after ERA40
;       markavg2d(i,k)=markavg2d(i,k)+total(marklev(index))
;    endif
;    endfor
;
endfor
print,'N anticyclones ',years(i),markavg(i)
endfor
;erase
;;index=where(markavg ne 0.)
;;plot,years(index),markavg(index),psym=1,color=0
;;plot,years(index),markavg(index),color=0
;      imin=10.
;      imax=max(markavg2d)
;nlvls=20
;level=imin+((imax-imin)/float(nlvls))*findgen(nlvls-1)
;        col1=1+indgen(nlvls)*mcolor/nlvls
;index=where(years ge 1991,nn)
;contour,markavg2d(index,*),years(index),th,yrange=[500.,2000.],levels=level,c_color=col1,/fill,$
;        color=0,/noeras,ytitle='Theta',title='# Anticyclone Points',xticks=nn-1,xtickv=years(index),$
;        xtickname=strcompress(years(index),/remove_all)
;contour,markavg2d(index,*),years(index),th,levels=level,color=0,/follow,/overplot,/noeras
;
;stop
;
; restore MetO data
;
restore,'MetO_tbar_ubar_alllat_1991-2009.sav
syyyymmdd_all=strcompress(yyyymmdd_all,/remove_all)
syear=strmid(syyyymmdd_all,0,4)
smon=strmid(syyyymmdd_all,4,2)
;
; interpolate MetO to ERA40 pressure surfaces
;
epress=[1000.,925.,850.,775.,700.,600.,500.,400.,300.,250.,200.,150.,$
        100.,70.,50.,30.,20.,10.,7.,5.,3.,2.,1.]
enl=n_elements(epress)
nday=n_elements(YYYYMMDD_ALL)
gbar_all_int=fltarr(nday,nr,enl)
tbar_all_int=fltarr(nday,nr,enl)
ubar_all_int=fltarr(nday,nr,enl)
for k=0L,enl-1L do begin
    p0=epress(k)
    for kk=0L,nl-2 do begin
        if press(kk) eq p0 then begin
           gbar_all_int(*,*,k)=gbar_all(*,*,kk)
           tbar_all_int(*,*,k)=tbar_all(*,*,kk)
           ubar_all_int(*,*,k)=ubar_all(*,*,kk)
        endif
        if press(kk) ne p0 then begin
        if press(kk) gt p0 and press(kk+1) lt p0 then begin
           zscale=(p0-press(kk+1))/(press(kk)-press(kk+1))
           gbar_all_int(*,*,k)=gbar_all(*,*,kk+1)+zscale*(gbar_all(*,*,kk)-gbar_all(*,*,kk+1))
           tbar_all_int(*,*,k)=tbar_all(*,*,kk+1)+zscale*(tbar_all(*,*,kk)-tbar_all(*,*,kk+1))
           ubar_all_int(*,*,k)=ubar_all(*,*,kk+1)+zscale*(ubar_all(*,*,kk)-ubar_all(*,*,kk+1))
        endif
        endif
    endfor
endfor
gbar_all=gbar_all_int
tbar_all=tbar_all_int
ubar_all=ubar_all_int
press=epress
;
; compute february monthly averages at 30 hPa
;
ilev50=where(press eq 50.)
ilev50=ilev50(0)
ilev30=where(press eq 30.)
ilev30=ilev30(0)
ilev20=where(press eq 20.)
ilev20=ilev20(0)
ilev10=where(press eq 10.)
ilev10=ilev10(0)
ilev1=where(press eq 1.)
ilev1=ilev1(0)
mgplev30=reform(GBAR_ALL(*,nr-1,ilev30)) ; daily GPH at NP and 30 hPa
mtplev30=reform(TBAR_ALL(*,nr-1,ilev30)) ; daily temp at NP and 30 hPa
ilat=where(alat eq 60.)
ilat=ilat(0)
ilat0=where(alat eq 0.)
ilat0=ilat0(0)
mulev60N10=reform(UBAR_ALL(*,ilat,ilev10))       ; daily Ubar at 60N and 10 hPa
mulev60N1=reform(UBAR_ALL(*,ilat,ilev1)) ; daily Ubar at 60N and 1 hPa
mulev30=reform(UBAR_ALL(*,ilat0,ilev30))    ; daily Ubar at Equator and 30 hPa
mulev50=reform(UBAR_ALL(*,ilat0,ilev50))    ; daily Ubar at Equator and 50 hPa
mulev10=reform(UBAR_ALL(*,ilat0,ilev10))    ; daily Ubar at Equator and 10 hPa
mulev20=reform(UBAR_ALL(*,ilat0,ilev20))    ; daily Ubar at Equator and 10 hPa
mgpavg30=-99.+0.*fltarr(nyear)                   ; Feb monthly average GPH at the NP and 30 hPa
mtpavg30=-99.+0.*fltarr(nyear)                   ; Feb monthly average temp at the NP and 30 hPa
muavgeq30=-99.+0.*fltarr(nyear)                  ; Feb monthly average zonal wind at the equator
muavgeq50=-99.+0.*fltarr(nyear)                  ; Feb monthly average zonal wind at the equator
muavgeq10=-99.+0.*fltarr(nyear)                  ; Feb monthly average zonal wind at the equator
muavgeq20=-99.+0.*fltarr(nyear)                  ; Feb monthly average zonal wind at the equator
musigeq50=-99.+0.*fltarr(nyear)                  ; Feb monthly standard deviation of zonal wind at the equator
muchgeq50=-99.+0.*fltarr(nyear)                  ; Feb monthly standard deviation of zonal wind at the equator
muchgeq30=-99.+0.*fltarr(nyear)                  ; Feb monthly standard deviation of zonal wind at the equator
muchgeq10=-99.+0.*fltarr(nyear)                  ; Feb monthly standard deviation of zonal wind at the equator
muavg60N10=-99.+0.*fltarr(nyear)                 ; Feb monthly average zonal wind at 60 N
muavg60N1=-99.+0.*fltarr(nyear)                  ; Feb monthly average zonal wind at 60 N
urad=-99.+0.*fltarr(nyear)                  ; Feb monthly average zonal wind at radiosonde locations
nssw=-99.+0.*fltarr(nyear)
for i=0,nyear-1L do begin
    index=where(syear eq strcompress(years(i),/remove_all) and smon eq '02',nn)
    index2=where(radyear eq strcompress(years(i),/remove_all) and radmon eq '02')
    if index(0) ne -1 then begin        ; there are years ini f10.7 record before and after ERA40
       mtpavg30(i)=total(mtplev30(index))/float(nn)
       mgpavg30(i)=total(mgplev30(index))/float(nn)
       muavgeq30(i)=total(mulev30(index))/float(nn)
       muavgeq50(i)=total(mulev50(index))/float(nn)
       muavgeq10(i)=total(mulev10(index))/float(nn)
       muavgeq20(i)=total(mulev20(index))/float(nn)
       musigeq50(i)=stdev(mulev50(index))
       if min(mulev50(index)) lt 0. and max(mulev50(index)) gt 0. then muchgeq50(i)=1
       if min(mulev10(index)) lt 0. and max(mulev10(index)) gt 0. then muchgeq10(i)=1
       if min(mulev30(index)) lt 0. and max(mulev30(index)) gt 0. then muchgeq30(i)=1
       muavg60N1(i)=min(mulev60N1(index))
       muavg60N10(i)=min(mulev60N10(index))

urad(i)=up(index2(0))
uvals=reform(mulev60N10(index))
mmwindex=where(uvals lt 0.,nmmw)
nssw(i)=nmmw
print,years(i),mtpavg30(i),urad(i),muavgeq50(i),musigeq50(i),muchgeq50(i),muavg60N10(i),nmmw
    endif
endfor
;
; restore ERA40 files
;
; ALAT            FLOAT     = Array[73]
; COMMENT         STRING    = 'dT/dy=Tbar85-Tbar60 and Ubar65 is zonal mean zonal wind at 65N/S'
; NDAY            LONG      =        16436
; NH_DTDY_ALL     FLOAT     = Array[16436, 23]
; NH_UBAR65_ALL   FLOAT     = Array[16436, 23]
; NL              LONG      =           23
; NR              LONG      =           73
; PRESS           FLOAT     = Array[23]
; SH_DTDY_ALL     FLOAT     = Array[16436, 23]
; SH_UBAR65_ALL   FLOAT     = Array[16436, 23]
; TBAR_ALL        FLOAT     = Array[16436, 73, 23]
; UBAR_ALL        FLOAT     = Array[16436, 73, 23]
; UBAR_EQ_ALL     FLOAT     = Array[16436, 23]
; YYYYMMDD_ALL    LONG      = Array[16436]
;
restore,'ERA40_dTdy_Ubar_QBO_1957-2002.sav'
restore,file='ERA40_tbar_ubar_alllat_era40_1957-2002.sav'	; no alat?
nr=73
alat=-90.+2.5*findgen(nr)
syyyymmdd_all=strcompress(yyyymmdd_all,/remove_all)
syear=strmid(syyyymmdd_all,0,4)
smon=strmid(syyyymmdd_all,4,2)
;
; compute february monthly averages at 30 hPa
;
;     years(iyear)=long(vals(0))
;     f10_feb(iyear)=vals(2)/10.
; TBAR_ALL        FLOAT     = Array[16436, 73, 23]
; GBAR_ALL        FLOAT     = Array[16436, 73, 23]
; UBAR_EQ_ALL     FLOAT     = Array[16436, 23]
; YYYYMMDD_ALL    LONG      = Array[16436]

ilev50=where(press eq 50.)
ilev50=ilev50(0)
ilev30=where(press eq 30.)
ilev30=ilev30(0)
ilev10=where(press eq 10.)
ilev10=ilev10(0)
ilev1=where(press eq 1.)
ilev1=ilev1(0)
gplev30=reform(GBAR_ALL(*,nr-1,ilev30))	; daily GPH at NP and 30 hPa
tplev30=reform(TBAR_ALL(*,nr-1,ilev30))	; daily temp at NP and 30 hPa
ilat=where(alat eq 60.)
ilat=ilat(0)
ulev60N10=reform(UBAR_ALL(*,ilat,ilev10))	; daily Ubar at 60N and 10 hPa
ulev60N1=reform(UBAR_ALL(*,ilat,ilev1))	; daily Ubar at 60N and 1 hPa
ulev30=reform(UBAR_EQ_ALL(*,ilev30))	; daily Ubar at Equator and 30 hPa
ulev50=reform(UBAR_EQ_ALL(*,ilev50))	; daily Ubar at Equator and 50 hPa
ulev10=reform(UBAR_EQ_ALL(*,ilev10))	; daily Ubar at Equator and 10 hPa
gpavg30=-99.+0.*fltarr(nyear)			; Feb monthly average GPH at the NP and 30 hPa 
tpavg30=-99.+0.*fltarr(nyear)			; Feb monthly average temp at the NP and 30 hPa 
uavgeq30=-99.+0.*fltarr(nyear)			; Feb monthly average zonal wind at the equator
uavgeq50=-99.+0.*fltarr(nyear)			; Feb monthly average zonal wind at the equator
uavgeq10=-99.+0.*fltarr(nyear)			; Feb monthly average zonal wind at the equator
usigeq50=-99.+0.*fltarr(nyear)			; Feb monthly standard deviation of zonal wind at the equator
uchgeq50=-99.+0.*fltarr(nyear)			; Feb monthly standard deviation of zonal wind at the equator
uchgeq30=-99.+0.*fltarr(nyear)			; Feb monthly standard deviation of zonal wind at the equator
uchgeq10=-99.+0.*fltarr(nyear)			; Feb monthly standard deviation of zonal wind at the equator
uavg60N10=-99.+0.*fltarr(nyear)			; Feb monthly average zonal wind at 60 N
uavg60N1=-99.+0.*fltarr(nyear)			; Feb monthly average zonal wind at 60 N
for i=0,nyear-1L do begin
    index=where(syear eq strcompress(years(i),/remove_all) and smon eq '02',nn)
    index2=where(radyear eq strcompress(years(i),/remove_all) and radmon eq '02')
    if index(0) ne -1 then begin	; there are years ini f10.7 record before and after ERA40
       tpavg30(i)=total(tplev30(index))/float(nn)
       gpavg30(i)=total(gplev30(index))/float(nn)
       uavgeq30(i)=total(ulev30(index))/float(nn)
       uavgeq50(i)=total(ulev50(index))/float(nn)
       uavgeq10(i)=total(ulev10(index))/float(nn)
       usigeq50(i)=stdev(ulev50(index))
       if min(ulev50(index)) lt 0. and max(ulev50(index)) gt 0. then uchgeq50(i)=1
       if min(ulev10(index)) lt 0. and max(ulev10(index)) gt 0. then uchgeq10(i)=1
       if min(ulev30(index)) lt 0. and max(ulev30(index)) gt 0. then uchgeq30(i)=1
       uavg60N1(i)=min(ulev60N1(index))
       uavg60N10(i)=min(ulev60N10(index))

urad(i)=up(index2(0))

uvals=reform(ulev60N10(index))
mmwindex=where(uvals lt 0.,nmmw)
nssw(i)=nmmw
print,years(i),tpavg30(i),urad(i),uavgeq50(i),usigeq50(i),uchgeq50(i),uavg60N10(i),nmmw
    endif
endfor
;
; plot
;
erase
!type=2^2+2^3+2^7
xmn=xorig(0)
xmx=xorig(0)+xlen
ymn=yorig(0)
ymx=yorig(0)+ylen
set_viewport,xmn,xmx,ymn,ymx
tmin=-30.
tmax=20.
tpavg30=gpavg30/10000.
index=where(tpavg30 ne -99.,nn)
plot,f10_feb(index),urad(index),xtitle='Solar Flux 10.7 cm',ytitle='Equatorial Wind at '+slev,psym=8,color=0,$
     xrange=[60.,260.],yrange=[tmin,tmax]
plots,60.,0.
plots,260.,0.,/continue,color=0
plots,140.,tmin
plots,140.,tmax,/continue,color=0
for i=0L,nn-1L do xyouts,f10_feb(index(i)),urad(index(i)),strmid(strcompress(years(index(i)),/remove_all),2,2),$
    charsize=1,/data,color=0
;index=where(tpavg30 ne -99. and uavg60N10 lt 0.,nn)
;oplot,f10_feb(index),urad(index),psym=8,color=mcolor*.9
;for i=0L,nn-1L do xyouts,f10_feb(index(i)),urad(index(i)),strmid(strcompress(years(index(i)),/remove_all),2,2),$
;    charsize=1,/data,color=mcolor*.9
;index=where(mtpavg30 ne -99. and muavg60N10 lt 0. and years gt 2002,nn)
;oplot,f10_feb(index),urad(index),psym=8,color=mcolor*.9
;for i=0L,nn-1L do xyouts,f10_feb(index(i)),urad(index(i)),strmid(strcompress(years(index(i)),/remove_all),2,2),$
;    charsize=1,/data,color=mcolor*.9
;
; 5 or more days of major warming 
;
index=where(tpavg30 ne -99. and nssw ge 2.,nn)
oplot,f10_feb(index),urad(index),psym=8,color=mcolor*.9
for i=0L,nn-1L do xyouts,f10_feb(index(i)),urad(index(i)),strmid(strcompress(years(index(i)),/remove_all),2,2),$
    charsize=1,/data,color=mcolor*.9
;
; color points by 10 hPa zonal mean wind at 60 N
;
        level=-30.+10.*findgen(9)
      imin=min(level)
      imax=max(level)
        nlvls=n_elements(level)
        col1=1+indgen(nlvls)*mcolor/nlvls
muavg60n10=muavg60n1
uavg60n10=uavg60n1
loadct,0
for i=0L,nyear-1L do if muavg60n10(i) ne -99. then oplot,[f10_feb(i),f10_feb(i)],$
    [urad(i),urad(i)],psym=8,color=170	;((muavg60N10(i)-imin)/(imax-imin))*mcolor
for i=0L,nyear-1L do if uavg60n10(i) ne -99. then oplot,[f10_feb(i),f10_feb(i)],$
    [urad(i),urad(i)],psym=8,color=170	;((uavg60N10(i)-imin)/(imax-imin))*mcolor
loadct,39
;
; color points by number of anticyclone points
;
imin=100.
imax=1200.
for i=0L,nyear-1L do if uavg60n10(i) ne -99. and markavg(i) ne 0. then oplot,[f10_feb(i),f10_feb(i)],$
    [urad(i),urad(i)],psym=8,color=((markavg(i)-imin)/(imax-imin))*mcolor

      ymnb=yorig(0)-cbaryoff
      ymxb=ymnb+cbarydel
      set_viewport,xorig(0)+0.1,xorig(0)+xlen-0.1,ymnb,ymxb
      !type=2^2+2^3+2^6
      plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],color=0,$
           xtitle='# Anticyclone points >40N at 500 K'
      ybox=[0,10,10,0,0]
      x1=imin
      dx=(imax-imin)/float(nlvls)
      for j=0,nlvls-1 do begin
          xbox=[x1,x1,x1+dx,x1+dx,x1]
          polyfill,xbox,ybox,color=col1(j)
          x1=x1+dx
      endfor

;
; save jpg
;
if setplot eq 'ps' then begin
device,/close
spawn,'convert -trim solar_flux_vs_singapore_winds+ssw_'+slev+'_antic.ps -rotate -90 solar_flux_vs_singapore_winds+ssw_'+slev+'_antic.jpg'
endif
end
