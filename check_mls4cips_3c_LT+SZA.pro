;
; check
; instead of separating ASC/DES by SZA, use LT
; save MLS T and H2O for each CIPS level 3c summary file
; VLH 10/10/2011
;
@stddat
@kgmt
@ckday
@kdate
@mkltime

re=40000./2./!pi
rad=double(180./!pi)
dtr=double(!pi/180.)

loadct,39
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
device,decompose=0
!p.background=icolmax
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
setplot='ps'
read,'setplot=',setplot
nxdim=750
nydim=750
xorig=[0.15,0.15,0.55,0.55]
yorig=[0.55,0.15,0.55,0.15]
xlen=0.225
ylen=0.3
cbaryoff=0.02
cbarydel=0.01
!NOERAS=-1
;if setplot ne 'ps' then window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smonth=['J','F','M','A','M','J','J','A','S','O','N','D']
mdir='/atmos/aura6/data/MLS_data/Datfiles_SOSST/'
cdir='/atmos/harvey/CIPS_data/Datfiles/Level_3c_Summary/cips_3c_'
odir='/Volumes/Data/CIPS_data/Datfiles_MLS_DMP/cips_3c_mls_'
version='v04.20_r05'
;version='v05.10_r01'

    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !p.font=0
       device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
              /bold,/color,bits_per_pixel=8,/helvetica,filename='check_mls4cips_3c_LT_SZA_'+version+'.ps'
       !p.charsize=1.25
       !p.thick=2
       !p.charthick=5
       !p.charthick=5
       !y.thick=2
       !x.thick=2
    endif
;
; Ask interactive questions- get starting/ending date and p surface
;
for ihem=0,1 do begin
if ihem eq 0 then shem='n'
if ihem eq 1 then shem='s'
lstyr=15
ialb=2
;read,' Enter hemisphere (n,s) ',shem
;read,' Enter PMC onset year (yy) ',lstyr
;read,' Enter Albedo threshhold (1,2,5,10) ',ialb
salb=string(format='(i2.2)',ialb)+'G'
if shem eq 'n' then begin
   shem='north'
   syr=string(format='(i2.2)',lstyr)
endif
if shem eq 's' then begin
   shem='south'
   syr=string(format='(i2.2)',lstyr)+string(format='(i2.2)',lstyr+1)
endif
;
; restore CIPS season
;
; ALB             FLOAT     = Array[589, 70]
; ALB_STD         FLOAT     = Array[589, 70]
; DFS             LONG      = Array[589]
; DOY             INT       = Array[589]
; IWC             FLOAT     = Array[589, 70]
; IWC_STD         FLOAT     = Array[589, 70]
; LATHI           INT       = Array[70]
; LATLO           INT       = Array[70]
; LON             FLOAT     = Array[589, 70]
; LTIME           FLOAT     = Array[589, 70]
; NBIN            INT       =       70
; NREV            LONG      =          589
; NUM_CLD         INT       = Array[589, 70]
; NUM_OBS         INT       = Array[589, 70]
; RAD             FLOAT     = Array[589, 70]
; RAD_STD         FLOAT     = Array[589, 70]
; REV             INT       = Array[589]
; SZA             FLOAT     = Array[589, 70]
; UT              FLOAT     = Array[589, 70]
; YEAR            INT       = Array[589]
;
restore,cdir+shem+'_'+syr+'_'+version+'_'+salb+'_all.sav'
print,'restored '+cdir+shem+'_'+syr+'_'+version+'_'+salb+'_all.sav'
mls_tp=-999.+0.*alb
mls_h2o=-999.+0.*alb
if shem eq 'south' then latlo=-1.*latlo
if shem eq 'south' then lathi=-1.*lathi
;
LIM=-999.
;SZALIM=91.	; DATA WITH SZA > SZALIM ARE BAD (IN NH THIS CAN ONLY HAPPEN ON THE ASCENDING NODE)
SZALIM=180.	; DON'T GET RID OF ANY DATA BASED ON SZA.
;
; r05 does not have doy array
;
doy=fltarr(nrev)
for i=0L,nrev-1L do begin
    iyr=long(strmid(strcompress(date(i),/remove_all),0,4))
    imn=long(strmid(strcompress(date(i),/remove_all),4,2))
    idy=long(strmid(strcompress(date(i),/remove_all),6,2))
    z=kgmt(imn,idy,iyr,iday)
    doy(i) = iday
endfor
year=long(strmid(strcompress(date,/remove_all),0,4))
;
; loop over DOYs in the cips summary file
;
n=findgen(nrev)
n1=1+findgen(nrev)
index=where(doy(n)-doy(n1) ne 0)
days=[doy(index),doy(nrev-1)]
years=[year(index),year(nrev-1)]
nday=n_elements(days)
;
; days array has to be monotonically increasing
;
index=where(days lt 100.)
if index(0) ne -1L then days(index)=days(index)+max(days)
;
; compute DFS
;
dfs=0.*doy
if shem eq 'north' then dfs=doy-172.			; June 21
if shem eq 'south' then begin
   index=where(doy gt 100.)
   dfs(index)=doy(index)-355.				; Dec 21
   index=where(doy lt 100.)
   dfs(index)=doy(index)+max(doy)-355. 
endif
;
; loop over days
;
icount=0L
;for iday=days(0),days(nday-1) do begin
for iday=days(50),days(50) do begin
;
; adjust iday back to doy
;
    kday=iday
    if year(0) mod 4 eq 0 then begin
       if iday gt 365 then kday=iday-365L
    endif
    if year(0) mod 4 ne 0 then begin
       if iday gt 366 then kday=iday-366L
    endif
;
; check for data on this day
;
    today=where(doy eq kday,nprof)
    if today(0) eq -1L then goto,jumpday
    lt_today=ltime(today,*)

    kdate,float(kday),years(icount),imn,idy
    ckday,kday,years(icount)
    z = stddat(imn,idy,years(icount),ndays)
    syrr=string(FORMAT='(I4)',years(icount))
    smn=string(FORMAT='(I2.2)',imn)
    sdy=string(FORMAT='(I2.2)',idy)
    sday=string(FORMAT='(I3.3)',kday)
    sdate=syrr+smn+sdy
    print,sdate,' ',doy(today(0)),dfs(today(0))
;
; read MLS data today and apply mask
;
    dum=findfile(mdir+'cat_mls_v4.2_'+sdate+'.sav')
    if dum(0) eq '' then goto,jumpday
    restore,mdir+'cat_mls_v4.2_'+sdate+'.sav'
    restore,mdir+'h2o_mls_v4.2_'+sdate+'.sav'
    restore,mdir+'tpd_mls_v4.2_'+sdate+'.sav'
    index=where(mask eq -99.)
    if index(0) ne -1L then mix(index)=-99.
    mh2o=mix
    good=where(mh2o ne -99.)
    mh2o(good)=mh2o(good)*1.e6
    index=where(temperature_mask eq -99.)
    if index(0) ne -1L then temperature(index)=-99.
    mtemp=temperature
    mpress=pressure
    mprof=n_elements(longitude)
    mlev=n_elements(altitude)
    muttime=time
    mlat=latitude
    mlon=longitude
    mltime=0.*time
    mldoy=0.*time
    mkltime,muttime,mlon,mltime,mldoy
;
; compute solar zenith angle for each MLS profile
;
    pi=3.14159265
    dtor=pi/180.
    earinc=23.5
    msza=-99.+0*fltarr(mprof)
    for ii=0L,mprof-1 do begin
        rlat=mlat(ii)
        rlon=mlon(ii)
        gmt=muttime(ii)
        sinlat=sin(rlat*dtor)
        coslat=sqrt(1.-sinlat^2.)
        sinlon=sin(rlon*dtor)
        coslon=cos(rlon*dtor)
        soya=(kday-81.25)*pi/182.5           ; day angle
        soha=2.*pi*(gmt-12.)/24.            ; hour angle
        soha=-soha
        sininc=sin(earinc*dtor)
        sindec=sininc*sin(soya)
        cosdec= sqrt(1.-sindec^2.)
        coszen=cos(soha)*coslon+sin(soha)*sinlon
        coszen=coszen*cosdec*coslat
        coszen=sindec*sinlat+coszen
        coszen=min([max([coszen,-1.]),1.])
        chi = acos(coszen)
        msza(ii) = chi/dtor
    endfor
;
; check
;
;if ihem eq 0 and iday eq 0 then erase
set_viewport,0.2,0.8,0.3,0.7
!type=2^2+2^3
plot,mltime,mlat,psym=2,color=0,xrange=[0,24],yrange=[-90,90],ytitle='Geographic Latitude',xtitle='Local Time',charsize=2,charthick=2,title=version

latlo2d=0.*ltime
for i=0,nrev-1L do latlo2d(i,*)=latlo
LATLO2D1=reform(LATLO2D(*,0:34))
LATLO2D2=reform(LATLO2D(*,35:69))
ltime1=reform(ltime(*,0:34))
ltime2=reform(ltime(*,35:69))

oplot,ltime2,latlo2d2,psym=1,color=250
oplot,ltime1,latlo2d1,psym=1,color=100
xyouts,18,0,'MLS',/data,color=0,charsize=2,charthick=2
xyouts,18,-15,'CIPS ASC',/data,color=100,charsize=2,charthick=2
xyouts,18,-30,'CIPS DES',/data,color=250,charsize=2,charthick=2

;if shem eq 'north' then begin
;   index=where(latitude ge min(latlo))
;   mltime=mltime(index)
;   msza=msza(index)
;   set_viewport,0.2,0.8,0.3,0.45
;   xyouts,.5,.4,'NH',/normal,charsize=2,charthick=2,color=0
;endif
;if shem eq 'south' then begin
;   index=where(latitude le max(latlo))
;   mltime=mltime(index)
;   msza=msza(index)
;   set_viewport,0.2,0.8,0.1,0.25
;   xyouts,.5,.2,'SH',/normal,charsize=2,charthick=2,color=0
;endif
;!type=2^2+2^3
;plot,mltime,msza,psym=2,color=0,xrange=[0,24],yrange=[40,120],ytitle='SZA',charsize=2,charthick=2   ;,title=sdate
;sza1=reform(sza(*,0:34))
;sza2=reform(sza(*,35:69))
;oplot,ltime2,sza2,psym=1,color=250
;oplot,ltime1,sza1,psym=1,color=100
;
    icount=icount+1L
    jumpday:
endfor	; loop over days
endfor	; loop over hemispheres
;
if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim check_mls4cips_3c_LT_SZA_'+version+'.ps -rotate -90 check_mls4cips_3c_LT_SZA_'+version+'.jpg'
;  spawn,'rm -f check_mls4cips_3c_LT_SZA_'+version+'.ps'
endif

end
