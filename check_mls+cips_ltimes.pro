;
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
if setplot ne 'ps' then window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smonth=['J','F','M','A','M','J','J','A','S','O','N','D']
ndir='/Volumes/earth/harvey/NOGAPS_Alpha/Datfiles/'
mdir='/Volumes/earth/aura6/data/MLS_data/Datfiles_SOSST/'
cdir='/Volumes/Data/CIPS_data/Datfiles/cips_3c_'
;
; Ask interactive questions- get starting/ending date and p surface
;
;
; postscript file
;
    if setplot eq 'ps' then begin
       set_plot,'ps'
       xsize=nxdim/100.
       ysize=nydim/100.
       !p.font=0
       device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
              /bold,/color,bits_per_pixel=8,/helvetica,filename='check_mls+cips_ltimes.ps'
       !p.charsize=1.5
       !p.thick=2
       !p.charthick=5
       !p.charthick=5
       !y.thick=2
       !x.thick=2
    endif

icount=0
shems=['n','s']
icount=0L
for ihem=0,1 do begin
shem=shems(ihem)
for lstyr=7,7 do begin
;lstyr=7
ialb=1
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
restore,cdir+shem+'_'+syr+'_v04.20_r04_'+salb+'_all.sav'
print,'restored '+cdir+shem+'_'+syr+'_v04.20_r04_'+salb+'_all.sav'
mls_tp=-999.+0.*alb
mls_h2o=-999.+0.*alb
if shem eq 'south' then latlo=-1.*latlo
if shem eq 'south' then lathi=-1.*lathi
;
LIM=-999.
;SZALIM=91.	; DATA WITH SZA > SZALIM ARE BAD (IN NH THIS CAN ONLY HAPPEN ON THE ASCENDING NODE)
SZALIM=180.	; DON'T GET RID OF ANY DATA BASED ON SZA.
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
for iday=days(0),days(1) do begin
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
; read data on first day
;
    if icount eq 0L then begin
;
; read NOGAPS data today
;
       filename=ndir+'H2O_'+sdate+'12_360x181x60_aim9c.NAV'
       data=read_NAV_macd88(filename,lat=alat,lon=alon,p=p,yyyymmddhh=yyyymmddhh)
;
; reform into 2d arrays as MLS
;
       nr=n_elements(alat)
       nc=n_elements(alon)
       nl=n_elements(p)
       nprof=nr*nc
       nlat=fltarr(nprof)
       nlon=fltarr(nprof)
       ncount=0L
       for i=0L,nc-1L do begin
           for j=0L,nr-1L do begin
               nlon(ncount)=alon(i)
               nlat(ncount)=alat(j)
               ncount=ncount+1L
           endfor
       endfor
;
; compute local time
;
       nuttime=12.+0.*nlon
       nltime=0.*nlon
       nldoy=kday+0.*nlon
       mkltime,nuttime,nlon,nltime,nldoy
;
; compute solar zenith angle for each NOGAPS profile
;
       pi=3.14159265
       dtor=pi/180.
       earinc=23.5
       nsza=-99.+0*fltarr(nprof)
       for ii=0L,nprof-1 do begin
           rlat=nlat(ii)
           rlon=nlon(ii)
           gmt=nuttime(ii)
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
           nsza(ii) = chi/dtor
       endfor
;
; read MLS data
;
       dum=findfile(mdir+'cat_mls_v3.3_'+sdate+'.sav')
       if dum(0) eq '' then goto,jumpday
       restore,mdir+'cat_mls_v3.3_'+sdate+'.sav'
       restore,mdir+'h2o_mls_v3.3_'+sdate+'.sav'
       restore,mdir+'tpd_mls_v3.3_'+sdate+'.sav'
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
       erase
       set_viewport,0.2,0.8,0.6,0.9
       !type=2^2+2^3
       plot,mltime,mlat,psym=1,color=0,xrange=[0,24],yrange=[-90,90],/noeras,xtitle='Local Time',ytitle='Latitude',title='MLS (black) NOGAPS (grey)'
       loadct,0
       oplot,nltime,nlat,psym=1,color=150
       loadct,39
       latlo2d=0.*ltime
       for i=0,nrev-1L do latlo2d(i,*)=latlo
       LATLO2D1=reform(LATLO2D(*,0:34))
       LATLO2D2=reform(LATLO2D(*,35:69))
       ltime1=reform(ltime(*,0:34))
       ltime2=reform(ltime(*,35:69))
       oplot,ltime2,latlo2d2,psym=1,color=250
       oplot,ltime1,latlo2d1,psym=1,color=90
       plots,6,-90
       plots,6,90,thick=5,color=0,/continue
       plots,18,-90
       plots,18,90,thick=5,color=0,/continue
       xyouts,.35,.45,'ASC',color=90,/normal,charsize=2,charthick=2
       xyouts,.55,.45,'DES',color=250,/normal,charsize=2,charthick=2
       set_viewport,0.2,0.8,0.1,0.4
       plot,mltime,msza,psym=1,color=0,xrange=[0,24],yrange=[0,180],/noeras,xtitle='Local Time',ytitle='SZA'
       loadct,0
       oplot,nltime,nsza,psym=1,color=150
       loadct,39
       plots,0,90
       plots,24,90,/continue,color=0,thick=5
       day=where(msza lt 90.)
       oplot,ltime(*,0:34),sza(*,0:34),psym=1,color=90
       oplot,ltime(*,35:69),sza(*,35:69),psym=1,color=250
       oplot,mltime,msza,psym=1,color=0
       plots,6,0
       plots,6,180,/continue,color=0,thick=5
       plots,18,0
       plots,18,180,/continue,color=0,thick=5
       icount=1
    endif	; if first day

    set_viewport,0.2,0.8,0.6,0.9
    latlo2d=0.*ltime
    for i=0,nrev-1L do latlo2d(i,*)=latlo
    LATLO2D1=reform(LATLO2D(*,0:34))
    LATLO2D2=reform(LATLO2D(*,35:69))
    ltime1=reform(ltime(*,0:34))
    ltime2=reform(ltime(*,35:69))
    plot,mltime,mlat,psym=1,color=0,xrange=[0,24],yrange=[-90,90],/noeras
    oplot,ltime2,latlo2d2,psym=1,color=250
    oplot,ltime1,latlo2d1,psym=1,color=90

    set_viewport,0.2,0.8,0.1,0.4
    plot,mltime,msza,psym=1,color=0,xrange=[0,24],yrange=[0,180],/noeras
    oplot,ltime(*,0:34),sza(*,0:34),psym=1,color=90
    oplot,ltime(*,35:69),sza(*,35:69),psym=1,color=250
    oplot,mltime,msza,psym=1,color=0
;
; loop over CIPS latitude bins
; indices 0-34: MLS local time should be < 6 or > 18
; indices 35-69: MLS local time should be between 6 and 18
; 
jumpday:
endfor	; loop over days
endfor	; loop over years
endfor	; loop over hemispheres

if setplot ne 'ps' then stop
if setplot eq 'ps' then begin
   device, /close
   spawn,'convert -trim check_mls+cips_ltimes.ps -rotate -90 check_mls+cips_ltimes.jpg'
;  spawn,'rm -f check_mls+cips_ltimes.ps'
endif
end
