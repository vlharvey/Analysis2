;
; MLS and SABER local time evaluation

@stddat
@kgmt
@ckday
@kdate
@rd_geos5_nc3_meto
@compvort
@drawvectors

sver='v3.3'

a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill

loadct,39
mcolor=!p.color
icolmax=byte(!p.color)
mcolor=icolmax
icmm1=icolmax-1B
icmm2=icolmax-2B
device,decompose=0
nlvls=19
col1=1+indgen(nlvls)*icolmax/nlvls
!NOERAS=-1
SETPLOT='ps'
read,'setplot',setplot
nxdim=750
nydim=750
xorig=[0.15,0.15,0.55,0.55]
yorig=[0.55,0.15,0.55,0.15]
xlen=0.3
ylen=0.3
cbaryoff=0.02
cbarydel=0.01
if setplot ne 'ps' then begin
   !p.background=mcolor
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
dirm='/atmos/aura6/data/MLS_data/Datfiles_SOSST/'
sdir='/atmos/harvey/SABER_data/Datfiles/'

lstmn=1
lstdy=1
lstyr=2013
ledmn=1
leddy=31
ledyr=2013
lstday=0
ledday=0
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
;
; Compute initial Julian date
;
iyr = lstyr
idy = lstdy
imn = lstmn
z = kgmt(imn,idy,iyr,iday)
iday = iday - 1
icount=0L

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr
;
; --- Test for end condition and close windows.
;
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then stop,' normal termination condition '
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sdate=syr+smn+sdy
      print,sdate
;
; read MLS temperature and water vapor
;
    dum=findfile(dirm+'cat_mls_'+sver+'_'+sdate+'.sav')
    if dum(0) eq '' then goto,jumpday
    restore,dirm+'cat_mls_'+sver+'_'+sdate+'.sav'             ; altitude
    mprof=n_elements(longitude)
    mlat=latitude
    mlon=longitude
;
; compute local time from UT time and longitude
;
lon=longitude
MLTIME=TIME
mdoy=iday+0*time

   L=LON
   X=WHERE(L GT 180,NX)
   IF NX GT 0 THEN L(X)=360-L(X)
   X=WHERE(LON LE 180,NX)
   IF NX GT 0 THEN MLTIME(X)=TIME(X)+L(X)*24./360.
   X=WHERE(LON GT 180,NX)
   IF NX GT 0 THEN MLTIME(X)=TIME(X)-L(X)*24./360.
   X=WHERE(MLTIME GT 24,NX)
   IF NX GT 0 THEN BEGIN
      MLTIME(X)=MLTIME(X)-24
      mdoy(x)=mdoy(x)+1
   ENDIF
   X=WHERE(MLTIME LT 0,NX)
   IF NX GT 0 THEN BEGIN
      MLTIME(X)=MLTIME(X)+24
      mdoy(x)=mdoy(x)-1
   ENDIF
;
; eliminate bad uttimes
;
    index=where(time gt 0.,mprof)
    if index(0) eq -1L then goto,jump
    time=reform(time(index))
    mltime=reform(mltime(index))
    mlat=reform(mlat(index))
    mlon=reform(mlon(index))
    mdoy=reform(mdoy(index))
;
; SABER
;
    dum=findfile(sdir+'SABER_TPZ_'+sdate+'.sav')
    restore,dum
    sprof=n_elements(mode)
    slat=reform(latitude(*,100))	; extract lat, lon, time at 100 km
    slon=reform(longitude(*,100))
    time=reform(time(*,100))
;
; compute local time from UT time and longitude
;
lon=slon
SLTIME=TIME
sdoy=iday+0*time

   L=LON
   X=WHERE(L GT 180,NX)
   IF NX GT 0 THEN L(X)=360-L(X)
   X=WHERE(LON LE 180,NX)
   IF NX GT 0 THEN SLTIME(X)=TIME(X)+L(X)*24./360.
   X=WHERE(LON GT 180,NX)
   IF NX GT 0 THEN SLTIME(X)=TIME(X)-L(X)*24./360.
   X=WHERE(SLTIME GT 24,NX)
   IF NX GT 0 THEN BEGIN
      SLTIME(X)=SLTIME(X)-24
      sdoy(x)=sdoy(x)+1
   ENDIF
   X=WHERE(SLTIME LT 0,NX)
   IF NX GT 0 THEN BEGIN
      SLTIME(X)=SLTIME(X)+24
      sdoy(x)=sdoy(x)-1
   ENDIF
;
; eliminate bad uttimes
;
    index=where(time gt 0.,sprof)
    if index(0) eq -1L then goto,jump
    time=reform(time(index))
    sltime=reform(sltime(index))
    slat=reform(slat(index))
    slon=reform(slon(index))
    sdoy=reform(sdoy(index))
;
; postscript file
;
    if setplot eq 'ps' then begin
       lc=0
       xsize=nxdim/100.
       ysize=nydim/100.
       set_plot,'ps'
       !p.font=0
       device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
              /bold,/color,bits_per_pixel=8,/helvetica,filename='ylt_mls_saber_'+sdate+'.ps'
       !p.charsize=1.25
       !p.thick=2
       !p.charthick=5
       !p.charthick=5
       !y.thick=2
       !x.thick=2
    endif
;
; plot
;
    erase
    !type=2^2+2^3
    xmn=xorig(0)
    xmx=xorig(0)+xlen
    ymn=yorig(0)
    ymx=yorig(0)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    map_set,0,0,0,/contin,/grid,title=sdate,color=0
    imin=0.
    imax=24.
    for i=0L,mprof-1 do oplot,[mlon(i),mlon(i)],[mlat(i),mlat(i)],psym=8,color=((mltime(i)-imin)/(imax-imin))*mcolor
    ymnb=ymn -cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1,xtitle='MLS Local Time'
    ybox=[0,10,10,0,0]
    x2=imin
    dx=(imax-imin)/(float(nlvls)-1)
    for j=1,nlvls-1 do begin
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j)
        x2=x2+dx
    endfor
    !type=2^2+2^3
    xmn=xorig(1)
    xmx=xorig(1)+xlen
    ymn=yorig(1)
    ymx=yorig(1)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    plot,mltime,mlat,psym=1,color=0,yticks=6,ytitle='Latitude',xtitle='Local Time (hours)',yrange=[-90,90],xrange=[0,24]

    xmn=xorig(2)
    xmx=xorig(2)+xlen
    ymn=yorig(2)
    ymx=yorig(2)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    map_set,0,0,0,/contin,/grid,title=sdate,color=0,/noeras
    for i=0L,sprof-1 do oplot,[slon(i),slon(i)],[slat(i),slat(i)],psym=8,color=((sltime(i)-imin)/(imax-imin))*mcolor
    ymnb=ymn -cbaryoff
    ymxb=ymnb+cbarydel
    set_viewport,xmn+0.01,xmx-0.01,ymnb,ymxb
    !type=2^2+2^3+2^6
    plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1,xtitle='SABER Local Time'
    ybox=[0,10,10,0,0]
    x2=imin
    dx=(imax-imin)/(float(nlvls)-1)
    for j=1,nlvls-1 do begin
        xbox=[x2,x2,x2+dx,x2+dx,x2]
        polyfill,xbox,ybox,color=col1(j)
        x2=x2+dx
    endfor
    !type=2^2+2^3
    xmn=xorig(3)
    xmx=xorig(3)+xlen
    ymn=yorig(3)
    ymx=yorig(3)+ylen
    set_viewport,xmn,xmx,ymn,ymx
    plot,sltime,slat,psym=1,color=0,yticks=6,ytitle='Latitude',xtitle='Local Time (hours)',yrange=[-90,90],xrange=[0,24],/noeras

    icount=icount+1

; Close PostScript file and return control to X-windows
     if setplot ne 'ps' then stop
     if setplot eq 'ps' then begin
        device, /close
        spawn,'convert -trim ylt_mls_saber_'+sdate+'.ps -rotate -90 ylt_mls_saber_'+sdate+'.jpg'
     endif
     jumpday:
goto,jump
end
