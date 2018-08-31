;
; compare CIPS level 4 v3.20 and v3.22 data
; VLH 4/26/2010
;
@stddat
@kgmt
@ckday
@kdate
@mkltime
@range_ring
@ks_stats

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
xorig=[0.15]
yorig=[0.15]
xlen=0.7
ylen=0.7
cbaryoff=0.02
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then begin
   lc=icolmax
   window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162
endif
mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
smonth=['J','F','M','A','M','J','J','A','S','O','N','D']
;
; restore CIPS procedures and functions
;
restore,'read_cips_file.sav
;pth='/aura7/harvey/CIPS_data/Datfiles/cips_sci_4_orbit_'
;pth='/ptmp/harvey/CIPS_data/Datfiles/cips_sci_4_orbit_'
pth='/Volumes/earth/harvey/CIPS_data/Datfiles/cips_sci_4_orbit_'

lstmn=7
lstdy=31
lstyr=2007
ledmn=7
leddy=31
ledyr=2007
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date and p surface
;
;read,' Enter starting date ',lstmn,lstdy,lstyr
;read,' Enter ending date ',ledmn,leddy,ledyr
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
;
; USE THE CLOUD PRESENCE MAP ARRAY TO CALCULATE FREQUENCIES. AND USES THE CPM=1 VALUE TO GET THE ALBEDOS.
; LOWER LIMIT FOR ALBEDO -- ANYTHING SMALLER THAN THIS IS ASSUMED NOT TO BE A CLOUD.
; USING LIM=-99 ESSENTIALLY INCLUDES ALL POINTS THAT ARE FOUND WITH CLOUD_PRESENCE_MAP,
; EVEN IF THE ALBEDO IS NEGATIVE (WHICH DOES HAPPEN) -- BUT THEN THE ALB/ALB_ERR TEST MIGHT CATCH IT.
LIM=1.
ERRLIM=1.0      ;MAXIMUM ALLOWED RATIO OF ALBEDO_ERR/ALBEDO
SZALIM=91.      ;DATA WITH SZA > SZALIM ARE BAD (IN NH THIS CAN ONLY HAPPEN ON THE ASCENDING NODE)
szamin=50.
;SZALIM=180.    ;DON'T GET RID OF ANY DATA BASED ON SZA.

; --- Loop here --------
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      if ndays gt ledday then stop,' Normal Termination Condition '
      syr=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sday=string(FORMAT='(I3.3)',iday)
      sdate=syr+smn+sdy
;
; get nc filenames on this day (I get an error if I try to read the gzipped file)
;
      spawn,'ls '+pth+'*'+syr+'-'+sday+'*.nc',fnames
      if fnames(0) eq '' then goto,jump
      norbit=n_elements(fnames)
;
; average albedo, IWC, radii in each latitude bin
;
;norbit=1L	; testing purposes
      npts=300000L
      CLOUD_INDEX_ALL=-99L+0.*lonarr(norbit,npts)
      ALB_ALL=-99.+0.*fltarr(norbit,npts)
      ALB_ERR_ALL=-99.+0.*fltarr(norbit,npts)
      SZA_ALL=-99.+0.*fltarr(norbit,npts)
      IWC_ALL=-99.+0.*fltarr(norbit,npts)
      IWC_ERR_ALL=-99.+0.*fltarr(norbit,npts)
      RAD_ALL=-99.+0.*fltarr(norbit,npts)
      RAD_ERR_ALL=-99.+0.*fltarr(norbit,npts)
      LAT_ALL=-99.+0.*fltarr(norbit,npts)
      LON_ALL=-99.+0.*fltarr(norbit,npts)
      UT_TIME_ALL=-99.+0.*fltarr(norbit,npts)
      UT_DATE_ALL=-99L+0L*intarr(norbit,npts)
;
; loop over orbits
;
      FOR iorbit = 0,norbit-1 DO BEGIN
          FNAME=FNAMES(iorbit)
          print,fname
          data=read_cips_file(fname,/full_path,attributes=attributes)
stop
;
; extract contents from data structure
;
          AIM_ORBIT_NUMBER=data.AIM_ORBIT_NUMBER                ;INT Cumulative mission orbit number
          VERSION=data.VERSION                                  ;STRING    '03.20'
          PRODUCT_CREATION_TIME=data.PRODUCT_CREATION_TIME      ;STRING    '2009/040-13:23:03' Version number of data product
          DEPENDANT1BVERSION=data.DEPENDANT1BVERSION            ;STRING    '03.20'
          UTDATE=DATA.UT_DATE                                  ;LONG       2009001 UTC date of this orbit
          HEMISPHERE=DATA.HEMISPHERE                            ;STRING    'S'
          ORBIT_START_TIME=data.ORBIT_START_TIME                ;DOUBLE 9.1480689e+14 Orbit start time in gps microseconds
          ORBIT_START_TIME_UT=data.ORBIT_START_TIME_UT          ;DOUBLE   2.0548680 UTC time of the start of the orbit
          ORBIT_END_TIME=data.ORBIT_END_TIME                    ;DOUBLE 9.1481268e+14 Orbit end time in gps microseconds
          STACK_ID=data.STACK_ID                                ;INT        0 uniquely identify the Level 1B data
          XDIM=data.XDIM                                        ;INT        X dimension of data. Average of 600
          YDIM=data.YDIM                                        ;INT        Y dimension of data. Average is 150
          QUALITY_FLAGS=data.QUALITY_FLAGS                      ;LONG      TBD
          X_TILE_DIM=data.X_TILE_DIM                            ;INT       Array size (columns) of tiles is x direction
          Y_TILE_DIM=data.Y_TILE_DIM                            ;INT       Array size (rows) of tiles in the y direction
          KM_PER_PIXEL=data.KM_PER_PIXEL                        ;INT              5
          BBOX=data.BBOX                                        ;LONG      Array[4] {x, y} bounding box of map projected image
          CENTER_LON=data.CENTER_LON    ;Center longitude of map proj, NOT data. Used for orienting the data horizontally.
          NLAYERS=(*data[0].nlayers)    ;# data layers corresponding to a pixel at [xDim, yDim] a diff SCA. Avg 8.
          UTTIME=(*data[0].ut_time)                            ;POINTER  Number of seconds elapsed since orbit_start_time_ut
          PERCENT_CLOUDS=data.PERCENT_CLOUDS                    ;FLOAT    Percentage of tiles that have an identified cloud
          CLOUD_INDEX=(*DATA[0].CLOUD_PRESENCE_MAP)             ;1 FOR CLOUD, 0 FOR NO CLOUD
;         COMMON_VOLUME_MAP=(*DATA[0].COMMON_VOLUME_MAP)        ;1=tile in the CVO, 0=otherwise
          LATITUDE=(*DATA[0].LATITUDE)                          ;latitude for every pixel
          LONGITUDE=(*DATA[0].LONGITUDE)                        ;longitude for every pixel
          X=WHERE(LATITUDE GT 90,NX)
          IF NX GT 0 THEN LATITUDE(X)=180-LATITUDE(X)           ;correct latitude for crossing over the NP
          X=WHERE(LATITUDE lt -90.,nx)
          if nx gt 0L then latitude(x)=-90.-(latitude(x)+90.)   ;correct latitude for crossing over the SP
          X=WHERE(LONGITUDE LT 0,NX)
          IF NX GT 0 THEN LONGITUDE(X)=LONGITUDE(X)+360
;         CLD_PHASE_ALBEDO=(*data[0].cld_phase_albedo)          ;The cloud phase albedo for each tile
;         CLD_PHASE_ALBEDO_ERR=(*data[0].cld_phase_albedo_UNC)  ;The cloud phase albedo uncertainty for each tile
          SZA = (*DATA[0].ZENITH_ANGLE_RAY_PEAK)                ;Zenith angle at the peak of the Rayleigh contribution
;         VIEW = (*DATA[0].VIEW_ANGLE_RAY_PEAK)                 ;View angle at the peak of the Rayleigh contribution
;         SCA = (*DATA[0].SCATTERING_ANGLE)                     ;Scattering angle at the peak of the Rayleigh contribution
          ALB = (*data[0].cld_albedo)                           ;Cloud albedo in Garys (10^-6 sr^-1) i.e., alb x 1.e6
          ALB_ERR = (*DATA[0].CLD_ALBEDO_UNC)                   ;1 sigma formal uncertainty of cld_albedo
          IWC = (*data[0].ICE_WATER_CONTENT)
          IWC_ERR = (*DATA[0].ICE_WATER_CONTENT_UNC)
          RAD=(*DATA[0].PARTICLE_RADIUS)                        ;Particle radius for each tile
          RAD_ERR=(*DATA[0].PARTICLE_RADIUS_UNC)                ;Particle radius uncertainty for each tile
;
; not in this version?
;
;         OZONE_COL_DENSITY: Ozone column density for each tile.
;         OZONE_COL_DENSITY_UNC: Ozone column density uncertainty for each tile.
;         SCALE_HEIGHT_RATIO: Scale height ratio for each tile.
;         SCALE_HEIGHT_RATIO_UNC: Scale height ratio uncertainty for each tile.
;         PRE_PHASE_CORR_ALB: fillValue = 1
;         PRE_PHASE_CORR_RAD: fillValue = 1
;
; free memory
;
          HEAP_GC
;         RESULT=MEMORY(/CURRENT)
;         PRINT,'MEMORY IS: ',RESULT
;         PRINT,' '
;
;---DO NOT--- GET RID OF ANY INFINITE DATA AND ANY negative DATA.
; Retain all data prior to search for coincidences (to compute cloud frequency)
;
;         good=WHERE(FINITE(ALB) EQ 1 and ALB gt 5. and sza gt szamin and SZA le SZALIM,ngood)
          good=WHERE(finite(uttime) eq 1 and sza gt szamin and SZA le SZALIM,ngood)	; all data
          IF NGOOD GT 0 THEN BEGIN
             if ngood gt npts then stop,'increase npts'
;    
; save petals into daisy
;
             CLOUD_INDEX_ALL(iorbit,0:ngood-1L)=cloud_index(good)
             ALB_ALL(iorbit,0:ngood-1L)=alb(good)
             ALB_ERR_ALL(iorbit,0:ngood-1L)=alb_err(good)
             SZA_ALL(iorbit,0:ngood-1L)=sza(good)
             IWC_ALL(iorbit,0:ngood-1L)=iwc(good)
             IWC_ERR_ALL(iorbit,0:ngood-1L)=iwc_err(good)
             RAD_ALL(iorbit,0:ngood-1L)=rad(good)
             RAD_ERR_ALL(iorbit,0:ngood-1L)=rad_err(good)
             LAT_ALL(iorbit,0:ngood-1L)=latitude(good)
             LON_ALL(iorbit,0:ngood-1L)=longitude(good)
             UT_TIME_ALL(iorbit,0:ngood-1L)=uttime(good)
             UT_DATE_ALL(iorbit,0:ngood-1L)=utdate(good)
;print,'UT time ',min(uttime(good)),max(uttime(good)),ngood,' points'
print,'min/max lat ',min(latitude(good)),max(latitude(good)),ngood,' points'
          ENDIF
      endfor  ; loop over orbits
;
; make arrays 1-d
;
      good=where(lat_all ge -90.,ngood)
      cloud_index=cloud_index_all(good)
      alb=alb_all(good)
      ALB_ERR=alb_err_all(good)
      SZA=sza_all(good)
      IWC=iwc_all(good)
      IWC_ERR=iwc_err_all(good)
      RAD=rad_all(good)
      RAD_ERR=rad_err_all(good)
      clat=lat_all(good)
      clon=lon_all(good)
      cuttime=ut_time_all(good)
      cutdate=ut_date_all(good)
;
; clear memory
;
      ALB_ALL=0 & ALB_ERR_ALL=0 & SZA_ALL=0 & VIEW_ALL=0 & SCA_ALL=0 & IWC_ALL=0 & IWC_ERR_ALL=0
      RAD_ALL=0 & RAD_ERR_ALL=0 & LAT_ALL=0 & LON_ALL=0 & UT_TIME_ALL=0 & UT_DATE_ALL=0 & cloud_index_all=0
;
; Compute LOCAL TIME and LOCAL DAY from UT time and longitude
;
      lday=iday+0L*indgen(ngood)
      MKLTIME,CUTTIME,CLON,LTIME,lday
      clttime=ltime
      clday=lday
      print,'CLTTIME ',min(clttime),max(clttime)
      print,'CLDAY ',min(clday),max(clday)
;;loadct,39
;;erase
;;plot,cuttime,clon,psym=1,color=0,xrange=[0.,24.]
;;oplot,clttime,clon,psym=1,color=.3*mcolor
;
; plot CIPS locations colored by local time
;
;
; postscript file
; 
      if setplot eq 'ps' then begin
         lc=0
         set_plot,'ps'
         xsize=nxdim/100.
         ysize=nydim/100.
         !p.font=0
         device,/landscape,/inch,xoff=4.25-ysize/2.,yoff=5.5+xsize/2.,xsize=xsize,ysize=ysize,$
                /bold,/color,bits_per_pixel=8,/helvetica,filename='compare_cips_l4_v3.20_v3.22_rawfiles_'+sdate+'.ps'
         !p.charsize=1.25
         !p.thick=2
         !p.charthick=5
         !p.charthick=5
         !y.thick=2
         !x.thick=2
      endif

      erase
      loadct,39
      xmn=xorig(0)
      xmx=xorig(0)+xlen
      ymn=yorig(0)
      ymx=yorig(0)+ylen
      set_viewport,xmn,xmx,ymn,ymx
      !type=2^2+2^3
;     if min(clat) gt 0. then map_set,90,0,0,/ortho,/contin,/grid,title='CIPS '+sdate,color=0,limit=[50.,0.,80.,360.],/noerase
;     if min(clat) lt 0. then map_set,-90,0,0,/ortho,/contin,/grid,title='CIPS '+sdate,color=0,limit=[-50.,0.,-80.,360.],/noerase
;     for ii=0L,ngood-1L do $
;         oplot,[clon(ii),clon(ii)],[clat(ii),clat(ii)],psym=8,color=((clttime(ii)-utmin)/(utmax-utmin))*mcolor,symsize=0.1
;     ymxb=ymn-cbaryoff
;     ymnb=ymxb+cbarydel
;     set_viewport,xmn,xmx,ymnb,ymxb
;     imin=utmin
;     imax=utmax
;     !type=2^2+2^3+2^6
;     plot,[imin,imax],[0,0],yrange=[0,10],xrange=[imin,imax],/noeras,color=0,charsize=1,xtitle='Local Time (hours)'
;     ybox=[0,10,10,0,0]
;     x2=imin
;     nlvls=11
;     col1=1+indgen(nlvls)*mcolor/nlvls
;     dx=(imax-imin)/(float(nlvls)-1)
;     for j=1,nlvls-1 do begin
;         xbox=[x2,x2,x2+dx,x2+dx,x2]
;         polyfill,xbox,ybox,color=col1(j)
;         x2=x2+dx
;     endfor
;
    if setplot ne 'ps' then stop
    if setplot eq 'ps' then begin
       device, /close
       spawn,'convert -trim compare_cips_l4_v3.20_v3.22_rawfiles_'+sdate+'.ps -rotate -90 compare_cips_l4_v3.20_v3.22_rawfiles_'+sdate+'.jpg'
       spawn,'rm -f compare_cips_l4_v3.20_v3.22_rawfiles_'+sdate+'.ps'
    endif
      skipmls:
      icount=icount+1L
goto,jump
end
