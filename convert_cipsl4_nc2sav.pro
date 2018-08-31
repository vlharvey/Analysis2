;
; "cvo" data
; read CIPS level 4 data and store in daily IDL save files
;
restore,'read_cips_file.sav
pth='/aura7/harvey/CIPS_data/Datfiles/'

lstmn=1
lstdy=1
lstyr=2009
ledmn=1
leddy=1
ledyr=2009
lstday=0
ledday=0
;
; Ask interactive questions- get starting/ending date
;
;print, ' '
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

; USE THE CLOUD PRESENCE MAP ARRAY TO CALCULATE FREQUENCIES. AND USES THE CPM=1 VALUE TO GET THE ALBEDOS.

;LIM=1   ;LOWER LIMIT FOR ALBEDO -- ANYTHING SMALLER THAN THIS IS ASSUMED NOT TO BE A CLOUD.
;;USING LIM=-99 ESSENTIALLY INCLUDES ALL POINTS THAT ARE FOUND WITH CLOUD_PRESENCE_MAP,
;;EVEN IF THE ALBEDO IS NEGATIVE (WHICH DOES HAPPEN) -- BUT THEN THE ALB/ALB_ERR TEST
;;MIGHT CATCH IT.
;;
;ERRLIM=1.0   ;MAXIMUM ALLOWED RATIO OF ALBEDO_ERR/ALBEDO
;SZALIM=91    ;DATA WITH SZA > SZALIM ARE BAD (IN NH THIS CAN ONLY HAPPEN ON THE ASCENDING NODE)
;SZALIM=180   ;DON'T GET RID OF ANY DATA BASED ON SZA.
;
; loop over days
;
jump: iday = iday + 1
      kdate,float(iday),iyr,imn,idy
      ckday,iday,iyr

; --- Test for end condition and close windows.
      z = stddat(imn,idy,iyr,ndays)
      if ndays lt lstday then stop,' starting day outside range '
      syear=string(FORMAT='(I4)',iyr)
      smn=string(FORMAT='(I2.2)',imn)
      sdy=string(FORMAT='(I2.2)',idy)
      sday=string(FORMAT='(I3.3)',iday)
      sdate=syear+smn+sdy
;
; if save file exists then skip
;
      dum=findfile('/aura7/harvey/CIPS_data/Datfiles/*cvo*'+syear+'-'+sday+'_v03.20.sav')
      if dum(0) ne '' then goto,jump
;
; loop over orbits
;
      spawn,'ls /aura7/harvey/CIPS_data/Datfiles/*cvo*'+syear+'-'+sday+'_v03.20.nc',fnames
      norbit=n_elements(fnames)
      FOR i = 0,norbit-1 DO BEGIN
          FNAME=FNAMES(i)
          print,fname
          data=read_cips_file(fname,/full_path,attributes=attributes)
;
; extract contents from data structure
;
          32UT_TIME:
          33    description = Number of seconds that have elapsed from orbit_start_time_ut for each tile in each layer.
          34PERCENT_CLOUDS:
          35    description = Percentage of tiles that have an identified cloud
          36CLOUD_PRESENCE_MAP:
          37    description = '1' indicates a tile that has an identified cloud, '0' otherwise.
          38COMMON_VOLUME_MAP:
          39    description = '1' indicates a tile in the CVO, '0' otherwise.
          40LATITUDE:
          41
    description = Downsampled by geometry_downsample quantity. If this file is read by read_cips_1a.pro, this field will be interpolated to create a value for every pixel.
          42LONGITUDE:
          43
    description = Downsampled by geometry_downsample quantity to reduce file size. If this file is read by read_cips_1a.pro, this field will be interpolated to create a value for every pixel.
          44CLD_PHASE_ALBEDO:
          45    description = The cloud phase albedo for each tile.
          46CLD_PHASE_ALBEDO_UNC:
          47    description = The cloud phase albedo uncertainty for each tile.
          48    fillValue = 1
          49ZENITH_ANGLE_RAY_PEAK:
          50
    description = Zenith angle at the peak of the Rayleigh contribution. Downsampled by geometry_downsample quantity to reduce file size. If this file is read by read_cips_1a.pro, this field will be interpolated to create a value for every pixel.
          51    units = degrees
          52VIEW_ANGLE_RAY_PEAK:
          53
    description = Ziew angle at the peak of the Rayleigh contribution. Downsampled by geometry_downsample quantity to reduce file size. If this file is read by read_cips_1a.pro, this field will be interpolated to create a value for every pixel.
          54    units = degrees
          55SCATTERING_ANGLE:
          56
    description = Downsampled by geometry_downsample quantity to reduce file size. If this file is read by read_cips_1a.pro, this field will be interpolated to create a value for every pixel.
          57    units = degrees
          58CLD_ALBEDO:
          59    description = Cloud albedo in units of Garys (10^-6 sr^-1) i.e., albedo multiplied by 1.0e6
          60CLD_ALBEDO_UNC:
          61    description = 1 sigma formal uncertainty of cld_albedo.
          62PARTICLE_RADIUS:
          63    description = Particle radius for each tile.
          64PARTICLE_RADIUS_UNC:
          65    description = Particle radius uncertainty for each tile.
          66OZONE_COL_DENSITY:
          67    description = Ozone column density for each tile.
          68OZONE_COL_DENSITY_UNC:
          69    description = Ozone column density uncertainty for each tile.
          70SCALE_HEIGHT_RATIO:
          71    description = Scale height ratio for each tile.
          72SCALE_HEIGHT_RATIO_UNC:
          73    description = Scale height ratio uncertainty for each tile.
          74PRE_PHASE_CORR_ALB:
          75    fillValue = 1
          76PRE_PHASE_CORR_RAD:
          77    fillValue = 1

          AIM_ORBIT_NUMBER=data.AIM_ORBIT_NUMBER		;INT Cumulative mission orbit number
          VERSION=data.VERSION					;STRING    '03.20'
          PRODUCT_CREATION_TIME=data.PRODUCT_CREATION_TIME	;STRING    '2009/040-13:23:03' Version number of data product
          DEPENDANT1BVERSION=data.DEPENDANT1BVERSION		;STRING    '03.20'
          UT_DATE=DATA.UT_DATE					;LONG       2009001 UTC date of this orbit
          HEMISPHERE=DATA.HEMISPHERE				;STRING    'S'
          ORBIT_START_TIME=data.ORBIT_START_TIME		;DOUBLE 9.1480689e+14 Orbit start time in gps microseconds
          ORBIT_START_TIME_UT=data.ORBIT_START_TIME_UT		;DOUBLE   2.0548680 UTC time of the start of the orbit
          ORBIT_END_TIME=data.ORBIT_END_TIME			;DOUBLE 9.1481268e+14 Orbit end time in gps microseconds
          STACK_ID=data.STACK_ID				;INT        0 uniquely identify the Level 1B data
          XDIM=data.XDIM					;INT        X dimension of data. Average of 600
          YDIM=data.YDIM					;INT        Y dimension of data. Average is 150
          QUALITY_FLAGS=data.QUALITY_FLAGS			;LONG      TBD
          X_TILE_DIM=data.X_TILE_DIM				;INT       Array size (columns) of tiles is x direction
          Y_TILE_DIM=data.Y_TILE_DIM				;INT       Array size (rows) of tiles in the y direction
          KM_PER_PIXEL=data.KM_PER_PIXEL			;INT              5
          BBOX=data.BBOX					;LONG      Array[4] {x, y} bounding box of map projected image
          CENTER_LON=data.CENTER_LON				;DOUBLE    Center longitude of map projection, NOT data. Used for orienting the data horizontally.
          NLAYERS=(*data[0].nlayers)				;POINTER   Number of data layers where each layer corresponds to a pixel at [xDim, yDim] with a unique scattering angle (nominally, 8).

          UT_TIME=(*data[0].ut_time)				;POINTER  Number of seconds elapsed since orbit_start_time_ut
          PERCENT_CLOUDS=data.PERCENT_CLOUDS			;FLOAT    Percentage of tiles that have an identified cloud
          CLOUD_INDEX=(*DATA[0].CLOUD_PRESENCE_MAP)		;1 FOR CLOUD, 0 FOR NO CLOUD
          COMMON_VOLUME_MAP=(*DATA[0].COMMON_VOLUME_MAP)	;1=tile in the CVO, 0=otherwise
          LATITUDE=(*DATA[0].LATITUDE)				;latitude for every pixel
          LONGITUDE=(*DATA[0].LONGITUDE)			; longitude for every pixel
          X=WHERE(LATITUDE GT 90,NX)
          IF NX GT 0 THEN LATITUDE(X)=180-LATITUDE(X)
          X=WHERE(LATITUDE lt -90.)
          if nx gt 0L then latitude(x)=-90.-(latitude(x)+90.)
          X=WHERE(LONGITUDE LT 0,NX)
          IF NX GT 0 THEN LONGITUDE(X)=LONGITUDE(X)+360

          CLD_PHASE_ALBEDO=(*data[0].cld_phase_albedo)
          CLD_PHASE_ALBEDO_ERR=(*data[0].cld_phase_albedo_UNC)
          SZA = (*DATA[0].ZENITH_ANGLE_RAY_PEAK)
          VIEW = (*DATA[0].VIEW_ANGLE_RAY_PEAK)
          SZA = (*DATA[0].SCATTERING_ANGLE)
          ALB = (*data[0].cld_albedo)
          ALB_ERR = (*DATA[0].CLD_ALBEDO_UNC)
          IWC = (*data[0].ICE_WATER_CONTENT)
          IWC_ERR = (*DATA[0].ICE_WATER_CONTENT_UNC)
          RAD=(*DATA[0].PARTICLE_RADIUS)
          RAD_ERR=(*DATA[0].PARTICLE_RADIUS_UNC)
 
 HEAP_GC

 RESULT=MEMORY(/CURRENT)
 PRINT,'MEMORY IS: ',RESULT
 PRINT,' '
;
;;OMIT ANY ALBEDOS THAT ARE SMALLER THAN LIM.
;;DO NOT USE ASCENDING NODE DATA WHERE SZA GT ~91 DEG
;;OMIT ANY DATA WHERE ALB_ERR/ALB GT ERRLIM.
;
;;FIRST GET RID OF ANY INFINITE DATA AND ANY DATA ON THE ASCENDING NODE WHERE SZA GT SZA_LIM (SZA IS NEVER
;;   GT SZA_LIM ON THE DESCENDING NODE).
;;
stop
;;     BAD=WHERE(FINITE(ALB) EQ 0 OR FINITE(ALB_ERR) EQ 0 OR SZA GT SZALIM,NBAD)
;;     IF NBAD GT 0 THEN BEGIN
;;        ALB(BAD)=-99 & ALB_ERR(BAD)=-99
;;        RAD(BAD)=-99 & RAD_ERR(BAD)=-99
;;     ENDIF
;
      ENDFOR   ;*** END LOOP OVER ORBITS
;
;;  fout=pth+strdt+'.sav'
;  SAVE,LATGRID,LONGRID,ALBEDO,ALBEDO_ERR,NCLOUDS,NPOINTS,RADIUS,RADIUS_ERR,DT,file=fout

   ;IF MAX(NCLOUDS) NE 0 THEN stop

;
; save daily file
;



goto,jump
end
