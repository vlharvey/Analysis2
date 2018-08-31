@get_day_string
;========================================================================
;  pro SABER_CALC_BT03
;========================================================================
;
; Caution: Do not surround a yaw period
;
;  2003:  Calendar days: Jan 15, Mar 19, May 22, July 17, Sept 19, Nov 22
;         Julian days:   15,         78,     142,     198,    262,   326
;
COMMON   ZPARM, ZGRID, NUMZ, ZINC
COMMON LATPARM, LATITUDE,  NUMLAT, LATINC
COMMON LONPARM, LONGITUDE, NUMLON, LONINC
COMMON CORVALS,CORPARM, TANLATA, LATRAD, ACOSPHI

;close,/all

;var='X' ; change this variable to 'X' to save fields during the execution of
        ; the program or to 'saber_fluxes' if using fields that have already
        ; been calculated

PLOT_SUMMARY = 0

max_orbs_per_day     = 16
max_number_of_events = 100L
maxlevs              = 400
sec_per_hour         = 60.*60.

IF N_ELEMENTS(YEAR) NE 1 THEN READ,'ENTER YEAR (yyyy): ', YEAR

lineyear:
;
; http://saber.gats-inc.com/operations/Yearly/index.php for yaw dates
;
IF YEAR EQ 2002 THEN STDAY  = 25
IF YEAR EQ 2002 THEN ENDDAY = 78
  
IF YEAR EQ 2003 THEN STDAY  = 16
IF YEAR EQ 2003 THEN ENDDAY = 77

IF YEAR EQ 2004 THEN STDAY  = 16
IF YEAR EQ 2004 THEN ENDDAY = 76
 
IF YEAR EQ 2005 THEN STDAY  = 15
IF YEAR EQ 2005 THEN ENDDAY = 76
 
IF YEAR EQ 2006 THEN STDAY  = 14
IF YEAR EQ 2006 THEN ENDDAY = 76
 
IF YEAR EQ 2007 THEN STDAY  = 13
IF YEAR EQ 2007 THEN ENDDAY = 74
 
;IF YEAR EQ 2008 THEN STDAY  = 16	; 9-17 missing
IF YEAR EQ 2008 THEN STDAY  = 18
IF YEAR EQ 2008 THEN ENDDAY = 74
 
IF YEAR EQ 2009 THEN STDAY  = 12
IF YEAR EQ 2009 THEN ENDDAY = 73
 
IF YEAR EQ 2010 THEN STDAY  = 12
IF YEAR EQ 2010 THEN ENDDAY = 73
 
IF YEAR EQ 2011 THEN STDAY  = 11
IF YEAR EQ 2011 THEN ENDDAY = 72

IF YEAR EQ 2012 THEN STDAY  = 10
IF YEAR EQ 2012 THEN ENDDAY = 71

IF YEAR EQ 2013 THEN STDAY  = 8
IF YEAR EQ 2013 THEN ENDDAY = 69

IF YEAR EQ 2014 THEN STDAY  = 6
IF YEAR EQ 2014 THEN ENDDAY = 69

SOUTH    = -85
NORTH    =  85
LATINC   =   5.
NUMLAT   = FIX( (NORTH-SOUTH)/LATINC) + 1
LATITUDE = SOUTH + FINDGEN(NUMLAT)*LATINC
LONINC   =  30.
NUMLON    = FIX( 360./LONINC)
LONGITUDE = FINDGEN(NUMLON)*LONINC
ERAD        = 6.37E6
SEC_PER_DAY = 60.*60.* 24.
PI2         = 2.*!pi
OMEGA       = PI2/SEC_PER_DAY
LATRAD      = LATITUDE*!dtor
CORPARM     = 2.*OMEGA*SIN(LATRAD)
TANLATA     = TAN(LATRAD)/ERAD
ACOSPHI     = ERAD * COS(LATRAD)
MAXWAVES   = NUMLON/2
NUMZ      = 120   ; Maximum number of levels on a fixed altitude grid
ZINC      = .15
ZBASE     = 2.
ZGRID     = ZBASE + FINDGEN(NUMZ)*ZINC
numdays  = endday - stday + 1
MAXDAYS = NUMDAYS
;
; Build file name
;
lteheader     = '_nlte'
wildcard      = '*'
yearstring    = STRTRIM( STRING(year,'(I4)') )
;dirheader     = '/Volumes/MacD68-2/france/data/SABER_data/'
dirheader='/atmos/harvey/SABER_data/Datfiles/'
print, dirheader
;
; Allow for 16 orbits in a day
; Originally set for 15 but has problems in 2002
; 2013 has 17 orbits?
;
MAX_IN_BIN    = 16
MAX_IN_BIN    = 17
KPLOT = numz/2

AV_AGEO   = FLTARR(numdays,NUMLON,NUMLAT,NUMZ)
ZM_AGEO   = FLTARR(numdays,NUMLAT,NUMZ)
AV_ATEMP  = FLTARR(numdays,NUMLON,NUMLAT,NUMZ)
ZM_ATEMP  = FLTARR(numdays,NUMLAT,NUMZ)
 
UWIND          = FLTARR(NUMLAT,numdays, NUMZ)
UPERTWIND      = COMPLEXARR(NUMLAT,numdays,MAXWAVES,NUMZ)
VPERTWIND      = COMPLEXARR(NUMLAT,numdays,MAXWAVES,NUMZ)

ZONAL_TSPECTRA  = COMPLEXARR(NUMLAT,numdays,MAXWAVES,NUMZ)
ZONAL_GSPECTRA  = COMPLEXARR(NUMLAT,numdays,MAXWAVES,NUMZ)
STAWAVES_G       = FLTARR(NUMLAT,NUMZ, NUMLON/2+1, 2,maxdays)
STAWAVES_T       = FLTARR(NUMLAT,NUMZ, NUMLON/2+1, 2,maxdays)
 
UPROFILES      = FLTARR(numdays,NUMLON,NUMLAT,NUMZ)
VPROFILES      = FLTARR(numdays,NUMLON,NUMLAT,NUMZ)
GPROFILES      = FLTARR(numdays,NUMLON,NUMLAT,NUMZ)
TPROFILES      = FLTARR(numdays,NUMLON,NUMLAT,NUMZ)

EPDIV_ALL = FLTARR(NUMLAT,numdays, NUMZ)
EPY_ALL = FLTARR(NUMLAT,numdays, NUMZ)
EPZ_ALL = FLTARR(NUMLAT,numdays, NUMZ)

NUM_IN_BINS     = INTARR(NUMLON,NUMLAT)
AGEO_REC        = FLTARR(NUMZ, MAX_IN_BIN,NUMLON,NUMLAT)
ATEMP_REC        = FLTARR(NUMZ, MAX_IN_BIN,NUMLON,NUMLAT)

first_orbit_set_today = 0
last_orb_today  = 0

daystring = GET_DAY_STRING(stday)
yrdayheader = yearstring + daystring
if stday le 31L then saber_header = 'SABER_Temp_O3_January'+yearstring
if stday gt 31L then saber_header = 'SABER_Temp_O3_February'+yearstring
prefix=dirheader+saber_header+'*'
length_of_prefix = strlen(prefix)
;
; Find files for today
;
this_days_files = prefix + wildcard
print, this_days_files
;
; Do we have a SABER file for this day?
;
todays_files = FINDFILE(this_days_files, COUNT=cvar )

filename           = todays_files(0)
length_of_filename = strlen(filename)
;
; Get the data (netcdf reader)
;
print, ' Reading ', filename
read_netCDF, data, filename, attributes, status
;
names_of_fields = tag_names(data)
FOR I = 0, N_TAGS(data) - 1 do print, i, names_of_fields(i)

index = where(names_of_fields eq 'ORBIT')
orbit1       = data.(index(0))

index = where(names_of_fields eq 'TPAD')
;adflag1       = data.(index(0))

index = where(names_of_fields eq 'DATE')
date1       = data.(index(0))

index = where(names_of_fields eq 'TIME')
utimes1       = data.(index(0))*1.e-3/sec_per_hour

index = where(names_of_fields eq 'TPALTITUDE')
tpaltitudes1  = data.(index(0))

index = where(names_of_fields eq 'TPLATITUDE')
tplatitudes1  = data.(index(0))

index = where(names_of_fields eq 'TPLONGITUDE')
tplongitudes1 = data.(index(0))

index = where(names_of_fields eq 'TPSOLARLT')
tploct1       = data.(index(0))*1.e-3/sec_per_hour

index = where(names_of_fields eq 'PRESSURE')
pressures1    = data.(index(0))

index = where(names_of_fields eq 'KTEMP')         ;Temperature
temps1        = data.(index(0))

index = where(names_of_fields eq 'TPGPALTITUDE')   ;Geopotential
geops1 = data.(index(0))

dates = uniq(date1)
ndays = n_elements(dates)
for doy = stday, endday do begin
    dind = doy - stday
;
;For monthly files, restore February on day 32
;
    if doy eq 32L then begin
       saber_header  = 'SABER_Temp_O3_February'+yearstring
       prefix=dirheader+saber_header+'*'
       length_of_prefix = strlen(prefix)
;
; Find files for today
;
       this_days_files = prefix + wildcard
       print, this_days_files
;
; Do we have a SABER file for this day?
;
       todays_files = FINDFILE(this_days_files, COUNT=cvar )
       filename           = todays_files(0)
       length_of_filename = strlen(filename)
;
; Get the data (netcdf reader)
;
       print, ' Reading ', filename
       read_netCDF, data, filename, attributes, status
       names_of_fields = tag_names(data)
       FOR I = 0, N_TAGS(data) - 1 do print, i, names_of_fields(i)

       index = where(names_of_fields eq 'ORBIT')
       orbit1       = data.(index(0))

       index = where(names_of_fields eq 'TPAD')
;       adflag1       = data.(index(0))

       index = where(names_of_fields eq 'DATE')
       date1       = data.(index(0))

       index = where(names_of_fields eq 'TIME')
       utimes1       = data.(index(0))*1.e-3/sec_per_hour

       index = where(names_of_fields eq 'TPALTITUDE')
       tpaltitudes1  = data.(index(0))

       index = where(names_of_fields eq 'TPLATITUDE')
       tplatitudes1  = data.(index(0))

       index = where(names_of_fields eq 'TPLONGITUDE')
       tplongitudes1 = data.(index(0))

       index = where(names_of_fields eq 'TPSOLARLT')
       tploct1       = data.(index(0))*1.e-3/sec_per_hour

       index = where(names_of_fields eq 'PRESSURE')
       pressures1    = data.(index(0))

       index = where(names_of_fields eq 'KTEMP')         ;Temperature
       temps1        = data.(index(0))

       index = where(names_of_fields eq 'TPGPALTITUDE')   ;Geopotential
       geops1 = data.(index(0))

    endif

;For monthly files, restore March

    if (doy eq 60L and year mod 4 gt 0) or (doy eq 61L and year mod 4 eq 0) then begin

       saber_header  = 'SABER_Temp_O3_March' +yearstring
       prefix=dirheader+saber_header+'*'
       length_of_prefix = strlen(prefix)
;
; Find files for today
;
       this_days_files = prefix + wildcard
       print, this_days_files
;
; Do we have a SABER file for this day?
;
       todays_files = FINDFILE(this_days_files, COUNT=cvar )

       filename           = todays_files(0)
       length_of_filename = strlen(filename)
;
; Get the data (netcdf reader)
;
       print, ' Reading ', filename
       read_netCDF, data, filename, attributes, status
       names_of_fields = tag_names(data)
       FOR I = 0, N_TAGS(data) - 1 do print, i, names_of_fields(i)

       index = where(names_of_fields eq 'ORBIT')
       orbit1       = data.(index(0))

       index = where(names_of_fields eq 'TPAD')
;       adflag1       = data.(index(0))

       index = where(names_of_fields eq 'DATE')
       date1       = data.(index(0))

       index = where(names_of_fields eq 'TIME')
       utimes1       = data.(index(0))*1.e-3/sec_per_hour

       index = where(names_of_fields eq 'TPALTITUDE')
       tpaltitudes1  = data.(index(0))

       index = where(names_of_fields eq 'TPLATITUDE')
       tplatitudes1  = data.(index(0))

       index = where(names_of_fields eq 'TPLONGITUDE')
       tplongitudes1 = data.(index(0))

       index = where(names_of_fields eq 'TPSOLARLT')
       tploct1       = data.(index(0))*1.e-3/sec_per_hour

       index = where(names_of_fields eq 'PRESSURE')
       pressures1    = data.(index(0))

       index = where(names_of_fields eq 'KTEMP')         ;Temperature
       temps1        = data.(index(0))

       index = where(names_of_fields eq 'TPGPALTITUDE')   ;Geopotential
       geops1 = data.(index(0))
    endif
;
;Select profiles (events) for each date and each orbit strip
;
    dateevents = where(strmid(strtrim(strcompress(date1),2),4,3) eq doy)

    if doy eq stday then yyyyddd = strarr(endday - stday + 1L)
    yyyyddd[doy - stday] = yearstring+string(doy, '(I3.3)')

    dailyorbits = reform(orbit1[dateevents])
    orbits2 = uniq(dailyorbits)
    norbits = n_elements(orbits2)

    print, 'DOY: ',doy,',   Daily orbits: ',norbits

    uts_this_day     = fltarr(max_number_of_events)
    loct_this_day    = fltarr(max_number_of_events)
    tplats_this_day  = fltarr(max_number_of_events)
    tplons_this_day  = fltarr(max_number_of_events)
    num_in_bins      = num_in_bins*0L

    for iorbit = 0L, norbits - 1L do begin
        events = where(orbit1 eq dailyorbits[orbits2[iorbit]] and strmid(strtrim(strcompress(date1),2),4,3) eq doy,number_of_events)
        orbit        = orbit1[events]
;        adflag       = adflag1[events]
        date         = date1[events]
        utimes       = utimes1[*,events]
        tpaltitudes  = tpaltitudes1[*,events]
        tplatitudes  = tplatitudes1[*,events]
        tplongitudes = tplongitudes1[*,events]
        tploct       = tploct1[*,events]
        pressures    = pressures1[*,events]
        temps        = temps1[*,events]
        geops        = geops1[*,events]
        print, 'number of events: ', number_of_events 
 
        for nev = 0, number_of_events-1, 1 do begin
            nalts = n_elements(tpaltitudes(*,nev) )
;
;  Coordinates set at nalts/2 level value.
;
            uts_this_day(nev)  = utimes(nalts/2,nev)
            loct_this_day(nev)  = tploct(nalts/2,nev)
            tplats_this_day(nev) = tplatitudes(nalts/2, nev)
            tplons_this_day(nev)  = tplongitudes(nalts/2, nev)

            GRID_IN_LAT,nev,loct_this_day,tplats_this_day,tplons_this_day,uts_this_day,$
                        pressures,geops,temps,NUM_IN_BINS,AGEO_REC,ATEMP_REC,PLOT_SUMMARY
        endfor ; nev loop
    endfor ; cvar loop

    AV_BINS, NUM_IN_BINS, AGEO_REC, AV_AGEO, ZM_AGEO, zonal_gspectra, STAWAVES_G, DIND
    AV_BINS, NUM_IN_BINS, ATEMP_REC, AV_ATEMP, ZM_ATEMP, zonal_tspectra, STAWAVES_T, DIND

    GET_MEAN_WIND,  DIND, ZM_AGEO, UWIND
    GET_PERTU_WIND, DIND, MAXWAVES, zonal_gspectra, UWIND, UPERTWIND
    GET_PERTV_WIND, DIND, MAXWAVES, zonal_gspectra, UWIND, VPERTWIND

    BUILD_WIND_PROFILES, UPERTWIND,             DIND, UPROFILES
    BUILD_WIND_PROFILES, VPERTWIND,             DIND, VPROFILES
    BUILD_WIND_PROFILES, zonal_gspectra*1.E3,   DIND, GPROFILES
    BUILD_WIND_PROFILES, zonal_tspectra,        DIND, TPROFILES
;
;Convert temperature to Potential Temperature
;
    GET_THETA, TPROFILES, ZM_ATEMP, THETA, THETA_Z, DIND
    GET_EP_FLUXES, ZM_ATEMP, THETA_Z, UPROFILES, VPROFILES, THETA, EPDIV, EPY, EPZ, DIND

    EPDIV_ALL(*,DIND,*) = EPDIV(*,*)
    EPY_ALL(*,DIND,*) = EPY(*,*)
    EPZ_ALL(*,DIND,*) = EPZ(*,*)

nextday:
endfor ; doy loop
;
;.SAV FILES USED IN SABER_PLOTS_BT01.PRO
;YYYY = YEAR
;0XX = START DAY (WITH PRECEEDING ZERO)
;0XX = END DAY (WITH PRECEEDING ZERO)
;
SAVE, FILENAME=DIRHEADER+'/SABER_TEMP_'+yearstring+'_'+strtrim(strcompress( STRING(STDAY,'(I3)')),2)+$
'_'+strtrim(strcompress( STRING(ENDDAY,'(I3)')),2)+'.SAV',AV_ATEMP, ZM_ATEMP, ZONAL_TSPECTRA,latitude, longitude,zgrid,yyyyddd
SAVE, FILENAME=DIRHEADER+'/SABER_GEOP_'+yearstring+'_'+strtrim(strcompress( STRING(STDAY,'(I3)')),2)+$
'_'+strtrim(strcompress( STRING(ENDDAY,'(I3)')),2)+'.SAV', AV_AGEO, ZM_AGEO, ZONAL_GSPECTRA,latitude, longitude,zgrid,yyyyddd
SAVE, FILENAME=DIRHEADER+'/SABER_WINDS_'+yearstring+'_'+strtrim(strcompress( STRING(STDAY,'(I3)')),2)+$
'_'+strtrim(strcompress( STRING(ENDDAY,'(I3)')),2)+'.SAV', UWIND, UPERTWIND, VPERTWIND,latitude, longitude,zgrid,yyyyddd
SAVE, FILENAME=DIRHEADER+'/SABER_PROFS_'+yearstring+'_'+strtrim(strcompress( STRING(STDAY,'(I3)')),2)+$
'_'+strtrim(strcompress( STRING(ENDDAY,'(I3)')),2)+'.SAV', UPROFILES, VPROFILES, GPROFILES, TPROFILES,latitude, longitude,zgrid,yyyyddd
SAVE, FILENAME=DIRHEADER+'/SABER_FLUXES_'+yearstring+'_'+strtrim(strcompress( STRING(STDAY,'(I3)')),2)+$
'_'+strtrim(strcompress( STRING(ENDDAY,'(I3)')),2)+'.SAV', EPDIV_ALL, EPY_ALL, EPZ_ALL,latitude, longitude,zgrid,yyyyddd
SAVE, FILENAME=DIRHEADER+'/SABER_STAWAVES_'+yearstring+'_'+strtrim(strcompress( STRING(STDAY,'(I3)')),2)+$
'_'+strtrim(strcompress( STRING(ENDDAY,'(I3)')),2)+'.SAV', STAWAVES_G, STAWAVES_T,latitude, longitude,zgrid,yyyyddd

;return
end
;===================================================================
PRO GET_THETA, TPROFILES, ZM_TEMP, THETA, THETA_Z, DIND
;===================================================================

COMMON   ZPARM, ZGRID, NUMZ, ZINC
COMMON LATPARM, LATITUDE,  NUMLAT, LATINC
COMMON LONPARM, LONGITUDE, NUMLON, LONINC

RAIR =  287.
CP = 1004.
KAPPA = RAIR/CP
HSC = 6.5E3
KGRID       = EXP(KAPPA*ZGRID)
;THETA       = TPROFILES
;THETA_Z     = ZM_TEMP
THETA       = FLTARR(NUMZ,NUMLON,NUMLAT)
THETA_Z     = FLTARR(NUMLAT, NUMZ)

FOR J = 0, NUMLAT-1 DO BEGIN
    THETA_MEAN = FLTARR(NUMZ)
    FOR K = 0, NUMZ-1 DO BEGIN
        THETA_MEAN(K) = ZM_TEMP(DIND,J,K)*KGRID(K)
        FOR I = 0, NUMLON-1 DO BEGIN
            THETA(K,I,J) = TPROFILES(DIND,I,J,K)*KGRID(K)
        ENDFOR  ;I LOOP
    ENDFOR ; K LOOP
    K = 0
    THETA1  = THETA_MEAN(K+1)
    THETA0  = THETA_MEAN(K)
    IF THETA1 NE 0. AND THETA0 NE 0. THEN THETA_Z(J,K) = (THETA1 - THETA0)/ZINC/HSC
    FOR K = 1, NUMZ-2 DO BEGIN
        THETA1   = THETA_MEAN(K+1)
        THETA_1  = THETA_MEAN(K-1)
        IF THETA1 NE 0. AND THETA_1 NE 0. THEN THETA_Z(J,K) = (THETA1 - THETA_1)/ZINC/HSC/2.
    ENDFOR ; K LOOP
    K = NUMZ-1
    THETA0   = THETA_MEAN(K)
    THETA_1  = THETA_MEAN(K-1)
    IF THETA0 NE 0. AND THETA_1 NE 0. THEN THETA_Z(J,K) = (THETA0 - THETA_1)/ZINC/HSC
ENDFOR ; J LOOP
RETURN
END
;=====================================================================
 PRO GET_EP_FLUXES, ZM_TEMPS, THETA_Z, UPROFILES, VPROFILES, $
                    THETA, EPDIV, EPY_TEMP, EPZ_TEMP, DIND
;=====================================================================

 COMMON   ZPARM, ZGRID, NUMZ, ZINC
 COMMON LATPARM, LATITUDE,  NUMLAT, LATINC
 COMMON LONPARM, LONGITUDE, NUMLON, LONINC

 sec_per_Day = 60. *60.*24.
 omega = !pi * 2./sec_per_Day
 erad = 6.37e6                         ;Earth's Radius
 dlat = latinc *!dtor                  ;Degrees to Radians
  dy  = dlat*erad                      ;Radians to Meters 

 RAIR =  287.
   CP = 1004.
 KAPPA = RAIR/CP
    P0 = 1.E5
    F0 = 2. * OMEGA * SIN(latitude*!dtor)
    HSC = 6.5E3

  PRESS = P0*EXP(-ZGRID)
 
     RHO  = FLTARR(NUMLAT, NUMZ)
    EPDIV = FLTARR(NUMLAT, NUMZ)
      EPY = FLTARR(NUMLAT, NUMZ)
      EPZ = FLTARR(NUMLAT, NUMZ)
  MOMFLUX = FLTARR(NUMLAT, NUMZ)
 HEATFLUX = FLTARR(NUMLAT, NUMZ)

 EPZ_TEMP = FLTARR(NUMLAT, NUMZ)
 EPY_TEMP = FLTARR(NUMLAT, NUMZ)

 FOR J = 0, NUMLAT-1 DO BEGIN 

  FOR K = 0, NUMZ-1 DO BEGIN 
 
    TAIR = ZM_TEMPS(DIND,J,K)
    MOM_AT_LONGITUDE = FLTARR(NUMLON)
    HEAT_AT_LONGITUDE = FLTARR(NUMLON)

    IF TAIR NE 0. THEN BEGIN 

     RHO(J,K)  = PRESS(K)/RAIR/TAIR
   
     FOR I = 0, NUMLON-1 DO BEGIN

      MOM_AT_LONGITUDE(I)  = -VPROFILES(DIND,I,J,K)*UPROFILES(DIND,I,J,K)
      HEAT_AT_LONGITUDE(I)  = VPROFILES(DIND,I,J,K)*THETA(K,I,J)

     ENDFOR ; I LOOP 
;
;    Zonal average
;
     FULL = WHERE(MOM_AT_LONGITUDE NE 0. )

     IF FULL(0) NE -1 THEN  BEGIN 

       MOMFLUX(J,K)  = RHO(J,K)*MEAN(MOM_AT_LONGITUDE(FULL) )
       HEATFLUX(J,K) = RHO(J,K)*F0(J) * $
             MEAN(HEAT_AT_LONGITUDE(FULL) )/THETA_Z(J,K)
     
     ENDIF ; FULL(0) NE -1 

    ENDIF ; ZM_TEMPS(K,J) NE 0. 

  ENDFOR ; K LOOP 
;
; Vertical flux divergence
;
  K = 0
  THETA1  = HEATFLUX(J,K+1)
  THETA0  = HEATFLUX(J,K)

  IF THETA1 NE 0. AND THETA0 NE 0. THEN $
      EPZ(J,K) = (THETA1 - THETA0)/ZINC/HSC

  FOR K = 1, NUMZ-2 DO BEGIN 

   THETA1   = HEATFLUX(J,K+1)
   THETA_1  = HEATFLUX(J,K-1)

   IF THETA1 NE 0. AND THETA_1 NE 0. THEN $
      EPZ(J,K) = (THETA1 - THETA_1)/ZINC/HSC/2.

  ENDFOR ; K LOOP 

  K        = NUMZ-1
  THETA0    = HEATFLUX(J,K)
  THETA_1  = HEATFLUX(J, K-1)

  IF THETA0 NE 0. AND THETA_1 NE 0. THEN $
      EPZ(J,K) = (THETA0 - THETA_1)/ZINC/HSC
 ;
 ; Horizontal flux divergence
 ;

 K = 0
  THETA1  = MOMFLUX(J,K+1)
  THETA0  = MOMFLUX(J,K)

  IF THETA1 NE 0. AND THETA0 NE 0. THEN $
      EPY(J,K) = (THETA1 - THETA0)/(dy)    

  FOR K = 1, NUMZ-2 DO BEGIN 

   THETA1   = MOMFLUX(J, K+1)
   THETA_1  = MOMFLUX(J, K-1)

   IF THETA1 NE 0. AND THETA_1 NE 0. THEN $
      EPY(J,K) = (THETA1 - THETA_1)/(dy)/2.

  ENDFOR ; K LOOP 

  K        = NUMZ-1
  THETA0    = MOMFLUX(J,K)
  THETA_1  = MOMFLUX(J,K-1)

  IF THETA0 NE 0. AND THETA_1 NE 0. THEN $
      EPY(J,K) = (THETA0 - THETA_1)/(dy)

 
  FOR K = 0, NUMZ-1 DO BEGIN
  IF RHO(J,K) NE 0. THEN $
       
       EPZ_TEMP(J,K) = EPZ(J,K)*sec_per_Day/RHO(J,K)
       EPY_TEMP(J,K) = EPY(J,K)*sec_per_Day/RHO(J,K)
       EPDIV(J,K) =  (EPZ(J,K)  + EPY(J,K))*sec_per_Day/RHO(J,K)

 ENDFOR ;K LOOP
 ENDFOR ; J LOOP 
 
 RETURN
 END


;==================================================================
 PRO  AV_BINS, NUM_IN_BINS, AGEO_REC, AV_AGEO_REC, ZM_AGEO_REC, $
                          zonal_spectra, STAWAVES, DIND
;===========================================================

COMMON   ZPARM, ZGRID, NUMZ, ZINC
COMMON LATPARM, LATITUDE,  NUMLAT, LATINC
COMMON LONPARM, LONGITUDE, NUMLON, LONINC

IF NUMLON MOD 2 EQ 0 THEN KNYQ = NUMLON/2-1 ELSE $
                          KNYQ = NUMLON/2
                           
 MINLON       = NUMLON - 2

 MAXWAVES   = NUMLON/2 

 FOR J = 0, NUMLAT - 1 DO BEGIN 

     FOR I = 0, NUMLON - 1 DO BEGIN
       
       N1 = NUM_IN_BINS(I,J)
                
       FOR A = 0,NUMZ-1 DO BEGIN                                    
      
        WORKTEMPS =  AGEO_REC(A,0:N1,I,J) 
        FULL = WHERE(WORKTEMPS NE 0.) 
        IF FULL(0) NE -1 THEN $
         AV_AGEO_REC(dind,I,J,A) = MEAN(WORKTEMPS(FULL) )  ; avg in lat-lon
         
      ENDFOR   ; A LOOP

     ENDFOR ; I loop 

 ENDFOR ; J loop 

   FOR J = 0, NUMLAT - 1 DO BEGIN 

     FOR A = 0, NUMZ -1 DO BEGIN

         WORKARR       = FLTARR(NUMLON)
         WORKARR(*)    = av_ageo_rec(DIND,0:NUMLON-1,J,A)
        FULL_THIS_A = WHERE( WORKARR NE 0.)

        IF FULL_THIS_A(0) NE -1 THEN BEGIN
;
;        Zonal mean geopotential
;
         ZM_AGEO_REC(dind,J,A) = MEAN(WORKARR(FULL_THIS_A) )
;        
;        FFT of geopotential
;
         INTERP, FULL_THIS_A, WORKARR, MINLON, J
;
;        Note: FULL_THIS_A gets changed in INTERP
;
         IF N_ELEMENTS(FULL_THIS_A) EQ NUMLON THEN BEGIN

           PXFORM = FFT(workarr,-1)
           zonal_spectra(J,dind,0:MAXWAVES-1,A)= PXFORM(0:MAXWAVES-1)
           
            FOR I = 0, KNYQ DO BEGIN

             STAWAVES(J,A,I,0,DIND) =  FLOAT(     PXFORM(I) )
             STAWAVES(J,A,I,1,DIND) = -IMAGINARY( PXFORM(I) )

           ENDFOR ; I LOOP
           
;           print, LATITUDE(J), zonal_spectra(J,A,1,dind)
         
         ENDIF ; N_ELEMENTS(FULL_THIS_A) EQ NUMLON 

        ENDIF ; FULL_THIS_A(0) NE -1 
                           
     ENDFOR  ; A loop

 ENDFOR ; J loop 
 RETURN
 END

;===========================================================================
 PRO GRID_IN_LAT, nev, loct_this_day, tplats_this_day, tplons_this_day,   $
                  uts_this_day, press_this_day, geo_this_day,             $
                  temps_this_day, NUM_IN_BINS, AGEO_REC, ATEMP_REC,       $ 
                  PLOT_SUMMARY
;==========================================================================

 COMMON   ZPARM, ZGRID, NUMZ, ZINC
 COMMON LATPARM, LATITUDE,  NUMLAT, LATINC
 COMMON LONPARM, LONGITUDE, NUMLON, LONINC

 PI2    = 2. * !pi
 SOUTH    = -85

 MINDIFF     = LONINC
 MAXDIFF     = 360. - MINDIFF 

  
 PRINTSCREEN = 1
 PRINTDIAG   = 0
 PAUSE       = ' ' 

 LOCT        = loct_this_day(nev)
 UNVT        = uts_this_day(nev)
 LATI        = tplats_this_day(nev)
 LONI        = tplons_this_day(nev)

; step through in latitude
 FOR J = 0, NUMLAT-1 DO BEGIN 
      UPPERBOUND = LATITUDE(J) + LATINC/2.
      LOWERBOUND = LATITUDE(J) - LATINC/2.
; Choose latitude bounds
      IF (LOWERBOUND LE LATI AND LATI LT UPPERBOUND) THEN BEGIN
; Choose longitude bounds and set a flag
        i        = 0
        rightlon =        loninc/2.
        leftlon  = 360. - loninc/2.
        IF (loni LT rightlon OR loni GE leftlon) THEN I_INDEX = I 
        for i = 1, numlon-1 do begin 
          rightlon = longitude(i) + loninc/2.
          leftlon  = longitude(i) - loninc/2.
          IF (leftlon LE loni AND loni LT rightlon) THEN I_INDEX = I 
        endfor ; i loop

         AGEO   =   geo_this_day(*,nev)
         ATEMPS = temps_this_day(*,nev)
         AVERT = press_this_day(*,nev)

         grid_in_z, ATEMPS, AVERT, tgridded
         grid_in_z, AGEO,   AVERT, ggridded
          
;        NUM_IN_BINS(I_INDEX,J)           = NUM_IN_BINS(I_INDEX,J) + 0
        NUM_IN_BINS(I_INDEX,J)           = NUM_IN_BINS(I_INDEX,J) + 1
;Subscript range values of the form low:high must be >= 0, < size, with low <= high: ATEMP_REC.       
         NIB                              = NUM_IN_BINS(I_INDEX,J)
         ATEMP_REC (*, NIB,I_INDEX,J)     = tgridded
         AGEO_REC (*, NIB,I_INDEX,J)      = ggridded
         if PLOT_SUMMARY then PLOTS,  loni, lati, PSYM = 1, color = 3
    
     ENDIF ; Latitude matching
    
   ENDFOR ; J LOOP

 RETURN
 END

;===============================================================
  PRO grid_in_z,  TVALS, PVALS, AVZ
;===============================================================

  COMMON   ZPARM, ZGRID, NUMZ, ZINC

  PLOTTING    = 0
  PAUSE       = ' '

  MISSING      = -999.
  PREF         = 1000.
  AVZ          = FLTARR(NUMZ)
  NUM_IN_KBIN  = INTARR(NUMZ)
  
  ALTI           = PVALS
  TEMPI          = TVALS
  nt		 = n_elements(ALTI)

  NONZERO       = WHERE (PVALS NE MISSING) 
  ALTI(NONZERO) = -ALOG(PVALS(NONZERO)/PREF) 
;
;  Interpolate to ZGRID
;
   FOR K = 0, NUMZ-1 DO BEGIN

      ALTK          = ZGRID(K)
      UPPERBOUND    = ALTK + ZINC
      LOWERBOUND    = ALTK - ZINC

 ;     print
 ;     print, K, ALTK
 ;     prnit

      FOR n = 0, nt-1 DO BEGIN

	  IF(LOWERBOUND LE ALTI(n) AND ALTI(n) LT UPPERBOUND) THEN BEGIN

          IF(TEMPI(n) NE MISSING) THEN BEGIN
           
             newval	    = TEMPI(n)
             AVZ(K)         = AVZ(K) + newval
             NUM_IN_KBIN(K) = NUM_IN_KBIN(K) + 1
             
            ; print, ALTK, ALTI(n), newval 
            
          ENDIF ; TEMPI(n) NE MISSING) 

        ENDIF ; LOWERBOUND LE ALTI(n) AND ALTI(n) LT UPPERBOUND) 

      ENDFOR ; LOOP OVER ALTI(n)

   ENDFOR ; K LOOP

   FOR K = 0, NUMZ - 1 DO BEGIN
 
     IF NUM_IN_KBIN(K) GT 0 THEN AVZ(K) = AVZ(K)/NUM_IN_KBIN(K) 
;     PRINT, K, ZGRID(K), AVZ(K), NUM_IN_KBIN(K)
       
    ENDFOR ; K LOOP

 RETURN
END

;==================================================================
 PRO GET_MEAN_WIND, DIND, SUM_ZMEAN, UWIND
;==================================================================

 COMMON   ZPARM, ZGRID, NUMZ, ZINC
 COMMON LATPARM, LATITUDE,  NUMLAT, LATINC
 COMMON CORVALS,CORPARM, TANLATA, LATRAD, ACOSPHI

 PRINTSCREEN = 0

 SOUTH    = -85
 JEQ = FIX( (  0. - SOUTH)/LATINC )

 ERAD        = 6.37E6
 GRAV        = 9.8
 SEC_PER_DAY = 60.*60.* 24.
 PI2         = 2.*!pi
 OMEGA       = PI2/SEC_PER_DAY
 OMEGA2      = 2.*OMEGA

 GEWIND      = FLTARR(NUMZ,NUMLAT)
 GRWIND      = FLTARR(NUMZ,NUMLAT)

  K          = NUMZ/2
;
; Locate first and last elements of  SUM_ZMEAN(J,K,DIND) 
;
  FOR J = 0, NUMLAT-1 DO BEGIN 

    IF SUM_ZMEAN(DIND,J,K) NE 0. THEN BEGIN 

      FIRSTJ = J
      GOTO, JUMPOUT1

    ENDIF ; SUM_ZMEAN(J,K,DIND) NE 0. 

  ENDFOR ; J LOOP 
;
; If we're here, no FIRSTJ located
;
  PRINT, ' No FIRSTJ '
  goto,noj
  JUMPOUT1:

  FOR J = NUMLAT-1, 0, -1 DO BEGIN 

    IF SUM_ZMEAN(DIND,J,K) NE 0. THEN BEGIN 

      LASTJ = J
      GOTO, JUMPOUT2

    ENDIF ; SUM_ZMEAN(J,K,DIND) NE 0. 

  ENDFOR ; J LOOP 
;
; If we're here, no LASTJ located
;
  PRINT, ' No LASTJ '
  STOP
  JUMPOUT2:

  FOR K = 0, NUMZ-1 DO BEGIN 
;
;  Get meridional gradient of phi 
;
   J    = FIRSTJ
   U0   = SUM_ZMEAN(DIND,J,K)   ; units are in km. 
   U1   = SUM_ZMEAN(DIND,J+1,K)   ; units are in km. 
   DPHI = LATRAD(J+1) - LATRAD(J)

   IF U1 NE 0. AND U0 NE 0. THEN BEGIN 

     DU   = (U1 - U0)
    GEWIND(K,J) = -DU*GRAV*1.E3/DPHI/ERAD   

   ENDIF ; U1 NE 0. AND U0 NE 0. 

   FOR J = FIRSTJ+1, LASTJ-1 DO BEGIN    

     U_1   = SUM_ZMEAN(DIND,J-1,K)
     U1    = SUM_ZMEAN(DIND,J+1,K)
     DPHI  = LATRAD(J+1) - LATRAD(J-1)

     IF U1 NE 0. AND U_1 NE 0. THEN BEGIN 

      DU   = (U1 - U_1)
      GEWIND(K,J) = -DU*GRAV*1.E3/DPHI/ERAD

     ENDIF ;  U1 NE 0. AND U_1 NE 0. 

   ENDFOR ; J LOOP

   J    = LASTJ
   U1   = SUM_ZMEAN(DIND,J,K)
   U0   = SUM_ZMEAN(DIND,J-1,K)
   DPHI = LATRAD(J) - LATRAD(J-1)

   IF U1 NE 0. AND U0 NE 0. THEN BEGIN 

     DU   = (U1 - U0)
    GEWIND(K,J) = -DU*GRAV*1.E3/DPHI/ERAD

   ENDIF ; U1 NE 0. AND U0 NE 0. 
;
;  Compute geostrophic winds 
;
   FOR J = FIRSTJ, JEQ-1 DO BEGIN

       GEWIND(K,J) = GEWIND(K,J)/CORPARM(J)
       GRWIND(K,J) = CORPARM(J) * GEWIND(K,J)/  $
                     (CORPARM(J) + GEWIND(K,J)*TANLATA(J) )
       UWIND(J,DIND,K) = GRWIND(K,J)

   ENDFOR ; J LOOP

   FOR J = JEQ, JEQ  DO BEGIN 

     U_2   = SUM_ZMEAN(DIND,J-2,K)
      U0   = SUM_ZMEAN(DIND,J,  K)
      U2   = SUM_ZMEAN(DIND,J+2,K)
     DPHI  = LATRAD(J+1) - LATRAD(J-1)

     IF U2 NE 0. AND U_2 NE 0. AND U0 NE 0. THEN BEGIN 

      DU              = U2 + U_2 - 2.*U0
      UG              = -DU*GRAV*1.E3/DPHI/DPHI/ERAD/OMEGA2
      GEWIND(K,J)     = UG
      GRWIND(K,J)     = -DU*GRAV*1.E3/DPHI/DPHI/ERAD/(OMEGA2 + UG/ERAD)
      UWIND(J,DIND,K) = GRWIND(K,J)

     ENDIF ; U2 NE 0. AND U_2 NE 0. AND U0 NE 0.

   ENDFOR ; J LOOP

   FOR J = JEQ+1, LASTJ DO BEGIN

       GEWIND(K,J) = GEWIND(K,J)/CORPARM(J)
       GRWIND(K,J) = CORPARM(J) * GEWIND(K,J)/  $
                     (CORPARM(J) + GEWIND(K,J)*TANLATA(J) )
       UWIND(J,DIND,K) = GRWIND(K,J)

   ENDFOR ; J LOOP

  ENDFOR ; K LOOP
noj:
 RETURN
 END

;==================================================================
 PRO GET_PERTU_WIND, DIND, MAXWAVES, SPECTRA, UWIND, PERTWIND
;==================================================================

 COMMON   ZPARM, ZGRID, NUMZ, ZINC
 COMMON LATPARM, LATITUDE,  NUMLAT, LATINC
 COMMON CORVALS,CORPARM, TANLATA, LATRAD, ACOSPHI

 PRINTSCREEN = 0
 IROOT       = COMPLEX(0,1)

 SOUTH    = -85
 J_60 = FIX( (SOUTH - SOUTH)/LATINC )
  J80 = FIX( (  85. - SOUTH)/LATINC )
  JEQ = FIX( (  0. - SOUTH)/LATINC )

 ERAD        = 6.37E6
 GRAV        = 9.8
 SEC_PER_DAY = 60.*60.* 24.
 PI2         = 2.*!pi
 OMEGA       = PI2/SEC_PER_DAY
 OMEGA2      = 2.*OMEGA

 K          = NUMZ/2
 I          = 1
;
; Locate first and last elements of  SPECTRA
;
  FOR J = 0, NUMLAT-1 DO BEGIN 

    IF SPECTRA(J,DIND,I,K) NE 0. THEN BEGIN 

      FIRSTJ = J
      GOTO, JUMPOUT1

    ENDIF ; SUM_ZMEAN(J,K,DIND) NE 0. 

  ENDFOR ; J LOOP 
;
; If we're here, no FIRSTJ located
;
  PRINT, ' No FIRSTJ '
  goto,noj
  JUMPOUT1:

  FOR J = NUMLAT-1, 0, -1 DO BEGIN 

    IF SPECTRA(J,DIND,I,K)  NE 0. THEN BEGIN 

      LASTJ = J
      GOTO, JUMPOUT2

    ENDIF ; SPECTRA(J,K,I,DIND)  NE 0. 

  ENDFOR ; J LOOP 
;
; If we're here, no LASTJ located
;
  PRINT, ' No LASTJ '
  STOP
  JUMPOUT2:

 FOR I = 1,  MAXWAVES -1 DO BEGIN 

 GEWIND      = COMPLEXARR(NUMZ,NUMLAT)
 GRWIND      = COMPLEXARR(NUMZ,NUMLAT)

 FOR K = 0, NUMZ-1 DO BEGIN 
;
;  Get meridional gradient of phi'
;
   J    = FIRSTJ
   U0   = SPECTRA(J,DIND,I,K)   ; units are in km. 
   U1   = SPECTRA(J+1,DIND,I,K)   ; units are in km. 
   DPHI = LATRAD(J+1) - LATRAD(J)

   IF U1 NE 0. AND U0 NE 0. THEN BEGIN 

     DU   = (U1 - U0)
    GEWIND(K,J) = -DU*GRAV*1.E3/DPHI/ERAD   

   ENDIF ; U1 NE 0. AND U0 NE 0. 

   FOR J = FIRSTJ+1, LASTJ-1 DO BEGIN    

     U_1   = SPECTRA(J-1,DIND,I,K)
     U1    = SPECTRA(J+1,DIND,I,K)
     DPHI  = LATRAD(J+1) - LATRAD(J-1)

     IF U1 NE 0. AND U_1 NE 0. THEN BEGIN 

      DU   = (U1 - U_1)
      GEWIND(K,J) = -DU*GRAV*1.E3/DPHI/ERAD

     ENDIF ;  U1 NE 0. AND U_1 NE 0. 

   ENDFOR ; J LOOP

   J    = LASTJ
   U1   = SPECTRA(J,DIND,I,K)
   U0   = SPECTRA(J-1,DIND,I,K)
   DPHI = LATRAD(J) - LATRAD(J-1)

   IF U1 NE 0. AND U0 NE 0. THEN BEGIN 

     DU   = (U1 - U0)
    GEWIND(K,J) = -DU*GRAV*1.E3/DPHI/ERAD

   ENDIF ; U1 NE 0. AND U0 NE 0. 
;
;  Compute geostrophic winds 
;
   FOR J = FIRSTJ, JEQ-1 DO BEGIN

       UBAR  = UWIND(J,DIND, K)
       DENOM = CORPARM(J) + 2.*UBAR*TANLATA(J) 
       PERTWIND(J,DIND,I,K) =  GEWIND(K,J)/DENOM
;       PRINT, LATITUDE(J), PERTWIND(J,K,I,DIND) 

   ENDFOR ; J LOOP

   FOR J = JEQ, JEQ  DO BEGIN 

      UBAR  = UWIND(J,DIND, K) 
      DENOM = OMEGA + UBAR/ERAD
     U_2   = SPECTRA(J-2,DIND,I,K)
      U0   = SPECTRA(J,DIND,I,K)
      U2   = SPECTRA(J+2,DIND,I,K)
     DPHI  = LATRAD(J+1) - LATRAD(J-1)

     IF U2 NE 0. AND U_2 NE 0. AND U0 NE 0. THEN BEGIN 

       DU                   = U2 + U_2 - 2.*U0
       PERTWIND(J,DIND,I,K) = -DU*GRAV*1.E3/DPHI/DPHI/ERAD/2./DENOM
;       PRINT, LATITUDE(J), PERTWIND(J,K,I,DIND) 

     ENDIF ; U2 NE 0. AND U_2 NE 0. AND U0 NE 0.

   ENDFOR ; J LOOP

   FOR J = JEQ+1, LASTJ DO BEGIN

       UBAR  = UWIND(J,DIND, K)
       DENOM = CORPARM(J) + 2.*UBAR*TANLATA(J) 
       PERTWIND(J,DIND,I,K) =  GEWIND(K,J)/DENOM
;       PRINT, LATITUDE(J), PERTWIND(J,K,I,DIND) 

   ENDFOR ; J LOOP

  ENDFOR ; K LOOP

 ENDFOR ; I LOOP
noj:
 RETURN
 END

;==================================================================
 PRO GET_PERTV_WIND, DIND, MAXWAVES, SPECTRA, UWIND, PERTWIND
;==================================================================

 COMMON   ZPARM, ZGRID, NUMZ, ZINC
 COMMON LATPARM, LATITUDE,  NUMLAT, LATINC
 COMMON CORVALS,CORPARM, TANLATA, LATRAD, ACOSPHI

 PRINTSCREEN = 0
 IROOT       = COMPLEX(0,1)

 SOUTH    = -85
 J_60 = FIX( (SOUTH - SOUTH)/LATINC )
  J80 = FIX( (  85. - SOUTH)/LATINC )
  JEQ = FIX( (  0. - SOUTH)/LATINC )

 ERAD        = 6.37E6
 GRAV        = 9.8
 SEC_PER_DAY = 60.*60.* 24.
 PI2         = 2.*!pi
 OMEGA       = PI2/SEC_PER_DAY
 OMEGA2      = 2.*OMEGA

 K          = NUMZ/2
 I          = 1
;
; Locate first and last elements of  SPECTRA
;
  FOR J = 0, NUMLAT-1 DO BEGIN 

    IF SPECTRA(J,DIND,I,K) NE 0. THEN BEGIN 

      FIRSTJ = J
      GOTO, JUMPOUT1

    ENDIF ; SPECTRA(J,K,I,DIND) NE 0. 

  ENDFOR ; J LOOP 
;
; If we're here, no FIRSTJ located
;
  PRINT, ' No FIRSTJ '
goto,noj
  JUMPOUT1:

  FOR J = NUMLAT-1, 0, -1 DO BEGIN 

    IF SPECTRA(J,DIND,I,K)  NE 0. THEN BEGIN 

      LASTJ = J
      GOTO, JUMPOUT2

    ENDIF ; SPECTRA(J,K,I,DIND)  NE 0. 

  ENDFOR ; J LOOP 
;
; If we're here, no LASTJ located
;
  PRINT, ' No LASTJ '
  STOP
  JUMPOUT2:

 FOR I = 1,  MAXWAVES -1 DO BEGIN 

 GEWIND      = COMPLEXARR(NUMZ,NUMLAT)
 YSHEAR      = COMPLEXARR(NUMZ,NUMLAT)

 FOR K = 0, NUMZ-1 DO BEGIN 
;
;  Get longitudinal gradient of phi'
;
   FOR J = 0, NUMLAT-1 DO BEGIN    

     U1    = SPECTRA(J,DIND,I,K)
     IF U1 NE 0. THEN $
     GEWIND(K,J) = IROOT*FLOAT(I)*U1*GRAV*1.E3/ACOSPHI(J)

   ENDFOR ; J LOOP
;
; Get meridional gradient of Ubar
;
   J    = FIRSTJ
   U0   = UWIND(J,DIND, K)   
   U1   = UWIND(J+1,DIND, K)   
   DPHI = LATRAD(J+1) - LATRAD(J)

   IF U1 NE 0. AND U0 NE 0. THEN $
    YSHEAR(K,J) = (U1 - U0)/DPHI/ERAD   

   FOR J = FIRSTJ+1, LASTJ-1 DO BEGIN    

     U_1   = UWIND(J-1,DIND, K)
     U1    = UWIND(J+1,DIND, K)
     DPHI  = LATRAD(J+1) - LATRAD(J-1)
     IF U1 NE 0. AND U_1 NE 0. THEN $
      YSHEAR(K,J) = (U1 - U_1)/DPHI/ERAD

   ENDFOR ; J LOOP

   J    = LASTJ
   U1   = UWIND(J, DIND, K)
   U0   = UWIND(J-1,DIND, K)
   DPHI = LATRAD(J) - LATRAD(J-1)
   IF U1 NE 0. AND U0 NE 0. THEN $
     YSHEAR(K,J) = (U1 - U0)/DPHI/ERAD
;
;  Compute gradient winds 
;
   FOR J = FIRSTJ, JEQ-1 DO BEGIN

       UBAR   = UWIND(J,DIND,K)
       UBARY  = YSHEAR(K,J)
       DENOM = -( CORPARM(J) - UBARY + UBAR*TANLATA(J) )
       PERTWIND(J,DIND,I,K) =  -GEWIND(K,J)/DENOM

   ENDFOR ; J LOOP

   FOR J = JEQ+1, LASTJ DO BEGIN

       UBAR  = UWIND(J,DIND,K)
       UBARY  = YSHEAR(K,J)
       DENOM = -( CORPARM(J) - UBARY + UBAR*TANLATA(J) )
       PERTWIND(J,DIND,I,K) =  -GEWIND(K,J)/DENOM

   ENDFOR ; J LOOP

  ENDFOR ; K LOOP

 ENDFOR ; I LOOP
noj:
 RETURN
 END

;========================================================================
;PRO LAT_HT, DATARR, BASEVAL, CONINC, ZEROVAL, PLOTLABEL,DIND,daystring,yearstring,plotlabel_date
PRO LAT_HT, DATARR, BASEVAL, CONINC, ZEROVAL, PLOTLABEL
;========================================================================

COMMON   ZPARM, ZGRID, NUMZ, ZINC
COMMON LATPARM, LATITUDE,  NUMLAT, LATINC
COMMON CORVALS,CORPARM, TANLATA, LATRAD, ACOSPHI

SPECIAL_VALUE = 999.
!P.FONT     = 17
!P.CHARSIZE = 1.5
ZTOP       =   90.   ; 100.
ZBOT       =   20.
H = 6.5
KLOWER     = ROUND(( ZBOT/H - 2.)/ZINC)
KUPPER     = ROUND(( ZTOP/H - 2.)/ZINC)
LOWER_Z    = ZGRID(KLOWER)
UPPER_Z    = ZGRID(KUPPER)
NUM_CON  =   21
CONVALS  =   FLTARR(num_con)
CONLINES =   INTARR(num_con)
CONTHICK =   INTARR(num_con)
CONLABS  =   INTARR(num_con)
CVEC     =   INTARR(num_con)
THICKVAL = 50.

FOR I = 0, NUM_CON-1L DO BEGIN
    CONVALS(I)  = BASEVAL + I*CONINC
    CONLINES(I) = 1
    CONLABS(I)  = 1
    CONTHICK(I) = 1
ENDFOR

CVEC(*)  =   255
ZERINDEX = WHERE(CONVALS EQ ZEROVAL)

IF ZERINDEX(0) NE -1 THEN BEGIN
   CONTHICK( ZERINDEX) = 3
   FOR I = ZERINDEX(0), NUM_CON-1L DO CONLINES(I) = 0
ENDIF

THICKINDEX = WHERE(CONVALS EQ THICKVAL)
IF THICKINDEX(0) NE -1 THEN FOR I = THICKINDEX(0), NUM_CON-1 DO CVEC(I) = 220

THICKINDEX = WHERE(CONVALS EQ -50.)
IF THICKINDEX(0) NE -1 THEN FOR I = 0, THICKINDEX(0) -1 DO CVEC(I) = 220

CONARR = DATARR
MISSING = WHERE(CONARR EQ 0.)
IF MISSING(0) NE -1 THEN CONARR(MISSING) = SPECIAL_VALUE

;date_label=label_date(date_format=['%M!C%D'],offset=-1)

loadct,13
set_plot,'ps'
device, filename='EP_fluxes_'+yearstring+daystring+'.ps',/color

nlev=10
bin=(-200)/nlev
userlev=indgen(nlev)*bin + 60
CONLAB = ['400','800','1200','1600','2000']
CONTOUR, CONARR, LATITUDE, ZGRID, LEVELS = convals, C_LINESTYLE = CONLINES, C_THICK = CONTHICK, $
         C_LABELS = CONLABS, C_CHARSIZE = 1.25, /FOLLOW, MAX_V = 120, MIN_VALUE=-200, thick=3, $
         yrange = [LOWER_Z, UPPER_Z], Xrange = [18, 62], $
         ; Xrange = [25., 85.],        $
         ;TITLE  = PLOTLABEL + plotlabel_date+daystring+', '+yearstring, $
         TITLE  = PLOTLABEL , XTITLE = 'DAY NUMBER', YTITLE = 'SCALED HEIGHT', $
         XSTYLE=1, YSTYLE=1, YTYPE = 0
OPLOT, [0,0], [3.,14.], LINESTYLE = 0

device,/close
RETURN
END

;=================================================================
 PRO LAT_LON_WIND, UWINDS, VWINDS, GEOPS, DIND, KPLOT,daystring,yearstring,plotlabel_date
;=================================================================

COMMON   ZPARM, ZGRID, NUMZ, ZINC
COMMON LATPARM, LATITUDE,  NUMLAT, LATINC
COMMON LONPARM, LONGITUDE, NUMLON, LONINC

PAUSE    = ' '
SPECIAL_VALUE = 9999.
STWAV  = 1
ENDWAV = 4
y0  =  25.
y1  =  85.
x0  =   0.
x1  =  360.
SOUTH    = -85.
JEQ     = FIX( ( y0 - SOUTH)/LATINC)
J84     = FIX( ( y1 - SOUTH)/LATINC)
K = kplot
;
;  Making a color table
;     0     1  2  3  4   5
;   Black   G  B  R  Y  Magenta
;
Red = [0, 0, 0, 1, 1, 1]
Green=[0, 1, 0, 0, 1, 0]
Blue =[0, 0, 1, 0, 0, 1]

Tvlct, 255*Red, 255*Green, 255*Blue

YASP=20.
XASP=YASP*!d.x_vsize*(!x.window(1)-!x.window(0)) / $
	  (!d.y_vsize*(!y.window(1)-!y.window(0)))
;
; Arrow heads
;
r=0.3				; arrow head length
angle = 22.5 * !dtor		; Angle of arrowhead
st = r * sin(angle)		; sin 22.5 degs * length of head
ct = r * cos(angle)
;
;  Define plotting space and contour the data
;
!P.FONT = 17
!P.CHARSIZE = 1.5
BASEVAL  = -50.  ; Temperature
CONINC   =   5.
;BASEVAL  = -2000.  ;Geopotential
;CONINC   =   200.
ZEROVAL  =     0.
NUM_CON  =   20
CONVALS  =   FLTARR(NUM_CON)
CONLINES =   INTARR(NUM_CON)
CONTHICK =   INTARR(NUM_CON)
CONLABS  =   INTARR(NUM_CON)

FOR I = 0, NUM_CON - 1 DO BEGIN
    CONVALS(I)  = BASEVAL + I*CONINC
    CONLINES(I) = 1
    CONLABS(I)  = 1
    CONTHICK(I) = 1
ENDFOR ; I LOOP

ZERINDEX = WHERE(CONVALS EQ ZEROVAL)
IF ZERINDEX(0) NE -1 THEN BEGIN
   CONTHICK( ZERINDEX) = 3
   FOR I = ZERINDEX(0), NUM_CON-1 DO CONLINES(I) = 0
ENDIF

ucon = fltarr(numlon,numlat)
vcon = fltarr(numlon,numlat)
ucon(*,jeq:j84) = uwinds(K,*,jeq:j84,dind)
vcon(*,jeq:j84) = vwinds(K,*,jeq:j84,dind)

full      = where(vcon ne 0.)
if full(0) ne -1 then begin
   maxmag = max( sqrt(ucon(full)^2. + vcon(full)^2. ) )
   length = 1.0
   xfact = length * 360. / XASP / maxmag
   yfact = length * 180. / YASP / maxmag
   gcon      =  FLTARR(NUMLON, NUMLAT)
   gcon(*,*) = geops(K,*,*,dind)
   avg_vcon=fltarr(numlon,numlat)
   avg_gcon=fltarr(numlon,numlat)

   openw,100,yearstring+daystring+'_zgrid20_averages.txt'
   printf,100,' latitude','      varv','      varg','      Cxy'
   vgcon = fltarr(numlon)
   for lat_int=0,numlat-1 do begin
; Sample zonal means
     avg_vcon = (total(vcon(0:numlon-1,lat_int)))/numlon
     avg_gcon = (total(gcon(0:numlon-1,lat_int)))/numlon
; VG Product and zonal mean
     vgcon(0:numlon-1) = (vcon(0:numlon-1,lat_int)-avg_vcon)*(gcon(0:numlon-1,lat_int)-avg_gcon)
     avg_vgcon = total(vgcon)/float(numlon)
; Sample zonal variances
     varv = total((vcon(0:numlon-1,lat_int)-avg_vcon)^2)/float(numlon-1)
     varg = total((gcon(0:numlon-1,lat_int)-avg_gcon)^2)/float(numlon-1)
; Zonal mean cross-correlation
     Cxy=avg_vgcon/(sqrt(varv*varg))
     print, lat_int, latitude(lat_int), avg_vcon, avg_gcon, avg_vgcon, varv, varg, cxy
     printf,100,latitude(lat_int),' ',varv,' ',varg,' ',Cxy
   endfor
   close,100

   MISSING = WHERE(gcon eq 0.)
   IF MISSING(0) NE -1 THEN GCON(MISSING) =   SPECIAL_VALUE

   loadct,13
   set_plot,'ps'
   device, filename='zgrid20_tempwinds_'+yearsting+daystring+'.ps',/color

   CONTOUR, GCON(*,JEQ:J84), LONGITUDE, LATITUDE(JEQ:J84), XSTYLE=1,YSTYLE=1,$
            YRANGE = [y0, y1], XRANGE = [x0, x1], $
            LEVELS = CONVALS, C_LINESTYLE = CONLINES, C_THICK = CONTHICK, $
            C_LABELS = CONLABS, MAX_VALUE=9998., /FOLLOW ,$
            XTITLE = 'LONGITUDE', YTITLE = 'LATITUDE',$
            TITLE = 'TEMP & WINDS'+' 5 SCALE HEIGHTS:'+plotlabel_date+daystring+', '+yearstring

   FOR J = JEQ, J84 DO BEGIN
       lat_to_plot   = LATITUDE(J)
       PRINT, LATITUDE(J), GCON(*,J)
       FOR I = 0, NUMLON-1 DO BEGIN
           lon_to_plot   = LONGITUDE(I)
           uval      = UWINDS(K,I,J,DIND)
           vval      = VWINDS(K,I,J,DIND)
           if vval ne 0. and uval ne 0. then begin
              xx0    =  lon_to_plot
              dx     =  xfact * uval
              xx1    =  xx0 + dx
              yy0    = lat_to_plot
              dy     = yfact * vval
              yy1    = yy0 + dy

; plots,[xx0,xx1, xx1-(ct*dx-st*dy),xx1,xx1-(ct*dx+st*dy)], $
;       [yy0,yy1,yy1-(ct*dy+st*dx),yy1,yy1-(ct*dy-st*dx)], color = 3
; print, [xx0,xx1, xx1-(ct*dx-st*dy),xx1,xx1-(ct*dx+st*dy)],$
;        [yy0,yy1,yy1-(ct*dy+st*dx),yy1,yy1-(ct*dy-st*dx)]

              arrow, xx0,yy0,xx1,yy1,/data,color=255,thick=2	 	;hsize=2, color=3, thick=3.0
           endif ; vval ne 0.
       ENDFOR ; I LOOP
   ENDFOR ; J LOOP
endif ; full(0) ne -1
device,/close
RETURN
END

;===========================================================================
  PRO BUILD_WIND_PROFILES, PERTWIND, DIND, WIND_PROFILES 
;===========================================================================

COMMON   ZPARM, ZGRID, NUMZ, ZINC
COMMON LATPARM, LATITUDE,  NUMLAT, LATINC
COMMON LONPARM, LONGITUDE, NUMLON, LONINC

PAUSE    = ' '
STWAV  = 1
ENDWAV = 5
IROOT     = COMPLEX(0,1)
LONRAD    = LONGITUDE*!dtor

for J = 0, NUMLAT-1 DO BEGIN
    for K = 0, NUMZ-1 DO BEGIN
    IF PERTWIND(J,DIND,1,K) NE 0. THEN BEGIN
       FOR I = 0, NUMLON-1 DO BEGIN
           FOR IK = STWAV, ENDWAV DO BEGIN
               CARG = IROOT*IK*LONRAD(I)
               C1 = PERTWIND(J,DIND,IK,K)*EXP(CARG)
               C_1  = CONJ( C1 )
               WIND_PROFILES(DIND,I,J,K) = WIND_PROFILES(DIND,I,J,K) + FLOAT( C1 + C_1 )
           ENDFOR ; IK LOOP
       ENDFOR ; I LOOP
    ENDIF ; PERTWIND(J,K,1,DIND) NE 0.
    ENDFOR ; K LOOP
endfor ; j loop

RETURN
END

;========================================================================
 PRO TIME_HT, stday, endday, DATARR, M, J, BASEVAL, CONINC, PLOTLABEL
;========================================================================

COMMON   ZPARM, ZGRID, NUMZ, ZINC
COMMON LATPARM, LATITUDE,  NUMLAT, LATINC
COMMON CORVALS,CORPARM, TANLATA, LATRAD, ACOSPHI
SPECIAL_VALUE = 9999.
!P.FONT     = 17
!P.CHARSIZE = 1.5
ZTOP       =   70.
ZBOT       =   20.
H = 6.5
KLOWER     = ROUND(( ZBOT/H - 2.)/ZINC)
KUPPER     = ROUND(( ZTOP/H - 2.)/ZINC)
LOWER_Z    = ZGRID(KLOWER)
UPPER_Z    = ZGRID(KUPPER)
NUM_CON  =   21
CONVALS  =   FLTARR(NUM_CON)
CONLINES =   INTARR(NUM_CON)
CONTHICK =   INTARR(NUM_CON)
CONLABS  =   INTARR(NUM_CON)
CVEC     =   INTARR(NUM_CON)
ZEROVAL =  0.
THICKVAL = 1500.

FOR I = 0, NUM_CON-1L DO BEGIN
    CONVALS(I)  = BASEVAL + I*CONINC
    CONLINES(I) = 1
    CONLABS(I)  = 1
    CONTHICK(I) = 1
ENDFOR

CVEC(*)  =   255
ZERINDEX = WHERE(CONVALS EQ ZEROVAL)

IF ZERINDEX(0) NE -1 THEN BEGIN
   CONTHICK( ZERINDEX) = 3
   FOR I = ZERINDEX(0), NUM_CON-1L DO CONLINES(I) = 0
ENDIF

THICKINDEX = WHERE(CONVALS EQ THICKVAL)
IF THICKINDEX(0) NE -1 THEN FOR I = THICKINDEX(0), NUM_CON-1 DO CVEC(I) = 220
THICKINDEX = WHERE(CONVALS EQ -50.)
IF THICKINDEX(0) NE -1 THEN FOR I = 0, THICKINDEX(0) -1 DO CVEC(I) = 220 

numdays = endday - stday + 1
days = findgen(numdays)+stday
CONARR = FLTARR(NUMDAYS,NUMZ)

FOR DIND = 0, NUMDAYS-1 DO CONARR(DIND,*) = 2.*ABS( DATARR(*,M,J,DIND) )
MISSING = WHERE(CONARR EQ 0.)
IF MISSING(0) NE -1 THEN CONARR(MISSING) = SPECIAL_VALUE

CONTOUR, CONARR, DAYS, ZGRID, LEVELS = CONVALS, C_LINESTYLE = CONLINES, C_THICK = CONTHICK, $
         C_LABELS = CONLABS, C_CHARSIZE = 1.25, /FOLLOW, MAX_V = SPECIAL_VALUE-1, thick=3, $
         yrange = [LOWER_Z, UPPER_Z], Xrange = [-80., 80.], TITLE  = PLOTLABEL, YTITLE = 'Scaled height', $
         XSTYLE=1, YSTYLE=1, YTYPE = 0

OPLOT, [0,0], [3.,12.], LINESTYLE = 0

PAUSE    = ' '
READ, PAUSE

RETURN
END

;====================================================
 PRO INTERP, FULL, AUX, MINFULL, J
;====================================================

  COMMON LONPARM, LONGITUDE, NUMLON, LONINC
  COMMON LATPARM, LATITUDE,  NUMLAT, LATINC

  FIRST_LON = 0
   LAST_LON = NUMLON-1

  NUMFULL = N_ELEMENTS(FULL)

  IF NUMFULL GE MINFULL THEN BEGIN 

   FOR N = 1, NUMFULL-1 DO BEGIN

     LAST_FULL_INDEX = FULL(N-1)
     THIS_FULL_INDEX = FULL(N)

     LAST_FULL_VALUE = AUX(LAST_FULL_INDEX)
     THIS_FULL_VALUE = AUX(THIS_FULL_INDEX)

     DIFF    = THIS_FULL_INDEX - LAST_FULL_INDEX
     DIFFVAL = THIS_FULL_VALUE - LAST_FULL_VALUE

     IF DIFF GT 1 AND DIFF LT 3 THEN BEGIN 

      SLOPE = DIFFVAL/DIFF

      FOR NN = LAST_FULL_INDEX+1, THIS_FULL_INDEX-1 DO $ 
       AUX(NN) = LAST_FULL_VALUE + SLOPE*(NN-LAST_FULL_INDEX)

     ENDIF ; DIFF GT 1

   ENDFOR ; N LOOP
;
;  What if the missing element is the final longitude? 
; 
   MISSING  = WHERE(AUX EQ 0.) 

   IF MISSING(0) NE -1 THEN BEGIN 

    LAST_IS_MISSING = WHERE(MISSING EQ LAST_LON)

    IF LAST_IS_MISSING(0) NE -1 THEN BEGIN
;
;    Need to "wrap" the interpolation 
;
     LAST_FULL_INDEX = LAST_LON - 1
     THIS_FULL_INDEX = FIRST_LON

     LAST_FULL_VALUE = AUX(LAST_FULL_INDEX)
     THIS_FULL_VALUE = AUX(THIS_FULL_INDEX)

     IF LAST_FULL_VALUE NE 0. AND THIS_FULL_VALUE NE 0. THEN BEGIN 

      DIFF    = 2
      DIFFVAL = THIS_FULL_VALUE - LAST_FULL_VALUE
      SLOPE = DIFFVAL/DIFF
      AUX(LAST_LON) = LAST_FULL_VALUE + SLOPE * 1. 
      PRINT, LATITUDE(J), ' Last longitude interpolated'    

     ENDIF ; IF LAST_FULL_VALUE NE 0. AND THIS_FULL_VALUE NE 0. 

    ENDIF ; LAST_IS_MISSING(0) NE -1 

   ENDIF ;  MISSING(0) NE -1 
;
;  What if the missing element is the first longitude? 
; 
   MISSING  = WHERE(AUX EQ 0.) 

   IF MISSING(0) NE -1 THEN BEGIN 

    FIRST_IS_MISSING = WHERE(MISSING EQ FIRST_LON)

    IF FIRST_IS_MISSING(0) NE -1 THEN BEGIN
;
;    Need to "wrap" the interpolation 
;
     LAST_FULL_INDEX = LAST_LON 
     THIS_FULL_INDEX = FIRST_LON + 1

     LAST_FULL_VALUE = AUX(LAST_FULL_INDEX)
     THIS_FULL_VALUE = AUX(THIS_FULL_INDEX)

     IF LAST_FULL_VALUE NE 0. AND THIS_FULL_VALUE NE 0. THEN BEGIN 

      DIFF    = 2
      DIFFVAL = THIS_FULL_VALUE - LAST_FULL_VALUE
      SLOPE = DIFFVAL/DIFF
      AUX(FIRST_LON) = LAST_FULL_VALUE + SLOPE * 1. 
      PRINT, LATITUDE(J), ' First longitude interpolated'    

     ENDIF ; IF LAST_FULL_VALUE NE 0. AND THIS_FULL_VALUE NE 0. 

    ENDIF ; FIRST_IS_MISSING(0) NE -1 

   ENDIF ;  MISSING(0) NE -1 

  ENDIF ; NUMFULL GT 20

  FULL  = WHERE(AUX NE 0.) 

 RETURN
 END

;========================================================================
 pro read_netCDF, data, filename, attributes, status
;========================================================================
;
; NAME:
;	read_netCDF.pro
;
; PURPOSE:
;	Read netCDF file into structure variable
;
; CATEGORY:
;	All levels of processing
;
; CALLING SEQUENCE:  =

;	read_netCDF, filename, data, attributes, status
;
; INPUTS:
;	filename = filename for existing netCDF file
;
; OUTPUTS:  =
;
;	data       = structure variable for data read from netCDF file
;	attributes = array of strings of the attributes from the netCDF file
;	status     = result status:      0 = OK_STATUS
;                                       -1 = BAD_PARAMS
;                                       -2 = BAD_FILE
;			                -3 = BAD_FILE_DATA
;                                       -4 = FILE_ALREADY_OPENED
;
; COMMON BLOCKS:
;	None
;
; PROCEDURE:
;	Check for valid input parameters
;
;	Open the netCDF file
;
;	Create structures based on the netCDF definitions
;
;	Once structures are defined, then read the netCDF variables 
;       into the structure's data
;
;	Read the attributes into a string array
;
;	Close the netCDF file
;
;	NetCDF IDL Procedures / Process:
;
;	1.	NCDF_OPEN: Open an existing netCDF file.
;
;	2.	NCDF_INQUIRE: Call this function to find the format of the netCDF file.
;
;	3.	NCDF_DIMINQ: Retrieve the names and sizes of dimensions in the file.
;
;	4.	NCDF_VARINQ: Retrieve the names, types, and sizes of variables in the file.
;
;	5.	NCDF_ATTINQ: Optionally, retrieve the types and lengths of attributes.
;
;	6.	NCDF_ATTNAME: Optionally, retrieve attribute names.
;
;	7.	NCDF_ATTGET: Optionally, retrieve the attributes.
;
;	8.	NCDF_VARGET: Read the data from the variables.
;
;	9.	NCDF_CLOSE: Close the file.
;
; MODIFICATION HISTORY:
;
;	9/20/99		Tom Woods		Original release of code, Version 1.00
;
;
;	Generic "status" values
;
;
OK_STATUS           =  0
BAD_PARAMS          = -1
BAD_FILE            = -2
BAD_FILE_DATA       = -3
FILE_ALREADY_OPENED = -4

debug_mode =  0 ; set to 1 if want to debug this procedure
;
; Check for valid parameters
;
status = OK_STATUS

if (debug_mode gt 2) and ( !d.name eq 'MAC' ) then begin
	SEE_MAC_CODE = !dir + ':SEE DPS =C4:'
	full_file = SEE_MAC_CODE + 'see_data:' + filename
endif else begin
	full_file = filename
endelse
;
; Open the netCDF file
;
; 1.	NCDF_OPEN: Open an existing netCDF file.
;
if (debug_mode gt 0) then print, 'Opening ', filename, ' ...'
fid = NCDF_OPEN( full_file, /NOWRITE )
;
; Create structures based on the netCDF definitions
;
; 2.	NCDF_INQUIRE: Call this function to find the format of the netCDF file.
;
; 3.	NCDF_DIMINQ: Retrieve the names and sizes of dimensions in the file.
;
; 4.	NCDF_VARINQ: Retrieve the names, types, and sizes of variables in the file.
;
finq = NCDF_INQUIRE( fid )	; finq /str = ndims, nvars, ngatts, recdim
;
; get dimension definitions first
; get unlimited dimension (finq.recdim) =

dim_unlimited = finq.recdim	 ; = -1 if undefined, otherwise index into dim array

if ( finq.ndims gt 0 ) then begin

   dimstr = ' '
   dimsize = 0L
   dim_name = strarr( finq.ndims )
   dim_size = lonarr( finq.ndims )

   for k=0,finq.ndims-1 do begin

    NCDF_DIMINQ, fid, k, dimstr, dimsize
    dim_name[k] = dimstr
    dim_size[k] = dimsize

   endfor

endif
;
; Get variable definitions next
; Also determine nested structure levels, max. dimension, and command dim=ension value
;
; LIMITATION: 6 dimensions allowed per variable
;	      netCDF does not really define unsigned variable types
;
;  Have internal structure definition for tracking variables / structures
;
;	name       = name from netCDF file
;	var_name   = name from structure definition (last word after last '.')
;	type       = data type value (same values as used by size())
;	natts      = number of attributes for this variable
;	ndims      = number of dimensions in "dim"
;	dim        = dimension index into dim_size[]
;	nest_level = nest level of structures (number of '.' in name)
;	nest_name  = structure name (nested)
;	nest_id    = index to first case of structure name (nested)
;	nest_cnt   = index of variable within a single structure (nested)
;	ptr        = data variable pointer
;	str_ptr    = structure pointer (if first case of new structure)
;
;
var_inq1 = { name : " ", var_name : " ", type : 0, natts : 0L, ndims : 0L, $
               dim: lonarr(8), nest_level : 0, nest_name: strarr(6),           $
               nest_id : lonarr(6), nest_cnt : lonarr(6), ptr : PTR_NEW(),     $
               str_ptr : PTRARR(6) }

var_inq   = replicate( var_inq1, finq.nvars )
max_level = 0		; track max structure nest level while getting variable definitions
max_dim = 1		; track max base structure dimension required
has_common_dim = 1	; assume TRUE to start out, any conflict makes it FALSE
;
; Sort out first the dimensions and attribute numbers
; Check for max. dim needed for base structure
;       and if should have base structure array (if all the same last dim)
;
for k=0, finq.nvars-1 do begin

  var_def = NCDF_VARINQ( fid, k )
  var_inq[k].ndims = var_def.ndims
  var_inq[k].natts = var_def.natts

  if (var_def.ndims gt 0) then begin

   for j=0, var_def.ndims-1 do var_inq[k].dim[j] = var_def.dim[j]

  endif

  if (var_def.ndims gt 0) then begin

    lastdim = dim_size[ var_def.dim[var_def.ndims-1] ]
    if (lastdim gt max_dim) then max_dim = lastdim
    if (var_inq[k].dim[var_inq[k].ndims-1] ne $
      var_inq[0].dim[var_inq[0].ndims-1]) then $
       has_common_dim = 0

  endif else has_common_dim = 0

endfor  ; k loop

if (debug_mode gt 0) then begin

  print, ' '
  if (has_common_dim) then $
    print, 'Array dimension for base structure = ', strtrim(max_dim, 2) $
  else print, $
   'Single structure element will be defined - max dim. seen though is ', $
    strtrim(max_dim, 2)

endif

if (has_common_dim eq 0) then $
   max_dim    = 1    ;  make single-element structure only
str_dim_limit = 1    ; define limit for converting BYTE array into STRING

if (has_common_dim) then str_dim_limit = 2
;
; Now define variables
;
for k=0, finq.nvars-1 do begin

  var_def = NCDF_VARINQ( fid, k )
  var_inq[k].name = var_def.name

  case strupcase(var_def.datatype) of

    'BYTE': begin

	    theType = 1		; use size() definitions for data type numbers
	    if (var_def.ndims ge str_dim_limit) then begin

 	     if (debug_mode gt 0) then $
                print, 'Forcing STRING type for ', var_def.name
	     theType = 7
			
            endif
			
            end
		
     'CHAR': begin
			
             theType = 7		; expect STRING type
	     if (debug_mode gt 0) then print, 'STRING type for ', var_def.name
			
             end
		
      'SHORT': theType = 2
      'LONG': theType = 3
      'DOUBLE': theType = 5
		
   else: theType = 4		; default is FLOAT
	
  endcase
;
; Set up structure variable definitions, assume nest level 0 before looking for '.'
; Increase nest_level for each '.' found and fill in nest_name, nest_id[], nest_cnt[]
;
  var_inq[k].type = theType
  var_inq[k].nest_level = 0

  for ii=0,5 do begin

    var_inq[k].nest_name[ii] = ''
    var_inq[k].nest_id[ii]   = 0
    var_inq[k].nest_cnt[ii]  = 0

  endfor

  var_inq[k].nest_id[0] = 0
  if (k eq 0) then var_inq[k].nest_cnt[0] = 0 $
  else var_inq[k].nest_cnt[0] =  var_inq[k-1].nest_cnt[0] + 1
  dotpos = 0
  
  while (dotpos ge 0) do begin

   lastpos = dotpos
   dotpos = strpos( var_def.name, '.', lastpos )
    
   if (dotpos ge 0) then begin

    var_inq[k].nest_level = var_inq[k].nest_level + 1
    nn = var_inq[k].nest_level
    if (nn gt max_level) then max_level = nn

    if (nn gt 5) then begin

     print, 'ERROR: write_netCDF can not handle more than 4 nested structures !'
     print, 'Aborting...'
     NCDF_CONTROL, fid, /ABORT
     status = BAD_FILE_DATA
     return
			
    endif
			
    newname = strmid(var_def.name, lastpos, dotpos-lastpos)
    var_inq[k].nest_name[nn] = newname
    if (k eq 0) then k1=0 else k1 = k - 1
   
    if (k ne 0) and ( var_inq[k1].nest_level ge nn ) and $
                     (var_inq[k1].nest_name[nn] eq newname) then begin

      var_inq[k].nest_cnt[nn-1] = var_inq[k].nest_cnt[nn-1] - 1
      var_inq[k].nest_id[nn] = var_inq[k1].nest_id[nn]
      var_inq[k].nest_cnt[nn] = var_inq[k1].nest_cnt[nn] + 1

    endif else begin
				
     var_inq[k].nest_id[nn] = k
     var_inq[k].nest_cnt[nn] = 0

    endelse
			
    dotpos = dotpos + 1
		
   endif
	
  endwhile

  var_inq[k].var_name = strmid( var_def.name, lastpos, strlen(var_def.name) - lastpos )
;
;  Now define variable and save as PTR
;  uses dumb dimension rules : 
;  ndim_var = ndim_total - 1				for base structure being an array
;  if (CHAR) then ndim_var = ndim_var - 1		for string definitions
;
  ndim_array = var_inq[k].ndims
  if (has_common_dim) then ndim_array = ndim_array - 1
  if (var_inq[k].type eq 7) then ndim_array = ndim_array - 1
  if (ndim_array lt 0) then ndim_array = 0

  case ndim_array of

    0:	begin

        case var_inq[k].type of 

	 1: theData = 0B
	 2: theData = 0
	 3: theData = 0L
	 5: theData = 0.0D0
	 7: theData = ''
      else: theData = 0.0

        endcase

       end

    1: begin

         case var_inq[k].type of 

	   1: theData = bytarr( dim_size[ var_inq[k].dim[0] ] )
	   2: theData = intarr( dim_size[ var_inq[k].dim[0] ] )
	   3: theData = lonarr( dim_size[ var_inq[k].dim[0] ] )
	   5: theData = dblarr( dim_size[ var_inq[k].dim[0] ] )
	   7: theData = strarr( dim_size[ var_inq[k].dim[1] ] )	; offset 1 Dim for char array
        else: theData = fltarr( dim_size[ var_inq[k].dim[0] ] )

         endcase
			
         end

    2:  begin
	
         case var_inq[k].type of 

           1: theData = bytarr( dim_size[ var_inq[k].dim[0] ], dim_size[ var_inq[k].dim[1] ] )
	   2: theData = intarr( dim_size[ var_inq[k].dim[0] ], dim_size[ var_inq[k].dim[1] ] )
	   3: theData = lonarr( dim_size[ var_inq[k].dim[0] ], dim_size[ var_inq[k].dim[1] ] )
	   5: theData = dblarr( dim_size[ var_inq[k].dim[0] ], dim_size[ var_inq[k].dim[1] ] )
	   7: theData = strarr( dim_size[ var_inq[k].dim[1] ], dim_size[ var_inq[k].dim[2] ] )
        else: theData = fltarr( dim_size[ var_inq[k].dim[0] ], dim_size[ var_inq[k].dim[1] ] )
			
         endcase
			
          end
	
     3:  begin

          case var_inq[k].type of 

	   1: theData = bytarr( dim_size[ var_inq[k].dim[0] ], dim_size[ var_inq[k].dim[1] ], $
				dim_size[ var_inq[k].dim[2] ] )
	   2: theData = intarr( dim_size[ var_inq[k].dim[0] ], dim_size[ var_inq[k].dim[1] ], $
  				dim_size[ var_inq[k].dim[2] ]  )
	   3: theData = lonarr( dim_size[ var_inq[k].dim[0] ], dim_size[ var_inq[k].dim[1] ], $
				dim_size[ var_inq[k].dim[2] ]  )
	   5: theData = dblarr( dim_size[ var_inq[k].dim[0] ], dim_size[ var_inq[k].dim[1] ], $
				dim_size[ var_inq[k].dim[2] ]  )
	   7: theData = strarr( dim_size[ var_inq[k].dim[1] ], dim_size[ var_inq[k].dim[2] ], $
				dim_size[ var_inq[k].dim[3] ]  ) ; offset 1 Dim for char array
        else: theData = fltarr( dim_size[ var_inq[k].dim[0] ], dim_size[ var_inq[k].dim[1] ], $
				dim_size[ var_inq[k].dim[2] ]  )
	  endcase
	
            end
		
     4:  begin
	
            case var_inq[k].type of 
	
            1: theData = bytarr( dim_size[ var_inq[k].dim[0] ], dim_size[ var_inq[k].dim[1] ], $
            	                 dim_size[ var_inq[k].dim[2] ], dim_size[ var_inq[k].dim[3] ] )
	    2: theData = intarr( dim_size[ var_inq[k].dim[0] ], dim_size[ var_inq[k].dim[1] ], $
                                 dim_size[ var_inq[k].dim[2] ], dim_size[ var_inq[k].dim[3] ]  )
	    3: theData = lonarr( dim_size[ var_inq[k].dim[0] ], dim_size[ var_inq[k].dim[1] ], $
                                 dim_size[ var_inq[k].dim[2] ], dim_size[ var_inq[k].dim[3] ]  )
	    5: theData = dblarr( dim_size[ var_inq[k].dim[0] ], dim_size[ var_inq[k].dim[1] ], $
                                 dim_size[ var_inq[k].dim[2] ], dim_size[ var_inq[k].dim[3] ]  )
	    7: theData = strarr( dim_size[ var_inq[k].dim[1] ], dim_size[ var_inq[k].dim[2] ], $
                                 dim_size[ var_inq[k].dim[3] ], dim_size[ var_inq[k].dim[4] ]  )
         else: theData = fltarr( dim_size[ var_inq[k].dim[0] ], dim_size[ var_inq[k].dim[1] ], $
                                 dim_size[ var_inq[k].dim[2] ], dim_size[ var_inq[k].dim[3] ]  )
         endcase

      end

      else: begin

            print, 'ERROR: read_netCDF can only handle 4 dimensions for arrays'
            print, 'Aborting...'
            NCDF_CONTROL, fid, /ABORT
            status = BAD_FILE_DATA
            return
            
      end

  endcase

  var_inq[k].ptr = PTR_NEW( theData )

endfor

if (debug_mode gt 0) then begin

  print, ' '
  nvar = n_elements( var_inq )
  print, 'Indx Lvl -- 0  1 ID 2  3--< 0  1 CT 2  3 >  NAME'
 for jj=0,nvar-1 do print, $
     jj, var_inq[jj].nest_level, var_inq[jj].nest_id[0:3], var_inq[jj].nest_cnt[0:3], $
     var_inq[jj].name, form="(10I4,'   ',A)"
 stop, 'Check out var_inq and dim_name, dim_size...'

endif
;
; Define structures based on var and dim definitions from netCDF file
;  using anonymous structure name with CREATE_STRUCT()
;
; Start with largest nest level and work down to zero level
; Store higher level structures as PTR (in var_inq[XX].str_ptr)
;
; Search backwards in variables for structure definitions
; Assume structure variables are grouped together
;
for nn=max_level,0,-1 do begin

  for k=0, finq.nvars-1 do begin
;
;  Check if new structure found (same nest level as "nn" and cnt = 0)
;  If new, then ss = CREATE_STRUCT( tag, value ) for first parameter and
;    then ss = CREATE_STRUCT( ss, tag, value ) for other parameters
;
    if (k eq 0) then firstzero = var_inq[k].nest_cnt[nn] eq 0 $
    else             firstzero = (var_inq[k].nest_cnt[nn] eq 0) and $
     ( (var_inq[k-1].nest_cnt[nn] ne 0) or (var_inq[k-1].nest_id[nn] ne var_inq[k].nest_id[nn]) )

    if (var_inq[k].nest_level ge nn) and (firstzero) then begin

     if (nn lt var_inq[k].nest_level) then begin

	ss = CREATE_STRUCT( var_inq[k].nest_name[nn+1], *(var_inq[k].str_ptr[nn+1]) )

     endif else begin

        ss = CREATE_STRUCT( var_inq[k].var_name, *(var_inq[k].ptr) )

     endelse

     k1 = k
			
     for kk=k+1, finq.nvars-1 do begin

	k2 = kk
	if ( var_inq[k2].nest_level ge nn ) and $
           ( var_inq[k2].nest_id[nn] eq var_inq[k].nest_id[nn] ) and $
           ( var_inq[k2].nest_cnt[nn] eq (var_inq[k1].nest_cnt[nn] + 1) ) then begin

	   if (nn lt var_inq[kk].nest_level) then begin

            ss = CREATE_STRUCT( ss, var_inq[kk].nest_name[nn+1], *(var_inq[kk].str_ptr[nn+1]) )

           endif else begin

            ss = CREATE_STRUCT( ss, var_inq[kk].var_name, *(var_inq[kk].ptr) )

           endelse

           k1 = k2

        endif

     endfor
;
;    Store new structure as PTR
;    If BASE structure, then replicate for all data reading later
;
     var_inq[k].str_ptr[nn] = PTR_NEW( ss )

     if (nn eq 0) then begin

	data = replicate( ss, max_dim )

     endif

     if (debug_mode gt 0) then begin

      if (nn gt 0) then print, k, nn, '  Structure defined for ', var_inq[k].nest_name[nn] $
      else print, k, nn, '  Base Structure defined as ' 
      help, ss, /struct

     endif

   endif

 endfor

endfor

if (debug_mode gt 0) then begin

  print, ' '
  print, '"data" array size is ', strtrim(max_dim,2)
  stop, 'Check out structure definitions in data...'

endif

;
; Once structures are defined, then read the netCDF variables into "data"
;
; 8.	NCDF_VARGET: Read the data from the variables.
;
for k=0, finq.nvars-1 do begin

  case var_inq[k].nest_level of

    0:  begin

        NCDF_VARGET, fid, k, value
	if ( var_inq[k].type eq 7 ) then $
         data.(var_inq[k].nest_cnt[0]) = string( value ) $
        else data.(var_inq[k].nest_cnt[0]) = value

        end

    1:  begin

        NCDF_VARGET, fid, k, value
        if ( var_inq[k].type eq 7 ) then $
         data.(var_inq[k].nest_cnt[0]).(var_inq[k].nest_cnt[1]) = string( value ) $
        else data.(var_inq[k].nest_cnt[0]).(var_inq[k].nest_cnt[1]) = value

        end

    2:  begin

        NCDF_VARGET, fid, k, value
        if ( var_inq[k].type eq 7 ) then $
          data.(var_inq[k].nest_cnt[0]).(var_inq[k].nest_cnt[1]).(var_inq[k].nest_cnt[2]) = string( value ) $
	else data.(var_inq[k].nest_cnt[0]).(var_inq[k].nest_cnt[1]).(var_inq[k].nest_cnt[2]) = value

       end

    3:  begin
 
        NCDF_VARGET, fid, k, value
        if ( var_inq[k].type eq 7 ) then $
         data.(var_inq[k].nest_cnt[0]).(var_inq[k].nest_cnt[1]).(var_inq[k].nest_cnt[2]).(var_inq[k].nest_cnt[3]) = string( value ) $
	else data.(var_inq[k].nest_cnt[0]).(var_inq[k].nest_cnt[1]).(var_inq[k].nest_cnt[2]).(var_inq[k].nest_cnt[3]) = value

        end

  else: begin

         print, 'ERROR: read_netCDF can only process 4 nested structures'
         print, '       data is lost for ', var_inq[k].name

	end

  endcase

endfor
;
;  Now define "attributes" as string array and read attributes from the netCDF file
;
;   5.	NCDF_ATTINQ: Optionally, retrieve the types and lengths of attributes.
;   6.	NCDF_ATTNAME: Optionally, retrieve attribute names.
;   7.	NCDF_ATTGET: Optionally, retrieve the attributes.
;
;   LIMITATION: limit attributes with more than 1 parameter are compressed into single string
;
CR = string( [ 13B ] )
num_att = 0L
;
; finq.ngatts	= number of GLOBAL attributes from NCDF_INQUIRE earlier
;
if (finq.ngatts gt 0) then num_att = finq.ngatts + 1
for k=0, finq.nvars-1 do if (var_inq[k].natts gt 0) then num_att = num_att + var_inq[k].natts + 1

if ( num_att gt 0 ) then begin

  attributes = strarr( num_att )
  acnt = 0L
;
; Do global variables first
;
  if ( finq.ngatts gt 0) then begin

    attributes[acnt] = 'GLOBAL:' ;	+ CR
    acnt = acnt + 1

    for jj=0,finq.ngatts-1 do begin

     att_name = NCDF_ATTNAME( fid, /GLOBAL, jj )
     NCDF_ATTGET, fid, /GLOBAL, att_name, att_value
     att_str = string( att_value )
     n_str = n_elements(att_str)

     if (n_str gt 1) then begin

      new_str = ''
      for ii=0,n_str-1 do new_str = new_str + ' ' + strtrim(att_str[ii],2)
      att_str = new_str

     endif

     attributes[acnt] = '    ' + att_name + ' = ' + att_str ; + CR
     acnt = acnt + 1

    endfor

  endif
	
  for k=0, finq.nvars-1 do begin

   if (var_inq[k].natts gt 0) then begin

    attributes[acnt] = var_inq[k].name + ':' ;  + CR
    acnt = acnt + 1

    for jj=0,var_inq[k].natts-1 do begin

     att_name = NCDF_ATTNAME( fid, k, jj )
     NCDF_ATTGET, fid, k, att_name, att_value
     att_str = string( att_value )
     n_str = n_elements(att_str)

     if (n_str gt 1) then begin

       new_str = ''
       for ii=0,n_str-1 do new_str = new_str + ' ' + strtrim(att_str[ii],2)
       att_str = new_str

     endif

     attributes[acnt] = '    ' + att_name + ' = ' + att_str ; + CR
     acnt = acnt + 1

    endfor

   endif

  endfor

endif else begin
  attributes = "NONE"
endelse
;
; Close the netCDF file
; 9. NCDF_CLOSE: Close the file.
;
NCDF_CLOSE, fid
;
; Free up Pointers before exiting
;
for k=0, finq.nvars-1 do begin

   if PTR_VALID( var_inq[k].ptr ) then PTR_FREE, var_inq[k].ptr
   for jj=0,5 do if PTR_VALID( var_inq[k].str_ptr[jj] ) then PTR_FREE, var_inq[k].str_ptr[jj]

endfor

return
end
