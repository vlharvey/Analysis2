
; assign instrument number 40 for SAGE III/ISS
;
;produce satellite position/time data file

; This file was modified by Josh Pettit on 3/9/2018. My goal was to output each instrument separately as its own netcdf sathist file

output_path='/Users/harvey/Brakebusch/'			;W:\Desktop\SatHist\'

;;;;; IMPORTANT: UPDATE THOSE LINES BEFORE RUNNING !!!!!

;live instruments (or instruments that need a last update)

ACE_FTS_period = [20170607,20180430] 		; SAGE III/ISS data range

ace_obs_path='/atmos/aura3/data/SAGE_3ISS_data/Datfiles_Sathist/'
period_array = [ace_fts_period]
start_date = strcompress(min(period_array),/r)
end_date   = strcompress(max(period_array),/r)

ndays=date2julday(end_date) - date2julday(start_date) +1

start_julday = date2julday(start_date)
start_doy = strmid(start_date,0,4)+strcompress(string(date2doy(start_date),format='(i3.3)'),/r)
end_doy = strmid(end_date,0,4)+strcompress(string(date2doy(end_date),format='(i3.3)'),/r)

date_reference = 19700101
julian_reference = date2julday(date_reference)
if julian_reference gt start_julday then stop

str_ace_fts_period=strcompress(string(ace_fts_period,format='(i8)'),/r) & str_ace_fts_period = str_ace_fts_period[0]+' - '+str_ace_fts_period[1]

;open NetCDF file for sequential writing
  ncid = ncdf_create(output_path+'sathist_40_SAGE3ISS_v5_'+start_date+'-'+end_date+'_c'+today()+'.nc',/clobber)
  ncdf_attput, ncid, 'Description', 'Satellite profile positions and times from '+start_date+' through '+end_date, /global
  ncdf_attput, ncid, 'Author', 'V. Lynn Harvey, file was created on '+today()+' using satellite_locations_vlh_sage3iss_c20180309_update.pro', /global
  ncdf_attput, ncid, 'Comment', '(40) SAGE3ISS geolocations in snoe_mission_profile_locs.sav from Scott Bailey were used, '+str_ace_fts_period+'; ', /global
  ncdf_attput, ncid, 'prof_num', '(40) SAGE3ISS: prof_num ranges from 1 to nprof each day ', /global
  ncdf_attput, ncid, 'orbit_num', '(40) SAGE3ISS: Orbit Number is the original orbit number and uniquely identifies profiles ', /global

  tid = ncdf_dimdef(ncid, 'profs', /unlimited)

  vid = ncdf_vardef(ncid, 'date', tid, /LONG)
  vid = ncdf_vardef(ncid, 'time', tid, /LONG)
  vid = ncdf_vardef(ncid, 'lat', tid, /FLOAT)
  vid = ncdf_vardef(ncid, 'lon', tid, /FLOAT)
  vid = ncdf_vardef(ncid, 'orbit_num', tid, /LONG)
  vid = ncdf_vardef(ncid, 'prof_num', tid, /LONG)
  vid = ncdf_vardef(ncid, 'instr_num', tid, /LONG)

  ncdf_attput, ncid, 'date', 'long_name', 'date[yyyymmdd]' & ncdf_attput, ncid, 'date', 'units', 'yyyymmdd'
  ncdf_attput, ncid, 'time', 'long_name', 'time of day' & ncdf_attput, ncid, 'time', 'units', 's'
  ncdf_attput, ncid, 'lat', 'long_name', 'latitude' & ncdf_attput, ncid, 'lat', 'units', 'degrees'
  ncdf_attput, ncid, 'lon', 'long_name', 'longitude' & ncdf_attput, ncid, 'lon', 'units', 'degrees'
  ncdf_attput, ncid, 'orbit_num', 'long_name', 'orbit number' & ncdf_attput, ncid, 'orbit_num', 'units', '1'
  ncdf_attput, ncid, 'prof_num', 'long_name', 'profile number, ... see ''prof_num'' comment' & ncdf_attput, ncid, 'prof_num', 'units', '1'
  ncdf_attput, ncid, 'instr_num', 'long_name', 'MLS=1, ACE-FTS=2, HIRDLS=3, ... see global comment' & ncdf_attput, ncid, 'instr_num', 'units', '1'

  ncdf_control, ncid, /endef

  ;sequentially read data and write to NetCDF
  offset=0l
  ace_old_file=''

  close,2
  openw,2,output_path+'satellite_locations_per_day_c'+today()+'.txt'
  printf,2,'date  nprofs_orig nprofs_incl SAGE3ISS'
;
; loop over days (this is currently unnecessary but done for historical purposes)
;
  for iday=0l,ndays-1l do begin
;   print, iday+1, ndays

    ;produce date variables
    current_julday = start_julday + iday
    caldat, current_julday, month, day, year
    current_doy = current_julday-julday(1,1,year)+1
    current_date = long(year*10000l+month*100l+day)
    str_year = strcompress(string(year),/r)
    str_date_code = date2datecode(julday2date(start_julday + iday))
;   seconds_since_19930101=86400*(current_julday-julday(1,1,1993)) ;for MLS
    ;print, str_date_code

    ;read entire dataset
    file=file_search(ace_obs_path+'SAGE3ISS_v5_sathist_2017-2018.sav',count=nfiles)
    if iday eq 0L then begin
;
; DATE            LONG      = Array[6846441]
; LATITUDE        FLOAT     = Array[6846441]
; LONGITUDE       FLOAT     = Array[6846441]
; ORBIT           LONG      = Array[6846441]
; PROF_NUM        LONG      = Array[6846441]
; TIME            DOUBLE    = Array[6846441]
;
       restore,file
profile_number=prof_num
orbit_number=orbit
;
; ace_orig_occ_num is profile number and ultimately event
;
ace_orig_occ_num=PROFILE_NUMBER
;
; set ace_orig_date to date
;
ace_orig_date=date
;
; set longitude, latitude, time, orbit number
;
ace_orig_longitudes=longitude
ace_orig_latitudes=latitude
ace_orig_time=(time/24.d)*86400.d
ace_orig_orbit_num=ORBIT_NUMBER
      endif	; if first day

      index=where(ace_orig_date eq current_date,nindex)
      if nindex ne 0 then begin
        ace_nprofs = nindex
        ace_longitudes = ace_orig_longitudes[index]
        ace_latitudes = ace_orig_latitudes[index]
        ace_time = ace_orig_time[index]
        ace_occ_num = ace_orig_occ_num[index]
        ace_orbit_num = ace_orig_orbit_num[index]
      endif else ace_nprofs = 0
    sat_num=0 ; needed for log file

    ;concatenate
    instrument=0 & latitudes=0 & longitudes=0 & daytime=0 & ltime=0
    orbit=0l & event=0L & occ_type=0 & sza=0 & julian=0
    nprofs = long(ace_nprofs)
    if nprofs eq 0 then goto, nodata
    if ace_nprofs ne 0 then begin
      instrument=[instrument,replicate(40,ace_nprofs)]			; 40=SAGE3ISS
      latitudes=[latitudes,ace_latitudes]
      longitudes=[longitudes,ace_longitudes]
      ace_sid=round(ace_time)
      daytime=[daytime,ace_sid] ;make it time of day in seconds, i.e. starts each day at 0 and goes up to 86400
      orbit=[orbit,ace_orbit_num]
      event=[event,ace_occ_num]
    endif

    instrument=instrument[1:*]
    latitudes=latitudes[1:*]
    longitudes=longitudes[1:*]
    daytime=daytime[1:*]
    orbit=orbit[1:*]
    event=event[1:*]

    ;data preparation
    yyyymmdd = replicate(current_date,nprofs)
    yyyyddd = replicate(long(year*1000l+current_doy),nprofs)
    longitudes = (longitudes lt 0)*(360+longitudes) + (longitudes ge 0)*longitudes
;   index = where(longitudes ge 359.99, nindex) & if nindex ne 0 then longitudes[index]=0.
    isort=sort(daytime)
    instrument=instrument[isort]
    latitudes=latitudes[isort]
    longitudes=longitudes[isort]
    daytime=daytime[isort]
    orbit=orbit[isort]
    event=event[isort]
    yyyymmdd=yyyymmdd[isort]
    yyyyddd=yyyyddd[isort]

    ;intercept error codes (-999.99) and values out of range
    nprofs_old = nprofs
    good = where(yyyymmdd ge start_date and yyyymmdd le end_date and $
                 daytime ge 0 and daytime lt 86400 and $
                 latitudes ge -90 and latitudes le 90 and $
                 longitudes ge 0 and longitudes lt 360 and $
                 yyyyddd ge start_doy and yyyyddd le end_doy and $
                 ltime ge 0 and ltime le 24.5 and $
                 event ne -999.99,ngood,complement=bad,ncomplement=nbad)
    if ngood eq 0 then goto, nodata else begin
      yyyymmdd = yyyymmdd[good]
      daytime = daytime[good]
print,current_date,max(daytime)
if max(daytime) gt 86399 then stop
      latitudes = latitudes[good]
      longitudes = longitudes[good]
      yyyyddd = yyyyddd[good]
      instrument = instrument[good]
      orbit=orbit[good]
      event=event[good]
    endelse
    nprofs=long(ngood)

    ace_incl_nprofs = total(instrument eq 40,/int)
 
    ;write satellite_locations_per_day.txt

    printf,2,str_date_code+$
             '  '+strcompress(string(ace_nprofs),/r)+$
             '  '+strcompress(string(ace_incl_nprofs),/r)

;    print,str_date_code+$
;         '  '+strcompress(string(mls_nprofs),/r)+$
;         '  '+strcompress(string(ace_nprofs),/r)+$
;         '  '+strcompress(string(hirdls_nprofs),/r)+$
;         '  '+strcompress(string(saber_nprofs),/r)+' done.' ;MLS, ACE, HIRDLS, SABER

    ;write NetCDF file
    ncdf_varput, ncid, 'date', yyyymmdd, offset=offset, count=nprofs
    ncdf_varput, ncid, 'time', daytime, offset=offset, count=nprofs
    ncdf_varput, ncid, 'lat', latitudes, offset=offset, count=nprofs
    ncdf_varput, ncid, 'lon', longitudes, offset=offset, count=nprofs
    ncdf_varput, ncid, 'orbit_num', orbit, offset=offset, count=nprofs
    ncdf_varput, ncid, 'prof_num', event, offset=offset, count=nprofs
    ncdf_varput, ncid, 'instr_num', instrument, offset=offset, count=nprofs

    ;increase offset for next day
    offset+=nprofs

    ;jump here if no data on this day
    nodata:
  endfor

ncdf_close, ncid
close,2
print,'...done.'
end
