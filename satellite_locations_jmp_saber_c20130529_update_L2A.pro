;
; profile number is now an integer ranging from 1 to nprofiles each day
; orbit number is now date. 
; use profile number with orbit number to uniquely identify a profile.
;
;produce satellite position/time data file
;
output_path='/Users/harvey/Brakebusch/'		;C:\Users\jope9587\Desktop\Sathist\'

;;;;; IMPORTANT: UPDATE Date Range LINES BEFORE RUNNING !!!!!

;live instruments (or instruments that need a last update)

MLS_period = [20020125,20180703] ;based on daily L2A files

mls_obs_path='/atmos/harvey/SABER_data/Datfiles_Sathist/'	;Z:\Brakebusch\MLS3\original he5\' ;this path contains cat !!!


;WACCM grid
waccm_lon = indgen(144)/144.*360. & waccm_lon_step = 2.5
waccm_lat = indgen(96)/95.*180.-90. & waccm_lat_step = 1.89475

period_array = [mls_period] 
start_date = strcompress(min(period_array),/r)
end_date   = strcompress(max(period_array),/r)

ndays=date2julday(end_date) - date2julday(start_date) +1
      
start_julday = date2julday(start_date)
start_doy = strmid(start_date,0,4)+strcompress(string(date2doy(start_date),format='(i3.3)'),/r)
end_doy = strmid(end_date,0,4)+strcompress(string(date2doy(end_date),format='(i3.3)'),/r)

date_reference = 19700101
julian_reference = date2julday(date_reference)
if julian_reference gt start_julday then stop

str_mls_period=strcompress(string(mls_period,format='(i8)'),/r) & str_mls_period = str_mls_period[0]+' - '+str_mls_period[1]

;open NetCDF file for sequential writing
ncid = ncdf_create(output_path+'sathist_04_SABER_v2_'+start_date+'-'+end_date+'_c'+today()+'.nc',/clobber)
ncdf_attput, ncid, 'Description', 'Satellite profile positions and times from '+start_date+' through '+end_date, /global
ncdf_attput, ncid, 'Author', 'V. Lynn Harvey, file was created on '+today()+' using satellite_locations_jmp_mls_c20130326_update_L2A.pro', /global
ncdf_attput, ncid, 'Comment', '(4) SABER daily L2A files version v2 were used, '+str_mls_period, /global
ncdf_attput, ncid, 'prof_num', '(4) SABER: Orbit profile number', /global
tid = ncdf_dimdef(ncid, 'profs', /unlimited)
  
  vid = ncdf_vardef(ncid, 'date', tid, /LONG)
  vid = ncdf_vardef(ncid, 'time', tid, /LONG)
  vid = ncdf_vardef(ncid, 'lat', tid, /FLOAT)
  vid = ncdf_vardef(ncid, 'lon', tid, /FLOAT)
  vid = ncdf_vardef(ncid, 'orbit_num', tid, /LONG)
  vid = ncdf_vardef(ncid, 'prof_num', tid, /LONG)
  vid = ncdf_vardef(ncid, 'instr_num', tid, /LONG)
  
  ncdf_attput, ncid, 'date', 'long_name', 'date[yyyymmdd]' & ncdf_attput, ncid, 'date', 'units', 'yyyymmdd'
  ncdf_attput, ncid, 'time', 'long_name', 'time of day at 80 km' & ncdf_attput, ncid, 'time', 'units', 's'
  ncdf_attput, ncid, 'lat', 'long_name', 'latitude at 80 km' & ncdf_attput, ncid, 'lat', 'units', 'degrees'
  ncdf_attput, ncid, 'lon', 'long_name', 'longitude at 80 km' & ncdf_attput, ncid, 'lon', 'units', 'degrees'
  ncdf_attput, ncid, 'orbit_num', 'long_name', 'Orbit number' & ncdf_attput, ncid, 'orbit_num', 'units', '1'
  ncdf_attput, ncid, 'prof_num', 'long_name', 'profile number, ... see ''prof_num'' comment' & ncdf_attput, ncid, 'prof_num', 'units', '1'
  ncdf_attput, ncid, 'instr_num', 'long_name', 'MLS=1, ACE-FTS=2, HIRDLS=3, ... see global comment' & ncdf_attput, ncid, 'instr_num', 'units', '1'
  
  ncdf_control, ncid, /endef

  ;sequentially read data and write to NetCDF
  offset=0l
  uars_mls_old_file=''

  close,2
  openw,2,output_path+'satellite_locations_per_day_c'+today()+'.txt'
  printf,2,'date  nprofs_orig nprofs_incl SABER_orig SABER_incl '
;
;read SABER and retain original geolocation information
;
; DATE            LONG      = Array[6561263]
; LATITUDE        FLOAT     = Array[6561263]
; LONGITUDE       FLOAT     = Array[6561263]
; ORBIT           LONG      = Array[6561263]
; TIME            FLOAT     = Array[6561263]

file=file_search(mls_obs_path+'SABER_L2A_v2_sathist_2002-2018.sav',count=nfiles)
restore,file
print,file

good=where(date ge MLS_period(0) and date le MLS_period(1))
date=date(good)
latitude=latitude(good)
longitude=longitude(good)
orbit=long(orbit(good))
time=time(good)
;
; define prof_num
;
nevent=n_elements(date)
n0=findgen(nevent)
n1=1+findgen(nevent)
index=where(orbit(n1)-orbit(n0) gt 0L)
orbit_numbers=[orbit(index),orbit(-1)]
prof_num=0L*orbit
for i=0L,n_elements(orbit_numbers)-1L do begin
    index=where(orbit eq orbit_numbers(i))
    prof_num(index)=long(1+indgen(n_elements(index)))
endfor

date_orig=date
longitude_orig=longitude
latitude_orig=latitude
time_orig=time
orbit_orig=orbit
profnum_orig=prof_num
instr_num=4L+0*date		; SABER is instrument 4
;
; loop over days to count daily profiles
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
    str_date_mmdd = string(FORMAT='(I2.2)',month)+string(FORMAT='(I2.2)',day)

    seconds_since_19930101=86400*(current_julday-julday(1,1,1993)) ;for MLS
    ;print, str_date_code

    index=where(date_orig eq current_date,nindex)
    print,current_date,nindex
    if nindex ne 0 then begin
      mls_nprofs = nindex
      mls_dates = date_orig(index)
      mls_latitudes = latitude_orig(index)
      mls_longitudes = longitude_orig(index)
      mls_time = time_orig(index)
      mls_orbit = orbit_orig(index)
      mls_profnum = profnum_orig(index)
    endif else mls_nprofs = 0
    
    ;concatenate
    yyyymmdd=0 & instrument=0 & latitudes=0 & longitudes=0 & daytime=0 & orbit=0l & event=0
    nprofs = long(mls_nprofs*1l)
    if nprofs eq 0 then goto, nodata
    if mls_nprofs ne 0 then begin
      yyyymmdd=[yyyymmdd,mls_dates]
      instrument=[instrument,replicate(4,mls_nprofs)]
      latitudes=[latitudes,mls_latitudes]
      longitudes=[longitudes,mls_longitudes]
      mls_sid=(mls_time/24.)*86400.	;-seconds_since_19930101)
      daytime=[daytime,mls_sid] ;make it time of day in seconds, i.e. starts each day at 0 and goes up to 86400
      orbit=[orbit,mls_orbit]				;replicate(-1l,mls_nprofs)]
      event=[event,mls_profnum]
    endif
;
; remove first profile artificially generated by starting to concatenate without declaring arrays
;
     yyyymmdd=yyyymmdd[1:*]
     instrument=instrument[1:*]
     latitudes=latitudes[1:*]
     longitudes=longitudes[1:*]
     daytime=daytime[1:*]
     orbit=orbit[1:*]
     event=event[1:*]
        
    ;data preparation
; data is already sorted!
;   yyyymmdd = replicate(current_date,nprofs)
;   longitudes = (longitudes lt 0)*(360+longitudes) + (longitudes ge 0)*longitudes		; add 360 to any negatives
;   index = where(longitudes ge 359.99, nindex) & if nindex ne 0 then longitudes[index]=0.
;   isort=sort(daytime)
;   yyyymmdd=yyyymmdd[isort]
;   instrument=instrument[isort]
;   latitudes=latitudes[isort]
;   longitudes=longitudes[isort]
;   daytime=daytime[isort]
;   orbit=orbit[isort]
;   event=event[isort]

    ;intercept error codes (-999.99) and values out of range   
    nprofs_old = nprofs
    good = where(yyyymmdd ge start_date and yyyymmdd le end_date and $
                 daytime ge 0 and daytime le 86400 and $
                 latitudes ge -90 and latitudes le 90 and $
                 longitudes ge 0 and longitudes lt 360 and $
                 event ne -999.99,ngood,complement=bad,ncomplement=nbad)
    if ngood eq 0 then goto, nodata else begin
      yyyymmdd = yyyymmdd[good]
      daytime = daytime[good]
      latitudes = latitudes[good]
      longitudes = longitudes[good]
      instrument = instrument[good]
      orbit=orbit[good]
      event=event[good]
    endelse
    nprofs=long(ngood)
    
    mls_incl_nprofs = total(instrument eq 4,/int)

    printf,2,str_date_code+$
             '  '+strcompress(string(mls_nprofs),/r)+$
             '  '+strcompress(string(mls_incl_nprofs),/r)

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
