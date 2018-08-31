;
; read all MLS cat files from /atmos/aura6/data/MLS_data/Datfiles_SOSST/
;
; profile number is now an integer ranging from 1 to nprofiles each day
; orbit number is now date. 
; use profile number with orbit number to uniquely identify a profile.
;
;produce satellite position/time data file
;
output_path='/Users/harvey/Brakebusch/'		;C:\Users\jope9587\Desktop\Sathist\'

;;;;; available data (as of 20130227) - IMPORTANT: UPDATE THOSE LINES BEFORE RUNNING !!!!!

;live instruments (or instruments that need a last update)

MLS_period = [20040808,20180704] ;uses original MLS he5 O3 files [last update 20180704]

mls_obs_path='/atmos/aura6/data/MLS_data/Datfiles_SOSST/'	;Z:\Brakebusch\MLS3\original he5\' ;this path contains cat !!!


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
ncid = ncdf_create(output_path+'sathist_01_EOSMLS_v4_'+start_date+'-'+end_date+'_c'+today()+'.nc',/clobber)
  ncdf_attput, ncid, 'Description', 'Satellite profile positions and times from '+start_date+' through '+end_date, /global
  ncdf_attput, ncid, 'Author', 'V. Lynn Harvey, file was created on '+today()+' using satellite_locations_jmp_mls_c20130529_update_v4.pro', /global
  ncdf_attput, ncid, 'Comment', '(1) EOS MLS O3 he5 files version v04-20 were used (leap seconds not corrected!), '+str_mls_period, /global
  ncdf_attput, ncid, 'prof_num', '(1) MLS: Profile number ranging from 1 to nprofiles each day', /global
  tid = ncdf_dimdef(ncid, 'profs', /unlimited)
  
  vid = ncdf_vardef(ncid, 'date', tid, /LONG)
  vid = ncdf_vardef(ncid, 'time', tid, /LONG)
  vid = ncdf_vardef(ncid, 'lat', tid, /FLOAT)
  vid = ncdf_vardef(ncid, 'lon', tid, /FLOAT)
  vid = ncdf_vardef(ncid, 'orbit_num', tid, /LONG)
  vid = ncdf_vardef(ncid, 'prof_num', tid, /LONG)
  vid = ncdf_vardef(ncid, 'instr_num', tid, /LONG)
  
; vid = ncdf_vardef(ncid, 'doy', tid, /LONG)
; vid = ncdf_vardef(ncid, 'local_time', tid, /FLOAT)
; vid = ncdf_vardef(ncid, 'occ_type', tid, /SHORT)
; vid = ncdf_vardef(ncid, 'instr_sza', tid, /FLOAT)
; vid = ncdf_vardef(ncid, 'julian', tid, /DOUBLE)
  
  ncdf_attput, ncid, 'date', 'long_name', 'date[yyyymmdd]' & ncdf_attput, ncid, 'date', 'units', 'yyyymmdd'
  ncdf_attput, ncid, 'time', 'long_name', 'time of day' & ncdf_attput, ncid, 'time', 'units', 's'
  ncdf_attput, ncid, 'lat', 'long_name', 'latitude' & ncdf_attput, ncid, 'lat', 'units', 'degrees'
  ncdf_attput, ncid, 'lon', 'long_name', 'longitude' & ncdf_attput, ncid, 'lon', 'units', 'degrees'
  ncdf_attput, ncid, 'orbit_num', 'long_name', 'date[yyyymmdd]' & ncdf_attput, ncid, 'orbit_num', 'units', 'yyyymmdd'
  ncdf_attput, ncid, 'prof_num', 'long_name', 'profile number, ... see ''prof_num'' comment' & ncdf_attput, ncid, 'prof_num', 'units', '1'
  ncdf_attput, ncid, 'instr_num', 'long_name', 'MLS=1, ACE-FTS=2, HIRDLS=3, ... see global comment' & ncdf_attput, ncid, 'instr_num', 'units', '1'
  
; ncdf_attput, ncid, 'doy', 'long_name', 'year, day of year' & ncdf_attput, ncid, 'doy', 'units', 'yyyyddd'
; ncdf_attput, ncid, 'local_time', 'long_name', 'local solar time' & ncdf_attput, ncid, 'local_time', 'units', 'hours'
; ncdf_attput, ncid, 'occ_type', 'long_name', 'type of occultation' & ncdf_attput, ncid, 'occ_type', 'units', '1 = sunrise, -1 = sunset, 0 = N/A'
; ncdf_attput, ncid, 'instr_sza', 'long_name', 'solar zenith angle' & ncdf_attput, ncid, 'instr_sza', 'units', 'degrees'
; ncdf_attput, ncid, 'julian', 'long_name', 'julian date' & ncdf_attput, ncid, 'julian', 'units', 'days since '+strcompress(date_reference,/r)
  
  ncdf_control, ncid, /endef

  ;sequentially read data and write to NetCDF
  offset=0l
  uars_mls_old_file=''

  close,2
  openw,2,output_path+'satellite_locations_per_day_c'+today()+'.txt'
  printf,2,'date  nprofs_orig nprofs_incl MLS_orig MLS_incl '

  for iday=0l,ndays-1l do begin
    print, iday+1, ndays
  
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
        
    ;read MLS
    file=file_search(mls_obs_path+'cat_mls_v4.2_'+str_year+str_date_mmdd+'.sav',count=nfiles)
    print,file
    if nfiles gt 1 then file=(file[sort(file)])[nfiles-1] ;pick highest version
    if nfiles ne 0 then begin
;      data=readl2gp_std(file[0])
;      data=read_mls_v3(file[0])
      restore,file
      mls_nprofs = n_elements(fdoy)	; data.ntimes
      mls_latitudes = latitude		; data.latitude
      mls_longitudes = longitude	; data.longitude
      mls_time = time			; data.time
;     mls_ltime = 0.*time		; data.localsolartime
;     mls_sza = 0.*time			; data.solarzenithangle
;     mls_oga = 0.*time			; (data.orbitgeodeticangle)*10
      
;      ;read and apply time correction
;      hdfid = h5f_open(file[0])
;        groupid = h5g_open(hdfid, '/HDFEOS/ADDITIONAL/FILE_ATTRIBUTES')
;          attributeid = h5a_open_name(groupid,'TAI93At0zOfGranule')
;          tai93atozofgranule = (h5a_read(attributeid))[0]
;          h5a_close, attributeid
;        h5g_close,groupid
;      h5f_close, hdfid
;      mls_time = mls_time - tai93atozofgranule
    endif else mls_nprofs = 0
   
    
    ;concatenate
    instrument=0 & latitudes=0 & longitudes=0 & daytime=0 & ltime=0
    orbit=0l & event=0 & occ_type=0 & sza=0 & julian=0
    nprofs = long(mls_nprofs*1l)
    if nprofs eq 0 then goto, nodata
    if mls_nprofs ne 0 then begin
      instrument=[instrument,replicate(1,mls_nprofs)]
      latitudes=[latitudes,mls_latitudes]
      longitudes=[longitudes,mls_longitudes]
      mls_sid=(mls_time/24.)*86400.	;-seconds_since_19930101)
      daytime=[daytime,mls_sid] ;make it time of day in seconds, i.e. starts each day at 0 and goes up to 86400
;     ltime=[ltime,mls_ltime]
      orbit=[orbit,date]				;replicate(-1l,mls_nprofs)]
      event=[event,long(1L+indgen(mls_nprofs))]		;mls_oga]
;     occ_type=[occ_type,replicate(0,mls_nprofs)]
;     sza=[sza,mls_sza]
;     mls_h=floor(mls_sid/3600.)
;     mls_min=floor((mls_sid/3600.-mls_h)*60.)
;     mls_s=round(mls_sid-mls_h*3600.-mls_min*60.)
;     julian=[julian,julday(month,day,year,mls_h,mls_min,mls_s)-julian_reference]
    endif
;
; remove first profile artificially generated by starting to concatenate without declaring arrays
;
     instrument=instrument[1:*]
     latitudes=latitudes[1:*]
     longitudes=longitudes[1:*]
     daytime=daytime[1:*]
;    ltime=ltime[1:*]
     orbit=orbit[1:*]
     event=event[1:*]
;    occ_type=occ_type[1:*]
;    sza=sza[1:*]
;    julian=julian[1:*]
        
    ;data preparation
    yyyymmdd = replicate(current_date,nprofs)
    yyyyddd = replicate(long(year*1000l+current_doy),nprofs)
    longitudes = (longitudes lt 0)*(360+longitudes) + (longitudes ge 0)*longitudes		; add 360 to any negatives
    index = where(longitudes ge 359.99, nindex) & if nindex ne 0 then longitudes[index]=0.
    isort=sort(daytime)
    instrument=instrument[isort]
    latitudes=latitudes[isort]
    longitudes=longitudes[isort]
    daytime=daytime[isort]
;   ltime=ltime[isort]
    orbit=orbit[isort]
    event=event[isort]
;   occ_type=occ_type[isort]
    yyyymmdd=yyyymmdd[isort]
    yyyyddd=yyyyddd[isort]
;   sza=sza[isort]
;   julian=julian[isort]

    ;intercept error codes (-999.99) and values out of range   
    nprofs_old = nprofs
    good = where(yyyymmdd ge start_date and yyyymmdd le end_date and $
                 daytime ge 0 and daytime le 86400 and $
                 latitudes ge -90 and latitudes le 90 and $
                 longitudes ge 0 and longitudes lt 360 and $
                 yyyyddd ge start_doy and yyyyddd le end_doy and $
                 ltime ge 0 and ltime le 24.5 and $
                 event ne -999.99,ngood,complement=bad,ncomplement=nbad)
    if ngood eq 0 then goto, nodata else begin
      yyyymmdd = yyyymmdd[good]
      daytime = daytime[good]
      latitudes = latitudes[good]
      longitudes = longitudes[good]
      yyyyddd = yyyyddd[good]
;     ltime = ltime[good]
      instrument = instrument[good]
      orbit=orbit[good]
      event=event[good]
;     occ_type=occ_type[good]
;     sza=sza[good]
;     julian=julian[good]
    endelse
    nprofs=long(ngood)
    
    mls_incl_nprofs = total(instrument eq 1,/int)

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
    
;   ncdf_varput, ncid, 'doy', yyyyddd, offset=offset, count=nprofs
;   ncdf_varput, ncid, 'local_time', ltime, offset=offset, count=nprofs
;   ncdf_varput, ncid, 'occ_type', occ_type, offset=offset, count=nprofs
;   ncdf_varput, ncid, 'instr_sza', sza, offset=offset, count=nprofs
;   ncdf_varput, ncid, 'julian', julian, offset=offset, count=nprofs
    
    ;increase offset for next day
    offset+=nprofs
    
    ;jump here if no data on this day
    nodata: 
  endfor

ncdf_close, ncid
close,2
print,'...done.'
end
