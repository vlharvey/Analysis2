
;produce satellite position/time data file

; This file was modified by Josh Pettit on 4/6/2015. My goal was to output each instrument separately as its own netcdf sathist file
 
output_path='W:\Desktop\SatHist\'

;;;;; available data (as of 20130227) - IMPORTANT: UPDATE THOSE LINES BEFORE RUNNING !!!!!

;live instruments (or instruments that need a last update)

ACE_FTS_period = [20040110,20120818] ;uses yearly ACE-FTS nc files (produced using ACE_raw2netcdf.pro) [last updated 20130226]

ace_obs_path='Z:\Brakebusch\SOSST\ACE original nc\'


;WACCM grid
waccm_lon = indgen(144)/144.*360. & waccm_lon_step = 2.5
waccm_lat = indgen(96)/95.*180.-90. & waccm_lat_step = 1.89475

period_array = [ace_fts_period] 
start_date = strcompress(min(period_array),/r)
end_date   = strcompress(max(period_array),/r)

;start_date = '20070514'
;end_date = '20130129'


ndays=date2julday(end_date) - date2julday(start_date) +1
      
start_julday = date2julday(start_date)
start_doy = strmid(start_date,0,4)+strcompress(string(date2doy(start_date),format='(i3.3)'),/r)
end_doy = strmid(end_date,0,4)+strcompress(string(date2doy(end_date),format='(i3.3)'),/r)

date_reference = 19700101
julian_reference = date2julday(date_reference)
if julian_reference gt start_julday then stop

;;read Aura time correction file (TAI93At0zOfGranule)
;file='Z:\Brakebusch\MLS\TAIat0z.txt'
;data=read_ascii(file,template=aura_tpl)
;aura_strdatecode = strmid(data.date,0,4)+'d'+strmid(data.date,5,3)
;aura_time_corr = data.TAI93At0zOfGranule


str_ace_fts_period=strcompress(string(ace_fts_period,format='(i8)'),/r) & str_ace_fts_period = str_ace_fts_period[0]+' - '+str_ace_fts_period[1]



;open NetCDF file for sequential writing
ncid = ncdf_create(output_path+'sat_hist_'+start_date+'-'+end_date+'_ace_c'+today()+'.nc',/clobber)
  ncdf_attput, ncid, 'Description', 'Satellite profile positions and times from '+start_date+' through '+end_date, /global
  ncdf_attput, ncid, 'Author', 'Matthias Brakebusch, file was created on '+today()+' using satellite_locations_c20130326.pro', /global
  ncdf_attput, ncid, 'Comment', '(2) ACE-FTS ascii files version v3.0 were used, '+str_ace_fts_period+'; ', /global
                                
  ncdf_attput, ncid, 'prof_num', '(2) ACE-FTS: OccultationNumber; ', /global

  tid = ncdf_dimdef(ncid, 'profs', /unlimited)
  
  vid = ncdf_vardef(ncid, 'date', tid, /LONG)
  vid = ncdf_vardef(ncid, 'time', tid, /LONG)
  vid = ncdf_vardef(ncid, 'lat', tid, /FLOAT)
  vid = ncdf_vardef(ncid, 'lon', tid, /FLOAT)
  vid = ncdf_vardef(ncid, 'orbit_num', tid, /LONG)
  vid = ncdf_vardef(ncid, 'prof_num', tid, /LONG)
  vid = ncdf_vardef(ncid, 'instr_num', tid, /LONG)
  
  vid = ncdf_vardef(ncid, 'doy', tid, /LONG)
  vid = ncdf_vardef(ncid, 'local_time', tid, /FLOAT)
  vid = ncdf_vardef(ncid, 'occ_type', tid, /SHORT)
  vid = ncdf_vardef(ncid, 'instr_sza', tid, /FLOAT)
  vid = ncdf_vardef(ncid, 'julian', tid, /DOUBLE)
  
  ncdf_attput, ncid, 'date', 'long_name', 'date[yyyymmdd]' & ncdf_attput, ncid, 'date', 'units', 'yyyymmdd'
  ncdf_attput, ncid, 'time', 'long_name', 'time of day' & ncdf_attput, ncid, 'time', 'units', 's'
  ncdf_attput, ncid, 'lat', 'long_name', 'latitude' & ncdf_attput, ncid, 'lat', 'units', 'degrees'
  ncdf_attput, ncid, 'lon', 'long_name', 'longitude' & ncdf_attput, ncid, 'lon', 'units', 'degrees'
  ncdf_attput, ncid, 'orbit_num', 'long_name', 'orbit number' & ncdf_attput, ncid, 'orbit_num', 'units', '1'
  ncdf_attput, ncid, 'prof_num', 'long_name', 'profile number, ... see ''prof_num'' comment' & ncdf_attput, ncid, 'prof_num', 'units', '1'
  ncdf_attput, ncid, 'instr_num', 'long_name', 'MLS=1, ACE-FTS=2, HIRDLS=3, ... see global comment' & ncdf_attput, ncid, 'instr_num', 'units', '1'
  
  ncdf_attput, ncid, 'doy', 'long_name', 'year, day of year' & ncdf_attput, ncid, 'doy', 'units', 'yyyyddd'
  ncdf_attput, ncid, 'local_time', 'long_name', 'local solar time' & ncdf_attput, ncid, 'local_time', 'units', 'hours'
  ncdf_attput, ncid, 'occ_type', 'long_name', 'type of occultation' & ncdf_attput, ncid, 'occ_type', 'units', '1 = sunrise, -1 = sunset, 0 = N/A'
  ncdf_attput, ncid, 'instr_sza', 'long_name', 'solar zenith angle' & ncdf_attput, ncid, 'instr_sza', 'units', 'degrees'
  ncdf_attput, ncid, 'julian', 'long_name', 'julian date' & ncdf_attput, ncid, 'julian', 'units', 'days since '+strcompress(date_reference,/r)
  
  ncdf_control, ncid, /endef

  ;sequentially read data and write to NetCDF
  offset=0l
  ace_old_file=''

  openw,2,output_path+'satellite_locations_per_day_c'+today()+'.txt'
  printf,2,'date  nprofs_orig nprofs_incl ACE-FTS_orig ACE-FTS_incl'


  for iday=0l,ndays-1l do begin
    print, iday+1, ndays
  
    ;produce date variables
    current_julday = start_julday + iday
    caldat, current_julday, month, day, year
    current_doy = current_julday-julday(1,1,year)+1
    current_date = long(year*10000l+month*100l+day)
    str_year = strcompress(string(year),/r)
    str_date_code = date2datecode(julday2date(start_julday + iday))
    seconds_since_19930101=86400*(current_julday-julday(1,1,1993)) ;for MLS
    ;print, str_date_code

    ;read ACE-FTS
    file=file_search(ace_obs_path+'ACE-FTS_v3.0_'+str_year+'.nc',count=nfiles)
    if nfiles gt 1 then file=(file[sort(file)])[nfiles-1] ;pick highest version - shouldn't happen with ACE data anyway
    if nfiles ne 0 then begin
      if file ne ace_old_file then begin 
        ncid2 = ncdf_open(file[0])
          ncdf_varget,ncid2,ncdf_varid(ncid2,'longitude'), ace_orig_longitudes
          ncdf_varget,ncid2,ncdf_varid(ncid2,'latitude'), ace_orig_latitudes
          ncdf_varget,ncid2,ncdf_varid(ncid2,'date'), ace_orig_date
          ncdf_varget,ncid2,ncdf_varid(ncid2,'time'), ace_orig_time
          ncdf_varget,ncid2,ncdf_varid(ncid2,'occultation_number'), ace_orig_occ_num
          ncdf_varget,ncid2,ncdf_varid(ncid2,'occ_type'), ace_orig_occ_type
        ncdf_close,ncid2
        ace_old_file=file
        ace_orig_nprofs=n_elements(ace_orig_longitudes)
        print,'starting ACE-FTS '+str_year
      endif
      index=where(ace_orig_date eq current_date,nindex)
      if nindex ne 0 then begin
        ace_nprofs = nindex
        ace_longitudes = ace_orig_longitudes[index]
        ace_latitudes = ace_orig_latitudes[index]
        ace_time = ace_orig_time[index]
        ace_occ_num = ace_orig_occ_num[index] 
        ace_occ_type = ace_orig_occ_type[index]
        ace_occ_type = fix(1*(ace_occ_type eq 114)-1*(ace_occ_type eq 115)) ;change r (114) to 1 and s (115) to -1
      endif else ace_nprofs = 0
    endif else ace_nprofs = 0
    sat_num=0 ; needed for log file
              
    ;concatenate
    instrument=0 & latitudes=0 & longitudes=0 & daytime=0 & ltime=0
    orbit=0l & event=0 & occ_type=0 & sza=0 & julian=0
    nprofs = long(ace_nprofs)
    if nprofs eq 0 then goto, nodata
    if ace_nprofs ne 0 then begin
      instrument=[instrument,replicate(2,ace_nprofs)]
      latitudes=[latitudes,ace_latitudes]
      longitudes=[longitudes,ace_longitudes]
      ace_sid=round(ace_time*3600)
      daytime=[daytime,ace_sid] ;make it time of day in seconds, i.e. starts each day at 0 and goes up to 86400
      ltime=[ltime,utc2ltime(ace_time,ace_longitudes)]
      orbit=[orbit,replicate(-1l,ace_nprofs)]
      event=[event,ace_occ_num]
      occ_type=[occ_type,ace_occ_type]
      sza=[sza,replicate(-1,ace_nprofs)]
      ace_h=floor(ace_sid/3600.)
      ace_min=floor((ace_sid/3600.-ace_h)*60.)
      ace_s=round(ace_sid-ace_h*3600.-ace_min*60.)
      julian=[julian,julday(month,day,year,ace_h,ace_min,ace_s)-julian_reference]
    endif


    instrument=instrument[1:*]
    latitudes=latitudes[1:*]
    longitudes=longitudes[1:*]
    daytime=daytime[1:*]
    ltime=ltime[1:*]
    orbit=orbit[1:*]
    event=event[1:*]
    occ_type=occ_type[1:*]
    sza=sza[1:*]
    julian=julian[1:*]
        
    ;data preparation
    yyyymmdd = replicate(current_date,nprofs)
    yyyyddd = replicate(long(year*1000l+current_doy),nprofs)
    longitudes = (longitudes lt 0)*(360+longitudes) + (longitudes ge 0)*longitudes
    index = where(longitudes ge 359.99, nindex) & if nindex ne 0 then longitudes[index]=0.
    isort=sort(daytime)
    instrument=instrument[isort]
    latitudes=latitudes[isort]
    longitudes=longitudes[isort]
    daytime=daytime[isort]
    ltime=ltime[isort]
    orbit=orbit[isort]
    event=event[isort]
    occ_type=occ_type[isort]
    yyyymmdd=yyyymmdd[isort]
    yyyyddd=yyyyddd[isort]
    sza=sza[isort]
    julian=julian[isort]
        
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
      ltime = ltime[good]
      instrument = instrument[good]
      orbit=orbit[good]
      event=event[good]
      occ_type=occ_type[good]
      sza=sza[good]
      julian=julian[good]
    endelse
    nprofs=long(ngood)
    
    ace_incl_nprofs = total(instrument eq 2,/int)
 
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
    
    ncdf_varput, ncid, 'doy', yyyyddd, offset=offset, count=nprofs
    ncdf_varput, ncid, 'local_time', ltime, offset=offset, count=nprofs
    ncdf_varput, ncid, 'occ_type', occ_type, offset=offset, count=nprofs
    ncdf_varput, ncid, 'instr_sza', sza, offset=offset, count=nprofs
    ncdf_varput, ncid, 'julian', julian, offset=offset, count=nprofs
    
    ;increase offset for next day
    offset+=nprofs
    
    ;jump here if no data on this day
    nodata: 
  endfor

ncdf_close, ncid
close,2
print,'...done.'
end