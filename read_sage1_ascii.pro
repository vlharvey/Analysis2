obs_path = 'D:\SOSST\SAGE1\original (HDF)\ascii\'

yymm = ['7902','7903','7904','7905','7906','7907','7908','7909','7910','7911','7912',$
        '8001','8002','8003','8004','8005','8006','8007','8008','8009','8010','8011','8012',$
        '8101','8102','8103','8104','8105','8106','8107','8108','8109','8110','8111']
        
datafields = ['5','7','9','11','13','19']         

ncid = ncdf_create(obs_path+'SAGE1_catalogfile_c'+today()+'.nc',/clobber)
  ncdf_attput, ncid, 'Description', 'SAGE 1 positions and times', /global
  ncdf_attput, ncid, 'Author', 'Matthias Brakebusch, file was created on '+today()+' using read_sage1_ascii.pro', /global
  ncdf_attput, ncid, 'Comment', 'original HDF files were used, only date, time, lat, lon, occultation type, and fdoy were extracted.', /global
  
  tid = ncdf_dimdef(ncid, 'profs', /unlimited)
  
  vid = ncdf_vardef(ncid, 'date', tid, /LONG)
  vid = ncdf_vardef(ncid, 'time', tid, /LONG)
  vid = ncdf_vardef(ncid, 'lat', tid, /FLOAT)
  vid = ncdf_vardef(ncid, 'lon', tid, /FLOAT)
  vid = ncdf_vardef(ncid, 'occ_type', tid, /LONG)
  vid = ncdf_vardef(ncid, 'fdoy', tid, /FLOAT)
  
  ncdf_attput, ncid, 'date', 'long_name', 'date[yyyymmdd]' & ncdf_attput, ncid, 'date', 'units', 'yyyymmdd'
  ncdf_attput, ncid, 'time', 'long_name', 'time of day' & ncdf_attput, ncid, 'time', 'units', 's'
  ncdf_attput, ncid, 'lat', 'long_name', 'latitude' & ncdf_attput, ncid, 'lat', 'units', 'degrees'
  ncdf_attput, ncid, 'lon', 'long_name', 'longitude' & ncdf_attput, ncid, 'lon', 'units', 'degrees'
  
  ncdf_attput, ncid, 'occ_type', 'long_name', 'type of occultation' & ncdf_attput, ncid, 'occ_type', 'units', '1 = sunrise, -1 = sunset, 0 = N/A'
  ncdf_attput, ncid, 'fdoy', 'long_name', 'floating point day of year' & ncdf_attput, ncid, 'fdoy', 'units', 'd'
  
  ncdf_control, ncid, /endef
  
  offset = 0
  for imonth=0,n_elements(yymm)-1 do begin
    nprofs = -1
    for idata=0,n_elements(datafields)-1 do begin
      filename = obs_path+yymm[imonth]+'Data-Set-'+datafields[idata]+'.txt'
      raw_data = (read_ascii(filename)).field001
      
      case datafields[idata] of
        '5': date = 19000000l+long(raw_data)
        '7': begin
               h = floor(raw_data/10000.)
               min = floor((raw_data - h*10000.)/100.)
               sec = round(raw_data-h*10000.-min*100.)
               sid = h*3600.+min*60.+sec 
             end 
        '9': lat = raw_data
        '11': lon = raw_data
        '13': occ_type = 1*(raw_data eq 0.)-1*(raw_data eq 1.)  ;1 is rise, -1 is set
        '19': fdoy = raw_data
      endcase
      if nprofs eq -1 then nprofs = n_elements(raw_data) $
                      else begin
                        if nprofs ne n_elements(raw_data) then stop
                      endelse
    endfor ;end of idata loop
    
    ;write NetCDF file
    ncdf_varput, ncid, 'date', date, offset=offset, count=nprofs
    ncdf_varput, ncid, 'time', sid, offset=offset, count=nprofs
    ncdf_varput, ncid, 'lat', lat, offset=offset, count=nprofs
    ncdf_varput, ncid, 'lon', lon, offset=offset, count=nprofs
    ncdf_varput, ncid, 'occ_type', occ_type, offset=offset, count=nprofs
    ncdf_varput, ncid, 'fdoy', fdoy, offset=offset, count=nprofs
    
    ;increase offset for next day
    offset += nprofs
    
  endfor ;end of yymm loop

ncdf_close,ncid


print, '...done.'
end