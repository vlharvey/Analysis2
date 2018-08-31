function read_mls_v3, filename
  ;extract species
  str_species = strmid(filename,strpos(filename,'L2GP-')+5)
  str_species = strmid(str_species,0,strpos(str_species,'_v0'))
  str_datecode = strmid(filename,11,8,/rev)
  str_julday = date2julday(doy2date(strmid(str_datecode,5),strmid(str_datecode,0,4)))
    
  hdfid = h5f_open(filename)
    
    ;read attribute
    groupid = h5g_open(hdfid, '/HDFEOS/ADDITIONAL/FILE_ATTRIBUTES')
      if str_julday lt 2455778 then begin ;after 2011d216 (=20110804, =2455778) there is no TAI93At0zOfGranule anymore
        attributeid = h5a_open_name(groupid,'TAI93At0zOfGranule')
        tai93atozofgranule = h5a_read(attributeid)
        h5a_close, attributeid
      endif else tai93atozofgranule = 0.0
    h5g_close, groupid
    
    ;read general parameters
    general=['nLevels','nTimes']
    for i=0,n_elements(general)-1 do begin
      datasetid = h5d_open(hdfid, '/HDFEOS/SWATHS/'+str_species+'/'+general[i])
      case general[i] of
        'nLevels': nlev = h5d_read(datasetid)
        'nTimes':  ntimes = h5d_read(datasetid)
      endcase
    endfor
    h5d_close, datasetid
    
    ;read Geolocation Fields
    geolocation=['Latitude','Longitude','LocalSolarTime','Pressure','SolarZenithAngle','Time','OrbitGeodeticAngle']
    for i=0,n_elements(geolocation)-1 do begin
      datasetid = h5d_open(hdfid, '/HDFEOS/SWATHS/'+str_species+'/Geolocation Fields/'+geolocation[i])
      case geolocation[i] of
        'Latitude': lat = h5d_read(datasetid)
        'Longitude': lon = h5d_read(datasetid)
        'LocalSolarTime': ltime = h5d_read(datasetid)
        'Pressure': p = h5d_read(datasetid)
        'SolarZenithAngle': sza = h5d_read(datasetid)
        'Time': time = h5d_read(datasetid)
        'OrbitGeodeticAngle': orbitgeodeticangle = h5d_read(datasetid) 
      endcase
    endfor      
    h5d_close, datasetid

    ;read Data Fields
    datafields=['Convergence','L2gpPrecision','L2gpValue','Quality','Status']
    for i=0,n_elements(datafields)-1 do begin
      datasetid = h5d_open(hdfid, '/HDFEOS/SWATHS/'+str_species+'/Data Fields/'+datafields[i])
      case datafields[i] of
        'Convergence': convergence = h5d_read(datasetid)
        'L2gpPrecision': precision = h5d_read(datasetid)
        'L2gpValue': value = h5d_read(datasetid)
        'Quality': quality = h5d_read(datasetid)
        'Status': status = h5d_read(datasetid)
      endcase
    endfor      
    h5d_close, datasetid

  h5f_close, hdfid
  
  data = {swathname:str_species, $
          nlevels:n_elements(nlev), ntimes:n_elements(ntimes), tai93atozofgranule:tai93atozofgranule,$
          latitude:lat, longitude:lon, localsolartime:ltime, pressure:p, solarzenithangle:sza, time:time, orbitgeodeticangle:orbitgeodeticangle, $
          convergence:convergence, l2gpprecision:precision, l2gpvalue:value, quality:quality, status:status}
  
  return, data
end