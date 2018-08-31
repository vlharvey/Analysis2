function utc2ltime, utc, lon
  utc_tmp = utc
  lon_tmp = lon

  ;correct utc values out of range, force to cyclic 0-24 interval
  index = where(utc_tmp le -24 or utc_tmp ge 24, nindex)
  if nindex ne 0 then utc_tmp[index] = utc_tmp[index] mod 24
  index = where(utc_tmp lt 0, nindex)
  if nindex ne 0 then utc_tmp[index] = 24 + utc_tmp[index]
  
  ;correct lon values out of range, here a -180 - 180 range is needed
  index = where(lon_tmp gt 180, nindex)
  if nindex ne 0 then lon_tmp[index] = lon_tmp[index] - 360 
  
  ;calculate local time
  ltime = utc_tmp + lon_tmp / 15.
  
  ;compensate day change
  index = where(ltime ge 24, nindex)
  if nindex ne 0 then ltime[index] = ltime[index] mod 24  
  index = where(ltime lt  0, nindex)
  if nindex ne 0 then ltime[index] = 24 + ltime[index]
  
  return, ltime
end