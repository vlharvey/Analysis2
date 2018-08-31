function date2datecode, date
  doy=date2doy(date)
  year=strmid(strcompress(string(date),/r),0,4) 
  return, year+'d'+string(doy,format='(i3.3)')
end