function doy2date, doy, year
  doy = long(doy)
  caldat,doy+julday(1,1,year)-1,m,d,y
  yyyymmdd=year*10000l+m*100l+d*1l
  return,yyyymmdd
end