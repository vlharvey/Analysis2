function date2doy, date
  ldate = long(date)
  year=floor(ldate/10000.)
  month=floor((ldate-year*10000)/100.)
  day=ldate-year*10000-month*100
  doy=julday(month, day, year)-julday(1,1,year)+1
  return,doy
end