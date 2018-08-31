function date2julday, yyyymmdd
  date=yyyymmdd
  return, julday(strmid(strcompress(date,/r),4,2),strmid(strcompress(date,/r),6,2),strmid(strcompress(date,/r),0,4))
end