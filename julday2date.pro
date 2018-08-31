function julday2date, julian
  jul=round(julian)
  caldat, jul, month, day, year
  date=year*10000+month*100+day
  return, date
end