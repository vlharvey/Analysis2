function today, dummy
  ;returns todays date in the format yyyymmdd
  ;dummy is not to be used, call function like in "print, today()"
  return, strmid(systime(),3,4,/reverse_offset)+$
         strcompress(string(where(['','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'] eq strmid(systime(),4,3)),format='(i2.2)'),/remove_all)+$
         string(strcompress(strmid(systime(),8,2),/r),format='(i2.2)')
end