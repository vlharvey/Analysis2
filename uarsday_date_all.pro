month=[           'sep_91','oct_91','nov_91','dec_91',$
'jan_92','feb_92','mar_92','apr_92','may_92','jun_92',$
'jul_92','aug_92','sep_92','oct_92','nov_92','dec_92',$
'jan_93','feb_93','mar_93','apr_93','may_93','jun_93',$
'jul_93','aug_93','sep_93','oct_93','nov_93','dec_93',$
'jan_94','feb_94','mar_94','apr_94','may_94','jun_94',$
'jul_94','aug_94','sep_94','oct_94','nov_94','dec_94',$
'jan_95','feb_95','mar_95','apr_95','may_95','jun_95',$
'jul_95','aug_95','sep_95','oct_95','nov_95','dec_95',$
'jan_96','feb_96','mar_96','apr_96','may_96','jun_96',$
'jul_96','aug_96','sep_96','oct_96','nov_96','dec_96',$
'jan_97','feb_97','mar_97','apr_97','may_97','jun_97',$
'jul_97','aug_97','sep_97','oct_97','nov_97','dec_97',$
'jan_98','feb_98','mar_98','apr_98','may_98','jun_98',$
'jul_98','aug_98','sep_98','oct_98','nov_98','dec_98']
days=[                 20 ,     31 ,     30 ,     31 ,$
     31 ,     29 ,     31 ,     30 ,     31 ,     30 ,$
     31 ,     31 ,     30 ,     31 ,     30 ,     31,$
     31 ,     28 ,     31 ,     30 ,     31 ,     30 ,$
     31 ,     31 ,     30 ,     31 ,     30 ,     31,$
     31 ,     28 ,     31 ,     30 ,     31 ,     30 ,$
     31 ,     31 ,     30 ,     31 ,     30 ,     31,$
     31 ,     28 ,     31 ,     30 ,     31 ,     30 ,$
     31 ,     31 ,     30 ,     31 ,     30 ,     31,$
     31 ,     29 ,     31 ,     30 ,     31 ,     30 ,$
     31 ,     31 ,     30 ,     31 ,     30 ,     31,$
     31 ,     28 ,     31 ,     30 ,     31 ,     30 ,$
     31 ,     31 ,     30 ,     31 ,     30 ,     31,$
     31 ,     28 ,     31 ,     30 ,     31 ,     30 ,$
     31 ,     31 ,     30 ,     31 ,     30 ,     31]
uars_mon=[         20.0000, 51.0000, 81.0000, 112.000,$
 143.000, 172.000, 203.000, 233.000, 264.000, 294.000,$
 325.000, 356.000, 386.000, 417.000, 447.000, 478.000,$
 509.000, 537.000, 568.000, 598.000, 629.000, 659.000,$
 690.000, 721.000, 751.000, 782.000, 812.000, 843.000,$
 874.000, 902.000, 933.000, 963.000, 994.000, 1024.00,$
 1055.00, 1086.00, 1116.00, 1147.00, 1177.00, 1208.00,$
 1239.00, 1267.00, 1298.00, 1328.00, 1359.00, 1389.00,$
 1420.00, 1451.00, 1481.00, 1512.00, 1542.00, 1573.00,$
 1604.00, 1633.00, 1664.00, 1694.00, 1725.00, 1755.00,$
 1786.00, 1817.00, 1847.00, 1878.00, 1908.00, 1939.00,$
 1970.00, 1998.00, 2029.00, 2059.00, 2090.00, 2120.00,$
 2151.00, 2182.00, 2212.00, 2243.00, 2273.00, 2304.00,$
 2335.00, 2363.00, 2394.00, 2424.00, 2455.00, 2485.00,$
 2516.00, 2547.00, 2577.00, 2608.00, 2638.00, 2669.00]

;nmonth=n_elements(days)
;uars_mon=fltarr(nmonth)
;;
;; uars day 1 is sep 11 1991
;;
;for n=0,nmonth-1 do begin
;uars_mon(n)=total(days(0:n))
;endfor
;ndays=total(days)
day=0
udays=36+findgen(2669) 
openw,99,'dates_all.fil'
for mu=0,2668 do begin
day1=udays(mu)
;read,'Input UARS day',day
day=day1
index=where(uars_mon gt day)
m=index(0)
date=days(m)-uars_mon(m)+day+1
mon1=strmid(month(m),0,3)
year1=strmid(month(m),4,2)
if date lt 10. then sdate1='0'+strmid(string(date),6,1)
if date ge 10. then sdate1=strmid(string(date),6,2)
if sdate1 eq '01' then begin
print,'UARS day ',day1,'=',sdate1,' ',mon1,' ',year1 
printf,99,'UARS day ',day1,'=',sdate1,' ',mon1,' ',year1 
endif
endfor
end
