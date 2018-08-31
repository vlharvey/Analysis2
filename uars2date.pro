pro uars2date,imn,idy,iyr,uday
;
; return (day,month,year) information given UARS day
; --- UARS day 1 = September 12, 1991 (jday=255)
;
month=[31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.]
ibase=1991
leapdays=[171.,1632.,3093.,4554.,6015.]
year1=[1992.,1993.,1994.,1995.,1996.,1997.,1998.,1999.,2000.,2001.,$
       2002.,2003.,2004.,2005.,2006.,2007.,2008.,2009.,2010.,2011.]
jday1=[112.0,478.0,843.0,1208.,1573.,1939.,2304.,2669.,3034.,3400.,$
       3765.,4130.,4495.,4861.,5226.,5591.,5956.,6322.,6687.,7052.]
nyear=n_elements(jday1)
if uday lt jday1(0) then begin
   jday=uday+254.
   iyr=1991.
endif
for i=0L,nyear-1L do begin
    if uday ge jday1(i) then begin
       jday=uday-jday1(i)+1.0
       iyr=year1(i)
    endif
endfor
if (iyr mod 4) eq 0 then month(1)=29.
idy=jday
imn=1.0
if jday le month(0) then return
i=0
for j=0,11 do begin
    i=i*1.0+month(j)
    if i gt jday then begin
       imn=j+1.0
       idy=jday+month(j)-i*1.0
       return
    endif
    if i eq jday then begin
       imn=j+1.0
       idy=month(j)
       return
    endif
endfor
end
