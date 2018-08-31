;
; read in monthly tem wbarstar and vbarstar from sdwaccm continuation run and save out yearly files
;
dir='/Volumes/cloud/data/WACCM_data/Datfiles_SD/f_1975-2010_2deg_refc1sd_wa4_tsmlt.002.cam.h4tem.'
spawn,'ls '+dir+'*-00000.sav',ifiles
for ifile=0L,n_elements(ifiles)-1L do begin
    restore,ifiles(ifile)
    print,ifiles(ifile)
    if ifile eq 0L then wstar_all=transpose(wstar)
    if ifile eq 0L then vstar_all=transpose(vstar)
    if ifile eq 0L then date_all=date
    if ifile gt 0L then begin
        wstar_all=[wstar_all,transpose(wstar)]
        vstar_all=[vstar_all,transpose(vstar)]
        date_all=[date_all,date]
    endif
help,wstar_all,date_all
endfor
;
; remove times after April 30th 2012 (end is 1 July)
;
good=where(date_all le 20120430L)
date_all=date_all(good)
wstar_all=wstar_all(good,*,*)
vstar_all=vstar_all(good,*,*)
;
; append VLH continuation of 2012 and 2013 at May 1 2012 so that NH summer of 2012 is entirely update
;
restore,'/Volumes/cloud/data/WACCM_data/Datfiles_SD/sdwaccm2012-2014_1_2_2.cam.h4tem.2012.sav
good=where(date ge 20120501L)
date_all=[date_all,date(good)]
wstar_all=[wstar_all,transpose(wstar(*,*,good))]
vstar_all=[vstar_all,transpose(vstar(*,*,good))]

restore,'/Volumes/cloud/data/WACCM_data/Datfiles_SD/sdwaccm2012-2014_1_2_2.cam.h4tem.2013.sav
date_all=[date_all,date]
wstar_all=[wstar_all,transpose(wstar)]
vstar_all=[vstar_all,transpose(vstar)]
;
; check monotonically increasing jday
;
year=long(strmid( strcompress(date_all,/r),0,4))
mon=long(strmid( strcompress(date_all,/r),4,2))
day=long(strmid( strcompress(date_all,/r),6,2))
jday=julday(mon,day,year)

date=date_all
vstar=vstar_all
wstar=wstar_all
save,file='/Volumes/cloud/data/WACCM_data/Datfiles_SD/f_1975-2010_2deg_refc1sd_wa4_tsmlt.002.cam.h4tem.1975-2013.sav',date,ilev,lat,vstar,wstar
end
