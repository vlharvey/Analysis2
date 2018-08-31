;
; read in monthly tem wbarstar and vbarstar from sdwaccm continuation run and save out yearly files
;
dir='/Volumes/cloud/data/WACCM_data/Datfiles_SD/sdwaccm2012-2014_1_2_2.cam.h4tem.'
spawn,'ls '+dir+'2013*.sav',ifiles
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
endfor
date=date_all
vstar=transpose(vstar_all)
wstar=transpose(wstar_all)
savfile = dir+'2013.sav
save,file=savfile,date,ilev,lat,vstar,wstar
end
