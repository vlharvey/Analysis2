;
; read in all tem wbarstar and vbarstar from sdwaccm and all daily h5 T, U, V, Z, H2O and save daily zonal means
;
dir='/Volumes/cloud/data/WACCM_data/Datfiles_SD/f_1975-2010_2deg_refc1sd_wa4_tsmlt.002.cam.h5.'         ;20120131.sav
;
; read in all TEM diagnostics from 1 file
;
restore,'/Volumes/cloud/data/WACCM_data/Datfiles_SD/f_1975-2010_2deg_refc1sd_wa4_tsmlt.002.cam.h4tem.1975-2013.sav'
good=where(date lt 20000101L,ndays)
date_tem=date(good)
wstar_all=wstar(good,*,*)
vstar_all=vstar(good,*,*)
lev_tem=ilev	; level edges. T, U, V are at midpoint levels. be mindful that tem levels have 1 extra level
sdate_tem=strcompress(date_tem,/r)
;
; loop over dates with TEM diagnostics and save daily zonal means
;
for ii=0L,ndays-1L do begin
    result=file_search(dir+sdate_tem(ii)+'.sav')
    if result(0) eq '' then goto,skipday
    print,result
    restore,result
;
; zonal means
;
    tbar=mean(t,dim=1)
    ubar=mean(u,dim=1)
    vbar=mean(v,dim=1)
    zbar=mean(z,dim=1)
    h2obar=mean(h2o,dim=1)
    wstar=transpose(reform(wstar_all(ii,*,*)))
    vstar=transpose(reform(vstar_all(ii,*,*)))
;
; save daily zonal means
;
    save,file='/Volumes/cloud/data/WACCM_data/Datfiles_SD/f_1975-2010_2deg_refc1sd_wa4_tsmlt.002.cam.zm.'+sdate_tem(ii)+'.sav',lev,ilev,lat,vstar,wstar,tbar,ubar,vbar,zbar,h2obar

skipday:
endfor
end
