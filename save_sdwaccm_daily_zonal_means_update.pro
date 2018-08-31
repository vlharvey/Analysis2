;
; read in all tem wbarstar and vbarstar from sdwaccm and all daily h1 T, U, V, Z, H2O and save daily zonal means
;
dir='/Volumes/cloud/data/WACCM_data/Datfiles_SD/sdwaccm2012-2014_1_2_2.cam.h1.'         ;YYYY-MM-DD-00000.nc
;
; read in all TEM diagnostics from 2 files
; DATE            LONG      = Array[390]
; ILEV            DOUBLE    = Array[89]
; LAT             DOUBLE    = Array[96]
; VSTAR           FLOAT     = Array[96, 89, 390]
; WSTAR           FLOAT     = Array[96, 89, 390]
;
restore,'/Volumes/cloud/data/WACCM_data/Datfiles_SD/sdwaccm2012-2014_1_2_2.cam.h4tem.2012.sav'
date_tem=date
result=size(VSTAR)
wstar1=reform(wstar,result(3),result(1),result(2))
vstar1=reform(vstar,result(3),result(1),result(2))
restore,'/Volumes/cloud/data/WACCM_data/Datfiles_SD/sdwaccm2012-2014_1_2_2.cam.h4tem.2013.sav'
date_tem=[date_tem,date]
result=size(VSTAR)
wstar2=reform(wstar,result(3),result(1),result(2))
vstar2=reform(vstar,result(3),result(1),result(2))
wstar_all=[wstar1,wstar2]
vstar_all=[vstar1,vstar2]
lev_tem=ilev	; level edges. T, U, V are at midpoint levels. be mindful that tem levels have 1 extra level
sdate_tem=strcompress(date_tem,/r)
ndays=n_elements(date_tem)
syear=strmid(sdate_tem,0,4)
smon=strmid(sdate_tem,4,2)
sday=strmid(sdate_tem,6,2)
;
; loop over dates with TEM diagnostics and save daily zonal means
;
for ii=0L,ndays-1L do begin
    result=file_search(dir+syear(ii)+'-'+smon(ii)+'-'+sday(ii)+'-00000.nc')
    if result(0) eq '' then goto,skipday
    print,result
;
; read 
;
    fname=result
    ncid = ncdf_open(fname)
    ncdf_varget, ncid, 'lat',     lat
    ncdf_varget, ncid, 'lev',     lev
    ncdf_varget, ncid, 'ilev',    ilev
    ncdf_varget, ncid, 'T',   t
    ncdf_varget, ncid, 'U',   u
    ncdf_varget, ncid, 'V',   v
    ncdf_varget, ncid, 'Z3',  z
    ncdf_varget, ncid, 'H2O', h2o
    ncdf_close, ncid
;
; zonal means
;
    tbar=mean(t,dim=1)
    ubar=mean(u,dim=1)
    vbar=mean(v,dim=1)
    zbar=mean(z,dim=1)
    h2obar=mean(h2o,dim=1)
    wstar=reform(wstar_all(ii,*,*))
    vstar=reform(vstar_all(ii,*,*))
;
; save daily zonal means
;
    save,file='/Volumes/cloud/data/WACCM_data/Datfiles_SD/sdwaccm2012-2014_1_2_2.cam.zm.'+sdate_tem(ii)+'.sav',lev,ilev,lat,vstar,wstar,tbar,ubar,vbar,zbar,h2obar

skipday:
endfor
end
