;
; print missing rad files
; Input data:  IDL> restore,'/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_press_YYYYMMDD.sav
; Input data:  IDL> restore,'/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_rad_YYYYMMDD.sav
; 
dirw='/Volumes/Data/MERRA_data/Datfiles/MERRA-on-WACCM_'
ifiles=file_search(dirw+'press_????????.sav')
nfile=n_elements(ifiles)
;
; loop over files
;
FOR n=0l,nfile-1l DO BEGIN
    result=strsplit(ifiles(n),'_',/extract)
    result2=strsplit(result(3),'.',/extract)
    sdate=result2(0)
    dum=findfile(dirw+'rad_'+sdate+'.sav')
    if dum(0) eq '' then print,sdate
ENDFOR
end
