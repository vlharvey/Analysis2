;
; read MLS gridded IDL save files.
; Input data:  IDL> restore,'/Volumes/earth/aura6/data/MLS_data/Datfiles_Grid/MLS_grid5_ALL_v3.3_YYYYMMDD.sav
; Input data:  IDL> restore,'/Volumes/earth/aura6/data/MLS_data/Datfiles_Grid/MLS_grid5_U_V_v3.3_YYYYMMDD.sav
; Check for missing files
; 
sver = 'v3.3'
dirw='/Volumes/earth/aura6/data/MLS_data/Datfiles_Grid/MLS_grid5_'
ifiles=file_search(dirw+'ALL_'+sver+'_????????.sav')
ufiles=file_search(dirw+'U_V_'+sver+'_????????.sav')
nfile=n_elements(ifiles)
nfile2=n_elements(ufiles)
if nfile ne nfile2 then stop,'Check number of all and uv files'
;
; loop over files
;
FOR n=0l,nfile-1l DO BEGIN

result=strsplit(ifiles(n),'_',/extract)
result2=strsplit(result(-1),'.',/extract)
sdate=result2(0)
dum=findfile('/Volumes/earth/aura6/data/MLS_data/Datfiles_Grid/MLS_grid5_U_V_v3.3_'+sdate+'.sav')
if dum(0) eq '' then print,'missing UV '+sdate
endfor
end
