;
; check L3 daisy sizes
;
; restore CIPS procedures and functions
;
restore,'read_cips_file.sav
pth='/Volumes/Data/CIPS_data/Datfiles/cips_sci_3a_'
spawn,'ls '+pth+'*'+'_v04.20_r04.nc',ifiles
nfiles=n_elements(ifiles)
;
; loop over days
;
for ifile=0,nfiles-1 do begin
;
    fname=ifiles(ifile)
    data = read_cips_file(fname)
;
; extract variables from data structure
;
    alat=data.latitude
    alon=data.longitude
    ut_date=data.ut_date
    result=size(alat)
    nc=result(1)
    nr=result(2)
if nc ne 1953L or nr ne 1953L then print,fname,' ',ut_date,' ',nc,' ',nr
ENDFOR   ; loop over days
end
