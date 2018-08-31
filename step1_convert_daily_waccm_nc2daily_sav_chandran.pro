;
; read SD-WACCM netcdf data from Amal Chandran with the following variables from ncdump -h
; convert monthly to daily files
;
; aura% pwd
; /aura7/harvey/WACCM_data/Datfiles/Datfiles_SD
; -rw-r--r--@ 1 harvey  staff   2.0G Sep 26 12:20 wa4_cesm1_1_b02_geos5_2x_sim153f.cam2.h3.1998-10-17-00000.nc.nc
; -rw-r--r--@ 1 harvey  staff   2.0G Sep 26 12:24 wa4_cesm1_1_b02_geos5_2x_sim153f.cam2.h3.1998-12-29-00000.nc.nc
; -rw-r--r--@ 1 harvey  staff   2.0G Sep 26 15:05 wa4_cesm1_1_b02_geos5_2x_sim153f.cam2.h3.2002-08-04-00000.nc.nc
; -rw-r--r--@ 1 harvey  staff   2.0G Sep 26 13:39 wa4_cesm1_1_b02_geos5_2x_sim153f.cam2.h3.2005-12-27-00000.nc.nc
; -rw-r--r--@ 1 harvey  staff   2.0G Sep 26 13:55 wa4_cesm1_1_b02_geos5_2x_sim153f.cam2.h3.2008-12-26-00000.nc.nc
;
dir='/Volumes/earth/harvey/WACCM_data/Datfiles/Datfiles_SD/'
spawn,'ls '+dir+'wa4_cesm1_1_b02_geos5_2x_sim153f.cam2.h3.*-00000.nc.nc',ncfiles
nfile=n_elements(ncfiles)
for ifile=0L,nfile-1L do begin

ncfile=ncfiles(ifile)
print,'opening '+ncfile
ncid=ncdf_open(ncfile)
result0=ncdf_inquire(ncid)
for idim=0,result0.ndims-1 do begin
    ncdf_diminq,ncid,idim,name,dim
    if name eq 'lon' then nc=dim
    if name eq 'lat' then nr=dim
    if name eq 'lev' then nl=dim
    if name eq 'time' then nt=dim
;   print,'read ',name,' dimension ',dim
endfor
for ivar=0,result0.nvars-1 do begin
    result=ncdf_varinq(ncid,ivar)
;   print,result.name
    ncdf_varget,ncid,ncdf_varid(ncid,result.name),data

    if result.name eq 'P0' then p0=data
    if result.name eq 'PS' then ps=data
    if result.name eq 'lat' then alat=data
    if result.name eq 'lon' then alon=data
    if result.name eq 'lev' then lev=data
    if result.name eq 'ilev' then ilev=data
    if result.name eq 'time' then time=data
    if result.name eq 'hyai' then hyai=data
    if result.name eq 'hybi' then hybi=data
    if result.name eq 'hyam' then hyam=data
    if result.name eq 'hybm' then hybm=data
    if result.name eq 'date' then begin
       date=data
       sdate=strcompress(date,/remove_all)
    endif
    if result.name eq 'T' then tgrd=data
    if result.name eq 'U' then ugrd=data
    if result.name eq 'V' then vgrd=data
    if result.name eq 'Z3' then zgrd=data/1000.
endfor	; loop over variables
ncdf_close,ncid
;
; only 00Z data
;
time=0.
stime=string(format='(i2.2)',long(time))
;
; Calculate 3d Pressure: p(i,j,k) = A(k)*PO + B(k)*PS(i,j) in Pascals
;
    pgrd=fltarr(nc,nr,nl)
    Pzero=P0
    FOR ilon=0,nc-1 DO $
        FOR ilat=0,nr-1 DO $
            FOR ialt=0,nl-1 DO $
                pgrd(ilon,ilat,ialt)=(hyam(ialt)*Pzero + hybm(ialt)*PS(ilon,ilat)) / 100.
;
; IDL save file for each output time
;
    for n=0L,nt-1L do begin
        ofile=dir+'wa4_cesm1_1_b02_geos5_2x_sim153f.cam2.h3.'+sdate(n)+'.sav'
        print,ofile
stop
;       save,file=ofile,alon,alat,pgrd,tgrd,ugrd,vgrd
    endfor	; loop over days
endfor	; loop over files
end
