;
; read yearly netcdf data one day at a time
;
ncfile='netcdf_filename.nc'
ncid=ncdf_open(ncfile)
result0=ncdf_inquire(ncid)
for idim=0,result0.ndims-1 do begin
    ncdf_diminq,ncid,idim,name,dim
    if name eq 'lon' then nc=dim
    if name eq 'lat' then nr=dim
    if name eq 'lev' then nl=dim
    if name eq 'time' then nt=dim
    print,'read ',name,' dimension ',dim
endfor
for ivar=0,result0.nvars-1 do begin
    result=ncdf_varinq(ncid,ivar)
    help,result.name
    if result.name eq 'H2O_VMR_inst' then begin
       count = [nc,nr,nl,1]
       for n=0L,nt-1L do begin
           offset = [0,0,0,n]
           ncdf_varget,ncid,ncdf_varid(ncid,result.name),h2o,count=count,offset=offset
help,h2o
       endfor
    endif
endfor		; loop over variables
ncdf_close,ncid
end
