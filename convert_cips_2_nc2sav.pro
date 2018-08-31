;
; read CIPS 2 netcdf files and archive in IDL save format
;
comment=strarr(25)
comment(0)='dimensions: structure_elements=1, dim1_BBOX=4, nCols=1738, nRows=289, dim1_VIEW_ANGLE=1738, dim2_VIEW_ANGLE=289, dim1_SCATTERING_ANGLE=1738, dim2_SCATTERING_ANGLE=289, dim1_CLD_ALBEDO=1738, dim2_CLD_ALBEDO=289'
comment(1)='Variables:'
comment(2)='short AIM_ORBIT_NUMBER(structure_elements) cumulative mission orbit number
comment(3)='byte VERSION(structure_elements) the version number of this data
comment(4)='int UT_DATE(structure_elements) ut date of this orbit
comment(5)='double ORBIT_START_TIME(structure_elements) start timestamp in gps microseconds of the orbit corresponding with this data
comment(6)='double ORBIT_START_TIME_UT(structure_elements) ut timestamp of the start of the orbit corresponding with this data
comment(7)='double ORBIT_END_TIME(structure_elements) end timestamp in gps microseconds of the orbit corresponding with this data
comment(8)='short XDIM(structure_elements) 
comment(9)='short YDIM(structure_elements)
comment(10)='byte NLAYERS(structure_elements) fillValue=1
comment(11)='byte UT_TIME(structure_elements) fillValue=1
comment(12)='int QUALITY_FLAGS(structure_elements)
comment(13)='short X_TILE_DIM(structure_elements)
comment(14)='short Y_TILE_DIM(structure_elements)
comment(15)='int BBOX(structure_elements, dim1_BBOX)
comment(16)='double CENTER_LON(structure_elements)
comment(17)='byte LATITUDE(structure_elements) fillValue=1
comment(18)='byte LONGITUDE(structure_elements) fillValue=1
comment(19)='double ZENITH_ANGLE(structure_elements, nRows, nCols)
comment(20)='double VIEW_ANGLE(structure_elements, dim2_VIEW_ANGLE, dim1_VIEW_ANGLE)
comment(21)='double SCATTERING_ANGLE(structure_elements, dim2_SCATTERING_ANGLE, dim1_SCATTERING_ANGLE)
comment(22)='double CLD_ALBEDO(structure_elements, dim2_CLD_ALBEDO, dim1_CLD_ALBEDO)

dir='/aura7/harvey/CIPS_data/Datfiles/'
spawn,'ls '+dir+'cips_sci_2*nc',ncfiles
nfiles=n_elements(ncfiles)
for n=0L,nfiles-1L do begin
;
; build filenames
;
ncfile=ncfiles(n)
result=strsplit(ncfile,'/',/extract)
result2=strsplit(result(4),'.',/extract)
ofile=dir+result2(0)+'.'+result2(1)+'.sav'
print,' '
print,ncfile
print,ofile
;
; open and read netcdf data
;
ncid=ncdf_open(ncfile)
result=ncdf_inquire(ncid)
for idim=0,result.ndims-1 do begin
    ncdf_diminq,ncid,idim,name,dim
    if name eq 'structure_elements' then structure_elements=dim
    if name eq 'dim1_BBOX' then dim1_BBOX=dim
    if name eq 'nCols' then ncols=dim
    if name eq 'nRows' then nrows=dim
    if name eq 'dim1_VIEW_ANGLE' then dim1_VIEW_ANGLE=dim
    if name eq 'dim2_VIEW_ANGLE' then dim2_VIEW_ANGLE=dim
    if name eq 'dim1_SCATTERING_ANGLE' then dim1_SCATTERING_ANGLE=dim
    if name eq 'dim2_SCATTERING_ANGLE' then dim2_SCATTERING_ANGLE=dim
    if name eq 'dim1_CLD_ANGLE' then dim1_CLD_ANGLE=dim
    if name eq 'dim2_CLD_ANGLE' then dim2_CLD_ANGLE=dim
    if name eq 'dim1_CLD_ALBEDO' then dim1_CLD_ALBEDO=dim
    if name eq 'dim2_CLD_ALBEDO' then dim2_CLD_ALBEDO=dim
    print,'read ',name,' dimension ',dim
endfor
for ivar=0,result.nvars-1 do begin
    result=ncdf_varinq(ncid,ivar)
    ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
    if result.name eq 'AIM_ORBIT_NUMBER' then AIM_ORBIT_NUMBER=data
    if result.name eq 'VERSION' then VERSION=data
    if result.name eq 'UT_DATE' then UT_DATE=data
    if result.name eq 'ORBIT_START_TIME' then ORBIT_START_TIME=data
    if result.name eq 'ORBIT_START_TIME_UT' then ORBIT_START_TIME_UT=data
    if result.name eq 'ORBIT_END_TIME' then ORBIT_END_TIME=data
    if result.name eq 'XDIM' then XDIM=data
    if result.name eq 'YDIM' then YDIM=data
    if result.name eq 'NLAYERS' then NLAYERS=data
    if result.name eq 'UT_TIME' then UT_TIME=data
    if result.name eq 'QUALITY_FLAGS' then QUALITY_FLAGS=data
    if result.name eq 'X_TILE_DIM' then X_TILE_DIM=data
    if result.name eq 'Y_TILE_DIM' then Y_TILE_DIM=data
    if result.name eq 'BBOX' then BBOX=data
    if result.name eq 'CENTER_LON' then CENTER_LON=data
    if result.name eq 'LATITUDE' then LATITUDE=data
    if result.name eq 'LONGITUDE' then LONGITUDE=data
    if result.name eq 'ZENITH_ANGLE' then ZENITH_ANGLE=data
    if result.name eq 'VIEW_ANGLE' then VIEW_ANGLE=data
    if result.name eq 'SCATTERING_ANGLE' then SCATTERING_ANGLE=data
    if result.name eq 'CLD_ALBEDO' then CLD_ALBEDO=data
    print,'read ',result.name,' variable'
endfor	; loop over variables
;
; write CIPS data in IDL save format
;
save,file=ofile,/variables
stop
endfor	; loop over files
end
