;
; read CIPS 1b netcdf files and archive in IDL save format
;
comment=strarr(25)
comment(0)='dimensions: nLayers_3=1, string=2, nLayers=108, dim2_IMAGE_TLM_TIMESTAMP=4, dim3_IMAGE_TLM_TIMESTAMP=27, nCols=216, nRows=208, nCols_2=2060, nRows_2=319, dim1_VIEW_ANGLE=216, dim2_VIEW_ANGLE=208, dim3_VIEW_ANGLE=108, dim1_CAMERA_ID=108'
comment(1)='Variables:'
comment(2)='short AIM_ORBIT_NUMBER(nLayers_3) cumulative mission orbit number
comment(3)='byte VERSION(nLayers_3) the version number of this data
comment(4)='short STACK_ID(nLayers_3) used by the CIPS science processing system to uniquely identify this Level 1B data
comment(5)='char UT_DATE(nLayers_3, string) ut date of this orbit
comment(6)='short NROWS(nLayers_3) the number of rows in image stack, nominally 250
comment(7)='short NCOLS(nLayers_3) the number of cols in image stack, nominally 450
comment(8)='short NLAYERS(nLayers_3) the number of layers in this image stack
comment(9)='double STACK_START_TIME(nLayers_3) timestamp in gps microseconds of the first image used to create this data file
comment(10)='double UT_TIME(nLayers_3) ut time in fractional hours for each layer
comment(11)='char LA_TIME(nLayers_3, string) ascii time per layer (e.g. 2000/001-00:00:00), provided for convenience
comment(12)='double IMAGE_TLM_TIMESTAMP(nLayers_3, dim3_IMAGE_TLM_TIMESTAMP, dim2_IMAGE_TLM_TIMESTAMP, nLayers) gps timestamp of each scene (in gps microseconds)
comment(13)='double ALBEDO(nLayers_3, nLayers, nRows, nCols) fillValue = -1
comment(14)='double ALBEDO_UNCERTAINTY(nLayers_3, nLayers, nRows, nCols) fillValue = -1
comment(15)='byte COMMON_VOLUME_MAP(nLayers_3, nRows_2, nCols_2) Ô1Õ indicates a tile in the CVO, Ô0Õ otherwise
comment(16)='int QUALITY_FLAGS(nLayers_3) TBD
comment(17)='double CENTER_LON(nLayers_3) center longitude of map projection, NOT data. used for orienting the data horizontally
comment(18)='int BBOX(nLayers_3, 4) {x, y} bounding box of image stack
comment(19)='short SCENE_BBOXES(nLayers_3, 4, nLayers) {x, y} bounding boxes for each stack layer
comment(20)='double SCATTERING_ANGLE(nLayers_3, nLayers, nRows, nCols) units=degrees, fillValue=-1
comment(21)='double VIEW_ANGLE(nLayers_3, dim3_VIEW_ANGLE, dim2_VIEW_ANGLE, dim1_VIEW_ANGLE) fillValue=-1
comment(22)='double ZENITH_ANGLE(nLayers_3, nLayers, nRows, nCols) units =degrees, fillValue=-1
comment(23)='byte COMMON_VOLUME_FLAG(nLayers_3) COMMON_VOLUME_FLAG:fillValue=1
comment(24)='char CAMERA_ID(nLayers_3, dim1_CAMERA_ID, string)

dir='/aura7/harvey/CIPS_data/Datfiles/'
spawn,'ls '+dir+'cips_sci_1b_*nc',ncfiles
nfiles=n_elements(ncfiles)
for n=0L,nfiles-1L do begin
;
; build filenames
;
ncfile=ncfiles(n)
result=strsplit(ncfile,'/',/extract)
result2=strsplit(result(4),'.',/extract)
ofile=dir+result2(0)+'.sav'
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
    if name eq 'nLayers_3' then nlayers_3=dim
    if name eq 'string' then string_dimension=dim
    if name eq 'nLayers' then nlayers=dim
    if name eq 'dim2_IMAGE_TLM_TIMESTAMP' then dim2_IMAGE_TLM_TIMESTAMP=dim
    if name eq 'dim3_IMAGE_TLM_TIMESTAMP' then dim3_IMAGE_TLM_TIMESTAMP=dim
    if name eq 'nCols' then ncols=dim
    if name eq 'nRows' then nrows=dim
    if name eq 'nCols_2' then ncols_2=dim
    if name eq 'nRows_2' then nrows_2=dim
    if name eq 'dim1_VIEW_ANGLE' then dim1_VIEW_ANGLE=dim
    if name eq 'dim2_VIEW_ANGLE' then dim2_VIEW_ANGLE=dim
    if name eq 'dim3_VIEW_ANGLE' then dim3_VIEW_ANGLE=dim
    if name eq 'dim1_CAMERA_ID' then dim1_CAMERA_ID=dim
    print,'read ',name,' dimension ',dim
endfor
for ivar=0,result.nvars-1 do begin
    result=ncdf_varinq(ncid,ivar)
    ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
    if result.name eq 'AIM_ORBIT_NUMBER' then AIM_ORBIT_NUMBER=data
    if result.name eq 'VERSION' then VERSION=data
    if result.name eq 'STACK_ID' then STACK_ID=data
    if result.name eq 'UT_DATE' then UT_DATE=data
    if result.name eq 'NROWS' then NROWS=data
    if result.name eq 'NCOLS' then NCOLS=data
    if result.name eq 'NLAYERS' then NLAYERS=data
    if result.name eq 'STACK_START_TIME' then STACK_START_TIME=data
    if result.name eq 'UT_TIME' then UT_TIME=data
    if result.name eq 'LA_TIME' then LA_TIME=data
    if result.name eq 'IMAGE_TLM_TIMESTAMP' then IMAGE_TLM_TIMESTAMP=data
    if result.name eq 'ALBEDO' then ALBEDO=data
    if result.name eq 'ALBEDO_UNCERTAINTY' then ALBEDO_UNCERTAINTY=data
    if result.name eq 'COMMON_VOLUME_MAP' then COMMON_VOLUME_MAP=data
    if result.name eq 'QUALITY_FLAGS' then QUALITY_FLAGS=data
    if result.name eq 'CENTER_LON' then CENTER_LON=data
    if result.name eq 'BBOX' then BBOX=data
    if result.name eq 'SCENE_BBOXES' then SCENE_BBOXES=data
    if result.name eq 'SCATTERING_ANGLE' then SCATTERING_ANGLE=data
    if result.name eq 'VIEW_ANGLE' then VIEW_ANGLE=data
    if result.name eq 'ZENITH_ANGLE' then ZENITH_ANGLE=data
    if result.name eq 'COMMON_VOLUME_FLAG' then COMMON_VOLUME_FLAG=data
    if result.name eq 'CAMERA_ID' then CAMERA_ID=data
    print,'read ',result.name,' variable'
endfor	; loop over variables
;
; write CIPS data in IDL save format
;
save,file=ofile,/variables

endfor	; loop over files
end
