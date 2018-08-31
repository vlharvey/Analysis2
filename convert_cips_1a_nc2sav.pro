;
; read CIPS 1a netcdf files and archive in IDL save format
;
comment=strarr(34)
comment(0)='Dimensions: structure_elements=27, string=17, nCols=170, nRows=340'
comment(1)='Variables:'
comment(2)='short AIM_ORBIT_NUMBER(structure_elements) cumulative mission orbit number'
comment(3)='byte VERSION(structure_elements) the version number of this data'
comment(4)='char UT_DATE(structure_elements, string) ut date of this orbit'
comment(5)='char CAMERA_ID(structure_elements, string) values=PX:+X, MX:-X, PY:+Y Nadir, MY:-Y Nadir'
comment(6)='short NROWS(structure_elements) the max number of rows found in albedo images in current file'
comment(7)='short NCOLS(structure_elements) the max number of columns found in albedo images in current file'
comment(8)='byte OBSERVATION_TYPE(structure_elements) 1=Dark, 2=First Light, 3=Primary Science, 4=Mapping, 5=Flat Field, 6=Last Light'
comment(9)='double UT_TIME(structure_elements) ut time in fractional hours for each image'
comment(10)='char LA_TIME(structure_elements, string) ascii time per image (e.g. 2000/001-00:00:00), provided for convenience'
comment(11)='double IMAGE_TLM_TIMESTAMP(structure_elements) time used to associate images together as part of a scene (in gps microsecs)'
comment(12)='double SHUTTER_OPEN_TIMESTAMP(structure_elements) time the shutter actually opened (in gps seconds)'
comment(13)='double SHUTTER_CLOSE_TIMESTAMP(structure_elements) time the shutter actually closed (in gps seconds)'
comment(14)='double INTEGRATION_TIME(structure_elements) number of seconds that the shutter was open'
comment(15)='double MID_IMAGE_TIMESTAMP(structure_elements) mid point of the integration (in gps seconds)'
comment(16)='double ALBEDO(structure_elements, nRows, nCols)'
comment(17)='double ALBEDO_UNCERTAINTY(structure_elements, nRows, nCols)'
comment(18)='byte COMMON_VOLUME_FLAG(structure_elements) does an image have any pixels that are in the common volume region? 1=yes, 0=no'
comment(19)='byte COMMON_VOLUME_MAP(structure_elements) Ô2Õ indicates a tile in the CVO at the cloud deck altitude'
comment(20)='int QUALITY_FLAGS(structure_elements) TBD'
comment(21)='double SC_LATITUDE(structure_elements) spacecraft geodetic latitude'
comment(22)='double SC_LONGITUDE(structure_elements) spacecraft geodetic latitude'
comment(23)='double SC_ALTITUDE(structure_elements) spacecraft altitude'
comment(24)='float ALONGTRACK_UNCERTAINTY(structure_elements) units=km'
comment(25)='float CROSSTRACK_UNCERTAINTY(structure_elements) units=km'
comment(26)='float RADIAL_UNCERTAINTY(structure_elements) units=km'
comment(27)='float TRACK_AZIMUTH(structure_elements) units=degrees'
comment(28)='short GEOMETRY_DOWNSAMPLE(structure_elements) Amount of downsampling in lat, lon, and viewing angles. Nominally=10'
comment(29)='double LATITUDE(structure_elements, nRows, nCols) downsampled by geometry_downsample quantity to reduce file size'
comment(30)='double LONGITUDE(structure_elements, nRows, nCols) downsampled by geometry_downsample quantity to reduce file size'
comment(31)='double SCATTERING_ANGLE(structure_elements, nRows, nCols) downsampled by geometry_downsample quantity to reduce file size'
comment(32)='double VIEWING_ANGLE(structure_elements, nRows, nCols) downsampled by geometry_downsample quantity to reduce file size'
comment(33)='double ZENITH_ANGLE(structure_elements, nRows, nCols) downsampled by geometry_downsample quantity to reduce file size'

dir='/aura7/harvey/CIPS_data/Datfiles/'
spawn,'ls '+dir+'cips_sci_1a_*nc',ncfiles
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
    if name eq 'nCols' then ncols=dim
    if name eq 'nRows' then nrows=dim
    if name eq 'structure_elements' then structure_elements=dim
    if name eq 'string' then string_dimension=dim
    print,'read ',name,' dimension ',dim
endfor
for ivar=0,result.nvars-1 do begin
    result=ncdf_varinq(ncid,ivar)
    ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
    if result.name eq 'AIM_ORBIT_NUMBER' then AIM_ORBIT_NUMBER=data
    if result.name eq 'VERSION' then VERSION=data
    if result.name eq 'UT_DATE' then UT_DATE=data
    if result.name eq 'CAMERA_ID' then CAMERA_ID=data
    if result.name eq 'NROWS' then NROWS=data
    if result.name eq 'NCOLS' then NCOLS=data
    if result.name eq 'OBSERVATION_TYPE' then OBSERVATION_TYPE=data
    if result.name eq 'UT_TIME' then UT_TIME=data
    if result.name eq 'LA_TIME' then LA_TIME=data
    if result.name eq 'IMAGE_TLM_TIMESTAMP' then IMAGE_TLM_TIMESTAMP=data
    if result.name eq 'SHUTTER_OPEN_TIMESTAMP' then SHUTTER_OPEN_TIMESTAMP=data
    if result.name eq 'SHUTTER_CLOSE_TIMESTAMP' then SHUTTER_CLOSE_TIMESTAMP=data
    if result.name eq 'INTEGRATION_TIME' then INTEGRATION_TIME=data
    if result.name eq 'MID_IMAGE_TIMESTAMP' then MID_IMAGE_TIMESTAMP=data
    if result.name eq 'ALBEDO' then ALBEDO=data
    if result.name eq 'ALBEDO_UNCERTAINTY' then ALBEDO_UNCERTAINTY=data
    if result.name eq 'COMMON_VOLUME_FLAG' then COMMON_VOLUME_FLAG=data
    if result.name eq 'COMMON_VOLUME_MAP' then COMMON_VOLUME_MAP=data
    if result.name eq 'QUALITY_FLAGS' then QUALITY_FLAGS=data
    if result.name eq 'SC_LATITUDE' then SC_LATITUDE=data
    if result.name eq 'SC_LONGITUDE' then SC_LONGITUDE=data
    if result.name eq 'SC_ALTITUDE' then SC_ALTITUDE=data
    if result.name eq 'ALONGTRACK_UNCERTAINTY' then ALONGTRACK_UNCERTAINTY=data
    if result.name eq 'CROSSTRACK_UNCERTAINTY' then CROSSTRACK_UNCERTAINTY=data
    if result.name eq 'RADIAL_UNCERTAINTY' then RADIAL_UNCERTAINTY=data
    if result.name eq 'TRACK_AZIMUTH' then TRACK_AZIMUTH=data
    if result.name eq 'GEOMETRY_DOWNSAMPLE' then GEOMETRY_DOWNSAMPLE=data
    if result.name eq 'LATITUDE' then LATITUDE=data
    if result.name eq 'LONGITUDE' then LONGITUDE=data
    if result.name eq 'SCATTERING_ANGLE' then SCATTERING_ANGLE=data
    if result.name eq 'VIEWING_ANGLE' then VIEWING_ANGLE=data
    if result.name eq 'ZENITH_ANGLE' then ZENITH_ANGLE=data
    print,'read ',result.name,' variable'
endfor	; loop over variables
;
; write CIPS data in IDL save format
;
save,file=ofile,/variables

endfor	; loop over files
end
