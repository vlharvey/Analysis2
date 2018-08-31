pro write_waccm3_nc3,ofile,nlg,nlat,nth,alon,alat,thlev,$
          ipv,prs,u,v,qdf,mark,sf,o3,ch4,no2,h2o

; Create netCDF file and erase existing netCDF file if it exists
ncid = ncdf_create(ofile,/CLOBBER)
latdimid=ncdf_dimdef(ncid, 'number_of_latitudes' , nlat)
londimid=ncdf_dimdef(ncid, 'number_of_longitudes',  nlg)
levdimid=ncdf_dimdef(ncid, 'number_of_levels'    ,  nth)
lonsid = ncdf_vardef(ncid, 'longitudes',  londimid)
latsid = ncdf_vardef(ncid, 'latitudes' ,  latdimid)
levsid = ncdf_vardef(ncid, 'th_levels' ,  levdimid)
ipvid  = ncdf_vardef(ncid, 'ipv'       , [latdimid,londimid,levdimid])
prsid  = ncdf_vardef(ncid, 'press'     , [latdimid,londimid,levdimid])
uuuid  = ncdf_vardef(ncid, 'u_wind'    , [latdimid,londimid,levdimid])
vvvid  = ncdf_vardef(ncid, 'v_wind'    , [latdimid,londimid,levdimid])
qdfid  = ncdf_vardef(ncid, 'qdf'       , [latdimid,londimid,levdimid])
mksid  = ncdf_vardef(ncid, 'mark'      , [latdimid,londimid,levdimid])
stfid  = ncdf_vardef(ncid, 'strm'      , [latdimid,londimid,levdimid])
o3fid  = ncdf_vardef(ncid, 'ozone'     , [latdimid,londimid,levdimid])
ch4id  = ncdf_vardef(ncid, 'methane'   , [latdimid,londimid,levdimid])
no2id  = ncdf_vardef(ncid, 'no2'       , [latdimid,londimid,levdimid])
h2oid  = ncdf_vardef(ncid, 'h2o'       , [latdimid,londimid,levdimid])

ncdf_control,ncid,/ENDEF
ncdf_varput, ncid, lonsid, alon , COUNT=[nlg]
ncdf_varput, ncid, latsid, alat , COUNT=[nlat]
ncdf_varput, ncid, levsid, thlev, COUNT=[nth]
ncdf_varput, ncid, ipvid , ipv  , COUNT=[nlat,nlg,nth]
ncdf_varput, ncid, prsid , prs  , COUNT=[nlat,nlg,nth]
ncdf_varput, ncid, uuuid , u    , COUNT=[nlat,nlg,nth]
ncdf_varput, ncid, vvvid , v    , COUNT=[nlat,nlg,nth]
ncdf_varput, ncid, qdfid , qdf  , COUNT=[nlat,nlg,nth]
ncdf_varput, ncid, mksid , mark , COUNT=[nlat,nlg,nth]
ncdf_varput, ncid, stfid , sf   , COUNT=[nlat,nlg,nth]
ncdf_varput, ncid, o3fid , o3   , COUNT=[nlat,nlg,nth]
ncdf_varput, ncid, ch4id , ch4  , COUNT=[nlat,nlg,nth]
ncdf_varput, ncid, no2id , no2  , COUNT=[nlat,nlg,nth]
ncdf_varput, ncid, h2oid , h2o  , COUNT=[nlat,nlg,nth]
ncdf_close,ncid
return
end
