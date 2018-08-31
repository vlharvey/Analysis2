pro write_waccm3_nc3_bardeen,ofile,nlg,nlat,nth,alon,alat,thlev,$
          ipv,prs,u,v,qdf,mark,vp,sf

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
vptid  = ncdf_vardef(ncid, 'vel_pot'   , [latdimid,londimid,levdimid])
stfid  = ncdf_vardef(ncid, 'strm'      , [latdimid,londimid,levdimid])

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
ncdf_varput, ncid, vptid , vp   , COUNT=[nlat,nlg,nth]
ncdf_varput, ncid, stfid , sf   , COUNT=[nlat,nlg,nth]
ncdf_close,ncid
return
end
