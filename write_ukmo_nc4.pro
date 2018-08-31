pro write_ukmo_nc4,ofile,nlg,nlat,nth,alon,alat,thlev,mark

; Create netCDF file and erase existing netCDF file if it exists
ncid = ncdf_create(ofile,/CLOBBER)
latdimid=ncdf_dimdef(ncid, 'number_of_latitudes' , nlat)
londimid=ncdf_dimdef(ncid, 'number_of_longitudes',  nlg)
levdimid=ncdf_dimdef(ncid, 'number_of_levels'    ,  nth)
lonsid = ncdf_vardef(ncid, 'longitudes',  londimid)
latsid = ncdf_vardef(ncid, 'latitudes' ,  latdimid)
levsid = ncdf_vardef(ncid, 'th_levels' ,  levdimid)
mksid  = ncdf_vardef(ncid, 'mark'      , [latdimid,londimid,levdimid])
ncdf_control,ncid,/ENDEF
ncdf_varput, ncid, lonsid, alon , COUNT=[nlg]
ncdf_varput, ncid, latsid, alat , COUNT=[nlat]
ncdf_varput, ncid, levsid, thlev, COUNT=[nth]
ncdf_varput, ncid, mksid , mark , COUNT=[nlat,nlg,nth]
ncdf_close,ncid
return
end
