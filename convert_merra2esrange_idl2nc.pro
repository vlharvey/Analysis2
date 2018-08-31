;
; convert merra_at_esrange_20150228.sav to netcdf
; ouput MERRA theta data at Esrange 67.8939° N, 21.1069° E
;
;longitude=slon
;latitude=slat
;potential_temperature_profile=th
;date=long(sdate)
;comment=strarr(11)
;comment(0)='Profiles based on MERRA daily averaged data'
;comment(1)='dates in YYYYMMDD'
;comment(2)='longitude in degrees east'
;comment(3)='latitude in degrees'
;comment(4)='potential_temperature = potential temperature profile (K)'
;comment(5)='vortex_marker_profiles = positive (negative) values in vortex (anticyclones)'
;comment(6)='pressure_profiles = Pressure profiles (hPa)'
;comment(7)='temperature_profiles = Temperature profiles (K)'
;comment(8)='altitude_profiles = Geometric Altitude profiles (km)'
;comment(9)='zonal_wind_profiles = Zonal Wind profiles (km)'
;comment(10)='meridional_wind_profiles = Meridional Wind profiles (km)'

pthout='/Volumes/Data/MERRA_data/Datfiles/'
sdate='20150228'
restore,pthout+'merra_at_esrange_'+sdate+'.sav'	;,dates,time,longitude,latitude,$
;     potential_temperature_profile,vortex_marker_profiles,pressure_profiles,$
;     temperature_profiles,altitude_profiles,zonal_wind_profiles,meridional_wind_profiles,comment
;print,'saved DMP '+sdate

ntime=n_elements(dates)
ntheta=n_elements(potential_temperature_profile)

;open NetCDF file for sequential writing
ncid = ncdf_create(pthout+'merra_at_esrange_1997-2015.nc',/clobber)
ncdf_attput, ncid, 'Comment', 'Profiles based on MERRA daily averaged data', /global
ncdf_attput, ncid, 'Author', 'V. Lynn Harvey, file was created on Mon Apr 13 06:58:39 MDT 2015', /global

tid = ncdf_dimdef(ncid, 'ntime', ntime)
zid = ncdf_dimdef(ncid, 'ntheta', ntheta)

vid = ncdf_vardef(ncid, 'dates', tid, /LONG)
vid = ncdf_vardef(ncid, 'latitude', /FLOAT)
vid = ncdf_vardef(ncid, 'longitude',/FLOAT)
vid = ncdf_vardef(ncid, 'POTENTIAL_TEMPERATURE_PROFILE', zid, /FLOAT)
vid = ncdf_vardef(ncid, 'ALTITUDE_PROFILES', [tid,zid] , /FLOAT)
vid = ncdf_vardef(ncid, 'PRESSURE_PROFILES', [tid,zid] , /FLOAT)
vid = ncdf_vardef(ncid, 'TEMPERATURE_PROFILES', [tid,zid] , /FLOAT)
vid = ncdf_vardef(ncid, 'VORTEX_MARKER_PROFILES', [tid,zid] , /FLOAT)
vid = ncdf_vardef(ncid, 'ZONAL_WIND_PROFILES', [tid,zid] , /FLOAT)
vid = ncdf_vardef(ncid, 'MERIDIONAL_WIND_PROFILES', [tid,zid] , /FLOAT)

ncdf_attput, ncid, 'dates', 'long_name', 'dates[yyyymmdd]' & ncdf_attput, ncid, 'dates', 'units', 'YYYYMMDD'
ncdf_attput, ncid, 'latitude', 'long_name', 'latitude' & ncdf_attput, ncid, 'latitude', 'units', 'degrees'
ncdf_attput, ncid, 'longitude', 'long_name', 'longitude' & ncdf_attput, ncid, 'longitude', 'units', 'degrees'
ncdf_attput, ncid, 'POTENTIAL_TEMPERATURE_PROFILE', 'long_name', 'Potential Temperature Profile' & ncdf_attput, ncid, 'POTENTIAL_TEMPERATURE_PROFILE', 'units', 'K'
ncdf_attput, ncid, 'ALTITUDE_PROFILES', 'long_name', 'Altitude Profiles' & ncdf_attput, ncid, 'ALTITUDE_PROFILES', 'units', 'km'
ncdf_attput, ncid, 'PRESSURE_PROFILES', 'long_name', 'Pressure Profiles' & ncdf_attput, ncid, 'PRESSURE_PROFILES', 'units', 'hPa'
ncdf_attput, ncid, 'TEMPERATURE_PROFILES', 'long_name', 'Temperature Profiles' & ncdf_attput, ncid, 'TEMPERATURE_PROFILES', 'units', 'K'
ncdf_attput, ncid, 'VORTEX_MARKER_PROFILES', 'long_name', 'Vortex Marker Profiles' & ncdf_attput, ncid, 'VORTEX_MARKER_PROFILES', 'units', '-1 in anticyclones; zero in neither; between 0 and 1 in vortex edge; +1 inside vortex'
ncdf_attput, ncid, 'ZONAL_WIND_PROFILES', 'long_name', 'Zonal Wind Profiles' & ncdf_attput, ncid, 'ZONAL_WIND_PROFILES', 'units', 'm/s'
ncdf_attput, ncid, 'MERIDIONAL_WIND_PROFILES', 'long_name', 'Meridional Wind Profiles' & ncdf_attput, ncid, 'MERIDIONAL_WIND_PROFILES', 'units', 'm/s'

ncdf_control, ncid, /endef

ncdf_varput, ncid, 'dates', dates, count=[ntime]
ncdf_varput, ncid, 'latitude', latitude, count=1	;[ntime,ntheta]
ncdf_varput, ncid, 'longitude', longitude, count=1	;[ntime,ntheta]
ncdf_varput, ncid, 'POTENTIAL_TEMPERATURE_PROFILE', POTENTIAL_TEMPERATURE_PROFILE, count=[ntheta]
ncdf_varput, ncid, 'ALTITUDE_PROFILES', ALTITUDE_PROFILES, count=[ntime,ntheta]
ncdf_varput, ncid, 'PRESSURE_PROFILES', PRESSURE_PROFILES, count=[ntime,ntheta]
ncdf_varput, ncid, 'TEMPERATURE_PROFILES', TEMPERATURE_PROFILES, count=[ntime,ntheta]
ncdf_varput, ncid, 'VORTEX_MARKER_PROFILES', VORTEX_MARKER_PROFILES, count=[ntime,ntheta]
ncdf_varput, ncid, 'ZONAL_WIND_PROFILES', ZONAL_WIND_PROFILES, count=[ntime,ntheta]
ncdf_varput, ncid, 'MERIDIONAL_WIND_PROFILES', MERIDIONAL_WIND_PROFILES, count=[ntime,ntheta]

ncdf_close, ncid




end
