   pro sofie_read_l2_netcdf,file,idump,$
       dp_version,nevent,nz,nwave,wave,miss,orbit,event,mode,date,time83,lat83,lon83,$
       z,p,t,t_pr,h2o,h2o_pr,o3,o3_pr,ch4,ch4_pr,no,no_pr,co2,co2_pr,ext,ext_pr,ext1037dv,ext1037dv_pr,$
       LOS_bearing,FOV_lock,       ext2_unc,ext3_unc,ext4_unc,extdv_unc

;-----------------------------------------------------------------------------------------------------------
; Routine reads a SOFIE level2 NetCDF file and returns named variables.  
;
; Input:
;
;   file.......path + name of the SOFIE L2 netcdf file
;   idump......1 = echo the variable names,  0 = be quiet
;
; Output:
;
;   Scalar header info:
;
;   dp_version..SOFIE data product version, string, scalar
;   nevent......number of events in this file, scalar
;   nz..........number of altitude points in profiles, scalar
;   nwave.......number of SOFIE wavelengths, scalar
;   miss........number that represents missing data (-1e24)
;
;   Vectors with 1 value per event:
;
;   orbit........AIM orbit number
;   event........SOFIE event number
;   mode.........occultation mode, 0 = sunrise, 1 = sunset
;   date.........date as yyyyddd
;   time83.......time at 83 km tangent point, seconds since UNIX epoch 
;   lat83........latitude (deg, -90 to 90), @ 83 km tangent point
;   lon83........longitude (deg E, 0 to 360), @ 83 km tangent point
;   LOS_bearing..line of sight bearing to Sun (deg E)
;   FOV_lock.....Field of view lockdown angle relative to the top edge of the sun (arcmin)
;
;   The common altitude grid for all profile data:
;
;   z.........altitude vector (km), dblarr(nz)
;
;   Profile data, all are dblarr(nz,nevent):
;
;   p.........pressure (hPa)
;   t.........temperature (K)
;   t_pr......temperature precision (K)
;   h2o.......water vapor mixing ratio (ppmv)
;   h2o_pr....water vapor precision (ppmv)
;   o3........ozone mixing ratio (ppmv)
;   o3_pr.....ozone precision (ppmv)
;   no........nitric oxide mixing ratio (ppbv)
;   no_pr.....nitric oxide precision (ppbv)
;   ch4.......methane mixing ratio (ppmv)
;   ch4_pr....methane precision (ppmv)
;   co2.......carbon dioxide mixing ratio (ppmv)
;   co2_pr....carbon dioxide precision (ppmv)
;
;   below are the aerosol extinction profile data:
;
;   wave..........SOFIE band center wavelengths (microns), fltarr(nwave)
;   ext...........aerosol extinction (1/km), dblarr(nz,nevent,nwave), the radiometer ("V") rerievals
;   ext_pr........aerosol extinction precision (1/km), dblarr(nz,nevent,nwave)
;   ext1037dv.....aerosol extinction from chan. 2 dV, 1.037 micron wavelength (1/km), dblarr(nz,nevent)
;   ext1037dv_pr..aerosol extinction precission from chan. 2 dV, 1.037 micron wavelength (1/km), dblarr(nz,nevent)
;
; Source: Mark Hervig, GATS Inc.
;
; Revision:  Mar 16, 2008 (MH)
;            Dec 12, 2008 (MH), added read statements for los_bearing, fov_lock
;            Dec 29, 2008 (MH), added more read statements for precision variables
;
;----------------------------------------------------------------------------------------------------------------

;- miscelaneous

   miss = -1d24  ; missing data value
   
;- open the file

   id   = ncdf_open(file,/nowrite)
   glob = ncdf_inquire(id)
   
;- echo the file contents

   glob = ncdf_inquire(id)
     
   varname = strarr(glob.nvars)   ; this will contain the variable names
     
   for i = 0,glob.nvars-1 do begin     
     info = ncdf_varinq(id, i)       
     varname(i) = info.name       
     if (idump eq 1) then  print,info.name       
   endfor
      
;- get the global attributes

   ncdf_diminq,id,'altitude',name,nz   ; number of points in profiles
   
   ncdf_attget,id,/global,'Title',out           &   title       = string(out)
   ncdf_attget,id,/global,'DP_Type',out         &   dp_type     = string(out)
   ncdf_attget,id,/global,'Source',out          &   source      = string(out)
   ncdf_attget,id,/global,'Mission',out         &   mission     = string(out)
   ncdf_attget,id,/global,'DP_Version',out      &   dp_version  = string(out)
   ncdf_attget,id,/global,'PF_Version',out      &   pf_version  = string(out)
   ncdf_attget,id,/global,'SW_Version',out      &   sw_version  = string(out)
   ncdf_attget,id,/global,'SW_Name',out         &   sw_name     = string(out)
   ncdf_attget,id,/global,'Calib_Version',out   &   cal_version = string(out)
   ncdf_attget,id,/global,'Description',out     &   description = string(out)
   ncdf_attget,id,/global,'History',out         &   history     = string(out)
   ncdf_attget,id,/global,'Gen_Date',out        &   gen_date    = string(out)

   if (idump eq 1) then begin   
     print,'Title:         ',title
     print,'DP_Type:       ',dp_type
     print,'Source:        ',source
     print,'DP_Version:    ',dp_version
     print,'PF_Version:    ',pf_version
     print,'SW_Version:    ',sw_version
     print,'SW_Name:       ',sw_name
     print,'Calib_Version: ',cal_version
     print,'Gen_Date:      ',gen_date
   endif
   
;- Read the header data

   ncdf_varget,id,'event'         , event
   ncdf_varget,id,'Orbit'         , orbit
   ncdf_varget,id,'Date'          , date
   ncdf_varget,id,'Mode'          , mode   
   ncdf_varget,id,'Time_83km'     , time83
   ncdf_varget,id,'Longitude_83km', lon83
   ncdf_varget,id,'Latitude_83km' , lat83
   
   LOS_bearing = lat83 * 0. + miss
   FOV_lock    = lat83 * 0. + miss
   
   k=where(varname eq 'LOS_Bearing',nk) & if nk gt 0 then ncdf_varget,id,'LOS_Bearing',LOS_bearing
   k=where(varname eq 'FOV_Lock'   ,nk) & if nk gt 0 then ncdf_varget,id,'FOV_Lock'   ,FOV_lock
   
   nevent = n_elements(event)
   
;- Read the profile data
      
   ncdf_varget,id,'Altitude'   , z
   ncdf_varget,id,'Pressure'   , p
   ncdf_varget,id,'Temperature', t
   
   ncdf_varget,id,'H2O_vmr', h2o   &  h2o = h2o *1d6
   ncdf_varget,id,'O3_vmr' , o3    &  o3  = o3  *1d6
   ncdf_varget,id,'CH4_vmr', ch4   &  ch4 = ch4 *1d6
   
   no = t * 0. + miss
   k = where(varname eq 'NO_vmr',nk) 
   if nk gt 0 then begin
     ncdf_varget,id,'NO_vmr' , no    &  no  = no  *1d9 ; ppbv
   endif
   
   co2 = t * 0. + miss
   k = where(varname eq 'CO2_vmr',nk) 
   if nk gt 0 then begin
     ncdf_varget,id,'CO2_vmr', co2   &  co2 = co2 *1d6 ; ppmv
   endif

   t_pr = t * 0. + miss
   k = where(varname eq 'Temperature_precision',nk) 
   if nk gt 0 then begin
     ncdf_varget,id,'Temperature_precision' , t_pr
   endif
      
   h2o_pr = t * 0. + miss
   k = where(varname eq 'H2O_vmr_precision',nk) 
   if nk gt 0 then begin
    ncdf_varget,id,'H2O_vmr_precision' , h2o_pr   &  h2o_pr = h2o_pr *1d6 ; ppmv
   endif
   
   o3_pr = t * 0. + miss
   k = where(varname eq 'O3_vmr_precision',nk) 
   if nk gt 0 then begin
    ncdf_varget,id,'O3_vmr_precision'  , o3_pr    &  o3_pr  = o3_pr  *1d6 ; ppmv
   endif
   
   ch4_pr = t * 0. + miss
   k = where(varname eq 'CH4_vmr_precision',nk) 
   if nk gt 0 then begin
    ncdf_varget,id,'CH4_vmr_precision' , ch4_pr   &  ch4_pr = ch4_pr *1d6 ; ppmv
   endif
   
   no_pr = t * 0. + miss
   k = where(varname eq 'NO_vmr_precision',nk) 
   if nk gt 0 then begin
     ncdf_varget,id,'NO_vmr_precision'  , no_pr    &  no_pr  = no_pr  *1d9 ; ppbv
   endif
   
   co2_pr = t * 0. + miss
   k = where(varname eq 'CO2_vmr_precision',nk) 
   if nk gt 0 then begin
     ncdf_varget,id,'CO2_vmr_precision' , co2_pr   &  co2_pr = co2_pr *1d6 ; ppmv
   endif
   
   ncdf_varget,id,'Extinction_0330' , e0330
   ncdf_varget,id,'Extinction_0867' , e0867
   ncdf_varget,id,'Extinction_1037' , e1037
   ncdf_varget,id,'Extinction_2462' , e2462
   ncdf_varget,id,'Extinction_2939' , e2939
   ncdf_varget,id,'Extinction_3064' , e3064
   ncdf_varget,id,'Extinction_3186' , e3186
   ncdf_varget,id,'Extinction_3479' , e3479
   ncdf_varget,id,'Extinction_4646' , e4646
   ncdf_varget,id,'Extinction_5006' , e5006
   
   k=where(varname eq 'Extinction_0330_precision',nk2 ) & if nk2  gt 0 then ncdf_varget,id,'Extinction_0330_precision' , e0330_pr
   k=where(varname eq 'Extinction_0867_precision',nk3 ) & if nk3  gt 0 then ncdf_varget,id,'Extinction_0867_precision' , e0867_pr
   k=where(varname eq 'Extinction_1037_precision',nk4 ) & if nk4  gt 0 then ncdf_varget,id,'Extinction_1037_precision' , e1037_pr
   k=where(varname eq 'Extinction_2462_precision',nk5 ) & if nk5  gt 0 then ncdf_varget,id,'Extinction_2462_precision' , e2462_pr
   k=where(varname eq 'Extinction_2939_precision',nk8 ) & if nk8  gt 0 then ncdf_varget,id,'Extinction_2939_precision' , e2939_pr
   k=where(varname eq 'Extinction_3064_precision',nk9 ) & if nk9  gt 0 then ncdf_varget,id,'Extinction_3064_precision' , e3064_pr
   k=where(varname eq 'Extinction_3186_precision',nk10) & if nk10 gt 0 then ncdf_varget,id,'Extinction_3186_precision' , e3186_pr
   k=where(varname eq 'Extinction_3479_precision',nk12) & if nk12 gt 0 then ncdf_varget,id,'Extinction_3479_precision' , e3479_pr
   k=where(varname eq 'Extinction_4646_precision',nk14) & if nk14 gt 0 then ncdf_varget,id,'Extinction_4646_precision' , e4646_pr
   k=where(varname eq 'Extinction_5006_precision',nk15) & if nk15 gt 0 then ncdf_varget,id,'Extinction_5006_precision' , e5006_pr
   
   k=where(varname eq 'Extinction_1037dv',nk) & if (nk gt 0) then ncdf_varget,id,'Extinction_1037dv', ext1037dv 
   
   ext1037dv_pr = t * 0. + miss
   k=where(varname eq 'Extinction_1037dv_precision',nk) & if (nk gt 0) then ncdf_varget,id,'Extinction_1037dv_precision', ext1037dv_pr
      
;- Order the extinctions into one array

   nwave = 16L
   
;  band     1     2     3     4     5     6     7     8     9     10    11    12    13    14    15    16
   wave  = [0.292,0.330,0.867,1.037,2.462,2.618,2.785,2.939,3.064,3.186,3.384,3.479,4.324,4.646,5.006,5.316] ; SOFIE band center wavelengths, microns

   ext    = dblarr(nz,nevent,nwave)
   ext_pr = dblarr(nz,nevent,nwave)
   
  ;ext(*,*, 0) = 
   ext(*,*, 1) = e0330
   ext(*,*, 2) = e0867
   ext(*,*, 3) = e1037
   ext(*,*, 4) = e2462
  ;ext(*,*, 5) = 
  ;ext(*,*, 6) = 
   ext(*,*, 7) = e2939
   ext(*,*, 8) = e3064
   ext(*,*, 9) = e3186
  ;ext(*,*,10) = 
   ext(*,*,11) = e3479
  ;ext(*,*,12) = 
   ext(*,*,13) = e4646
   ext(*,*,14) = e5006
  ;ext(*,*,15) = 

   if nk2  gt 0 then ext_pr(*,*, 1) = e0330_pr
   if nk3  gt 0 then ext_pr(*,*, 2) = e0867_pr
   if nk4  gt 0 then ext_pr(*,*, 3) = e1037_pr
   if nk5  gt 0 then ext_pr(*,*, 4) = e2462_pr
   if nk8  gt 0 then ext_pr(*,*, 7) = e2939_pr
   if nk9  gt 0 then ext_pr(*,*, 8) = e3064_pr
   if nk10 gt 0 then ext_pr(*,*, 9) = e3186_pr
   if nk12 gt 0 then ext_pr(*,*,11) = e3479_pr
   if nk14 gt 0 then ext_pr(*,*,13) = e4646_pr
   if nk15 gt 0 then ext_pr(*,*,14) = e5006_pr

;- Read some in-house research variables,  these are not in public release data

   ext2_unc = t * 0. + miss  ; uncorrected extinctions
   ext3_unc = t * 0. + miss
   ext4_unc = t * 0. + miss
   extdv_unc = t * 0. + miss
   
   k = where(varname eq 'Band2_uncorrected',nk) & if nk gt 0 then ncdf_varget,id,'Band2_uncorrected' , ext2_unc   
   k = where(varname eq 'Band3_uncorrected',nk) & if nk gt 0 then ncdf_varget,id,'Band3_uncorrected' , ext3_unc   
   k = where(varname eq 'Band4_uncorrected',nk) & if nk gt 0 then ncdf_varget,id,'Band4_uncorrected' , ext4_unc   
   k = where(varname eq 'Ch2_DV_uncorrected',nk) & if nk gt 0 then ncdf_varget,id,'Ch2_DV_uncorrected' , extdv_unc
   
;- Close the file

   ncdf_close,id  

;- Done

   return
   end