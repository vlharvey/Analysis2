;
; convert CIPS Level 3c IDL save files to NetCDF
;
; VLH 20180821
;
; Input files:
; cips_3c_north_18_v05.20_r02_all.sav
; cips_3c_north_18_v05.20_r02_cld.sav
; cips_3c_north_18_v05.20_r02_nocld.sav
;
; get Today's date
;
spawn,'date',today
;
; Define Input and Output directories
;
input_path='/atmos/harvey/CIPS_data/Datfiles/Level_3c_Summary/'
output_path='/atmos/harvey/CIPS_data/Datfiles/Level_3c_Summary/'
;
; 3-dimensional CIPS 3c files to convert
;
ifiles=[$
'cips_3c_north_18_v05.20_r02_all.sav',$
'cips_3c_north_18_v05.20_r02_cld.sav',$
'cips_3c_north_18_v05.20_r02_nocld.sav']
nfile=n_elements(ifiles)
;
; loop over files
;
for ii=0L,nfile-1L do begin
;
; restore file (contents are identical for all, cld, nocld)
; ALB             FLOAT     = Array[35, 1459, 120]
; ALB_AIR         FLOAT     = Array[35, 1459, 120]
; ALB_AIR_STD     FLOAT     = Array[35, 1459, 120]
; ALB_STD         FLOAT     = Array[35, 1459, 120]
; DATE            LONG      = Array[1459]
; IWC             FLOAT     = Array[35, 1459, 120]
; IWC_AIR         FLOAT     = Array[35, 1459, 120]
; IWC_AIR_STD     FLOAT     = Array[35, 1459, 120]
; IWC_STD         FLOAT     = Array[35, 1459, 120]
; LATHI           INT       = Array[120]
; LATLO           INT       = Array[120]
; LON             FLOAT     = Array[35, 1459, 120]
; LTIME           FLOAT     = Array[35, 1459, 120]
; NBIN            INT       =      120
; NREV            LONG      =         1459
; NTHRESH         INT       =       35
; NUM_CLD         INT       = Array[35, 1459, 120]
; NUM_OBS         INT       = Array[35, 1459, 120]
; RAD             FLOAT     = Array[35, 1459, 120]
; RAD_STD         FLOAT     = Array[35, 1459, 120]
; REV             LONG      = Array[1459]
; SZA             FLOAT     = Array[35, 1459, 120]
; THRESHOLD       FLOAT     = Array[35]
; UT              FLOAT     = Array[35, 1459, 120]
;
    restore,input_path+ifiles(ii)
;
; build output filename
;
    result=strsplit(ifiles(ii),'.',/extract)
    ofile=output_path+result(0)+'.'+result(1)+'.nc'
;
; create NetCDF file and describe global attributes and dimensions
;
    ncid = ncdf_create(ofile,/clobber)
    ncdf_attput, ncid, 'Description', 'CIPS Level 3c data with NTHRESH Gary thresholds', /global
    ncdf_attput, ncid, 'Author',      'V. Lynn Harvey, file was created on '+today(0)+' using convert_cips_3c_save2nc.pro', /global
    ncdf_attput, ncid, 'Comment',     'Input files used were cips_3c v05.20_r02 IDL save format', /global

    ncdf_attput, ncid, 'NTHRESH', 'Number of albedo threshold values (fixed at 35)', /global
    ncdf_attput, ncid, 'NREV', 'Total number of orbits in the season', /global
    ncdf_attput, ncid, 'NBIN', 'Number of latitude bins (fixed at 120)', /global
;
; define dimensions
;
    ntid = ncdf_dimdef(ncid, 'NTHRESH', NTHRESH) 
    nrid = ncdf_dimdef(ncid, 'NREV', NREV)
    nbid = ncdf_dimdef(ncid, 'NBIN', NBIN)
;
; define variable dimensions and type
;
    vid = ncdf_vardef(ncid, 'THRESHOLD'  , ntid, /FLOAT)
    vid = ncdf_vardef(ncid, 'DATE'       , nrid, /LONG)
    vid = ncdf_vardef(ncid, 'REV'        , nrid, /LONG)
    vid = ncdf_vardef(ncid, 'LATHI'      , nbid, /SHORT)
    vid = ncdf_vardef(ncid, 'LATLO'      , nbid, /SHORT)
    vid = ncdf_vardef(ncid, 'UT'         , [ntid,nrid,nbid], /FLOAT)
    vid = ncdf_vardef(ncid, 'LTIME'      , [ntid,nrid,nbid], /FLOAT)
    vid = ncdf_vardef(ncid, 'LON'        , [ntid,nrid,nbid], /FLOAT)
    vid = ncdf_vardef(ncid, 'SZA'        , [ntid,nrid,nbid], /FLOAT)
    vid = ncdf_vardef(ncid, 'NUM_CLD'    , [ntid,nrid,nbid], /FLOAT)
    vid = ncdf_vardef(ncid, 'NUM_OBS'    , [ntid,nrid,nbid], /FLOAT)
    vid = ncdf_vardef(ncid, 'RAD'        , [ntid,nrid,nbid], /FLOAT)
    vid = ncdf_vardef(ncid, 'RAD_STD'    , [ntid,nrid,nbid], /FLOAT)
    vid = ncdf_vardef(ncid, 'ALB'        , [ntid,nrid,nbid], /FLOAT)
    vid = ncdf_vardef(ncid, 'ALB_STD'    , [ntid,nrid,nbid], /FLOAT)
    vid = ncdf_vardef(ncid, 'IWC'        , [ntid,nrid,nbid], /FLOAT)
    vid = ncdf_vardef(ncid, 'IWC_STD'    , [ntid,nrid,nbid], /FLOAT)
    vid = ncdf_vardef(ncid, 'ALB_AIR'    , [ntid,nrid,nbid], /FLOAT)
    vid = ncdf_vardef(ncid, 'ALB_AIR_STD', [ntid,nrid,nbid], /FLOAT)
    vid = ncdf_vardef(ncid, 'IWC_AIR'    , [ntid,nrid,nbid], /FLOAT)
    vid = ncdf_vardef(ncid, 'IWC_AIR_STD', [ntid,nrid,nbid], /FLOAT)
;
; define variable names and units
;
    ncdf_attput, ncid, 'THRESHOLD',   'long_name', 'Albedo threshold'                                & ncdf_attput, ncid, 'THRESHOLD', 'units', '10-6 sr-1'
    ncdf_attput, ncid, 'DATE',        'long_name', 'DATE[YYYYMMDD]'                                  & ncdf_attput, ncid, 'DATE', 'units', 'YYYYMMDD'
    ncdf_attput, ncid, 'REV',         'long_name', 'AIM orbit number'                                & ncdf_attput, ncid, 'REV', 'units', 'N/A'
    ncdf_attput, ncid, 'LATHI',       'long_name', 'Upper latitude of bin'                           & ncdf_attput, ncid, 'LATHI', 'units', 'Degrees'
    ncdf_attput, ncid, 'LATLO',       'long_name', 'Lower latitude of bin'                           & ncdf_attput, ncid, 'LATLO', 'units', 'Degrees'
    ncdf_attput, ncid, 'UT',          'long_name', 'Mean UT time in bin'                             & ncdf_attput, ncid, 'UT', 'units', 'Hours'
    ncdf_attput, ncid, 'LTIME',       'long_name', 'Mean local time in bin'                          & ncdf_attput, ncid, 'LTIME', 'units', 'Hours'
    ncdf_attput, ncid, 'LON',         'long_name', 'Mean longitude in bin'                           & ncdf_attput, ncid, 'LON', 'units', 'Degrees'
    ncdf_attput, ncid, 'SZA',         'long_name', 'Mean solar zenith angle in bin'                  & ncdf_attput, ncid, 'SZA', 'units', 'Degrees'
    ncdf_attput, ncid, 'NUM_CLD',     'long_name', 'Number of pixels in bin containing clouds'       & ncdf_attput, ncid, 'NUM_CLD', 'units', 'N/A'
    ncdf_attput, ncid, 'NUM_OBS',     'long_name', 'Total number of data points in bin'              & ncdf_attput, ncid, 'NUM_OBS', 'units', 'N/A'
    ncdf_attput, ncid, 'RAD',         'long_name', 'Mean particle radius in bin'                     & ncdf_attput, ncid, 'RAD', 'units', 'nm'
    ncdf_attput, ncid, 'RAD_STD',     'long_name', 'Particle radius standard deviation in bin'       & ncdf_attput, ncid, 'RAD_STD', 'units', 'nm'
    ncdf_attput, ncid, 'ALB',         'long_name', 'Mean cloud albedo in bin'                        & ncdf_attput, ncid, 'ALB', 'units', '10-6 sr-1'
    ncdf_attput, ncid, 'ALB_STD',     'long_name', 'Cloud albedo standard deviation in bin'          & ncdf_attput, ncid, 'ALB_STD', 'units', '10-6 sr-1'
    ncdf_attput, ncid, 'IWC',         'long_name', 'Mean ice water content in bin'                   & ncdf_attput, ncid, 'IWC', 'units', 'g km-2'
    ncdf_attput, ncid, 'IWC_STD',     'long_name', 'Ice water content standard deviation in bin'     & ncdf_attput, ncid, 'IWC_STD', 'units', 'ug m-2'
    ncdf_attput, ncid, 'ALB_AIR',     'long_name', 'Mean AIR cloud albedo in bin'                    & ncdf_attput, ncid, 'ALB_AIR', 'units', '10-6 sr-1'
    ncdf_attput, ncid, 'ALB_AIR_STD', 'long_name', 'AIR cloud albedo standard deviation in bin'      & ncdf_attput, ncid, 'ALB_AIR_STD', 'units', '10-6 sr-1'
    ncdf_attput, ncid, 'IWC_AIR',     'long_name', 'Mean AIR ice water content in bin'               & ncdf_attput, ncid, 'IWC_AIR', 'units', 'g km-2'
    ncdf_attput, ncid, 'IWC_AIR_STD', 'long_name', 'AIR ice water content standard deviation in bin' & ncdf_attput, ncid, 'IWC_AIR_STD', 'units', 'ug m-2'

    ncdf_control, ncid, /endef
;
; write NetCDF file
;
    ncdf_varput, ncid, 'THRESHOLD',   THRESHOLD,   count=[NTHRESH]
    ncdf_varput, ncid, 'DATE',        DATE,        count=[NREV]
    ncdf_varput, ncid, 'REV',         REV,         count=[NREV]
    ncdf_varput, ncid, 'LATHI',       LATHI,       count=[NBIN]
    ncdf_varput, ncid, 'LATLO',       LATLO,       count=[NBIN]
    ncdf_varput, ncid, 'UT',          UT,          count=[NTHRESH,NREV,NBIN]
    ncdf_varput, ncid, 'LTIME',       LTIME,       count=[NTHRESH,NREV,NBIN]
    ncdf_varput, ncid, 'LON',         LON,         count=[NTHRESH,NREV,NBIN]
    ncdf_varput, ncid, 'SZA',         SZA,         count=[NTHRESH,NREV,NBIN]
    ncdf_varput, ncid, 'NUM_CLD',     NUM_CLD,     count=[NTHRESH,NREV,NBIN]
    ncdf_varput, ncid, 'NUM_OBS',     NUM_OBS,     count=[NTHRESH,NREV,NBIN]
    ncdf_varput, ncid, 'RAD',         RAD,         count=[NTHRESH,NREV,NBIN]
    ncdf_varput, ncid, 'RAD_STD',     RAD_STD,     count=[NTHRESH,NREV,NBIN]
    ncdf_varput, ncid, 'ALB',         ALB,         count=[NTHRESH,NREV,NBIN]
    ncdf_varput, ncid, 'ALB_STD',     ALB_STD,     count=[NTHRESH,NREV,NBIN]
    ncdf_varput, ncid, 'IWC',         IWC,         count=[NTHRESH,NREV,NBIN]
    ncdf_varput, ncid, 'IWC_STD',     IWC_STD,     count=[NTHRESH,NREV,NBIN]
    ncdf_varput, ncid, 'ALB_AIR',     ALB_AIR,     count=[NTHRESH,NREV,NBIN]
    ncdf_varput, ncid, 'ALB_AIR_STD', ALB_AIR_STD, count=[NTHRESH,NREV,NBIN]
    ncdf_varput, ncid, 'IWC_AIR',     IWC_AIR,     count=[NTHRESH,NREV,NBIN]
    ncdf_varput, ncid, 'IWC_AIR_STD', IWC_AIR_STD, count=[NTHRESH,NREV,NBIN]
;
; close NetCDF file
;
    ncdf_close, ncid

endfor	; end loop over input files
end
