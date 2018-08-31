;
; read CIPS RAA L2 scene and orbit data and check that Nan values are consistent among arrays
;
; /atmos/harvey/CIPS_data/Datfiles/RAA/cips_raa_2a_orbit_?????_YYYY-DOY_v01.00_r02_alb.nc (2A=Scene data: alb, ang, cat)
; /atmos/harvey/CIPS_data/Datfiles/RAA/cips_raa_2b_orbit_?????_YYYY-DOY_v01.00_r02_alb.nc (2B=Orbit data: alb, cat)
;
loadct,39
icolmax=byte(!p.color)
icolmax=fix(icolmax)
if icolmax eq 0 then icolmax=255
mcolor=icolmax
device,decompose=0
!p.background=icolmax
a=findgen(8)*(2*!pi/8.)
usersym,cos(a),sin(a),/fill
setplot='ps'
read,'setplot=',setplot
nxdim=750
nydim=750
xorig=[0.15]
yorig=[0.15]
xlen=0.8
ylen=0.8
cbaryoff=0.07
cbarydel=0.01
!NOERAS=-1
if setplot ne 'ps' then window,4,xsize=nxdim,ysize=nydim,retain=2,colors=162

dir='/atmos/harvey/CIPS_data/Datfiles/RAA/'
ifiles=file_search(dir+'cips_raa_2a*_alb.nc')

nfiles=n_elements(ifiles)
orbitnum_all=lonarr(nfiles)
fyr_all=fltarr(nfiles)
fracgood_all=fltarr(nfiles)

for ifile=0L,nfiles-1L do begin
    ncfile0=ifiles(ifile)

    ncid=ncdf_open(ncfile0)
    result=ncdf_inquire(ncid)
    nvars=result.nvars        ;# variables in the file
    for ivar=0,nvars-1 do begin
        result=ncdf_varinq(ncid,ivar) ;get the data name, type, dimensions
        ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
        if ( ( size( data, /n_dimensions )  EQ 1 ) && $
             ( size( data, /type ) EQ 1 ) ) then $
               data = string( data )
        if Execute(result.name + ' = data') eq 0 then $
           Print, ' "Execute" command failed -- are you in Virtual Machine mode?'
;       print,result.name
    endfor
    ncdf_close,ncid

index=strpos(ncfile0,'alb')
catfile=strmid(ncfile0,0,index)+'cat.nc'
angfile=strmid(ncfile0,0,index)+'ang.nc'

index=strpos(ncfile0,'orbit')

orbitnum_all(ifile)=long(strmid(ncfile0,index(0)+6,5))
if long(strmid(ncfile0,index(0)+12,4)) mod 4 ne 0 then fyr_all(ifile)=float(strmid(ncfile0,index(0)+12,4))+float(strmid(ncfile0,index(0)+17,3))/365.
if long(strmid(ncfile0,index(0)+12,4)) mod 4 eq 0 then fyr_all(ifile)=float(strmid(ncfile0,index(0)+12,4))+float(strmid(ncfile0,index(0)+17,3))/366.

    ncid=ncdf_open(catfile)
    result=ncdf_inquire(ncid)
    nvars=result.nvars        ;# variables in the file
    for ivar=0,nvars-1 do begin
        result=ncdf_varinq(ncid,ivar) ;get the data name, type, dimensions
        ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
        if ( ( size( data, /n_dimensions )  EQ 1 ) && $
             ( size( data, /type ) EQ 1 ) ) then $
               data = string( data )
        if Execute(result.name + ' = data') eq 0 then $
           Print, ' "Execute" command failed -- are you in Virtual Machine mode?'
;       print,result.name
    endfor
    ncdf_close,ncid

    ncid=ncdf_open(angfile)
    result=ncdf_inquire(ncid)
    nvars=result.nvars        ;# variables in the file
    for ivar=0,nvars-1 do begin
        result=ncdf_varinq(ncid,ivar) ;get the data name, type, dimensions
        ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
        if ( ( size( data, /n_dimensions )  EQ 1 ) && $
             ( size( data, /type ) EQ 1 ) ) then $
               data = string( data )
        if Execute(result.name + ' = data') eq 0 then $
           Print, ' "Execute" command failed -- are you in Virtual Machine mode?'
;       print,result.name
    endfor
    ncdf_close,ncid

;help,longitude,latitude,RAYLEIGH_ALBEDO,RAYLEIGH_ALBEDO_ANOMALY,RAYLEIGH_ALBEDO_ANOMALY_UNC,SCATTERING_ANGLE,VIEW_ANGLE,VIEW_ANGLE_DERIVATIVE,ZENITH_ANGLE,ZENITH_ANGLE_DERIVATIVE
result=size(longitude)
npoints=result(-1)

good1=where(finite(longitude) eq 1)
good2=where(finite(latitude) eq 1)
good3=where(finite(RAYLEIGH_ALBEDO) eq 1)
good4=where(finite(RAYLEIGH_ALBEDO_ANOMALY) eq 1)
good5=where(finite(RAYLEIGH_ALBEDO_ANOMALY_UNC) eq 1)
good6=where(finite(SCATTERING_ANGLE) eq 1)
good7=where(finite(VIEW_ANGLE) eq 1)
good8=where(finite(VIEW_ANGLE_DERIVATIVE) eq 1)
good9=where(finite(ZENITH_ANGLE) eq 1)
good10=where(finite(ZENITH_ANGLE_DERIVATIVE) eq 1)

print,ncfile0,100.*n_elements(good10)/float(npoints),' good data'
print,100.*n_elements(good1)/float(npoints),' good lons'

fracgood_all(ifile)=100.*n_elements(good10)/float(npoints)

if n_elements(good3) ne n_elements(good4) or $
   n_elements(good3) ne n_elements(good5) or $
   n_elements(good3) ne n_elements(good6) or $
   n_elements(good3) ne n_elements(good7) or $
   n_elements(good3) ne n_elements(good8) or $
   n_elements(good3) ne n_elements(good9) or $
   n_elements(good3) ne n_elements(good10) or $
   n_elements(good4) ne n_elements(good5) or $
   n_elements(good4) ne n_elements(good6) or $
   n_elements(good4) ne n_elements(good7) or $
   n_elements(good4) ne n_elements(good8) or $
   n_elements(good4) ne n_elements(good9) or $
   n_elements(good4) ne n_elements(good10) or $
   n_elements(good5) ne n_elements(good6) or $
   n_elements(good5) ne n_elements(good7) or $
   n_elements(good5) ne n_elements(good8) or $
   n_elements(good5) ne n_elements(good9) or $
   n_elements(good5) ne n_elements(good10) or $
   n_elements(good6) ne n_elements(good7) or $
   n_elements(good6) ne n_elements(good8) or $
   n_elements(good6) ne n_elements(good9) or $
   n_elements(good6) ne n_elements(good10) or $
   n_elements(good7) ne n_elements(good8) or $
   n_elements(good7) ne n_elements(good9) or $
   n_elements(good7) ne n_elements(good10) or $
   n_elements(good8) ne n_elements(good9) or $
   n_elements(good8) ne n_elements(good10) or $
   n_elements(good9) ne n_elements(good10) then stop,'different number of Nans in 2 or more arrays'

endfor

erase
loadct,39
plot,fyr_all,fracgood_all,psym=8,color=0,title=ncfile0

end
