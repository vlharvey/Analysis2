;
; check contents of SABER sathist file
;
output_path='/Users/harvey/Brakebusch/'		;C:\Users\jope9587\Desktop\Sathist\'
;
; open NetCDF file 
;
ncid=ncdf_open(output_path+'sathist_04_SABER_v2_20020125-20150301_c20150502.nc')
result=ncdf_inquire(ncid)
for ivar=0,result.nvars-1 do begin
    result=ncdf_varinq(ncid,ivar)
    ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
    if result.name eq 'date' then date=data
    if result.name eq 'time' then time=data
    if result.name eq 'latitude' then latitude=data
    if result.name eq 'longitude' then longitude=data
    if result.name eq 'orbit_num' then orbit_num=data
    if result.name eq 'prof_num' then prof_num=data
endfor
;
; check that orbit number is monotonically increasing with adjusted date/time
;
nevent=n_elements(date)
n0=findgen(nevent)
n1=1+findgen(nevent)
sdate=strcompress(date,/r)
y=long(strmid(sdate,0,4))
m=long(strmid(sdate,4,2))
d=long(strmid(sdate,6,2))
jday=julday(m,d,y,24.*time/86400.)

if min(orbit_num(n1)-orbit_num(n0)) lt 0 then stop,'check orbit numbers not monotonically increasing'
if min(jday(n1)-jday(n0)) lt 0 then stop,'check times not monotonically increasing'
end
