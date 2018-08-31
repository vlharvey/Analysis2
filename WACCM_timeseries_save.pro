@/Volumes/MacD68-1/france/idl_files/stddat			; Determines the number of days since Jan 1, 1956
@/Volumes/MacD68-1/france/idl_files/kgmt			; This function computes the Julian day number (GMT) from the
								;    day, month, and year information.
@/Volumes/MacD68-1/france/idl_files/ckday			; This routine changes the Julian day from 365(6 if leap yr)
								;    to 1 and increases the year, if necessary.
@/Volumes/MacD68-1/france/idl_files/kdate			; gives back kmn,kdy information from the Julian day #.
@/Volumes/MacD68-1/france/idl_files/rd_ukmo_nc3			; reads the data from nc3 files
@/Volumes/MacD68-1/france/idl_files/date2uars			; This code returns the UARS day given (jday,year) information.
@/Volumes/MacD68-1/france/idl_files/plotPosition		; defines positions for n number of plots
@/Volumes/MacD68-1/france/idl_files/rd_GEOS5_nc3



restore, '/Volumes/MacD68-1/france/WACCM_paper/Post_process/Post_process_1_height_temps_mark.sav'



;-----------Determine the range of days and set up the counter for the date----------------------

lstdy = 1L
lstmn = 1L
lstyr = 3L
read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
days = 14600L
;read,' Enter number of days  ', days
;  set endday equal to startday so for each plot, only one day is plotted, for each elements of 'days'
ledmn = lstmn
leddy = lstdy
ledyr = lstyr

doy = fltarr(days)
idoy = 1L

meandailyProfile = fltarr(n_elements(dates),200L)

mon=['jan_','feb_','mar_','apr_','may_','jun_',$
    'jul_','aug_','sep_','oct_','nov_','dec_']
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1991 then stop,'Year out of range '
;if ledyr lt 1991 then stop,'Year out of range '
stddat,lstmn,lstdy,lstyr,lstday
stddat,ledmn,leddy,ledyr,ledday
if ledday lt lstday then stop,' Wrong dates! '
iyr = lstyr
idy = lstdy
imn = lstmn
kgmt,imn,idy,iyr,iday
startday = julday(lstmn, lstdy, lstyr)
endday = julday(ledmn, leddy, ledyr)


dir='/Volumes/MacD68-2/france/data/WACCM_4/Datfiles/mee00fpl_FW2.cam2.h3.dyns.20'

djfdate = 0L

FOR iii = 0L, days - 1L DO BEGIN


sys1 = SYSTIME( /JULIAN)


if iii gt 0L then begin

  startday = startday + 1L
  endday = endday + 1L
  caldat, startday, lstmn, lstdy, lstyr
  lstyr = lstyr - 2000L
  ledmn = lstmn
  leddy = lstdy
  ledyr = lstyr
  print, lstyr, ledmn, leddy
  mon=['jan_','feb_','mar_','apr_','may_','jun_',$
       'jul_','aug_','sep_','oct_','nov_','dec_']
  if lstyr lt 91 then lstyr=lstyr+2000
  if ledyr lt 91 then ledyr=ledyr+2000
  if lstyr lt 1900 then lstyr=lstyr+1900
  if ledyr lt 1900 then ledyr=ledyr+1900
  if lstyr lt 1991 then stop,'Year out of range '
  ;if ledyr lt 1991 then stop,'Year out of range '
  stddat,lstmn,lstdy,lstyr,lstday
  stddat,ledmn,leddy,ledyr,ledday
  if ledday lt lstday then stop,' Wrong dates! '
  iyr = lstyr
  idy = lstdy
  imn = lstmn
  kgmt,imn,idy,iyr,iday

endif


  ;----------------------------------------------
  ; determine the date
  ckday,iday,iyr
  kdate,float(iday),iyr,imn,idy
  stddat, imn,idy,iyr,ndays
  if iyr lt 2000 then iyr1=iyr-1900
  if iyr ge 2000 then iyr1=iyr-2000
  syr   = string(FORMAT='(I2.2)',iyr1)
  smn   = string(FORMAT='(I2.2)',imn)
  sdy   = string(FORMAT='(I2.2)',idy)
  doy = string(Format='(I03)',iday)
  ifile = string(syr+smn+sdy)
  ;----------------------------------------------


nc = 0L
ifiles=file_search(dir+ifile+'_3D_dyn.nc3',count=nfile)
	if ifiles[0] eq '' then continue
    result=strsplit(ifiles(0),'.',/extract)
    result2=strsplit(result(4),'_',/extract)
    sdate=result2(0)
    iflag=0
;
; read daily file
;
    	ncfile0=ifiles(0)

       	ncid=ncdf_open(ncfile0)
    	result0=ncdf_inquire(ncid)
    	for idim=0,result0.ndims-1 do begin
        	ncdf_diminq,ncid,idim,name,dim
        	if name eq 'number_of_latitudes' then nr=dim
        	if name eq 'number_of_longitudes' then nc=dim
        	if name eq 'number_of_levels' then nth=dim
;       	print,'read ',name,' dimension ',dim
    endfor
    	for ivar=0,result0.nvars-1 do begin
        	result=ncdf_varinq(ncid,ivar)
        	ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
        	if result.name eq 'latitude' then alat=data
        	if result.name eq 'longitude' then alon=data
        	if result.name eq 'theta' then th=data
;        	if result.name eq 'IPV' then pv2=data
        	if result.name eq 'P' then p2=data
;        	if result.name eq 'U' then u2=data
;        	if result.name eq 'V' then v2=data
;        	if result.name eq 'QDF' then qdf2=data
;        	if result.name eq 'Q' then q2=data
        	if result.name eq 'GPH' then gph2=data
;        	if result.name eq 'TTGW' then ttgw2=data
;        	if result.name eq 'SF' then sf2=data
        	if result.name eq 'MARK' then marksf2=data
;        	print,ivar,result.name,min(data),max(data)
	endfor
    		ncdf_close,ncid


n_lat = n_elements(alat)
n_lon = n_elements(alon)
n_levels    = 200L
n_theta = nth
		if nc eq 0L then continue
		WACCMtemp=!values.f_nan*p2
zgrd=fltarr(nr,nc,nth)   ; 3-d geopotential height
altitude_grid = findgen(200L)+1L
altitude = altitude_grid
intWACCMtemps = fltarr(n_lat,n_lon,200L)
intWACCMmark = fltarr(n_lat,n_lon,200L)

k = 8.314
g=9.81
M = .02896
rtd=double(180./!pi)
dtr=1./rtd
ks=1.931853d-3
ecc=0.081819
gamma45=9.80
gph2 = gph2/1000.


        ; if file exists, then determine daily mean T and P

lat = alat
lon = alon
      
WACCMtemp     = fltarr(n_lat, n_lon, n_theta) * !values.f_nan	; averaged temperatures over n_dayAvg and longitude, fit to lat grid
zgrd     = fltarr(n_lat, n_lon, n_theta) * !values.f_nan

for ii=0L,nr-1L do begin

	; CONVERT GEOPOTENTIAL TO GEOMETRIC HEIGHT

	sin2=sin( (alat(ii)*dtr)^2.0 )
	numerator=1.0+ks*sin2
	denominator=sqrt( 1.0 - (ecc^2.0)*sin2 )
	gammas=gamma45*(numerator/denominator)
	r=6378.137D/(1.006803-(0.006706*sin2))
	for jj = 0L, nc - 1L do begin
		waccmtemp[ii,jj,*] = th[*]*( (p2[ii,jj,*]/1000.D)^(.286) )
	    	z=(r*reform(gph2(ii,jj,*)))/ ( (gammas/gamma45)*r - reform(gph2(ii,jj,*) ))
		zgrd[ii,jj,*] = z

		x =interpol(reverse(reform(waccmtemp[ii,jj,*])),reverse(reform(zgrd[ii,jj,*])),altitude_grid)
		intWACCMtemps[ii,jj,*] = x    
	endfor
endfor		


if iyr eq 2003L and imn eq 1L and idy eq 1L then begin
	strdates = strarr(120L)
	dateSeries = findgen(120L) + 1L
	tempSeries = fltarr(n_elements(dateSeries), n_elements(altitude))
	tday = 0L
	mode = 0L
endif




if imn eq 11L and idy eq 1L then begin

	smoothtemp = tempseries*!values.f_nan
	x = where(tempseries le 0. or tempseries gt 1000.,nx)
	if nx gt 0L then tempseries[x] = !values.f_nan

	for iiii = 4L, n_elements(strdates) - 5L do begin
		for iz = 0L, n_elements(altitude) - 1L do begin
			smoothtemp[iiii] = mean(tempseries[iiii-4L:iiii+4L,iz],/nan)
		endfor
	endfor

	save, filename = '/Volumes/MacD68-1/france/WACCM_paper/Post_process/timeseries_SH_20'+syr+'.sav', tempseries, dateseries, altitude, strdates, smoothtemp
	strdates = strarr(181L)
	dateSeries = findgen(181L) + 1L
	tempSeries = fltarr(n_elements(dateSeries), n_elements(altitude))
	tday = 0L
	mode = 0L
endif
if imn eq 5L and idy eq 1L then begin

	smoothtemp = tempseries*!values.f_nan
	x = where(tempseries le 0. or tempseries gt 1000.,nx)
	if nx gt 0L then tempseries[x] = !values.f_nan

	for iiii = 4L, n_elements(strdates) - 5L do begin
		for iz = 0L, n_elements(altitude) - 1L do begin
			smoothtemp[iiii] = mean(tempseries[iiii-4L:iiii+4L,iz],/nan)
		endfor
	endfor

	save, filename = '/Volumes/MacD68-1/france/WACCM_paper/Post_process/timeseries_NH_20'+syr+'.sav', tempseries, dateseries, altitude, strdates, smoothtemp
	strdates = strarr(184L)
	dateSeries = findgen(184L) + 1L
	tempSeries = fltarr(n_elements(dateSeries), n_elements(altitude))
	tday = 0L
	mode = 1L
endif

x = where(intwaccmtemps le 100. and intwaccmtemps gt 1500.,nx)
if nx gt 0L then intwaccmtemps[x] = !values.f_nan
if mode eq 0L then lat70 = where(lat ge 70L)
if mode eq 1L then lat70 = where(lat le -70L)
latweight = SQRT(COS(2*3.141592*abs(lat70)/360.))
for kk = 0L, n_elements(altitude) - 1L do begin
		temp = reform(intwaccmtemps[lat70,*,kk])
		latweights = temp*0.
		for i = 0L, n_elements(lat70) - 1L do latweights[i,*] = latweight[i]
		x = where(temp gt 100. and temp lt 600.,nx)
		if nx gt 0L then tempSeries[tday,kk] = TOTAL(latweights[x]*temp[x],/nan)/total(latweights[x],/nan)
endfor
strdates[tday] = '20' + ifile
meandailyProfile[iii,*] = tempseries[tday,*]

tday = tday + 1L
print, iii
endfor

save, filename = '/Volumes/MacD68-1/france/WACCM_paper/Post_process/Dailymean_lat_weight_profile_70-90N.sav',dates, dailymeanprofile, altitude

end