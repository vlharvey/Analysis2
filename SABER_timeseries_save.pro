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






;-----------Determine the range of days and set up the counter for the date----------------------

lstdy = 1L
lstmn = 1L
lstyr = 12L
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
days = 200L
;read,' Enter number of days  ', days
;  set endday equal to startday so for each plot, only one day is plotted, for each elements of 'days'
ledmn = lstmn
leddy = lstdy
ledyr = lstyr

doy = fltarr(days)
idoy = 25L
dates = strarr(days)
meandailyProfile = fltarr(days,201L)

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


dir='/Volumes/MacD68-1/france/SABER_data/Datfiles/saber_T_grid5_20'

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
dates[iii] = '20'+ifile

nc = 0L
ifiles=file_search(dir+ifile+'.sav',count=nfile)
if imn ne 5L or idy ne 1L then begin

	if ifiles[0] eq '' then continue
   	 restore, ifiles[0]
endif

	sabertemps = t_grid


if iyr eq 2012L and imn eq 1L and idy eq 1L then begin
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

	save, filename = '/Volumes/MacD68-1/france/WACCM_paper/Post_process/SABER_timeseries_SH_20'+syr+'.sav', tempseries, dateseries, altitude, strdates, smoothtemp
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

	save, filename = '/Volumes/MacD68-1/france/ES_paper/Post_process/SABER_timeseries_NH_20'+syr+'.sav', tempseries, dateseries, altitude, strdates, smoothtemp
	strdates = strarr(184L)
	dateSeries = findgen(184L) + 1L
	tempSeries = fltarr(n_elements(dateSeries), n_elements(altitude))
	tday = 0L
	mode = 1L
endif

x = where(SABERtemps le 100. and SABERtemps gt 1500.,nx)
if nx gt 0L then SABERtemps[x] = !values.f_nan
if mode eq 0L then lat70 = where(lat ge 70L)
if mode eq 1L then lat70 = where(lat le -70L)
latweight = SQRT(COS(2*3.141592*abs(lat70)/360.))
for kk = 0L, n_elements(altitude) - 1L do begin
		temp = reform(SABERtemps[*,lat70,kk])
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

save, filename = '/Volumes/MacD68-1/france/ES_paper/Post_process/SABER_Dailymean_lat_weight_profile_70-90N.sav',dates, meandailyprofile, altitude

end