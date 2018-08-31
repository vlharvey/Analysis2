;-----------------------------------------------------------------------------------------------------------------------------
;SAVES CODE FOR ELEVATED STRATOPAUSE DAYS FOR WACCM.
;   	 -------------------------------
;       |         Jeff France           |
;       |         LASP, ATOC            |
;       |    University of Colorado     |
;       |     modified: 04/11/2012      |
;   	 -------------------------------
;
;
;
;-----------------------------------------------------------------------------------------------------------------------------
;
@/Volumes/MacD68-1/france/idl_files/bsort

;------------RESTORE SAVE FILE WITH DAYS, LON, LAT FOR WACCM AND GEOS----------


restore, '/Volumes/MacD68-1/france/WACCM_paper/Post_process/area_weighted_strat_daily_mean.sav'
restore,'/Volumes/MacD68-1/france/WACCM_paper/Post_process/Dailymean_lat_weight_profile_70-90N.sav'



minjump = 25.


;elevatedstrat,dates,meanstrat
dum = indgen(n_elements(dates))

xx = where(meanstrat le 0.)
meanstrat[xx] = !values.f_nan
y = smooth(meanstrat,9,/nan)

dayzeroindex = fltarr(n_elements(dates))
ESindex = fltarr(n_elements(dates))

totdays = n_elements(dates)
dayzeros = fltarr(200)
stratjump = fltarr(200)
tempstddev = fltarr(200)
nd = 0L
for ii = 7L, n_elements(y) - 8L do begin

	jump = mean(meanstrat[ii+3L:ii+7L]) - mean(meanstrat[ii-7L:ii-3L])
	if jump gt 10. then begin ; require 
		x = where(meanstrat[ii-7L:ii+7L] eq max( meanstrat[ii-7L:ii+7L]))
		dayzeros[nd] = x + ii - 7L
		stratjump[nd] = max(meanstrat[ii+3L:ii+7L]) - min(meanstrat[ii-7L:ii-3L])
		
		tempstddev[nd] = stddev(meandailyprofile[x + ii - 7L, where((altitude - min(meanstrat[ii-7L:ii-3L])) eq min((altitude - min(meanstrat[ii-7L:ii-3L]))))$
		:where((altitude - max(altitude - meanstrat[ii+3L:ii+7L])) eq min((altitude - max(meanstrat[ii+3L:ii+7L]))))],/nan)
		
		if stratjump[nd] ge minjump then dayzeroindex[ii] = 1L
		
		nd = nd + 1L
		ii = ii+20L
	endif
endfor


x = where(stratjump ge minjump,nx)

dayzeros = reform(dayzeros[x])

x = where(dayzeroindex eq 1L,nx)

for ii = 0L, nx - 1L do ESindex[x[ii]:x[ii] + 5L] = 1L


; Determine ES period following each day zero-- ES lasts until stratopause descends to within one standard deviation above the daily mean

meanstrat = smooth(meanstrat,5,/nan)
x = where(dayzeroindex eq 1L,nx)
for ii = 0L, nx - 1L do begin
	xdate = where(strmatch(strmid(dates,4,4), strmid(dates[x[ii]],4,4)) eq 1L) + 6L ; Add 6 because the first 5 days are assumed to be elevated

	for jj = 0L, 90L do begin ; loop over days following ES day zero to determine how long it's elevated
		xdate2 = xdate + jj
		esday = x[ii] + 6L + jj 			; Add 6 because the first 5 days are assumed to be elevated
		y = where(meanstrat[xdate2] gt 20. and meanstrat[xdate2] lt 100.,ny)
		dailymean = mean(meanstrat[xdate2[y]],/nan)
		dailystddev = stddev(meanstrat[xdate2[y]],/nan)
		if meanstrat[esday] gt (dailymean + dailystddev) then begin
			ESindex[esday] = 1L
		endif else break
	endfor
endfor
;for ii = 0L, n_elements(esindex) - 1L do begin;
;	if esindex[ii] eq 0L then continue
;	xdate = where(strmatch(strmid(dates,4,4), strmid(dates[ii],4,4)) eq 1L) 
;		y = where(meanstrat[xdate] gt 20. and meanstrat[xdate] lt 100.,ny)
;		dailymean = mean(meanstrat[xdate[y]],/nan)
;		dailystddev = stddev(meanstrat[xdate[y]],/nan)
;
;		if meanstrat[ii] lt (dailymean + dailystddev) then begin
;		print, meanstrat[ii],dailymean, dailystddev
;			ESindex[ii] = 0L
;		endif
;endfor

x = where(esindex eq 1L, comp = y)
esindex = x
noesindex = y

save,tempstddev,dayzeros, dayzeroindex, ESindex, noESindex, stratjump,dates, filename='/Volumes/MacD68-1/france/ES_paper/Post_process/elevated_strat.sav'
END	