
restore,'/Volumes/MacD68-1/france/ES_paper/Post_process/SABER_Dailymean_lat_weight_profile_70-90N.sav'
;,dates, meandailyprofile, altitude
x = where(meandailyprofile eq 0.)
meandailyprofile[x] = !values.f_nan


SABERstratopausepoly = fltarr(n_elements(dates))
SABERstratopause = fltarr(n_elements(dates))
for ii = 0l, n_elements(dates) - 1L do begin
;Result = POLY_FIT( X, Y, Degree [, CHISQ=variable] [, COVAR=variable] [, /DOUBLE] [, MEASURE_ERRORS=vector]
; [, SIGMA=variable] [, STATUS=variable] [, YBAND=variable] [, YERROR=variable] [, YFIT=variable] )

	st = reform(meandailyprofile[ii,15:150])
	z1 = reform(altitude[15:150])
	x = where(st lt 0. or finite(st) eq 0.,nx, comp = y)
	if nx gt 0L then st[x] = !values.f_nan
	if max(y) ge 0 then begin
	x = poly_fit(z1[y],st[y],6)
	z = findgen(400)*.35 + 10
	SABERt = x[0,6]*z^6 + x[0,5]*z^5 + x[0,4]*z^4 + x[0,3]*z^3+x[0,2]*z^2+x[0,1]*z+x[0,0]
	
	yz = where(finite(st))
	yzalt = max(z1[yz])
	top = where(z gt yzalt, nt)
	if nt gt 0L then SABERt[top] = !values.f_nan
	x = where(z gt 20. and z le 85.)	
	xx = where(SABERt[x] eq max(SABERt[x], /nan), nxx)
	if nxx ge 1L then begin
	if SABERt[x[xx[0]]] gt 0. and SABERt[x[xx[0]]]lt 1000. then SABERstratopausepoly[ii] = z[x[xx[0]]]
	endif
	endif

	SABERt = reform(meandailyprofile[ii,*])
	x = where(altitude gt 20. and altitude le 85.)	
	xx = where(SABERt[x] eq max(SABERt[x], /nan), nxx)
	if nxx ge 1L then begin
		if SABERt[x[xx[0]]] gt 0. and SABERt[x[xx[0]]]lt 1000. then SABERstratopause[ii] = altitude[x[xx[0]]]
	endif
endfor


meanstrat = SABERstratopause
meanstratpoly = SABERstratopausepoly

save, dates,meanstratpoly,meanstrat,filename = '/Volumes/MacD68-1/france/ES_paper/Post_process/SABER_area_weighted_strat_daily_mean.sav'


END	