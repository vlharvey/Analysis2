;-----------------------------------------------------------------------------------------------------------------------------
; Reads in WACCM data and returns stratopause height and temperature
;	 -------------------------------
;       |         Jeff France           |
;       |         LASP, ATOC            |
;       |    University of Colorado     |
;       |     modified: 05/24/2010      |
;	 -------------------------------
;
;
;
;-----------------------------------------------------------------------------------------------------------------------------
;
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


;-----------------------------------------------------
;The following determines the date range

restore, '/Volumes/MacD68-1/france/WACCM_paper/Post_process/elevated_strat.sav'

ndates = 30.*n_elements(dayzeros)
dir='/Volumes/MacD68-2/france/data/WACCM_4/Datfiles/mee00fpl_FW2.cam2.h3.dyns.'

ESday = fltarr(n_elements(dayzeros)*30L)
Date = strarr(n_elements(dayzeros)*30L)
nday = 0L
dayofES = 1L

for iES = 0L, n_elements(dayzeros) - 1L do begin
	ydate = where(dates eq dayzerodates[iES])
	for iday = 0L, 29L do begin
		ifile = dates[ydate+iday]
		print, ifile

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
;       		print,'read ',name,' dimension ',dim
    		endfor
    		for ivar=0,result0.nvars-1 do begin
        		result=ncdf_varinq(ncid,ivar)
        		ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
        		if result.name eq 'latitude' then alat=data
        		if result.name eq 'longitude' then alon=data
        		if result.name eq 'theta' then th=data
        		if result.name eq 'IPV' then pv2=data
        		if result.name eq 'P' then p2=data
        		if result.name eq 'U' then u2=data
        		if result.name eq 'V' then v2=data
        		if result.name eq 'QDF' then qdf2=data
        		if result.name eq 'Q' then q2=data
        		if result.name eq 'GPH' then gph2=data
        		if result.name eq 'TTGW' then ttgw2=data
        		if result.name eq 'SF' then sf2=data
        		if result.name eq 'MARK' then marksf2=data
        		print,ivar,result.name,min(data),max(data)
		endfor
    		ncdf_close,ncid

	; CONVERT GEOPOTENTIAL TO GEOMETRIC HEIGHT

	k = 8.314
	g=9.81
	M = .02896
	rtd=double(180./!pi)
	dtr=1./rtd
	ks=1.931853d-3
	ecc=0.081819
	gamma45=9.80
	gph2 = gph2/1000.

zgrd=fltarr(nr,nc,nth)   ; 3-d geopotential height
altitude_grid = findgen(100L)+1L
altitude =altitude_grid
n_lat = n_elements(alat)
n_lon = n_elements(alon)
intWACCMgph = fltarr(n_lat,n_lon,46)
waccmtemp = gph2*0.
intwaccmmark = intwaccmsf*0.


pressure_grid = exp(findgen(46)*.32855-9.11)

for ii = 0L, n_elements(alat) - 1L do begin
	sin2=sin( (alat(ii)*dtr)^2.0 )
	numerator=1.0+ks*sin2
	denominator=sqrt( 1.0 - (ecc^2.0)*sin2 )
	gammas=gamma45*(numerator/denominator)
	r=6378.137D/(1.006803-(0.006706*sin2))
	for jj = 0L, n_elements(alon) - 1L do begin
		x =interpol(reform(gph2[ii,jj,*]),reform(p2[ii,jj,*]),pressure_grid)
		intWACCMgph[ii,jj,*] = x    
	endfor
endfor



 		if ies eq 0L and iday eq 0L then begin
			phase60	= fltarr(n_elements(dayzeros)*30L, n_elements(th))
			phase70	= fltarr(n_elements(dayzeros)*30L, n_elements(th))
			phase75	= fltarr(n_elements(dayzeros)*30L, n_elements(th))
			amp60	= fltarr(n_elements(dayzeros)*30L, n_elements(th))
			amp70	= fltarr(n_elements(dayzeros)*30L, n_elements(th))
			amp75	= fltarr(n_elements(dayzeros)*30L, n_elements(th))
		endif

			x = where(abs(alat - 75.) eq min(abs(alat - 75.)))
			Date[nday] = ifile

	   		for kk = 0L, n_elements(pressure_grid) - 1L do begin
				z = where(abs(alat - 60.) eq min(abs(alat - 60.)),nz)
	    	 	if nz ge 1L then begin
					level = reform(intWACCMgph[z,*,kk])
					X = alon 
					Y = level
					weights = Y*0. + 1. ; Define a vector of weights.
					levelmean = mean(Y,/nan)
					A2guess = where(abs(Y - levelmean) eq min(abs(Y - levelmean),/nan))
					A2guess = X[A2guess[0]]
					A = [1000.0,A2guess, levelmean] ; Provide an initial guess of the function's parameters.
					yfit = CURVEFIT(X, Y, weights, A, SIGMA, FUNCTION_NAME='gfunct') ; fit data to sine function : y = a[0]*sin(2*pi*x+a[1]) + a[2]
;					PRINT, 'Function parameters: ', A ; Print the parameters returned in A. 
					Amp60[nday,kk] = abs(A[0])
					Phase60[nday,kk] = A[1]
				endif

				z = where(abs(alat - 70.) eq min(abs(alat - 70.)),nz)
	    	 	if nz ge 1L then begin
					level = reform(intWACCMgph[z,*,kk])
					X = alon
					Y = level
					weights = Y*0. + 1. ; Define a vector of weights.
					levelmean = mean(Y,/nan)
					A2guess = where(abs(Y - levelmean) eq min(abs(Y - levelmean),/nan))
					A2guess = X[A2guess[0]]
					A = [1000.0,A2guess, levelmean] ; Provide an initial guess of the function's parameters.
					yfit = CURVEFIT(X, Y, weights, A, SIGMA, FUNCTION_NAME='gfunct') ; fit data to sine function : y = a[0]*sin(2*pi*x+a[1]) + a[2]
;					PRINT, 'Function parameters: ', A ; Print the parameters returned in A. 
					Amp70[nday,kk] = abs(A[0])
					Phase70[nday,kk] = A[1]
				endif

				z = where(abs(alat - 75.) eq min(abs(alat - 75.)),nz)
	    	 	if nz ge 1L then begin
					
					level = reform(intWACCMgph[z,*,kk])
					X = alon 
					Y = level
					weights = Y*0. + 1. ; Define a vector of weights.
					levelmean = mean(Y,/nan)
					A2guess = where(abs(Y - levelmean) eq min(abs(Y - levelmean),/nan))
					A2guess = X[A2guess[0]]
					A = [1000.0,A2guess, levelmean] ; Provide an initial guess of the function's parameters.
					yfit = CURVEFIT(X, Y, weights, A, SIGMA, FUNCTION_NAME='gfunct') ; fit data to sine function : y = a[0]*sin(2*pi*x+a[1]) + a[2]
					PRINT, 'Function parameters: ', A ; Print the parameters returned in A. 
					Amp75[nday,kk] = abs(A[0])
					Phase75[nday,kk] = A[1]
			endif
		endfor

nday = nday+1L
dayofes = dayofes+1L
if dayofES eq 31L then dayofES = 1L
endfor
endfor
Save, filename = '/Volumes/MacD68-1/france/ES_paper/Post_process/WACCM_ES_wave-1_phase_pre_process_5km.sav',amp60,amp70,amp75,phase60,phase70,phase75,Date, pressure_grid


end