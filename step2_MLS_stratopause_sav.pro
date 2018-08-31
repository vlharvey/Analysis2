;-----------------------------------------------------------------------------------------------------------------------------
; Reads in MLS data and returns stratopause Z, T at all latitude and longitudes.
;	 -------------------------------
;       |         Jeff France           |
;       |         LASP, ATOC            |
;       |    University of Colorado     |
;	 -------------------------------
;
;
;
;-----------------------------------------------------------------------------------------------------------------------------
;
@stddat			; Determines the number of days since Jan 1, 1956
@kgmt			; This function computes the Julian day number (GMT) from the
								;    day, month, and year information.
@ckday			; This routine changes the Julian day from 365(6 if leap yr)
								;    to 1 and increases the year, if necessary.
@kdate			; gives back kmn,kdy information from the Julian day #.


;-----------------------------------------------------
;The following determines the date range

lstdy = 1L
lstmn = 3L
lstyr = 2013L

;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
days = 1L
;read,' Enter number of days  ', days
; set endday equal to startday so for each plot, only one day is plotted, for each elements of 'days'
ledmn = lstmn
leddy = lstdy
ledyr = lstyr
n_dayAvg = 1L    ; 1 day running average



mon=['jan_','feb_','mar_','apr_','may_','jun_',$
     'jul_','aug_','sep_','oct_','nov_','dec_']
if lstyr lt 91 then lstyr=lstyr+2000
if ledyr lt 91 then ledyr=ledyr+2000
if lstyr lt 1900 then lstyr=lstyr+1900
if ledyr lt 1900 then ledyr=ledyr+1900
if lstyr lt 1991 then stop,'Year out of range '
;if ledyr lt 1991 then stop,'Year out of range '
z = stddat(lstmn,lstdy,lstyr,lstday)
z = stddat(ledmn,leddy,ledyr,ledday)
if ledday lt lstday then stop,' Wrong dates! '
iyr = lstyr
idy = lstdy
imn = lstmn
mn1 = ledmn
dy1 = leddy
yr1 = ledyr

z = kgmt(imn,idy,iyr,iday)
startday = julday(lstmn, lstdy, lstyr)
endday = julday(ledmn, leddy, ledyr)

idir='/Volumes/earth/aura6/data/MLS_data/Datfiles_Grid/'
odir='/Volumes/earth/aura6/data/MLS_data/Datfiles_Strat/'

; -----------------------------
; loop over each day
nth_day = 10
for day = 0L, days - 1L do begin

	nday = iday
;-------------------------------------------------------------------------

	; define variables
	startday      = julday(lstmn,lstdy,2000 + lstyr)
	endday        = julday(ledmn,leddy,2000 + ledyr)

	comment = 'Stratopause algorithm: 1)apply an 11 km boxcar smoothing 2)find tmax_smooth between 20-85 km 3)find tmax within 15 km of tmax_smooth '$
	+'4)lapse rate must be positive (negative) for 5 km below (above) tmax to define as stratopause.'



  	;----------------------------------------------
  	; determine the date
  	ckday,iday,iyr
  	kdate,float(iday),iyr,imn,idy
        z = stddat(imn,idy,iyr,ndays)
  	if iyr lt 2000 then iyr1=iyr-1900
  	if iyr ge 2000 then iyr1=iyr-2000
  	syr   = string(FORMAT='(I4.4)',iyr)
  	smn   = string(FORMAT='(I2.2)',imn)
  	sdy   = string(FORMAT='(I2.2)',idy)
  	doy = string(Format='(I03)',iday)
  	ifile = string(syr+smn+sdy)
  	;----------------------------------------------

    	file =  idir + 'MLS_T_grid5_v3.3_' + ifile + '.sav'
    	spawn, 'ls ' + file, nfile
    	if nfile ne '' then begin
      		version = 'v3.3'
      		restore, file

       		; if file exists, then determine daily mean T and P
      		MLStempAvg = t_grid
      		MLSpressAvg = p_grid

		;Define variables
		mlslat = lat
		mlslon = lon
		n_lat = n_elements(lat)
		n_lon = n_elements(lon)
		n_levels    = n_elements(altitude)
		dlon        = 360./(n_lon - 1L)
		dlat        = 180./(n_lat - 1L)
		tempAvg     = fltarr(n_dayAvg, n_lon, n_lat, n_levels) - 99.
		pressAvg    = fltarr(n_dayAvg, n_lon, n_lat, n_levels) - 99.
		MLSgphtave = fltarr(n_lon, n_lat, n_levels) * !values.f_nan
		MLStheta       = fltarr(n_lon, n_lat, n_levels) *!values.f_nan
		MLSGPHT     = fltarr(n_lon, n_lat, n_levels)*!values.f_nan
		MLSGPheightstrat = fltarr(10L, n_lon, n_lat)*!values.f_nan
		MLSthetastrat = fltarr(10L, n_lon, n_lat)*!values.f_nan
		MLSpressstrat = fltarr(10L, n_lon, n_lat)*!values.f_nan
		MLSTstrat   = fltarr(10L, n_lon, n_lat)*!values.f_nan
		MLSzstrat   = fltarr(10L, n_lon, n_lat)*!values.f_nan
		MLSzstratconvert= fltarr(10L, n_lon, n_lat)*!values.f_nan
		Alt         = fltarr(10L, n_lon, n_lat)*!values.f_nan





    		for jj = 0L, n_lat - 1L do begin
   			for kk = 0L, n_lon - 1L do begin
     				x = where(MLStempAvg[kk,jj,*] gt -99. and MLSpressAvg[kk,jj,*] gt -99., nx)
      				if nx gt 0L then MLStheta[kk,jj,*]  = MLStempAvg[kk,jj,*] * ((1000. / MLSpressAvg[kk,jj,*])^(.286))

    			endfor
    		endfor

		;DEFINE STRATOPAUSE 
		;DEFINITION: Stratopause algorithm:
		;1) apply 11 km boxcar smoothing
		;2) find tmax_smooth between 20-85 km 
		;3) find tmax from raw profile within 15 km of tmax_smooth
		;4) lapse rate must be positive (negative) for 5 km below (above) tmax to define as stratopause.
		; Multiple stratopauses can be flagged
		for ii = 0L, n_lon - 1L do begin

			for jj = 0L, n_lat - 1L do begin
  				nstrat = 0L

				x = where (MLStempAvg[ii,jj,*] le 0., nx)
				if nx gt 0L then MLStempAvg[ii,jj,x] = !values.f_nan
	
    				MLSsmoothT = reform(MLStempAvg[ii,jj,*])
  				if max(MLStempAvg[ii, jj, *], /nan) gt 0. then begin

					nsmooth = 0L
					tsmooth1 = 0.
					;1) apply 11 km boxcar smoothing
					tsmooth = smooth(mlssmootht,11,/nan)

					;2) find tmax_smooth between 20-85 km 
      					for kk = 10L, n_levels - 10L do begin
        					;Stratopause definition using the static stability
        					; Determines if ii,jj,kk location is a local max temperature

          					if (altitude[kk] gt 20. and altitude[kk] lt 85. and tsmooth[kk] gt -99. and $
              					(tsmooth[kk] gt tsmooth[kk+1L]) and (tsmooth[kk] gt tsmooth[kk-1L]) and $
	              				(tsmooth[kk] gt tsmooth[kk+3L]) and (tsmooth[kk] gt tsmooth[kk-3L]) and $
        	      				(tsmooth[kk] gt tsmooth[kk+4L]) and (tsmooth[kk] gt tsmooth[kk-4L]) and $
              					(tsmooth[kk] gt tsmooth[kk+5L]) and (tsmooth[kk] gt tsmooth[kk-5L]) and $
             					(tsmooth[  kk] gt tsmooth[ kk+2L]) and (tsmooth[ kk] gt tsmooth[kk-2L])) then begin
            					nsmooth = nsmooth + 1L
							if tsmooth[kk] gt tsmooth1 then begin
					 			stratsmooth = altitude[kk]
								tsmooth1 = tsmooth[kk]
	 						endif
      	  					endif
        				endfor


					;3) find tmax from raw profile within 15 km of tmax_smooth
					;4) lapse rate must be positive (negative) for 5 km below (above) tmax to define as stratopause.
					for kk = 10L, n_levels - 10L do begin	
						IF NSTRAT lt 10L then begin
	
	      						if (altitude[kk] gt (stratsmooth - 15.) and altitude[kk] lt (stratsmooth + 15.) and MLSsmoothT[kk] gt -99. and $
        	     					(MLSsmoothT[kk] gt MLSsmoothT[kk+1L]) and (MLSsmoothT[kk] gt MLSsmoothT[kk-1L]) and $
              						(MLSsmoothT[kk+1] gt MLSsmoothT[kk+2L]) and (MLSsmoothT[kk+2] gt MLSsmoothT[kk+3L]) and $
              						(MLSsmoothT[kk+3] gt MLSsmoothT[kk+4L]) and (MLSsmoothT[kk+4] gt MLSsmoothT[kk+5L]) and $
              						(MLSsmoothT[kk-1] gt MLSsmoothT[kk-2L]) and (MLSsmoothT[kk-2] gt MLSsmoothT[kk-3L]) and $
              						(MLSsmoothT[kk-3] gt MLSsmoothT[kk-4L]) and (MLSsmoothT[kk-4] gt MLSsmoothT[kk-5L]) and $
              						(MLSsmoothT[kk] gt MLSsmoothT[kk+3L]) and (MLSsmoothT[kk] gt MLSsmoothT[kk-3L]) and $
              						(MLSsmoothT[kk] gt MLSsmoothT[kk+4L]) and (MLSsmoothT[kk] gt MLSsmoothT[kk-4L]) and $
              						(MLSsmoothT[kk] gt MLSsmoothT[kk+5L]) and (MLSsmoothT[kk] gt MLSsmoothT[kk-5L]) and $
              						(MLSsmootht[  kk] gt MLSsmootht[ kk+2L]) and (MLSsmootht[ kk] gt MLSsmoothT[kk-2L])) and nstrat lt 10L then begin



               							MLSthetastrat[nstrat,ii, jj] = MLStheta[ii,jj,kk]
                						MLSTstrat[nstrat,ii,jj] = MLStempAvg[ii,jj,kk]
                						MLSpressstrat[nstrat,ii,jj] = MLSpressAvg[ii,jj,kk]
                						MLSGPheightstrat[nstrat,ii,jj] = MLSgpht[kk]
                						MLSzstrat[nstrat,ii,jj] = altitude[kk]
                                				nstrat = nstrat + 1L
	   						endif
						endif ;NSTRAT lt 10L

					endfor ;n_levels


				endif ;max(MLStempAvg[ii, jj, *], /nan) gt 0.
			ENDFOR ;n_lat
		ENDFOR ;n_lon
	
		if max(MLSzstrat,/nan) gt 0. then begin
			print, 'file saved'
			x = where (lat gt 82. or lat lt -82.)
			MLSTstrat[*,*,x] = !values.f_nan
			MLSthetastrat[*,*,x] = !values.f_nan
			MLSGPheightstrat[*,*,x] = !values.f_nan
			MLSpressstrat[*,*,x] = !values.f_nan
			MLSzstrat[*,*,x] = !values.f_nan
			date = ifile



    			SAVE, MLSTstrat,MLSZstrat,MLSthetastrat,MLSGPheightstrat,lat,lon,date,MLSpressstrat, $
   			FILENAME = odir + 'MLS_stratopause_' + ifile +'.sav'
		endif
	endif
	; ---------------------------------------------------------------------
	; ---------------------------------------------------------------------

	; CHANGE DATE TO NEXT DAY
	iday = nday
	startday = startday + 1L
	endday = endday + 1L
	caldat, startday, lstmn, lstdy, lstyr
	lstyr = lstyr - 2000L
	ledmn = lstmn
	leddy = lstdy
	ledyr = lstyr
	print, lstyr, ledmn, leddy
	nday  = nday + 1L
	mon = ['jan_','feb_','mar_','apr_','may_','jun_',$
	       'jul_','aug_','sep_','oct_','nov_','dec_']
	if lstyr lt 91 then lstyr=lstyr+2000
	if ledyr lt 91 then ledyr=ledyr+2000
	if lstyr lt 1900 then lstyr=lstyr+1900
	if ledyr lt 1900 then ledyr=ledyr+1900
	if lstyr lt 1991 then stop,'Year out of range '
	;if ledyr lt 1991 then stop,'Year out of range '
        z = stddat(lstmn,lstdy,lstyr,lstday)
        z = stddat(ledmn,leddy,ledyr,ledday)
	if ledday lt lstday then stop,' Wrong dates! '
	iyr = lstyr
	idy = lstdy
	imn = lstmn
        z = kgmt(imn,idy,iyr,iday)

endfor; day


end
