
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
@/Volumes/MacD68-1/france/idl_files/rd_waccm_nc3
@/Volumes/MacD68-1/france/idl_files/frac_index


SABERheights = 0
GEOS5heights = 0
WACCMheights = 0
SABERheights = 0
marker = 0.

;-----------Determine the range of days and set up the counter for the date----------------------

lstmn = 1L
lstdy = 25L
lstyr = 2L
;read,' Enter starting date (month, day, year) ',lstmn,lstdy,lstyr
days = 3600L
;read,' Enter number of days  ', days
;  set endday equal to startday so for each plot, only one day is plotted, for each elements of 'days'
ledmn = lstmn
leddy = lstdy
ledyr = lstyr

doy = fltarr(days)
idoy = 25L

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

;  ----------- define variables and set plot window------------------
; loop over each day
nth_day = 10
nday = iday
ifile1 = 0L
ifile2 = 0L
ifile3 = 0L
height1 = 0L
height2 = 0L
height3 = 0L
temp1 = 0L
temp2 = 0L
temp3 = 0L
alon = findgen(96)
alat = findgen(72)
nth=0L
saberyaw = 0L
ngeos=0L
FOR iii = 0L, days-1 DO BEGIN          	; n_days is set in plot code and it determines how many days will be averaged together in each plot.
  	nsave = 0L

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
  	ifile = string(syr+smn+sdy)

  	ifile3 = ifile2
  	ifile2 = ifile1
  	ifile1 = ifile

	if iii eq 0L then dates = strarr(days)
	dates[iii] = '20' + ifile


if imn eq 1l and idy eq 1l then idoy = 0l
idoy = idoy + 1L
doy[iii] = idoy

print, doy[iii]






	;--------RESTORE DATA-------------------------


marker = marker*0.
		;------------restore GEOS5 data from IDL save files----------------------------------
		geos5file = 0L
		mark3 = 0
  		print ,ifile
 		; restore GEOS5 sav files
  		@/Volumes/MacD68-1/france/GEOS5_data/Pre_process/read_GEOS5  ; read in GEOS5 data
	if max(mark3,/nan) gt 0L then begin
  		GEOS5lat = alat
       		GEOS5lon = alon
       		GEOS5lat = alat
       		GEOS5mark4 = mark4
		latgrid = alat
		longrid = alon
		ngeos = ngeos+1L

	x = where(mark4 le -10.,nx)
	if nx gt 0L then mark4[x] = 0.
	x = where(mark4 gt 0. or mark4 lt 0.,nx)
	if nx gt 0L then mark4[x] = mark4[x]/abs(mark4[x])
endif
	
 	;------------restore SABER data from IDL save files----------------------------------
  	print ,ifile
  	directory='/Volumes/MacD68-2/france/data/SABER_data/'

    	;------------restore SABER data from IDL save files----------------------------------
  	print ,ifile
  	directory='/Volumes/MacD68-1/france/SABER_data/Datfiles/stratopause/'
SABERfile = 0
 	; restore SABER sav files
  	file =  directory + 'SABER_pause_height_grid_smooth_for_20' + ifile + '.sav'
  	spawn, 'ls ' + file, nfile
  	if nfile ne '' then begin
    		restore, file
    		SABERlat = lat
    		SABERlon = lon
    		SABERT = Tempstrat
		SABERfile = 1L
	endif

	
	
;-------------------------------------DEFINE VARIABLES--------------------------------------------------------------


		if iii eq 0L then day1 = date
		if iii eq days - 1L then day2 = date
	    	
		if ngeos eq 1L then begin
    			GEOS5heights = fltarr(days, n_elements(GEOS5lon), n_elements(GEOS5lat))*0. - 999.
    			GEOS5temps = fltarr(days,n_elements(GEOS5lon), n_elements(GEOS5lat))*0. - 999.
    			GEOS5markers = fltarr(days, n_elements(GEOS5lon), n_elements(GEOS5lat))*0. - 999.
  				GEOS5double = fltarr(days, n_elements(geos5lon), n_elements(geos5lat))*0.
		endif

	if iii eq 0L then day1 = date
	if iii eq days - 1L then day2 = date
    	





		if iii eq days - 1L then endday = date
	    	

		if n_elements(SABERheights) eq 1L then begin
    			SABERheights = fltarr(days, n_elements(SABERlon), n_elements(SABERlat))*0. - 999.
    			SABERtemps = fltarr(days,n_elements(SABERlon), n_elements(SABERlat))*0. - 999.
    			SABERmarks = fltarr(days,n_elements(SABERlon), n_elements(SABERlat))*0.
    			SABERdouble = fltarr(days, n_elements(SABERlon), n_elements(SABERlat))*0.
		endif
		
icount = 0

;x = where(geos5lat gt .,nx)
;if nx gt 0L then mark4[x,*,*] = 0.

;-------------------------------LOOP OVER LAT/LON TO DETERMINE STRATOPAUSE Z/T---------------------------------------------------
	for ii = 0L, n_elements(geos5lon)-1L do begin
		for jj = 0L, n_elements(geos5lat)-1L do begin
 
      			if SABERfile eq 1L then begin 

      			

				y = where(SABERzstrat[*,ii,jj] gt 20. and SABERzstrat[*,ii,jj] lt 85.,ny)
				if ny gt 0L then begin
					x = where(SABERt[y,ii,jj] eq max(SABERt[y,ii,jj],/nan),nx)
					if nx eq 1 then begin
 						SABERheights[iii,ii,jj] = SABERzstrat[y[x],ii,jj]
						SABERtemps[iii,ii,jj] = SABERt[y[x],ii,jj]
						x = where(abs(SABERheights[iii,ii,jj] - zgrd[ii,jj,*]) eq min(abs(SABERheights[iii,ii,jj] - zgrd[ii,jj,*]),/nan),nx)

						if nx eq 2L then begin
							if SABERheights[iii,ii,jj] lt 80. and mark4[jj,ii,x[0]] eq mark4[jj,ii,x[1]]$
							 then SABERmarks[iii,ii,jj] = mark4[jj,ii,x[0]]
						endif							
						if x[0] eq 0L and nx eq 1L then SABERmarks[iii,ii,jj] = mark4[jj,ii,x]
						
						if x[0] gt 0L and x[0] lt n_elements(zgrd[0,0,*]) and nx eq 1L then begin
						if (SABERheights[iii,ii,jj] - zgrd[ii,jj,x-1]) lt (SABERheights[iii,ii,jj] - zgrd[ii,jj,x+1]) then begin
							if SABERheights[iii,ii,jj] lt 80. and mark4[jj,ii,x] eq mark4[jj,ii,x-1]$
						 	then SABERmarks[iii,ii,jj] = mark4[jj,ii,x]
						endif
						if ngeos ge 1L then begin
 						if (SABERheights[iii,ii,jj] - zgrd[ii,jj,x-1]) gt (SABERheights[iii,ii,jj] - zgrd[ii,jj,x+1]) then begin
							if SABERheights[iii,ii,jj] lt 80. and mark4[jj,ii,x] eq mark4[jj,ii,x+1]$
							 then SABERmarks[iii,ii,jj] = mark4[jj,ii,x]
						endif
						endif
						endif
					endif
				endif

				if ii gt 0L and jj gt 0L then begin
					diff = abs(SABERheights[iii,ii,jj] - SABERheights[iii,ii,jj-1L])
					diff2 = abs(SABERheights[iii,ii,jj] - SABERheights[iii,ii-1L,jj])
					diff3 = abs(SABERheights[iii,ii,jj] - SABERheights[iii,ii-1L,jj-1L])
					if diff ge 10. and diff le 100. or diff2 ge 10. and diff2 le 100. $
					and diff3 ge 10. and diff3 le 100. then SABERdouble[iii,ii,jj] = 1L
				endif
			ENDIF; SABERFILE=1

		    	
		    	
			if ngeos eq 1L then begin
    				GEOS5heights = fltarr(days, n_elements(GEOS5lon), n_elements(GEOS5lat))*0. - 999.
    				GEOS5temps = fltarr(days,n_elements(GEOS5lon), n_elements(GEOS5lat))*0. - 999.
    				GEOS5markers = fltarr(days,n_elements(GEOS5lon), n_elements(GEOS5lat))*0. - 999.
			endif
		

		    	
		    	

		    	
					; STATIC STABILITY AS IN MANNEY
		    	

			endfor
		endfor

  	startday = startday + 1L
  	endday = endday + 1L
  	caldat, startday, lstmn, lstdy, lstyr
  	lstyr = lstyr - 2000L
  	ledmn = lstmn
  	leddy = lstdy
  	ledyr = lstyr
 	print, lstyr, ledmn, leddy
  	nday  = nday + 1L
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
endfor


save, filename = '/Volumes/MacD68-1/france/ES_paper/Post_process/SABER_temp_height_mark.sav',$
	SABERtemps,SABERheights,SABERmarks, dates, lat, lon, doy


end
