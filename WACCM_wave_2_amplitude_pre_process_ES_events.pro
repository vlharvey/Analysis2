;-----------------------------------------------------------------------------------------------------------------------------
;Using the database of GEOS5, WACCM, MLS, this creates a save file of stratopause height and T on geos5 grid and original strat savefile grid
;	 -------------------------------
;       |         Jeff France           |
;       |         LASP, ATOC            |
;       |    University of Colorado     |
;       |     modified: 05/11/2010      |
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



;;----------------DEFINE THE ELEVATED STRATOPAUSE--------------------------------------;;
nnx = 0L
;remove elevated stratopause days
restore, '/Volumes/MacD68-1/france/WACCM_paper/Post_process/elevated_strat.sav'
xdayzeros = dayzeros

;Define pressure grid
Pgrid = reverse(exp(findgen(100)*.163-9.21)) ;.003 to 1018. hPa
waccm_levels=n_elements(Pgrid)







;-----------Determine the range of days and set up the counter for the date----------------------


dir='/Volumes/MacD68-2/france/data/WACCM_4/Datfiles/mee00fpl_FW2.cam2.h3.dyns.'

iii = 0L
xes = 0L
iesday = -30.


days = 61L*n_elements(dayzeros)
esdates = strarr(days)
FOR iii = 0L, days-1 DO BEGIN
	if iesday eq 31. then begin
		xes = xes+1L
		iesday = -30.
	endif

	if iii eq 0L then ESbefore0after1 = fltarr(n_elements(dayzeros)*61.)
	if iesday lt 0L then ESbefore0after1[iii] = 0L	
	if iesday ge 0L then ESbefore0after1[iii] = 1L	


	sys1 = SYSTIME( /JULIAN)


	nc = 0L
	ifile = dates[xdayzeros[xes]+iesday]
	print, ifile, xes, iesday,xdayzeros[xes]
	iesday = iesday+1L

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
 ;       	if result.name eq 'U' then u2=data
 ;       	if result.name eq 'V' then v2=data
 ;       	if result.name eq 'QDF' then qdf2=data
 ;       	if result.name eq 'Q' then q2=data
        	if result.name eq 'GPH' then gph2=data
 ;       	if result.name eq 'TTGW' then ttgw2=data
 ;       	if result.name eq 'SF' then sf2=data
        	if result.name eq 'MARK' then marksf2=data
        	print,ivar,result.name,min(data),max(data)
	endfor
    		ncdf_close,ncid



waccmlat = alat
waccmlon = alon


latgrid = fltarr(nr,nc)

;-------------Put WACCM on a pressure grid----------------------

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
intWACCMtemps = fltarr(n_lat,n_lon,100L)
intWACCMmark = fltarr(n_lat,n_lon,100L)
intWACCMgph = fltarr(n_lat,n_lon,100L)
intpgrid = fltarr(n_lat,n_lon,100L)
intWACCMq = fltarr(n_lat,n_lon,100L)
intWACCMttgw = fltarr(n_lat,n_lon,100L)
intWACCMsf = fltarr(n_lat,n_lon,100L)
intWACCMipv = fltarr(n_lat,n_lon,100L)
intWACCMu = fltarr(n_lat,n_lon,100L)
intWACCMv = fltarr(n_lat,n_lon,100L)
intWACCMqdf = fltarr(n_lat,n_lon,100L)
intWACCMtheta = fltarr(n_lat,n_lon,100L)
waccmtemp = gph2*0.
intwaccmmark = intwaccmsf*0.


;put GPH on Pressure grid
waccmgph = fltarr(nr,nc,waccm_levels)
waccmtemp = fltarr(nr,nc,waccm_levels)
waccmt = fltarr(nr,nc,nth)
for ii = 0L, n_elements(alat) - 1L do begin
	sin2=sin( (alat(ii)*dtr)^2.0 )
	numerator=1.0+ks*sin2
	denominator=sqrt( 1.0 - (ecc^2.0)*sin2 )
	gammas=gamma45*(numerator/denominator)
	r=6378.137D/(1.006803-(0.006706*sin2))

	for jj = 0L, n_elements(alon) - 1L do begin
		waccmGPH[ii,jj,*] =interpol(reverse(reform(gph2[ii,jj,*])),reverse(reform(p2[ii,jj,*])),pgrid)

		WACCMt[ii,jj,*] = th[*]*( (p2[ii,jj,*]/1000.D)^(.286) )
		x =interpol(reverse(reform(waccmt[ii,jj,*])),reverse(reform(p2[ii,jj,*])),pgrid)
		WACCMtemp[ii,jj,*] = x    

	endfor
endfor
		
waccmgph = waccmgph*1000.
		
 	if iii eq 0L then begin
 		WACCMtemps = fltarr(days, WACCM_levels)-999.
 		WACCMGPHgrid50 = fltarr(nc, WACCM_levels)-999.
 		WACCMGPHgrid60 = fltarr(nc, WACCM_levels)-999.
 		WACCMGPHgrid70 = fltarr(nc, WACCM_levels)-999.
		WACCMtemps70 = fltarr(days, WACCM_levels)-999.
		WACCMtemps60 = fltarr(days, WACCM_levels)-999.
		WACCMtemps50 = fltarr(days, WACCM_levels)-999.
		WACCMheight70 = fltarr(days, WACCM_levels)-999.
		WACCMheight60 = fltarr(days, WACCM_levels)-999.
		WACCMheight50 = fltarr(days, WACCM_levels)-999.
	endif

;---WACCM----
	x = where(WACCMgph le 0.,nx)
	if nx gt 0L then WACCMgph[x] = !values.f_nan
	xlat = where(WACCMlat gt 47. and WACCMlat lt 53.)
	for ii = 0L, nc - 1L do begin
			for iz = 0L, WACCM_levels - 1L do WACCMGPHgrid50[ii,iz] = mean(WACCMgph[xlat,ii,iz],/nan)
	endfor
	x = where(WACCMgph le 0.,nx)
	if nx gt 0L then WACCMgph[x] = !values.f_nan
	xlat = where(WACCMlat gt 57. and WACCMlat lt 63.)
	for ii = 0L, nc - 1L do begin
			for iz = 0L, WACCM_levels - 1L do WACCMGPHgrid60[ii,iz] = mean(WACCMgph[xlat,ii,iz],/nan)
	endfor
	x = where(WACCMgph le 0.,nx)
	if nx gt 0L then WACCMgph[x] = !values.f_nan
	xlat = where(WACCMlat gt 67. and WACCMlat lt 73.)
	for ii = 0L, nc - 1L do begin
			for iz = 0L, WACCM_levels - 1L do WACCMGPHgrid70[ii,iz] = mean(WACCMgph[xlat,ii,iz],/nan)
	endfor

	g = 9.81
	p = pgrid


	HGPH50 = smooth(WACCMgphgrid50,5,/nan)
	dGPHdLon50 = fltarr(nc,WACCM_levels) *!values.f_nan
	for ii = 1L, nc - 2L do begin
		for iz = 0L, WACCM_levels - 1L do dGPHdLon50[ii,iz] = (HGPH50[ii+1L,iz] - HGPH50[ii-1L,iz])  /  (waccmlon[ii+1L] - waccmlon[ii-1L]) 
	endfor

	WACCMepflux50 = fltarr(WACCM_levels) *!values.f_nan
	for iz = 0L, WACCM_levels - 1L do WACCMepflux50[iz] = (1/g)*p[iz]*mean(dGPHdLon50[*,iz],/nan)


	HGPH60 = smooth(WACCMgphgrid60,5,/nan)
	dGPHdLon60 = fltarr(nc,WACCM_levels) *!values.f_nan
	for ii = 1L, nc - 2L do begin
		for iz = 0L, WACCM_levels - 1L do dGPHdLon60[ii,iz] = (HGPH60[ii+1L,iz] - HGPH60[ii-1L,iz])  /  (waccmlon[ii+1L] - waccmlon[ii-1L]) 
	endfor

	WACCMepflux60 = fltarr(WACCM_levels) *!values.f_nan
	for iz = 0L, WACCM_levels - 1L do WACCMepflux60[iz] = (1/g)*p[iz]*mean(dGPHdLon60[*,iz],/nan)


	HGPH70 = smooth(WACCMgphgrid70,5,/nan)
	dGPHdLon70 = fltarr(nc,WACCM_levels) *!values.f_nan
	for ii = 1L, nc - 2L do begin
		for iz = 0L, WACCM_levels - 1L do dGPHdLon70[ii,iz] = (HGPH70[ii+1L,iz] - HGPH70[ii-1L,iz])  /  (waccmlon[ii+1L] - waccmlon[ii-1L]) 
	endfor

	WACCMepflux70 = fltarr(WACCM_levels) *!values.f_nan
	for iz = 0L, WACCM_levels - 1L do WACCMepflux70[iz] = (1/g)*p[iz]*mean(dGPHdLon70[*,iz],/nan)


 ; if file exists, then determine daily mean T and P
		
    for kk = 0L, WACCM_levels-1L do begin

    

		level  = reform(WACCMGPHgrid50[*,kk])
		x = where(level le 0.,nx)
		if nx gt 0L then level[x] = !values.f_nan
      	; bin data into 2 degree lat bins

      	z = where(level gt 0. and level lt 90000., nz)
      	if nz gt 10. then begin
	  		X = WACCMlon[z] 
			Y = level[z]
			weights = Y*0. + 1. ; Define a vector of weights.
			levelmean = mean(Y,/nan)
			A2guess = where(abs(Y - levelmean) eq min(abs(Y - levelmean),/nan))
			A2guess = X[A2guess[0]]
			A = [1000.0,A2guess, levelmean] ; Provide an initial guess of the function's parameters.
			yfit = CURVEFIT(X, Y, weights, A, SIGMA, FUNCTION_NAME='wave2')
;			PRINT, 'Function parameters: ', A ; Print the parameters returned in A. 
			WACCMheight50[iii,kk] = abs(A[0])
		endif

		level  = reform(WACCMGPHgrid60[*,kk])
      	z = where(level gt 0. and level lt 90000., nz)
   		if nz gt 10. then begin
			X = WACCMlon[z] 
			Y = level[z]
			weights = Y*0. + 1. ; Define a vector of weights.
			levelmean = mean(Y,/nan)
			A2guess = where(abs(Y - levelmean) eq min(abs(Y - levelmean),/nan))
			A2guess = X[A2guess[0]]
			A = [1000.0,A2guess, levelmean] ; Provide an initial guess of the function's parameters.
			yfit = CURVEFIT(X, Y, weights, A, SIGMA, FUNCTION_NAME='wave2') ; fit data to sine function : y = a[0]*sin(a[1]*x+a[2]) + a[3]
;			PRINT, 'Function parameters: ', A ; Print the parameters returned in A. 
			WACCMheight60[iii,kk] = abs(A[0])
		endif
		      
		level  = reform(WACCMGPHgrid70[*,kk])
      	z = where(level gt 0. and level lt 90000., nz)
   		if nz gt 10. then begin
			X = WACCMlon[z] 
			Y = level[z]
			weights = Y*0. + 1. ; Define a vector of weights.
			levelmean = mean(Y,/nan)
			A2guess = where(abs(Y - levelmean) eq min(abs(Y - levelmean),/nan))
			A2guess = X[A2guess[0]]
			A = [1000.0,A2guess, levelmean] ; Provide an initial guess of the function's parameters.
			yfit = CURVEFIT(X, Y, weights, A, SIGMA, FUNCTION_NAME='wave2') ; fit data to sine function : y = a[0]*sin(a[1]*x+a[2]) + a[3]
;			PRINT, 'Function parameters: ', A ; Print the parameters returned in A. 
			WACCMheight70[iii,kk] = abs(A[0])
		endif
      
      
      
   	endfor ; kk

meanwaccmtemp50 = fltarr(waccm_levels)
meanwaccmtemp60 = fltarr(waccm_levels)
meanwaccmtemp70 = fltarr(waccm_levels)
meanwaccmtemppole = fltarr(waccm_levels)
if iii eq 0L then begin
	WACCMstratopause50 = fltarr(days)
	WACCMstratopause60 = fltarr(days)
	WACCMstratopause70 = fltarr(days)
	WACCMstratopausepole = fltarr(days)
endif
	
x = where(alat gt 47. and alat lt 53.)
for ii = 0L, waccm_levels - 1L do meanwaccmtemp50[ii] = mean(waccmtemp[x,*,ii],/nan)
x = where(alat gt 57. and alat lt 63.)
for ii = 0L, waccm_levels - 1L do meanwaccmtemp60[ii] = mean(waccmtemp[x,*,ii],/nan)
x = where(alat gt 67. and alat lt 73.)
for ii = 0L, waccm_levels - 1L do meanwaccmtemp70[ii] = mean(waccmtemp[x,*,ii],/nan)
x = where(alat ge 70.)
for ii = 0L, waccm_levels - 1L do meanwaccmtemppole[ii] = mean(waccmtemp[x,*,ii],/nan)
	

x = where(pgrid lt 50. and pgrid gt .001 and max(meanWACCMtemp50 gt 100.), nx)
if nx gt 0L then begin
	y = where(max(meanWACCMtemp50[x]) eq meanWACCMtemp50[x],ny)
	if ny gt 0L then WACCMstratopause50[iii] = pgrid[x[y[0]]]
endif

x = where(pgrid lt 50. and pgrid gt .001 and max(meanWACCMtemp60 gt 100.), nx)
if nx gt 0L then begin
	y = where(max(meanWACCMtemp60[x]) eq meanWACCMtemp60[x],ny)
	if ny gt 0L then WACCMstratopause60[iii] = pgrid[x[y[0]]]
endif

x = where(pgrid lt 50. and pgrid gt .001 and max(meanWACCMtemp70 gt 100.), nx)
if nx gt 0L then begin
	y = where(max(meanWACCMtemp70[x]) eq meanWACCMtemp70[x],ny)
	if ny gt 0L then WACCMstratopause70[iii] = pgrid[x[y[0]]]
endif
x = where(pgrid lt 50. and pgrid gt .001 and max(meanWACCMtemppole gt 100.), nx)
if nx gt 0L then begin
	y = where(max(meanWACCMtemppole[x]) eq meanWACCMtemppole[x],ny)
	if ny gt 0L then WACCMstratopausepole[iii] = pgrid[x[y[0]]]
endif
print, WACCMstratopause50[iii]

esdates[iii] = ifile
endfor

save, filename = '/Volumes/MacD68-1/france/ES_paper/Post_process/Platentary_Wave_2_amplitudes.sav',$
	esdates,WACCMheight50,WACCMheight60,WACCMheight70,pgrid,WACCMstratopause50,WACCMstratopause60,WACCMstratopause70,waccmstratopausepole



end
