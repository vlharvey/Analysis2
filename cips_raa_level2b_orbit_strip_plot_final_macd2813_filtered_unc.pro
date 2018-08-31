; Laboratory for Atmospheric and Space Physics
; University of Colorado, Boulder, Colorado, USA
;
; FILENAME:
;   CIPS_raa_orbit_strip_plot.pro
;
; AUTHOR:
;    Jeff France
;
; DATE:    Feb 19, 2018
;
; PURPOSE:
;   This program reads in Level 2b RAA files and creates orbit strip plots and corresponding thumbnails
;  
; NOTES:
;    

;Identify orbits to loop over
startorbit = 48633
endorbit = 50729

RAA_Level2_path='/atmos/harvey/CIPS_data/Datfiles/RAA/'
outpath = 'Figures/'


;Read in MERRA2 file for determining land/water using surface GPHT
file=RAA_Level2_path+'MERRA2_300.inst3_3d_asm_Np.20040101.SUB.nc4'
ncfile = file_search(file,count = count);Open the *.nc file
ncid=ncdf_open(ncfile)
result=ncdf_inquire(ncid)   ;Inquire about the data
nvars=result.nvars        ;# variables in the file
;Read in the data
for ivar=0,nvars-1 do begin
		result=ncdf_varinq(ncid,ivar) ;get the data name, type, dimensions
		;Puts data into array called "data" and variable name into "result.name":
		ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
		
	 ; Invoke some magic to check for string data
	 ; masquerading as byte data, but don't convert
	 ; byte data blindly, i.e., quality_flags is a 2-dimensional
	 ; array of byte data. 
	 ; (This is done to account for a bug in the ncdf write routine.)
		  if ( ( size( data, /n_dimensions )  EQ 1 ) && $
			( size( data, /type ) EQ 1 ) ) then $
			data = string( data )

	 ;Extract each variable from the "data" structure and name it
	;  the corresponding "name" from "result.name":
		   if Execute(result.name + ' = data') eq 0 then $
		   Print, ' "Execute" command failed -- are you in Virtual Machine mode?'            
endfor
NCDF_CLOSE, ncid

;Assume surface GPHT less than 1m and greater than -1m is water
x = where(abs(mean(phis ,dim=3)) gt .1)
continent_mask = reform(phis[*,*,0])*0.
continent_mask[x] = 1.
continent_mask = smooth(continent_mask,5,/nan)
lon_mask = lon
x = where(lon_mask lt 0.)
lon_mask[x] = lon_mask[x] + 360.
lat_mask = lat


a=findgen(8)*(2*!pi/8.)
usersym,1.5*cos(a),1.5*sin(a),/fill


orbitList = findgen(endorbit-startorbit)+startorbit
strOrbits = string(orbitlist, format = '(i5.5)')
;LOOP OVER ORBITS
for iorbit = 0, n_elements(orbitlist) - 1L do begin ; loop over orbits

;set_plot,'ps'
;device,/inch,xoff=0,yoff=0.,xsize=xsize,ysize=ysize,/encapsulate,$;
; /bold,/color,bits_per_pixel=8,/helvetica,filename=RAA_Plot_filename_ps,decompose=0

	print, strorbits[iorbit]
	;---------------------READ IN CIPS DATA USING READ_CIPS_NCDF-------------------------------------------
	file=RAA_Level2_path+'cips_raa_2b_orbit_'+strOrbits[iorbit]+'*v01.10_r04_alb.nc'
	ncfile = file_search(file,count = count);Open the *.nc file
	if count ne 1 then continue
	 ncid=ncdf_open(ncfile)
	 result=ncdf_inquire(ncid)   ;Inquire about the data
	 nvars=result.nvars        ;# variables in the file
	 ;Read in the data
	for ivar=0,nvars-1 do begin
			result=ncdf_varinq(ncid,ivar) ;get the data name, type, dimensions
			;Puts data into array called "data" and variable name into "result.name":
			ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
			
		 ; Invoke some magic to check for string data
		 ; masquerading as byte data, but don't convert
		 ; byte data blindly, i.e., quality_flags is a 2-dimensional
		 ; array of byte data. 
		 ; (This is done to account for a bug in the ncdf write routine.)
			  if ( ( size( data, /n_dimensions )  EQ 1 ) && $
				( size( data, /type ) EQ 1 ) ) then $
				data = string( data )

		 ;Extract each variable from the "data" structure and name it
		;  the corresponding "name" from "result.name":
			   if Execute(result.name + ' = data') eq 0 then $
			   Print, ' "Execute" command failed -- are you in Virtual Machine mode?'            
	endfor
	NCDF_CLOSE, ncid

	file=RAA_Level2_path+'cips_raa_2b_orbit_'+strOrbits[iorbit]+'*v01.10_r04_cat.nc'
	ncfile = file_search(file,count = count);Open the *.nc file
	if count ne 1 then continue
	ncid=ncdf_open(ncfile)
	result=ncdf_inquire(ncid)   ;Inquire about the data
	nvars=result.nvars        ;# variables in the file
	;Read in the data
	for ivar=0,nvars-1 do begin
		result=ncdf_varinq(ncid,ivar) ;get the data name, type, dimensions
		;Puts data into array called "data" and variable name into "result.name":
		ncdf_varget,ncid,ncdf_varid(ncid,result.name),data
		
		; Invoke some magic to check for string data
		; masquerading as byte data, but don't convert
		; byte data blindly, i.e., quality_flags is a 2-dimensional
		; array of byte data.
		; (This is done to account for a bug in the ncdf write routine.)
		if ( ( size( data, /n_dimensions )  EQ 1 ) && $
			( size( data, /type ) EQ 1 ) ) then $
			data = string( data )
			
		;Extract each variable from the "data" structure and name it
		;  the corresponding "name" from "result.name":
		if Execute(result.name + ' = data') eq 0 then $
			Print, ' "Execute" command failed -- are you in Virtual Machine mode?'
	endfor
	NCDF_CLOSE, ncid
	
	;-----------------------------------------------------------------------

raa_unc=RAYLEIGH_ALBEDO_ANOMALY_unc
err_frac=abs(raa_unc)

	;Filter out high SZA
	x=where(ZENITH_ANGLE ge 95. or finite(RAYLEIGH_ALBEDO_ANOMALY) eq 0., nx)
	RAA_outline = RAYLEIGH_ALBEDO_ANOMALY
	RAYLEIGH_ALBEDO_ANOMALY[x] = !values.f_nan
	x=where(finite(raa_outline))
	raa_outline[x] = 1.

	

	;redefine plot variables and shift longitude to 0-360
	plot_RAA = RAYLEIGH_ALBEDO_ANOMALY
	plot_lon = longitude
	plot_lat = latitude
        x = where(plot_lon lt 0.)
        plot_lon[x] = plot_lon[x] +360.
        x = where(longitude lt 0.)
        longitude[x] = longitude[x] +360.

	
	;Determine orbit start time in UT hours:minutes
	x=where(finite(ut_time))                  
	UT = ut_time[x[0]] - floor(ut_time[x[0]])
	hours = strmid(ORBIT_START_TIME_UT,9,2);  STRING    = '2017/266-06:59:44'
	minutes = strmid(ORBIT_START_TIME_UT,12,2)

	plot_title = 'CIPS RAA (%) -- '+' Orbit Number = '+strOrbits[iorbit]+', Date = '+UT_DATE_ORBIT_START+' '+hours+':'+minutes+'Z'
	RAA_Plot_filename = outpath+'RAA_Level2b_v01.1_rev04_Orbit-'+strOrbits[iorbit]+'_'+UT_DATE_ORBIT_START+'_contour_filtered.png'
	RAA_Plot_filename_60 = outpath+'small_RAA_Level2b_v01.1_rev04_Orbit-'+strOrbits[iorbit]+'_'+UT_DATE_ORBIT_START+'_contour_filtered.png'
	RAA_Plot_filename_thumb = outpath+'RAA_Level2b_Orbit-'+strOrbits[iorbit]+'_'+UT_DATE_ORBIT_START+'.thumb_filtered.png'

	
	; determine whether each CIPS pixel is over land or water using MERRA2 surface GPHT for plotting continents
	continent_plot = RAYLEIGH_ALBEDO_ANOMALY*0.	
	for ii = 0, n_elements(RAYLEIGH_ALBEDO_ANOMALY) - 1L do begin
		xlat = latitude[ii]
		xlon = longitude[ii]
		x =min(abs(lat_mask - xlat),latindex)
		x =min(abs(lon_mask - xlon),lonindex)
		continent_plot[ii] = continent_mask[lonindex[0],latindex[0]]
	endfor


; add white box around title text

	;Define linear colortable	

  ;create generic plotting arrays
	xarray = findgen(n_elements(RAYLEIGH_ALBEDO_ANOMALY[*,0]))
	yarray = findgen(n_elements(RAYLEIGH_ALBEDO_ANOMALY[0,*]))

	x_y_ratio = (n_elements(xarray)*1.)/n_elements(yarray)
	if x_y_ratio gt 25. then x_y_ratio = 25.
 
  thisDevice = !D.Name
	Set_Plot, 'Z'
	;device,/inch,xoff=0,yoff=0.,xsize=xsize,ysize=ysize,$;
	;   filename=RAA_Plot_filename,decompose=1
	xsize = n_elements(xarray)
	if xsize lt 1000. then xsize = 1000.
  charsizefactor =(1.9 - .4*2500./xsize) * .6
  charthickfactor =(1.9 - .6*2500./xsize) * .43


  xsize = xsize*1.5
	ysize = xsize/x_y_ratio *1.1
	if ysize lt 200 then ysize = 200.

  xsize = xsize*.6
  ysize = ysize*.6
  
 
  ;hardwire for wide orbits
	if n_elements(xarray) gt 3000 and n_elements(yarray) gt 500 then charsizefactor = 1.5
	if n_elements(xarray) gt 3000 and n_elements(yarray) gt 500 then charthickfactor = 1.3

	Device, Decomposed=0, Set_Pixel_Depth=24, Set_Resolution=[xsize,ysize]
;	Set_Plot, thisDevice
;	!p.charsize=1.4    ; test with 1.5
;	!p.thick=8
;	!p.charthick=3.1
;	!y.thick=2
;	!x.thick=2
;	!p.font=0
	erase

loadct,0
  contour,continent_plot,xarray, yarray,/fill,level = [0,.85], c_color = [0,150],$
  position = [.01,.01,.93,.84], xstyle=5,ystyle=5, color = 0, /nodata,background = 0
  x = where(finite(continent_plot) eq 0,comp = y)
  continent_plot1 =continent_plot
  continent_plot1[x] = 1.
  continent_plot1[y] = 0.
  contour,continent_plot1,xarray, yarray,/fill,level = [0], c_color = [0],/overplot

	;Determine if orbit is during PMC season and add red shading poleward of 60
	year = strmid(UT_DATE_ORBIT_START,0,4)
	month = strmid(UT_DATE_ORBIT_START,4,2)
	day = strmid(UT_DATE_ORBIT_START,6,2)
	current_doy = julday(month,day,year)-julday('01','01',year)+1
	if current_doy ge 321 or current_doy le 54 then SH_PMC = 1 else SH_PMC = 0
	if current_doy ge 135 and current_doy le 241 then NH_PMC = 1 else NH_PMC = 0

;	loadct, 70
;	if NH_PMC eq 1 then begin
;		x = where(latitude ge 60.,comp = y)
;		continent_plot2 = continent_plot*0. + 1.
;		continent_plot2[y] = 0
;		contour,continent_plot2,xarray, yarray,/fill,level = [1], c_color = [80],/overplot
;	endif
;	if SH_PMC eq 1 then begin
;		x = where(latitude le -60.,comp = y)
;		continent_plot2 = continent_plot*0. + 1.
;		continent_plot2[y] = 0
;		contour,continent_plot2,xarray, yarray,/fill,level = [1], c_color = [80],/overplot
;	endif
	
	
  ;add continents	
	loadct,0
	continent_plot = smooth(continent_plot,5,/edge_truncate)
	contour,continent_plot,xarray, yarray,/fill,level = [.85], c_color = [150],/overplot

	;Add lat/lon contours	
	x = where(longitude le 5.)
	longitude[x] = !values.f_nan
	x=where(longitude gt 359.0)
	longitude[x] = 360.
	
	x = where(abs(latitude) gt 70.)
	xlons = longitude
	xlons[x] = !values.f_nan
	
	x=where(finite(rayleigh_albedo_anomaly))
	latitude[x]=!values.f_nan
	longitude[x]=!values.f_nan
	xlons[x] = !values.f_nan
  contour, latitude, xarray,yarray,levels = [-80,-60,-50,-40,-30,-20,-10,10,20,30,40,50,60,80],c_linestyle = 0, /follow, /overplot,color = 250,$
     c_charsize = charsizefactor*1.5, c_thick = 1.5,C_Labels=Replicate(1, 13),c_charthick = charsizefactor*2
  contour, latitude, xarray,yarray,levels = [-90,-80,-70,-60,-50,-40,-30,-20,-10,10,20,30,40,50,60,70,80,90],c_linestyle = 0, /overplot,color = 250,$
       c_charsize = charsizefactor*1.5, c_thick = 1.5,C_Labels=Replicate(0, 19),c_charthick = charsizefactor*2
  contour, latitude, xarray,yarray,levels = [0],c_linestyle = 0, /follow, /overplot,color = 250, c_charsize = charsizefactor*1.5, c_thick = 3,C_Labels=[0],C_annot = ['Eq'],c_charthick = charsizefactor*2
  contour, xlons, xarray,yarray,levels = findgen(37)*10.,c_linestyle = 0, /follow, /overplot,color = 250, c_charsize = charsizefactor*1.5, c_thick =1.5,c_charthick = 2
  contour, longitude, xarray,yarray,levels = findgen(37)*10.,c_linestyle = 0, /overplot,color = 250, c_charsize = charsizefactor*1.5, c_thick =1.5,C_Labels=Replicate(0, 36),c_charthick = charsizefactor*2
  contour, longitude, xarray,yarray,levels = [360],c_linestyle = 0, /overplot,color = 250, c_charsize = charsizefactor*1.5, c_thick =1.5,C_Labels=Replicate(0, 1),c_charthick = charsizefactor*2
  contour, longitude, xarray,yarray,levels = [0],c_linestyle = 0, /overplot,color = 250, c_charsize = charsizefactor*1.5, c_thick =1.5,C_Labels=Replicate(0, 1),c_charthick = charsizefactor*2
  ;contour,continent_plot,xarray, yarray,/fill,level = [0,.85], c_color = [0,150],position = [.01,.01,.9,.9], xstyle=5,ystyle=5, color = 0, /nodata, /noerase

	
	;loop over CIPS pixels and plot using psym
;	loadct,60
;	for i=0, n_elements(rayleigh_albedo_anomaly[*,0]) - 1L do begin
;	for j=0, n_elements(rayleigh_albedo_anomaly[0,*]) - 1L do begin
;		dum = min(abs(levels1 - rayleigh_albedo_anomaly[i,j]),icol)
; testing... tweak symbol size/thickness.....
;		if finite(dum) then oplot, [xarray[i],xarray[i]],[yarray[j],yarray[j]],color = col1[icol], psym=8,symsize = .3, thick = 2
;		if finite(dum) then oplot, [xarray[i],xarray[i]],[yarray[j],yarray[j]],color = col1[icol], psym=8,symsize = .23, thick = 2
;		if finite(dum) then oplot, [xarray[i],xarray[i]],[yarray[j],yarray[j]],color = col1[icol], psym=8,symsize = .14, thick = 1
; testing....

;	endfor
;	endfor

;Contour image
raa_2std = 2.*stddev(rayleigh_albedo_anomaly,/nan)
if sh_pmc eq 1 then begin
	x = where(finite(plot_lat) and plot_lat ge -60.,comp=y)
	raa_2std = 2.*stddev(rayleigh_albedo_anomaly[x],/nan)
	pmc_2std = 2.*stddev(rayleigh_albedo_anomaly[y],/nan)
endif	
if nh_pmc eq 1 then begin
	x = where(finite(plot_lat) and plot_lat le 60.,comp=y)
	raa_2std = 2.*stddev(rayleigh_albedo_anomaly[x],/nan)
	pmc_2std = 2.*stddev(rayleigh_albedo_anomaly[y],/nan)
endif	

	if raa_2std le 3. and raa_2std gt 0.5 then begin
	levels1 = (findgen(62)*(2.*raa_2std/60.)-raa_2std)
	endif else 	levels1 = (findgen(62)*.1-3.)
nlvls  = n_elements(levels1)
col1 = reverse((1 + indgen(nlvls)) * 250. / nlvls)	; define colors
col1[0] = 255.
col1[-1] = 0.


x=where(rayleigh_albedo_anomaly lt min(levels1) and finite(rayleigh_albedo_anomaly),nx)
if nx gt 0 then rayleigh_albedo_anomaly[x] = min(levels1)

loadct,0
		x = where(finite(raa_outline))
		datalocations = fltarr(n_elements(xarray),n_elements(yarray))
		datalocations[x] =10.
				;Plot outline of CIPS scene
		contour, smooth(datalocations,1,/nan),xarray,yarray,level=[.1], /overplot,/noerase,/close,c_color = 250, c_thick = 2*charsizefactor

;contour, smooth(raa_outline,3,/nan), xarray,yarray,/overplot, levels = [.99],c_color= 0,c_thick=charsizefactor*3

raa_pmc = rayleigh_albedo_anomaly

err_thresh=.7

if sh_pmc eq 1 then begin
  x = where(finite(plot_lat) and plot_lat le -60.,comp=y)
    raa_pmc[y] = !values.f_nan
    loadct,62
	if pmc_2std le 3. and pmc_2std gt 0.5 then begin
	levels2 = (findgen(62)*(2.*pmc_2std/60.)-pmc_2std)
	endif else 	levels2 = (findgen(62)*.1-3.)
	nlvls2  = n_elements(levels2)
	col2 = reverse((1 + indgen(nlvls2)) * 250. / nlvls2)	; define colors
	col2[0] = 255.
	col2[-1] = 0.
index=where(err_frac gt err_thresh)    ; where uncertainty is larger than |RAA|
if index(0) ne -1L then raa_pmc(index)=0./0.
if index(0) ne -1L then rayleigh_albedo_anomaly(index)=0./0.

    contour, raa_pmc, xarray,yarray,/overplot, levels = levels2, c_color = col2, /cell_fill
    rayleigh_albedo_anomaly[x] = !values.f_nan
    ;    xarray2d[x] = !values.f_nan
;    yarray2d[x] = !values.f_nan
    loadct,60
  contour, rayleigh_albedo_anomaly, xarray,yarray,/overplot, levels = levels1, c_color = col1, /cell_fill 

endif else begin
  if nh_pmc eq 1 then begin
  x = where(finite(plot_lat) and plot_lat ge 60., comp = y)
    raa_pmc[y] = !values.f_nan
    loadct,62
	if pmc_2std le 3. and pmc_2std gt 0.5 then begin
	levels2 = (findgen(62)*(2.*pmc_2std/60.)-pmc_2std)
	endif else 	levels2 = (findgen(62)*.1-3.)
	nlvls2  = n_elements(levels2)
	col2 = reverse((1 + indgen(nlvls2)) * 250. / nlvls2)	; define colors
	col2[0] = 255.
	col2[-1] = 0.

index=where(err_frac gt err_thresh)    ; where uncertainty is larger than |RAA|
if index(0) ne -1L then raa_pmc(index)=0./0.
if index(0) ne -1L then rayleigh_albedo_anomaly(index)=0./0.

    contour, raa_pmc, xarray,yarray,/overplot, levels = levels2, c_color = col2, /cell_fill
    rayleigh_albedo_anomaly[x] = !values.f_nan

    loadct,60
  contour, rayleigh_albedo_anomaly, xarray,yarray,/overplot, levels = levels1, c_color = col1, /cell_fill 

  endif else begin
      loadct,60

index=where(err_frac gt err_thresh)    ; where uncertainty is larger than |RAA|
if index(0) ne -1L then raa_pmc(index)=0./0.
if index(0) ne -1L then rayleigh_albedo_anomaly(index)=0./0.

      contour, rayleigh_albedo_anomaly, xarray,yarray,/overplot, levels = levels1, c_color = col1, /cell_fill
  endelse
endelse
loadct,0
	plot,[0,0],[0,1],xrange=[0,1],yrange=[0,1],/noeras,/nodata,xstyle=4,ystyle = 4,position = [.01,.84, .9,.99]
	xbox=[0,1,1,0,0]
	ybox=[0,0,1,1,0]
	polyfill,xbox,ybox,color=250
	!type=2^2+2^3+2^6     ; no y title or ticsks
	plot, [0,0],[0,0],xstyle = 4, ystyle =  4,position = [.01,.84,.93,.99], charsize = 1.125,thick = 4,/noerase ; style= 4 supresses axis
	xyouts, .5,.2,plot_title, align = .5, charsize = charsizefactor*2.5, color = 0, charthick = charthickfactor*5.
	; -----------------plot the color bar-----------------------	
	slab=' '+strarr(n_elements(levels1))

	loadct,0
	!p.title = ' '
	plot,[0,0],[0,1],xrange=[0,1],yrange=[0,1],/noeras,/nodata,$
	xstyle=4,ystyle = 4,position = [.93,0.01, .999,.84]
		xbox=[0,1,1,0,0]
		ybox=[0,0,1,1,0]
		polyfill,xbox,ybox,color=250
		!type=2^2+2^3+2^6     ; no y title or ticsks
		plot, [0,0],[0,0],xstyle = 4, ystyle =  4,$
		position = [.01,.01,.99,.84], charsize = charsizefactor*1.125,thick = 4,/noerase ; style= 4 supresses axis

	;Add the red color table if necessary
	if nh_pmc eq 1 or sh_pmc eq 1 then begin
		xyouts, .99,.5,'Rayleigh Albedo Anomaly (%)', orientation = 90, align = .5, charsize =charsizefactor*1., color = 0, charthick = charthickfactor*1.7
		position = [.93,0.01, .95,.84]
		plot,[0,0],[0,1],xrange=[0,10],yrange=[0,1],/noeras,$
		position = [.95,.01,.97,.84],xstyle=1, color = 250
		loadct,62
		xbox=[0,10,10,0,0]
		x2=0
		for j=1,n_elements(col2)-1 do begin
			dx= 1./(n_elements(levels2)-1.)
			ybox=[x2,x2,x2+dx,x2+dx,x2]
			polyfill,xbox,ybox,color=col2[j-1]
			x2=x2+dx
		endfor
		loadct,0
		slab=strcompress(string(format='(f8.3)',levels2),/remove_all)
		slabcolor = fltarr(n_elements(levels2))*0.
		;slabcolor[0:2] = 255	
	;	slabcolor[-4:-1] = 255
		x1=dx/2+.01
	for i=0L,n_elements(slab)-2L,3 do begin
		slab0=slab[i]
		flab0=float(slab[i])
		slab0=strcompress(string(format='(f5.1)',flab0),/remove_all)
		xyouts,5,x1-dx/2.,slab0,charsize=charsizefactor*.8,/data,color=slabcolor[i], align = .5,charthick=charthickfactor*1.3
		x1=x1+.96*dx*3
	endfor

	endif

	if nh_pmc eq 0 and sh_pmc eq 0 then begin
		xyouts, .98,.5,'Rayleigh Albedo Anomaly (%)', orientation = 90, align = .5, charsize =charsizefactor*1., color = 0, charthick = charthickfactor*2
		position = [.93,0.01, .96,.84]
	endif
	;Add the purple color table
	plot,[0,0],[0,1],xrange=[0,10],yrange=[0,1],/noeras,$
		position = position,xstyle=1, color = 250
	loadct,60
	xbox=[0,10,10,0,0]
 	x2=0
	for j=1,n_elements(col1)-1 do begin
		dx= 1./(n_elements(levels1)-1.)
		ybox=[x2,x2,x2+dx,x2+dx,x2]
		polyfill,xbox,ybox,color=col1[j-1]
		x2=x2+dx
	endfor
	loadct,0
	slab=strcompress(string(format='(f8.3)',levels1),/remove_all)
	slabcolor = fltarr(n_elements(levels1))*0.
	;slabcolor[0:2] = 255	
;	slabcolor[-4:-1] = 255
	x1=dx/2+.01
	
	for i=0L,n_elements(slab)-2L,3 do begin
		slab0=slab[i]
		flab0=float(slab[i])
		slab0=strcompress(string(format='(f5.1)',flab0),/remove_all)
		xyouts,5,x1-dx/2.,slab0,charsize=charsizefactor*.8,/data,color=slabcolor[i], align = .5,charthick=charthickfactor*1.3
		x1=x1+.96*dx*3
	endfor




	; convert .ps to png and save image and thumbnail
	;device, /close
  write_png,RAA_Plot_filename ,tvrd(channel=1,/true)
 ; spawn,'convert '+ RAA_Plot_filename+' -resize 18% '+RAA_Plot_filename_thumb
 ; spawn,'convert '+ RAA_Plot_filename+' -resize 50% '+RAA_Plot_filename_60
endfor
end
